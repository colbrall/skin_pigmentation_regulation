# regulation_over_time.jl
# given ensembl id(s) and a set of individuals, calculates and plots regression for predictions over time
# @author Laura Colbran
# julia 1.4

using ArgParse
using GZip, DataFrames, CSV
using GLM, Distributions
using Seaborn, Plots

SKIN_GENES = ["ENSG00000161091","ENSG00000167986","ENSG00000104044",
"ENSG00000128731","ENSG00000188467","ENSG00000077498","ENSG00000107165",
"ENSG00000164175","ENSG00000258839","ENSG00000101440","ENSG00000115138",
"ENSG00000080166","ENSG00000049130","ENSG00000157404","ENSG00000137265",
"ENSG00000176797","ENSG00000185664","ENSG00000162341","ENSG00000187098",
"ENSG00000197535","ENSG00000047579","ENSG00000069974","ENSG00000088812",
"ENSG00000143669","ENSG00000115648","ENSG00000166189","ENSG00000134160",
"ENSG00000151694","ENSG00000173157","ENSG00000146648","ENSG00000145244",
"ENSG00000112038","ENSG00000157168","ENSG00000173068","ENSG00000040531",
"ENSG00000136160","ENSG00000124205"]

# parses command-line arguments
function parseCommandline()
	s = ArgParseSettings()
	# @add_arg_table! s begin
	@add_arg_table! s begin
		"--regulation", "-r"
			nargs='*'
			help="path(s) to predicted regulation files"
			arg_type = String
			required = true
		"--genes","-g"
			help="file with gene IDs to track"
			arg_type=String
			required = true
		"--column","-c"
			help = "1-indexed column containing gene IDs to search for"
			arg_type=Int64
			required = true
        "--samples","-s"
            help="path to file with samples to track over time."
            arg_type=String
			required = true
        "--date_col","-d"
            help = "1-indexed column in sample file that contains date"
            arg_type=Int64
			required = true
		"--modern_pop","-p"
			nargs='*'
			help="path(s) to predicted regulation files for a modern population to include in plotting."
			arg_type = String
			default=["None"]
		"--modern_samples","-m"
			help="path to file with modern samples to include, if you're filtering them. Ignored for regression calculation"
			arg_type=String
			default="None"
		"--covariates","-v"
			help = "file with covariates to regress out (e.g. PCs). expects all individuals in one file."
			arg_type=String
			default="None"
		"--cov_cols","-l"
			help = "columns from covariates file to add. assumes first given is ind_id"
			nargs = '*'
			arg_type =Int64
		"--out_dir", "-o"
            help = "directory path to write output files"
            arg_type = String
            default = "./"
		"--plot"
			help = "if you want to plot the scatterplot and regression for each model"
			action=:store_true
	end
	return parse_args(s)
end

function sampDF(s_path::String,col::Int64)
    df = CSV.read(s_path, DataFrame; delim='\t')
	return df[:,[1,col]]
end

function geneArray(gene_path::String,gene_col::Int64)
	genes = String[]
	open(gene_path) do f
		for line in eachline(f)
			if startswith(line, "#") continue end
			push!(genes,split(chomp(line),"\t")[gene_col])
		end
	end
	return genes
end

# returns tissue name from prediction file
function getTissue(file::String)::String
    return join(split(basename(file),"")[1:end-24])
end

#returns list of IDs, and predictions for particular gene
function pullModern(pop_path::Array{String,1},s_path::String,gene_ids::Array{String,1},tiss::String)
	preds = DataFrames.DataFrame("ind_id" => String[],"date"=>Int64[])
    indices = Dict{SubString,Array{Int64,1}}()
	for t in pop_path
		t_name = getTissue(tiss)
		if getTissue(t) == t_name
			GZip.open("$(t)") do f
		        for line in eachline(f)
		            l = split(chomp(line),'\t')
		            if startswith(line,"gene")
		                if s_path != "None"
							tmp = CSV.read(s_path,DataFrame;delim='\t',header=[:sample,:pop])
							preds = DataFrames.DataFrame("ind_id" => tmp[!,:sample], "date"=>fill(0,nrow(tmp)))
						else
							preds = DataFrames.DataFrame("ind_id" => l[2:end], "date"=>fill(0,length(l)-1))
						end
						indices = [findfirst(x->x==i,l) for i in preds[:,1]]
		                continue
		            end
		            if !in(l[1], gene_ids) continue end
		            preds[!,Symbol("$(l[1])_$t_name")] = parse.(Float64,l[indices])
				end
			end
		end
	end
	return preds
end
# calculates linear regression across time
function timeSeries(pred_path::Array{String,1},gene_path::String,gene_col::Int64,s_path::String,mod_pop::Array{String,1},mod_s::String,s_col::Int64,covariates::String,cov_cols::Array{Int64,1},out_path::String,plot::Bool)
	var = sampDF(s_path,s_col) # ends up with cols sample, date, [genes x tissue]
	id_name = names(var)[1]
	date_name = names(var)[2]
	gene_ids = geneArray(gene_path,gene_col)
    indices = Dict{SubString,Array{Int64,1}}()
	mod_preds = DataFrames.DataFrame("ind_id" => String[],"date"=>Int64[])
	covs = DataFrames.DataFrame("ind_id" => String[])
	results = DataFrames.DataFrame("model"=> String[],"coeff" => Float64[],"stderr"=> Float64[],"t"=> Float64[],"pval"=>Float64[],"low95"=> Float64[],"up95"=> Float64[])
	for tiss in pred_path
		t_name = getTissue(tiss)
		if mod_pop != ["None"] # add modern individuals to var if present
			if nrow(mod_preds) > 0
				tmp = pullModern(mod_pop,mod_s,gene_ids,tiss)
				for id in gene_ids
					mod_preds[Symbol("$(id)_$t_name")] = tmp[Symbol("$(id)_$t_name")] # just add new column(s)
				end
			else
				mod_preds = pullModern(mod_pop,mod_s,gene_ids,tiss)
			end
		end
		if covariates != "None"
			# read in covariates file
		    covs = CSV.read(covariates, DataFrame; delim='\t')[!,cov_cols]
			new = ["ind_id"]
			for i in 1:(ncol(covs)-1) push!(new,"cov$(i)") end
			rename!(covs,new)
		end
	    GZip.open("$(tiss)") do f
	        for line in eachline(f)
	            l = split(chomp(line),'\t')
	            if startswith(line,"gene")
	                indices = [findfirst(x->x==i,l) for i in var[:,1]]
	                continue
	            end
	            if !in(l[1], gene_ids) continue end
	            var[!,Symbol("$(l[1])_$t_name")] = parse.(Float64,l[indices])
				# calculate regression
				if mod_pop != ["None"]
					tmp = DataFrames.DataFrame(:ind_id => vcat(var[!,Symbol("$id_name")],mod_preds[!,:date]),
							:expr => vcat(var[!,Symbol("$(l[1])_$t_name")],mod_preds[!,Symbol("$(l[1])_$t_name")]),
							:date => vcat(var[!,Symbol("$date_name")],mod_preds[!,:date]))
				else
					tmp = DataFrames.DataFrame(:ind_id => vcat(var[!,Symbol("$id_name")],mod_preds[!,:date]),
							:expr => var[!,Symbol("$(l[1])_$t_name")],
							:date => var[!,Symbol("$date_name")])
				end
				if covariates != "None"
					# join covariates to tmp
					tmp = innerjoin(tmp,covs, on = :ind_id)
				end
				fm = @eval (Term(:expr) ~ sum(term.($(Symbol.(names(tmp)[3:end])))))
				model = coeftable(lm(fm,tmp))#fit regression
				tmp = tmp[tmp[:,:date] .< 15000,:]
				model_15ky = coeftable(lm(fm,tmp))#fit regression
				# println(model)
				# push!(results,["$(l[1])_$(t_name):",model.cols[1][2],model.cols[2][2],model.cols[3][2],model.cols[4][2],model.cols[5][2],model.cols[6][2]])
				push!(results,["$(l[1])",model.cols[1][2],model.cols[2][2],model.cols[3][2],model.cols[4][2],model.cols[5][2],model.cols[6][2]])

				#plot date vs pred. expression
				if plot
					s_plot = Seaborn.regplot(x=vcat(var[!,Symbol("$date_name")],mod_preds[!,:date]),
								y=vcat(var[!,Symbol("$(l[1])_$t_name")],mod_preds[!,Symbol("$(l[1])_$t_name")]),
								fit_reg=:false)
	                s_plot.set_title("$(l[1])_$t_name")
	                s_plot.set_ylabel("Pred. Norm. Expr.")
					s_plot.set_xlabel("Age (yBP)")

					x_vals = collect(1:s_plot.get_xlim()[2])
					x_vals_15k = collect(1:15000)
				    y_vals = x_vals .* model.cols[1][2] .+ model.cols[1][1]
					y_vals_15ky = x_vals_15k .* model_15ky.cols[1][2] .+ model_15ky.cols[1][1]
					Seaborn.regplot(x=x_vals_15k,y=y_vals_15ky,color=:black,scatter = :false)
				    Seaborn.regplot(x=x_vals, y=y_vals,color = :red,scatter =:false)
	                Seaborn.savefig("$(out_path)$(l[1])_$(t_name)_time_series.png") #N.B. on ancient EUR pdf has so many oints inkscape won't open it
	                clf()
				end
	        end
	    end
	end
	CSV.write("$(out_path)regr_stats.txt",results;delim = "\t") #qqplot can be made from this output using qqplot.jl

end

function main()
    parsed_args = parseCommandline()

    timeSeries(parsed_args["regulation"],parsed_args["genes"],parsed_args["column"],
			parsed_args["samples"],parsed_args["modern_pop"],parsed_args["modern_samples"],
			parsed_args["date_col"],parsed_args["covariates"],parsed_args["cov_cols"],
			parsed_args["out_dir"],parsed_args["plot"])
end

main()

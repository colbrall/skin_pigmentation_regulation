# regulation_over_time.jl
# given ensembl id(s) and a set of individuals, calculates and plots regression for predictions over time
# @author Laura Colbran
# julia 1.4

using ArgParse
using GZip, DataFrames, CSV
using GLM
using Seaborn

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
			help = "column(s) in sample file to regress out (e.g. PCs)"
			nargs='*'
			arg_type=Int64
		"--out_dir", "-o"
            help = "directory path to write output files"
            arg_type = String
            default = "./"
	end
	return parse_args(s)
end

function sampDF(s_path::String,col::Int64,cov_cols::Array{Int64,1})
    df = CSV.read(s_path, DataFrame; delim='\t')
	if length(cov_cols) > 0
		return df[:,[1,col,cov_cols]]
	else
    	return df[:,[1,col]]
	end
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
function timeSeries(pred_path::Array{String,1},gene_path::String,gene_col::Int64,s_path::String,mod_pop::Array{String,1},mod_s::String,s_col::Int64,cov_cols::Array{Int64,1},out_path::String)
	var = sampDF(s_path,s_col,cov_cols) # ends up with cols sample, date, [covs], [genes x tissue]
	date_name = names(var)[2]
	gene_ids = geneArray(gene_path,gene_col)
    indices = Dict{SubString,Array{Int64,1}}()
	mod_preds = DataFrames.DataFrame("ind_id" => String[],"date"=>Int64[])
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
	    GZip.open("$(tiss)") do f
	        for line in eachline(f)
	            l = split(chomp(line),'\t')
	            if startswith(line,"gene")
	                indices = [findfirst(x->x==i,l) for i in var[:,1]]
	                continue
	            end
	            if !in(l[1], gene_ids) continue end
	            var[!,Symbol("$(l[1])_$t_name")] = parse.(Float64,l[indices])
				#plot date vs pred. expression
				s_plot = Seaborn.regplot(vcat(var[!,Symbol("$date_name")],mod_preds[!,:date]),
							vcat(var[!,Symbol("$(l[1])_$t_name")],mod_preds[!,Symbol("$(l[1])_$t_name")]))
                s_plot.set_title("$(l[1])_$t_name")
                s_plot.set_ylabel("Pred. Norm. Expr.")
				s_plot.set_xlabel("Age (yBP)")
                Seaborn.savefig("$(out_path)$(l[1])_$(t_name)_time_series.pdf")
                clf()
				if length(cov_cols) > 0
					#fit regression with covariates as well as date
					println("covariate regression not implemented yet")
					exit()
				else
					tmp = DataFrames.DataFrame(:x => vcat(var[!,Symbol("$date_name")],mod_preds[!,:date]),
							:y => vcat(var[!,Symbol("$(l[1])_$t_name")],mod_preds[!,Symbol("$(l[1])_$t_name")]))
					fm = @eval @formula($(:y) ~ $(:x))
					model = lm(fm,tmp)#fit regression
					println("$(l[1])_$(t_name):")
					println(model)
				end
	        end
	    end
	end
end

function main()
    parsed_args = parseCommandline()

    # println("Parsed args:")
    # for (arg,val) in parsed_args
    #     println("  $arg  =>  $val")
    # end
    timeSeries(parsed_args["regulation"],parsed_args["genes"],parsed_args["column"],parsed_args["samples"],parsed_args["modern_pop"],parsed_args["modern_samples"],parsed_args["date_col"],parsed_args["covariates"],parsed_args["out_dir"])
end

main()

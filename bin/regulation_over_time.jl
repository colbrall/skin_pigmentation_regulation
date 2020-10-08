# regulation_over_time.jl
# given ensembl id(s) and a set of individuals, calculates and plots regression for predictions over time
# @author Laura Colbran
# julia 1.4

using ArgParse
using GZip, DataFrames, CSV
using GLM
using Plots

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
            help = "column in sample file that contains date"
            arg_type=Int64
			required = true
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
    df = CSV.read(s_path; delim='\t',allowmissing=:none)
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

# calculates linear regression across time
function timeSeries(pred_path::Array{String,1},gene_path::String,gene_col::Int64,s_path::String,s_col::Int64,cov_cols::Array{Int64,1},out_path::String)
	var = sampDF(s_path,s_col,cov_cols) # ends up with cols sample, date, [covs], [genes x tissue]
	date_name = names(var)[2]
	gene_ids = geneArray(gene_path,gene_col)
    indices = Dict{SubString,Array{Int64,1}}()
	for tiss in pred_path
		t_name = getTissue(tiss)
	    GZip.open("$(tiss)") do f
	        for line in eachline(f)
	            l = split(chomp(line),'\t')
	            if startswith(line,"gene")
	                indices = [findfirst(x->x==i,l) for i in var[:,1]]
	                continue
	            end
	            if !in(l[1], gene_ids) continue end
	            var[Symbol("$(l[1])_$t_name")] = parse.(Float64,l[indices])
				#plot date vs pred. expression
				s_plot = Plots.scatter(var[:,2],var[Symbol("$(l[1])_$t_name")],xlabel="Age (yBP)",ylabel="Pred. Expr.",title ="$(l[1])_$t_name",margin=10Plots.mm)
	            Plots.savefig(s_plot,"$(out_path)$(l[1])_$(t_name)_time_series.pdf")
				if length(cov_cols) > 0
					#fit regression with covariates as well as date
					println("covariate regression not implemented yet")
					exit()
				else
					fm = @eval @formula($(Symbol("$(l[1])_$t_name")) ~ $(date_name))
					model = lm(fm,var)#fit regression
					println(model)
				end
	        end
	    end
	end
end

function main()
    parsed_args = parseCommandline()
    timeSeries(parsed_args["regulation"],parsed_args["genes"],parsed_args["column"],parsed_args["samples"],parsed_args["date_col"],parsed_args["covariates"],parsed_args["out_dir"])
end

main()

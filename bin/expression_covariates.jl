# expression_covariates.jl
# @author Laura Colbran
#
# converts zhang expression text file into correct format, and calculates PEER covariates.
# julia 1.4
# dependencies:

using ArgParse
using GZip

NUM_SAMPLES = 10
RSEM_THRESH = 0.5
READ_THRESH = 6

# parses command-line arguments
function parseCommandline()
	s = ArgParseSettings()
	# @add_arg_table! s begin
	@add_arg_table! s begin
		"--expression", "-e"
			help="path(s) to text file with zhang expression file"
			arg_type = String
		"--pca", "-c"
			help="path(s) to pca eigenvectors"
			arg_type = String
		"--out_name", "-o"
            help = "path and prefix to write output files to"
            arg_type = String
            default = "./out"
	end
	return parse_args(s)
end

# reformats zhang expression file to tsv matrix for use in calculating PEER factors and normalization. goes ahead and filters by rsem and read counts
function expMatrix(expr_path::String,out_path::String)
    fpkm = Dict{String,Array{Float64,1}}() #gene => expr
    rsem = Dict{String,Array{Float64,1}}() #gene => expr
    counts = Dict{String,Array{Float64,1}}() #gene => expr
    ind_ids = Set{String}()
    GZip.open(expr_path) do f
        for line in eachline(f)
            if startswith(line,"Sample") continue end
            l = split(chomp(line),"\t")
            push!(ind_ids,l[1])
            if haskey(fpkm,l[5])
                push!(fpkm[l[5]],parse(Float64,l[8]))
                push!(rsem[l[5]],parse(Float64,l[9]))
                push!(counts[l[5]],parse(Float64,l[6]))
            else
                fpkm[l[5]] = [parse(Float64,l[8])]
                rsem[l[5]] = [parse(Float64,l[9])]
                counts[l[5]] = [parse(Float64,l[6])]
            end
        end
    end
    open("$out_path.txt","w") do f
        write(f,"#gene_id\t$(join(ind_ids,'\t'))\n")
        for gene in keys(fpkm)
            if length(findall(i -> i > RSEM_THRESH,rsem[gene])) < NUM_SAMPLES continue end
            if length(findall(i -> i >= READ_THRESH,counts[gene])) < NUM_SAMPLES continue end
            write(f,"$(gene)\t$(join(fpkm[gene],'\t'))\n")
        end
    end
    run(`gzip "$out_path.txt"`)
end

# wrapper function to calculate peer factors
function peerFactors()
end

#combines PCA and PEER factor output into a covariate class for predixcan processing
function combineFiles()
end

function main()
    parsed_args = parseCommandline()
    expMatrix(parsed_args["expression"],parsed_args["out_name"])
    peerFactors(parsed_args["out_name"])
    combineFiles()
end

main()

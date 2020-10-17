# expression_covariates.jl
# @author Laura Colbran
#
# converts zhang expression text file into correct format, and calculates PEER covariates.
# julia 1.4
# dependencies: R >3.3 with PEER installed

using ArgParse
using GZip
using RCall
using DataFrames, CSV

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
            help = "path and prefix to write output files to. (or prefix for matrix to normalize)"
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
            if haskey(fpkm,l[4])
                push!(fpkm[l[4]],parse(Float64,l[8]))
                push!(rsem[l[4]],parse(Float64,l[9]))
                push!(counts[l[4]],parse(Float64,l[6]))
            else
                fpkm[l[4]] = [parse(Float64,l[8])]
                rsem[l[4]] = [parse(Float64,l[9])]
                counts[l[4]] = [parse(Float64,l[6])]
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
    run(`gzip -f "$out_path.txt"`)
end

# normalizes expression in prep for PEER values and PrediXcan
function normExpr(matrix_path::String)
	expr =  CSV.read(GZip.open("$(matrix_path).txt.gz");delim='\t')
	gene_ids = expr[Symbol("#gene_id")]
	deletecols!(expr,Symbol("#gene_id"))
	@rput expr
	@rput gene_ids
	@rput matrix_path
	R"""
		# borrowed code from here: https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/
		# tweaked tie handling to "average" to match bioconductor
		quantileNormalisation <- function(df){
		    df_rank <- apply(df,2,rank,ties.method="average")
		    df_sorted <- data.frame(apply(df, 2, sort))
		    df_mean <- apply(df_sorted, 1, mean)

		    index_to_mean <- function(my_index, my_mean){
		    return(my_mean[my_index])
		    }

		    df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
		    rownames(df_final) <- rownames(df)
		    return(df_final)
		}
		# normalize across samples
		expr <- quantileNormalisation(expr)

	    # "For each gene, expression values were inverse quantile normalized to a standard normal distribution across samples."
	    # command from eric's link (pg18): https://media.nature.com/original/nature-assets/nature/journal/v490/n7419/extref/nature11401-s1.pdf
	    expr <- t(expr)
		# normalize each gene
	    for (i in 1:ncol(expr)) {
	        norm_dist <- qnorm((rank(expr[,i],na.last="keep")-0.5)/sum(!is.na(expr[,i])))
	        expr[,i] <- norm_dist
	    }
		expr = t(expr)
		num_col = ncol(expr)
		expr = cbind(expr,gene_ids)
		expr = expr[,c(ncol(expr),1:num_col)]
		write.table(expr,paste(matrix_path,"_norm.txt",sep=""),quote=FALSE,row.names = FALSE,sep='\t')
	"""
	run(`gzip -f $(matrix_path)_norm.txt`)
end

# wrapper function to calculate peer factors and get residuals
function peerFactors(matrix_path::String,pca_path::String)
	expr =  CSV.read(GZip.open("$(matrix_path)_norm.txt.gz");delim='\t')
	gene_ids = expr[:gene_ids]
	deletecols!(expr,:gene_ids)
	ind_ids = names(expr)

	covs = CSV.read(pca_path;delim='\t')
	deletecols!(covs,:Column1)
	col_names = ["id","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","pop"]
	rename!(covs,Symbol.(col_names))
	for i in 1:nrow(covs)
		covs[i,:id] = split(covs[i,:id],":")[2]
	end
	# sort covs to match expr
	covs = covs[[:pc1,:pc2,:pc3]] #match Zhang predixcan run

	@rput expr
	@rput covs
	@rput gene_ids
	@rput matrix_path
	@rput ind_ids
	R"""
	suppressPackageStartupMessages(library(peer))
	model = PEER()
	PEER_setPhenoMean(model,as.matrix(t(expr)))
	PEER_setNk(model,15) #match Zhang eQTL
	PEER_update(model)
	PEER_setCovariates(model, as.matrix(covs))
	residuals = PEER_getResiduals(model)

	rownames(residuals) <- ind_ids
	residuals = t(residuals)
	rownames(residuals) <- gene_ids
	write.table(residuals,paste(matrix_path,"_residuals.txt",sep=""),quote=FALSE,row.names = TRUE,sep='\t')
	"""
	run(`gzip -f $(matrix_path)_residuals.txt`)
end

function main()
    parsed_args = parseCommandline()
    # expMatrix(parsed_args["expression"],parsed_args["out_name"])
    # normExpr(parsed_args["out_name"])
    peerFactors(parsed_args["out_name"],parsed_args["pca"])
end

main()

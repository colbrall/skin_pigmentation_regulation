# predixcan_search.jl
# given ensembl id(s), searches predixcan model database(s) for a model. returns SNP count and r2
# @author Laura Colbran
# julia 1.4

using SQLite
using ArgParse
using DataFrames
using DBInterface

# parses command-line arguments
function parseCommandLine()
	s = ArgParseSettings()
	@add_arg_table! s begin
		"--databases", "-d"
			nargs='*'
			help="path to model database(s)"
			arg_type = String
			required = true
		"--genes","-g"
			help="file with gene IDs to search for"
			arg_type=String
		"--column","-c"
			help = "1-indexed column containing gene IDs to search for"
			arg_type=Int64
	end
	return parse_args(s)
end

function dbSearch(dbs::Array{String,1},gene_file::String,col::Int64)
	println("gene_id\tdb\tnum_SNPs\tr2")
	open(gene_file) do f
		for line in eachline(f)
			if startswith(line, "#") continue end
			l = split(chomp(line),"\t")
			q = "SELECT * FROM 'extra' WHERE gene='$(l[col])'"
			for db in dbs
				db_name = join(split(basename(db),".")[1:(end-1)],".")
				model = DataFrame(DBInterface.execute(SQLite.DB(db),q))
				for r in 1:nrow(model)
					println("$(l[col])\t$(db_name)\t$(model[r,Symbol("n.snps.in.model")])\t$(model[r,Symbol("pred.perf.R2")])")
				end
			end
		end
	end
end

function main()
	parsed_args = parseCommandLine()
	dbSearch(parsed_args["databases"],parsed_args["genes"],parsed_args["column"])
end

main()

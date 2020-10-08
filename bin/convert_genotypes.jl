# convert_genotypes.jl
# @author Laura Colbran
#
# convert illumina A/B genotype text file to plink with TOP genotypes. then wraps strand conversion.
# works with Zhang data. not necessarily transferable
# julia 1.4
# dependencies: transposeFile.awk

using ArgParse
using GZip

TRANS_PATH = "/project/mathilab/colbranl/bin/transposeFile.awk"
ALLELES = Dict("A" => 1, "C" => 2, "G" => 3, "T" => 4,"a" => 1, "c" => 2, "g" => 3, "t" => 4)
# parses command-line arguments
function parseCommandline()
	s = ArgParseSettings()
	# @add_arg_table! s begin
	@add_arg_table! s begin
		"--array_map", "-a"
			help="path(s) to text file with A/B to TOP allele matching by rsID"
			arg_type = String
			required = true
		"--genotypes","-g"
			help = "1-indexed column containing gene IDs to search for"
			arg_type=String
			required = true
		"--out_name", "-o"
            help = "path and prefix to write output files to"
            arg_type = String
            default = "./out"
	end
	return parse_args(s)
end

# read A/B -> TOP file into dictionary. returns {rsID -> (A,B)}
function topDict(map_file::String)
    map = Dict{String,Array{String,1}}()
    open(map_file) do f
        for line in eachline(f)
            l = split(chomp(line),"\t")
            map[l[1]] = split(l[3]," ")
        end
    end
    return map
end

# given an array of indices, pulls IDs out of another array
function pullIDs(l::Array{SubString{String},1},inds::Array{Int64,1})
    ids = String[]
    for index in inds
        push!(ids,split(l[index],".")[1])
    end
    return ids
end

# converts genotypes of a snp
function convertGenos(genos::Array{SubString{String}, 1},alleles::Array{String,1})
    new_genos = String[]
    for g in genos
        if g == "NC"
            push!(new_genos,"0 0")
        elseif g == "AA"
            push!(new_genos,"$(ALLELES[alleles[1]]) $(ALLELES[alleles[1]])")
        elseif g == "AB"
            push!(new_genos,"$(ALLELES[alleles[1]]) $(ALLELES[alleles[2]])")
        elseif g == "BB"
            push!(new_genos,"$(ALLELES[alleles[2]]) $(ALLELES[alleles[2]])")
        else
            println("ERROR-- weird genotype present: $g")
            exit()
        end
    end
    return new_genos
end

# converts A/B alleles to TOP alleles, and writes to plink PED and MAP files
function abToTOP(map_file::String,geno_file::String,out_path::String)
    allele_map = topDict(map_file)
    gtype_indices = Int64[]
    GZip.open(geno_file) do genf
        open("tmp.txt","w") do tped
            open("tmp.map","w") do map
                for line in eachline(genf)
                    l = split(chomp(line),"\t")
                    if startswith(line,"Index")
                        gtype_indices = findall(i->endswith(i,"GType"),l)
                        write(tped,"$(join(repeat([0],length(gtype_indices)),' '))\n")
                        write(tped,"$(join(pullIDs(l,gtype_indices),' '))\n")
                        write(tped,"$(join(repeat([0],length(gtype_indices)),' '))\n")
                        write(tped,"$(join(repeat([0],length(gtype_indices)),' '))\n")
                        write(tped,"$(join(repeat([0],length(gtype_indices)),' '))\n")
                        write(tped,"$(join(repeat(["-9"],length(gtype_indices)),' '))\n")
                        continue
                    end
                    rsid = l[2]
                    chr = l[4]
                    pos = l[5]
                    write(tped,"$(join(convertGenos(l[gtype_indices],allele_map[rsid]),' '))\n")
                    write(map,"$chr\t$rsid\t0\t$pos\n")
                end
            end
        end
    end
    run(`awk -v delim=" " -f $TRANS_PATH tmp.txt '>' $(out_path).ped`)
    run(`rm tmp.txt`)
end

function main()
    parsed_args = parseCommandline()
    abToTOP(parsed_args["array_map"],parsed_args["genotypes"],parsed_args["out_name"])
end

main()

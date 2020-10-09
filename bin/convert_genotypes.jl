# convert_genotypes.jl
# @author Laura Colbran
#
# convert illumina A/B genotype text file to plink with TOP genotypes. then wraps strand conversion script, and outputs some QC things
# works with Zhang data. not necessarily transferable
# julia 1.4
# dependencies: transposeFile.awk, plink, update_build.sh

using ArgParse
using GZip
using MultivariateStats
using Plots,CSV,DataFrames

TRANS_PATH = "/project/mathilab/colbranl/bin/transposeFile.awk"
STRAND_PATH = "/project/mathilab/colbranl/pigmentation/bin/update_build.sh"
ALLELES = Dict("A" => 1, "C" => 2, "G" => 3, "T" => 4,"a" => 1, "c" => 2, "g" => 3, "t" => 4) #not currently, but if i end up needing number versions instead

# parses command-line arguments
function parseCommandline()
	s = ArgParseSettings()
	# @add_arg_table! s begin
	@add_arg_table! s begin
		"--array_map", "-a"
			help="path(s) to text file with A/B to TOP allele matching by rsID"
			arg_type = String
		"--genotypes","-g"
			help = "1-indexed column containing gene IDs to search for"
			arg_type=String
		"--strand_ref","-s"
			help = "path to strand file for update_build.sh to use"
			arg_type=String
		"--out_name", "-o"
            help = "path and prefix to write output files to"
            arg_type = String
            default = "./out"
		"--pop_path","-p"
			help="path to file with sample->pop mapping for the PCA."
			arg_type=String
		"--eigen_vecs","-e"
			help="output from simplepca to plot"
			arg_type=String
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
    new_genos_1 = String[]
	new_genos_2 = String[]
    for g in genos
        if g == "NC"
            push!(new_genos_1,"0")
			push!(new_genos_2,"0")
        elseif g == "AA"
            push!(new_genos_1,"$(alleles[1])")
			push!(new_genos_2,"$(alleles[1])")
        elseif g == "AB"
            push!(new_genos_1,"$(alleles[1])")
			push!(new_genos_2,"$(alleles[2])")
        elseif g == "BB"
            push!(new_genos_1,"$(alleles[2])")
			push!(new_genos_2,"$(alleles[2])")
        else
            println("ERROR-- weird genotype present: $g")
            exit()
        end
    end
    return new_genos_1,new_genos_2
end

# converts A/B alleles to TOP alleles, and writes to tmp plink PED and MAP files
function abToTOP(map_file::String,geno_file::String)
    allele_map = topDict(map_file)
    gtype_indices = Int64[]
    GZip.open(geno_file) do genf
        open("tmp.txt","w") do tped
            open("tmp.map","w") do map
                for line in eachline(genf)
                    l = split(chomp(line),"\t")
                    if startswith(line,"Index")
                        gtype_indices = findall(i->endswith(i,"GType"),l)
						println(gtype_indices)
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
					if !haskey(allele_map,rsid)
						println("Skipped: $rsid\t$chr\t$pos")
						continue
					end
					new_genos_1,new_genos_2 = convertGenos(l[gtype_indices],allele_map[rsid])
                    write(tped,"$(join(new_genos_1,' '))\n")
                    write(tped,"$(join(new_genos_2,' '))\n")
                    write(map,"$chr\t$rsid\t0\t$pos\n")
                end
            end
        end
    end
	# make ped, then convert to bed
    run(pipeline(`awk -v delim=" " -f $TRANS_PATH tmp.txt`,stdout="tmp.ped"))
    run(`rm tmp.txt`)
end

# wrapper for update_build.sh (converts to binary files first)
function strandFlip(strand_ref::String,out_path::String)
	run(`plink --file tmp --make-bed --out tmp`)
	# use update_build to flip strands to match hg19
	run(`$(STRAND_PATH) tmp $(strand_ref) $(out_path)`)
	run(`rm tmp."*"`)
end

# makes PCA plots (and other QC things?) using only autosomes
function qc(eigen_path::String)
	eigen_vecs = CSV.read(eigen_path;delim='\t')
	deletecols!(eigen_vecs,:Column1)
	col_names = ["id","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","pop"]
	rename!(eigen_vecs,Symbol.(col_names))
	p = Plots.scatter(eigen_vecs[:pc1],eigen_vecs[:pc2],xlabel="PC1",ylabel="PC2",group=eigen_vecs[:pop],margin=10Plots.mm,legend=:best)
	Plots.savefig(p,"pc1_pc2.pdf")
end

function main()
    parsed_args = parseCommandline()
    # abToTOP(parsed_args["array_map"],parsed_args["genotypes"])
	# strandFlip(parsed_args["strand_ref"],parsed_args["out_name"])
	qc(parsed_args["eigen_vecs"])
end

main()

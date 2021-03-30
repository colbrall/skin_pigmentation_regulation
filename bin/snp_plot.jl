-# snp_plot.jl
# @author Laura Colbran
# given a list of snps and populations, makes a plot of allele dosages in individuals
# assumes you're giving SNPs on a specific chromosome
# also just assumes you want the alternate allele dosage.
# Julia 1.5. reqs R 4.0 with pegas installed to call hap_network.R

using ArgParse
using Seaborn
using GZip,DataFrames,CSV

Seaborn.set(style="white", palette="muted")
set_style(Dict("font.family" =>["DejaVu Sans"]))

SCRIPT_DIR = dirname(Base.source_path())

# columns to parse dosage file
RSID_COL = 2
REF_COL = 4
ALT_COL = 5
FIRST_DOS_IND = 7
FIRST_VCF_IND = 10
POS_COL = 2

function parseCommandLine()
    s = ArgParseSettings()

        @add_arg_table s begin
            "--pop_file","-f"
                help = "path to genotype file. required. expects dosage for heatmap and vcf for network."
                arg_type = String
                required = true
            "--snps","-s"
                help = "file with location and ids of SNPs to plot"
                arg_type = String
                required = true
            "--target","-t"
                help = "id for target SNP to sort by"
                arg_type = String
                default = "none"
            "--pops","-p"
                help = "file with individual -> population assignments"
                arg_type = String
            "--delim","-d"
                help = "delimiter for pop_file"
                arg_type = String
                default = "\t"
            "--heatmap","-m"
                help = "to plot a heatmap of dosages in individuals"
                action=:store_true
            "--network","-n"
                help = "to build a haplotype network"
                action = :store_true
            "--window"
                help = "if pulling all SNPs in window. if not grab specific snps in --snps"
                action = :store_true
            "--targ_pop", "-r"
                help = "list of individuals to include in hap plot"
                default = ""
                arg_type = String
            "--ind_col", "-i"
                help = "column in targ_pop file with the ids"
                default = 1
                arg_type = Int64
            "--min_freq","-q"
                help = "minimum MAF of SNP to include"
                default = 0.0
                arg_type = Float64
    end
    return parse_args(s)
end

function readSNPs(file::String)
    # read in target snps by coordinate
    snps = Array{Any,1}[] #[[loc,rsid,ref,alt]]
    open(file) do f
        for line in eachline(f)
            if startswith(line,"#") continue end
            l = split(chomp(line),"\t")
            snps = vcat(snps,[[parse(Int64,l[2]),l[6],l[4],l[5]]])
        end
    end
    # sort snps by location (bc it's first entry in arrays)
    sort!(snps)
    return snps
end

#return array of individuals to filter
function readInds(path::String,col::Int64)
    inds = String[]
    open(path) do f
        for line in eachline(f)
            inds = vcat(inds,[split(chomp(line),"\t")[col]])
        end
    end
    return inds
end

function filterAF(loci::Array{AbstractString,2},min_freq::Float64)
    snps_to_keep = Int64[1]
    for i in 2:size(loci)[2]
        freq = length(findall(x->isequal("0",x),loci[1:end,i]))/(size(loci)[1]-1)
        if (1 - freq) >= min_freq
            snps_to_keep = vcat(snps_to_keep,[i])
        end
    end
    return loci[1:end,snps_to_keep]
end

# grabs all variants in the range of target snp set, and writes them to a temporary file for pegas to read
function writeLociWindow(snps::Array{Array{Any,1},1},fpath::String,targpop::String,col::Int64,min_freq::Float64)
    min,max = snps[1][1],snps[end][1]
    loci = Array{String,2}
    indices = Int64[]
    GZip.open(fpath) do inf
        for line in eachline(inf)
            if startswith(line,"##") continue end
            if startswith(line,"#")
                l = split(chomp(line),"\t")
                # filter individuals if necessary
                if targpop == ""
                    indices = collect(FIRST_VCF_IND:length(l))
                else
                    targ = readInds(targpop,col)
                    indices = [findfirst(x->x==i,l) for i in targ]
                end
                loci = vcat("inds",l[indices])
                continue
            end
            line = split(chomp(line),"\t")
            pos = parse(Int64,line[POS_COL])
            if pos < min continue end
            if pos > max break end #assumes vcf is sorted
            loci =hcat(loci,vcat(line[POS_COL],vcat([split(split(i,":")[1],"/")[2] for i in line[indices]]))) # add genotypes
        end
    end
    println(size(loci))
    if min_freq > 0
        loci = filterAF(loci,min_freq)
    end
    println(size(loci))
    CSV.write("loci.txt",  DataFrame(loci); writeheader=false,delim = "\t")
    println("Loci written...")
    return (size(loci)[2] -1)
end

# grabs all specific variants in target snp set, and writes them to a temporary file for pegas to read
function writeLociSpec(snps::Array{Array{Any,1},1},fpath::String,targpop::String,col::Int64)
    snp_pos = [entry[1] for entry in snps]
    loci = Array{String,2}
    indices = Array{Int64,1}
    GZip.open(fpath) do inf
        for line in eachline(inf)
            if startswith(line,"##") continue end
            if startswith(line,"#")
                l = split(chomp(line),"\t")
                # filter individuals if necessary
                if targpop == ""
                    indices = collect(FIRST_VCF_IND:length(l))
                else
                    targ = readInds(targpop,col)
                    indices = [findfirst(x->x==i,l) for i in targ]
                end
                loci = vcat("inds",l[indices])
                continue
            end
            line = split(chomp(line),"\t")
            pos = parse(Int64,line[POS_COL])
            if pos in snp_pos
                println(pos)
                loci =hcat(loci,vcat(line[POS_COL],vcat([split(split(i,":")[1],"/")[2] for i in line[indices]]))) # add genotypes
            end
        end
    end
    println(size(loci))
    println(loci[1:6])
    CSV.write("loci.txt",  DataFrame(loci); writeheader=false,delim = "\t")
    println("Loci written...")
    return (size(loci)[2] -1)
end

# plots heatmap of dosages of SNPs of interest
function plotSNPs(snp_file::String,target_id::String,dos_path::String,pop_path,delim::String)
    snps = readSNPs(snp_file)
    # build array of snps by dosage in individuals
    inds = split(chomp(read(pipeline(`zcat $dos_path`,`head -1`),String)),delim)[FIRST_DOS_IND:end] #currently doesn't make sense if the first line isn't a header, but also doesn't break
    println("Num individuals: $(length(inds))")
    dosages = fill(0.0,length(inds),length(snps)) #allocate array of correct size
    targ_ind = 1 #sorts by first if no target specified
    for i in 1:length(snps)
        if snps[i][2] == target_id
            targ_ind = i
        end
        cmd = `zgrep "$(snps[i][2])" $dos_path`
        try
            line = split(chomp(read(cmd,String)),delim)
            snps[i][3] = line[REF_COL]
            snps[i][4] = line[ALT_COL]
            dosages[:,i] .= [parse(Float64,x) for x in line[FIRST_DOS_IND:end]]
        catch y
            println("$(snps[i][2]) missing from population. dosage set to 0 for all.")
       end
    end
    # sort individuals by specific SNP dosage
    dosages = sortslices(dosages,dims = 1,by=x->x[targ_ind],rev=false)
    # plot heatmap dosage
    h_plot = heatmap(dosages, yticklabels=false, xticklabels=["$(x[2]) $(x[3])/$(x[4])" for x in snps],cmap="Purples",cbar=false)
    h_plot.set_xticklabels(h_plot.get_xticklabels(), rotation=45)
    Seaborn.savefig("$(target_id)_heatmap.pdf")
    clf()
end

# plots haplotype network for SNPs of interest
function hapNet(snp_file::String,vcf_path::String,target::String,window::Bool,targpop::String,ind_col::Int64,min_freq::Float64)
    snps = readSNPs(snp_file)
    if window
        nloci = writeLociWindow(snps,vcf_path,targpop,ind_col,min_freq)
    else
        nloci = writeLociSpec(snps,vcf_path,targpop,ind_col)
    end
    println("N Loci: $nloci")
    # run(`Rscript $(SCRIPT_DIR)/hap_network.R loci.txt $nloci $target`)
end

function main()
    parsed_args = parseCommandLine()
    if parsed_args["heatmap"]
        plotSNPs(parsed_args["snps"],parsed_args["target"],parsed_args["pop_file"],parsed_args["pops"],parsed_args["delim"])
    end
    if parsed_args["network"]
        hapNet(parsed_args["snps"],parsed_args["pop_file"],parsed_args["target"],parsed_args["window"],parsed_args["targ_pop"],parsed_args["ind_col"],parsed_args["min_freq"])
    end
end

main()

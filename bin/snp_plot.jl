# snp_plot.jl
# @author Laura Colbran
# given a list of snps and populations, makes a plot of allele dosages in individuals
# assumes you're giving SNPs on a specific chromosome
# also just assumes you want the alternate allele dosage.
# Julia 1.5. reqs R 4.0 with pegas installed

using ArgParse
using Seaborn
using RCall
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
    end
    return parse_args(s)
end

function readSNPs(file::String)
    # read in target snps by coordinate
    snps = Array{Any,1}[] #[[loc,rsid,ref,alt]]
    chr = ""
    open(file) do f
        for line in eachline(f)
            if startswith(line,"#") continue end
            l = split(chomp(line),"\t")
            chr = l[1]
            snps = vcat(snps,[[parse(Int64,l[2]),l[6],l[4],l[5]]])
        end
    end
    # sort snps by location (bc it's first entry in arrays)
    sort!(snps)
    return chr,snps
end

# grabs all variants in the range of target snp set, and writes them to a temporary file for pegas to read
function writeLoci(min::Int64,max::Int64,fpath::String)
    loci = Array{String,2}
    N = 0
    GZip.open(fpath) do inf
        for line in readlines(inf)
            if startswith(line,"##") continue end
            if startswith(line,"#")
                loci = vcat("inds",split(chomp(line),"\t")[FIRST_VCF_IND:end])
                continue
            end
            line = split(chomp(line),"\t")
            pos = parse(Int64,line[POS_COL])
            if pos < min continue end
            if pos > max break end #assumes vcf is sorted
            loci =hcat(loci,vcat(line[POS_COL],vcat([split(i,":")[1] for i in line[FIRST_VCF_IND:end]]))) # add genotypes
            N += 1
        end
    end
    # println(loci)
    CSV.write("loci.txt",  DataFrame(loci); writeheader=false,delim = "\t")
    return N
end

# plots heatmap of dosages of SNPs of interest
function plotSNPs(snp_file::String,target_id::String,dos_path::String,pop_path,delim::String)
    chr,snps = readSNPs(snp_file)
    # build array of snps by dosage in individuals
    inds = split(chomp(read(pipeline(`zcat $dos_path`,`head -1`),String)),delim)[FIRST_DOS_IND:end] #currently doesn't make sense if the first line isn't a header, but also doens't break
    println("Num individuals: $(length(inds))")
    dosages = fill(0.0,length(inds),length(snps)) #allocate array of correct size
    targ_ind = 0
    for i in 1:length(snps)
        if snps[i][2] == target_id
            targ_ind = i
        end
        cmd = `zgrep "$(snps[i][2])" $dos_path`
        try
            line = split(chomp(read(cmd,String)),delim)
            # println(line)
            snps[i][3] = line[REF_COL]
            snps[i][4] = line[ALT_COL]
            dosages[:,i] .= [parse(Float64,x) for x in line[FIRST_DOS_IND:end]]
            # exit()
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
function hapNet(snp_file::String,vcf_path::String)
    chr,snps = readSNPs(snp_file)
    # println(chr,snps)
    nloci = writeLoci(snps[1][1],snps[end][1],vcf_path)
    # println(nloci)
    run(`Rscript $(SCRIPT_DIR)/hap_network.R $nloci`)
end

function main()
    parsed_args = parseCommandLine()
    if parsed_args["heatmap"]
        plotSNPs(parsed_args["snps"],parsed_args["target"],parsed_args["pop_file"],parsed_args["pops"],parsed_args["delim"])
    end
    if parsed_args["network"]
        hapNet(parsed_args["snps"],parsed_args["pop_file"])
    end
end

main()

# snp_plot.jl
# @author Laura Colbran
# given a list of snps and populations, makes a plot of allele dosages in individuals
# assumes you're giving SNPs on a specific chromosome
# also just assumes you want the alternate allele dosage.

using ArgParse
using Seaborn

Seaborn.set(style="white", palette="muted")
set_style(Dict("font.family" =>["DejaVu Sans"]))

# columns to parse dosage file
RSID_COL = 2
REF_COL = 4
ALT_COL = 5
FIRST_IND = 7

function parseCommandLine()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--pop_file","-f"
            help = "path to dosage file containing populations. required."
            arg_type = String
        "--snps","-s"
            help = "file with location and ids of SNPs to plot"
            arg_type = String
        "--target","-t"
            help = "id for target SNP to sort by"
            arg_type = String
        "--pops","-p"
            help = "file with individual -> population assignments"
            arg_type = String
        "--delim","-d"
            help = "delimiter for dosage file"
            arg_type = String
            default = "\t"
    end
    return parse_args(s)
end

function plotSNPs(snp_file::String,target_id::String,dos_path::String,pop_path,delim::String)
    # read in target snps by coordinate
    snps = Array{Any,1}[] #[[loc,rsid,ref,alt]]
    chr = ""
    open(snp_file) do f
        for line in eachline(f)
            if startswith(line,"#") continue end
            l = split(chomp(line),"\t")
            chr = l[1]
            snps = vcat(snps,[[parse(Int64,l[2]),l[6],l[4],l[5]]])
        end
    end
    # sort snps by location (bc it's first entry in arrays)
    sort!(snps)
    # build array of snps by dosage in individuals
    inds = split(chomp(read(pipeline(`zcat $dos_path`,`head -1`),String)),delim)[FIRST_IND:end] #currently doesn't make sense if the first line isn't a header, but also doens't break
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
            dosages[:,i] .= [parse(Float64,x) for x in line[FIRST_IND:end]]
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

function main()
    parsed_args = parseCommandLine()
    plotSNPs(parsed_args["snps"],parsed_args["target"],parsed_args["pop_file"],parsed_args["pops"],parsed_args["delim"])
end

main()

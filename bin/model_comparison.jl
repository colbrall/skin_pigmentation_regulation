# @author Laura Colbran
# 2019-02-20
# compares how well two sets of PrediXcan models match
#
# julia 1.4

using ArgParse
using DataFrames
using StatsBase
using HypothesisTests
using StatsPlots
using SQLite
using GZip
using CSV

default(color=:blue,leg=false,grid=false,fontfamily="arial",alpha=0.5)

# parses command-line arguments
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--databases","-d"
            nargs='*'
            help = "paths to the model databases"
            arg_type = String
            required = true
        "--performance"
            help = "to run comparison between certain models' performances and SNP counts across different training instances"
            action = :store_true
        "--reference","-r"
            nargs=1
            help = "paths to the reference database"
            arg_type = String
        "--targets","-t"
            help = "models to include, if you want to not include all"
            arg_type = String
        "--box"
            help = "to make plots of performance metrics"
            action = :store_true
        # "--snps"
        #     help = "to compare specific SNPs and their weights"
        #     action = :store_true
        "--genomes","-g"
            help = "path to directory containing genomes"
            arg_type = String
        "--missingness"
            help = "to do analyses on missing SNPs by model"
            action = :store_true
        "--summarize"
            help="summarize genes across tissues"
            action=:store_true
    end
    return parse_args(s)
end

function readGenes(gene_file::String,tiss::SubString{String}) #return array of gene IDs to keep in tissue
    genes = CSV.read(gene_file;delim="\t")
    return genes[genes[:tissue] .== tiss,:gene]
end

# for each SNP, count the number of people missing it, and for each individual,
function missingDict(path::String)
    ind_dict = Dict{String,Array{String,1}}()
    snp_dict = Dict{String,Array{String,1}}()
    for item in readdir(path)
        if !endswith(item,"dos.gz") continue end
        println(item)
        GZip.open("$(realpath(path))/$item") do f
            for line in eachline(f)
                l = split(chomp(line),"\t")
                ind_ids = [string(i) for i in collect(1:length(l)-6)] #if ids aren't present, ids will just be index of individual from 1
                if startswith(line,"#")
                    ind_ids = l[7:length(l)]
                    continue
                end
                for i in collect(7:length(l))
                    if l[i] == "NA"
                        #println(line)
                        if in(l[2], keys(snp_dict))
                            append!(snp_dict[l[2]],[ind_ids[i-6]])
                        else
                            snp_dict[l[2]] = [ind_ids[i-6]]
                        end
                        if in(ind_ids[i-6], keys(ind_dict))
                            append!(ind_dict[ind_ids[i-6]],[l[2]])
                        else
                            ind_dict[ind_ids[i-6]] = [l[2]]
                        end
                    end
                end
            end
        end
    end
    return ind_dict,snp_dict
end

#compares R2 and Number of SNPs in models from both sets
function performance(db_files::Array{String,1})
    q = "SELECT * FROM 'extra'"
    models_1 = DataFrame(DBInterface.execute(SQLite.DB(db_files[1]),q))
    for i in 1:nrow(models_1)
        models_1[i,:gene] = split(models_1[i,:gene],".")[1]
    end
    models_2 = DataFrame(DBInterface.execute(SQLite.DB(db_files[2]),q))
    for i in 1:nrow(models_2)
        models_2[i,:gene] = split(models_2[i,:gene],".")[1]
    end
    if :R2 in names(models_1)
        rename!(models_1, :R2 => Symbol("pred.perf.R2"))
        rename!(models_1, Symbol("n.snps") => Symbol("n.snps.in.model"))
    end
    if :R2 in names(models_2)
        rename!(models_2, :R2 => Symbol("pred.perf.R2"))
        rename!(models_2, Symbol("n.snps") => Symbol("n.snps.in.model"))
    end
    println("Number of models in $(db_files[1]): $(nrow(models_1))")
    println("R2:")
    describe(models_1[Symbol("pred.perf.R2")])
    println("NumSNPs:")
    describe(models_1[Symbol("n.snps.in.model")])
    println("Number of models in $(db_files[2]): $(nrow(models_2))")
    println("R2:")
    describe(models_2[Symbol("pred.perf.R2")])
    println("NumSNPs:")
    describe(models_2[Symbol("n.snps.in.model")])
    CSV.write("models_1_all.txt",models_1;delim='\t')
    CSV.write("models_2_all.txt",models_2;delim='\t')
    CSV.write("models_1_unique.txt",join(models_1,models_2,on=:gene,kind=:anti);delim='\t')
    CSV.write("models_2_unique.txt",join(models_2,models_1,on=:gene,kind=:anti);delim='\t')
    all_models = join(models_1,models_2,on=:gene,kind=:inner,makeunique=true)
    println("Number models shared: $(nrow(all_models))")

    #println(join(models_2,models_1,on=:gene,kind=:anti,makeunique=true))
    #println(first(all_models,6))
	println("In order of arguments, models that are shared:")
    println("R2:")
	describe(all_models[Symbol("pred.perf.R2")])
	describe(all_models[Symbol("pred.perf.R2_1")])
    describe(all_models[Symbol("n.snps.in.model")])
    describe(all_models[Symbol("n.snps.in.model_1")])

    rho = corspearman(convert(Array{Float64,1},all_models[Symbol("pred.perf.R2")]),convert(Array{Float64,1},all_models[Symbol("pred.perf.R2_1")]))
    #rho = corspearman(convert(Array{Float64,1},all_models[Symbol("pred.perf.R2")]),convert(Array{Float64,1},all_models[Symbol("R2")]))
    println("Spearman rho: $(rho)")
    println(OneSampleZTest(atanh(rho), 1, nrow(all_models)))
    println("Num SNPs:")

    rho = corspearman(convert(Array{Float64,1},all_models[Symbol("n.snps.in.model")]),convert(Array{Float64,1},all_models[Symbol("n.snps.in.model_1")]))
    #rho = corspearman(convert(Array{Float64,1},all_models[Symbol("n.snps.in.model")]),convert(Array{Float64,1},all_models[Symbol("n.snps")]))
    println("Spearman rho: $(rho)")
    println(OneSampleZTest(atanh(rho), 1, nrow(all_models)))

    #plots
    r2_1_plot = histogram(models_1[Symbol("pred.perf.R2")],xlabel="R2",ylabel="Num. Models",bins=50,margin=10Plots.mm)
    r2_1_plot = histogram!(all_models[Symbol("pred.perf.R2")],bins=50,color=:red)
    Plots.savefig(r2_1_plot,"arg1_r2.pdf")
    r2_2_plot = histogram(models_2[Symbol("pred.perf.R2")],xlabel="R2",ylabel="Num. Models",bins=50,margin=10Plots.mm)
    r2_2_plot = histogram!(all_models[Symbol("pred.perf.R2_1")],bins=50,color=:red)
    Plots.savefig(r2_2_plot,"arg2_r2.pdf")
    nsnps_1_plot = histogram(models_1[Symbol("n.snps.in.model")],xlabel="Num_SNPs",ylabel="Num. Models",bins=50,margin=10Plots.mm)
    nsnps_1_plot = histogram!(all_models[Symbol("n.snps.in.model")],bins=50,color=:red)
    Plots.savefig(nsnps_1_plot,"arg1_nsnps.pdf")
    nsnps_2_plot = histogram(models_2[Symbol("n.snps.in.model")],xlabel="Num_SNPs",ylabel="Num. Models",bins=50,margin=10Plots.mm)
    nsnps_2_plot = histogram!(all_models[Symbol("n.snps.in.model_1")],bins=50,color=:red)
    Plots.savefig(nsnps_2_plot,"arg2_nsnps.pdf")
    r2_corr_plot = scatter(all_models[Symbol("pred.perf.R2")],all_models[Symbol("pred.perf.R2_1")],xlabel="v6p r2",ylabel="v8 r2",margin=10Plots.mm)
    Plots.savefig(r2_corr_plot,"r2_corrplot.pdf")
end

# analyses of missing SNPs in models
function missingSNPs(db_files::Array{String,1},genome_path::String)
    ind_dict,snp_dict = missingDict(genome_path)
    for db in db_files
        name = splitext(splitpath(db)[end])[1]
        println(name)
        gene_snps = Dict{String,Array{String,1}}()
        q = "SELECT * FROM 'weights'"
        for snp in DBInterface.execute(SQLite.DB(db),q)
            #println(typeof(snp[:gene]))
            #println(typeof(snp[:rsid]))
            if !in(snp[:gene], keys(gene_snps))
                gene_snps[snp[:gene]] = [snp[:rsid]]
            else
                append!(gene_snps[snp[:gene]],[snp[:rsid]])
            end
        end
        gene_info = DataFrame(gene_id=String[],num_ind_any=Int64[],max_ind_per_snp=Int64[],num_snps_any=Int64[],max_snp_per_ind=Int64[])
        for gene in keys(gene_snps)
            ind_set = Set{String}()
            max = 0
            missing_snps = length(gene_snps[gene])
            for snp in gene_snps[gene]
                if !in(snp,keys(snp_dict)) # number of SNPs missing in at least one person for each model
                    missing_snps -= 1
                    continue
                end
                for ind in snp_dict[snp] # number of people missing a SNP in each model
                    push!(ind_set,ind)
                end
                if length(snp_dict[snp]) > max # max number of people missing a single SNP per model
                    max = length(snp_dict[snp])
                end
        # max number of missing SNPs in a single person for each model
            end
        end
    end
end

# plots boxplots of performance metrics
function boxplots(ref_file::Array{String,1},db_files::Array{String,1})
    q = "SELECT * FROM 'extra'"
    ref_models = DataFrame(DBInterface.execute(SQLite.DB(ref_file[1]),q))
    for i in 1:nrow(ref_models)
        ref_models[i,:gene] = split(ref_models[i,:gene],".")[1]
    end
    if :R2 in names(ref_models)
        rename!(ref_models, :R2 => Symbol("pred.perf.R2"))
        rename!(ref_models, Symbol("n.snps") => Symbol("n.snps.in.model"))
    end
    r2_boxes = DataFrames.DataFrame(R2 = ref_models[Symbol("pred.perf.R2")], set = repeat(["ref_$(nrow(ref_models))"],nrow(ref_models)),cat= repeat(["ref"],nrow(ref_models)))
    nsnps_boxes = DataFrames.DataFrame(nsnps = ref_models[Symbol("n.snps.in.model")], set = repeat(["ref_$(nrow(ref_models))"],nrow(ref_models)),cat= repeat(["ref"],nrow(ref_models)))
    i = 1
    for file in db_files
        println("$(i)_comp = $file")
        models = DataFrame(DBInterface.execute(SQLite.DB(file),q))
        for j in 1:nrow(models)
            models[j,:gene] = split(models[j,:gene],".")[1]
        end
        if :R2 in names(models)
            rename!(models, :R2 => Symbol("pred.perf.R2"))
            rename!(models, Symbol("n.snps") => Symbol("n.snps.in.model"))
        end
        append!(r2_boxes,DataFrames.DataFrame(R2 = models[Symbol("pred.perf.R2")], set = repeat(["$(i)_comp_$(nrow(models))"],nrow(models)),cat= repeat(["comp"],nrow(models))))
        append!(nsnps_boxes,DataFrames.DataFrame(nsnps = models[Symbol("n.snps.in.model")], set = repeat(["$(i)_comp_$(nrow(models))"],nrow(models)),cat= repeat(["comp"],nrow(models))))
        models = join(ref_models,models, on=:gene,kind=:inner,makeunique=true)
        append!(r2_boxes,DataFrames.DataFrame(R2 = models[Symbol("pred.perf.R2")], set = repeat(["$(i)_shared_$(nrow(models))"],nrow(models)),cat= repeat(["shared"],nrow(models))))
        append!(nsnps_boxes,DataFrames.DataFrame(nsnps = models[Symbol("n.snps.in.model")], set = repeat(["$(i)_shared_$(nrow(models))"],nrow(models)),cat= repeat(["shared"],nrow(models))))
        # scatter plot R2 ref vs. R2 file
        s_plot = scatter(models[Symbol("pred.perf.R2")],models[Symbol("pred.perf.R2_1")],xlabel="ref_models R2",ylabel="$(i)_comp R2",margin=10Plots.mm)
        Plots.savefig(s_plot,"scatter_$(i)_comp.pdf")
        i+=1
    end
    r2_plot = boxplot(r2_boxes[:set], r2_boxes[:R2],notch = true, ylabel="R2",xlabel= "Training Set",margin=10Plots.mm,outliers=false)
    Plots.savefig(r2_plot,"boxes_r2.pdf")
    nsnps_plot = boxplot(nsnps_boxes[:set], nsnps_boxes[:nsnps],notch = true, ylabel="N_SNPs",xlabel= "Training Set",margin=10Plots.mm,outliers=false)
    Plots.savefig(nsnps_plot,"boxes_nsnps.pdf")
    for lab in unique(r2_boxes[:set])
        println("$lab R2:")
        describe(r2_boxes[r2_boxes[:set].==lab,:R2])
    end
    for lab in unique(nsnps_boxes[:set])
        println("$lab num_snps:")
        describe(nsnps_boxes[nsnps_boxes[:set].==lab,:nsnps])
    end
end

#compares R2 and Number of SNPs in models from both sets
function modelSumm(db_files::Array{String,1},gene_file)
    q = "SELECT * FROM 'extra'"
    tiss_data = DataFrames.DataFrame(tissue=String[],num_models=Int64[])
    genes = Dict{String,Int64}()
    for file in db_files
        tiss = split(splitpath(file)[end],"_1240k_alpha0.5")[1]
        if occursin("_BA",tiss) tiss = split(tiss,"_BA")[1] end #to catch ACC  and frontal cortex
        models = DataFrame(DBInterface.execute(SQLite.DB(file),q))
        if typeof(gene_file) != Nothing
            targets = readGenes(gene_file,tiss)
            models = models[[in(i,targets) for i in models[:gene]],:]
        end
        push!(tiss_data,[tiss,nrow(models)])
        for gene in models[:gene]
            if in(gene,keys(genes))
                genes[gene] += 1
            else
                genes[gene] = 1
            end
        end
        if :R2 in names(models)
            rename!(models, :R2 => Symbol("pred.perf.R2"))
            rename!(models, Symbol("n.snps") => Symbol("n.snps.in.model"))
        end
        println("Number of models in $(file): $(nrow(models))")
        println("R2:")
        describe(models[Symbol("pred.perf.R2")])
        println("NumSNPs:")
        describe(models[Symbol("n.snps.in.model")])
        r2_1_plot = histogram(models[Symbol("pred.perf.R2")],xlabel="R2",ylabel="Num. Models",bins=50,margin=10Plots.mm)
        Plots.savefig(r2_1_plot,"$(tiss)_r2.pdf")
        nsnps_1_plot = histogram(models[Symbol("n.snps.in.model")],xlabel="Num_SNPs",ylabel="Num. Models",bins=50,margin=10Plots.mm)
        Plots.savefig(nsnps_1_plot,"$(tiss)_nsnps.pdf")
    end
    sort!(tiss_data,:num_models)
    tiss_plot = bar(tiss_data[:tissue],tiss_data[:num_models],xlabel="",ylabel="# Genes",margin=10Plots.mm, xtickfontcolor=:white)
    Plots.savefig(tiss_plot,"tissue_numModels.pdf")
    describe(tiss_data[:num_models])
    for t in tiss_data[:tissue]
        println(t)
    end

    gene_data = DataFrames.DataFrame(gene=String[],num_tiss=Int64[])
    for g in keys(genes)
        push!(gene_data,[g,genes[g]])
    end
    gene_plot = histogram(gene_data[:num_tiss],xlabel="# Tissues",ylabel="# Genes",margin=10Plots.mm,bins= 50)
    Plots.savefig(gene_plot,"gene_numTissues.pdf")
    describe(gene_data[:num_tiss])
end

function main()
    parsed_args = parse_commandline()
    if parsed_args["performance"]
        performance(parsed_args["databases"])
    end
    if parsed_args["missingness"]
        missingSNPs(parsed_args["databases"],parsed_args["genomes"])
    end
    if parsed_args["box"]
        boxplots(parsed_args["reference"],parsed_args["databases"])
    end
    if parsed_args["summarize"]
        modelSumm(parsed_args["databases"],parsed_args["targets"])
    end
end

main()

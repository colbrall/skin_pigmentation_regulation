# quick script to get sample summary by num_snps threshold.
# adapted from functions in snp_choice_analyses.jl

using CSV
using DataFrames
using Plots
using Seaborn
using StatsBase
using RCall

#set plotting defaults
default(color=:black,leg=false,grid=false,fontfamily="arial")

Seaborn.set(style="white", palette="muted")
set_style(Dict("font.family" =>["DejaVu Sans"]))

const N_SNPS_CHIP = 1233013

EUROPE = ["Austria","Belgium","Bulgaria","Croatia","Czech Republic","Denmark",
            "Estonia","Finland","France","","Germany","Great Britain",
            "Greece","Hungary","Iceland","Ireland","Italy","Kosovo","Latvia",
            "Liechtenstein","Lithuania","Luxembourg","Macedonia","Malta","Moldova",
            "Monaco","Montenegro","Netherlands","Norway","Poland","Portugal",
            "Romania","San Marino","Serbia","Slovakia","Slovenia","Spain","Sweden",
            "Switzerland","Ukraine","Vatican","Gibraltar"] #countries in europe in annotation file
ASIA = ["Armenia", "Cambodia", "China","Georgia","India","Iran","Israel","Japan",
        "Jordan", "Kazakhstan", "Kyrgyzstan", "Laos", "Lebanon","Mongolia",
        "Myanmar", "Nepal", "Russia","Thailand", "Turkey", "Turkmenistan","Vietnam",
        "Pakistan","Uzbekistan","Tajikistan","Afghanistan","Crimea"]
        #countries in asia in annotation file
N_AMERICA = ["Bahamas", "Belize", "Canada","Greenland","Mexico", "USA"] #countries in north america in annotation file
S_AMERICA = ["Argentina","Brazil","Chile", "Peru"] #countries in south america in annotation file
OCEANIA = ["French Polynesia", "Indonesia","Malaysia", "Philippines", "Solomon Islands",
            "Vanuatu"] #countries in oceania in annotation file
AFRICA = ["Canary Islands", "Egypt", "Ethiopia","Kenya","Malawi", "Morocco",
        "South Africa", "Tanzania","Tonga","Cameroon"] #countries in africa in annotation file
EXCLUDE = ["Ancestor.REF","Href.REF","Chimp.REF", "Gorilla.REF","Altai_snpAD.DG",
            "Denisova_snpAD.DG","Altai_published.DG","Denisova_published.DG",
            "Denisova11.SG","Vindija_snpAD.DG","VindijaG1_final.SG","Goyet_final_provisional.SG",
            "Les_Cottes_final_provisional.SG","Mezmaiskaya1_final_provisional.SG","VindijaG1_final_provisional.SG",
            "Spy_final_provisional.SG","Mezmaiskaya2_final_provisional.SG"] #non-human ids in annotation file

# parses Reich Lab sample annotation file
function parseAnnoFile(anno_file::String)
    inds = CSV.read(anno_file; delim='\t',allowmissing=:none)
    inds = inds[[2,3,4,7,9,13,18,19,33],]
#    inds = inds[[1,2,3,7,8,9],]
    names!(inds,[:ind_id,:master_id,:skeleton_id,:publication,:date,:country,:coverage,:num_snps,:assessment])
#    inds = inds[inds[:date] .!= "0",:]
    inds = inds[inds[:date] .!= 0,:] #remove present-day samples
    inds[:continent] = ".."
    inds[:flt_cov] = 0.0
    for i in 1:nrow(inds)
        # if inds[i,:date] == ".."
        #     inds[i,:ind_date] = -1
        # else
        #     inds[i,:ind_date] = parse(Int,inds[i,:date])
        # end

        if inds[i,:coverage] == ".."
            inds[i,:flt_cov] = -1.0
        else
            inds[i,:flt_cov] = parse(Float64,inds[i,:coverage])
        end

        if inds[i,:country] in EUROPE
            inds[i,:continent] = "Europe"
        elseif inds[i,:country] in ASIA
            inds[i,:continent] = "Asia"
        elseif inds[i,:country] in N_AMERICA
            inds[i,:continent] = "N_America"
        elseif inds[i,:country] in S_AMERICA
            inds[i,:continent] = "S_America"
        elseif inds[i,:country] in AFRICA
            inds[i,:continent] = "Africa"
        elseif inds[i,:country] in OCEANIA
            inds[i,:continent] = "Oceania"
        elseif inds[i,:country] == ".." #for 2 raghavan 2014 science samples.
            inds[i,:country] = "Canada"
            inds[i,:continent] = "N_America"
        end

        if inds[i,:ind_id] in EXCLUDE
            inds[i,:continent] = "NON-HUMAN"
        end

        if startswith(inds[i,:assessment],"PASS")
            inds[i,:assessment] = "yes"
        elseif startswith(inds[i,:assessment],"QUESTIONABLE_CRITICAL")
            inds[i,:assessment] = "no"
        elseif startswith(inds[i,:assessment],"QUESTIONABLE")
            inds[i,:assessment] = "maybe"
        end
    end
    # inds[:date] = inds[:ind_date]
    inds[:coverage] = inds[:flt_cov]
    # deletecols!(inds,[:ind_date,:flt_cov])
    deletecols!(inds,[:flt_cov])
    inds = inds[inds[:continent] .!== "NON-HUMAN",:] #remove archaic hominins, reference genomes
    return inds
end

function mapSampleIDs(df::DataFrames.DataFrame,map_file::String,info_file::String)
    map_ids = CSV.read(map_file)
    df = join(df,map_ids,on=:master_id => :reich_id,kind=:inner)
    df = join(df,CSV.read(info_file,allowmissing=:none),kind=:inner,on=:colbran_id => :sample,makeunique=true)
    return df
end

# builds a swarmplot with up to 6 subplots
function sPlot(df::DataFrames.DataFrame,sub_var::Symbol, x_var::Symbol, x_or::Array{String,1}, y_var::Symbol,y_lab::String,out_n::String)
    f, axes = Seaborn.subplots(2, 3, figsize=(10,7), sharey="all")
    x = y = 2
    plot_num = 1
    # order=["HG","Past.","Agri.","NA"],
    for loc in unique(df[sub_var])
        life_plot = swarmplot(df[df[sub_var] .== loc,x_var],df[df[sub_var] .== loc,y_var],
                        order=x_or,ax=axes[x, y])
        life_plot.set_title("$(loc) (N = $(nrow(df[df[sub_var] .== loc,:])))")
        life_plot.set_ylabel(y_lab)
        plot_num += 1
        x = plot_num%2 + 1
        y = plot_num%3 + 1
    end
    Seaborn.savefig("$out_n.pdf")
    clf()
end

function summaryMaps(df::DataFrames.DataFrame)
    R"""
    library(readr)
    library(ggplot2)
    library(maps)
    sources <- $df
    world <- map_data("world")
    world <- world[world$region != "Antarctica",]
    europe <- c("Albania","Andorra","Austria","Belarus","Belgium","Bosnia and Herzegovina","Bulgaria","Croatia","Cyprus","Czech Republic",
            "Denmark","Estonia","Finland","France","","Germany","Greece","Hungary","Iceland","Ireland","Italy","Kosovo","Latvia",
            "Liechtenstein","Lithuania","Luxembourg","Macedonia","Malta","Moldova","Monaco","Montenegro","Netherlands","Norway","Poland",
            "Portugal","Romania","San Marino","Serbia","Slovakia","Slovenia","Spain","Sweden","Switzerland","Ukraine","UK",
            "Vatican")

    # world, all
    ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group),fill = "grey", color = "black") + theme_void() +
               coord_fixed(1.3) + scale_color_gradient(low="blue",high="red",limits=c(min(sources$date),max(sources$date))) +
               geom_point(data=sources[sources$continent != "Europe"& sources$type != "WGS",], aes(x=longitude,y=latitude,color=date),position = "jitter",shape=1,size=2.5)+
               geom_point(data=sources[sources$continent != "Europe"& sources$type == "WGS",], aes(x=longitude,y=latitude,color=date),position = "jitter",size=1.5)
    ggsave("world_map.pdf",width = 8, height = 6,units="in",device="pdf")

    #europe, all
    ggplot() + geom_polygon(data = subset(world, region %in% europe), aes(x=long, y = lat, group = group),fill = "grey", color = "black") +
        theme_void() + coord_fixed(1.3) + scale_color_gradient(low="blue",high="red",limits=c(min(sources$date),max(sources$date))) +
            geom_point(data=sources[sources$continent == "Europe"& sources$type != "WGS",], aes(x=longitude,y=latitude,color=date),position = "jitter",shape=1,size=2.5)+
            geom_point(data=sources[sources$continent == "Europe"& sources$type == "WGS",], aes(x=longitude,y=latitude,color=date),position = "jitter",size=1.5)
    ggsave("europe_map.pdf",width = 8, height = 6,units="in",device="pdf")
    """
end

function eurasiaMaps(df::DataFrames.DataFrame)
    R"""
    library(readr)
    library(ggplot2)
    library(maps)
    sources <- $df
    world <- map_data("world")
    world <- world[world$region != "Antarctica",]
    eurasia <- c("Albania","Andorra","Austria","Belarus","Belgium","Bosnia and Herzegovina","Bulgaria","Croatia","Cyprus","Czech Republic",
        "Denmark","Estonia","Finland","France","","Germany","Greece","Hungary","Iceland","Ireland","Italy","Kosovo","Latvia",
        "Liechtenstein","Lithuania","Luxembourg","Macedonia","Malta","Moldova","Monaco","Montenegro","Netherlands","Norway","Poland",
        "Portugal","Romania","San Marino","Serbia","Slovakia","Slovenia","Spain","Sweden","Switzerland","Ukraine","UK",
        "Vatican","Afghanistan","Armenia","Azerbaijan","Bahrain","Bangladesh","Bhutan","Brunei","Cambodia","China","Georgia","India","Iran","Iraq","Israel","Japan",
         "Jordan","Kazakhstan","Kuwait","Kyrgyzstan","Laos","Lebanon","Mongolia","Myanmar","Nepal","North Korea","Oman","Pakistan","Palestine","Philippines","Qatar","Russia",
         "Saudi Arabia","South Korea","Sri Lanka","Syria","Tajikistan","Thailand","Turkey","Turkmenistan","United Arab Emirates","Uzbekistan","Vietnam","Yemen")
     cbPalette <- c("#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#E69F00", "#56B4E9", "#999999")

     #eurasia, by lifestyle
     ggplot() + geom_polygon(data = subset(world, region %in% eurasia), aes(x=long, y = lat, group = group),color = "black",fill="grey") + theme_void() +
         scale_colour_manual(values=cbPalette) + coord_fixed(1.3) +
         geom_point(data=sources[(sources$continent == "Europe" | sources$continent == "Asia") & sources$lifestyle == "Agri.",], aes(x=longitude,y=latitude,color=lifestyle),position = "jitter",size=1.5,shape=3) +
         geom_point(data=sources[(sources$continent == "Europe" | sources$continent == "Asia") & sources$lifestyle == "Past.",], aes(x=longitude,y=latitude,color=lifestyle),position = "jitter",size=1.5,shape=3) +
         geom_point(data=sources[(sources$continent == "Europe" | sources$continent == "Asia") & sources$lifestyle == "HG",], aes(x=longitude,y=latitude,color=lifestyle),position = "jitter",size=1.5,shape=3)
     ggsave("reich_3q_eurasia_lifestyle.pdf",width = 8, height = 6,units="in",device="pdf")
    """
end

function sampleSummary()
    anno_file = "../../data/ancient_dna/reich_compilation/v42/v42.4.1240K.anno"
    map_file =  "data/reich_v42/mapped_sample_ids.csv"
    info_file = "data/aDNA_sources.csv"
    anno = parseAnnoFile(anno_file)
    anno = anno[anno[:assessment] .!= "no",:]
    println("Number of Samples at Start: $(nrow(anno))")
    anno = mapSampleIDs(anno,map_file,info_file)
    println("Number of Samples successfully mapped: $(nrow(anno))")

    anno[:prop_miss] = 1.0 .- (anno[:num_snps]/N_SNPS_CHIP)
    CSV.write("reich_all_QCpass.txt",anno[[:ind_id,:master_id,:date,:country,:sex,:coverage,:num_snps,:prop_miss,:continent,:colbran_id,:type,:location,:population,:lifestyle,:latitude,:longitude]];delim='\t')

    anno[anno[:lifestyle] .== "nomadic/pastoral",:lifestyle] = "Past."
    anno[anno[:lifestyle] .== "pastoral",:lifestyle] = "Past."
    anno[anno[:lifestyle] .== "Agriculture",:lifestyle] = "Agri."
    anno[anno[:lifestyle] .== "Hunter-gatherer",:lifestyle] = "HG"
    describe(anno[:num_snps])
    nsnps_plot = histogram(anno[:num_snps],xlabel="Number SNPs Genotyped",ylabel="Number Samples",bins=100,margin=10Plots.mm)
    Plots.savefig(nsnps_plot,"num_snps_hist_all.pdf")
    prop_plot = histogram(anno[:prop_miss],xlabel="Proportion SNPs Missing",ylabel="Number Samples",bins=100,margin=10Plots.mm)
    Plots.savefig(prop_plot,"prop_miss_hist_all.pdf")
    describe(anno[:date])
    summaryMaps(anno)
    println(countmap(anno[:lifestyle]))
    println(countmap(anno[findall([x == "Europe" || x == "Asia" for x in anno[:continent]]),:lifestyle]))

    sPlot(anno,:continent,:lifestyle,["HG","Past.","Agri.","NA"],:date,"Age (yBP)","lifestyle_all")
    # for ones with >771029 SNPs called (3rd quartile)
    anno = anno[anno[:num_snps].>771029,:]
    println(countmap(anno[:lifestyle]))
    println(countmap(anno[findall([x == "Europe" || x == "Asia" for x in anno[:continent]]),:lifestyle]))

    describe(anno[:num_snps])
    describe(anno[:date])
    eurasiaMaps(anno)
    sPlot(anno,:continent,:lifestyle,["HG","Past.","Agri.","NA"],:date,"Age (yBP)","lifestyle_3q")
    anno = anno[findall([x == "Europe" || x == "Asia" for x in anno[:continent]]),:]
    nsnps_plot = histogram(anno[:num_snps],xlabel="Number SNPs Genotyped",ylabel="Number Samples",bins=100,margin=10Plots.mm)
    Plots.savefig(nsnps_plot,"num_snps_hist_3rdqeurasia.pdf")
    date_plot = histogram(anno[:date],xlabel="Date (yBP)",ylabel="Number Samples",bins=100,margin=10Plots.mm)
    Plots.savefig(date_plot,"date_hist_3rdqeurasia.pdf")
    CSV.write("reich_3rdq_eurasians.txt",anno[[:ind_id,:master_id,:publication,:date,:country,:sex,:coverage,:num_snps,:prop_miss,:continent,:colbran_id,:type,:location,:population,:lifestyle,:latitude,:longitude]];delim='\t')

end

sampleSummary()

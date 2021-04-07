# locus_zoom.jl
# given a file with a stat and a bunch of SNPs, makes a locuszoom plot
# @author Laura Colbran
# julia 1.4

using ArgParse
using DataFrames,CSV
using Plots

POS_COL = 2

# parses command-line arguments
function parseCommandline()
	s = ArgParseSettings()
	@add_arg_table! s begin
        "--data","-d"
            help = "path to file with stat. Assumes 2nd col is position."
            arg_type=String
            required = true
        "--column","-c"
            help = "1-indexed column containing the stat to plot"
            arg_type=Int64
            required = true
        "--window","-w"
            nargs=2
            help = "Min Max position to plot"
            arg_type = Float64
            default = [-Inf,Inf]
    	"--log", "-l"
            help = "if you want to plot the negative log of the stat"
            action=:store_true
        "--target","-t"
            help = "position of SNP to highlight"
            arg_type = Int64
            default = -1
    end
    return parse_args(s)
end

function zoomPlot(data_path::String,col::Int64,window::Array{Float64,1},l::Bool,target::Int64)
    data = CSV.read(data_path, DataFrame; delim='\t')
    data = data[(data[:,POS_COL] .>= window[1]) .& (data[:,POS_COL] .<= window[2]),:]
    data[!,:grp] .= "not"
    # data[!,col] = parse.(Float64, data[:,col])
    if l
        data[!,col] = -log.(data[:,col])
    end
    if target > 0
        data[(data[!,2] .== target),:grp] .= "target"
    end
    qq = Plots.scatter(data[!,POS_COL],data[!,col],xlabel="Coordinate",ylabel = "$(names(data)[col])",
                    group = data[!,:grp],legend = false,xformatter = :plain)
	Plots.savefig("zoomplot.pdf")
end

function main()
    parsed_args = parseCommandline()

    zoomPlot(parsed_args["data"],parsed_args["column"],parsed_args["window"],
                parsed_args["log"],parsed_args["target"])
end

main()

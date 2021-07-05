# qqplot.jl
# given a group of p-values, plots a qqplot, and runs a genomic control
# @author Laura Colbran
# julia 1.4

using ArgParse
using DataFrames, CSV
using RCall, MultipleTesting, HypothesisTests
using Plots

DEG_F = 1

GENES = ["ENSG00000161091","ENSG00000167986","ENSG00000104044",
"ENSG00000128731","ENSG00000188467","ENSG00000077498","ENSG00000107165",
"ENSG00000164175","ENSG00000258839","ENSG00000101440","ENSG00000115138",
"ENSG00000080166","ENSG00000049130","ENSG00000157404","ENSG00000137265",
"ENSG00000176797","ENSG00000185664","ENSG00000162341","ENSG00000187098",
"ENSG00000197535","ENSG00000047579","ENSG00000069974","ENSG00000088812",
"ENSG00000143669","ENSG00000115648","ENSG00000166189","ENSG00000134160",
"ENSG00000151694","ENSG00000173157","ENSG00000146648","ENSG00000145244",
"ENSG00000112038","ENSG00000157168","ENSG00000173068","ENSG00000040531",
"ENSG00000136160","ENSG00000124205"] #skin pigmentation genes


# parses command-line arguments
function parseCommandLine()
	s = ArgParseSettings()
	# @add_arg_table! s begin
	@add_arg_table! s begin
		"--p_file", "-f"
			help="path to file with p-values"
			arg_type = String
			required = true
		"--out", "-o"
            help = "path to write output files"
            arg_type = String
            default = "./"
		"--column","-c"
			help = "1-indexed column containing p-values to plot"
			arg_type=Int64
			required = true
        "--gc"
            help = "if you want to do genomic control"
            action=:store_true
		"--multi"
			help = "if you want to  FDR correct p-values"
			action=:store_true
		end
	return parse_args(s)
end

function genomicControl(p_values::Array{Float64,1})
    @rput p_values
    @rput DEG_F
    R"""
        chi_stat <- qchisq(p_values,df = DEG_F, lower.tail=F)
        lambda <- median(chi_stat)/qchisq(0.5, df=DEG_F) #inflation factor
        corrected_p <- pchisq(chi_stat/lambda, df=DEG_F, lower.tail=F)
    """
    @rget corrected_p
    return corrected_p
end

function qqPlot(path::String,col::Int64,out::String,gc::Bool,multi::Bool)
    results = CSV.read(path, DataFrame; delim='\t')
    rename!(results,col => :pval)
    results[isnan.(results.pval), :pval] .= 1 #only happens if literally all the predictions are the exact same value (4 genes in melanocytes)
    sort!(results,:pval)
    results[!,:exp] .= 1.0

    if gc
        results[!,:corr_pval] = genomicControl(results[!,:pval])
    end

	if multi
		if gc
			results[!,:bh_Q] = adjust(results[!,:corr_pval],BenjaminiHochberg())
		else
			results[!,:bh_Q] = adjust(results[!,:pval],BenjaminiHochberg())
		end
	end

    for i in 1:nrow(results)
        results[i,:exp] = -log10(i/nrow(results))
    end
    x = [0,maximum(results[!,:exp])]

	println(ApproximateTwoSampleKSTest(results[[x in GENES for x in results[:,1]],:corr_pval],results[[!(x in GENES) for x in results[:,1]],:corr_pval]))
	println(pvalue(ApproximateTwoSampleKSTest(results[[x in GENES for x in results[:,1]],:corr_pval],results[[!(x in GENES) for x in results[:,1]],:corr_pval])))

    qq = Plots.plot(x,x,color = :grey,xlabel="-log10(Expected P)",ylabel = "-log10(observed P)",margin=7Plots.mm,grid=false) #,lims=(0,maximum(-log10.(results[!,:corr_pval])))
	qq = Plots.scatter!(results[[!in(x,GENES) for x in results[:,1]],:exp],-log10.(results[[!in(x,GENES) for x in results[:,1]],:corr_pval]),
            legend = false,color = :black, alpha = 0.5,markersize=3)
    qq = Plots.scatter!(results[[x in GENES for x in results[:,1]],:exp],-log10.(results[[x in GENES for x in results[:,1]],:corr_pval]),
            legend = false,color = :dodgerblue, alpha = 1,markershape=:vline,markersize=9)
    Plots.savefig("$(out)qqplot.pdf")
	CSV.write("$(out)corrected.txt",results;delim="\t")
end

function main()
    parsed_args = parseCommandLine()
    # println(typeof(parsed_args["p_file"]),typeof(parsed_args["col"]),typeof(parsed_args["out"]),typeof(parsed_args["gc"]))
    qqPlot(parsed_args["p_file"],parsed_args["column"],parsed_args["out"],parsed_args["gc"],parsed_args["multi"])
end

main()

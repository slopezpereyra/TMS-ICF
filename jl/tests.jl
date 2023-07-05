using HypothesisTests
using CSV 
using DataFrames
using EffectSizes
using Statistics
using Plots
using StatsPlots

df = CSV.read("JLRA_DF.csv", DataFrame)
df = filter(:ISI => x -> x != 0, df)

function label_filter(label::String, df)
    filter([:Label, :Subject] => (x, y) -> x == label && 
                                                y != 8, 
                                                df)
end

function run_dependent_sample_tests(df::DataFrame)
    HC = [label_filter("HC BL", df), label_filter("HC SWD", df), "HC"]
    MDD = [label_filter("MDD BL", df), label_filter("MDD SWD", df), "MDD"]

    for group in [HC, MDD]
        println("~~~~~~~~~~~~~~~~~~~~~~  ", group[3], "  ~~~~~~~~~~~~~~~~~~~~~~~")
        println(SignTest(group[1].EMGPeakToPeak, group[2].EMGPeakToPeak))
        cohen_display(CohenD(group[1].EMGPeakToPeak, group[2].EMGPeakToPeak))
    end
end

function run_independent_sample_tests(df::DataFrame)
    A = [label_filter("HC BL", df), label_filter("MDD BL", df), "HC BL vs MDD BL"]
    B = [label_filter("HC SWD", df), label_filter("MDD SWD", df), "HC SWD vs MDD SWD"]
    display(A[1])
    display(B[1])

    for group in [A, B]
        println("~~~~~~~~~~~~~~~~~~~~~~  ", group[3], "  ~~~~~~~~~~~~~~~~~~~~~~~")
        println(MannWhitneyUTest(group[1].EMGPeakToPeak, group[2].EMGPeakToPeak))
    end
end

run_independent_sample_tests(mepdf)

function run_tests(df::DataFrame)
    HC_BL = label_filter("HC BL", df); HC_SWD = label_filter("HC SWD", df)
    MDD_BL = label_filter("MDD BL", df); MDD_SWD = label_filter("MDD SWD", df)
   
    println("----------------------------------------------------")
    println("Running tests on independent samples ~ HCs vs. MDDs Baseline")
    println("----------------------------------------------------\n\n")
    run_independent_sample_tests(df)
    println("----------------------------------------------------")
    println("Running tests on dependent samples ~ HC BL vs HC SWD")
    println("----------------------------------------------------\n")
    run_dependent_sample_tests(df)
end

run_independent_sample_tests(mepdf)

# Peculiarities: 
#
# On ISI 15: HC BL and MDD BL did differ.
# On ISI 8 HC BL and HC SWD did not differ.
#
#!/usr/bin/env julia


function run_tests_by_isi(df::DataFrame)
    res = DataFrame()
    res.pvalue = []; res.mean = []; res.median = []; res.test_level =[];
    res.isi = []; res.test = String[]; res.cohen = []; res.cohenci = [];

    
    for isi in sort(unique(df.ISI))
        idf = filter(:ISI => x -> x == isi, df)
        HC_BL = label_filter("HC BL", idf); HC_SWD = label_filter("HC SWD", idf)
        MDD_BL = label_filter("MDD BL", idf); MDD_SWD = label_filter("MDD SWD", idf)
        bl_wilcox = MannWhitneyUTest(HC_BL.EMGPeakToPeak, MDD_BL.EMGPeakToPeak)
        swd_wilcox = MannWhitneyUTest(HC_SWD.EMGPeakToPeak, MDD_SWD.EMGPeakToPeak)
        sign_hc = SignTest(HC_BL.EMGPeakToPeak, HC_SWD.EMGPeakToPeak)
        sign_mdd = SignTest(MDD_BL.EMGPeakToPeak, MDD_SWD.EMGPeakToPeak)
        cohen_hc = CohenD(HC_BL.EMGPeakToPeak, HC_SWD.EMGPeakToPeak)
        cohen_mdd = CohenD(MDD_BL.EMGPeakToPeak, MDD_SWD.EMGPeakToPeak)
        push!(res, ( 
            test = "HC BL vs MDD BL",
            test_level = 0,
            pvalue = pvalue(bl_wilcox), 
            isi = isi, 
            mean = (mean(HC_BL.EMGPeakToPeak), mean(MDD_BL.EMGPeakToPeak)),
            median = (median(HC_BL.EMGPeakToPeak), median(MDD_BL.EMGPeakToPeak)),
            cohen = missing, 
            cohenci = missing
            ))
        push!(res, ( 
            test = "HC SWD vs MDD SWD",
            pvalue = pvalue(swd_wilcox), 
            test_level = 1,
            isi = isi, 
            mean = (mean(HC_SWD.EMGPeakToPeak), mean(MDD_SWD.EMGPeakToPeak)),
            median = (median(HC_SWD.EMGPeakToPeak), median(MDD_SWD.EMGPeakToPeak)),
            cohen = missing,
            cohenci = missing
            ))
        push!(res, ( 
            test = "HC BL vs HC SWD",
            test_level = 2,
            pvalue = pvalue(sign_hc), 
            isi = isi, 
            mean = (mean(HC_BL.EMGPeakToPeak), mean(HC_SWD.EMGPeakToPeak)),
            median = (median(HC_BL.EMGPeakToPeak), median(HC_SWD.EMGPeakToPeak)),
            cohen = effectsize(cohen_hc), 
            cohenci = round.(ci(cohen_hc).ci, digits=3),
            ))
        push!(res, ( 
            test = "MDD BL vs MDD SWD",
            test_level = 3,
            pvalue = pvalue(sign_mdd), 
            isi = isi, 
            mean = (mean(MDD_BL.EMGPeakToPeak), mean(MDD_SWD.EMGPeakToPeak)),
            median = (median(MDD_BL.EMGPeakToPeak), median(MDD_SWD.EMGPeakToPeak)),
            cohen = effectsize(cohen_mdd), 
            cohenci = round.(ci(cohen_mdd).ci, digits=3)
            ))
    end
    return res
end

function get_grouped_summary(df::DataFrame)
    groupby(df, [:ISI, :Label])
    # Group the DataFrame by the 'group' column and compute summary statistics
    summary_df = combine(groupby(df, [:ISI, :Label])) do grouped_df
        DataFrame(
            mean_EMGPeakToPeak = mean(grouped_df.EMGPeakToPeak),
            median_EMGPeakToPeak = median(grouped_df.EMGPeakToPeak),
            σ²_EMGPeakToPeak = var(grouped_df.EMGPeakToPeak),
            σ_EMGPeakToPeak = sqrt(var(grouped_df.EMGPeakToPeak)),
            min_EMGPeakToPeak = minimum(grouped_df.EMGPeakToPeak),
            max_EMGPeakToPeak = maximum(grouped_df.EMGPeakToPeak)
        )
    end
end

function plot_effect_sizes(df::DataFrame)
    df = filter([:cohen, :pvalue] => (x, y) -> !ismissing(x) && y < 0.05, df)
    display(df)
    x = df[!, [:isi, :test_level, :cohen]]
    x = unstack(x, :test_level, :isi, :cohen)
    replace!(x[!, 5], missing => 0); replace!(x[!, 4], missing => 0);
    display(x)
    heatmap(Matrix(x[!, 2:end]), yflip=true)
    #scatter(df.isi, df.cohen, group=df.test, legend=:topleft, markersize=6)
end


df = read_swipf()
rename!(df, [:Sample, :Subject, :Group, :Session, :Type, :Label, :ISI, :EMGPeakToPeak])
mepdf = filter(:ISI => x -> ismissing(x), df)

run_tests(mepdf)

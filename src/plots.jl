using Plots
using DataFrames
using Measurements
using CSV
using OrderedCollections



units = OrderedDict(
    "ω_n" => "1/year",
    "ζ" => 1,
    "T_∞" => "K",
    "σ" => "K",
)

var_desc = Dict(
    "σ" => "noise level",
    "ζ" => "damping ratio",
    "ω_n" => "natural frequency",
    "T_∞" => "equilibrium surface temperature anomaly"
)


df_results = CSV.read("tas_fit_all.csv", DataFrame)
df_results = DataFrames.unstack(df_results, :gcm, :parameters, :value)
for v in keys(units)
    df_results[!,v] = measurement.(df_results[!,v])
end


function plot_all()
    # remove when CSV has been cleaned
    #df_results = df_results[df_results[!,:gcm] .!= "Year",:]

    subplots = []
    for var_name in keys(units)
        n_gcms = size(df_results, 1)
        p = scatter(
            1:n_gcms, df_results[!,Symbol(var_name)], 
            xrotation=90, title="$(var_name) ($(var_desc[var_name]))",
            xticks=(1:n_gcms, df_results[!,:gcm]),
            label=nothing, ylabel="$(var_name) [$(units[var_name])]"
        )
        if var_name == "ω_n"
            plot!(twinx(p), 1.0 ./ yticks(p)[1][1], ylabel="years", xticks=:none, label=:none, alpha=0.0)
        end
        push!(subplots, p)
    end
    plot(subplots..., layout=(2,2), size=(1400, 1000), margin=12Plots.mm, left_margin=10Plots.mm, bottom_margin=16Plots.mm, plot_title="Damped-spring fit to all GCM models")
end

p = plot_all()
savefig(p, "tas_fit_all.png")


function zeta_frequency_plot()
    xvar = "ω_n"
    yvar = "ζ"
    scatter(
        #Measurements.value.(df_results[!,:ω_n]), 
        #Measurements.value.(df_results[!, :σ]), 
        df_results[!,xvar],
        df_results[!,yvar], 
        xlabel="$xvar [$(units[xvar])]", ylabel="$(yvar) [$(units[yvar])]",
        label=:none
    )

    df_results_ann = df_results[df_results.ω_n .> 0.1, :]
    annotate!(
        Measurements.value.(df_results_ann[!,xvar]),
        Measurements.value.(df_results_ann[!, yvar]), Plots.text.(df_results_ann.gcm, 10, :bottom, )
    )

    df_results_ann = df_results[df_results.ζ .> 4.0, :]
    annotate!(
        Measurements.value.(df_results_ann[!,xvar]),
        Measurements.value.(df_results_ann[!, yvar]), Plots.text.(df_results_ann.gcm, 10, :left, :bottom)
    )
end

p = zeta_frequency_plot()
savefig(p, "tas_fit_all_zeta_vs_omega_n.png")

(df_results[!, :ω_n])
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



function plot_all()
    df_results = CSV.read("tas_fit_all.csv", DataFrame)
    df_results[!,:value] = measurement.(df_results[!,:value])
    # remove when CSV has been cleaned
    #df_results = df_results[df_results[!,:gcm] .!= "Year",:]

    subplots = []
    for var_name in keys(units)
        df_results_variable = df_results[df_results[!, :parameters] .== var_name,:]
        n_gcms = size(df_results_variable, 1)
        p = scatter(
            1:n_gcms, df_results_variable[!,:value], 
            xrotation=90, title="$(var_name) ($(var_desc[var_name]))", xticks=(1:n_gcms, df_results_variable[!,:gcm]),
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

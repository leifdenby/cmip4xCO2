using CSV
using DataFrames
using Plots
using Turing
using OrdinaryDiffEq
using StatsPlots
using Measurements


function damped_spring(du, u, p, t)
    # v = dT/dt, the "velocity" of the temperature
    (T, v) = u
    (ζ, ω_n, T_∞) = p
    # dx/dt
    du[1] = v
    # d^2x/dt^2 = dv/dt
    du[2] = -2*ζ*ω_n*v - ω_n^2 * (T - T_∞)
end


@model function temperature_equilibration(years, tas)
    # priors
    #ζ ~ Beta(1, 1)
    #ω0 ~ Beta(1, 1)
    ζ ~ TruncatedNormal(1.0, 1.0, 0.0, Inf)
    ω_n ~ TruncatedNormal(1.0, 1.0, 0.0, Inf)
    σ ~ TruncatedNormal(0.0, 1.0, 0.0, Inf)
    T_∞ ~ Normal(tas[end], 1.0)  # tas as t -> ∞

    # initial condition
    v0 = (tas[2] - tas[1])/(years[2] - years[1])
    u0 = [tas[1], v0]

    # solution
    tspan = (minimum(years), maximum(years))
    prob = ODEProblem(damped_spring, u0, tspan, (ζ, ω_n, T_∞))
    prediction = solve(prob, Tsit5())

    # posterior
    for i in 1:length(years)
        tas[i] ~ Normal(prediction(years[i])[1], σ)
    end
    return nothing
end

function fit_spring_model_and_plot(df, gcm_model_name)
    t_years = df[!,"Year"]
    tas = df[!, gcm_model_name]

    #chain = sample(temperature_equilibration(t_years, tas), MH(), 10)
    chain = sample(temperature_equilibration(t_years, tas), NUTS(), MCMCThreads(), 500, 12, progress=true);

    p1 = plot(chain)
    p2 = plot(; legend=false)
    posterior_samples = sample(chain[[:ζ, :ω_n, :T_∞]], 300; replace=false)
    for p in eachrow(Array(posterior_samples))
        v0 = (tas[2] - tas[1])/(t_years[2] - t_years[1])
        prob = ODEProblem(damped_spring, [tas[1], v0], (minimum(t_years), maximum(t_years)), p)
        sol_p = solve(prob, Tsit5(); saveat=0.1)
        plot!(p2, sol_p; alpha=0.1, color="#BBBBBB")
    end
    plot!(p2, t_years, tas, xlabel="year", ylabel="surface air temperature anomaly [K]")

    plot(p1, p2, plot_title=gcm_model_name, layout=grid(2,1,heights=[0.7,0.3]), size=(800, 1200), left_margin=10Plots.mm)

    savefig("plots/tas_oscillator_fit_$(gcm_model_name).png")

    df_chain_summary = DataFrame(summarize(chain))
    # make a column which encodes mean ± std as the estimate
    insertcols!(df_chain_summary, :est => df_chain_summary[!,:mean] .± df_chain_summary[!,:std])
    # and add the gcm model name as a column
    insertcols!(df_chain_summary, :gcm => gcm_model_name)
    df_results = DataFrames.stack(df_chain_summary[!, [:parameters, :est, :gcm]])
    return df_results
end

function main()
    df_results = nothing
    df = CSV.read("cmip6/delta_tas_abrupt-4xCO2_cmip6.csv", DataFrame)

    gcm_names = [name for name in names(df) if name != "Year"]

    for gcm_name in gcm_names
        df_results_gcm = fit_spring_model_and_plot(df, gcm_name)
        if df_results == nothing
            df_results = df_results_gcm
        else
            append!(df_results, df_results_gcm)
        end
    end
    df_results = DataFrames.unstack(df_results, :gcm, :parameters, :value)
    CSV.write("tas_fit_all.csv", df_results)
end

main()
using CSV
using DataFrames
using Plots
using Turing
using OrdinaryDiffEq
using StatsPlots

df = CSV.read("cmip6/delta_tas_abrupt-4xCO2_cmip6.csv", DataFrame)
names(df)


function damped_spring(du, u, p, t)
    (x, v) = u
    (ζ, ω_n, x0) = p
    # dx/dt
    du[1] = v
    # d^2x/dt^2 = dv/dt
    du[2] = -2*ζ*ω_n*v - ω_n^2 * (x - x0)
end


@model function temperature_equilibration(years, tas)
    # priors
    #ζ ~ Beta(1, 1)
    #ω0 ~ Beta(1, 1)
    ζ ~ TruncatedNormal(1.0, 1.0, 0.0, Inf)
    ω_n ~ TruncatedNormal(1.0, 1.0, 0.0, Inf)
    σ ~ TruncatedNormal(0.0, 1.0, 0.0, Inf)
    x0 ~ Normal(tas[end], 1.0)  # tas as t -> ∞

    # initial condition
    v0 = (tas[2] - tas[1])/(years[2] - years[1])
    u0 = [tas[1], v0]

    # solution
    tspan = (minimum(years), maximum(years))
    prob = ODEProblem(damped_spring, u0, tspan, (ζ, ω_n, x0))
    prediction = solve(prob, Tsit5())

    # posterior
    for i in 1:length(years)
        tas[i] ~ Normal(prediction(years[i])[1], σ)
    end
    return nothing
end

model_name = "INM-CM4-8"
t_years = df[!,"Year"]
tas = df[!, model_name]
#chain = sample(temperature_equilibration(t_years, tas), MH(), 10)
chain = sample(temperature_equilibration(t_years, tas), NUTS(), MCMCThreads(), 500, 4, progress=true);


p1 = plot(chain)
p2 = plot(; legend=false)
posterior_samples = sample(chain[[:ζ, :ω_n, :x0]], 300; replace=false)
for p in eachrow(Array(posterior_samples))
    v0 = (tas[2] - tas[1])/(t_years[2] - t_years[1])
    prob = ODEProblem(damped_spring, [tas[1], v0], (minimum(t_years), maximum(t_years)), p)
    sol_p = solve(prob, Tsit5(); saveat=0.1)
    plot!(p2, sol_p; alpha=0.1, color="#BBBBBB")
end
plot!(p2, t_years, tas, xlabel="year", ylabel="surface air temperature anomaly [K]")

plot(p1, p2, plot_title=model_name, layout=grid(2,1,heights=[0.7,0.3]), size=(800, 1200), left_margin=10Plots.mm)

savefig("tas_oscillator_fit_$(model_name).png")

## plot all
plot()
for name in names(df)[2:end]
    plot!(df[!,"Year"], df[!,name], label=name)
end
plot!()
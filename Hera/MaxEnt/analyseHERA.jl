using LensFactory
using JLD2
using CairoMakie
using LensFactory.LFUtils
using LaTexStrings

path = "/home/arnav/Parametric/Hera_outputs/MaxEnt/"
name = "Hera_MaxEnt_20261113"
N_constraints = 130     # 65 total images
println("N_constrain = $N_constraints")
clustername = "Hera"
if clustername == "SMACS0723"
    X_lim_plot = 50.0
    Y_lim_plot = 50.0
    kappa_range = (0, 2.0)
elseif clustername == "Hera"
    X_lim_plot = 70.0
    Y_lim_plot = 70.0
    kappa_range = (0, 2.70)
else
    X_lim_plot = 100.0
    Y_lim_plot = 100.0
    kappa_range = (0, 3.0)
end

last_flag = true

time_start = time()
# Read the JLD2 file
jldopen("$(path)$(name).jld2", "r") do data
    # Read data
    model = data["model"]
    results = data["optimizer"]
    chains = data["chains"]
    logL = data["logL"]
    last_kappa = data["optimizer"][:kappa]
    last_prior = data["optimizer"][:prior]
    runs = data["optimizer"][:runs]
    x_grid = data["optimizer"][:x]
    y_grid = data["optimizer"][:y]

    # Get the best fit lens model
    #=best_lens, best_F = LensModel.get_best_model("$(path)$(name).jld2", run = :best)
    println("Best target function is : $(best_F)")=#

    last_lens, last_F = LensModel.get_best_model("$(path)$(name).jld2", run = :last)
    println("Last target function is : $(last_F)")

    # Common quantities
    FOV = model.observation.FOV
    pixel_scale = model.observation.pixel_scale
    param_ref = Dict(p.key => p.refer for p in model.parameters)
    cosmo     = LensModel.LensModelUtils.current_cosmology(model, param_ref)
    zs_ = 9
    zd_ = model.observation.z_d
    cosmo = Cosmology.init_cosmology()      # default cosmo
    Dds_ = Cosmology.angular_diameter_distance(cosmo, zd_, zs_)
    Ds_ = Cosmology.angular_diameter_distance(cosmo, 0.0, zs_)
    adis_9 = Dds_ / Ds_

    last_kappa9 = copy(last_kappa) * adis_9
    last_prior9 = copy(last_prior) * adis_9

    # plots for last model
    cc_fig, cc_axes = Lenses.plot_image_plane(last_lens, x_grid, y_grid, adis_9, two_panel = true, figure_size = (1200,600))        # saving critical curves
    xlims!(cc_axes[1], -X_lim_plot, X_lim_plot)
    ylims!(cc_axes[1], -Y_lim_plot, Y_lim_plot)
    xlims!(cc_axes[2], -X_lim_plot, X_lim_plot)
    ylims!(cc_axes[2], -Y_lim_plot, Y_lim_plot)
    save("$(path)$(name)_last_critical_curves.png", cc_fig)

    sd_fig, sd_axes = Lenses.plot_sky(0.5*FOV[1], 0.5*FOV[2], figure_size = (500, 400))
    hm = heatmap!(sd_axes, x_grid[:,1], y_grid[1,:], last_prior9, colormap = :turbo, colorrange = kappa_range)
    cb = Colorbar(fig[1,2], hm; label = L"κ_{prior}")
    xlims!(sd_axes, -X_lim_plot, X_lim_plot)
    ylims!(sd_axes, -Y_lim_plot, Y_lim_plot)
    save("$(path)$(name)_last_prior_map.png", sd_fig)
    
    sd_fig, sd_axes = Lenses.plot_sky(0.5*FOV[1], 0.5*FOV[2], figure_size = (500, 400))
    hm = heatmap!(sd_axes, x_grid[:,1], y_grid[1,:], last_kappa9, colormap = :turbo, colorrange = kappa_range)
    cb = Colorbar(fig[1,2], hm; label = "κ")
    xlims!(sd_axes, -X_lim_plot, X_lim_plot)
    ylims!(sd_axes, -Y_lim_plot, Y_lim_plot)
    save("$(path)$(name)_last_kappa_map.png", sd_fig)

    mg_fig, mg_axes = Lenses.plot_magnification_map(last_lens, x_grid, y_grid, adis_9, heatmap_kws = (colormap=:turbo, colorrange=(0, 100)))
    xlims!(mg_axes, -X_lim_plot, X_lim_plot)
    ylims!(mg_axes, -Y_lim_plot, Y_lim_plot)
    save("$(path)$(name)_last_magnification_map.png", mg_fig)

    # Plot runs stats
    n = eachindex(runs)
    logpost = getproperty.(runs, :log_posterior)
    loglike = getproperty.(runs, :log_likelihood)
    logprior = getproperty.(runs, :log_prior)
    iterations = getproperty.(runs, :iterations)
    change = getproperty.(runs, :change)

    fig = Figure(size = (1400, 400))
    ax1 = Axis(fig[1,1],
        title = "Log probabilities",
        xlabel = "Run",
        ylabel = "Log value")
    ax2 = Axis(fig[1,2],
        title = "Optimizer iterations",
        xlabel = "Run",
        ylabel = "Iterations")
    ax3 = Axis(fig[1,3],
        title = "Relative map change",
        xlabel = "Run",
        ylabel = "Change")

    lines!(ax1, n, logpost, label = "Log-posterior")
    lines!(ax1, n, loglike, label = "Log-likelihood")
    lines!(ax1, n, logprior, label = "Log-prior")
    axislegend(ax1, position = :rt)

    scatterlines!(ax2, n, iterations)
    scatterlines!(ax3, n, change)

    save("$(path)$(name)_runstats.png", fig)

    println("plots saved.")
    println("Calculating best fit rms...")
    flush(stdout)
    t0 = time()
    LensModel.get_best_fit_rms(model, last_kappa)
    println("Time taken to calculate best fit rms: $(time() - t0) seconds")
    println("----------------------------------------")
    flush(stdout)

    # plots for best model

    # plots for first model

    # plotting each individual image system predictions
    n_sources = length(model.source_config.sources)
    plot_sids = collect(1:n_sources)
    for sid in plot_sids
        LensModel.plot_best_model(model, optimizer; two_panel = true, source=sid, plot_name="$(path)$(name)_last_model_$(sid).png")
    end
end

time_end = time()

println("Total time taken: $(time_end - time_start) seconds")

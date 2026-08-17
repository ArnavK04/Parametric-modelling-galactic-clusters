using LensFactory
using JLD2
using CairoMakie

path = "/home/arnav/Parametric/SMACSJ0723_outputs/2026-08-14/"
name = "SMACSJ0723_2026-08-14"
N_constraints = 98     # 49 total images

time_start = time()
# Read the JLD2 file
jldopen("$(path)$(name).jld2", "r") do data
    # Read data
    model = data["model"]
    results = data["optimizer"]
    chains = data["chains"]
    logL = data["logL"]
	
    free_params = LensModel.free_parameter_names(model)
    t0 = time()
    println("Free parameters loaded: $free_params")
    println("----------------------------------------")
    best_theta, log_post, lower_err, upper_err = LensModel.get_best_fit_parameters(logL, chains=chains, burn_in=0.2, with_errors=true, thin=1, print_table=true, free_parameter_names = free_params)
    println("Best fit parameters loaded: $best_theta")
    println("Time taken to load best fit parameters: $(time() - t0) seconds")
    println("----------------------------------------")
    flush(stdout)
    # Get rms
    println("Calculating best fit rms...")
    t0 = time()
    LensModel.get_best_fit_rms(model, chains, logL)
    println("Time taken to calculate best fit rms: $(time() - t0) seconds")
    println("----------------------------------------")
    flush(stdout)
    # get aic and bic
    AIC = LensModel.get_AIC(model, logL)
    BIC = LensModel.get_BIC(model, logL, N_constraints)
    println("AIC: $AIC")
    println("BIC: $BIC")

    # MCMC diagnostics
    println("----------------------------------------")
    t0 = time()
    println("MCMC diagnostics:")
    LensModel.print_gr_report(chains, free_param_names = free_params, burn_in=0.4)
    LensModel.autocorrelation(chains, free_param_names=free_params, burn_in=0.4)
    LensModel.acceptance_rate(chains, burn_in=0.4)
    println("Time taken for MCMC diagnostics: $(time() - t0) seconds")
    println("----------------------------------------")
    flush(stdout)

    # Model plots
    println("Plotting modelling results...")
    t0 = time()
    fig_corner = LensModel.plot_corner(chains, logL, free_parameter_names=free_params, burn_in=0.3, plot_name="$(path)$(name)_corner_plot.png")
    fig_trace = LensModel.plot_trace(chains, free_parameter_names=free_params, burn_in=0.3, plot_name="$(path)$(name)_trace_plot.png")
    fig_model = LensModel.plot_best_model(model, chains, logL, z_s = 1.5, plot_name="$(path)$(name)_best_model.png")
    fig_scatter, ax_scatter = LensModel.plot_image_scatter(model, chains, logL, plot_name="$(path)$(name)_image_scatter.png")
    println("Time taken to plot modelling results: $(time() - t0) seconds")
    println("Done plotting modelling results.")
    flush(stdout)

    # plotting each individual image system predictions
    n_sources = length(model.source_config.sources)
    plot_sids = collect(1:n_sources)
    for sid in plot_sids
        LensModel.plot_best_model(model, chains, logL, source=sid, plot_name="$(path)$(name)_best_model_$(sid).png")
    end

    cosmo_best = LensModel.get_cosmology(data; burn_in=0.2, thin=1, with_errors=false)
    LensModel.save_best_fits(data)

    println("fits saved")

end

time_end = time()

println("Total time taken: $(time_end - time_start) seconds")

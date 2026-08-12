using LensFactory
using JLD2
using CairoMakie

path = "/home/arnav/Parametric/SMACSJ0723_outputs/"
name = "SMACSJ0723_2026-08-11"
N_constraints = 98     # 49 total images

# Read the JLD2 file
jldopen("$(path)$(name).jld2", "r") do data
    # Read data
    model = data["model"]
    results = data["optimizer"]
    chains = data["chains"]
    logL = data["logL"]
	
    free_params = LensModel.free_parameter_names(model)
    println("Free parameters loaded: $free_params")
    println("----------------------------------------")
    best_theta, log_post, lower_err, upper_err = LensModel.get_best_fit_parameters(logL, chains=chains, burn_in=0.2, with_errors=true, thin=1, print_table=true, free_param_names = free_params)
    println("----------------------------------------")
    # Get rms
    LensModel.get_best_fit_rms(model, chains, logL)
    println("----------------------------------------")
    # get aic and bic
    AIC = get_AIC(model, chains, logL)
    BIC = get_BIC(model, chains, logL, N_constraints)
    println("AIC: $AIC")
    println("BIC: $BIC")

    # MCMC diagnostics
    println("----------------------------------------")
    println("MCMC diagnostics:")
    LensModel.print_gr_report(chains, free_param_names = free_params, burn_in=0.4)
    LensModel.autocorrelation(chains, free_param_names=free_params, burn_in=0.4)
    LensModel.acceptance_rate(chains, burn_in=0.4)
    println("----------------------------------------")

    # Model plots
    println("Plotting modelling results...")
    fig_corner = LensModel.plot_corner(chains, logL, free_parameter_names=free_params, burn_in=0.3, plot_name="$(path)$(name)_corner_plot.png")
    fig_trace = LensModel.plot_trace(chains, free_parameter_names=free_params, burn_in=0.3, plot_name="$(path)$(name)_trace_plot.png")
    fig_model = LensModel.plot_best_model(model, chains, logL, z_s = 1.5, plot_name="$(path)$(name)_best_model.png")
    fig_scatter, ax_scatter = LensModel.plot_image_scatter(model, chains, logL, plot_name="$(path)$(name)_image_scatter.png")

    println("Done plotting modelling results.")

    cosmo_best = LensModel.get_cosmology(data; burn_in=0.2, thin=1, with_errors=false)
    save_best_fits(data)

    println("fits saved")

end
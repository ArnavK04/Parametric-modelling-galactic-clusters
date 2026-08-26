using LensFactory
using JLD2
using CairoMakie
using LensFactory.LFUtils

path = "/home/arnav/Parametric/SMACSJ0723_outputs/2026-08-11/"
name = "SMACSJ0723_2026-08-11"
N_constraints = 98     # 49 total images
println("N_constrain = $N_constraints")

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

    # Get the best fit lens model
    best_lens, best_logL = LensModel.get_best_model(model, mcmc_chains = chains, mcmc_logL = logL)
    println("Best log-likelihood in source plane is : $(best_logL)")

    # Construct grid
    FOV = model.observation.FOV
    pixel_scale = model.observation.pixel_scale
    x_grid, y_grid = Lenses.get_meshgrid(0.5 * FOV[1], 0.5 * FOV[2], pixel_scale)

    param_ref = Dict(p.key => p.refer for p in model.parameters)
    pvals     = LensModel.LensModelUtils.param_dict(model, best_theta, param_ref)
    cosmo     = LensModel.LensModelUtils.current_cosmology(model, pvals)
    adis      = LensModel.LensModelUtils.adis_current(model, param_ref, cosmo)
    ax_all, ay_all    = LensModel.LensModelUtils.lens_quantities_def(model, best_lens)
    A_all     = LensModel.LensModelUtils.lens_quantities_jac(model, best_lens)

    logL_img, β_mod_s, θ_mod_s, all_converged = LensModel.Likelihood.logL_imageplane_fast(model, best_lens, adis, ax_all, ay_all, A_all)
    println("Log-likelihood for the best fit model in image plane is: $(logL_img)")
    println("All sources converged: $(all_converged)")
    flush(stdout)

    sid = 5
    kid = 1

    src = model.source_config.sources[sid]
    adis_value = adis[sid]
    z_s = 1.425400  # get this from the images files for source 5

    knot = src.knots[kid]
    x = knot.x
    y = knot.y
    images_obs = [(xi, yi) for (xi, yi) in zip(x, y)]
    if size(x, 1) > 1
        # Get deflection at the image positions
        αx, αy = Lenses.get_deflection(best_lens, x, y)

        # Get magnification at the image positions
        μ_obs = Lenses.get_magnification_image(best_lens, x, y, adis_value)

        # Calculate individual image source positions
        βx_indi = x - adis_value * αx
        βy_indi = y - adis_value * αy

        # Calculate barycenter source position
        βx_model = sum(βx_indi .* μ_obs.^2) / sum(μ_obs.^2)
        βy_model = sum(βy_indi .* μ_obs.^2) / sum(μ_obs.^2)
    else
        αx, αy = Lenses.get_deflection(best_lens, x, y)
        βx_model = x - adis_value * αx
        βy_model = y - adis_value * αy
    end
    images = Lenses.get_image(best_lens, x_grid, y_grid, adis_value, (βx_model, βy_model))
    # get the lens redshift and add
    z_d = model.observation.z_d
    D_d = Cosmology.angular_diameter_distance(cosmo, 0.0, z_d)

    td_obs = Lenses.get_time_delay(best_lens, x, y, adis_value, z_d, D_d, (βx_model, βy_model))
    td_model = Lenses.get_time_delay(best_lens, first.(images), last.(images), adis_value, z_d, D_d, (βx_model, βy_model))
    reltd_obs = td_obs .- minimum(td_obs)
    reltd_model = td_model .- minimum(td_model)
    println("Observed images are at $(images_obs)")
    println("Model predicted images are at $(images)")
    println("Observed time delays are $(reltd_obs) seconds or $((reltd_obs)./86400) days")
    println("Model time delays are $(reltd_model) seconds or $((reltd_model)./86400) days")
    println("source redshift is $z_s")
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
    fig_model = LensModel.plot_best_model(model, chains, logL, two_panel = true, z_s = 1.5, plot_name="$(path)$(name)_best_model.png")
    fig_scatter, ax_scatter = LensModel.plot_image_scatter(model, chains, logL, plot_name="$(path)$(name)_image_scatter.png")
    println("Time taken to plot modelling results: $(time() - t0) seconds")
    println("Done plotting modelling results.")
    flush(stdout)

    # plotting each individual image system predictions
    n_sources = length(model.source_config.sources)
    plot_sids = collect(1:n_sources)
    for sid in plot_sids
        LensModel.plot_best_model(model, chains, logL, two_panel = true, source=sid, plot_name="$(path)$(name)_best_model_$(sid).png")
    end

    cosmo_best = LensModel.get_cosmology(data; burn_in=0.2, thin=1, with_errors=false)
    LensModel.save_best_fits(data)

    println("fits saved")

end

time_end = time()

println("Total time taken: $(time_end - time_start) seconds")

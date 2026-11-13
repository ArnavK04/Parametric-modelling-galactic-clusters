using LensFactory
using JLD2
using CairoMakie
using LensFactory.LFUtils
using LensFactory.Constants
using LensFactory.LensModel.LensModelIO
using LaTeXStrings
using Interpolations
using FITSIO
using LinearAlgebra

path = "/home/arnav/Parametric/Hera_outputs/"
name = "Hera_Parametric_20261110"
N_constraints = 130     # 65 total images
println("N_constrain = $N_constraints")
clustername = "Hera"
if clustername == "SMACS0723"
    X_lim_plot = 50.0
    Y_lim_plot = 50.0
    kappa_range = (0, 2.0)
elseif clustername == "Hera"
    X_lim_plot = 50.0
    Y_lim_plot = 50.0
    kappa_range = (0, 2.70)
else
    X_lim_plot = 100.0
    Y_lim_plot = 100.0
    kappa_range = (0, 3.0)
end
kappa_dev_range = (-2.0, 2.0)
mag_dev_range = (-4.0, 4.0)

function load_fitsfile(filename::String)

    local gridx_fits, gridy_fits, kappa, gamma1, gamma2

    FITS("/home/arnav/FreeFormModelling/$(filename)_data/kappa_z9_0.fits") do f
        hdr = read_header(f[1])
        kappa = read(f[1])

        NX = hdr["NAXIS1"]
        NY = hdr["NAXIS2"]

        pixel_scale = hdr["CDELT1"] * 3600.0  # Convert degrees to arcseconds
        CRVAL1 = hdr["CRVAL1"] * 3600.0
        CRVAL2 = hdr["CRVAL2"] * 3600.0
        CRPIX1 = hdr["CRPIX1"]
        CRPIX2 = hdr["CRPIX2"]

        X_max = (NX - CRPIX1) * pixel_scale + CRVAL1
        Y_max = (NY - CRPIX2) * pixel_scale + CRVAL2

        gridx_fits, gridy_fits = Lenses.get_meshgrid(X_max, Y_max, pixel_scale)
    end

    FITS("/home/arnav/FreeFormModelling/$(filename)_data/gammax_z9_0.fits") do f
        gamma1 = read(f[1])
    end

    FITS("/home/arnav/FreeFormModelling/$(filename)_data/gammay_z9_0.fits") do f
        gamma2 = read(f[1])
    end

    return gridx_fits[2:end, 2:end], gridy_fits[2:end, 2:end], kappa, gamma1, gamma2    # because the fits arrays are 2048, but grid_fits are 2049
end

function refine_map(map_coarse::M, gridx::M, gridy::M, new_grid_limx::T, new_grid_limy::T, resolution::T, order::Int) where {T <:RV, M <: ROA}

    x_nodes = range(gridx[1, 1], stop=gridx[end, 1], length=size(gridx, 2))
    y_nodes = range(gridy[1, 1], stop=gridy[1, end], length=size(gridy, 1))

    if order == 1
        itp = interpolate(map_coarse, BSpline(Linear()))
    elseif order == 2
        itp = interpolate(map_coarse, BSpline(Quadratic(Line(OnGrid()))))
    elseif order == 3
        itp = interpolate(map_coarse, BSpline(Cubic(Line(OnGrid()))))
    else
        error("Unsupported interpolation order $order. upto cubic(3) is allowed.")
    end
    
    itp = scale(itp, x_nodes, y_nodes)

    x_fine, y_fine = Lenses.get_meshgrid(new_grid_limx, new_grid_limy, resolution)

    map_fine = itp.(x_fine, y_fine)

    return map_fine, x_fine, y_fine

end


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
    k = length(free_params)
    println("Number of free parameters: $k")
    println("Free parameters loaded: $free_params")
    println("----------------------------------------")
    best_theta, log_post, lower_err, upper_err = LensModel.get_best_fit_parameters(logL, chains=chains, burn_in=0.2, with_errors=true, thin=1)
    println("Best fit parameters loaded: $best_theta")
    println("Time taken to load best fit parameters: $(time() - t0) seconds")
    println("----------------------------------------")
    flush(stdout)

    # Get the best fit lens model
    best_lens, best_logL = LensModel.get_best_model(model, chains = chains, mcmc_logL = logL)
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

    zs_ = 9
    zd_ = model.observation.z_d
    cosmo = Cosmology.init_cosmology()      # default cosmo
    Dds_ = Cosmology.angular_diameter_distance(cosmo, zd_, zs_)
    Ds_ = Cosmology.angular_diameter_distance(cosmo, 0.0, zs_)
    adis_9 = Dds_ / Ds_

    cc_fig, cc_axes = Lenses.plot_image_plane(best_lens, x_grid, y_grid, adis_9, two_panel = true)        # saving critical curves
    xlims!(cc_axes[1], -X_lim_plot, X_lim_plot)
    ylims!(cc_axes[1], -Y_lim_plot, Y_lim_plot)
    xlims!(cc_axes[2], -X_lim_plot, X_lim_plot)
    ylims!(cc_axes[2], -Y_lim_plot, Y_lim_plot)
    save("$(path)$(name)_critical_curves.png", cc_fig)

    sd_fig, sd_axes = Lenses.plot_surface_density(best_lens, x_grid, y_grid, adis_9, unit = :convergence, heatmap_kws = (colormap=:turbo, colorrange=kappa_range))
    xlims!(sd_axes, -X_lim_plot, X_lim_plot)
    ylims!(sd_axes, -Y_lim_plot, Y_lim_plot)
    save("$(path)$(name)_kappa_map.png", sd_fig)   

    mg_fig, mg_axes = Lenses.plot_magnification_map(best_lens, x_grid, y_grid, adis_9, heatmap_kws = (colormap=:turbo, colorrange=(0, 100)))
    xlims!(mg_axes, -X_lim_plot, X_lim_plot)
    ylims!(mg_axes, -Y_lim_plot, Y_lim_plot)
    save("$(path)$(name)_magnification_map.png", mg_fig)

    logL_img, β_mod_s, θ_mod_s, all_converged = LensModel.Likelihood.logL_imageplane_fast(model, best_lens, adis, ax_all, ay_all, A_all)
    println("Log-likelihood for the best fit model in image plane is: $(logL_img)")
    println("All sources converged: $(all_converged)")
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
    println("# Free parameters: $k")
    println("# Constraints: $N_constraints")
    println("SOURCE PLANE AIC and BIC:")
    println("AIC: $AIC")
    println("BIC: $BIC")
    println("IMAGE PLANE  AIC and BIC:")
    img_AIC = -2.0 * logL_img + 2.0 * k
    img_BIC = -2.0 * logL_img + k * log(N_constraints)
    println("AIC: $img_AIC")
    println("BIC: $img_BIC")
    flush(stdout)

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
    fig_model = LensModel.plot_best_model(model, chains, logL, two_panel = true, z_s = 9.0, plot_name="$(path)$(name)_best_model.png")      # best model at 9
    fig_scatter, ax_scatter = LensModel.plot_image_scatter(model, chains, logL, plot_name="$(path)$(name)_image_scatter.png", point_kws = (markersize  = 5, 
                                                    marker      = :circle, 
                                                    color       = :transparent, 
                                                    strokecolor = :black, 
                                                    strokewidth = 2))
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

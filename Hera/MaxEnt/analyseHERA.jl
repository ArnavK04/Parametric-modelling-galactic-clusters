using LensFactory
using LensFactory.Constants
using LensFactory.LensModel.LensModelIO
using JLD2
using CairoMakie
using LensFactory.LFUtils
using LaTeXStrings
using Interpolations
using FITSIO
using LinearAlgebra

include("/home/arnav/FreeFormModelling/FreeForm-modelling-LensFactory/utility_functions.jl")

path = "/home/arnav/Parametric/Hera_outputs/Hera_MaxEnt_20261113/"
name = "Hera_MaxEnt_20261113"
N_constraints = 130     # 65 total images
println("N_constrain = $N_constraints")
clustername = "Hera"
plot_individual_src = false

if clustername == "SMACS0723"
    X_lim_plot = 50.0
    Y_lim_plot = 50.0
    kappa_range = (0, 2.0)
elseif clustername == "Hera"
    X_lim_plot = 70.0
    Y_lim_plot = 70.0
    X_lim = 100.0
    Y_lim = 100.0
    kappa_range = (0, 2.70)
else
    X_lim_plot = 100.0
    Y_lim_plot = 100.0
    kappa_range = (0, 3.0)
end
kappa_dev_range = (-2.0, 2.0)
mag_dev_range = (-4.0, 4.0)

last_flag = true

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

function get_time_delays_all(model, lens, adis, x_grid, y_grid, cosmo; match_tol = 2.0)

    zd = model.observation.z_d
    Dd = Cosmology.angular_diameter_distance(cosmo, 0.0, zd)

    all_delays = Vector{Vector{Float64}}()      # one entry per knot, in observed-image order

    for (sid, src) in enumerate(model.source_config.sources)
        adis_value = adis[sid]

        for knot in src.knots
            β_knot = LensModel._get_source_position(lens, knot, adis_value)
            images = Lenses.get_image(lens, x_grid, y_grid, adis_value, β_knot)

            delays = fill(NaN, length(knot.x))

            if !isempty(images)
                imx = [p[1] for p in images]
                imy = [p[2] for p in images]
                t = Lenses.get_time_delay(lens, imx, imy, adis_value, zd, Dd, β_knot)

                used = falses(length(imx))      # predicted images already matched

                for j in eachindex(knot.x)
                    d = hypot.(imx .- knot.x[j], imy .- knot.y[j])
                    d[used] .= Inf              # exclude already-matched images
                    k = argmin(d)
                    if d[k] <= match_tol        # also false if everything is Inf
                        delays[j] = t[k]
                        used[k] = true          # only consume it on a real match
                    else
                        @warn "No match found for knot $j of source $sid (min distance = $(d[k]))"
                    end
                end

                # relative to the earliest-arriving matched image
                valid = filter(!isnan, delays)
                isempty(valid) || (delays .-= minimum(valid))
            end

            push!(all_delays, delays)
        end
    end

    return all_delays
end

time_start = time()
# Read the JLD2 file
jldopen("$(path)$(name).jld2", "r") do data
    # Read data
    model = data["model"]
    optimizer = data["optimizer"]
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

    println("Calculating best fit rms...")
    flush(stdout)
    t0 = time()
    #LensModel.get_best_fit_rms(model, last_kappa, method = :adaptive)
    println("Time taken to calculate best fit rms: $(time() - t0) seconds")
    println("----------------------------------------")
    flush(stdout)

    last_lens, last_F = LensModel.get_best_model("$(path)$(name).jld2", run = :last)
    println("Last target function is : $(last_F)")

    # Common quantities
    FOV = model.observation.FOV
    pixel_scale = model.observation.pixel_scale
    X_lim = 0.5*FOV[1]
    Y_lim = 0.5*FOV[2]
    param_ref = Dict(p.key => p.refer for p in model.parameters)
    cosmo     = LensModel.LensModelUtils.current_cosmology(model, param_ref)
    adis      = LensModel.LensModelUtils.adis_current(model, param_ref, cosmo)
    zs_ = 9
    zd_ = model.observation.z_d
    Dds_ = Cosmology.angular_diameter_distance(cosmo, zd_, zs_)
    Ds_ = Cosmology.angular_diameter_distance(cosmo, 0.0, zs_)
    adis_9 = Dds_ / Ds_

    last_kappa9 = copy(last_kappa) * adis_9
    last_prior9 = copy(last_prior) * adis_9

    # plots for last model
    crt_curv_plot = time()
    cc_fig, cc_axes = Lenses.plot_image_plane(last_lens, x_grid, y_grid, adis_9, two_panel = true, figure_size = (1200,600))        # saving critical curves
    xlims!(cc_axes[1], -X_lim_plot, X_lim_plot)
    ylims!(cc_axes[1], -Y_lim_plot, Y_lim_plot)
    xlims!(cc_axes[2], -X_lim_plot, X_lim_plot)
    ylims!(cc_axes[2], -Y_lim_plot, Y_lim_plot)
    save("$(path)$(name)_last_critical_curves.png", cc_fig)
    println("Critical curves plotted in $(time() - crt_curv_plot) sec")

    sd_fig, sd_axes = Lenses.plot_sky(0.5*FOV[1], 0.5*FOV[2], figure_size = (500, 400))
    hm = heatmap!(sd_axes, x_grid[:,1], y_grid[1,:], last_prior9, colormap = :turbo, colorrange = kappa_range)
    cb = Colorbar(sd_fig[1,2], hm; label = L"κ_{prior}")
    xlims!(sd_axes, -X_lim_plot, X_lim_plot)
    ylims!(sd_axes, -Y_lim_plot, Y_lim_plot)
    save("$(path)$(name)_last_prior_map.png", sd_fig)
    
    sd_fig, sd_axes = Lenses.plot_sky(0.5*FOV[1], 0.5*FOV[2], figure_size = (500, 400))
    hm = heatmap!(sd_axes, x_grid[:,1], y_grid[1,:], last_kappa9, colormap = :turbo, colorrange = kappa_range)
    cb = Colorbar(sd_fig[1,2], hm; label = "κ")
    xlims!(sd_axes, -X_lim_plot, X_lim_plot)
    ylims!(sd_axes, -Y_lim_plot, Y_lim_plot)
    save("$(path)$(name)_last_kappa_map.png", sd_fig)

    mag_plot = time()
    last_mag9 = Lenses.get_magnification_image(last_lens, x_grid, y_grid, adis_9)

    mag_fig, mag_axes = Lenses.plot_sky(0.5*FOV[1], 0.5*FOV[2], figure_size = (500, 400))
    hm = heatmap!(mag_axes, x_grid[:,1], y_grid[1,:], abs.(last_mag9), colormap = :turbo, colorrange = (0,100))
    cb = Colorbar(mag_fig[1,2], hm; label = L"|\mu|")
    xlims!(mag_axes, -X_lim_plot, X_lim_plot)
    ylims!(mag_axes, -Y_lim_plot, Y_lim_plot)
    save("$(path)$(name)_last_magnification_map.png", mag_fig)
    println("magnification maps plotted in $(time() - mag_plot) sec")

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
        ylabel = "Change",
        yscale = log10)

    lines!(ax1, n, logpost, label = "Log-posterior")
    lines!(ax1, n, loglike, label = "Log-likelihood")
    lines!(ax1, n, logprior, label = "Log-prior")
    axislegend(ax1, position = :rb)

    scatterlines!(ax2, n, iterations)
    scatterlines!(ax3, n, change)

    save("$(path)$(name)_runstats.png", fig)

    # plot deviation from the truth

    if clustername == "Hera" || clustername == "Ares"

        res = abs(x_grid[1,1] - x_grid[2,1])
        order = 3

        gridx_fits, gridy_fits, kappa, gamma1, gamma2 = load_fitsfile(clustername)
        Hera_offset = 0.0      # Hera image coordinates were offset by 12 in x direction, so we need to correct for that
        if clustername == "Hera"
            gridx_fits .+= Hera_offset
        end
        kappa = Float64.(kappa)  # Ensure kappa is of type Float64
        gamma1 = Float64.(gamma1)  # Ensure gamma1 is of type Float64
        gamma2 = Float64.(gamma2)  # Ensure gamma2 is of type Float64

        flush(stdout)
        println("interpolating...")
        # interpolate to a particular grid

        refinestart = time()
        kappa_finefits, x_fine, y_fine = refine_map(kappa, gridx_fits, gridy_fits, X_lim, Y_lim, res, order)
        gamma1_finefits, _, _ = refine_map(gamma1, gridx_fits, gridy_fits, X_lim, Y_lim, res, order)
        gamma2_finefits, _, _ = refine_map(gamma2, gridx_fits, gridy_fits, X_lim, Y_lim, res, order)
        magstrt = time()
        println("interpolation done in ", magstrt - refinestart, " seconds.")
        mag_finefits = @. 1.0 / ((1.0 - kappa_finefits)^2 - (gamma1_finefits^2 + gamma2_finefits^2))
        println("magnification calculated in ", time() - magstrt, " seconds.")
        flush(stdout)

        sd_fig, sd_axes = Lenses.plot_sky(0.5*FOV[1], 0.5*FOV[2], figure_size = (500, 400))
        hm = heatmap!(sd_axes, x_grid[:,1], y_grid[1,:], (last_kappa9 .- kappa_finefits)./kappa_finefits, colormap = :BrBG, colorrange = kappa_dev_range)
        cb = Colorbar(sd_fig[1,2], hm; label = L"\delta \kappa/\kappa_{truth}")
        xlims!(sd_axes, -X_lim_plot, X_lim_plot)
        ylims!(sd_axes, -Y_lim_plot, Y_lim_plot)
        save("$(path)$(name)_last_kappareldev_map.png", sd_fig)

        mag_fig, mag_axes = Lenses.plot_sky(0.5*FOV[1], 0.5*FOV[2], figure_size = (500, 400))
        hm = heatmap!(mag_axes, x_grid[:,1], y_grid[1,:], (last_mag9 .- mag_finefits)./mag_finefits, colormap = :BrBG, colorrange = mag_dev_range)
        cb = Colorbar(mag_fig[1,2], hm; label = L"\delta \mu/\mu_{truth}")
        xlims!(mag_axes, -X_lim_plot, X_lim_plot)
        ylims!(mag_axes, -Y_lim_plot, Y_lim_plot)
        save("$(path)$(name)_last_magreldev_map.png", mag_fig)

        maghist_fig = Figure(size=(500,400))
        maghist_ax = Axis(maghist_fig[1,1], xlabel = L"\delta \mu/\mu_{truth}", ylabel = "PDF")
        hist!(maghist_ax, vec((last_mag9 .- mag_finefits)./mag_finefits), normalization = :pdf)
        save("$(path)$(name)_last_magreldev_hist.png", maghist_fig)

        kappahist_fig = Figure(size=(500,400))
        kappahist_ax = Axis(kappahist_fig[1,1], xlabel = L"\delta \kappa/\kappa_{truth}", ylabel = "PDF")
        hist!(kappahist_ax, vec((last_kappa9 .- kappa_finefits)./kappa_finefits), normalization = :pdf)
        save("$(path)$(name)_last_kappareldev_hist.png", kappahist_fig)

        # plot radial profile
        if clustername == "Ares"
            centers1, profile1, edges1 = Lenses.get_radial_profile(last_kappa9, x_grid, y_grid, origin = (-20.0, -32.0))    
            centers2, profile2, edges2 = Lenses.get_radial_profile(last_kappa9, x_grid, y_grid, origin = (40.0, 40.0))
        elseif clustername == "Hera"
            centers1, profile1, edges1 = Lenses.get_radial_profile(last_kappa9, x_grid, y_grid, origin = (0.0, 0.0))    
            centers2, profile2, edges2 = Lenses.get_radial_profile(last_kappa9, x_grid, y_grid, origin = (-19.5408, 2.3782))
        end

        fig1 = Figure(size=(500,400))
        ax1 = Axis(fig1[1,1], xlabel = "r (\")", ylabel = L"\kappa")
        lines!(ax1, centers1, profile1)
        save("$(path)$(name)_last_radialprofile1.png", fig1)
        fig2 = Figure(size=(500,400))
        ax2 = Axis(fig2[1,1], xlabel = "r (\")", ylabel = L"\kappa")
        lines!(ax2, centers2, profile2)
        save("$(path)$(name)_last_radialprofile2.png", fig2)

        time_delays = get_time_delays_all(model, last_lens, adis, x_grid, y_grid, cosmo)
        
        #save all important arrays.
        jldsave("$(path)$(name)_finalarrays.jld2";
            kappa9 = last_kappa9,
            prior9 = last_prior9,
            mag9 = last_mag9,
            kappa_finefits = kappa_finefits,
            mag_finefits = mag_finefits,
            time_delays = time_delays,
            centers1 = centers1,
            profile1 = profile1,
            centers2 = centers2,
            profile2 = profile2
        )
    end

    println("plots saved.")

    # plots for best model

    # plots for first model
    if plot_individual_src
        # plotting each individual image system predictions
        n_sources = length(model.source_config.sources)
        plot_sids = collect(1:n_sources)
        for sid in plot_sids
            LensModel.plot_best_model(model, optimizer; two_panel = true, source=sid, plot_name="$(path)$(name)_last_model_$(sid).png")
        end
    end
end

time_end = time()

println("Total time taken: $(time_end - time_start) seconds")

include(joinpath(@__DIR__, "common.jl"))
using Random
using RectiGrids
using AxisKeysExtra
import ScatteringOpticsExtra as SO

gaussfits_fullband = deserialize(joinpath(INTERMEDIATE_DIR, "gaussfits_fullband.jls"))

fig = Figure(backgroundcolor=:black, figure_padding=0)
map(Real[1.4, 2.3, 5]u"GHz" |> enumerate) do (i, ν)
	λ = u"c" / ν |> u"cm"

	fit = argmin(r -> abs(r.freq - ν), gaussfits_fullband)
	@info "Broadening model" requested=ν used=fit.freq
	sky = @modify(MCM.pmedian, $(fit.sky) |> RecursiveOfType(MCM.Particles))

	# intrinsic source model:
	sourcemodel = construct(CircularGaussian, flux => flux(sky), fwhm_max => intrinsic_source_fwhm(ν), coords => (0,0))
	broadening = construct(EllipticGaussian, flux => NaN, fwhm_max => fwhm_max(sky), fwhm_min => fwhm_min(sky), position_angle => position_angle(sky), coords => (0,0))
	# intrinsic source image:
	imggrid = let
		fov = 40u"mas"
		npix = 1024
		grid(SVector, X=LinRange(0u"mas"±fov, npix), Y=LinRange(0u"mas"±fov, npix))
	end

	screen_spec = SO.StochasticPhaseScreen()  # pure turbulence

	sm = SO.ScatteringModel(;
		target_distance=1u"Gpc", screen_distance=1.4u"kpc",
		r_in=1000u"km", α=5/3,
		broadening, λ₀=λ, screen_spec
	)

	Random.seed!(Random.default_rng(), round(Int, 123 * ustrip(ν)))

	sourceimg = intensity(sourcemodel).(imggrid) .|> ustrip;
	scatimg = SO.scatter_image(sm, sourceimg; λ)
	scatimg = @modify(ak -> u"mas".(ak), axiskeys(scatimg)[∗])

	w = Dict(
		1.4u"GHz" => 250,
		1.8u"GHz" => 190,
		2.3u"GHz" => 120,
		5u"GHz" => 40,
	)[ν] * 1.4
	lim = 35
	ax,plt = image(fig[1,i], reverse(scatimg; dims=1); colormap=:inferno, axis=(;width=w, height=350, limits=(0±lim*w/350, 0±lim)))
	hidedecorations!()
	text!((0.5,1),space=:relative,align=(0.5,1),offset=(0,-2), text="$ν", color=:white)

	p = poly!(InterferometricModels.convolve(sourcemodel, broadening)|>ustrip, color=:transparent, strokecolor=:white, strokewidth=2, label="Observed FWHM size")
	translate!(p,0,0,100)
end
scalebar!(content(fig[1,2]), 1u"mas"; color=:white, muls=[15], position=(0.25, 0.07))
colgap!(fig.layout, 7)
resize_to_layout!()
save(joinpath(FIGS_DIR, "simimage.pdf"), fig; backend=CairoMakie)
@info "Saved simimage.pdf"

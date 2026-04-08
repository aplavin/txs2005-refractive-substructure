include(joinpath(@__DIR__, "common.jl"))

using RectiGrids
using AxisKeysExtra
using StatsBase: quantile
import ScatteringOpticsExtra as SO
import VLBISkyModels as VSM

gaussfits_fullband = deserialize(joinpath(INTERMEDIATE_DIR, "gaussfits_fullband_post.jls"))

isfast = false

function generate_scatterimages(ν, n)
	λ = u"c" / ν |> u"cm"

	fit = argmin(r -> abs(r.freq - ν), gaussfits_fullband)
	@info "Broadening model" requested=ν used=fit.freq
	sky = @modify(MCM.pmedian, $(fit.sky) |> RecursiveOfType(MCM.Particles))

	sourcemodel = construct(CircularGaussian, flux => fit.model_flux, fwhm_max => intrinsic_source_fwhm(ν), coords => (0,0))
	broadening = construct(EllipticGaussian, flux => NaN, fwhm_max => fwhm_max(sky), fwhm_min => fwhm_min(sky), position_angle => position_angle(sky), coords => (0,0))

	mul = isfast ? 3 : 7
	screen_spec = SO.StochasticPhaseScreen(pad=@o mul*_)

	sm = SO.ScatteringModel(;
		target_distance=1u"Gpc", screen_distance=1.4u"kpc",
		r_in=1000u"km", α=5/3,
		broadening, λ₀=λ, screen_spec
	)

	fov = 120u"mas" / (ν/1.4u"GHz")
	npix = 1024
	imggrid = grid(SVector, X=LinRange(0u"mas"±fov, npix), Y=LinRange(0u"mas"±fov, npix))
	pixarea = abs(prod(step, axiskeys(imggrid)))

	sourceimg = intensity(sourcemodel).(imggrid) * pixarea .|> ustrip
	fovmas = ustrip(u"mas", fov)
	screen_shift = Observable(SVector(0., 0.)u"mas")
	scatimg_obs = SO.scatter_image(sm, sourceimg; λ, screen_shift)
	shifts = map(1:n) do _
		SVector(rand(0±((mul-1)*fovmas)), rand(0±((mul-1)*fovmas))) * u"mas"
	end
	obsmap(screen_shift, shifts, scatimg_obs)
end

function generate_scatteruvtbls(scatimgs, ν)
	λ = u"c" / ν |> u"cm"

	uvspecs = let maxuv = 0.8 * 12800u"km" / λ |> NoUnits
		@p grid(R=range(0..maxuv, length=100), θ=range(0..2π, length=100)) vec map((; spec=_.R * UV(sincos(_.θ)), freq_spec=ν))
	end

	map(scatimgs) do scatimg
		model = VSM.ContinuousImage(VSM.IntensityMap(scatimg), VSM.DeltaPulse())
		@p VSM.visibilitymap(model, uvspecs) sort(by=norm(VLBI.UV(_)))
	end
end

function uvtbls_to_uvtblband(uvtbls)
	maxuv = @p uvtbls[1] maximum(norm(UV(_)))
	uvmul = maxuv / 200
	@p let
		uvtbls
		stack
		eachslice(dims=1)
		map() do sl
			res = sl[1]
			vals = map(r -> abs(r.value), sl)
			@set res.value = Interval(abs.(quantile(vals, [0.16, 0.84]))...)
		end
		group_vg(round(norm(VLBI.UV(_)) / uvmul) * uvmul)
		sort(by=key)
		map() do gr
			absvalrng = @p gr map(_.value) minimum(minimum, __)..maximum(maximum, __)
			absvalrng = minimum(absvalrng) .. maximum(absvalrng)
			@p let
				gr[1]
				@set VLBI.UV(__) = VLBI.UV(key(gr), 0)
				@set __.value = absvalrng
			end
		end
	end
end

n = isfast ? 30 : 600
freqs_flux = [
	(ν=1.4u"GHz", color=ColorSchemes.seaborn_bright6[1]),
	(ν=2.3u"GHz", color=ColorSchemes.seaborn_bright6[2]),
	(ν=5u"GHz", color=ColorSchemes.seaborn_bright6[3]),
]

fig = Figure()
ax = Axis(fig[1,1]; 
		yscale=log10,
		limits=((0, 8e3), 0.002..4),
		yticks=BaseMulTicks([1,2,5]),
		xticks=WilkinsonTicks(10, k_min=6),
		# xtickformat=Makie.Formatters.plain,
		width=400, height=300,
		xlabel="Baseline length (km)",
		ylabel="Amplitude (Jy)")

# Remaining frequencies
for (; ν, color) in freqs_flux
	@info "Generating scattered images..." ν
	scatimgs = @showtime generate_scatterimages(ν, n)
	uvtbls = @showtime generate_scatteruvtbls(scatimgs, ν)
	uvtbl_band = uvtbls_to_uvtblband(uvtbls)
	fplt = FPlot(uvtbl_band, VLBIPlots.AxFuncs.UVdist_u(u"km"), @o _.value)
	bandstroke!(fplt, label="$ν", color=(color, 0.25), strokecolor=color)

	save(joinpath(FIGS_DIR, "scatmodel_uvamp.pdf"), current_figure(); backend=CairoMakie)
end

# VLBA baseline range
let bl_rng_km = 130..5200
	vspan!(bl_rng_km, ymax=0.06, color=(:black, 0.18))
	vlines!(collect(endpoints(bl_rng_km)), ymax=0.06, color=:gray)
	text!((@lift ($(MakieExtra.data2rel(current_axis(), :x, mean(bl_rng_km))), 0)), space=:relative, align=(:center,:bottom), offset=(0,2), text="Baselines VLBA covers in all directions")
end

axislegend("Scattering model\n(measured broadening,\nexpected substructure)", merge=true)
resize_to_layout!()

save(joinpath(FIGS_DIR, "scatmodel_uvamp.pdf"), current_figure(); backend=CairoMakie)
@info "Saved scatmodel_uvamp.pdf"

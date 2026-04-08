include(joinpath(@__DIR__, "common.jl"))

gaussfits_perif = deserialize(joinpath(INTERMEDIATE_DIR, "gaussfits_perif.jls"))
powerlaw_fits = deserialize(joinpath(INTERMEDIATE_DIR, "powerlaw_fits.jls"))

let
	fig = Figure()
	fplt = FPlot((@p gaussfits_perif sort(by=_.date)),
				 AxFunc(label="Frequency (GHz)", scale=log10, ticks=[1.4,1.8, 2.3, 5], @o ustrip(u"GHz", _.freq)),
				 AxFunc(label="Gaussian FWHM (mas)", scale=log10, ticks=BaseMulTicks(), @o fwhm_max(_.sky) |> ustrip);
				 color=AsCategorical(@o year(_.date)), markersize=12)

	multiplot(fig[1,1], (axplot(scatter), rangebars), fplt, marker=:utriangle; axis=(;height=300, width=400))
	multiplot!((scatter, rangebars), (@set fplt[2] = AxFunc(scale=log10, ticks=BaseMulTicks(), @o fwhm_min(_.sky) |> ustrip)), marker=:dtriangle)

	# Major axis power-law
	let (; fit, func) = powerlaw_fits.major
		# 2σ credible interval
		modvals = @p StructArray(freq=range(1..10, 100)u"GHz") mapinsert(fwhm=func(_, fit)[1][2] |> U.Value) map(@set U.uncertainty(_.fwhm) *= 2)
		multiplot!((lines, band => (;alpha=0.4)), FPlot(modvals, (@o ustrip(u"GHz", _.freq)), last); color=:black,linestyle=:solid, label=f"Major axis: θ ~ ν"*superscript(f"{U.Value(fit.k):.2f}"), to_xy_attrs(autolimits=false)...)
	end

	# Minor axis power-law
	let (; fit, func) = powerlaw_fits.minor
		modvals = @p StructArray(freq=range(1..10, 100)u"GHz") mapinsert(fwhm=func(_, fit)[1][2] |> U.Value) map(@set U.uncertainty(_.fwhm) *= 2)
		multiplot!((lines, band => (;alpha=0.4)), FPlot(modvals, (@o ustrip(u"GHz", _.freq)), last); color=:black,linestyle=:dash, label=f"Minor axis: θ ~ ν"*superscript(f"{U.Value(fit.k):.2f}"), to_xy_attrs(autolimits=false)...)
	end

	axislegend(merge=true, position=(:left,:bottom), orientation=:horizontal, nbanks=4)

	multiplot(fig[2,1], (axplot(scatter), rangebars), (@set fplt[2] = AxFunc(label="Major axis P.A.", ticks=Makie.AngularTicks(rad2deg(1), "°"), @o position_angle(_.sky))); axis=(;height=170))
	hlines!(NoUnits(gallat_PA), linestyle=:dash, color=:gray)
	text!((@lift ($(MakieExtra.rel2data(:x, 0)), NoUnits(gallat_PA))); text="Constant Galactic latitude", xautolimits=false, offset=(3,0))

	autohide_axlabels!(fig[1:2, 1])

	resize_to_layout!()
	save(joinpath(FIGS_DIR, "gaussfit_summary.pdf"), fig; backend=CairoMakie)
	@info "Saved figs/gaussfit_summary.pdf"
end

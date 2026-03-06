include(joinpath(@__DIR__, "common.jl"))

@info "Loading fit results..."
fits_perif = deserialize(joinpath(INTERMEDIATE_DIR, "gaussfits_perif.jls"))
fits_fullband_post = deserialize(joinpath(INTERMEDIATE_DIR, "gaussfits_fullband_post.jls"))

@info "Loading raw data..."
alldata = load_all_data()

perif_data = @p let
	alldata
	flatmap() do r
		@p r.uvtbl group_vg(_.freq_spec) map() do gr
			@p let
				r
				@set __.uvtbl = gr
				@insert __.freq_band = r.freq
				@reset __.freq = frequency(key(gr))
			end
		end
	end
end

@assert length(perif_data) == length(fits_perif) "Mismatch: $(length(perif_data)) per-IF data vs $(length(fits_perif)) fits"

gaussfits_perif = map(perif_data, fits_perif) do r, f
	@assert visfile_band(r.f) == f.band && visfile_date(r.f) == f.date "Mismatch: $(visfile_band(r.f))/$(visfile_date(r.f)) vs $(f.band)/$(f.date)"
	@insert r.sky = f.sky
end

@info "Generating radplot_gauss figures..."
@p let
	gaussfits_perif
	group_vg(_.f)
	map() do fgr
		fig = Figure()
		ax = Axis(fig[1,1]; height=270, width=350)
		model_flux = filterfirst(r -> r.band == visfile_band(key(fgr)) && r.date == visfile_date(key(fgr)), fits_fullband_post).model_flux
		uvtbl = @p fgr flatmap() do r
			mod = r.sky
			uvtbl_r = @p value(r.uvtbl) mapinsert(value_mod=U.Value(visibility(abs, mod, UV(_))))
			mod_vis = map(row -> U.value(row.value_mod), uvtbl_r)
			gains = VLBI.solve_gains(VLBI.AmplitudeAnalytic(), mod_vis, uvtbl_r)
			gains = @modify(g -> clamp(g, 1..4), gains.gain[∗])  # clamping shouldn't affect the results much...

			data = VLBI.apply_gains(VLBI.AmplitudeAnalytic(), gains, uvtbl_r) #; default_gain=1.)
			@p data filter(U.uncertainty(_.value) < 0.02)  # drop large errors, just for display
		end
		uvtbl = @modify(v -> v * model_flux, uvtbl.value[∗])
		freq = @p uvtbl mean(frequency(_.freq_spec))
		if freq < 1.6u"GHz"
			# for more familiar labeling
			freq = 1.4u"GHz"
		end
		λ = u"c"/freq
		maxuv = @p uvtbl maximum(norm(UV(_)))

		θ = @p fgr map(_.sky) mean(MCM.pmedian(fwhm_average(_)))
		srcmod = construct(CircularGaussian, flux => model_flux, fwhm_max => intrinsic_source_fwhm(freq), coords => (0,0))
		θ_scat = √(θ^2 - fwhm_average(srcmod)^2) |> u"mas"
		other_scat_params = (;D=1.4u"kpc", α=5/3, r_in=1_000u"km")
		qs = let p = 0.95
			Interval(quantile(Rayleigh(), [(1-p)/2, 1-(1-p)/2])...)
		end

		avgmodel = let
			props = @p fgr map(_.sky) StructArray getproperties() map(col -> mean(x -> x isa MCM.Particles ? MCM.pmedian(x) : x, col))
			model_flux * EllipticGaussian(;props...)
		end

		gausscolor = ColorSchemes.seaborn_bright[1]
		gausslabel = "Gaussian structure\nwith θ"*subscript("avg")*f" = {fwhm_average(avgmodel):.1f}" => (;markersize=15)

		multiplot!((axplot(scatter!) => (;markersize=3), rangebars => (;linewidth=1, color=(:black, 0.2))), RadPlot(uvtbl,); color=:black, axis=(;yscale=SymLog(0.02, linscale=0.5), limits=(nothing, 0..(1.1*model_flux)), ylabel="Amplitude (Jy)"), label="Data " * rich("(self-calibrated)"; font=:italic) => (;markersize=15))
		scatter!(FPlot(uvtbl, VLBIPlots.F.UVdist(), @o model_flux * abs(_.value_mod)); color=gausscolor, markersize=2.5, label=gausslabel)

		bandstroke!(RadPlot(range(0..100e6, 500), model=avgmodel); color=(gausscolor, 0.12), strokecolor=(gausscolor, 1), strokewidth=1, xautolimits=false, label=gausslabel)

		let
			color = ColorSchemes.seaborn_bright[2]
			uvdist_start = filterfirst(uv -> visibility(abs, avgmodel, UV(uv, 0)) < 0.25 * visibility(abs, avgmodel, UV(0, 0)), range(0, 1.5maxuv, 1000))
			bandstroke!((uvdist_start)..(1.1maxuv), uvdist -> Interval((endpoints(qs) .* scat_vis_noise(srcmod; uvdist, θ_scat, ν=freq, other_scat_params...))...); color=(color, 0.12), strokecolor=color, strokewidth=1, label=f"Substructure (expectations)\n95% range", xautolimits=false)
		end

		axislegend(f"{Unitful.GHz(freq):.1f} on {visfile_date(key(fgr))}", position=(:right,:top), merge=true)

		kmlims = @lift ustrip.(u"km", λ .* endpoints(intervals($(ax.finallimits))[1]))
		Axis(fig[1,1], xaxisposition=:top, xlabel="Baseline (km)", limits=(@lift ($kmlims, nothing)), xgridvisible=false)
		hideydecorations!()

		resize_to_layout!()
		outpath = joinpath(FIGS_DIR, f"radplot_gauss_{visfile_band(key(fgr))}_{visfile_date(key(fgr))}.pdf")
		save(outpath, fig; backend=CairoMakie)
		@info "Saved $outpath"
	end
	collect
end

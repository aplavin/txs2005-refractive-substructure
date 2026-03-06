include(joinpath(@__DIR__, "common.jl"))

all_fringefits = deserialize(joinpath(INTERMEDIATE_DIR, "all_fringefits.jls"))

# Derive per-band Gaussian models from fullband fits (for contour overlays)
gaussfits_fullband = deserialize(joinpath(INTERMEDIATE_DIR, "gaussfits_fullband.jls"))
gaussfits = Dict(
	letter => let
		sky = last(filter(r -> r.band == letter, gaussfits_fullband)).sky
		@modify(MCM.pmean, sky |> RecursiveOfType(MCM.Particles))
	end
	for letter in ["L", "S", "C"]
)

let
	fig = Figure()
	@p all_fringefits filter(any(r -> r.band isa AbstractVector, _.fringefits)) flatmap() do r
		@p r.fringefits group_vg(_.freqix) map() do gr
			@reset r.fringefits = gr
			@insert r.freqix = key(gr)
		end
	end filter(frequency(_.fringefits[1].band) < 6u"GHz" && !occursin("/BG246AK/", _.uvf_path) && !(_.freqix == 3 && occursin("/BG246AR/", _.uvf_path))) group_vg(_.freqix) map() do freqgr
		col = fig[key(freqgr), 1]
		bandletter = ["L", "S", "C"][key(freqgr)]

		freqint = @p freqgr flatmap(__ -> __.fringefits |> flatmap(frequency(_.band, Interval) |> endpoints)) extrema Interval(__...)
		dateint = @p freqgr map(_.date_obs |> year) extrema Interval(__...)

		freqstr = @p endpoints(freqint) map(f"{_ |> Unitful.GHz:.1f}") join(__, " – ")
		datestr = @p endpoints(dateint) map(string) join(__, " – ")
		Label(col[0,1:2], rich(f"""VLBA observations: {bandletter} band, {freqstr}, in {datestr}"""); tellwidth=false)

		# we have fewer than 1e4 scans across all experiments, so this corresponds to a false detection probability of less than 1 scan across all experiments
		isdetected = @o _.ntrials * exp(-_.peakval^2 / (2*_.noisedist.σ^2)) < 1e-4
		color = AsCategorical(isdetected, label=d -> (d ? "Detection" : "Non-detection") => (;markersize=15))

		radax = Axis(col[1,1], palette=(;color=[Makie.Gray(0.8), Makie.ColorSchemes.seaborn_bright[3]]), width=220, height=250)
		uvax = Axis(col[1,2], palette=(;color=[Makie.Gray(0.8), Makie.ColorSchemes.seaborn_bright[3]]), width=220, height=250)
		colgap!(col, 3)

		for (i, r) in enumerate(@p freqgr sort(by=_.date_obs) first(__, 20))
			fplt = FPlot(
				(@p r.fringefits filter(!allequal(_.ants)) filter(_.band isa VLBI.FrequencyWindow)),
				VLBIPlots.F.UVdist(),
				AxFunc(limit=(-0.01, nothing), scale=SymLog(5; linscale=0.2), label="Fringe S/N ratio", @o _.peakval/_.noisedist.σ);
				color, markersize=4,
			)
			combdata = @p r.fringefits filter(!allequal(_.ants)) filter(!(_.band isa VLBI.FrequencyWindow)) sort(by=_.peakval / _.noisedist.σ)
			axplot(scatter!)(radax, (@set fplt.data = combdata))

			axplot(scatter!)(uvax, UVPlot((@p VLBI.add_conjvis(combdata)); markersize=5, color); label=nothing)
		end

		axislegend(radax, "Fringe finding", merge=true, position=(:right, :top))
		gaussfit = gaussfits[bandletter]
		contour!(uvax, UVPlot(0±150e6; model=gaussfit); levels=[1e-3], color=:black, linestyle=:dash, linewidth=2, to_xy_attrs(autolimits=false)..., label="Scatter broadening:\n1/1000 of peak amplitude")
		axislegend(uvax; framevisible=false, position=(0.5, 1), padding=0)
	end
	resize_to_layout!()
	save(joinpath(FIGS_DIR, "fringeplots_summary.pdf"), fig; backend=CairoMakie)
	@info "Saved figs/fringeplots_summary.pdf"
end

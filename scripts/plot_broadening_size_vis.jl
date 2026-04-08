include(joinpath(@__DIR__, "common.jl"))

gaussfits_fullband = deserialize(joinpath(INTERMEDIATE_DIR, "gaussfits_fullband.jls"))

let
	fig = Figure()
	ax = Axis(fig[1,1]; width=300, height=300, aspect=DataAspect(), autolimitaspect=1, xlabel="Relative RA (mas⋅GHz²)", ylabel="Relative Dec (mas⋅GHz²)", to_xy_attrs(ticks=WilkinsonTicks(10, k_min=7, k_max=20))..., xreversed=true)

	scatter!((0,0), markersize=0, label=rich("Observed size\n(FWHM) at:", font=:bold))

	@p let
		gaussfits_fullband
		map() do (;freq, sky)
			sky = @modify(MCM.pmedian, sky |> RecursiveOfType(MCM.Particles))
			core_νnorm = @p let
				sky
				@modify(θ -> θ * u"GHz"(freq)^2, __.σ_major)
				@modify(θ -> zero(θ) * u"GHz"(freq)^2, __.coords[∗])  # just for proper units
			end

			color = get(ColorSchemes.turbo, shift_range(log(ustrip(u"GHz", freq)), (log(1.3)..log(5)) => 0..0.8))
			poly!(core_νnorm; color=(color, 0.05), strokecolor=color, strokewidth=2, label=f"{freq|>u\"GHz\":.1f}")
		end
	end

	scatter!((0,0), markersize=0, label="")
	lines!(x -> cot(gallat_PA)*x, linestyle=:dash, color=:gray, label="Line of constant\nGalactic latitude", yautolimits=false)

	Legend(fig[1,2], ax, merge=true)
	save(joinpath(FIGS_DIR, "broadening_size_vis.pdf"), fig; backend=CairoMakie)
	@info "Saved figs/broadening_size_vis.pdf"
end

include(joinpath(@__DIR__, "common.jl"))

@info "Loading fullband fits..."
gaussfits_fullband = deserialize(joinpath(INTERMEDIATE_DIR, "gaussfits_fullband.jls"))

@info "Loading raw data..."
alldata = load_all_data()

@assert length(alldata) == length(gaussfits_fullband)

result = map(alldata, gaussfits_fullband) do r, f
	@assert visfile_band(r.f) == f.band
	sky = @modify(MCM.pmedian, $(f.sky) |> RecursiveOfType(MCM.Particles))
	uvtbl_r = @p r.uvtbl mapinsert(value_mod=U.Value(visibility(abs, sky, UV(_))))
	mod_vis = map(row -> U.value(row.value_mod), uvtbl_r)
	gains = VLBI.solve_gains(VLBI.AmplitudeAnalytic(), mod_vis, uvtbl_r)
	gains = @modify(g -> clamp(g, 1..4), gains.gain[∗])
	scal_lavg = 2 * median(log.(gains.gain))
	@info "  $(f.band) $(visfile_date(r.f)): model_flux = $(exp(scal_lavg))"
	f = @insert f.date = visfile_date(r.f)
	@insert f.model_flux = exp(scal_lavg)
end

mkpath(INTERMEDIATE_DIR)
serialize(joinpath(INTERMEDIATE_DIR, "gaussfits_fullband_post.jls"), result)
@info "Saved $(length(result)) entries to gaussfits_fullband_post.jls"

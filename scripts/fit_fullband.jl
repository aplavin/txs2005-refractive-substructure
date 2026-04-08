include(joinpath(@__DIR__, "common.jl"))

@info "Loading all data..."
alldata = load_all_data()
@info "Loaded $(length(alldata)) visibility files"

gaussfits_fullband = @p let
	alldata
	map() do (; f, uvtbl, freq, maxuv)
		@info "Fitting $(visfile_band(f)) band ($(basename(f)))" freq=u"GHz"(freq) maxuv
		uvtbl = @p uvtbl filter(U.nσ(_.value) > 10)
		uvtbl_gauss = @p uvtbl filter(norm(UV(_)) < maxuv)
		uvtbl_fit = @p uvtbl_gauss VLBI.closures_all(VLBI.ClosureAmpSpec; rescale_redundancy=true)
		@info "  $(length(uvtbl_fit)) closure amplitudes"
		am = amodel(uvtbl_fit; flux_int=1, size_int=(0.2u"mas"..70u"mas"))
		pt = pigeons(;target=am, multithreaded=true, record=[traces; round_trip; record_default()], seed=rand(Int), n_rounds=7, n_chains=15)
		(;sky) = samples(MCM.Particles, pt)
		@info "  Done"
		(; band=visfile_band(f), freq, maxuv, sky)
	end
end

mkpath(INTERMEDIATE_DIR)
serialize(joinpath(INTERMEDIATE_DIR, "gaussfits_fullband.jls"), gaussfits_fullband)
@info "Saved $(length(gaussfits_fullband)) fullband fits to $INTERMEDIATE_DIR/gaussfits_fullband.jls"

include(joinpath(@__DIR__, "common.jl"))

@info "Loading all data..."
alldata = load_all_data()
@info "Loaded $(length(alldata)) visibility files"

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
@info "Split into $(length(perif_data)) per-IF entries"

gaussfits_perif = @p let
	perif_data
	map() do r
		@info "Fitting $(visfile_band(r.f)) band IF at $(u"GHz"(r.freq)) ($(basename(r.f)))"
		(; uvtbl, maxuv) = r
		uvtbl = @p uvtbl filter(U.nσ(_.value) > 10)
		uvtbl_gauss = @p uvtbl filter(norm(UV(_)) < maxuv) 
		uvtbl_fit = @p uvtbl_gauss VLBI.closures_all(VLBI.ClosureAmpSpec; rescale_redundancy=true)
		@info "  $(length(uvtbl_fit)) closure amplitudes"
		am = amodel(uvtbl_fit; flux_int=1, size_int=(0.2u"mas"..70u"mas"))
		pt = pigeons(;target=am, multithreaded=true, record=[traces; round_trip; record_default()], seed=rand(Int), n_rounds=8, n_chains=15)
		(;sky) = samples(MCM.Particles, pt)
		@info "  Done"
		(; band=visfile_band(r.f), date=visfile_date(r.f), freq=r.freq, freq_band=r.freq_band, sky)
	end
end

mkpath(INTERMEDIATE_DIR)
serialize(joinpath(INTERMEDIATE_DIR, "gaussfits_perif.jls"), gaussfits_perif)
@info "Saved $(length(gaussfits_perif)) per-IF fits to $INTERMEDIATE_DIR/gaussfits_perif.jls"

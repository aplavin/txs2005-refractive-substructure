include(joinpath(@__DIR__, "common.jl"))

gaussfits_perif = deserialize(joinpath(INTERMEDIATE_DIR, "gaussfits_perif.jls"))
data = @p gaussfits_perif map(@modify(U.Value, _.sky |> RecursiveOfType(MCM.Particles))) filter(_.freq < 4u"GHz")

@info "Fitting power-law for major axis ($(length(data)) points)..."
fit_major = let
	params=(θmaj₁=LogUniform(10, 1000), k=Uniform(1, 3), fracerr=LogUniform(1e-4, 1e-1))
	func=(x, p) -> ((@o ustrip(u"mas", fwhm_max(_.sky))) => p.θmaj₁ * ustrip(u"GHz", x.freq)^(-p.k),)
	am = accmodel(func, data, params; distf=(modval, unc, p) -> Normal(modval, unc + p.fracerr * modval))
	pt = pigeons(;target=am, multithreaded=true, record=[traces; round_trip; record_default()], seed=rand(Int), n_rounds=6, n_chains=15)
	@info "  Done"
	(; fit=samples(MCM.Particles, pt), func)
end

@info "Fitting power-law for minor axis ($(length(data)) points)..."
fit_minor = let
	params=(θmin₁=LogUniform(10, 1000), k=Uniform(1, 3), fracerr=LogUniform(1e-4, 1e-1))
	func=(x, p) -> ((@o ustrip(u"mas", fwhm_min(_.sky))) => p.θmin₁ * ustrip(u"GHz", x.freq)^(-p.k),)
	am = accmodel(func, data, params; distf=(modval, unc, p) -> Normal(modval, unc + p.fracerr * modval))
	pt = pigeons(;target=am, multithreaded=true, record=[traces; round_trip; record_default()], seed=rand(Int), n_rounds=6, n_chains=15)
	@info "  Done"
	(; fit=samples(MCM.Particles, pt), func)
end

mkpath(INTERMEDIATE_DIR)
serialize(joinpath(INTERMEDIATE_DIR, "powerlaw_fits.jls"), (; major=fit_major, minor=fit_minor))
@info "Saved power-law fits to $INTERMEDIATE_DIR/powerlaw_fits.jls"

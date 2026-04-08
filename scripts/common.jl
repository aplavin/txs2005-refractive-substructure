using VLBIFiles
using InterferometricModels
using Unitful, UnitfulAngles, UnitfulAstro
using DataManipulation
using Statistics
using Glob
using Uncertain
using AccessorsExtra
using LinearAlgebra: norm
using StaticArrays
using Dates
using IntervalSets
using StructArrays
using Distributions
using Pigeons
using AccessibleModels
import MonteCarloMeasurements as MCM
using Serialization

import CairoMakie
using MakieExtra
using ColorSchemes
using PyFormattedStrings
using DateFormats
using VLBIPlots

const PATHS = let
	d = Dict{String,String}()
	for line in eachline(joinpath(@__DIR__, "..", "paths.mk"))
		m = match(r"^(\w+)\s*=\s*(.+)$", line)
		m !== nothing && (d[m[1]] = normpath(joinpath(@__DIR__, "..", strip(m[2]))))
	end
	d
end
const DATA_DIR = PATHS["DATA_DIR"]
const FIGS_DIR = PATHS["FIGS_DIR"]
const INTERMEDIATE_DIR = PATHS["INTERMEDIATE_DIR"]
const gallat_PA = 32.4u"°"
intrinsic_source_fwhm(ν) = 4.2u"mas" / (ν/1u"GHz")  # from Koryukova+23; exact value doesn't matter, source is scattering-dominated

visfile_band(path) = String(split(first(splitext(basename(path))), "_")[2])
visfile_date(path) = let parts = split(first(splitext(basename(path))), "_")
	Date(parse(Int, parts[3]), parse(Int, parts[4]), parse(Int, parts[5]))
end

function get_visfiles()
	[
		readdir(glob"rfc/*_vis.fits", DATA_DIR);
		@p readdir(glob"BG246T/*/*.uvf", DATA_DIR) filter(!endswith(_, "_crc.uvf"));
	]
end

function load_uvtable(f)
	@p VLBI.load(f) uvtable filter(_.stokes == :I || VLBIData.is_parallel_hands(_.stokes)) VLBI.rescale_visibility_errors(VLBI.ConsecutiveDifferencesStandard()) VLBI.average_data(VLBI.GapBasedScans())
end

maxuv_forfit(uvtbl) = @p uvtbl sort(by=norm(UV(_))) filterfirst(U.nσ(_.value) < 5) norm(UV(__))

function load_all_data()
	visfiles = get_visfiles()
	@p let
		map(visfiles) do f
			uvtbl = load_uvtable(f)
			maxuv = maxuv_forfit(uvtbl)
			freqs = @p uvtbl extrema(frequency) Interval(__...)
			(; f, uvtbl, maxuv, freq=mean(freqs))
		end
		sort(by=_.freq)
	end
end

# --- Model / likelihood functions ---

function loglike(model, data::AbstractVector)
	vismodel = visibility(model.sky)
	sum(x -> VLBI.loglike(vismodel, x), data)
end

function amodel(uvtbl; flux_int, size_int, initsky=EllipticGaussian(flux=mean(flux_int), σ_major=0.5u"mas", ratio_minor_major=0.9, pa_major=0., coords=(0.,0.) .* u"mas"))
	AccessibleModel(Base.Fix2(loglike, uvtbl), (;sortby=norm∘coords, sky=initsky), (
		(flux_int isa Number ? () : ((@o flux(_.sky)) => LogUniform(endpoints(flux_int)...),))...,
		(@o ustrip(u"mas", fwhm_max(_.sky))) => LogUniform(ustrip.(endpoints(size_int))...),
		(@o _ |> RecursiveOfType(EllipticGaussian; recurse=Union{Tuple,NamedTuple,MultiComponentModel}) |> _.ratio_minor_major) => Uniform(0, 1),
		(@o _ |> RecursiveOfType(EllipticGaussian; recurse=Union{Tuple,NamedTuple,MultiComponentModel}) |> _.pa_major) => Uniform(-π/2, π/2),
	))
end


function accmodel(func, data, params; distf=(modval, unc, p) -> Normal(modval, unc))
	oparams = AccessorsExtra.flat_concatoptic(params, RecursiveOfType(Distribution))
	AccessibleModel(modify(median, params, oparams), AccessorsExtra._optics(oparams) .=> getall(params, oparams)) do p
		sum(data) do x
			modvals = func(x, p)
			sum(modvals) do (o, modval)
				obs = o(x)::U.Value
				logpdf(distf(modval, U.uncertainty(obs), p), U.value(obs))
			end
		end
	end
end

# RMS refractive substructure noise in visibility: Eq. 19 of Johnson & Gwinn (2015), for Kolmogorov (α=5/3)
function scat_vis_noise_mul(; uvdist, λ, θ_scat, θ_img, D)
	B = uvdist * λ
	return 0.0038 * NoUnits(λ/6u"cm") * NoUnits(B/(10^5)u"km")^(-5/6) * NoUnits(θ_scat/30u"μas")^(5/6) * NoUnits(θ_img/300u"μas")^(-2) * NoUnits(D/1u"kpc")^(-1/6)
end

function scat_vis_noise(srcmod; uvdist, ν, θ_scat, D, α, r_in)
	@assert isapprox(α, 5/3; rtol=1e-3)
	@assert r_in ≈ 1000u"km"
	λ = u"c"/ν |> u"m"
	compnoises = map(InterferometricModels.components(srcmod)) do comp
		θ_img = hypot(fwhm_min(comp), θ_scat)
		flux(comp) * scat_vis_noise_mul(; uvdist, λ, θ_scat, θ_img, D)
	end
	return hypot(compnoises...)
end

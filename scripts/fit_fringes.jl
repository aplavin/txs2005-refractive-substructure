include(joinpath(@__DIR__, "common.jl"))

import FFTW, FITSIO
using SkyCoords
using AxisKeysExtra
using RectiGrids
using Unitful: °

# --- Helper functions ---

function extract_raw_data_for_source(uvf::VLBI.UVData, source::Function)
	fitsfile = FITSIO.FITS(uvf.path)
	source_id = @p StructArrays.fromtable(fitsfile["SOURCE"]) filteronly(source(_)) @oget __.var"ID_NO." __.SOURCE_ID
	uvdata_src_raw = @p VLBI.read_data_raw(uvf) filterview(@o _.SOURCE == source_id)
end

UV_from_uvrow(r) =
	UV(
		(@oget r.var"UU---SIN" r.var"UU--SIN" r.var"UU-L"),
		(@oget r.var"VV---SIN" r.var"VV--SIN" r.var"VV-L")
	) .* u"c*s" .|> u"m"

function calculate_timesteps(times::AbstractVector)
	ts = (times .- first(times)) .|> u"s" .|> float
	dt = @p ts diff map(abs) median
	tns = @p ts enumerate() map() do (i, t)
		n, Δ = divrem(t, dt, RoundNearest)
		@assert abs(Δ) < 0.3*dt (dt, t/dt)
		Int(n)
	end
	@assert issorted(tns)
	return (;dt, tns)
end

function uvtable_to_visarray(uvtable::StructVector{<:NamedTuple})
	(;dt, tns) = calculate_timesteps(uvtable.datetime)
	ns = @p tns extrema() range(__...)
	Z = zero(uvtable[1].FLUX)
	alldata = @p let
		ns
		map() do n
			ix = searchsorted(tns, n)
			isempty(ix) && return Z
			return uvtable.FLUX[only(ix)]
		end
		stack(KeyedArray(__, time=ns .* dt))
	end
end

function is_adjacent(a::VLBI.FrequencyWindow, b::VLBI.FrequencyWindow)
	fa = VLBIFiles.frequencies(a)
	fb = VLBIFiles.frequencies(b)
	return a < b && step(fa) ≈ step(fb) && first(fb) - last(fa) < 1.5*step(fa)
end

function adjacent_groups(fws)
	groups = [
		[fws[1]],
	]
	for fw in fws[2:end]
		if is_adjacent(last(groups)[end], fw)
			push!(last(groups), fw)
		else
			push!(groups, [fw])
		end
	end
	return groups
end

function zeropad(A::AbstractArray; factor)
	P = zeros(eltype(A), size(A) .* factor)
	P[CartesianIndices(A)] .= A
	return P
end
function zeropad(A::KeyedArray; factor)
	KeyedArray(zeropad(AxisKeys.keyless_unname(A); factor); map(ak -> expand_range(ak; factor), named_axiskeys(A))...)
end

expand_range(rng::StepRangeLen; factor::Int) = @set rng.len *= factor
expand_range(rng::LinRange; factor::Int) = LinRange(first(rng), first(rng) + (last(rng) - first(rng)) * factor, length(rng) * factor)


function vec_to_range(x::AbstractVector)
	Δs = diff(x)
	if !isapprox(maximum(Δs), minimum(Δs); rtol=1e-2)
		error("Different time steps: $(minimum(Δs)) .. $(maximum(Δs))")
	end
	range(start=first(x), stop=last(x), length=length(x))
end

const freqbins = [
	(1u"GHz"..2u"GHz"),
	(2u"GHz"..3u"GHz"),
	(4u"GHz"..6u"GHz"),
	(8u"GHz"..9u"GHz"),
]

freqbinix(freq::Quantity) = findonly(gb -> freq ∈ gb, freqbins)
freqbinix(freq::AbstractVector{<:VLBIFiles.FrequencyWindow}) = uniqueonly(map(freqbinix, freq))
freqbinix(freq::VLBIFiles.FrequencyWindow) = freqbinix(frequency(freq))

# --- Data loading ---

@info "Collecting FITS files..."
fitsfiles = vcat(readdir(glob"*/*.*fits", joinpath(DATA_DIR, "archival", "VLBA")), readdir(glob"*/*/*.*fits", joinpath(DATA_DIR, "archival", "VLBA")))
@info "Found $(length(fitsfiles)) FITS files"

tgtcoords = ICRSCoords(301.94°, 40.50°)

@info "Loading raw data..."
alldata1 = @p fitsfiles filtermap() do f
	@info "  Loading $(basename(dirname(f)))/$(basename(f))"
	uvf = VLBI.load(VLBI.UVData, f)
	source_pred = s -> SkyCoords.separation(ICRSCoords(s.RAEPO*°, s.DECEPO*°), tgtcoords) < 0.01°
	(;file=f, uvf, raw_data=extract_raw_data_for_source(uvf, source_pred))
end
@info "Loaded $(length(alldata1)) observations"

@info "Processing visibility data..."
alldata2 = map(alldata1) do r
	(;uvf, raw_data) = r
	@insert r.uvdata = @p let
		raw_data
		mapinsert⁻(datetime=@o VLBIFiles.DateTime_from_DATE_TIME(_.DATE, _.TIME))
		mapinsert(spec=@o VisSpec(
			VLBIFiles.Baseline_from_fits(_.BASELINE, uvf.ant_arrays),
			UV_from_uvrow(_)))
		mapset(FLUX=r -> @p let
			r.FLUX
			__[RA=1, DEC=1]
			complex.(__(COMPLEX=:re), __(COMPLEX=:im))
			@set axiskeys($__, :BAND) = @p uvf.freq_windows filter(_.freqid == r.FREQID)
		end)
		delete(__, @o _.var"UU---SIN" _.var"VV---SIN" _.var"WW---SIN" _.var"UU--SIN" _.var"VV--SIN" _.var"WW--SIN" _.var"UU-L" _.var"VV-L" _.var"WW-L" _.BASELINE _.FREQID)
		filter(!allequal(antenna_names(_)))
		VLBI.add_scan_ids(VLBI.GapBasedScans(30u"s"))
	end
end

# --- Fringe fitting ---

@info "Running fringe fitting..."
all_fringefits = @p alldata2 enumerate() map() do (i, r)
	@info "  Fringe fitting observation $i/$(length(alldata2)): $(basename(dirname(r.file)))/$(basename(r.file))"
	(;uvf, uvdata) = r
	fringefits = @p uvdata filter(!allequal(antenna_names(_))) group_vg((;named_axiskeys(_.FLUX).BAND, _.scan_id, ants=antenna_names(_))) collect flatmap() do scan
		alldata = uvtable_to_visarray(value(scan))
		bandsets = vcat(
			axiskeys(alldata, :BAND),
			@p adjacent_groups(axiskeys(alldata, :BAND)) filter(length(_) > 1)
		)

		@p grid(band=bandsets, stokes=@p axiskeys(alldata, :STOKES) filter(VLBIData.is_parallel_hands(_))) vec filtermap() do (;band, stokes)
			freqs = band isa VLBI.FrequencyWindow ? VLBIFiles.frequencies(band) : @p band flatmap(VLBIFiles.frequencies) vec_to_range()
			data = @p alldata(;BAND=band, STOKES=stokes) |>
				reshape(__, (:, length(axiskeys(__, :time)))) |>
				KeyedArray(__, FREQ=freqs, time=axiskeys(alldata, :time))

			pad_factor = 2
			fabs = @p let
				FFTW.fft(zeropad(data; factor=pad_factor))
				FFTW.fftshift
				@set dimnames(__) = (:delay, :rate)
				@modify(d -> d .|> u"ns", __ |> axiskeys(_, :delay))
				@modify(d -> d .|> u"mHz", __ |> axiskeys(_, :rate))
				map(abs)
			end

			med = median(fabs)
			σ = med / median(Rayleigh(1))
			noisedist = Rayleigh(iszero(σ) ? √(eps(σ)) : σ)

			ntrials = length(fabs)  # these trials aren't independent, making downstream calculations a conservative estimate
			peakval, peakloc = with_axiskeys(findmax)(fabs)

			(; band, stokes, key(scan)..., uv=mean(UV, scan), peakval, peakloc, noisedist, ntrials)
		end
	end
	@insert r.fringefits = fringefits
end

all_fringefits2 = map(all_fringefits) do r
	r = @set r.fringefits = @p r.fringefits mapinsert(freqix=freqbinix(_.band))
	r = @modify(r -> (@insert r.spec=VisAmpSpec(VisSpec(Baseline(r.ants), NoUnits.(r.uv / (u"c"/mean(frequency(r.band, Interval))))))), r.fringefits[∗])
	r = @modify(r -> (@insert r.value=r.peakval), r.fringefits[∗])
end

# --- Serialize ---

results = map(all_fringefits2) do r
	(; fringefits=r.fringefits, uvf_path=r.uvf.path, date_obs=r.uvf.header.date_obs)
end

mkpath(INTERMEDIATE_DIR)
serialize(joinpath(INTERMEDIATE_DIR, "all_fringefits.jls"), results)
@info "Saved $(length(results)) fringe fit results to intermediate/all_fringefits.jls"

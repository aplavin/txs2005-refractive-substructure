include(joinpath(@__DIR__, "common.jl"))
using MakieBake
using Random
using RectiGrids
using AxisKeysExtra
import ScatteringOpticsExtra as SO
import VLBISkyModels as VSM

gaussfits_fullband = deserialize(joinpath(INTERMEDIATE_DIR, "gaussfits_fullband_post.jls"))

# --- Observables (for MakieBake animation) ---
p = Observable((freq=1.4u"GHz", α=5/3, r_in=1000u"km"))

fit = @lift argmin(r -> abs(r.freq - $p.freq), gaussfits_fullband)
sky = @lift modify(MCM.pmedian, $fit.sky, RecursiveOfType(MCM.Particles))
model_flux = @lift $fit.model_flux

sourcemodel = @lift construct(CircularGaussian,
	flux => $model_flux,
	fwhm_max => intrinsic_source_fwhm($p.freq),
	coords => (0, 0))

broadening = @lift construct(EllipticGaussian,
	flux => 1,
	fwhm_max => fwhm_max($sky),
	fwhm_min => fwhm_min($sky),
	position_angle => position_angle($sky),
	coords => (0, 0))

convmodel = @lift InterferometricModels.convolve($sourcemodel, $broadening)

# Image grid: FOV scales with frequency
npix = 512
fov = @lift 50u"mas" / ($p.freq / 1.4u"GHz")^2
imggrid = @lift grid(SVector,
	X=LinRange(0u"mas" ± $fov, npix),
	Y=LinRange(0u"mas" ± $fov, npix))

# --- Column 1: Intrinsic source image (analytic) ---
sourceimg = @lift intensity(sourcemodel[]).($imggrid) .|> ustrip

# --- Column 2: Broadened image (analytic convolution) ---
convimg = @lift intensity(convmodel[]).($imggrid) .|> ustrip

# --- Column 3: Scattered image (depends on α, r_in too) ---
scatimg_raw = @lift begin
	(; freq, α, r_in) = p[]
	λ = u"c" / freq |> u"cm"
	sm = SO.ScatteringModel(;
		target_distance=1u"Gpc", screen_distance=1.4u"kpc",
		r_in, α, broadening=broadening[], λ₀=λ,
		screen_spec=SO.StochasticPhaseScreen()
	)
	Random.seed!(Random.default_rng(), round(Int, 123 * ustrip(freq)))
	pixarea = abs(prod(step, axiskeys($imggrid)))
	src = intensity(sourcemodel[]).($imggrid) * pixarea .|> ustrip
	SO.scatter_image(sm, src; λ)
end

# For display: convert axis keys to mas
scatimg_display = @lift @modify(ak -> u"mas".(ak), axiskeys($scatimg_raw)[∗])

# --- Visibility: 1D slice along u-axis (v=0) ---
maxbl_km = 8000
nuvpts = 1000

baselines_km = range(0, maxbl_km, length=nuvpts) |> collect

uvs = @lift map(range(0, maxbl_km, length=nuvpts)) do bl
	λ = u"c" / $p.freq |> u"cm"
	bl * u"km" / λ |> NoUnits
end

# Column 1 visibility (analytic)
vis1 = @lift map($uvs) do uv
	abs(visibility($sourcemodel, UV(uv, 0)))
end

# Column 2 visibility (analytic)
vis2 = @lift map($uvs) do uv
	abs(visibility($convmodel, UV(uv, 0)))
end

# Column 3 visibility (from scattered image, using raw axis keys for NFFT)
vis3 = @lift begin
	model = VSM.ContinuousImage(VSM.IntensityMap(scatimg_raw[]), VSM.DeltaPulse())
	uvspecs = map($uvs) do uv
		(; spec=UV(uv, 0), freq_spec=p[].freq)
	end
	vmap = VSM.visibilitymap(model, uvspecs)
	map(r -> abs(r.value), vmap)
end

# --- Figure ---
fig = Figure(figure_padding=100)

# Row 1: Images
titles = [
	"Intrinsic blazar core\n" * rich("(constrained from multifrequency)", color=Makie.Colors.Gray(0.6)),
	"+ Scatter Broadening\n" * rich("(directly measured observed size)", color=Makie.Colors.Gray(0.6)),
	"+ Scatter Substructure\n" * rich("(detected in observations)", color=Makie.Colors.Gray(0.6))
]
img_obss = [sourceimg, convimg, scatimg_display]
img_axes = map(enumerate(zip(titles, img_obss))) do (i, (title, img_obs))
	ax = Axis(fig[1, i]; aspect=DataAspect(), title,
		width=250, height=250,
		xticks=[0], yticks=[0]  # just to keep the axis size constant, no matter the limits
	)
	image!(ax, @lift(reverse($img_obs; dims=1)); colormap=:inferno)
	scalebar!(1u"mas", color=:white)
	hidedecorations!()
	ax
end

# Broadening ellipse on columns 2, 3
for i in 2:3
	pl = poly!(img_axes[i], @lift(ustrip(InterferometricModels.convolve($sourcemodel, $broadening)));
		color=:transparent, strokecolor=:white, strokewidth=2)
	translate!(pl, 0, 0, 100)
end

# Color borders on image axes
vis_colors = [:royalblue, :seagreen, :orangered]
for (i, ax) in enumerate(img_axes)
	for s in [:leftspinecolor, :rightspinecolor, :topspinecolor, :bottomspinecolor]
		getproperty(ax, s)[] = to_color(vis_colors[i])
	end
	ax.spinewidth = 3
end

# Row 2: Visibility amplitude
vis_obss = [vis1, vis2, vis3]
vis_axes = map(enumerate(vis_obss)) do (i, vis_obs)
	ax = Axis(fig[2, i]; yscale=AsinhScale(1e-2),
		limits=((0, maxbl_km), (0, 4)),
		xlabel="Baseline (km)", ylabel="Amplitude (Jy)",
		width=250, height=200)
	hidespines!()
	# Previous columns as dashed lines
	for j in 1:i-1
		lines!(ax, baselines_km, vis_obss[j]; color=vis_colors[j], linewidth=2, linestyle=:dash)
	end
	# Current column as solid line
	lines!(ax, baselines_km, vis_obs; color=vis_colors[i], linewidth=2)
	ax
end

# Link visibility y-axes
linkyaxes!(vis_axes...)

colgap!(fig.layout, 50)
rowgap!(fig.layout, 50)
resize_to_layout!()
save(joinpath(FIGS_DIR, "stages.png"), fig)

outdir = joinpath(FIGS_DIR, "stages_html")
bake_html(
	p => (
		(@o ustrip(u"GHz", _.freq)) => [1.4, 1.8, 2.3, 5],
		(@o _.α) => [1.5, 1.6, 1.6666666, 1.7, 1.8],
		(@o ustrip(u"km", _.r_in)) => [200, 500, 1000, 2000, 5000, 10000],
	);
	blocks=vcat(img_axes, vis_axes),
	outdir,
)
write(joinpath(outdir, "layout.js"), """
const HEADER = 'Scattering in the ISM:<br>from <span style="color:royalblue">intrinsic structure</span> to <span style="color:orangered">final image</span>';
const TITLE = 'Scattered image formation';
const MAXWIDTH = '1500px';
const LAYOUT = ["A B C S", "D E F S"];
""")

@info "Saved"

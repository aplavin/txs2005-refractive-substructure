include(joinpath(@__DIR__, "common.jl"))

let
	fplt = FPlot([
		(date=Date(2010, 11, 5), freq=1.4, exp="BG196"),
		(date=Date(2010, 11, 5), freq=1.6, exp="BG196"),
		(date=Date(2013, 7, 4), freq=1.6, exp="BS224"),

		(date=Date(2017, 7, 12), freq=1.4, exp=""),
		(date=Date(2017, 7, 12), freq=1.8, exp=""),
		(date=Date(2017, 7, 12), freq=2.3, exp=""),
		(date=Date(2017, 7, 12), freq=5.0, exp=""),
		(date=Date(2017, 8, 1), freq=1.4, exp=""),
		(date=Date(2017, 8, 1), freq=1.8, exp=""),
		(date=Date(2017, 8, 1), freq=2.3, exp=""),
		(date=Date(2017, 8, 1), freq=5.0, exp=""),
		(date=Date(2017, 8, 15), freq=1.4, exp="BG246"),
		(date=Date(2017, 8, 15), freq=1.8, exp="BG246"),
		(date=Date(2017, 8, 15), freq=2.3, exp="BG246"),
		(date=Date(2017, 8, 15), freq=5.0, exp="BG246"),
		(date=Date(2017, 9, 9), freq=1.4, exp=""),
		(date=Date(2017, 9, 9), freq=1.8, exp=""),
		(date=Date(2017, 9, 9), freq=2.3, exp=""),
		(date=Date(2017, 9, 9), freq=5.0, exp=""),
		(date=Date(2017, 9, 10), freq=1.4, exp=""),
		(date=Date(2017, 9, 10), freq=1.8, exp=""),
		(date=Date(2017, 9, 10), freq=2.3, exp=""),
		(date=Date(2017, 9, 10), freq=5.0, exp=""),
		(date=Date(2017, 9, 20), freq=1.4, exp=""),
		(date=Date(2017, 9, 20), freq=1.8, exp=""),
		(date=Date(2017, 9, 20), freq=2.3, exp=""),
		(date=Date(2017, 9, 20), freq=5.0, exp=""),
		(date=Date(2018, 1, 8), freq=1.4, exp="BG246"),
		(date=Date(2018, 1, 8), freq=1.8, exp="BG246"),
		(date=Date(2018, 1, 8), freq=2.3, exp="BG246"),
		(date=Date(2018, 1, 8), freq=5.0, exp="BG246"),

		(date=Date(2018, 8, 11), freq=2.3, exp="UG002"),
		(date=Date(2019, 2, 20), freq=1.4, exp="BG258"),
		(date=Date(2019, 2, 20), freq=1.6, exp="BG258"),
		(date=Date(2019, 2, 20), freq=2.3, exp="BG258"),
		(date=Date(2019, 2, 20), freq=5., exp="BG258"),
	],
		  AxFunc(label="Date", ticks=1:5000, @o yeardecimal(_.date)),
		  AxFunc(label="Frequency (GHz)", scale=log10, @o _.freq),
			marker='∘', markersize=30,
			color=Makie.ColorSchemes.seaborn_bright[3])
	fig,ax,plt = axplot(scatter)(fplt; axis=(;limits=(nothing, (1,14)), yticks=0:5), figure=(;size=(400, 220)))
	texts = map(@p fplt.data group_vg((;_.exp, _.date))) do gr
		plt = vlines!(yeardecimal(key(gr).date); color=(:black, 0.5), linewidth=1, linestyle=:dash); translate!(plt, 0,0,-100)
		textwithbox!((@lift (yeardecimal(key(gr).date), $(MakieExtra.rel2data(:y, 1)))), text=key(gr).exp, rotation=90u"°", align=(:right, :center), offset=(0, -2), xautolimits=false, yautolimits=false, poly=(;color=(:white, 1)))
	end
	save(joinpath(FIGS_DIR, "detections_summary.pdf"), fig; backend=CairoMakie)
	@info "Saved figs/detections_summary.pdf"
end

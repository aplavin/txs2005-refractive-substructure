JULIA = julia --project=scripts

include paths.mk

COMMON_DEPS = scripts/.pkg.stamp scripts/common.jl
UVFITS_FILES = $(wildcard $(DATA_DIR)/rfc/*.fits) $(wildcard $(DATA_DIR)/BG246T/*/*.uvf)
FITSIDI_FILES = $(wildcard $(DATA_DIR)/archival/VLBA/*/*.fits) $(wildcard $(DATA_DIR)/archival/VLBA/*/*/*.fits)

.PHONY: all crop clean

all: $(FIGS_DIR)/broadening_size_vis.pdf $(FIGS_DIR)/gaussfit_summary.pdf $(FIGS_DIR)/detections_summary.pdf $(FIGS_DIR)/radplot_gauss_C_2017-07-12.pdf $(FIGS_DIR)/simimage.pdf $(FIGS_DIR)/fringeplots_summary.pdf $(FIGS_DIR)/scatmodel_uvamp.pdf $(FIGS_DIR)/stages.png

crop:
	for f in **/*.pdf; do uvx pdfcropmargins -p 0 -o "$$f.tmp" "$$f" && mv "$$f.tmp" "$$f"; done

# --- Fitting (expensive) ---

$(INTERMEDIATE_DIR)/gaussfits_fullband.jls: scripts/fit_fullband.jl $(COMMON_DEPS) $(UVFITS_FILES)
	@mkdir -p $(INTERMEDIATE_DIR)
	$(JULIA) scripts/fit_fullband.jl

$(INTERMEDIATE_DIR)/gaussfits_fullband_post.jls: scripts/postprocess_fullband.jl $(COMMON_DEPS) $(INTERMEDIATE_DIR)/gaussfits_fullband.jls $(UVFITS_FILES)
	$(JULIA) scripts/postprocess_fullband.jl

$(INTERMEDIATE_DIR)/gaussfits_perif.jls: scripts/fit_perif.jl $(COMMON_DEPS) $(UVFITS_FILES)
	@mkdir -p $(INTERMEDIATE_DIR)
	$(JULIA) scripts/fit_perif.jl

$(INTERMEDIATE_DIR)/powerlaw_fits.jls: scripts/fit_powerlaw.jl $(COMMON_DEPS) $(INTERMEDIATE_DIR)/gaussfits_perif.jls
	$(JULIA) scripts/fit_powerlaw.jl

$(INTERMEDIATE_DIR)/all_fringefits.jls: scripts/fit_fringes.jl $(COMMON_DEPS) $(FITSIDI_FILES)
	@mkdir -p $(INTERMEDIATE_DIR)
	$(JULIA) scripts/fit_fringes.jl

# --- Plotting (fast) ---

$(FIGS_DIR)/broadening_size_vis.pdf: scripts/plot_broadening_size_vis.jl $(COMMON_DEPS) $(INTERMEDIATE_DIR)/gaussfits_fullband.jls
	@mkdir -p $(FIGS_DIR)
	$(JULIA) scripts/plot_broadening_size_vis.jl

$(FIGS_DIR)/gaussfit_summary.pdf: scripts/plot_gaussfit_summary.jl $(COMMON_DEPS) $(INTERMEDIATE_DIR)/gaussfits_perif.jls $(INTERMEDIATE_DIR)/powerlaw_fits.jls
	@mkdir -p $(FIGS_DIR)
	$(JULIA) scripts/plot_gaussfit_summary.jl

$(FIGS_DIR)/detections_summary.pdf: scripts/plot_detections_summary.jl $(COMMON_DEPS)
	@mkdir -p $(FIGS_DIR)
	$(JULIA) scripts/plot_detections_summary.jl

$(FIGS_DIR)/radplot_gauss_C_2017-07-12.pdf: scripts/plot_radplot_gauss.jl $(COMMON_DEPS) $(INTERMEDIATE_DIR)/gaussfits_perif.jls $(INTERMEDIATE_DIR)/gaussfits_fullband_post.jls $(UVFITS_FILES)
	@mkdir -p $(FIGS_DIR)
	$(JULIA) scripts/plot_radplot_gauss.jl

$(FIGS_DIR)/simimage.pdf: scripts/plot_simimage.jl $(COMMON_DEPS) $(INTERMEDIATE_DIR)/gaussfits_fullband.jls
	@mkdir -p $(FIGS_DIR)
	$(JULIA) scripts/plot_simimage.jl

$(FIGS_DIR)/fringeplots_summary.pdf: scripts/plot_fringeplots_summary.jl $(COMMON_DEPS) $(INTERMEDIATE_DIR)/all_fringefits.jls $(INTERMEDIATE_DIR)/gaussfits_fullband.jls
	@mkdir -p $(FIGS_DIR)
	$(JULIA) scripts/plot_fringeplots_summary.jl

$(FIGS_DIR)/scatmodel_uvamp.pdf: scripts/plot_scatmodel_uvamp.jl $(COMMON_DEPS) $(INTERMEDIATE_DIR)/gaussfits_fullband_post.jls
	@mkdir -p $(FIGS_DIR)
	$(JULIA) scripts/plot_scatmodel_uvamp.jl

$(FIGS_DIR)/stages.png: scripts/plot_stages.jl $(COMMON_DEPS) $(INTERMEDIATE_DIR)/gaussfits_fullband_post.jls
	@mkdir -p $(FIGS_DIR)
	$(JULIA) scripts/plot_stages.jl

scripts/.pkg.stamp: scripts/Project.toml scripts/Manifest.toml
	$(JULIA) -e 'using Pkg; Pkg.instantiate()'
	@touch $@

clean:
	rm -rf $(INTERMEDIATE_DIR)
	rm -rf $(FIGS_DIR)
	rm -f scripts/.pkg.stamp

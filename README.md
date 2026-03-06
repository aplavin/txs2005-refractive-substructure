# Reproducible code for: Refractive Substructure in TXS 2005+403

Paper "**Direct VLBI Detection of Interstellar Turbulence Imprint on a Quasar: TXS 2005+403**", [arXiv:2602.24255](https://arxiv.org/abs/2602.24255).

Code here performs Bayesian Gaussian model fitting to visibility data, characterizes frequency-dependent angular broadening, and conducts fringe fitting on archival VLBA data.
Generates all publication figures: see [./figs/](./figs/).

## Interactive figure

Explore interstellar scattering stages — intrinsic source, scatter broadening, and refractive substructure — across frequencies and turbulence parameters:
[**Scattering Stages**](https://aplavin.github.io/txs2005-refractive-substructure/figs/stages_html/)

## Prerequisites

Install:
- GNU Make
- Julia 1.10 ([downloads page](https://julialang.org/downloads/manual-downloads/))

And download this repo content to your computer.

## Data

Three data sources are required. All paths are relative to the parent `2005+403/` directory.

### Archival VLBA FITSIDI files

Download from the [NRAO Science Data Archive](https://data.nrao.edu/).
Project codes: BG196, BG246, BG258, BS224B, UG002.

Files go into `data/archival/VLBA/{project}/{subproject}/`.

### L and S band visibility files

UVFITS visibility files for J2007+4029 from [astrogeo.org](http://astrogeo.org/) catalog (L and S bands). Place into `data/rfc/`:
```
data/rfc/
├── J2007+4029_L_2010_11_05_pus_vis.fits
├── J2007+4029_S_1997_01_10_fey_vis.fits
├── J2007+4029_S_1997_01_10_pus_vis.fits
└── J2007+4029_S_2018_08_10_pet_vis.fits
```

### BG246T visibility files

Reach out to me to get these files, or alternatively retrieve BG246T experiment from NRAO archive and process yourself.
```
data/BG246T/
├── J2007+4029_C_2017_07_12/*.uvf
├── J2007+4029_L1_2017_07_12/*.uvf
├── J2007+4029_L2_2017_07_12/*.uvf
└── J2007+4029_S_2017_07_12/*.uvf
```

### Expected directory structure

```
/
├── data/
│   ├── rfc/
│   │   └── J2007+4029_{band}_{date}_{author}_vis.fits
│   ├── BG246T/
│   │   └── J2007+4029_{band}_{date}/*.uvf
│   └── archival/VLBA/
│       └── {project}/{subproject}/*.idifits
└── txs2005-refractive-substructure/  # this repo
    ├── scripts/        # Julia project, all scripts
    └── Makefile
```

## Running

Execute
```sh
make all
```
in the root directory of this repository.

The pipeline automatically installs required Julia package, runs computational scripts, and then plotting scripts.
The produce PDFs in `figs/`.
All data and intermediate dependencies are tracked by the Makefile; individual targets can be built separately, e.g. `make figs/gaussfit_summary.pdf`.

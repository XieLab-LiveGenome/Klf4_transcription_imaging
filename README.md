# Klf4_transcription_imaging

MATLAB pipelines for live-cell SR-SIM: 3D tracking of enhancer, promoter and MS2 loci, HMM burst calling, and detection, sizing, lifetime and colocalization of nuclear hubs.

Imaging platform: Zeiss Elyra 7 Lattice SIM (31.3 nm/pixel lateral, up to 4 colours).
Input: Carl Zeiss Image (`.czi`), 6D (S, T, Z, C, Y, X). Output: MATLAB `.mat`.

Full documentation of every script, parameter, and output variable lives in
[`CODE_METADATA.md`](CODE_METADATA.md) (locus tracking) and
[`CODE_METADATA_hub_analysis.md`](CODE_METADATA_hub_analysis.md) (hub analysis).

---

## 1. Setup

Each pipeline script begins with a project directory setup block that adds its
sibling folders to the MATLAB path. Two of those folders need to be populated with contents that are third-party software that cannot be
redistributed here, and must be installed before anything will run.

### 1.1 Bio-Formats (required by every script)

CZI files are read through Bio-Formats. Download `bioformats_package.jar` from

  https://www.openmicroscopy.org/bio-formats/downloads/

and place it in `bioformats/`. The setup block
registers the JAR with `javaaddpath`.

Bio-Formats is distributed by the Open Microscopy Environment under the GPL and
is not redistributed here.

### 1.2 vbFRET (required for HMM burst calling)

vbFRET is required by `HMM_fit_fun.m` for two-state HMM fitting of MS2 traces.
It is not redistributed here. Download from

  https://vbfret.sourceforge.net/

and place the package contents in `HMM fitting/`. Please cite
Bronson et al., *Biophysical Journal* 97(12):3196–205 (2009).

The vbFRET authors ask that the package not be further distributed without
their prior permission, which is why it is absent from this repository.

**Keep the original folder structure.** The setup block uses `genpath`, so
subfolders are added recursively — there is no need to flatten the package.
After adding it, check for name collisions with the Statistics and Machine
Learning Toolbox, since vbFRET bundles Netlab code that shadows some built-ins:

```matlab
which -all kmeans
which -all pca
```

If a MATLAB function is shadowed, add and remove the vbFRET path around the HMM
call rather than leaving it on the path for the whole session.

`HMM_fit_fun.m` (the HMM fitting wrapper function) is included
in `HMM fitting/`.

### 1.3 Cellpose (required for anything that segments nuclei)

Nuclear segmentation calls Cellpose through MATLAB's Python interface.

> **Version matters.** These pipelines were developed and validated against
> Cellpose 3.x. The `models.Cellpose` class was removed in Cellpose 4
> (Cellpose-SAM), so v4 raises
> `Unrecognized method, property, or field 'Cellpose' for class 'py.module'`.
> Cellpose-SAM is also a different model with different channel and diameter
> handling, so it produces different masks — and every hub detection threshold
> in this repository is referenced to the mean intensity of the segmented
> nucleus. Use 3.x for results comparable to the published analyses.

```bash
conda env create -f environment.yml
conda activate klf4-imaging
```

Then point MATLAB at that interpreter:

```matlab
pyenv('Version', '/path/to/envs/klf4-imaging/bin/python')
```

Verify with:

```matlab
pyenv
char(py.importlib.metadata.version('cellpose'))
```

### 1.4 MATLAB requirements

R2019b or later (`matchpairs`), with the Image Processing, Computer Vision,
Optimization, and Statistics and Machine Learning toolboxes.

### 1.5 Other third-party functions

`functions/` contains helper code from several sources. `pkfnd` and `cntrd` are
the Crocker–Grier particle-tracking routines; `feretDiameters` and a few others
come from the MATLAB File Exchange. Original attributions are preserved in the
function headers.

---

## 2. Repository layout

Klf4_transcription_imaging/
├── README.md
├── LICENSE
├── .gitignore
├── environment.yml
├── CODE_METADATA.md                  # locus tracking: full documentation
├── CODE_METADATA_hub_analysis.md     # hub analysis: full documentation
│
├── Image Analysis/                   # all pipelines documented here
│   │
│   │   # --- locus tracking (§3.1) ---
│   ├── Enhancer-Promoter_tracking.m
│   ├── Enhancer_promoter_MS2_tracking.m
│   ├── Enhancer_Hub_MS2_tracking.m
│   ├── Enhancer_promoter_MS2_Hub_tracking.m
│   ├── MS2_tracking_burst_analysis.m
│   │
│   │   # --- hub analysis (§3.2) ---
│   ├── Hub_counts_per_nucleus.m
│   ├── Hub_size_2D.m
│   ├── Hub_tracking_lifetime.m
│   ├── BRD4_MED14_Hub_colocalization.m
│   ├── Enhancer_MS2_BRD4_MED14_Hub_colocalization.m
│   │
│   │   # --- fixed-timepoint variants and calibration ---
│   ├── Enhancer_Hub_MS2_singleImages.m
│   ├── Enhancer_promoter_MS2_Hub_singleImages.m
│   ├── Bead_tracking.m
│   │
│   ├── functions/                    # shared helper functions
│   ├── bioformats/                   # Add bioformats_package.jar before starting (1.1)
│   └── HMM fitting/                  # wrapper function only, add vbFRET package before starting (1.2)
│
└── Data Analysis/                    # downstream analysis and plotting

```
Scripts resolve these folders relative to their own location, so the layout must
be preserved. Run scripts from their own directory.

---

## 3. Pipelines

### 3.1 Locus tracking

| Script | Channels | Tracks | Burst state | Hubs |
|---|---|---|---|---|
| `Enhancer-Promoter_tracking.m` | 2 | enhancer, promoter | – | – |
| `Enhancer_promoter_MS2_tracking.m` | 3 | enhancer, promoter, MS2 | HMM | – |
| `Enhancer_Hub_MS2_tracking.m` | 3 | enhancer, MS2 | HMM | 3 nearest |
| `Enhancer_promoter_MS2_Hub_tracking.m` | 4 | enhancer, promoter, MS2 | HMM | nearest to E and to P |
| `MS2_tracking_burst_analysis.m` | 1 | MS2 sites (many per field) | HMM + burst calling | – |

All five share one architecture: 2D peak detection and Gaussian localisation on
the MIP, 3D Gaussian refinement of the corresponding z-stack subvolume, and
chained seeding — the enhancer is tracked first from an operator-supplied seed
and its position seeds the search for every other channel. Distances are
reported in nm.

### 3.2 Hub analysis

| Script | Input | Measures |
|---|---|---|
| `Hub_counts_per_nucleus.m` | single image | hub counts, area, peak/mean enrichment |
| `Hub_size_2D.m` | folder batch | FWHM (2D Gaussian and radial), Feret diameter, aspect ratio |
| `Hub_tracking_lifetime.m` | timecourse | structural lifetimes, Kaplan–Meier survival |
| `BRD4_MED14_Hub_colocalization.m` | 2-channel image | nearest-neighbour distances, Manders M1/M2 |
| `Enhancer_MS2_BRD4_MED14_Hub_colocalization.m` | 4-channel image | locus-to-hub distances, per-hub MED14 colocalization |

**Two detection regimes**, and the distinction matters when comparing numbers
across scripts:

- **DBSCAN** — used for raw hub counts and photometry
  (`Hub_counts_per_nucleus.m`). Clustering recovers whole hub objects with their
  true footprint, which is what a count and an area measurement require.
- **3× nuclear mean intensity** — used everywhere else. Peak detection plus
  centroid refinement gives one stable sub-pixel centre per hub, which is what
  sizing, tracking, and distance measurements need.


---

## 4. Running a pipeline

1. Install the dependencies in §1.
2. Open the script and set the input paths (`Input_zstack`, `MIP_filename`),
   `scene`, and `save_filename`.
3. Set the channel indices for your acquisition.
4. Set the seed coordinate `C_in` and the frame range (`start_f`, `stop_f`).
   `MS2_tracking_burst_analysis.m` instead opens an interactive picker
   (`pick_TS_spots`) — click each transcription site.
5. Set search radii and threshold multipliers; hub scripts also need `nuc_ave`,
   a per-movie background constant.
6. Run. Every pipeline draws a per-frame overlay — watch it and abort if the
   marker leaves the locus.
7. Results are written to `save_filename.mat`.

The parameter values used for each dataset in the manuscript are recorded
per script in the two metadata files.

**Note on paths.** The scripts ship with absolute paths from the acquisition
machines (`/Volumes/...`). Replace them with your own before running.

---

## 5. Outputs

Each pipeline saves one `.mat` file. Variable names, dimensions, units, and
conventions — including the `MS2_score` column key and the `C_cent` row offset,
where row `k+1` holds the frame-`k` result — are tabulated in the data
dictionary sections of the two metadata files.

Distances are in nm, positions in µm, and lateral pixel size is hard-coded to
0.0313 µm; axial spacing is read from the CZI metadata.

---

## 6. Citation

If you use this code, please cite the associated manuscript (details to follow)
and the third-party methods it depends on:

- **vbFRET** — Bronson et al., *Biophysical J.* 97(12):3196–205 (2009)
- **Cellpose** — Stringer et al., *Nature Methods* 18:100–106 (2021)
- **Bio-Formats** — Linkert et al., *J. Cell Biol.* 189(5):777–782 (2010)
- **Particle tracking** (`pkfnd`, `cntrd`) — Crocker & Grier,
  *J. Colloid Interface Sci.* 179:298–310 (1996)

---

## 7. License

MIT (see `LICENSE`). The license covers the code in this repository only.
Bio-Formats, vbFRET, and Cellpose carry their own terms and are not
redistributed here.

---

## 8. Contact

Aniket Jana — Xie Lab, Lerner Research Institute, Cleveland Clinic.

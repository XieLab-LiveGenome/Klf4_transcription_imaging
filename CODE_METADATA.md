---
title: Live-cell SR-SIM tracking pipelines for enhancer–promoter dynamics, MS2 transcriptional bursting, and BRD4 hub association
language: MATLAB
type: analysis-code
components: 5 top-level pipeline scripts + shared helper functions
imaging_platform: Zeiss Elyra 7 Lattice SIM
input_format: Carl Zeiss Image (.czi), 6D (S, T, Z, C, Y, X)
output_format: MATLAB .mat
lateral_pixel_size_um: 0.0313
axial_pixel_size_um: read from CZI metadata (metadata.ScaleZ)
distance_units_reported: nanometres
---

# Code metadata — tracking pipelines

This document describes the MATLAB code used to extract 3D positions, pairwise
distances, and transcriptional states from multicolour live-cell SR-SIM
timecourses. It is written to accompany a code deposit and to support the
Methods / Code Availability section: it records what each script does, what it
consumes, what it emits, and every parameter that was set by the user rather
than derived from the data.

---

## 1. Scope and design

All five pipelines share one architecture:

1. **Load** a CZI z-stack (and, in most scripts, a separately saved maximum
   intensity projection, MIP) through Bio-Formats.
2. **Detect and track in 2D** on the MIP — local maxima (`pkfnd`) inside a
   search radius centred on the previous frame's position, followed by a 2D
   Gaussian fit of the brightest candidate (`track_spot2D`).
3. **Refine in 3D** by fitting an anisotropic 3D Gaussian to a subvolume of the
   full z-stack seeded at the 2D result (`fit_Gaussian3D`).
4. **Chain the channels**: the enhancer is tracked first and its position seeds
   the search for the promoter, the MS2 site, and the protein hub. Only the
   enhancer is tracked de novo from a user-supplied seed coordinate.
5. **Post-process**: reject or interpolate implausible axial jumps, convert to
   physical units, compute pairwise 3D and 2D separations in nm.
6. **Save** all trajectory and distance arrays to a `.mat` file.

Tracking is single-locus and semi-automated by design: the operator supplies one
seed coordinate per allele (`C_in`) and the frame range, and the run is verified
frame by frame through an on-screen overlay.

---

## 2. Pipeline inventory

| Script | Channels | Tracks | Transcription state | Hub detection | Primary outputs |
|---|---|---|---|---|---|
| `Enhancer-Promoter_tracking.m` | 2 | enhancer, promoter | – | – | E–P distance |
| `Enhancer_promoter_MS2_tracking.m` | 3 | enhancer, promoter, MS2 | HMM | – | E–P distance + burst state |
| `Enhancer_Hub_MS2_tracking.m` | 3 | enhancer, MS2, hubs | HMM | 3 nearest hubs | E–hub, TSS–hub distances |
| `Enhancer_promoter_MS2_Hub_tracking.m` | 4 | enhancer, promoter, MS2, hubs | HMM | nearest-to-E and nearest-to-P | E–P, E–hub, P–hub distances |
| `MS2_tracking_burst_analysis.m` | 1 | MS2 transcription sites (multiple per field) | HMM + burst calling | – | ON/OFF durations, burst size and amplitude |

### 2.1 `Enhancer-Promoter_tracking.m`
Two-colour, fast-interval (5 s) E–P imaging. This is the only script that builds
the MIP internally (per-timepoint, per-channel `max` along z) rather than reading
a pre-computed MIP file. Channel 1 = enhancer, channel 2 = promoter.
Configuration used in the deposited copy: `C_in = [53 82]`, frames 1–103,
`R = 30`, `R_p = 20`, `R_fit = 6`, `fit_thr = 0.2`, `m_enh = m_prom = 0.2`,
`thr_z = 0.3` µm.
Saved: `enh_xyz`, `prom_xyz`, `C_cent_e`, `C_cent_p`, `E_P_dist_3d`,
`E_P_dist_2d`, `r2`.

### 2.2 `Enhancer_promoter_MS2_tracking.m`
Three-colour E–P imaging with a simultaneous MS2 readout. Channel 1 = MS2
(mStayGold), channel 2 = enhancer, channel 3 = promoter. Reads the z-stack and a
pre-computed MIP as two separate CZI files. Transcriptional state is inferred by
two-state HMM fitting of the MS2 intensity trace (`HMM_fit_fun`, binarised by
`binary`).
Configuration used: `C_in = [193 189]`, frames 1–96, `R = 50`, `R_p = 25`,
`R_m = 20`, `R_fit = 10`, `fit_thr = 0.2`, `m_enh = m_prom = 0.1`,
`m_MS2 = 0.2`, `thr_z = 0.5` µm.
Saved: `enh_xyz`, `prom_xyz`, `TSS_xyz`, `C_cent_e`, `C_cent_p`, `C_cent_TSS`,
`MS2_score`, `E_P_dist_3d`, `E_P_dist_2d`, `r2`.

### 2.3 `Enhancer_Hub_MS2_tracking.m`
Three-colour enhancer / MS2 / protein-hub imaging. Channel 1 = MS2, channel 2 =
enhancer (SNAP-tag JF552), channel 3 = hub (HaloTag JF642). Hubs are segmented
on the MIP by adaptive thresholding — the threshold is `nuc_ave × mult`, scaled
each frame by the ratio of that frame's whole-image intensity SD to the first
frame's, which compensates for photobleaching — then filtered to ≥10 px objects
and passed with the enhancer position to `condensate_search_v6`, which returns
the three hubs nearest the enhancer with their 3D position, size, integrated
brightness, and aspect ratio.
Configuration used: `C_in = [99 96]`, frames 1–50, `R = 40`, `R_m = 30`,
`R_fit = 10`, `m_enh = m_MS2 = 0.1`, `thr_z = 0.5` µm, `nuc_ave = 1500`,
`mult = 3.0`, `R_search_c = 10` px.
Saved: `MS2_score`, `enh_xyz`, `TSS_xyz`, `Cond_xyz1..3`, `Cond_s1..3`,
`Cond_int1..3`, `E_TSS_dist`, `E_C_dist1..3`, `TSS_cond_dist`.

### 2.4 `Enhancer_promoter_MS2_Hub_tracking.m`
Four-colour pipeline: channel 1 = MS2 (405), channel 2 = enhancer (488),
channel 3 = promoter (561), channel 4 = hub (642). Differs from the three-colour
versions in two respects: (i) frames with no detected punctum are written as
`NaN` and later filled by linear interpolation (`fillmissing`) instead of
carrying the previous position forward, and (ii) `condensate_search_v9` is called
twice per frame, once seeded at the enhancer and once at the promoter, so the
enhancer-proximal and promoter-proximal hubs are resolved separately. This yields
both matched distances (E to its own nearest hub, P to its own) and cross
distances (E to the promoter's hub and vice versa).
Configuration used: `C_in = [262 225]`, frames 1–91, `R = 40`, `R_p = 20`,
`R_m = 20`, `R_fit = 10`, `m_enh = m_prom = 0.2`, `m_MS2 = 0.1`,
`thr_z = 0.6` µm, `nuc_ave = 4000`, `mult = 3.0`, `R_search_c = 20` px (≈600 nm).
Saved: `enh_xyz`, `prom_xyz`, `TSS_xyz`, `C_cent_e`, `C_cent_p`, `C_cent_TSS`,
`MS2_score`, `E_P_dist_3d`, `E_P_dist_2d`, `r2`, `Cond_xyz1..2`, `Cond_s1..2`,
`Cond_int1..2`, `E_C_dist1..2`, `P_C_dist1..2`.

### 2.5 `MS2_tracking_burst_analysis.m`
Single-channel MS2 burst quantification over many transcription sites in one
field. Operates on the MIP only (no z-stack). Four stages, each in its own
function: interactive site selection (`pick_TS_spots`), nucleus motion estimation
(`nuc_motion_MS2` or `track_nuclei_MS2`), motion-corrected site tracking
(`track_TS_spots`), and HMM fitting plus burst quantification
(`burst_analysis_MS2`).

Two switches control behaviour:
- `run_nuclei` — compute per-nucleus motion.
- `motion_correct` — shift the search seed by the host nucleus's frame-to-frame
  displacement before searching. Measurement is always taken at the fitted spot,
  so correction changes *where the search is centred*, not *what is measured*.
- `nuc_method` — `'blocks'` (recommended; `nuc_motion_MS2`) or `'masks'`
  (`track_nuclei_MS2`, per-frame Cellpose centroids).

Intensities are rescaled before HMM fitting by the ratio of the first frame's
whole-image intensity SD to the current frame's, correcting for global bleaching.
Configuration used: `time_int = 3`, `nuc_block = 10`, `nuc_preproc = 'none'`,
`pad = 100`, `sigma = 4`, `max_step = 150`, `cell_D_seg = 160`, `Rp = 80`,
`Rn = 150`, `m_MS2 = 0.3`, `R_c = 5`.
Saved (19 variables): `C_in`, `C_cent`, `sd`, `spot_size`, `I_t`, `I_t_norm`,
`fitMS2_2s`, `binary_2s`, `on_times`, `off_times`, `on_cellavg`, `off_cellavg`,
`on_compiled`, `off_compiled`, `Burst`, `Burst_amp`, `Burst_compiled`,
`Burst_amp_compiled`, `burst_cellavg`, `burst_amp_cellavg`.

---

## 3. Functions included in this deposit

| Function | Role |
|---|---|
| `fit_Gaussian3D.m` | 3D Gaussian refinement of a seed coordinate in a z-stack |
| `track_spot2D.m` | Peak detection + 2D Gaussian localisation within a search radius |
| `pick_TS_spots.m` | Interactive selection of transcription sites on a block average |
| `nuc_motion_MS2.m` | Per-nucleus rigid translation from block averages by phase correlation |
| `track_TS_spots.m` | Frame-by-frame MS2 site tracking with optional motion correction |

### 3.1 `fit_Gaussian3D(stack, guess, radius, spacing)`
Fits `A·exp(−[(x−x₀)²/2σx² + (y−y₀)²/2σy² + (z−z₀)²/2σz²]) + B` (8 free
parameters) to a subvolume by `lsqcurvefit`. The lateral ROI is `±radius` about
the guess; **the full z range is always used**. Fitting is performed in physical
units (µm) and the returned coordinates are converted back to voxel indices.
Lateral centre bounds are constrained to `±radius` of the guess; axial centre is
bounded to the stack. Returns `NaN` for coordinates, parameters, and `R²` on any
failure (out-of-bounds guess, empty or constant subvolume, fewer than 10 valid
points, non-finite fit), so callers must test for `NaN`. `R²` is
`1 − SS_res/SS_tot` over the fitted subvolume and is used upstream as the fit
quality gate (`fit_thr = 0.2`).

### 3.2 `track_spot2D(image, x_guess, y_guess, R_search, m)`
Thresholds at `m × max(image)`, runs `pkfnd` with a minimum peak separation of
3 px, keeps peaks within `R_search` of the guess, selects the brightest, crops a
17×17 px ROI (`FIT_R = 8`) about it, and fits a 6-parameter 2D Gaussian
(`fit2DGaussian`, initial σ = 3 px). Returns
`[x_refined, y_refined, n_candidates]` in global image coordinates; if no
candidate is found it returns the input guess with `n_candidates = 0`, which the
callers use as the miss flag.

### 3.3 `pick_TS_spots(MIP_image6d, channel, ...)`
Averages `n_frames` (default 10) from `start_f`, displays the average
(optionally inverted, with 1.0–99.9 percentile display limits), and collects
clicks. Duplicate clicks within 2 px are dropped and each click is snapped to the
brightest pixel within `snap_R` (default 8 px) on the raw average. Returns
`C_in` as `[x y]`. Passing `'manual', C_in_manual` bypasses picking, which is how
a previous selection is replayed for reproducibility.

### 3.4 `nuc_motion_MS2(MIP_image6d, channel, C_in, ...)`
Estimates nucleus motion without per-frame segmentation. Cellpose is run **once**
on the whole-movie time average to define nucleus ROIs and to assign each
transcription site to a nucleus; the movie is then averaged in blocks
(`block_size`, default 10 frames) and each nucleus ROI in block *b* is
phase-correlated against block 1, with parabolic sub-pixel refinement, Hann
windowing, and rejection of block-to-block jumps exceeding `max_step`. Block
shifts are interpolated onto every frame (`pchip`). Field names match
`track_nuclei_MS2`, so it is a drop-in replacement.

*Rationale, recorded here because it motivates the choice:* on weak diffuse
nuclear SIM signal, Cellpose mask area fluctuates ~25% frame to frame, moving the
mask centroid by tens of pixels when the nucleus has not moved. Registration
compares image content and is insensitive to where the mask boundary landed.
Blur must stay small — σ ≈ 4 px gives ~0.1 px error, whereas heavy blur (σ ≈ 30)
biases the estimate by ~2 px.

### 3.5 `track_TS_spots(MIP_image6d, channel, C_in, nuc, motion_correct, ...)`
Tracks each site with `track_spot2D`, seeding frame *k* from the frame *k−1*
result plus, if `motion_correct` is true, the host nucleus's displacement between
those frames. On a hit, the next search radius is `Rp` (80 px); on a miss the
position is held and the radius widens to `Rn` (150 px). Intensity is the mean in
a fixed circular mask of radius `R_c` (5 px) centred on the fitted position.
Requires `start_f = 1` when motion correction is on, so nucleus shifts and tracked
frames are indexed identically.

---

## 4. External dependencies (not included in this deposit)

**Must be on the MATLAB path** (expected in `./functions`, `./bioformats`, and
`./HMM fitting` relative to the calling script):

| Dependency | Used by | Purpose |
|---|---|---|
| Bio-Formats (`bioformats_package.jar`), `ReadImage6D2` | all | CZI reading into a 6D array + metadata struct |
| `single_zstack` | all 3D scripts | extract one (timepoint, channel) z-stack |
| `SNR_inc2` | all | contrast-stretch denoising before detection |
| `pkfnd` | `track_spot2D` | local maximum detection (Crocker–Grier lineage) |
| `fit2DGaussian` | `track_spot2D` | 6-parameter 2D Gaussian fit |
| `maskavg` | MS2 scripts | mean intensity in a circular mask |
| `HMM_fit_fun`, `binary` | MS2 scripts | two-state HMM fit and binarisation |
| `condensate_search_v6` / `_v9` | hub scripts | nearest-hub search and hub descriptors |
| `cellpose_seg_MS2` | `nuc_motion_MS2` | one-off nuclear segmentation |
| `track_nuclei_MS2` | burst script (`'masks'` mode) | per-frame Cellpose segmentation and linking |
| `burst_analysis_MS2` | burst script | HMM fit, ON/OFF durations, burst size and amplitude |

**MATLAB toolboxes:** Image Processing (`imbinarize`, `bwareafilt`,
`regionprops`, `imgaussfilt`, `imshow`, `getpts`), Computer Vision
(`insertObjectAnnotation`), Optimization (`lsqcurvefit`). Cellpose access is
required for the burst pipeline only. `pick_TS_spots` deliberately avoids the
Statistics Toolbox (percentiles computed by sorting).

---

## 5. Data dictionary

Array dimensions are `(frames × components × loci)` unless noted. `zs` is the
number of analysed frames, `vis` the number of tracked loci (or transcription
sites).

| Variable | Size | Units | Meaning |
|---|---|---|---|
| `C_in` | `vis × 2` | px | user-supplied seed, `[x y]` |
| `enh_xyz`, `prom_xyz`, `TSS_xyz` | `zs × 3 × vis` | µm | fitted 3D positions |
| `Cond_xyz1..3` | `zs × 3 × vis` | µm | hub 3D positions, ranked by proximity |
| `C_cent_e`, `C_cent_p` | `(zs+1) × 3 × vis` | px / slice | voxel-index positions; **row `k+1` holds the frame-`k` result**, row 1 holds the seed |
| `C_cent_TSS` | `zs × 3 × vis` | px / slice | rounded MS2 voxel position, row `k` = frame `k` |
| `r2` | `zs × 2 × vis` | – | 3D fit R²; column 1 = enhancer, column 2 = promoter |
| `E_P_dist_3d`, `E_P_dist_2d` | `zs × 2 × vis` | nm | column 1 from the 3D-fit coordinates, column 2 the alternate estimate using the 2D-tracked lateral positions |
| `E_C_dist*`, `P_C_dist*`, `E_TSS_dist`, `TSS_cond_dist` | `zs × 2 × vis` | nm | column 1 = 3D, column 2 = 2D (lateral only) |
| `Cond_s*`, `Cond_int*`, `Cond_AR*` | `zs × vis` | px², a.u., – | hub size, integrated brightness, aspect ratio in the MIP |
| `I_t`, `I_t_norm` | `zs × vis` | a.u. | raw and SD-corrected MS2 intensity |
| `sd` | `zs × vis` | a.u. | whole-frame intensity SD (bleaching proxy) |
| `on_times`, `off_times` | cell | frames / `time_int` | burst start, end, duration |
| `Burst`, `Burst_amp` | – | a.u. | per-burst integrated size and amplitude |

### `MS2_score` (`zs × 5 × vis`)

| Column | Content |
|---|---|
| 1 | detection flag: `1` = MS2 spot found, `−1` = no spot (intensity then sampled at the promoter/enhancer position) |
| 2 | raw masked MS2 intensity (radius 5 px) |
| 3 | binarised ON/OFF state (`1`/`0`) |
| 4 | two-state HMM fit of column 2 |
| 5 | column 2 normalised to that trajectory's maximum |

Downstream analysis scripts read column 3 directly rather than re-running the HMM.

---

## 6. Conventions, thresholds, and known limitations

**Coordinates.** All image coordinates are `[x y]` (column, row), matching
`track_spot2D` output. Lateral physical positions are `x·xpixel`; axial positions
are `(z−1)·zpixel`, i.e. z is zero-based and xy is one-based. This offsets the
absolute origin but not any reported separation, since all distances are
differences within a frame.

**Pixel size.** `xpixel` and `ypixel` are hard-coded to 0.0313 µm (Elyra 7 SIM,
lattice reconstruction) rather than read from metadata; `zpixel` is taken from
`metadata.ScaleZ`. Any dataset acquired at a different lateral sampling requires
these constants to be edited.

**Detection thresholds** (`m_enh`, `m_prom`, `m_MS2`) are fractions of the
per-frame maximum, so they adapt to intensity but are sensitive to bright
out-of-locus objects in the field. `nuc_ave` in the hub scripts is a per-movie
background constant set by the operator (1500 and 4000 in the two deposited
configurations) and must be re-set for a new acquisition or laser setting.

**Axial spike filtering.** Frames whose axial displacement exceeds `thr_z`
(0.3–0.6 µm depending on the frame interval), or whose z is negative, are
replaced by the mean of their two neighbours. This is a bad-fit rejection step,
not a smoothing step; it is applied once, non-iteratively, over frames 2 to
`zs−1`, so the first and last frames are never corrected.

**Missing frames.** The three-colour scripts carry the previous position forward
when no punctum is detected; the four-colour script writes `NaN` and interpolates
linearly. Trajectories from the two families are therefore not identically
gap-handled, and comparisons should either be restricted to detected frames or
use the detection flag (`MS2_score` column 1, `NaN` runs in `enh_xyz`).

**Chained seeding.** Promoter, MS2, and hub searches are all seeded from the
enhancer position. A lost enhancer track propagates to every other channel in
that frame; the per-frame overlay exists so this is caught during the run.

**Frame interval** is recorded in each script's header comment (5 s for the
two- and three-colour E–P data; longer for the hub timecourses) and in
`time_int` for the burst pipeline. It is not read from the CZI, so it should be
verified against the acquisition metadata for each dataset.

**Interactive steps.** `pick_TS_spots` and the seed coordinate `C_in` require an
operator. Re-running an analysis reproducibly requires the printed
`C_in_manual` block (emitted by `pick_TS_spots`) or the `C_in` saved in the
output `.mat`.

---

## 7. Reproducing a run

1. Place `functions/`, `bioformats/`, and `HMM fitting/` beside the pipeline
   script; the script adds them to the path and registers the Bio-Formats JAR.
2. Set `Input_zstack`, `MIP_filename`, `scene`, and `save_filename`.
3. Set the channel indices, `C_in`, `start_f`, and `stop_f`.
4. Set the search radii, threshold multipliers, and (hub scripts) `nuc_ave`.
5. Run; watch the per-frame overlay and abort if the marker leaves the locus.
6. Output is written to `save_filename.mat` in the working directory.

---

## 8. Draft text for the Code Availability statement

> Image analysis was performed with custom MATLAB code. Loci were localised in
> three dimensions by 2D peak detection and Gaussian fitting on maximum
> intensity projections followed by 3D Gaussian refinement of the corresponding
> z-stack subvolume; transcriptional states were called by two-state HMM fitting
> of MS2 intensity traces; protein hubs were identified by adaptive intensity
> thresholding and a nearest-neighbour search in the locus neighbourhood. All
> pipeline scripts, helper functions, and the parameter values used for each
> dataset are available at [repository / DOI]. Bio-Formats and Cellpose are
> third-party dependencies and are not redistributed.

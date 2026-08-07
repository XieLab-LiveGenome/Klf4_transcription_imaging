---
title: Hub detection, sizing, lifetime, and colocalization pipelines (BRD4 / MED14)
language: MATLAB
type: analysis-code
components: 5 top-level pipeline scripts
imaging_platform: Zeiss Elyra 7 Lattice SIM
input_format: Carl Zeiss Image (.czi), 6D (S, T, Z, C, Y, X); single images, batches, and timecourses
output_format: MATLAB .mat
lateral_pixel_size_um: 0.0313
distance_units_reported: nanometres
detection_regimes: DBSCAN (hub counts and photometry) | 3x nuclear mean intensity (all other scripts)
status: semi-final
---

# Code metadata — hub analysis

Companion to `CODE_METADATA.md` (locus tracking). This document covers the code
used to detect nuclear protein hubs, measure their size and enrichment, track
their structural lifetimes, and quantify BRD4/MED14 colocalization.

---

## 1. Two detection regimes

All hub detection is performed inside Cellpose-segmented nuclei, so every
threshold is referenced to that nucleus's own signal rather than to a global
image value. Two distinct detectors are used, and which one applies depends on
the quantity being measured:

**Density-based clustering (DBSCAN) — used for raw hub counts and photometry.**
`Hub_counts_per_nucleus.m` calls `detectPunctaDBSCAN` on the raw nuclear image
with a top-hat radius of 30 px, `Epsilon = 8` px, and a minimum cluster area of
40 px, thresholded at the nuclear mean + 2 SD. Clustering recovers whole hub
*objects* with their true footprint, which is what a count and an
area/enrichment measurement require; a peak finder would collapse a large or
irregular hub into one point and would split an elongated one into two.

**Peak detection at 3× nuclear mean intensity — used everywhere else.**
`Hub_size_2D.m`, `Hub_tracking_lifetime.m`, `BRD4_MED14_Hub_colocalization.m`,
and `Enhancer_MS2_BRD4_MED14_Hub_colocalization.m` threshold at
`thresh_mult × mean(nuclear pixels)` with `thresh_mult = 3` and locate hubs with
`pkfnd` followed by Crocker–Grier centroid refinement (`cntrd`). Sizing,
tracking, and distance measurements all need a single sub-pixel centre per hub
and a stable, intensity-referenced criterion for what counts as a hub across
frames, nuclei, and conditions — which is what the 3× rule provides.

One deliberate exception sits inside the colocalization scripts: the **Manders
coefficients use a gentler 2× nuclear mean** threshold to define the "signal
present in the other channel" mask. The 3× rule defines hub identity; the 2×
rule defines where the partner protein is detectable at all, and using the
stricter value there would systematically depress M1/M2.

---

## 2. Pipeline inventory

| Script | Input | Detection | Nuclear segmentation | Primary outputs |
|---|---|---|---|---|
| `Hub_counts_per_nucleus.m` | single MIP image | DBSCAN | `cellpose_seg`, diameter 150, border nuclei removed | hub counts, areas, peak/mean intensity normalised to nuclear mean |
| `Hub_size_2D.m` | folder batch of MIP images | `pkfnd` @ 3× | `cellpose_seg`, diameter 150 | FWHM (2D Gaussian and radial), Feret diameter, aspect ratio, local background |
| `Hub_tracking_lifetime.m` | 2-channel timecourse (SIM + widefield) | `pkfnd` + `cntrd` @ 3× | `cellpose_seg2` on the widefield channel, per frame | structural lifetimes, Kaplan–Meier survival, tracks |
| `BRD4_MED14_Hub_colocalization.m` | 2-channel single image | `pkfnd` + `cntrd` @ 3× | `cellpose_seg` on the BRD4 channel | nearest-neighbour distances, Manders M1/M2, colocalized vs isolated classification |
| `Enhancer_MS2_BRD4_MED14_Hub_colocalization.m` | 4-channel single image + z-stack | `condensate_search_v6` + `pkfnd` @ 3× | `cellpose_seg` on the summed MED14+BRD4 image | enhancer/TSS-to-hub distances, per-hub MED14 colocalization flags, M2 |

---

## 3. Pipeline descriptions

### 3.1 `Hub_counts_per_nucleus.m` — hub counts and enrichment

Single-image pipeline. Nuclei are segmented with `cellpose_seg` (diameter 150)
and **border nuclei are removed**: a nucleus is rejected if the total number of
its pixels touching any image edge exceeds 5 µm worth of pixels
(`distance_threshold_um / xpixel`), which discards truncated nuclei while
tolerating a nucleus that merely grazes the frame.

Within each retained nucleus, `detectPunctaDBSCAN` is run on the raw masked
image with `Epsilon = 8` px, `MinArea = 40` px, `TophatRadius = 30`, and a
threshold of `mean + 2·SD` over the nuclear pixels. Detection output is
relabelled with `bwlabel` and intersected with the nuclear mask; a mismatch
between the detector's region count and `bwlabel`'s is reported as a warning
rather than an error.

**Detection and photometry are deliberately separated**: clusters are found on
the detection image, but every intensity is read from the raw image. Per hub the
script records area (µm²), centroid, peak and mean raw intensity, and both
normalised to that nucleus's mean intensity — the peak ratio is the enrichment
metric and the mean ratio is the partition coefficient. Per nucleus it records
the whole-nucleus mean, the nucleoplasmic mean (nuclear pixels excluding hub
pixels), and the median peak ratio.

Results are pooled across nuclei with `pooled_nucleus_id` retained so that
per-nucleus structure is recoverable for statistics.
Configuration used: `channel_cond = 3` (642), `eps_val = 8`, `min_pixels = 40`.

### 3.2 `Hub_size_2D.m` — hub size, batch

Loops over every `.czi` in `folderPath` and writes one `.mat` per input file,
named after the source. Nuclei are segmented (`cellpose_seg`, diameter 150) and
**eroded by 10 px** before peak finding, so hubs sitting on the segmentation
boundary — where the local background estimate would be unreliable — are
excluded. Peaks are found with `pkfnd` at `3 × nuclear mean` and a minimum
separation of 25 px, then refined with `cntrd` (window 21 px). Points rejected by
`cntrd` near image edges are matched back by nearest distance so that
`nucleus_id` stays aligned.

Each hub is then measured in a 41×41 px ROI (`R_search = 20`) after subtracting a
radially estimated local background (`local_background_2Dg`), with negatives
clipped to zero:

- **Feret diameter** — Otsu threshold (`graythresh`) of the normalised ROI, keep
  only the connected object containing the ROI centre, then
  `(MaxFeret + MinFeret)/2 × xpixel × 1000` in nm, with the max/min ratio stored
  as an aspect ratio. Falls back to `feretDiameters` if `regionprops` fails.
- **FWHM (radial)** — `radial_spot_size2D` on the radial intensity profile.
- **FWHM (2D Gaussian)** — `gaussian_spot_size2D`, which also returns the
  aspect ratio and a refined centre (converted back to global coordinates).

Configuration used: `channel_cond = 1`, `pkfnd_sz = 25`, `edge_margin = 10`,
`thresh_mult = 3`, `sz_cntrd = 21`, `R_search = 20`.

### 3.3 `Hub_tracking_lifetime.m` — structural lifetime

Two-channel timecourse: channel 1 is the SIM BRD4 signal used for detection,
channel 2 is the widefield BRD4 signal used for segmentation. Widefield is used
for segmentation because the nuclear signal there is more homogeneous than in the
SIM reconstruction; detection and all measurements stay on SIM. Ten steps:

1. **Segmentation with graded fallback** (`pick_nucleus`): Cellpose
   (`cellpose_seg2`, diameter 120 at 0.3× downsampling) on the widefield frame,
   candidates filtered at `MIN_NUC_AREA = 30000` px and the one nearest the
   previous centroid selected, rejected if it jumps more than
   `MAX_CENT_JUMP = 80` px. On failure the average of the current and previous
   widefield frames is retried; on second failure the previous frame's mask is
   held. Each frame's outcome is recorded in `seg_status`
   (`direct`/`combined`/`fallback`), so segmentation quality is auditable after
   the run.
2. **Defective-frame detection**: the 99th percentile of SIM intensity inside
   the mask is compared with its own 5-frame median trend; frames with a robust
   z (MAD-scaled) below `DEFECT_Z = −2.0` are flagged. These are SIM
   reconstruction dropouts, not biology.
3. **Multiplicative correction**: flagged frames are rescaled to the linearly
   interpolated trend of the good frames. Correction is applied to intensities
   only, and the scale factors are saved.
4. **Peak detection**: Gaussian blur (σ = 1.0), `pkfnd` at
   `THR_MULT = 3 × mean nuclear intensity` with `PEAK_SEPARATION = 25` px, then
   `cntrd` (feature size 5), keeping only peaks inside the mask.
5. **Rigid-body registration**: peaks are transformed into a common reference
   frame — translation from the median-filtered nuclear centroid to the median
   centroid, plus rotation to the median orientation. **Rotation is gated**: it
   is applied only when the frame's orientation differs by more than
   `ROT_DEG_GATE = 3°` *and* the nucleus is elongated enough for orientation to
   be meaningful (`ROT_ECC_GATE = 0.60` eccentricity). This prevents the
   near-arbitrary orientation of a round nucleus from injecting spurious
   rotation.
6. **Linking**: Hungarian assignment (`matchpairs`) with a maximum displacement
   of `MAX_DISP_PX = 45` px and `MEMORY = 1`, so a hub may vanish for one frame
   without terminating its track.
7. **Lifetimes with edge censoring**: a track's duration is last-seen minus
   first-seen plus one, converted with `INTERVAL_SEC = 120`. Tracks present in
   frame 1 are left-censored and tracks present in the final frame are
   right-censored.
8. **Two summaries with different censoring rules**: the histogram uses only
   fully uncensored tracks (both birth and death observed), while the
   Kaplan–Meier curve drops left-censored tracks and treats right-censored
   tracks as non-events. The KM estimator is computed inline, not via a
   toolbox function.

Output is one `-v7.3` struct containing lifetimes, full tracks, masks,
centroids, registration shifts and rotations, the defective-frame list and scale
factors, per-frame peaks in both raw and reference coordinates, `seg_status`,
and the complete parameter struct `p`.

### 3.4 `BRD4_MED14_Hub_colocalization.m` — two-colour colocalization

Single-image, two-channel (channel 1 = MED14/mScarlet, channel 2 =
BRD4/HaloTag). Nuclei are segmented once on the BRD4 channel and the **same
masks are used for both channels**, so per-nucleus statistics are paired by
construction.

A local function `detect_and_size` runs the identical detection and sizing
routine on each channel (the same measurement stack as `Hub_size_2D.m`:
`pkfnd` at 3× nuclear mean with separation 15 px → `cntrd` → local background →
Feret, radial FWHM, 2D Gaussian FWHM and aspect ratio). Note the mask handling
differs from `Hub_size_2D.m`: masks are **dilated** by 10 px rather than eroded,
so hubs at the nuclear rim are included, while the intensity threshold is still
computed from the undilated mask.

**Proximity classification.** All pairwise centroid distances are computed with
`pdist2`; each MED14 hub is assigned its nearest BRD4 hub and vice versa. A hub
is colocalized if that distance is at most `coloc_radius_nm = 300` nm
(≈9.6 px). The classification is computed twice — globally across the image, and
again restricted within each nucleus for the per-nucleus tables.

**Manders coefficients**, computed per hub in an 11×11 px ROI
(`R_manders = 5`) intersected with the host nuclear mask:

- `M1` (per MED14 hub) = fraction of MED14 intensity in the ROI that overlaps
  BRD4 pixels above `2 × that nucleus's BRD4 mean`.
- `M2` (per BRD4 hub) = fraction of BRD4 intensity in the ROI that overlaps
  MED14 pixels above `2 × that nucleus's MED14 mean`.

M1 and M2 are then reported split by the distance-based classification —
colocalized versus isolated, pooled and per nucleus — which is the comparison
the figures use.

Configuration used: `pkfnd_sz = 15`, `dilate_margin = 10`, `thresh_mult = 3`,
`sz_cntrd = 21`, `R_search = 20`, `R_manders = 5`, `coloc_radius_nm = 300`,
`manders threshold = 2 × nuclear mean`.

### 3.5 `Enhancer_MS2_BRD4_MED14_Hub_colocalization.m` — locus-anchored colocalization

Four-channel single-timepoint pipeline (1 = MS2, 2 = enhancer, 3 = MED14,
4 = BRD4) that asks the colocalization question *at the locus* rather than
across the whole nucleus. Both the MIP and the full z-stack are loaded, so all
distances are 3D.

Nuclei are segmented on the **sum of the MED14 and BRD4 images**, which gives a
more complete nuclear footprint than either channel alone. Per-nucleus means are
stored for both channels, together with global fallback means computed over the
union of all masks.

Per locus (seeded by operator-supplied enhancer coordinates in `C_in`):

1. The enhancer is refined in 3D (`spotPosition3D`, which wraps
   `fit_Gaussian3D`, `R_fit = 8`).
2. The MS2 site is found by `pkfnd` in a 30 px window about the enhancer at 5%
   of the global image maximum, the brightest peak is refined in 3D, and its
   intensity is read in a 5 px mask. If no peak is found, `MS2_score` column 1 is
   set to −1 and the intensity is sampled at the enhancer position instead.
3. The **host nucleus** is identified by looking up the enhancer position in the
   Cellpose label image (`point_location2D`). If the locus falls outside every
   mask, the script drops to a documented fallback: a synthetic circular ROI of
   `fallback_radius_px = 300` about the locus together with the global mean
   intensities, and `host_nuc_id` is set to 0 so those loci can be excluded
   downstream.
4. The three nearest BRD4 hubs are returned by `condensate_search_v6` seeded at
   the enhancer, with an **auto-expanding search radius**: the call is retried
   with `k × R_search_c` for `k = 2…20` while the primary hub is still NaN.
5. For each of those hubs, a **focused local MED14 search** runs `pkfnd` in a
   ±300 nm window at 3× the nuclear MED14 mean; every candidate is refined in 3D
   and the nearest is kept, with a colocalization flag set at 300 nm. Absence of
   any MED14 peak in the window is recorded as not colocalized, with distances
   left NaN — distinguishable from a hub that was never evaluated.
6. Manders M2 is computed per BRD4 hub exactly as in §3.4 (5 px ROI, 2× nuclear
   MED14 mean, undilated mask).
7. Enhancer–TSS, enhancer–hub, and TSS–hub distances are computed in 3D and 2D,
   in nm.

---

## 4. Dependencies (not included in this deposit)

| Dependency | Used by | Purpose |
|---|---|---|
| Bio-Formats, `ReadImage6D2` (and `ReadImage6D` in §3.4) | all | CZI reading |
| `cellpose_seg` / `cellpose_seg2` | all | nuclear segmentation; returns per-nucleus masks and a label image |
| `detectPunctaDBSCAN` | §3.1 | density-based punctum detection |
| `pkfnd`, `cntrd` | §3.2–3.5 | peak detection and Crocker–Grier centroid refinement |
| `local_background_2Dg` | §3.2, §3.4 | radial local background and its SD |
| `radial_spot_size2D`, `gaussian_spot_size2D` | §3.2, §3.4 | FWHM and aspect ratio |
| `feretDiameters` | §3.2, §3.4 | fallback Feret measurement |
| `matchpairs` | §3.3 | Hungarian linking (R2019b+) |
| `condensate_search_v6` | §3.5 | ranked nearest-hub search with descriptors |
| `single_zstack`, `fit_Gaussian3D`, `spotPosition3D` | §3.5 | z-stack extraction and 3D refinement |
| `point_location2D`, `imagemask`, `maskavg` | §3.5 | label lookup, ROI masking, masked mean intensity |

**MATLAB toolboxes:** Image Processing (`imbinarize`, `bwareafilt`, `bwlabel`,
`regionprops`, `imerode`, `imdilate`, `imgaussfilt`, `graythresh`,
`labeloverlay`, `visboundaries`, `medfilt1`), Computer Vision
(`insertObjectAnnotation`), Statistics (`pdist2`, `prctile`). Cellpose is
required by every script here.

---

## 5. Data dictionary

### Counts and enrichment (§3.1)

| Variable | Type | Units | Meaning |
|---|---|---|---|
| `all_condensate_counts` | `nNuc × 1` | – | hubs per nucleus |
| `all_condensate_Areas` | cell | µm² | per-hub area |
| `all_condensate_centroids` | cell | px, `[x y]` | per-hub centroid |
| `all_condensate_peak_raw`, `_mean_raw` | cell | a.u. | raw peak / mean hub intensity |
| `all_condensate_peak_int_norm` | cell | ratio | peak ÷ nuclear mean (enrichment) |
| `all_condensate_mean_int_norm` | cell | ratio | mean ÷ nuclear mean (partition coefficient) |
| `nucleus_mean_int`, `nucleus_nucleoplasm_int` | `nNuc × 1` | a.u. | whole-nucleus and hub-excluded means |
| `pooled_*` | vector | – | all hubs pooled across nuclei |
| `pooled_nucleus_id` | vector | – | nucleus index for each pooled hub |

### Size (§3.2, and the `med` / `brd` structs in §3.4)

| Variable | Units | Meaning |
|---|---|---|
| `condensate_cent` / `.cent` | px `[x y]` | Gaussian-refined centre |
| `condensate_FWHM` / `.FWHM` | nm | FWHM from the 2D Gaussian fit |
| `condensate_FWHM_r` / `.FWHM_r` | nm | FWHM from the radial profile |
| `condensate_feret` / `.feret` | nm | mean of max and min Feret diameter |
| `condensate_AR` / `.AR` | – | aspect ratio from the Gaussian fit |
| `condensate_AR_feret` / `.AR_feret` | – | max ÷ min Feret |
| `condensate_bg`, `condensate_sigma_r` | counts | local background and its SD |
| `nucleus_id` | – | host nucleus for each hub |

### Lifetime (§3.3)

| Variable | Units | Meaning |
|---|---|---|
| `structural_lifetime_sec` / `_frames` | s / frames | uncensored lifetimes only |
| `duration_sec` | s | all tracks, censored included |
| `first_seen`, `last_seen` | frame | track birth and death |
| `left_censored`, `right_censored` | logical | present in frame 1 / final frame |
| `tracks` | struct array | `.frames`, `.coords` (reference-frame `[y x]`) |
| `shifts_yx`, `dthetas_rad` | px, rad | applied registration |
| `defective_frames`, `scale_factors` | – | intensity correction record |
| `seg_status` | string | `direct` / `combined` / `fallback` per frame |
| `parameters` | struct | full parameter set `p` |

### Colocalization (§3.4)

Counts, colocalized and isolated fractions, and Feret diameters are stored per
nucleus and pooled; `dist_MED14_to_BRD4_nm` and `dist_BRD4_to_MED14_nm` hold
nearest-neighbour distances; `M1_coloc_*` / `M1_isol_*` and `M2_coloc_*` /
`M2_isol_*` hold the Manders values split by classification, each in both a
pooled (`_all`) and per-nucleus (`_per_nuc`) form. The full `med` and `brd`
detection structs are saved alongside.

### Locus-anchored colocalization (§3.5)

| Variable | Size | Units | Meaning |
|---|---|---|---|
| `enh_xyz`, `TSS_xyz` | `vis × 3` | µm | refined 3D positions; `TSS_xyz` is NaN when no MS2 spot |
| `MS2_score` | `vis × 2` | – | col 1 = spot flag (1 / −1), col 2 = masked intensity |
| `BRD4_xyz` | `vis × 3 × 3` | µm | top-3 hubs per locus, ranked by proximity |
| `BRD4_size`, `BRD4_int`, `BRD4_AR` | `vis × 3` | – | hub descriptors |
| `MED14_nearest_xyz`, `_d3D`, `_d2D` | `vis × 3 [× 3]` | µm, nm | nearest MED14 to each BRD4 hub |
| `MED14_coloc_flag` | `vis × 3` | logical | ≤ 300 nm |
| `M2_BRD4_per_cond` | `vis × 3` | – | per-hub Manders M2 |
| `E_TSS_dist` | `vis × 2` | nm | col 1 = 3D, col 2 = 2D |
| `E_BRD4_dist`, `TSS_BRD4_dist` | `vis × 2 × 3` | nm | to each of the three hubs |
| `host_nuc_id` | `vis × 1` | – | >0 real nucleus, 0 = synthetic fallback ROI |

---

## 6. Conventions and limitations

**Thresholds are per nucleus, not per image.** Every detection threshold is a
multiple of the host nucleus's own mean intensity, which makes counts and sizes
comparable across nuclei of differing expression level and across conditions,
but also means a nucleus with a high diffuse background will have a
correspondingly high bar for hub calling. `Hub_counts_per_nucleus.m` uses
`mean + 2·SD` instead, so its detections are not numerically interchangeable
with the 3× detections used elsewhere — counts from the two regimes should not
be pooled.

**Mask handling differs by intent.** `Hub_size_2D.m` *erodes* the nuclear mask
by 10 px so that local-background estimation is never contaminated by the
nuclear boundary; the colocalization scripts *dilate* by 10 px so that rim hubs
are not lost from a pairing analysis. Thresholds in both cases come from the
undilated mask. `Hub_counts_per_nucleus.m` does neither, but removes
border-touching nuclei outright.

**Colocalization is proximity-based, then intensity-based.** The 300 nm cutoff
classifies hub pairs by centroid distance; the Manders coefficients then measure
intensity overlap within that classification. The two are reported together
because neither alone distinguishes a genuinely shared hub from two adjacent
ones — the distance test can pair hubs that merely sit near each other, and
Manders alone carries no information about which hub the overlap belongs to.

**Censoring.** Structural lifetimes are truncated by the movie, not only by
biology. The histogram and the Kaplan–Meier curve use different censoring rules
(§3.3, step 8), so their medians are not expected to agree, and the histogram
median is biased short.

**Registration limits.** Peak coordinates are corrected for nuclear translation
and gated rotation only. Nuclear deformation is not corrected, so hubs far from
the centroid in a strongly deforming nucleus carry more residual apparent motion
than hubs near it.

**Single-timepoint scripts** (§3.1, §3.2, §3.4, §3.5) read frame 1 only. §3.2 is
the only batch script; the others process one file per run and the save path is
set by hand, so it must be checked against the input path before each run.

**Points to verify before deposit:**
- §3.4 calls `ReadImage6D`, whereas every other script calls `ReadImage6D2`.
- §3.5 binarises the BRD4 image at `mult × t_cond_BRD4 / 2`, i.e. 1.5× the
  nuclear mean, while passing `mult = 3` to `condensate_search_v6`; the mask is
  therefore permissive and the 3× criterion is enforced inside the search
  function.
- §3.5 sets the MS2 detection threshold at 5% of the global image maximum, the
  one threshold in this set that is not nucleus-referenced.
- §3.5's auto-expansion loop iterates `k = 2…20` and re-tests the NaN condition
  each pass, so it continues iterating after a successful retry without further
  effect; the final radius reached is not recorded.
- §3.1 reuses the loop variable `k` for both the frame index and the mask index.
- §3.1's declared dependency header is inherited from `Hub_size_2D.m` and lists
  `pkfnd`/`cntrd`/`local_background_2Dg` rather than `detectPunctaDBSCAN`.

---

## 7. Draft text for Methods / Code Availability

> Nuclei were segmented with Cellpose. Hub counts and per-hub enrichment were
> obtained by density-based clustering (DBSCAN) of nuclear puncta, with
> intensities measured on the raw image and normalised to the mean intensity of
> the host nucleus. For hub sizing, tracking, and colocalization, hubs were
> detected as local intensity maxima exceeding three times the mean intensity of
> their host nucleus and localised to sub-pixel precision by centroid refinement.
> Hub size was quantified as the Feret diameter and as the FWHM of a
> background-subtracted 2D Gaussian fit. Structural lifetimes were obtained by
> linking hubs across frames after rigid-body registration of the nucleus, with
> lifetimes edge-censored and summarised by Kaplan–Meier survival analysis.
> Colocalization between BRD4 and MED14 hubs was classified by a 300 nm
> centroid-proximity cutoff and quantified by per-hub Manders coefficients
> computed against a two-fold nuclear-mean threshold in the partner channel. All
> scripts and the parameter values used for each dataset are available at
> [repository / DOI].

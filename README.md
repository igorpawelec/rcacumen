# rcacumen

<img src="https://raw.githubusercontent.com/igorpawelec/rcacumen/main/www/rcacumen.png" alt="rcacumen logo" align="right" width="200"/>

[![R-CMD-check](https://github.com/igorpawelec/rcacumen/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/igorpawelec/rcacumen/actions/workflows/R-CMD-check.yaml)
[![Release](https://img.shields.io/github/v/release/igorpawelec/rcacumen)](https://github.com/igorpawelec/rcacumen/releases)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-%3E%3D%203.6-blue.svg)](https://www.r-project.org)

**Individual tree crown delineation from canopy height models, by Hierarchical Region Growing.**

A canopy height model goes in; a crown per tree comes out, as a raster of ids and as polygons. R and plain C — no Rcpp, no compiled dependencies, and `terra` only if you want to read and write files.

> **Python users:** the same algorithm is in [pycacumen](https://github.com/igorpawelec/pycacumen). The two are separate repositories because their tooling and idioms do not mix — but they are not merely similar. On the shared synthetic test suite they produce **identical output, pixel for pixel**. On real canopy height models they agree closely but not exactly, because the two watershed implementations break plateau ties differently. Both numbers are measured, not claimed — see [Agreement with pycacumen](#agreement-with-pycacumen).

## The problem it solves

Every crown delineation that starts from a canopy height model meets the same obstacle: **tree tops are over-detected**. The upper surface of a broad crown is ragged, so a local-maximum detector finds several peaks on it. The usual remedy is to smooth harder or to raise the detection threshold, and both trade one error for another: the surplus peaks disappear together with the small, suppressed trees that were real.

rcacumen does not try to get the tops right first. It accepts the surplus and corrects it afterwards, in the growing itself:

1. **Watershed.** Every tree top seeds a marker-based watershed on the inverted CHM, so there is exactly one region per top — the regions *are* the detected trees.
2. **Region adjacency graph.** Neighbouring regions are joined by a weighted edge, *w(a,b) = α·|Δμ| + β·|Δσ| + γ/(border+1)*, and the region statistics come from the same single pixel pass.
3. **Growing.** Each tree grows greedily, absorbing neighbours while the combined height variance stays under a threshold. When two trees absorb each other, they were one crown — this is where over-detection is undone.
4. **Arbitration.** Growing runs independently per seed, so two crowns can claim the same region. The claim is settled on the data — the taller tree, the nearer tree, or the better height match — never on iteration order, so the result does not depend on how the tops were sorted.

<img src="https://raw.githubusercontent.com/igorpawelec/rcacumen/main/www/pipeline.png" alt="A 0.5 m canopy height model of a circular sample plot; the same tile with 253 local maxima of which 213 are kept as tree tops; and the 212 crowns delineated from them, drawn as outlines over the height model" width="100%"/>

*A 100 × 100 m tile of a 0.5 m canopy height model (`test_data/chm_150_2014.tif`, a mature stand on a circular sample plot). Left: the input. Middle: 253 local maxima after median smoothing (`hmin` 7 m, 5 px window), of which 213 survive merging within 5 px and screening below 10 m — the filled markers. Right: the 212 crowns grown from them at the default `variance_thresh` of 2. Made by `www/figures.R`, from this package.*

## The one parameter that matters

`variance_thresh` is the height variance (σ², in m²) a single crown may contain. It sets how readily a tree absorbs its neighbours, and therefore how many of the surplus tops are merged back. Everything else has a sensible default; this one you choose per stand.

<img src="https://raw.githubusercontent.com/igorpawelec/rcacumen/main/www/variance_thresh.png" alt="The same tile delineated at variance_thresh 2, 8 and 20: 212, 202 and 123 crowns; the tree tops absorbed into a neighbouring crown are drawn as hollow rings and cluster where crowns merged" width="100%"/>

*The same 213 tree tops grown at `variance_thresh` 2, 8 and 20. Filled markers are tops that ended up as a crown of their own; hollow rings are tops absorbed into a neighbour — 1, 11 and 90 of them. At 2 almost nothing merges and every peak on a broad crown stays a separate tree; at 20 the merging reaches across real crown boundaries in the dense south-east of the plot. The right value lies where the absorbed tops are the ragged-surface duplicates and not the suppressed neighbours, which is a judgement made by looking, and the point of drawing them. pycacumen's README shows the same tile from Python: identical at 2 and 20, one crown apart at 8 (203 against 202) — the watershed tie-breaking measured [below](#agreement-with-pycacumen), at work.*

Two further controls shape the outcome:

- **`conflict_rule`** decides who gets contested canopy: the taller tree (`"height"`, the default — dominant trees overtop their neighbours), the nearer seed (`"distance"`, classic ITC behaviour, splits rather than merges) or the better height match (`"similarity"`).
- **`protect_seeds = TRUE`** switches merging off entirely: every top keeps its crown. Use it when the tops are trusted, for instance field-measured stem positions.

The number of regions claimed by more than one crown is reported after every run (`attr(crowns, "n_contested")`, or `res$n_contested`); it grows with the threshold and is the quickest sign that a setting is too loose.

## When to use it, and when not

Use rcacumen when you have a canopy height model — from ALS, from photogrammetric point clouds, at any resolution around 0.25–1 m — and want one crown per tree with a defensible answer to over-detection. Bring your own tree tops if you have better ones; the growing works from any set of seeds.

Do not reach for it when you need Dalponte & Coomes delineation, or want to work on the point cloud directly: [lidR](https://github.com/r-lidar/lidR) and [itcSegment](https://cran.r-project.org/package=itcSegment) do that well in R, and rcacumen deliberately does not duplicate them.

The pipeline shape — smooth the CHM, find tops as local maxima, delineate crowns — comes from [PyCrown](https://github.com/manaakiwhenua/pycrown) (Zörner et al. 2018), via pycacumen; hierarchical region growing is the only delineation here.

### The package family

rcacumen is one step of a longer chain; the other steps are separate packages, each with an R and a Python twin.

| Step | R | Python |
|---|---|---|
| Colour-space conversion of orthophotos | [rgeopalette](https://github.com/igorpawelec/rgeopalette) | [pygeopalette](https://github.com/igorpawelec/pygeopalette) |
| Adaptive superpixels and seeded growing on orthophotos | [rgeoadaptels](https://github.com/igorpawelec/rgeoadaptels) | [pygeoadaptels](https://github.com/igorpawelec/pygeoadaptels) |
| Crowns from a canopy height model | **rcacumen** | [pycacumen](https://github.com/igorpawelec/pycacumen) |
| Standing dead trees on orthophotos | — | [pygeosnag](https://github.com/igorpawelec/pygeosnag) |
| The same, inside QGIS | | [qgis-geoadaptels-geopalette](https://github.com/igorpawelec/qgis-geoadaptels-geopalette), [qgis-geosnag](https://github.com/igorpawelec/qgis-geosnag) |
| Polish national geodata (GUGiK, BDL) | [rgeopl](https://github.com/igorpawelec/rgeopl) | — |

## Installation

```r
# install.packages("remotes")
remotes::install_github("igorpawelec/rcacumen")
```

A C compiler is needed, which R already requires on Linux and macOS; on Windows install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

`terra` is optional and only used to read and write rasters:

```r
install.packages("terra")
```

## Quick start

```r
library(rcacumen)

res <- delineate("chm.tif",
                 smooth_ws      = 3,
                 hmin           = 7,
                 detect_ws      = 5,
                 merge_distance = 5,
                 screen_hmin    = 10,
                 variance_thresh = 2)

res$crowns        # SpatRaster of crown ids, 0 = background
res$tops          # tree tops actually used
res$n_contested   # regions claimed by more than one crown

write_crowns(res$crowns, "crowns.tif")
plot(crowns_to_polygons(res$crowns))
```

### Matrices, without touching the disk

Every stage is a plain function, and none of the algorithm needs `terra`:

```r
sm   <- smooth_chm(chm, ws = 3, method = "median")
tops <- detect_tops(sm, hmin = 7, ws = 5)      # (n, 2) subpixel row/col
tops <- merge_tops(tops, distance = 5)

crowns <- grow_crowns(sm, as_pixels(tops), variance_thresh = 2)
attr(crowns, "n_contested")
```

## Reference

### Parameters

**Smoothing** — `smooth_chm(chm, ws, method)`

| Parameter | Default | Description |
|---|---|---|
| `ws` | 3 | Window size (px). Larger = fewer false tops, but merges close crowns |
| `method` | `"median"` | `median`, `mean`, `gaussian`, `maximum`. Median keeps crown edges sharp |

**Tree tops** — `detect_tops`, `merge_tops`, `screen_tops`

| Parameter | Default | Description |
|---|---|---|
| `hmin` | 2 | Minimum height (m) for a pixel to be a candidate |
| `ws` | 3 | Local-maximum window (px) = minimum spacing between tops |
| `distance` | 5 | Merge radius (px). Grouping is transitive — keep below crown diameter |

**Growing** — `grow_crowns`

| Parameter | Default | Description |
|---|---|---|
| `variance_thresh` | 2 | Max height variance within a crown. **The main lever** |
| `mask_thresh` | 0 | Minimum CHM height treated as canopy (m) |
| `morpho_radius` | 0 | Radius for opening/closing the mask. 0 = off |
| `alpha`, `beta`, `gamma` | 1, 0.5, 0.1 | Edge weights: mean diff, σ diff, inverse border length |
| `anneal_lambda` | 1 | Per-iteration tightening of the threshold. 1 = constant |
| `max_iters` | `NULL` | Cap on grow iterations per seed. `NULL` grows to natural termination |
| `conflict_rule` | `"height"` | Who wins contested canopy — see above |
| `protect_seeds` | `FALSE` | If TRUE, no tree is absorbed; every top yields a crown |
| `retry_rejected` | `FALSE` | Reconsider regions rejected earlier in the same grow |

Ties in the conflict rules resolve to the lower crown id, so output is fully reproducible.

### Agreement with pycacumen

Two implementations of one method are only useful if they agree, so this is measured, not asserted.

| Test | Result |
|---|---|
| Watershed vs `scikit-image`, 200 random scenes incl. plateaus | **0 of 422,687 px differ** |
| Full pipeline vs pycacumen, 70 synthetic scenes, all 3 conflict rules | 60 small scenes **exact**; 10 large crowded ones within **0.032 %** |
| Crown counts, contested counts | identical throughout |
| Watershed vs `scikit-image`, 7 **real** CHMs | 375 of 147,368 masked px differ (**0.25 %**) |
| Full pipeline vs pycacumen, same 7 real CHMs | 5,142 of 281,602 px differ (**1.8 %**) |

Reproduce it yourself:

```bash
pip install pycacumen
python3 tools/generate_pycacumen_reference.py
Rscript tools/cross_validate_against_pycacumen.R
```

<details>
<summary><b>Where the last 0.25 % goes</b></summary>

The synthetic scenes agree exactly; real canopy height models do not, and the honest version is that the equality holds on the test suite rather than in general. Measured stage by stage, each stage fed pycacumen's own output for the previous one:

- `smooth_chm` — identical to the last bit
- `detect_tops` — agrees to 2.8e-14
- the mask — identical
- the watershed — **375 of 147,368 masked pixels differ**

`src/watershed.c` reimplements `skimage.segmentation.watershed` rather than calling it, and on the plateau structure of a real CHM the two break ties differently. Which of the three documented details diverges — heap order, label-on-push, or neighbour order — is not yet established.

Those 375 pixels become 5,142 differing crown pixels, and on `chm_33_2012` 20 watershed pixels become 3,437: the watershed regions are the units the graph grows over, so moving one boundary pixel shifts a region's mean and variance, flips a merge, and carries a whole region to another crown. Treat any watershed difference as significant however small it looks.

Two details make that agreement possible, and both are easy to get wrong in a port:

- **R stores matrices by column, Python by row.** The neighbour visiting order decides how plateaus resolve, so the C code is handed `t(matrix)`, whose column-major buffer is byte-identical to a row-major NumPy array. Skip that and the watershed silently transposes its tie-breaking.
- **Floating-point addition is not associative.** The pixel pass traverses in the same order as pycacumen, so the same heights are summed in the same sequence and the variances match to the last bit.

</details>

<details>
<summary><b>Notes on behaviour</b></summary>

**Tree detection is unstable on plateaus, and the rest of the pipeline is not.** This is a property of the method, not of either implementation, and it matters more than the 0.25 % above. Adding 1 nanometre of noise to a CHM — far below any physical or instrumental meaning, and far below the single precision the data is stored in — changes nothing about the canopy, but it does break plateau ties differently. Measured over the seven real CHMs, eight repeats each (`tools/measure_tie_sensitivity.py`):

| | worst case |
|---|---|
| tree tops detected | **3 to 30 fewer**, some moving up to **13 px** |
| canopy pixels changing crown, tops free to move | up to **62 %** |
| canopy pixels changing crown, tops held fixed | **0 %**, on every scene |

Read that last row carefully: with the tree tops pinned, the watershed, the region growing and the conflict arbitration are perfectly stable. All of the instability enters through `detect_tops`. The tree count falls consistently rather than randomly, which suggests exact ties inflate it — a plateau of equal maxima yields tops that a hair of noise removes — though that mechanism is a hypothesis and has not been confirmed.

The practical consequence: **a tree count from this method carries an implementation-defined component**, and so do crown boundaries. Report the smoothing and detection parameters, and do not compare counts across implementations or library versions without re-running both.

**The crown count can be lower than the tree-top count.** That is the point: merged trees leave gaps in the id sequence. Use `protect_seeds = TRUE` if you need one crown per top.

**`retry_rejected` only bites in the middle.** A rejected region can become admissible later, because absorbing a large homogeneous neighbour can *lower* a crown's variance. On a 1105-tree synthetic scene this changed nothing at `variance_thresh` 2 or 8, and nothing at 120, but changed ~17% of the raster at 20. The direction is not predictable.

**Where the C is, and why.** Two routines are compiled: the watershed, which is a sequential priority flood that no amount of vectorising fixes, and the pixel pass, which touches every cell. Everything else — the graph, the growing, the arbitration — is ordinary R, because it runs once per *tree*, not once per pixel, and that is the part you might actually want to read or change.

</details>

<details>
<summary><b>Testing and the figures</b></summary>

```r
remotes::install_github("igorpawelec/rcacumen")
# or, from a clone:
#   R CMD build rcacumen && R CMD check rcacumen_0.5.0.tar.gz
```

The README figures are remade with `Rscript www/figures.R` from the tile in `test_data/`, using the package in the checkout (`pkgload`) rather than an installed copy.

</details>

## Citation

If you use rcacumen in your research, please cite the software and the work it builds on:

> Pawelec, I. (2026). *rcacumen: individual tree crown delineation from canopy height models by Hierarchical Region Growing* [R package]. https://github.com/igorpawelec/rcacumen

Machine-readable metadata, including everything cited below, is in [`CITATION.cff`](CITATION.cff).

**Upstream project.** The pipeline structure derives from PyCrown:

> Zörner, J., Dymond, J., Shepherd, J., Jolly, B. (2018). *PyCrown — Fast raster-based individual tree segmentation for LiDAR data.* Landcare Research NZ Ltd. https://doi.org/10.7931/M0SR-DN55
>
> Zörner, J., Dymond, J.R., Shepherd, J.D., Wiser, S.K., Bunting, P., Jolly, B. (2018). LiDAR-based regional inventory of tall trees — Wellington, New Zealand. *Forests* 9(11), 702. https://doi.org/10.3390/f9110702

**Methods used**

> Chan, T.F., Golub, G.H., LeVeque, R.J. (1983). Algorithms for computing the sample variance: analysis and recommendations. *The American Statistician* 37(3), 242–247. — the pairwise formula that merges two regions' statistics in O(1). Note this is *not* Welford's update, which adds one sample at a time; a crown absorbs a whole region at once.
>
> Beucher, S., Meyer, F. (1993). The morphological approach to segmentation: the watershed transformation. In: *Mathematical Morphology in Image Processing*, 433–481.
>
> Popescu, S.C., Wynne, R.H. (2004). Seeing the trees in the forest. *Photogrammetric Engineering & Remote Sensing* 70(5), 589–604. — local-maxima tree detection on a smoothed CHM.

## License

GNU General Public License v3.0 — see [LICENSE](LICENSE).

rcacumen derives from PyCrown, which is published under GPLv3; the licence carries over.

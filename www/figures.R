# The README figures, generated from test_data so they can be remade.
#
#     Rscript www/figures.R           # from the repository root
#
# Writes www/pipeline.png (CHM, tree tops, crowns) and www/variance_thresh.png
# (the same scene at three settings of the main lever, with the tops the
# growing absorbed drawn hollow) and prints the numbers the README captions
# quote. Same tile, same parameters and same palette as pycacumen's
# www/figures.py, so the two READMEs show the same scene from either language.
# Uses the package in this checkout, not an installed copy.
suppressPackageStartupMessages({ library(terra) })
pkgload::load_all(".", quiet = TRUE)

tile   <- "test_data/chm_150_2014.tif"
blues  <- c("#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
            "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b")
ramp   <- colorRampPalette(blues)(256)
orange <- "#eb6834"; ink <- "#1f1f1f"; muted <- "#6b6b6b"
hmin <- 7; ws_tops <- 5L; merge_d <- 5; screen_h <- 10

r    <- read_chm(tile)
chm  <- r$chm; tmpl <- r$template          # NA outside the plot stays NA, as delineate() would have it
disp <- chm; disp[!is.finite(disp)] <- 0   # zeroed only for display
vmax <- as.numeric(quantile(disp, 0.995))
sm   <- smooth_chm(chm, ws = 3L, method = "median")
raw  <- detect_tops(sm, hmin = hmin, ws = ws_tops)
kept <- screen_tops(sm, merge_tops(raw, distance = merge_d), hmin = screen_h)
px   <- as_pixels(kept)

# pixel (row, col) -> map x, y at the pixel centre
res_m <- terra::res(tmpl)[1]; e <- terra::ext(tmpl)
to_xy <- function(rc) cbind(e$xmin + (rc[, 2] - 0.5) * res_m, e$ymax - (rc[, 1] - 0.5) * res_m)

absorbed <- function(crowns, px) {
  present <- logical(nrow(px)); present[unique(as.vector(crowns))[unique(as.vector(crowns)) > 0]] <- TRUE
  gone <- !present
  for (i in which(gone)) {                       # every absorbed top sits inside somebody else's crown
    id <- crowns[px[i, 1], px[i, 2]]
    stopifnot(id > 0, id != i)
  }
  gone
}

runs <- list()
for (vt in c(2, 8, 20)) {
  crowns <- grow_crowns(sm, px, variance_thresh = vt)
  gone <- absorbed(crowns, px)
  n <- length(setdiff(unique(as.vector(crowns)), 0L))
  runs[[as.character(vt)]] <- list(crowns = crowns, n = n, a = sum(gone), gone = gone)
}
cat(sprintf("tile %s: %d x %d px at %g m, heights to %.1f m\n", basename(tile), ncol(chm), nrow(chm), res_m, max(disp)))
cat(sprintf("tops: %d local maxima (hmin %g m, window %d px) -> %d after merging within %g px and screening below %g m\n",
            nrow(raw), hmin, ws_tops, nrow(kept), merge_d, screen_h))
for (vt in names(runs)) cat(sprintf("variance_thresh %4s: %d crowns from %d tops, %d tops absorbed\n",
                                    vt, runs[[vt]]$n, nrow(kept), runs[[vt]]$a))

panel <- function(m, label) {
  rr <- terra::rast(tmpl); terra::values(rr) <- pmin(pmax(as.numeric(t(m)), 0), vmax)   # clamped both ways, like imshow
  terra::plot(rr, col = ramp, range = c(0, vmax), legend = FALSE, axes = FALSE, mar = c(3.2, 0.4, 0.4, 0.4), maxcell = Inf)
  rect(e$xmin, e$ymin, e$xmax, e$ymax, border = "#d9d9d9")             # frame and label in map units: terra resets par on exit
  h <- e$ymax - e$ymin
  text((e$xmin + e$xmax) / 2, e$ymin - 0.05 * h, label, adj = c(0.5, 1), cex = 0.86, col = ink, xpd = NA)
}
draw_crowns <- function(crowns) {
  terra::plot(crowns_to_polygons(crowns, tmpl), border = orange, lwd = 1.1, add = TRUE)
}
draw_tops <- function(rc, cex, filled = TRUE) {
  if (nrow(rc) == 0) return(invisible())
  xy <- to_xy(rc)
  if (filled) points(xy, pch = 21, bg = ink, col = "white", cex = cex, lwd = 0.8)
  else points(xy, pch = 21, bg = NA, col = ink, cex = cex, lwd = 1.3)
}
colour_bar <- function() {
  par(mar = c(3.2, 0.6, 0.4, 2.6))
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, vmax))
  z <- seq(0, vmax, length.out = 256)
  rect(0, z[-256], 1, z[-1], col = ramp, border = NA)
  rect(0, 0, 1, vmax, border = "#d9d9d9")
  at <- pretty(c(0, vmax)); at <- at[at <= vmax]
  axis(4, at = at, las = 1, cex.axis = 0.7, col = "#d9d9d9", col.axis = muted, lwd = 0, lwd.ticks = 1, line = -0.2)
  mtext("height [m]", side = 4, line = 1.7, cex = 0.72, col = muted)
}

# ---- figure 1: the pipeline ----------------------------------------------------
png("www/pipeline.png", width = 11.4, height = 4.1, units = "in", res = 160, bg = "white")
layout(matrix(1:4, 1), widths = c(1, 1, 1, 0.12))
panel(disp, "canopy height model, 0.5 m")
panel(sm, sprintf("%d local maxima, %d kept as tree tops", nrow(raw), nrow(kept)))
draw_tops(as_pixels(raw), cex = 0.65, filled = FALSE); draw_tops(px, cex = 0.8)
panel(sm, sprintf("%d crowns", runs[["2"]]$n)); draw_crowns(runs[["2"]]$crowns); draw_tops(px, cex = 0.6)
colour_bar()
invisible(dev.off())

# ---- figure 2: the main lever --------------------------------------------------
png("www/variance_thresh.png", width = 11.4, height = 4.1, units = "in", res = 160, bg = "white")
layout(matrix(1:3, 1))
for (vt in names(runs)) {
  x <- runs[[vt]]
  panel(sm, sprintf("variance_thresh = %s\n%d crowns, %d top%s absorbed", vt, x$n, x$a, if (x$a == 1) "" else "s"))
  draw_crowns(x$crowns)
  draw_tops(px[!x$gone, , drop = FALSE], cex = 0.55)
  draw_tops(px[x$gone, , drop = FALSE], cex = 1.1, filled = FALSE)
}
invisible(dev.off())
cat("written: www/pipeline.png, www/variance_thresh.png\n")

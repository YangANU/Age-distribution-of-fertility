########
# plots
########

savefig("Fig_1a", width = 12, height = 12, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(AUS_fert_NSW_fts, xlab = "", ylab = "", ylim = c(0,9000), yaxt="none", main = "NSW", lwd = CURVE_LWD)
axis(2, seq(0, 9000, 1500))
add_ylab("Birth counts")
add_xlab("Age")
ny = ncol(as.matrix(AUS_fert_NSW_fts$y))
col_rainbow <- rainbow(min(1024, 1.25 * ny))
curve_label("1975", x0 = 18, y0 = 6000,
            x1 = 20, y1 = AUS_fert_NSW_fts$y[6,1],
            x_text = 0, y_text = 200,
            text_col = col_rainbow[1])
curve_label("1999", x0 = 20, y0 = 7500,
            x1 = 24, y1 = AUS_fert_NSW_fts$y[which(age == 24),ceiling(ny/2)],
            x_text = 0, y_text = 200,
            text_col = col_rainbow[ceiling(ny/2)])
curve_label("2024", x0 = 25, y0 = 8500,
            x1 = 29, y1 = AUS_fert_NSW_fts$y[which(age == 29),ny],
            x_text = 0, y_text = 200,
            text_col = col_rainbow[ny])
dev.off()

savefig("Fig_1b", width = 12, height = 12, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(age, NSW_den), xlab = "", ylab = "", ylim = c(0,0.09), yaxt="none", main = "NSW", lwd = CURVE_LWD)
axis(2, seq(0, 0.09, 0.015))
lines(age, NSW_den[,24], col = "black", lty = 2, lwd = CURVE_LWD)
add_ylab("Density")
curve_label("Medoid", x0 = 45, y0 = 0.065,
            x1 = age[20], y1 = NSW_den[20,24], x_text = 1.5, y_text = 0.001)
add_xlab("Age")
dev.off()

savefig("Fig_1c", width = 12, height = 12, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(AUS_fert_VIC_fts, xlab = "", ylab = "", ylim = c(0,7500), yaxt="none", main = "VIC", lwd = CURVE_LWD)
axis(2, seq(0,7500, 1500))
add_ylab("Birth counts")
add_xlab("Age")
dev.off()

savefig("Fig_1d", width = 12, height = 12, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(age, VIC_den), xlab = "", ylab = "", ylim = c(0,0.09), yaxt="none", main = "VIC", lwd = CURVE_LWD)
axis(2, seq(0, 0.09, 0.015))
lines(age, VIC_den[,24], col = "black", lty = 2, lwd = CURVE_LWD)
add_ylab("Density")
curve_label("Medoid", x0 = 45, y0 = 0.065,
            x1 = age[20], y1 = VIC_den[20,24], x_text = 1.5, y_text = 0.001)
add_xlab("Age")
dev.off()

savefig("Fig_1e", width = 12, height = 12, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(AUS_fert_QLD_fts, xlab = "", ylab = "", ylim = c(0,5000), yaxt="none", main = "QLD", lwd = CURVE_LWD)
axis(2, seq(0, 5000, 1000))
add_ylab("Birth counts")
add_xlab("Age")
dev.off()

savefig("Fig_1f", width = 12, height = 12, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(age, QLD_den), xlab = "", ylab = "", ylim = c(0,0.09), yaxt="none", main = "QLD", lwd = CURVE_LWD)
axis(2, seq(0, 0.09, 0.015))
lines(age, QLD_den[,24], col = "black", lty = 2, lwd = CURVE_LWD)
add_ylab("Density")
add_xlab("Age")
curve_label("Medoid", x0 = 45, y0 = 0.065,
            x1 = age[20], y1 = QLD_den[20,24], x_text = 1.5, y_text = 0.001)
dev.off()

# All eight jurisdictions are drawn, so the figure shows the full spread. The
# five discussed in the text -- the three most populous states plus the two
# extremes, NT and ACT -- are drawn in colour and labelled; the remaining
# jurisdictions are drawn in grey as context, which keeps the middle of the
# range from becoming an unreadable tangle. NT and ACT deliberately avoid
# colours 7 (yellow) and 8 (grey), which do not survive printing.
if (!exists("savefig")) {
  if (is.null(src) || !exists("savefig", envir = src, inherits = FALSE))
    stop("savefig() not found; it lives in .RData")
  savefig <- get("savefig", envir = src)
}

lab_state <- c("NSW", "VIC", "QLD", "NT", "ACT")   # placed first -> roomiest
bg_state  <- setdiff(state, lab_state)             # SA, WA, TAS
COL <- c(NSW = "black",     VIC = "red",       QLD = "green3",
         NT  = "blue",      ACT = "magenta",
         SA  = "darkorange", WA = "darkcyan",  TAS = "saddlebrown")
LTY <- c(NSW = 1, VIC = 2, QLD = 3, NT = 4, ACT = 5,
         SA  = 6, WA  = 1, TAS = 2)
ARROW <- 5                     # length of the label arrow, in years

all_lab <- c(lab_state, bg_state)
LAB_COL <- as.list(COL)

# Choose a label anchor whose text box clears every curve -- including the
# unlabelled grey ones -- and every label already placed, preferring the
# position with the largest clearance. Returns NULL if some label cannot be
# placed, so the caller can retry with a smaller margin.
place <- function(series, ylim, label_these, pad = 1, fixed = list()) {
  placed <- list()
  for (nm in names(fixed)) {
    f <- fixed[[nm]]
    placed[[nm]] <- list(
      x0 = f$cx + 1.5, cx = f$cx, y0 = f$y0, tip = f$tip,
      val = approx(year, series[[nm]], xout = f$tip)$y,
      hw = strwidth(nm) / 2 + 0.6 * pad,
      hh = strheight(nm) / 2 + 0.012 * diff(ylim) * pad,
      clear = NA_real_)
  }
  label_these <- setdiff(label_these, names(fixed))
  for (nm in label_these) {
    v <- series[[nm]]; best <- NULL
    for (tip in seq(1982, 2020, by = 1)) {
      val <- approx(year, v, xout = tip)$y
      for (sgn in c(1, -1)) for (mult in c(1, 1.5, 2, 2.6)) {
        y0 <- val + sgn * mult * 0.045 * diff(ylim)
        if (y0 < ylim[1] + 0.05 * diff(ylim) ||
            y0 > ylim[2] - 0.05 * diff(ylim)) next
        x0 <- tip - ARROW
        cx <- x0 - 1.5                                  # text centre
        hw <- strwidth(nm) / 2 + 0.6 * pad
        hh <- strheight(nm) / 2 + 0.012 * diff(ylim) * pad
        if (cx - hw < 1976) next
        xs <- seq(cx - hw, cx + hw, length.out = 9)
        clear <- Inf
        for (s2 in names(series))
          clear <- min(clear,
                       min(abs(approx(year, series[[s2]], xout = xs)$y - y0)) - hh)
        if (clear <= 0) next
        bad <- FALSE
        for (p in placed)
          if (abs(p$cx - cx) < (hw + p$hw) * 1.1 &&
              abs(p$y0 - y0) < (hh + p$hh) * 1.6) bad <- TRUE
        if (bad) next
        if (is.null(best) || clear > best$clear)
          best <- list(x0 = x0, cx = cx, y0 = y0, tip = tip, val = val,
                       hw = hw, hh = hh, clear = clear)
      }
    }
    if (is.null(best)) return(NULL)
    placed[[nm]] <- best
  }
  placed
}

panel <- function(file, M, ylim, yat, ylab, h0 = NULL, fixed = list()) {
  series <- lapply(state, function(s) as.numeric(M[s, ])); names(series) <- state
  savefig(file, width = 12, height = 12, toplines = 0.8, type = "png",
          pointsize = 12)
  std_par()
  plot(year, series[[bg_state[1]]], xlab = "", ylab = "", xaxt = "n",
       yaxt = "n", type = "l", lwd = CURVE_LWD, ylim = ylim,
       col = COL[[bg_state[1]]], lty = LTY[[bg_state[1]]])
  axis(side = 1, at = seq(1975, 2025, 10))
  axis(side = 2, at = yat)
  # reference line where the sign of the skewness changes
  if (!is.null(h0)) abline(h = h0, lty = 3, lwd = AXIS_LWD, col = "grey50")
  for (s in c(bg_state[-1], lab_state))
    lines(year, series[[s]], col = COL[[s]], lty = LTY[[s]], lwd = CURVE_LWD)
  add_ylab(ylab)
  add_xlab("Year")
  pos <- NULL
  for (pad in c(1, 0.8, 0.6, 0.45, 0.3, 0.2)) {
    pos <- place(series, ylim, all_lab, pad, fixed)
    if (!is.null(pos)) break
  }
  if (is.null(pos)) stop("could not place all labels in ", file)
  cat(sprintf("  %s: labels placed, pad = %.2f, min clearance = %.4f%s\n",
              file, pad, min(sapply(pos, function(p) p$clear), na.rm = TRUE),
              if (length(fixed)) paste0(" (pinned: ", paste(names(fixed), collapse = ", "), ")") else ""))
  for (nm in names(pos)) {
    p <- pos[[nm]]
    # note: text_col, not col -- curve_label() forwards ... into text(),
    # which already receives col = text_col
    curve_label(nm, x0 = p$x0, y0 = p$y0, x_text = -1.5,
                x1 = p$tip, y1 = p$val, text_col = LAB_COL[[nm]])
  }
  dev.off()
  invisible(pos)
}


panel("Fig_2a", Gini, c(0.385, 0.605), seq(0.40, 0.60, 0.04), "Gini coefficient")
panel("Fig_2b", SDv,  c(4.5, 6.5),     seq(4.5, 6.5, 0.4),    "Standard deviation (years)")
panel("Fig_2c", IQRv, c(5.7, 10.15),   seq(6.0, 10.0, 1.0),   "Interquartile range (years)")
panel("Fig_skewness", Skew, c(-0.30, 0.78),  seq(-0.2, 0.6, 0.2),   "Skewness", h0 = 0)

savefig("Fig_cont_mod", width = 24, height = 12, toplines = 0.8, type = "png", pointsize = 12)
std_par()
par(mfrow = c(2, 4), mar = c(4.4, 3.4, 3.4, 0.8), mgp = c(2.1, 0.7, 0))
for (s in state) {
  plot(year, M_continuous[s, ], type = "n", xlab = "", ylab = "", yaxt="none", 
       cex.lab = 1.5)
  title(s, cex.main = 1.5, line = 1)
  add_ylab("Modal age at birth")
  add_xlab("Year")
  axis(2, seq(20, 36, 5))
  polygon(c(year, rev(year)), c(CI[1, s, ], rev(CI[2, s, ])),
          col = "grey80", border = NA)
  lines(year, M_continuous[s, ], lwd = 2)
  legend("bottomleft", bty = "n", cex = 1.4,
         legend = sprintf("CI width %.2f yr", mean(W[s, ])))
}
dev.off()

savefig("Fig_3a", width = 12, height = 12, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, AUS_den), xlab = "", ylab = "", col = "gray",
     ylim = c(0, 0.1), main = "AUS", lwd = CURVE_LWD)
lines(fts(ages, t(AUS_den_fore)), lwd = CURVE_LWD)
add_ylab("Age distribution\nof birth counts")
add_xlab("Age")
curve_label("Historical\ndata", y_text = 0.003, x0 = 20, y0 = 0.077,
            x1 = 25, y1 = 0.050, text_col = "gray")
curve_label("10-steps-ahead\nforecasts", x0 = 42, y0 = 0.077, y_text = 0.003,
            x_text = 0.5, x1 = 36, y1 = AUS_den_fore[5,22], text_col = "red")
dev.off()

savefig("Fig_3b", width = 12, height = 12, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, NSW_den), xlab = "", ylab = "", col = "gray",
     ylim = c(0, 0.1), main = "NSW", lwd = CURVE_LWD)
lines(fts(ages, t(NSW_den_fore)), lwd = CURVE_LWD)
add_ylab("Age distribution\nof birth counts")
add_xlab("Age")
# curve_label("Historical\ndata",   y_text = 0.003,   x0 = 20, y0 = 0.077, x1 = 25, y1 = 0.050, col = "gray")
# curve_label("10-step\nforecasts",   y_text = 0.003,  x0 = 42, y0 = 0.077, x1 = 36, y1 = NSW_den_fore[5,22], col = "red")
dev.off()

savefig("Fig_3c", width = 12, height = 12, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, VIC_den), xlab = "", ylab = "", col = "gray",
     ylim = c(0, 0.1), main = "VIC", lwd = CURVE_LWD)
lines(fts(ages, t(VIC_den_fore)), lwd = CURVE_LWD)
add_ylab("Age distribution\nof birth counts")
add_xlab("Age")
# curve_label("Historical\ndata",    x0 = 20, y0 = 0.077, x1 = 25, y1 = 0.050, col = "gray")
# curve_label("10-step\nforecasts",  x0 = 42, y0 = 0.077, x1 = 33, y1 = 0.055, col = "red")
dev.off()

savefig("Fig_3d", width = 12, height = 12, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, QLD_den), xlab = "", ylab = "", col = "gray",
     ylim = c(0, 0.1), main = "QLD", lwd = CURVE_LWD)
lines(fts(ages, t(QLD_den_fore)), lwd = CURVE_LWD)
add_ylab("Age distribution\nof birth counts")
add_xlab("Age")
# curve_label("Historical\ndata",    x0 = 20, y0 = 0.077, x1 = 25, y1 = 0.050, col = "gray")
# curve_label("10-step\nforecasts",  x0 = 42, y0 = 0.077, x1 = 33, y1 = 0.055, col = "red")
dev.off()

























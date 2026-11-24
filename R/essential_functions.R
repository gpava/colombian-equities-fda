# Draw and validate EPS figures from the fitted objects in memory.
font_publication <- "Helvetica"

jet_publication <- c("#00007F", "#0000FF", "#007FFF", "#00FFFF", "#7FFF7F", "#FFFF00", "#FF7F00", "#FF0000", 
    "#7F0000")

pal_beta <- function(n = 256L) {
    stopifnot(length(n) == 1L, is.finite(n), n >= 1, n == as.integer(n))
    (grDevices::colorRampPalette(jet_publication))(n)
}

pal_map <- function(n = 256L) pal_beta(n)

colours_map <- function(x, lim) coefficient_colours(x, lim, pal_beta())

pal_correlation <- function(n = 201L) (grDevices::colorRampPalette(c("#2166AC", "#67A9CF", "#D1E5F0", 
    "#F7F7F7", "#FDDBC7", "#EF8A62", "#B2182B")))(n)

limits_color <- function(B) {
    r <- range(B, finite = TRUE)
    if (length(r) != 2L || !all(is.finite(r))) 
        stop("The plot has no finite values.")
    if (r[1] == r[2]) 
        return(r + c(-1, 1) * max(abs(r[1]), 1) * 1e-08)
    r
}

pal_count <- function() pal_beta(14L)

col_companies <- function(n = 13L) grDevices::hcl.colors(n, "Dark 3")

sci <- function(x) ifelse(x == 0, "0", formatC(x, format = "e", digits = 1))

base_style <- function(mar = c(4.2, 4.5, 1.4, 1.2), ps = 10) {
    par(mar = mar, mgp = c(2.6, 0.65, 0), tcl = -0.22, las = 1, family = font_publication, ps = ps, cex = 1, 
        cex.axis = 0.85, cex.lab = 1, lwd = 0.8, xaxs = "i", yaxs = "r", fg = "#202A33", col.axis = "#202A33", 
        col.lab = "#202A33", bg = "white")
}

image_eps <- function(x, y, z, ..., useRaster = TRUE) {
    edges <- function(v, n) {
        if (length(v) == n + 1L) 
            return(v)
        stopifnot(length(v) == n, n > 1L, all(diff(v) > 0))
        d <- diff(v)/2
        c(v[1] - d[1], v[-length(v)] + d, v[length(v)] + tail(d, 1))
    }
    ex <- edges(x, nrow(z))
    ey <- edges(y, ncol(z))
    stopifnot(max(abs(diff(ex) - mean(diff(ex)))) < 1e-07 * max(1, diff(range(ex))), max(abs(diff(ey) - 
        mean(diff(ey)))) < 1e-07 * max(1, diff(range(ey))))
    k <- pmax(1L, ceiling(par("pin") * 600/dim(z)))
    zz <- z[rep(seq_len(nrow(z)), each = k[1]), rep(seq_len(ncol(z)), each = k[2]), drop = FALSE]
    xx <- seq(ex[1], tail(ex, 1), length.out = nrow(zz) + 1L)
    yy <- seq(ey[1], tail(ey, 1), length.out = ncol(zz) + 1L)
    graphics::image(xx, yy, zz, ..., useRaster = TRUE)
}

locate_ghostscript <- function() {
    candidates <- c(Sys.getenv("R_GSCMD"), Sys.which(c("gs", "gswin64c", "gswin32c")), "/opt/homebrew/bin/gs", 
        "/usr/local/bin/gs")
    if (.Platform$OS.type == "windows") 
        candidates <- c(candidates, Sys.glob("C:/Program Files/gs/gs*/bin/gswin64c.exe"))
    candidates <- unique(candidates[nzchar(candidates)])
    found <- candidates[file.exists(candidates)]
    if (!length(found)) 
        stop("EPS export requires Ghostscript. Install it or set R_GSCMD to its executable.")
    normalizePath(found[1], winslash = "/", mustWork = TRUE)
}

embed_fonts_eps <- function(input, output) {
    gs <- locate_ghostscript()
    args <- c("-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=eps2write", "-dEPSCrop", "-dAutoRotatePages=/None", 
        "-dEmbedAllFonts=true", "-dSubsetFonts=true", "-dAutoFilterColorImages=false", "-dColorImageFilter=/LZWEncode", 
        "-dDownsampleColorImages=false", "-dDownsampleGrayImages=false", "-dDownsampleMonoImages=false", 
        shQuote(paste0("-sOutputFile=", output)), "-c", shQuote("<</NeverEmbed [] /AlwaysEmbed [/StandardSymbolsPS /Helvetica /Helvetica-Bold /Helvetica-Oblique /Helvetica-BoldOblique]>> setdistillerparams"), 
        "-f", shQuote(input))
    log <- suppressWarnings(system2(gs, args, stdout = TRUE, stderr = TRUE))
    status <- attr(log, "status")
    if (!is.null(status) && status != 0L) 
        stop("Ghostscript: ", paste(log, collapse = "\n"))
    if (!file.exists(output) || file.info(output)$size < 100) 
        stop("Empty EPS file: ", output)
    h <- readLines(output, n = 20, warn = FALSE)
    if (!grepl("EPSF", h[1], fixed = TRUE)) 
        stop("The output has no EPSF header.")
    invisible(output)
}

with_device <- function(path, format = "eps", width_mm, height_mm, draw, dpi = 600L) {
    if (!identical(format, "eps") || tolower(tools::file_ext(path)) != "eps") 
        stop("Only EPS figures can be exported.")
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    raw <- tempfile("eps-fonts-", fileext = ".eps")
    temp <- tempfile("eps-final-", fileext = ".eps")
    id <- NA_integer_
    on.exit({
        if (!is.na(id) && id %in% dev.list()) try(dev.off(id), silent = TRUE)
        unlink(c(raw, temp))
    }, add = TRUE)
    grDevices::postscript(file = raw, onefile = FALSE, horizontal = FALSE, paper = "special", width = width_mm/25.4, 
        height = height_mm/25.4, pointsize = 10, family = font_publication, encoding = "ISOLatin1.enc", 
        bg = "white", colormodel = "srgb", title = basename(path))
    id <- dev.cur()
    old_option <- options(fda.current_format = "eps")
    on.exit(options(old_option), add = TRUE)
    base_style()
    draw()
    dev.off(id)
    id <- NA_integer_
    ps <- readLines(raw, warn = FALSE)
    ps <- gsub("/Symbol", "/StandardSymbolsPS", ps, fixed = TRUE, useBytes = TRUE)
    ps <- gsub("font Symbol", "font StandardSymbolsPS", ps, fixed = TRUE, useBytes = TRUE)
    writeLines(ps, raw, useBytes = TRUE)
    embed_fonts_eps(raw, temp)
    if (!file.copy(temp, path, overwrite = TRUE)) 
        stop("Could not save ", path)
    invisible(path)
}

raster_panel_eps <- function(draw, size, mar) {
    temp <- tempfile(fileext = ".png")
    raw <- tempfile(fileext = ".ps")
    child <- NA_integer_
    on.exit({
        if (!is.na(child) && child %in% dev.list()) dev.off(child)
        unlink(c(temp, raw))
    }, add = TRUE)
    grDevices::postscript(file = raw, onefile = FALSE, horizontal = FALSE, paper = "special", width = size[1], 
        height = size[2], pointsize = 10, family = font_publication, encoding = "ISOLatin1.enc", bg = "white", 
        colormodel = "srgb")
    child <- dev.cur()
    base_style(mar)
    draw()
    dev.off(child)
    child <- NA_integer_
    gs <- locate_ghostscript()
    args <- c("-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m", "-r600", "-dTextAlphaBits=4", 
        "-dGraphicsAlphaBits=4", "-dAutoRotatePages=/None", paste0("-g", ceiling(size[1] * 600), "x", 
            ceiling(size[2] * 600)), shQuote(paste0("-sOutputFile=", temp)), shQuote(raw))
    log <- suppressWarnings(system2(gs, args, stdout = TRUE, stderr = TRUE))
    status <- attr(log, "status")
    if (!is.null(status) && status != 0L) 
        stop("Ghostscript (3D panel): ", paste(log, collapse = "\n"))
    if (!file.exists(temp) || file.info(temp)$size < 100) 
        stop("The temporary 3D panel was not generated.")
    img <- png::readPNG(temp, native = TRUE)
    if (!identical(as.integer(dim(img)), as.integer(ceiling(rev(size) * 600)))) 
        stop("Incorrect dimensions for the 600 dpi 3D panel.")
    par(mar = c(0, 0, 0, 0))
    plot.new()
    plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
    rasterImage(img, 0, 0, 1, 1, interpolate = FALSE)
}

coefficient_colours <- function(x, lim, pal = pal_beta()) {
    stopifnot(all(is.finite(x)), length(lim) == 2L, lim[1] < lim[2], all(x >= lim[1]), all(x <= lim[2]))
    k <- findInterval(x, seq(lim[1], lim[2], length.out = length(pal) + 1L), rightmost.closed = TRUE, 
        all.inside = TRUE)
    pal[k]
}

abbreviate_labels <- function(x) {
    x <- gsub("Bolsa de Valores de Colombia", "BVC", x, fixed = TRUE)
    gsub("Corporación Financiera de Colombia", "Corficolombiana", x, fixed = TRUE)
}

company_names <- function(keys) {
    labels <- c(bancolombia = "Bancolombia", bogota = "Banco de Bogotá", bvc = "BVC", ceargos = "Cementos Argos", 
        celsia = "Celsia", ecopetrol = "Ecopetrol", corpcolombia = "Corficolombiana", bolivar = "Grupo Bolívar", 
        mineros = "Mineros", nutresa = "Nutresa", terpel = "Terpel", promigas = "Promigas", suramericana = "Suramericana")
    y <- unname(labels[as.character(keys)])
    y[is.na(y)] <- as.character(keys)[is.na(y)]
    y
}

date_axis <- function(dates, side = 1, n = 5, indices = TRUE) {
    ix <- unique(round(seq(1, length(dates), length.out = n)))
    axis(side, at = if (indices) 
        ix
    else as.numeric(dates[ix]), labels = paste(month.abb[as.integer(format(dates[ix], "%m"))], format(dates[ix], 
        "%Y"), sep = "\n"), cex.axis = 0.8, las = 1)
}

draw_colorbar <- function(lim, pal, label, discrete = FALSE) {
    base_style(c(4.2, 0.3, 3, 3.8))
    par(mgp = c(2, 0.55, 0))
    if (discrete) {
        image_eps(c(0, 1), 0:13, matrix(rep(0:13, each = 2), nrow = 2), col = pal, breaks = seq(-0.5, 
            13.5, 1), axes = FALSE, xlab = "", ylab = "")
        axis(4, at = 0:13, las = 1, cex.axis = 0.78)
        abline(h = seq(0.5, 12.5, 1), col = "#D9D9D9", lwd = 0.6)
    }
    else {
        v <- seq(lim[1], lim[2], length.out = length(pal))
        image_eps(c(0, 1), v, matrix(rep(v, each = 2), nrow = 2), col = pal, zlim = lim, axes = FALSE, 
            xlab = "", ylab = "", useRaster = FALSE)
        ticks <- if (lim[1] < 0 && lim[2] > 0) 
            c(lim[1], 0, lim[2])
        else c(lim[1], mean(lim), lim[2])
        axis(4, at = ticks, labels = sci(ticks), las = 1, cex.axis = 0.78)
    }
    if (!discrete && lim[1] < 0 && lim[2] > 0) 
        segments(par("usr")[1], 0, par("usr")[2], 0, col = "#111111", lwd = 0.9)
    box()
    mtext(label, side = 3, line = 0.3, cex = 0.9)
}

draw_continuous_map <- function(B, s, t, lim) {
    colours <- colours_map(as.vector(B), lim)
    palette <- unique(colours)
    index <- matrix(match(colours, palette), nrow = nrow(B), ncol = ncol(B))
    image_eps(s, t, index, col = palette, breaks = seq(0.5, length(palette) + 0.5, by = 1), axes = FALSE, 
        xlab = "", ylab = "", useRaster = FALSE)
}

draw_map <- function(B, s, t, dates = NULL, events = NULL, count = FALSE, active = NULL, label = expression(hat(beta)(t, 
    s)), diagonal = FALSE, xlabel = "s (log-volume session)", ylabel = "t (log-return session)") {
    layout(matrix(c(1, 2), 1), widths = c(7.3, 1.25))
    base_style(c(4.7, 4.7, if (is.null(events)) 1.6 else 3, 1))
    pal <- if (count) 
        pal_count()
    else pal_map()
    lim <- if (count) 
        c(0, 13)
    else limits_color(B)
    if (count) 
        image_eps(s, t, B, col = pal, breaks = seq(-0.5, 13.5, 1), axes = FALSE, xlab = "", ylab = "", 
            useRaster = FALSE)
    else draw_continuous_map(B, s, t, lim)
    if (!count) {
        lev <- pretty(lim, n = 8)
        if (any(lev < 0)) 
            contour(s, t, B, levels = lev[lev < 0], add = TRUE, drawlabels = FALSE, col = "#666666", 
                lwd = 0.6, lty = 2)
        if (any(lev > 0)) 
            contour(s, t, B, levels = lev[lev > 0], add = TRUE, drawlabels = FALSE, col = "#666666", 
                lwd = 0.6, lty = 1)
        if (lim[1] < 0 && lim[2] > 0) 
            contour(s, t, B, levels = 0, add = TRUE, drawlabels = FALSE, col = "#111111", lwd = 1)
    }
    if (!is.null(dates)) {
        date_axis(dates, 1, 5)
        date_axis(dates, 2, 5)
    }
    else {
        axis(1, cex.axis = 0.8)
        axis(2, cex.axis = 0.8)
    }
    mtext(xlabel, 1, 3.2)
    mtext(ylabel, 2, 3.5, las = 0)
    if (diagonal) 
        abline(0, 1, lty = 2, col = "#444444", lwd = 0.6)
    if (!is.null(events)) {
        abline(v = events$Index_s, col = "#555555", lty = 3, lwd = 0.55)
        for (i in seq_len(nrow(events))) mtext(paste0("E", i), side = 3, at = events$Index_s[i], line = if (i%%2) 
            0.35
        else 1.35, cex = 0.76)
    }
    if (!is.null(active)) 
        contour(s, t, active, levels = 0.05 * max(active), add = TRUE, drawlabels = FALSE, col = "#111111", 
            lwd = 0.8, lty = 2)
    box()
    draw_colorbar(lim, pal, label, count)
}

draw_surface <- function(B, s, t, label = expression(hat(beta)(t, s)), xlabel = "s (log-volume session)", 
    ylabel = "t (log-return session)") {
    layout(matrix(c(1, 2), 1), widths = c(7.4, 1.3))
    base_style(c(2.5, 4.6, 1.8, 0.5))
    pal <- pal_beta()
    lim <- limits_color(B)
    f <- (B[-1, -1] + B[-nrow(B), -1] + B[-1, -ncol(B)] + B[-nrow(B), -ncol(B)])/4
    cols <- coefficient_colours(f, lim, pal)
    render_panel <- function() {
        persp(s, t, B, theta = 35, phi = 26, expand = 0.6, scale = TRUE, zlim = lim, col = cols, border = NA, 
            shade = NA, ticktype = "detailed", nticks = 3, xlab = xlabel, ylab = ylabel, zlab = "", cex.axis = 0.75, 
            cex.lab = 0.85, ltheta = -135)
    }
    if (identical(getOption("fda.current_format"), "eps")) {
        size <- dev.size("in")
        size[1] <- size[1] * 7.4/8.7
        raster_panel_eps(render_panel, size, c(2.5, 4.6, 1.8, 0.5))
    }
    else render_panel()
    draw_colorbar(lim, pal, label)
}

draw_event_curve <- function(e, what = "mean", panel = NULL) {
    v <- if (what == "mean") 
        e$summary_beta_s$Beta_mean
    else e$summary_beta_s$Beta_RMS
    base_style(c(4.6, 5.8, 2.6, 1))
    lim <- range(c(v, 0))
    pad <- diff(lim) * 0.1
    plot(seq_along(v), v, type = "l", lwd = 1.2, xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = lim + 
        c(-pad, pad))
    ticks <- pretty(lim, n = 4)
    axis(2, at = ticks, labels = sci(ticks), cex.axis = 0.8)
    date_axis(e$dates_s, 1, 5)
    mtext("s (log-volume session)", 1, 3.2)
    mtext(if (what == "mean") 
        expression(bar(beta)(s))
    else expression(RMS[beta](s)), 2, 3.5, las = 0)
    if (what == "mean") 
        abline(h = 0, lty = 2, col = "#777777", lwd = 0.6)
    idx <- e$events_plot$Index_s
    abline(v = idx, col = "#AAAAAA", lty = 3, lwd = 0.6)
    points(idx, v[idx], pch = 21, bg = "white", col = "#202A33", cex = 0.75)
    for (i in seq_along(idx)) mtext(paste0("E", i), 3, at = idx[i], line = if (i%%2) 
        0.25
    else 1.25, cex = 0.76)
    if (!is.null(panel)) 
        mtext(panel, 3, adj = 0, line = 1.4, font = 2)
}

draw_r2 <- function(x, y, dates = NULL) {
    base_style()
    plot(x, y, type = "l", lwd = 1.15, xlab = "", ylab = expression(R^2(t)), ylim = c(min(y, 0, na.rm = TRUE), 
        1), xaxt = if (is.null(dates)) 
        "s"
    else "n")
    if (!is.null(dates)) 
        date_axis(dates, 1, 5, indices = FALSE)
    mtext(if (is.null(dates)) 
        "t (trading-session index)"
    else "Date of log-return", 1, 3)
    abline(h = 0, col = "#999999", lwd = 0.6)
    abline(h = mean(y, na.rm = TRUE), lty = 2, col = "#2166AC", lwd = 1)
    legend("bottomright", legend = sprintf("Temporal mean = %.3f", mean(y, na.rm = TRUE)), lty = 2, col = "#2166AC", 
        bty = "n", cex = 0.85)
}

draw_extremes <- function(e, type) {
    d <- if (type == "return") 
        e$returns_long
    else e$volume_observed_long
    chosen <- switch(type, return = e$table_returns_extremes, high = e$table_volume_high, low = e$table_volume_low)
    y <- if (type == "return") 
        "Return"
    else "Volume_log"
    if (type != "return") 
        d$Volume_log <- log1p(d$Volume)
    stopifnot(y %in% names(d), y %in% names(chosen))
    groups <- unique(as.character(d$Company))
    pal <- col_companies(length(groups))
    layout(matrix(c(1, 2), ncol = 1), heights = c(5, 1.35))
    base_style(c(3.4, 4.4, 0.8, 1))
    plot(as.numeric(d$date), d[[y]], type = "n", xaxt = "n", xlab = "", ylab = if (type == "return") 
        "Log-return"
    else "log(1 + observed trading volume)")
    points(as.numeric(d$date), d[[y]], col = "#C5CAD0", pch = 16, cex = 0.33)
    for (i in seq_along(groups)) {
        ix <- as.character(chosen$Company) == groups[i]
        points(as.numeric(chosen$date[ix]), chosen[[y]][ix], pch = 21, bg = pal[i], col = "#202A33", 
            cex = 0.9)
    }
    date_axis(sort(unique(d$date)), 1, 5, indices = FALSE)
    base_style(c(0, 0, 0, 0))
    plot.new()
    legend("center", abbreviate_labels(groups), col = "#202A33", pt.bg = pal, pch = 21, ncol = 3, cex = 0.75, 
        bty = "n", x.intersp = 0.5, y.intersp = 1.1)
}

evaluate_fd_stored <- function(fd, x) {
    b <- fd$basis
    stopifnot(identical(b$type, "bspline"), is.null(b$dropind), all(is.finite(x)), all(x >= b$rangeval[1]), 
        all(x <= b$rangeval[2]))
    ord <- b$nbasis - length(b$params)
    knots <- c(rep(b$rangeval[1], ord), b$params, rep(b$rangeval[2], ord))
    des <- splines::splineDesign(knots, x, ord = ord)
    stopifnot(ncol(des) == nrow(fd$coefs))
    des %*% fd$coefs
}

validate_evaluation_fd <- function(e) {
    pairs <- list(list(e$fd_returns, e$time_return, e$matrix_return_smoothed), list(e$fd_volume, 
        e$time_volume, e$matrix_volume_smoothed), list(e$fd_error_model, e$grid_diag, e$matrix_residuals_diag))
    for (p in pairs) stopifnot(isTRUE(all.equal(as.numeric(evaluate_fd_stored(p[[1]], p[[2]])), as.numeric(p[[3]]), 
        tolerance = 1e-12)))
    stopifnot(isTRUE(all.equal(as.numeric(evaluate_fd_stored(e$mean_returns, e$time_return)), 
        rowMeans(e$matrix_return_smoothed), tolerance = 1e-12)), isTRUE(all.equal(as.numeric(evaluate_fd_stored(e$mean_volume, 
        e$time_volume)), rowMeans(e$matrix_volume_smoothed), tolerance = 1e-12)))
    invisible(TRUE)
}

draw_smoothing_panels <- function(e, volume = FALSE) {
    x <- if (volume) 
        e$time_volume
    else e$time_return
    observed <- if (volume) 
        e$matrix_volume
    else e$matrix_returns
    fitted <- if (volume) 
        e$matrix_volume_smoothed
    else e$matrix_return_smoothed
    indices <- if (volume) 
        e$companies_volume_plot
    else e$indices_return_plot
    par(mfrow = c(length(indices), 2))
    for (j in indices) for (resid in c(FALSE, TRUE)) {
        base_style(c(4, 4.7, 2.4, 0.8))
        y <- if (resid) 
            observed[, j] - fitted[, j]
        else observed[, j]
        plot(x, y, pch = 16, cex = 0.4, col = "#A8B0B8", xlab = "Trading-session index", ylab = if (resid) 
            "Residual"
        else if (volume) 
            "log(1 + volume)"
        else "Log-return")
        if (resid) 
            abline(h = 0, lty = 2)
        else lines(x, fitted[, j], col = "#2166AC", lwd = 1.1)
        mtext(paste(company_names(colnames(observed)[j]), if (resid) 
            "residuals"
        else "smoothing"), 3, line = 0.5, cex = 0.8)
    }
}

draw_loco_compact <- function(loco, poster = FALSE) {
    fields <- c("Company_omitted", "Distance_L2_relative", "Correlation_surfaces", "Agreement_sign_region_active_5pct")
    m <- loco$metrics
    stopifnot(is.data.frame(m), all(fields %in% names(m)), nrow(m) > 1L, !anyNA(m[, fields]), !anyDuplicated(m$Company_omitted), 
        all(vapply(m[, fields[-1]], is.numeric, logical(1))), all(is.finite(as.matrix(m[, fields[-1]]))), 
        all(m$Distance_L2_relative >= 0), all(abs(m$Correlation_surfaces) <= 1), all(m$Agreement_sign_region_active_5pct >= 
            0 & m$Agreement_sign_region_active_5pct <= 1))
    m <- m[order(-m$Distance_L2_relative, as.character(m$Company_omitted)), , drop = FALSE]
    vals <- list(m$Distance_L2_relative, m$Correlation_surfaces, 100 * m$Agreement_sign_region_active_5pct)
    labels <- list(sprintf("%.3f", vals[[1]]), sprintf("%.3f", vals[[2]]), sprintf("%.1f", vals[[3]]))
    layout(matrix(1L))
    size <- par("din") * 25.4
    scale <- if (poster) 
        min(size/c(143.4, 103))
    else 1
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)
    par(mai = rep(0, 4), family = font_publication, ps = 10 * scale, cex = 1, xaxs = "i", yaxs = "i", 
        bg = "white", col = "black", xpd = FALSE)
    plot.new()
    plot.window(xlim = c(0, size[1]), ylim = c(0, size[2]))
    mm <- function(x) x * scale
    names <- company_names(m$Company_omitted)
    label_cex <- 0.85
    left <- mm(2)
    right <- size[1] - mm(2)
    names_end <- left + max(strwidth(names, cex = label_cex))
    first <- names_end + mm(4)
    gutter <- mm(4)
    span <- (right - first - 2 * gutter)/3
    number_width <- max(vapply(labels, function(x) max(strwidth(x, cex = label_cex)), numeric(1)))
    label_gap <- mm(2.4)
    data_width <- span - number_width - label_gap
    if (data_width < mm(18)) 
        stop("The LOCO plot needs more horizontal space. Use a width of at least 143.4 mm.")
    axis_y <- mm(15)
    row_bottom <- mm(21)
    row_top <- size[2] - mm(18)
    y <- seq(row_top, row_bottom, length.out = nrow(m))
    if (min(abs(diff(y))) < mm(4)) 
        stop("Increase the height of the LOCO plot to separate the company labels.")
    lim <- list(c(0, max(0.6, ceiling(max(vals[[1]]) * 10)/10)), c(min(0.85, floor(min(vals[[2]]) * 20)/20), 
        1), c(min(80, floor(min(vals[[3]])/5) * 5), 100))
    refs <- c(0, 1, 100)
    headers <- c("(a) Surface change", "(b) Surface correlation", "(c) Sign agreement")
    axis_labels <- list(expression(paste("Relative ", L^2, " distance")), "Pearson correlation", "Agreement (%)")
    text(left, size[2] - mm(7), "Omitted company", adj = c(0, 0.5), cex = label_cex)
    text(left, y, names, adj = c(0, 0.5), cex = label_cex)
    geometry <- vector("list", 3)
    for (k in 1:3) {
        x0 <- first + (k - 1) * (span + gutter)
        x1 <- x0 + data_width
        number_right <- x0 + span
        to_x <- function(v) x0 + (v - lim[[k]][1])/diff(lim[[k]]) * data_width
        tick_values <- if (k == 1 && lim[[k]][2] == 0.6) 
            c(0, 0.2, 0.4, 0.6)
        else if (k == 2 && lim[[k]][1] == 0.85) 
            c(0.85, 0.9, 0.95, 1)
        else if (k == 3 && lim[[k]][1] == 80) 
            c(80, 90, 100)
        else pretty(lim[[k]], n = 3)
        tick_values <- tick_values[tick_values >= lim[[k]][1] & tick_values <= lim[[k]][2]]
        tick_labels <- if (k == 1) 
            sprintf("%.1f", tick_values)
        else if (k == 2) 
            sprintf("%.2f", tick_values)
        else sprintf("%.0f", tick_values)
        segments(x0, y, x1, y, col = "#D8DEE4", lwd = 0.6)
        ref_x <- to_x(refs[k])
        segments(ref_x, axis_y + mm(1.7), ref_x, row_top + mm(2.1), col = "#667887", lty = 2, lwd = 0.7)
        points(to_x(vals[[k]]), y, pch = 21, bg = "#2166AC", col = "#163E5A", cex = 0.85, lwd = 0.6)
        text(number_right, y, labels[[k]], adj = c(1, 0.5), cex = label_cex)
        segments(x0, axis_y, x1, axis_y, lwd = 0.7)
        segments(to_x(tick_values), axis_y, to_x(tick_values), axis_y - mm(1.2), lwd = 0.7)
        text(to_x(tick_values), axis_y - mm(3.4), tick_labels, cex = 0.8)
        text((x0 + x1)/2, mm(5), axis_labels[[k]], cex = 0.85)
        text((x0 + number_right)/2, size[2] - mm(7), headers[k], cex = 0.85)
        if (k == 3) 
            text((x0 + number_right)/2, size[2] - mm(12), "Active region: 5%", cex = 0.8)
        label_left <- number_right - strwidth(labels[[k]], cex = label_cex)
        point_clearance <- label_left - x1
        line_clearance <- label_left - ref_x
        stopifnot(all(point_clearance >= label_gap - 1e-08), all(line_clearance >= label_gap - 1e-08))
        geometry[[k]] <- data.frame(Company = m$Company_omitted, Panel = letters[k], Value = vals[[k]], 
            Label = labels[[k]], Clearance_line_mm = line_clearance, Clearance_data_mm = point_clearance, 
            Text_nominal_pt = 10 * scale * label_cex)
    }
    invisible(do.call(rbind, geometry))
}

draw_loco_full_final <- function(loco) {
    m <- loco$metrics
    m <- m[order(m$Distance_L2_relative, decreasing = TRUE), ]
    y <- rev(seq_len(nrow(m)))
    v <- list(m$Distance_L2_relative, m$Correlation_surfaces, 100 * m$Agreement_sign_global, 100 * 
        m$Agreement_sign_region_active_5pct)
    par(mfrow = c(2, 2))
    for (k in 1:4) {
        base_style(c(3.3, 8.5, 3.5, 1), ps = 10)
        lim <- if (k == 1) 
            c(0, 0.6)
        else if (k == 2) 
            c(0.85, 1)
        else c(75, 100)
        stopifnot(all(v[[k]] >= lim[1]), all(v[[k]] <= lim[2]))
        plot(v[[k]], y, type = "n", axes = FALSE, xlim = lim + c(0, 0.03 * diff(lim)), ylim = c(0.5, 
            length(y) + 0.5), xlab = "", ylab = "")
        abline(h = y, col = "#DEE4E9", lwd = 0.7)
        points(v[[k]], y, pch = 21, bg = "#2166AC", cex = 0.9)
        axis(2, at = y, labels = company_names(m$Company_omitted), tick = FALSE, las = 1, cex.axis = 0.8)
        axis(1, at = if (k == 1) 
            c(0, 0.3, 0.6)
        else if (k == 2) 
            c(0.85, 0.925, 1)
        else c(75, 87.5, 100), cex.axis = 0.8)
        mtext(switch(k, expression(paste("(a) Relative ", L^2)), "(b) Surface", "(c) Global sign", "(d) Active sign"), 
            3, line = 2.1, cex = 0.9)
        mtext(c("distance", "correlation", "agreement (%)", "agreement (%)")[k], 3, line = 0.8, cex = 0.9)
        box()
    }
}

draw_sensitivity_final <- function(e) {
    ds <- e$df_sensitivity_plot
    lim <- limits_color(ds$beta)
    ticks <- function(d, x) {
        ix <- unique(round(seq(1, length(d), length.out = 3)))
        data.frame(coordinate = x[ix], label = paste(month.abb[as.integer(format(d[ix], "%m"))], format(d[ix], 
            "%Y"), sep = "\n"))
    }
    tx <- ticks(e$dates_s, e$time_volume)
    ty <- ticks(e$dates_t, e$time_return)
    layout(matrix(1:4, 2, byrow = TRUE))
    for (i in 1:3) {
        d <- ds[ds$scenario == e$scenarios_sensitivity$scenario[i], ]
        z <- matrix(d$beta, nrow = length(e$grid_s_sensitivity))
        stopifnot(nrow(d) == length(e$grid_s_sensitivity) * length(e$grid_t_sensitivity))
        base_style(c(4.4, 4.4, 3.8, 1.5), ps = 10)
        draw_continuous_map(z, e$grid_s_sensitivity, e$grid_t_sensitivity, lim)
        axis(1, at = tx$coordinate, labels = tx$label, cex.axis = 0.8)
        axis(2, at = ty$coordinate, labels = ty$label, cex.axis = 0.8, las = 1)
        lev <- pretty(lim, n = 8)
        contour(e$grid_s_sensitivity, e$grid_t_sensitivity, z, levels = lev[lev < 0], add = TRUE, drawlabels = FALSE, 
            col = "#666666", lwd = 0.6, lty = 2)
        contour(e$grid_s_sensitivity, e$grid_t_sensitivity, z, levels = lev[lev > 0], add = TRUE, drawlabels = FALSE, 
            col = "#666666", lwd = 0.6)
        contour(e$grid_s_sensitivity, e$grid_t_sensitivity, z, levels = 0, add = TRUE, drawlabels = FALSE, lwd = 1)
        mtext("s (volume session)", 1, line = 3, cex = 0.85)
        mtext("t (return session)", 2, line = 3, cex = 0.85, las = 0)
        mtext(c("(a) Less smoothing", "(b) Selected model", "(c) More smoothing")[i], 3, line = 2.3, 
            cex = 0.85)
        mtext(bquote(lambda[s] == 10^.(round(log10(e$scenarios_sensitivity$lambda_s[i]))) ~ "," ~ lambda[t] == 
            10^.(round(log10(e$scenarios_sensitivity$lambda_t[i])))), 3, line = 0.7, cex = 0.9)
        box()
    }
    base_style(c(4.4, 1.2, 3.8, 3), ps = 10)
    par(plt = c(0.18, 0.35, 0.22, 0.83))
    b <- seq(lim[1], lim[2], length.out = 256)
    image_eps(c(0, 1), b, matrix(rep(b, each = 2), nrow = 2), col = pal_beta(), axes = FALSE, xlab = "", 
        ylab = "", useRaster = FALSE, zlim = lim)
    at <- sort(unique(c(lim, 0, pretty(lim, n = 4))))
    at <- at[at >= lim[1] & at <= lim[2]]
    axis(4, at = at, labels = sci(at), cex.axis = 0.8, las = 1)
    segments(par("usr")[1], 0, par("usr")[2], 0, lwd = 1)
    box()
    mtext(expression(hat(beta)(t, s)), 3, line = 0.8, cex = 1)
}

draw_residual_fit_final <- function(e) {
    d <- e$df_res_fit
    base_style(c(4.2, 4.7, 1.3, 1), ps = 10)
    plot(d$Fitted, d$Residual, pch = 16, cex = 0.35, col = "#A6AFB7", xlab = "Fitted value", ylab = "Residual")
    abline(h = 0, lty = 2, lwd = 0.8, col = "#555555")
    fit <- stats::loess(Residual ~ Fitted, data = d, span = 0.75, degree = 2, family = "gaussian")
    x <- seq(min(d$Fitted), max(d$Fitted), length.out = 80)
    lines(x, predict(fit, newdata = data.frame(Fitted = x)), col = "#2166AC", lwd = 1.5)
}

draw_corr_final <- function(e) {
    z <- e$cor_residual_companies
    n <- nrow(z)
    base_style(c(8.5, 8.5, 1.5, 2), ps = 10)
    image_eps(seq_len(n), seq_len(n), z, col = pal_correlation(), zlim = c(-1, 1), axes = FALSE, xlab = "", 
        ylab = "", useRaster = FALSE)
    axis(1, at = seq_len(n), labels = company_names(colnames(z)), las = 2, tick = FALSE, cex.axis = 0.8)
    axis(2, at = seq_len(n), labels = company_names(rownames(z)), las = 1, tick = FALSE, cex.axis = 0.8)
    for (i in seq_len(n)) for (j in seq_len(n)) text(i, j, sprintf("%.1f", z[i, j]), cex = 0.65, col = if (abs(z[i, 
        j]) > 0.6) 
        "white"
    else "#202A33")
    box()
}

validate_numbers_final <- function(obj) {
    e <- obj$e
    l <- obj$loco
    B <- e$beta_final_mat
    weights <- function(x) {
        d <- diff(x)
        c(d[1]/2, (head(d, -1) + tail(d, -1))/2, tail(d, 1)/2)
    }
    w <- outer(weights(l$grid_s), weights(l$grid_t))
    norm <- sqrt(sum(w * B^2))
    active <- abs(B) >= 0.05 * max(abs(B))
    tol <- 1e-11
    equal <- function(a, b) stopifnot(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol)))
    for (i in seq_along(l$companies)) {
        z <- l$surfaces_loco[, , i]
        same <- sign(B) == sign(z)
        r <- l$metrics[i, ]
        stopifnot(identical(as.character(r$Company_omitted), as.character(l$companies[i])))
        equal(sqrt(sum(w * (z - B)^2))/norm, r$Distance_L2_relative)
        equal(cor(as.vector(B), as.vector(z)), r$Correlation_surfaces)
        equal(sum(w * same)/sum(w), r$Agreement_sign_global)
        equal(sum(w[active] * same[active])/sum(w[active]), r$Agreement_sign_region_active_5pct)
        equal(max(abs(z - B)), r$Max_difference_absolute)
    }
    ev <- e$table_events_beta
    ix <- ev$Index_s
    stopifnot(all(as.Date(e$dates_s[ix]) == as.Date(ev$date)), all(ev$date == ev$date_in_data))
    equal(rowMeans(B)[ix], ev$Beta_mean)
    equal(sqrt(rowMeans(B^2))[ix], ev$Beta_RMS)
    equal(rowMeans(B > 0)[ix], ev$Proportion_beta_positive)
    equal(diag(B)[ix], ev$Beta_contemporaneous)
    equal(rowMeans(B), e$summary_beta_s$Beta_mean)
    equal(sqrt(rowMeans(B^2)), e$summary_beta_s$Beta_RMS)
    stopifnot(length(e$res_finite) == 6370L, all(is.finite(e$res_finite)), all(is.finite(B)))
    invisible(TRUE)
}

verify_eps <- function(path, width_mm = NULL, height_mm = NULL) {
    all_lines <- readLines(path, warn = FALSE)
    h <- head(all_lines, 30)
    fonts <- grep("^<</BaseFont/", all_lines)
    for (i in fonts) if (!any(grepl("/FontDescriptor ", all_lines[i:min(i + 5, length(all_lines))], fixed = TRUE))) 
        stop("External font not embedded: ", path, " / ", all_lines[i])
    stopifnot(length(h) > 1, grepl("EPSF", h[1], fixed = TRUE), any(grepl("^%%Pages: 1$", h)))
    bb <- grep("^%%BoundingBox:", h, value = TRUE)[1]
    box <- as.numeric(strsplit(trimws(sub("^%%BoundingBox:", "", bb)), "[[:space:]]+")[[1]])
    stopifnot(length(box) == 4L, all(is.finite(box)), box[3] > box[1], box[4] > box[2])
    if (!is.null(width_mm)) 
        stopifnot(abs((box[3] - box[1]) - width_mm/25.4 * 72) < 2, abs((box[4] - box[2]) - height_mm/25.4 * 
            72) < 2)
    log <- suppressWarnings(system2(locate_ghostscript(), c("-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", 
        "-sDEVICE=nullpage", shQuote(path)), stdout = TRUE, stderr = TRUE))
    status <- attr(log, "status")
    if (!is.null(status) && status != 0L) 
        stop("Invalid EPS file: ", path, "\n", paste(log, collapse = "\n"))
    data.frame(File = path, Bytes = file.info(path)$size, Width_pt = box[3] - box[1], Height_pt = box[4] - 
        box[2], EPSF = TRUE, One_page = TRUE, Fonts_self_contained = TRUE, Ghostscript = TRUE)
}


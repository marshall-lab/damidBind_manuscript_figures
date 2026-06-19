#' Draw a ggplot2-compatible area-proportional Venn diagram
#'
#' Draws an approximately area-proportional Venn diagram for two or three sets
#' using `ggplot2` and `ggforce`. The function computes exact set-membership
#' counts, converts set sizes to circle areas, solves pairwise circle distances
#' analytically via numerical root finding, and returns a standard `ggplot`
#' object.
#'
#' For two circles, the pairwise area relationship can be represented exactly.
#' For three circles, exact representation of all seven Venn regions is not
#' always geometrically possible using circles. This function makes circle areas
#' proportional to set sizes and pairwise overlaps proportional where possible.
#' If the requested pairwise overlaps cannot form a valid triangle of circle
#' centres, distances are adjusted and a warning is issued.
#'
#' Labels show the exact set-region counts calculated from the input vectors,
#' regardless of any geometric approximation.
#'
#' Ported from and inspired by BioVenn by Tim Hulsen.
#'
#' @param list_x Character, numeric, or factor vector. Identifiers in set X.
#' @param list_y Character, numeric, or factor vector. Identifiers in set Y.
#' @param list_z Character, numeric, or factor vector. Identifiers in set Z.
#'   For a two-set Venn diagram, provide `character(0)` or `NULL` for one set.
#' @param title Character. Plot title.
#' @param subtitle Character or `NULL`. Plot subtitle.
#' @param xtitle Character. Label for set X.
#' @param ytitle Character. Label for set Y.
#' @param ztitle Character. Label for set Z.
#' @param nrtype Character. Type of region labels. One of `"abs"`, `"pct"`,
#'   or `"none"`.
#' @param palette Colour palette specification. May be:
#'   \itemize{
#'     \item `NULL`, equivalent to `"okabe-ito"`;
#'     \item a character vector of at least three colours;
#'     \item a palette function accepting `n`, e.g. `ggsci::pal_npg()`;
#'     \item one of `"okabe-ito"`, `"ggplot2"`, `"hue"`, `"npg"`, `"nejm"`,
#'       `"aaas"`, `"jama"`, or `"lancet"`.
#'   }
#' @param fill_alpha Numeric between 0 and 1. Fill transparency for circles.
#' @param stroke_colour Circle border colour. If `NULL` and `stroke_width > 0`,
#'   each circle uses its own set colour with alpha given by `stroke_alpha`.
#'   Use `NA` to suppress borders explicitly. May also be a single colour or
#'   a vector of colours for active circles.
#' @param stroke_alpha Numeric between 0 and 1. Alpha used for default
#'   set-coloured strokes when `stroke_colour = NULL`.
#' @param stroke_width Numeric. Circle border width. Set to `0` for no border.
#' @param label_colour Character. Colour of region-count labels.
#' @param label_size Numeric. Size of region-count labels, in ggplot2 text units.
#' @param label_fontface Character. Font face for region-count labels.
#' @param set_label_colour Character. Colour of set labels.
#' @param set_label_size Numeric. Size of set labels, in ggplot2 text units.
#' @param set_label_fontface Character. Font face for set labels.
#' @param show_set_labels Logical. If `TRUE`, set names are drawn.
#' @param repel_labels Logical. If `TRUE`, uses `ggrepel` to reposition region
#'   labels and set labels to reduce text overlap.
#' @param repel_seed Integer or `NULL`. Random seed passed to `ggrepel` for
#'   reproducible label placement.
#' @param repel_box_padding Numeric. Padding around label boxes for `ggrepel`.
#' @param repel_point_padding Numeric. Padding around label anchor points for
#'   `ggrepel`.
#' @param repel_force Numeric. Repulsion force used by `ggrepel`.
#' @param repel_force_pull Numeric. Pull force back towards original label
#'   positions used by `ggrepel`.
#' @param repel_max_iter Integer. Maximum number of iterations used by
#'   `ggrepel`.
#' @param repel_max_time Numeric. Maximum time, in seconds, used by `ggrepel`.
#' @param repel_segment_colour Character or `NA`. Colour of leader-line segments
#'   drawn by `ggrepel`. Use `NA` to suppress leader lines.
#' @param label_grid Integer. Resolution of the grid used to estimate region
#'   label centroids. Larger values improve label placement at some computational
#'   cost.
#' @param output Character. Output mode. One of `"plot"`, `"png"`, `"pdf"`,
#'   `"svg"`, `"jpg"`, or `"tif"`. `"plot"` returns the plot without saving.
#' @param filename Character or `NULL`. Output filename when `output != "plot"`.
#' @param width Numeric. Saved plot width in inches.
#' @param height Numeric. Saved plot height in inches.
#' @param dpi Numeric. Resolution for raster outputs.
#' @param bg Character. Background colour used by `ggsave()`.
#' @param set_label_position Character. Positioning strategy for set labels.
#'   `"outside_unique"` places each set label outside its circle in the direction
#'   closest to the set's unique region where possible. `"centre"` reproduces
#'   the previous centre-of-circle placement.
#' @param set_label_offset Numeric. Distance of outside set labels from the
#'   circle boundary, expressed as a fraction of the plot span.
#' @param set_count_nrtype Character. Optional count format appended to set
#'   labels. One of `"none"`, `"abs"`, `"pct"`, or `"abs_pct"`.
#' @param set_pct_denominator Character. Denominator used for set-label
#'   percentages. `"union"` uses the number of unique identifiers across all
#'   active sets. `"set_sum"` uses the sum of active set sizes, so percentages
#'   may sum to more than 100 when sets overlap.
#' @param return_lists Logical. If `TRUE`, returns a list containing the plot and
#'   all derived set lists. If `FALSE`, returns only the ggplot object, with
#'   derived lists attached as attributes.
#'
#' @return If `return_lists = FALSE`, a `ggplot` object. The object has
#'   attributes:
#'   \itemize{
#'     \item `"venn_lists"`: derived set lists and counts;
#'     \item `"circle_data"`: circle centres, radii, and colours;
#'     \item `"label_data"`: region label positions and text;
#'     \item `"plot_label_data"`: combined region and set label data before
#'       `ggrepel` adjustment.
#'   }
#'
#'   If `return_lists = TRUE`, a list with element `plot` plus the derived set
#'   lists and counts.
#'
#' @details
#' This function was shamelessly derived from BioVenn, ported over to a ggplot2 compatible
#' framework, and then some internals streamlined and optimised.  A ggrepel-based
#' labelling structure (optional) has been added, and there's compatibility with
#' ggsci colour schemes.
#'
#' All the credit goes to Tim Hulsen's original code.
#'
#' (Sorry, Tim, but I just couldn't get BioVenn to behave with patchwork ...)
#'
#' Region codes used internally are:
#' \itemize{
#'   \item `"100"`: X only;
#'   \item `"010"`: Y only;
#'   \item `"001"`: Z only;
#'   \item `"110"`: X and Y only;
#'   \item `"101"`: X and Z only;
#'   \item `"011"`: Y and Z only;
#'   \item `"111"`: X, Y, and Z.
#' }
#'
#' @examples
#' \dontrun{
#' list_x <- c("1007_s_at", "1053_at", "117_at", "121_at",
#'             "1255_g_at", "1294_at")
#' list_y <- c("1255_g_at", "1294_at", "1316_at",
#'             "1320_at", "1405_i_at")
#' list_z <- c("1007_s_at", "1405_i_at", "1255_g_at",
#'             "1431_at", "1438_at", "1487_at")
#'
#' p <- draw_venn_gg(
#'   list_x,
#'   list_y,
#'   list_z,
#'   title = "Example Venn",
#'   subtitle = "ggplot2-compatible proportional Venn",
#'   palette = "okabe-ito",
#'   stroke_width = 0.7
#' )
#'
#' p
#'
#' # ggsci palette by name
#' draw_venn_gg(list_x, list_y, list_z, palette = "npg")
#'
#' # ggsci palette function
#' draw_venn_gg(list_x, list_y, list_z, palette = ggsci::pal_jama())
#'
#' # Retrieve derived lists
#' attr(p, "venn_lists")$xy_only
#' attr(p, "venn_lists")$xyz
#'
#' # Save directly
#' draw_venn_gg(
#'   list_x,
#'   list_y,
#'   list_z,
#'   output = "png",
#'   filename = "venn.png",
#'   width = 7,
#'   height = 7,
#'   dpi = 300
#' )
#' }
#'
#' @importFrom stats uniroot aggregate
#' @export
draw_venn_gg <- function(
        list_x,
        list_y,
        list_z = NULL,
        title = NULL,
        subtitle = NULL,
        xtitle = "ID Set X",
        ytitle = "ID Set Y",
        ztitle = "ID Set Z",
        nrtype = c("abs", "pct", "abs_pct", "none"),
        palette = NULL,
        fill_alpha = 0.45,
        stroke_colour = NULL,
        stroke_alpha = 0.8,
        stroke_width = 0.7,
        label_colour = "black",
        label_size = 4,
        label_fontface = "plain",
        set_label_colour = "black",
        set_label_size = 4,
        set_label_fontface = "bold",
        show_set_labels = TRUE,
        repel_labels = TRUE,
        repel_seed = 1,
        repel_box_padding = 0.35,
        repel_point_padding = 0.15,
        repel_force = 1,
        repel_force_pull = 0.5,
        repel_max_iter = 10000,
        repel_max_time = 1,
        repel_segment_colour = NA,
        label_grid = 350,
        output = c("plot", "png", "pdf", "svg", "jpg", "tif"),
        filename = NULL,
        width = 7,
        height = 7,
        dpi = 300,
        bg = "white",
        set_label_position = c("outside_unique", "centre"),
        set_label_offset = 0.07,
        set_count_nrtype = c("none", "abs", "pct", "abs_pct"),
        set_pct_denominator = c("union", "set_sum"),
        return_lists = FALSE
) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.", call. = FALSE)
    }

    if (!requireNamespace("ggforce", quietly = TRUE)) {
        stop("Package 'ggforce' is required.", call. = FALSE)
    }

    if (!requireNamespace("scales", quietly = TRUE)) {
        stop("Package 'scales' is required.", call. = FALSE)
    }

    nrtype <- match.arg(nrtype)
    output <- match.arg(output)
    set_label_position <- match.arg(set_label_position)
    set_count_nrtype <- match.arg(set_count_nrtype)
    set_pct_denominator <- match.arg(set_pct_denominator)

    if (isTRUE(repel_labels) && !requireNamespace("ggrepel", quietly = TRUE)) {
        stop(
            "Package 'ggrepel' is required when `repel_labels = TRUE`.",
            call. = FALSE
        )
    }

    if (is.null(list_x)) list_x <- character(0)
    if (is.null(list_y)) list_y <- character(0)
    if (is.null(list_z)) list_z <- character(0)

    label_grid <- as.integer(label_grid)

    if (length(label_grid) != 1 || is.na(label_grid) || label_grid <= 0) {
        stop("`label_grid` must be a positive integer.", call. = FALSE)
    }

    if (label_grid < 50) {
        warning(
            "`label_grid` is quite low; label placement may be poor.",
            call. = FALSE
        )
    }

    if (!is.numeric(fill_alpha) || length(fill_alpha) != 1 ||
        is.na(fill_alpha) || fill_alpha < 0 || fill_alpha > 1) {
        stop("`fill_alpha` must be a numeric value between 0 and 1.", call. = FALSE)
    }

    if (!is.numeric(stroke_alpha) || length(stroke_alpha) != 1 ||
        is.na(stroke_alpha) || stroke_alpha < 0 || stroke_alpha > 1) {
        stop("`stroke_alpha` must be a numeric value between 0 and 1.", call. = FALSE)
    }

    if (!is.numeric(stroke_width) || length(stroke_width) != 1 ||
        is.na(stroke_width) || stroke_width < 0) {
        stop("`stroke_width` must be a non-negative numeric value.", call. = FALSE)
    }

    if (!is.logical(show_set_labels) || length(show_set_labels) != 1 ||
        is.na(show_set_labels)) {
        stop("`show_set_labels` must be TRUE or FALSE.", call. = FALSE)
    }

    if (!is.logical(repel_labels) || length(repel_labels) != 1 ||
        is.na(repel_labels)) {
        stop("`repel_labels` must be TRUE or FALSE.", call. = FALSE)
    }

    if (!is.numeric(repel_box_padding) || length(repel_box_padding) != 1 ||
        is.na(repel_box_padding) || repel_box_padding < 0) {
        stop("`repel_box_padding` must be a non-negative numeric value.", call. = FALSE)
    }

    if (!is.numeric(repel_point_padding) || length(repel_point_padding) != 1 ||
        is.na(repel_point_padding) || repel_point_padding < 0) {
        stop("`repel_point_padding` must be a non-negative numeric value.", call. = FALSE)
    }

    if (!is.numeric(repel_force) || length(repel_force) != 1 ||
        is.na(repel_force) || repel_force < 0) {
        stop("`repel_force` must be a non-negative numeric value.", call. = FALSE)
    }

    if (!is.numeric(repel_force_pull) || length(repel_force_pull) != 1 ||
        is.na(repel_force_pull) || repel_force_pull < 0) {
        stop("`repel_force_pull` must be a non-negative numeric value.", call. = FALSE)
    }

    if (!is.numeric(repel_max_time) || length(repel_max_time) != 1 ||
        is.na(repel_max_time) || repel_max_time <= 0) {
        stop("`repel_max_time` must be a positive numeric value.", call. = FALSE)
    }

    repel_max_iter <- as.integer(repel_max_iter)

    if (length(repel_max_iter) != 1 || is.na(repel_max_iter) ||
        repel_max_iter <= 0) {
        stop("`repel_max_iter` must be a positive integer.", call. = FALSE)
    }

    if (!is.numeric(set_label_offset) || length(set_label_offset) != 1 ||
        is.na(set_label_offset) || set_label_offset < 0) {
        stop("`set_label_offset` must be a non-negative numeric value.", call. = FALSE)
    }

    # ---------------------------------------------------------------------------
    # Palette handling
    # ---------------------------------------------------------------------------

    resolve_palette <- function(palette, n = 3) {
        if (is.null(palette)) {
            return(c("#E69F00", "#56B4E9", "#009E73")[seq_len(n)])
        }

        if (is.function(palette)) {
            cols <- palette(n)

            if (length(cols) < n) {
                stop(
                    "The supplied palette function returned fewer than ",
                    n,
                    " colours.",
                    call. = FALSE
                )
            }

            return(cols[seq_len(n)])
        }

        if (is.character(palette)) {
            if (length(palette) >= n) {
                return(palette[seq_len(n)])
            }

            if (length(palette) == 1) {
                pal_name <- tolower(palette)

                if (pal_name %in% c("ggplot", "ggplot2", "hue")) {
                    return(scales::hue_pal()(n))
                }

                if (pal_name %in% c(
                    "okabe-ito",
                    "okabe_ito",
                    "okabe ito",
                    "colourblind",
                    "colorblind"
                )) {
                    return(c("#E69F00", "#56B4E9", "#009E73")[seq_len(n)])
                }

                if (pal_name %in% c("nejm", "npg", "aaas", "jama", "lancet")) {
                    if (!requireNamespace("ggsci", quietly = TRUE)) {
                        stop(
                            "Package 'ggsci' is required for palette = '",
                            palette,
                            "'.",
                            call. = FALSE
                        )
                    }

                    pal_fun <- switch(
                        pal_name,
                        "npg" = ggsci::pal_npg(),
                        "nejm" = ggsci::pal_nejm(),
                        "aaas" = ggsci::pal_aaas(),
                        "jama" = ggsci::pal_jama(),
                        "lancet" = ggsci::pal_lancet()
                    )

                    return(pal_fun(n))
                }
            }
        }

        stop(
            "`palette` must be NULL, a colour vector of length >= 3, ",
            "a palette function, or one of: 'okabe-ito', 'ggplot2', ",
            "'npg', 'nejm', 'aaas', 'jama', 'lancet'.",
            call. = FALSE
        )
    }

    set_cols <- resolve_palette(palette, 3)
    names(set_cols) <- c("X", "Y", "Z")

    # ---------------------------------------------------------------------------
    # Set logic
    # ---------------------------------------------------------------------------

    compute_venn_regions <- function(list_x, list_y, list_z) {
        list_x <- unique(list_x)
        list_y <- unique(list_y)
        list_z <- unique(list_z)

        universe <- unique(c(list_x, list_y, list_z))

        if (length(universe) == 0) {
            stop(
                "At least one of `list_x`, `list_y`, or `list_z` must contain an identifier.",
                call. = FALSE
            )
        }

        in_x <- universe %in% list_x
        in_y <- universe %in% list_y
        in_z <- universe %in% list_z

        region <- paste0(
            as.integer(in_x),
            as.integer(in_y),
            as.integer(in_z)
        )

        split_regions <- split(universe, region)

        get_region <- function(code) {
            if (code %in% names(split_regions)) {
                split_regions[[code]]
            } else {
                character(0)
            }
        }

        regions <- list(
            x_only = get_region("100"),
            y_only = get_region("010"),
            z_only = get_region("001"),
            xy_only = get_region("110"),
            xz_only = get_region("101"),
            yz_only = get_region("011"),
            xyz = get_region("111")
        )

        region_counts <- vapply(regions, length, integer(1))

        counts <- c(
            x = region_counts[["x_only"]] +
                region_counts[["xy_only"]] +
                region_counts[["xz_only"]] +
                region_counts[["xyz"]],
            y = region_counts[["y_only"]] +
                region_counts[["xy_only"]] +
                region_counts[["yz_only"]] +
                region_counts[["xyz"]],
            z = region_counts[["z_only"]] +
                region_counts[["xz_only"]] +
                region_counts[["yz_only"]] +
                region_counts[["xyz"]],
            xy = region_counts[["xy_only"]] +
                region_counts[["xyz"]],
            xz = region_counts[["xz_only"]] +
                region_counts[["xyz"]],
            yz = region_counts[["yz_only"]] +
                region_counts[["xyz"]],
            xyz = region_counts[["xyz"]],
            region_counts
        )

        lists <- c(
            list(
                x = list_x,
                y = list_y,
                z = list_z,
                xy = c(regions$xy_only, regions$xyz),
                xz = c(regions$xz_only, regions$xyz),
                yz = c(regions$yz_only, regions$xyz)
            ),
            regions
        )

        list(
            lists = lists,
            counts = counts,
            region_counts = region_counts
        )
    }

    venn <- compute_venn_regions(list_x, list_y, list_z)

    counts <- venn$counts
    total_unique <- sum(venn$region_counts)

    # ---------------------------------------------------------------------------
    # Geometry helpers
    # ---------------------------------------------------------------------------

    circle_overlap_area <- function(r1, r2, d) {
        if (r1 <= 0 || r2 <= 0) {
            return(0)
        }

        if (d >= r1 + r2) {
            return(0)
        }

        if (d <= abs(r1 - r2)) {
            return(pi * min(r1, r2)^2)
        }

        cos1 <- (d^2 + r1^2 - r2^2) / (2 * d * r1)
        cos2 <- (d^2 + r2^2 - r1^2) / (2 * d * r2)

        cos1 <- min(1, max(-1, cos1))
        cos2 <- min(1, max(-1, cos2))

        term1 <- r1^2 * acos(cos1)
        term2 <- r2^2 * acos(cos2)

        term3 <- 0.5 * sqrt(
            max(
                0,
                (-d + r1 + r2) *
                    (d + r1 - r2) *
                    (d - r1 + r2) *
                    (d + r1 + r2)
            )
        )

        term1 + term2 - term3
    }

    distance_for_overlap <- function(r1, r2, overlap, tol = 1e-10) {
        if (r1 <= 0 || r2 <= 0) {
            return(r1 + r2)
        }

        max_overlap <- pi * min(r1, r2)^2

        if (overlap <= tol) {
            return(r1 + r2)
        }

        if (overlap >= max_overlap - tol) {
            return(abs(r1 - r2))
        }

        stats::uniroot(
            function(d) circle_overlap_area(r1, r2, d) - overlap,
            lower = abs(r1 - r2),
            upper = r1 + r2,
            tol = tol
        )$root
    }

    place_circle_centres <- function(r, overlaps, eps = 1e-10) {
        active <- r > 0

        centres <- data.frame(
            set = names(r),
            x = 0,
            y = 0,
            r = as.numeric(r),
            stringsAsFactors = FALSE
        )

        if (sum(active) == 1) {
            return(centres)
        }

        d <- c(
            XY = distance_for_overlap(r[["X"]], r[["Y"]], overlaps[["XY"]]),
            XZ = distance_for_overlap(r[["X"]], r[["Z"]], overlaps[["XZ"]]),
            YZ = distance_for_overlap(r[["Y"]], r[["Z"]], overlaps[["YZ"]])
        )

        active_sets <- names(r)[active]

        if (sum(active) == 2) {
            pair <- paste(active_sets, collapse = "")

            pair_distance <- switch(
                pair,
                "XY" = d[["XY"]],
                "XZ" = d[["XZ"]],
                "YZ" = d[["YZ"]]
            )

            centres$x[centres$set == active_sets[1]] <- 0
            centres$y[centres$set == active_sets[1]] <- 0
            centres$x[centres$set == active_sets[2]] <- pair_distance
            centres$y[centres$set == active_sets[2]] <- 0

            return(centres)
        }

        original_d <- d

        d[["XY"]] <- min(d[["XY"]], d[["XZ"]] + d[["YZ"]])
        d[["XZ"]] <- min(d[["XZ"]], d[["XY"]] + d[["YZ"]])
        d[["YZ"]] <- min(d[["YZ"]], d[["XY"]] + d[["XZ"]])

        if (any(abs(original_d - d) > eps)) {
            warning(
                "Requested pairwise overlaps cannot be represented exactly by ",
                "three circles. Distances were adjusted to satisfy the triangle ",
                "inequality.",
                call. = FALSE
            )
        }

        if (d[["XY"]] < eps) {
            d[["XY"]] <- eps
        }

        z_x <- (d[["XZ"]]^2 + d[["XY"]]^2 - d[["YZ"]]^2) /
            (2 * d[["XY"]])

        z_y <- sqrt(max(0, d[["XZ"]]^2 - z_x^2))

        centres$x[centres$set == "X"] <- 0
        centres$y[centres$set == "X"] <- 0

        centres$x[centres$set == "Y"] <- d[["XY"]]
        centres$y[centres$set == "Y"] <- 0

        centres$x[centres$set == "Z"] <- z_x
        centres$y[centres$set == "Z"] <- z_y

        centres
    }

    # ---------------------------------------------------------------------------
    # Compute circle layout
    # ---------------------------------------------------------------------------

    r <- c(
        X = sqrt(counts[["x"]] / pi),
        Y = sqrt(counts[["y"]] / pi),
        Z = sqrt(counts[["z"]] / pi)
    )

    overlaps <- c(
        XY = counts[["xy"]],
        XZ = counts[["xz"]],
        YZ = counts[["yz"]]
    )

    centres <- place_circle_centres(r, overlaps)

    centres$title <- c(xtitle, ytitle, ztitle)
    centres$colour <- set_cols[centres$set]

    centres <- centres[centres$r > 0, , drop = FALSE]

    centres$x <- centres$x - mean(range(centres$x))
    centres$y <- centres$y - mean(range(centres$y))

    # ---------------------------------------------------------------------------
    # Label placement via region-aware anchors
    # ---------------------------------------------------------------------------

    bbox <- data.frame(
        xmin = min(centres$x - centres$r),
        xmax = max(centres$x + centres$r),
        ymin = min(centres$y - centres$r),
        ymax = max(centres$y + centres$r)
    )

    pad_fraction <- if (isTRUE(repel_labels)) 0.16 else 0.08

    pad <- pad_fraction * max(
        bbox$xmax - bbox$xmin,
        bbox$ymax - bbox$ymin
    )

    bbox$xmin <- bbox$xmin - pad
    bbox$xmax <- bbox$xmax + pad
    bbox$ymin <- bbox$ymin - pad
    bbox$ymax <- bbox$ymax + pad

    xs <- seq(bbox$xmin, bbox$xmax, length.out = label_grid)
    ys <- seq(bbox$ymin, bbox$ymax, length.out = label_grid)

    grid <- expand.grid(
        x = xs,
        y = ys
    )

    get_centre_value <- function(set_name, column, missing_value) {
        value <- centres[centres$set == set_name, column]

        if (length(value) == 0) {
            missing_value
        } else {
            value
        }
    }

    cx <- c(
        X = get_centre_value("X", "x", Inf),
        Y = get_centre_value("Y", "x", Inf),
        Z = get_centre_value("Z", "x", Inf)
    )

    cy <- c(
        X = get_centre_value("X", "y", Inf),
        Y = get_centre_value("Y", "y", Inf),
        Z = get_centre_value("Z", "y", Inf)
    )

    rr <- c(
        X = get_centre_value("X", "r", 0),
        Y = get_centre_value("Y", "r", 0),
        Z = get_centre_value("Z", "r", 0)
    )

    active_circle_names <- names(rr)[rr > 0 & is.finite(cx) & is.finite(cy)]

    inside_x <- (grid$x - cx[["X"]])^2 +
        (grid$y - cy[["X"]])^2 <= rr[["X"]]^2

    inside_y <- (grid$x - cx[["Y"]])^2 +
        (grid$y - cy[["Y"]])^2 <= rr[["Y"]]^2

    inside_z <- (grid$x - cx[["Z"]])^2 +
        (grid$y - cy[["Z"]])^2 <= rr[["Z"]]^2

    grid$region <- paste0(
        as.integer(inside_x),
        as.integer(inside_y),
        as.integer(inside_z)
    )

    region_counts <- c(
        "100" = counts[["x_only"]],
        "010" = counts[["y_only"]],
        "001" = counts[["z_only"]],
        "110" = counts[["xy_only"]],
        "101" = counts[["xz_only"]],
        "011" = counts[["yz_only"]],
        "111" = counts[["xyz"]]
    )

    format_count_label <- function(count, type, denominator) {
        if (type == "none") {
            return("")
        }

        pct <- if (denominator > 0) {
            round(count / denominator * 100, 2)
        } else {
            NA_real_
        }

        if (type == "abs") {
            as.character(count)
        } else if (type == "pct") {
            paste0(pct, "%")
        } else if (type == "abs_pct") {
            paste0(count, " (", pct, "%)")
        } else {
            ""
        }
    }

    fallback_centroid <- function(region_code) {
        bits <- strsplit(region_code, "", fixed = TRUE)[[1]] == "1"

        centre_mat <- data.frame(
            set = c("X", "Y", "Z"),
            x = c(cx[["X"]], cx[["Y"]], cx[["Z"]]),
            y = c(cy[["X"]], cy[["Y"]], cy[["Z"]]),
            stringsAsFactors = FALSE
        )

        active_centres <- centre_mat[bits, , drop = FALSE]
        active_centres <- active_centres[
            is.finite(active_centres$x) & is.finite(active_centres$y),
            ,
            drop = FALSE
        ]

        if (nrow(active_centres) == 0) {
            return(c(x = 0, y = 0))
        }

        c(
            x = mean(active_centres$x),
            y = mean(active_centres$y)
        )
    }

    boundary_clearance <- function(points) {
        if (nrow(points) == 0 || length(active_circle_names) == 0) {
            return(numeric(0))
        }

        dmat <- vapply(
            active_circle_names,
            function(set_name) {
                abs(
                    sqrt(
                        (points$x - cx[[set_name]])^2 +
                            (points$y - cy[[set_name]])^2
                    ) - rr[[set_name]]
                )
            },
            numeric(nrow(points))
        )

        if (is.null(dim(dmat))) {
            dmat
        } else {
            apply(dmat, 1, min)
        }
    }

    region_anchor <- function(region_code) {
        points <- grid[grid$region == region_code, c("x", "y"), drop = FALSE]
        fallback <- fallback_centroid(region_code)

        if (nrow(points) == 0) {
            return(fallback)
        }

        clearance <- boundary_clearance(points)

        fallback_distance <- sqrt(
            (points$x - fallback[["x"]])^2 +
                (points$y - fallback[["y"]])^2
        )

        score <- clearance - 1e-6 * fallback_distance

        best <- which.max(score)

        c(
            x = points$x[best],
            y = points$y[best]
        )
    }

    label_regions <- names(region_counts)[region_counts > 0]

    if (length(label_regions) > 0) {
        label_df <- do.call(
            rbind,
            lapply(label_regions, function(region_code) {
                xy <- region_anchor(region_code)

                data.frame(
                    region = region_code,
                    x = xy[["x"]],
                    y = xy[["y"]],
                    count = unname(region_counts[[region_code]]),
                    stringsAsFactors = FALSE
                )
            })
        )
    } else {
        label_df <- data.frame(
            region = character(0),
            x = numeric(0),
            y = numeric(0),
            count = integer(0),
            stringsAsFactors = FALSE
        )
    }

    if (nrow(label_df) > 0 && nrtype != "none") {
        label_df$label <- vapply(
            label_df$count,
            format_count_label,
            character(1),
            type = nrtype,
            denominator = total_unique
        )
    } else {
        label_df <- label_df[0, , drop = FALSE]
        label_df$label <- character(0)
    }

    # ---------------------------------------------------------------------------
    # Set-label placement
    # ---------------------------------------------------------------------------

    outward_direction_for_set <- function(set_name) {
        centre_x <- cx[[set_name]]
        centre_y <- cy[[set_name]]

        unique_code <- switch(
            set_name,
            X = "100",
            Y = "010",
            Z = "001"
        )

        if (!is.null(region_counts[[unique_code]]) &&
            region_counts[[unique_code]] > 0 &&
            any(grid$region == unique_code)) {
            unique_xy <- region_anchor(unique_code)

            direction <- c(
                x = unique_xy[["x"]] - centre_x,
                y = unique_xy[["y"]] - centre_y
            )

            norm <- sqrt(sum(direction^2))

            if (is.finite(norm) && norm > 1e-12) {
                return(direction / norm)
            }
        }

        other_centres <- centres[
            centres$set != set_name &
                is.finite(centres$x) &
                is.finite(centres$y),
            ,
            drop = FALSE
        ]

        if (nrow(other_centres) > 0) {
            other_mean <- c(
                x = mean(other_centres$x),
                y = mean(other_centres$y)
            )

            direction <- c(
                x = centre_x - other_mean[["x"]],
                y = centre_y - other_mean[["y"]]
            )

            norm <- sqrt(sum(direction^2))

            if (is.finite(norm) && norm > 1e-12) {
                return(direction / norm)
            }
        }

        fallback_angle <- switch(
            set_name,
            X = pi,
            Y = 0,
            Z = pi / 2
        )

        c(
            x = cos(fallback_angle),
            y = sin(fallback_angle)
        )
    }

    plot_span <- max(
        bbox$xmax - bbox$xmin,
        bbox$ymax - bbox$ymin
    )

    offset_abs <- set_label_offset * plot_span

    set_label_pos_df <- centres[, c("set", "x", "y", "r", "title"), drop = FALSE]

    if (set_label_position == "outside_unique" && nrow(set_label_pos_df) > 0) {
        for (i in seq_len(nrow(set_label_pos_df))) {
            set_name <- set_label_pos_df$set[i]
            direction <- outward_direction_for_set(set_name)

            set_label_pos_df$x[i] <- cx[[set_name]] +
                direction[["x"]] * (rr[[set_name]] + offset_abs)

            set_label_pos_df$y[i] <- cy[[set_name]] +
                direction[["y"]] * (rr[[set_name]] + offset_abs)
        }

        label_pad <- 0.04 * plot_span

        bbox$xmin <- min(bbox$xmin, set_label_pos_df$x - label_pad)
        bbox$xmax <- max(bbox$xmax, set_label_pos_df$x + label_pad)
        bbox$ymin <- min(bbox$ymin, set_label_pos_df$y - label_pad)
        bbox$ymax <- max(bbox$ymax, set_label_pos_df$y + label_pad)
    }

    # ---------------------------------------------------------------------------
    # Circle fill and stroke aesthetics
    # ---------------------------------------------------------------------------

    circle_df <- centres
    circle_df$fill <- scales::alpha(circle_df$colour, fill_alpha)

    if (is.null(stroke_colour)) {
        circle_df$stroke <- if (stroke_width > 0) {
            scales::alpha(circle_df$colour, stroke_alpha)
        } else {
            NA_character_
        }
    } else if (length(stroke_colour) == 1) {
        circle_df$stroke <- stroke_colour
    } else if (length(stroke_colour) >= nrow(circle_df)) {
        circle_df$stroke <- stroke_colour[seq_len(nrow(circle_df))]
    } else {
        stop(
            "`stroke_colour` must be NULL, NA, a single colour, or a vector ",
            "with at least as many colours as active circles.",
            call. = FALSE
        )
    }

    # ---------------------------------------------------------------------------
    # Combined label data
    # ---------------------------------------------------------------------------

    plot_label_parts <- list()

    if (nrow(label_df) > 0) {
        region_label_df <- data.frame(
            x = label_df$x,
            y = label_df$y,
            label = label_df$label,
            label_type = "region",
            colour = label_colour,
            size = label_size,
            fontface = label_fontface,
            stringsAsFactors = FALSE
        )

        plot_label_parts[["region"]] <- region_label_df
    }

    if (isTRUE(show_set_labels) && nrow(circle_df) > 0) {
        set_label_source <- set_label_pos_df[
            match(circle_df$set, set_label_pos_df$set),
            ,
            drop = FALSE
        ]

        set_counts <- c(
            X = counts[["x"]],
            Y = counts[["y"]],
            Z = counts[["z"]]
        )

        active_set_counts <- set_counts[circle_df$set]

        set_pct_denominator_value <- if (set_pct_denominator == "union") {
            total_unique
        } else {
            sum(set_counts[circle_df$set])
        }

        set_count_labels <- vapply(
            active_set_counts,
            format_count_label,
            character(1),
            type = set_count_nrtype,
            denominator = set_pct_denominator_value
        )

        set_labels <- set_label_source$title

        if (set_count_nrtype != "none") {
            set_labels <- paste0(set_labels, "\n", set_count_labels)
        }

        set_label_df <- data.frame(
            x = set_label_source$x,
            y = set_label_source$y,
            label = set_labels,
            label_type = "set",
            colour = set_label_colour,
            size = set_label_size,
            fontface = set_label_fontface,
            stringsAsFactors = FALSE
        )

        plot_label_parts[["set"]] <- set_label_df
    }

    if (length(plot_label_parts) > 0) {
        plot_label_df <- do.call(rbind, plot_label_parts)
        rownames(plot_label_df) <- NULL
    } else {
        plot_label_df <- data.frame(
            x = numeric(0),
            y = numeric(0),
            label = character(0),
            label_type = character(0),
            colour = character(0),
            size = numeric(0),
            fontface = character(0),
            stringsAsFactors = FALSE
        )
    }

    # ---------------------------------------------------------------------------
    # Plot construction
    # ---------------------------------------------------------------------------

    p <- ggplot2::ggplot() +
        ggforce::geom_circle(
            data = circle_df,
            ggplot2::aes(
                x0 = x,
                y0 = y,
                r = r,
                fill = fill,
                colour = stroke
            ),
            linewidth = stroke_width,
            inherit.aes = FALSE
        ) +
        ggplot2::scale_fill_identity() +
        ggplot2::scale_colour_identity() +
        ggplot2::scale_size_identity() +
        ggplot2::coord_equal(
            xlim = c(bbox$xmin, bbox$xmax),
            ylim = c(bbox$ymin, bbox$ymax),
            expand = FALSE,
            clip = "off"
        ) +
        ggplot2::labs(
            title = title,
            subtitle = subtitle
        ) +
        ggplot2::theme_void() +
        ggplot2::theme(
            plot.background = ggplot2::element_rect(fill = bg, colour = NA),
            panel.background = ggplot2::element_rect(fill = bg, colour = NA),
            plot.title = ggplot2::element_text(
                hjust = 0.5,
                face = "bold"
            ),
            plot.subtitle = ggplot2::element_text(
                hjust = 0.5
            ),
            legend.position = "none",
            plot.margin = ggplot2::margin(12, 12, 12, 12)
        )

    if (nrow(plot_label_df) > 0) {
        if (isTRUE(repel_labels)) {
            p <- p +
                ggrepel::geom_text_repel(
                    data = plot_label_df,
                    ggplot2::aes(
                        x = x,
                        y = y,
                        label = label,
                        colour = colour,
                        size = size,
                        fontface = fontface
                    ),
                    box.padding = repel_box_padding,
                    point.padding = repel_point_padding,
                    force = repel_force,
                    force_pull = repel_force_pull,
                    max.iter = repel_max_iter,
                    max.time = repel_max_time,
                    max.overlaps = Inf,
                    seed = repel_seed,
                    segment.colour = repel_segment_colour,
                    inherit.aes = FALSE
                )
        } else {
            p <- p +
                ggplot2::geom_text(
                    data = plot_label_df,
                    ggplot2::aes(
                        x = x,
                        y = y,
                        label = label,
                        colour = colour,
                        size = size,
                        fontface = fontface
                    ),
                    inherit.aes = FALSE
                )
        }
    }

    # ---------------------------------------------------------------------------
    # Optional saving
    # ---------------------------------------------------------------------------

    if (output != "plot") {
        if (is.null(filename)) {
            filename <- paste0("biovenn.", output)
        }

        device <- switch(
            output,
            jpg = "jpeg",
            tif = "tiff",
            output
        )

        ggplot2::ggsave(
            filename = filename,
            plot = p,
            device = device,
            width = width,
            height = height,
            dpi = dpi,
            bg = bg
        )
    }

    # ---------------------------------------------------------------------------
    # Return
    # ---------------------------------------------------------------------------

    venn_lists <- c(
        venn$lists,
        list(
            counts = counts,
            region_counts = venn$region_counts
        )
    )

    attr(p, "venn_lists") <- venn_lists
    attr(p, "circle_data") <- circle_df
    attr(p, "label_data") <- label_df
    attr(p, "plot_label_data") <- plot_label_df
    attr(p, "set_label_data") <- set_label_pos_df

    if (isTRUE(return_lists)) {
        return(c(list(plot = p), venn_lists))
    }

    p
}

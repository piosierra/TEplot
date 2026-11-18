#!/usr/bin/env Rscript
### Tool to draw selfdotplots or pairwise dotplots, with ORF tracks
### Libs
if (suppressPackageStartupMessages(!require("data.table", quietly = TRUE))) {
  install.packages("data.table")
}
if (suppressPackageStartupMessages(!require("stringi", quietly = TRUE))) {
  install.packages("stringi")
}
if (suppressPackageStartupMessages(!require("ggplot2", quietly = TRUE))) {
  install.packages("ggplot2")
}

### Appearance settings (edit here)
SELF_STYLE  <- "triangle"          # self-dotplots: "triangle" (rotated upper half) or "square"
DOT_SIZE_MM <- 1.2                 # minimum dot diameter in mm
COL_FWD     <- "#4F7A1A"           # direct repeats / forward matches (blue)
COL_RC      <- "#D63384"           # inverted repeats / reverse-complement (orange)
COL_DIAG    <- "grey75"            # main diagonal in self-dotplots
COL_ORF_P   <- "#2F5D8A"
COL_ORF_M   <- "#B8700F" 
COL_TRACK   <- "grey96"            # ORF lane background
PANEL_IN    <- 10                  # length of the longer axis, in inches
DPI         <- 300

### ORF settings
ORF_MIN_AA      <- 200             # report ORFs strictly longer than this
ORF_REQUIRE_ATG <- TRUE            # TRUE: ATG..stop; FALSE: stop..stop
LANE_FRAC       <- 0.025           # height of one ORF lane, fraction of longer axis

### Functions
ffasta <- function(f) {
  check <- file.size(f)
  if (check == 0 | is.na(check)) {
    return(data.table(name = as.character(), seq = as.character()))
  } else {
    fa_raw <- fread(f, header = F, fill = T, sep = "\n")
    fa_raw[, h := grepl(">", V1)]
    fa_oneline <- fa_raw[, .(paste0(V1, collapse = "")), by = rleid(h)]
    return(data.table(
      name = fa_oneline[rep(c(TRUE, FALSE),  length = .N)]$V1,
      seq  = fa_oneline[rep(c(FALSE, TRUE),  length = .N)]$V1
    ))
  }
}

rc <- function(x) {
  return(stri_reverse(chartr("ATGCNatgcn", "TACGNTACGN", x)))
}

get_kmers <- function(data, k) {
  l <- nchar(data)
  if (l < k) return(data.table(pos = integer(), kmer = character()))
  pos <- seq_len(l - k + 1)   # include the last k-mer
  data.table(pos = pos, kmer = stri_sub(data, pos, pos + (k - 1)))
}

# Safe file name (very long FASTA headers can exceed filesystem limits)
safe_name <- function(x, max_len = 100) {
  substr(make.names(x), 1, max_len)
}

# Wrap long labels so titles never run off the image
wrap <- function(x, width = 70) paste(strwrap(x, width = width), collapse = "\n")

# Matches are placed at the CENTRE of each k-mer (start + (k - 1) / 2).
# Start positions are fine for direct repeats, but for inverted repeats the
# k-mer starting at i pairs with the k-mer ENDING at the mirror position.
kmer_centre <- function(pos, kmers) pos + (nchar(kmers$kmer[1]) - 1) / 2

match_pairs <- function(kmers, key) {
  a <- kmers[, .(V1 = pos, kmer)]
  b <- kmers[, .(V2 = pos, kmer = get(key))]
  m <- merge(a, b, by = "kmer", allow.cartesian = TRUE)
  m[V1 < V2, .(V1 = kmer_centre(as.numeric(V1), kmers),
               V2 = kmer_centre(as.numeric(V2), kmers))]
}

cross_pairs <- function(kmersA, kmersB, keyA) {
  a <- kmersA[, .(V1 = pos, kmer = get(keyA))]
  b <- kmersB[, .(V2 = pos, kmer)]
  m <- merge(a, b, by = "kmer", allow.cartesian = TRUE)
  m[, .(V1 = kmer_centre(as.numeric(V1), kmersA),
        V2 = kmer_centre(as.numeric(V2), kmersB))]
}
### ORF finding
# ORFs on the given strand of s (already upper case). Coordinates are 1-based
# on s, end includes the stop codon. ORFs that run off the end of the
# sequence without a stop are kept (complete = FALSE): consensus TEs are
# often truncated.
orfs_one_strand <- function(s, min_aa, require_atg) {
  L <- nchar(s)
  out <- list()
  for (f in 1:3) {
    if (L - f + 1 < 3) next
    st  <- seq(f, L - 2, by = 3)
    cod <- stri_sub(s, st, length = 3)
    is_stop <- cod %chin% c("TAA", "TAG", "TGA")
    # segment = codons up to and including the next stop
    seg <- c(0L, cumsum(is_stop)[-length(is_stop)])
    dt  <- data.table(st = st, seg = seg, stop = is_stop, atg = cod == "ATG")
    o <- dt[, {
      last_coding <- if (stop[.N]) .N - 1L else .N
      first <- if (require_atg) which(atg)[1] else 1L
      if (is.na(first) || first > last_coding) NULL else
        list(start = st[first], end = st[.N] + 2L,
             aa = last_coding - first + 1L, complete = stop[.N])
    }, by = seg]
    if (nrow(o)) out[[length(out) + 1]] <- o[aa > min_aa, .(start, end, aa, complete, frame = f)]
  }
  r <- rbindlist(out)
  if (!nrow(r)) r <- data.table(start = integer(), end = integer(), aa = integer(),
                                complete = logical(), frame = integer())
  r
}

# ORFs on both strands, in forward coordinates (start < end always)
find_orfs <- function(s, min_aa = ORF_MIN_AA, require_atg = ORF_REQUIRE_ATG) {
  s <- toupper(s)
  L <- nchar(s)
  plus  <- orfs_one_strand(s, min_aa, require_atg)[, strand := "+"]
  minus <- orfs_one_strand(rc(s), min_aa, require_atg)
  minus[, c("start", "end") := .(L - end + 1L, L - start + 1L)]
  minus[, strand := "-"]
  rbind(plus, minus)
}

### ORF track geometry
# Build arrow polygons for ORFs drawn in a track. `along` is the axis the
# sequence runs on (x for the top track, y for the right track); `base` is
# where the track starts on the other axis. Each ORF has a `lane` (1 = next
# to the dotplot). Arrows point in the direction of translation.
# side: "top" / "bottom" (horizontal tracks) or "right" (vertical track);
# "bottom" stacks lanes downwards from `base`.
orf_layers <- function(orfs, base, lane_h, seq_len_bp, side) {
  dir <- if (side == "bottom") -1 else 1
  if (!nrow(orfs)) return(list())
  o <- copy(orfs)
  pad <- lane_h * 0.15
  o[, `:=`(id = .I,
           b0 = base + dir * ((lane - 1) * lane_h + pad),
           b1 = base + dir * (lane * lane_h - pad),
           h  = pmin(0.012 * seq_len_bp, (end - start) * 0.4))]
  o[, bm := (b0 + b1) / 2]
  poly <- o[, {
    if (strand == "+") {
      a <- c(start, end - h, end, end - h, start)
    } else {
      a <- c(end, start + h, start, start + h, end)
    }
    list(a = a, b = c(b0, b0, bm, b1, b1))
  }, by = .(id, complete, strand)]
  poly[, fill := ifelse(strand == "+", COL_ORF_P, COL_ORF_M)]
  lab <- o[(end - start) > 0.05 * seq_len_bp,
           .(a = (start + end) / 2, b = bm, txt = paste0(aa, " aa"))]

  if (side != "right") {
    poly[, `:=`(x = a, y = b)]; lab[, `:=`(x = a, y = b)]; ang <- 0
  } else {
    poly[, `:=`(x = b, y = a)]; lab[, `:=`(x = b, y = a)]; ang <- -90
  }
  layers <- list(
    geom_polygon(data = poly, aes(x = x, y = y, group = id),
                 inherit.aes = FALSE, fill = poly$fill, colour = NA)
  )
  if (nrow(lab)) {
    layers <- c(layers, list(
      geom_text(data = lab, aes(x = x, y = y, label = txt), inherit.aes = FALSE,
                colour = "white", size = 2.9, angle = ang, fontface = "bold")
    ))
  }
  layers
}

# Lane backgrounds + lane labels ("+1", "-2", ...) for one track
track_frame <- function(lane_names, base, lane_h, lo, hi, side, label_at) {
  n <- length(lane_names)
  if (!n) return(list())
  dir <- if (side == "bottom") -1 else 1
  bg <- data.table(l = seq_len(n), b0 = base + dir * (seq_len(n) - 1) * lane_h,
                   b1 = base + dir * seq_len(n) * lane_h)
  bg[, lab := lane_names]
  bg[, col := ifelse(startsWith(lab, "+"), COL_ORF_P, COL_ORF_M)]
  if (side != "right") {
    list(
      annotate("rect", xmin = lo, xmax = hi, ymin = bg$b0, ymax = bg$b1,
               fill = COL_TRACK, colour = "white", linewidth = 0.6),
      annotate("text", x = label_at, y = (bg$b0 + bg$b1) / 2, label = bg$lab,
               hjust = 1, size = 3.5, colour = bg$col, fontface = "bold")
    )
  } else {
    list(
      annotate("rect", ymin = lo, ymax = hi, xmin = bg$b0, xmax = bg$b1,
               fill = COL_TRACK, colour = "white", linewidth = 0.6),
      annotate("text", y = label_at, x = (bg$b0 + bg$b1) / 2, label = bg$lab,
               vjust = 0, size = 3.5, colour = bg$col, fontface = "bold")
    )
  }
}

# Assign lanes: frames 1-3 of one strand, or + then - for a 6-lane track
assign_lanes <- function(orfs, strands) {
  o <- orfs[strand %in% strands]
  o[, lane := frame + 3L * (match(strand, strands) - 1L)]
  o
}
lane_names_for <- function(strands) {
  as.vector(t(outer(strands, 1:3, paste0)))
}

### Layout helpers
# Nice axis breaks (1, 2, 2.5 or 5 x 10^n), aiming for ~n_ticks
nice_breaks <- function(l, n_ticks = 5) {
  b <- pretty(c(0, l), n = n_ticks)
  b[b >= 0 & b <= l]
}

# Output size in inches for a plotting region of xr x yr data units, with the
# same physical scale on both axes and the longer side = PANEL_IN.
plot_dims <- function(xr, yr) {
  s <- PANEL_IN / max(xr, yr)
  list(width  = max(xr * s + 1.4, 6),
       height = yr * s + 2.0)
}

# Shared plotting function for self and pairwise dotplots.
# top_orfs / right_orfs: ORF tables with a `lane` column; top_lanes /
# right_lanes: lane labels (track height follows the number of lanes).
draw_dotplot <- function(dots, lenX, lenY, title, subtitle, xlab, ylab, fname,
                         top_orfs, top_lanes, right_orfs, right_lanes,
                         diagonal = FALSE) {
  dots[, type := factor(type, levels = c("fwd", "rc"))]
  setorder(dots, -type)            # draw forward matches on top

  M      <- max(lenX, lenY)
  pad    <- 0.01 * M               # breathing room around the dotplot
  gap    <- 0.012 * M              # space between dotplot and tracks
  lane_h <- LANE_FRAC * M
  loX <- 1 - pad; hiX <- lenX + pad
  loY <- 1 - pad; hiY <- lenY + pad
  topBase   <- hiY + gap
  rightBase <- hiX + gap
  topEnd    <- topBase   + length(top_lanes)   * lane_h
  rightEnd  <- rightBase + length(right_lanes) * lane_h

  bp_mm    <- PANEL_IN * 25.4 / max(rightEnd - loX, topEnd - loY)
  dot_size <- min(max(DOT_SIZE_MM, bp_mm), 3)
  dims     <- plot_dims(rightEnd - loX, topEnd - loY)

  brksX <- nice_breaks(lenX); brksY <- nice_breaks(lenY)
  grid  <- rbind(
    data.table(x = brksX, xend = brksX, y = loY, yend = hiY),
    data.table(x = loX, xend = hiX, y = brksY, yend = brksY)
  )

  p <- ggplot(dots, aes(x = V1, y = V2, colour = type)) +
    # dotplot frame and grid (drawn by hand so they stop at the tracks)
    annotate("rect", xmin = loX, xmax = hiX, ymin = loY, ymax = hiY,
             fill = "white", colour = NA) +
    geom_segment(data = grid, aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE, colour = "grey92", linewidth = 0.3)
  if (diagonal) {
    p <- p + annotate("segment", x = 1, y = 1, xend = lenX, yend = lenY,
                      colour = COL_DIAG, linewidth = dot_size * 0.6)
  }
  p <- p +
    geom_point(shape = 16, size = dot_size, stroke = 0) +
    annotate("rect", xmin = loX, xmax = hiX, ymin = loY, ymax = hiY,
             fill = NA, colour = "grey40", linewidth = 0.6) +
    # ORF tracks
    track_frame(top_lanes,   topBase,   lane_h, loX, hiX, "top",   loX - pad) +
    track_frame(right_lanes, rightBase, lane_h, loY, hiY, "right", hiY + pad * 0.3) +
    orf_layers(top_orfs,   topBase,   lane_h, lenX, "top") +
    orf_layers(right_orfs, rightBase, lane_h, lenY, "right") +
    scale_colour_manual(
      values = c(fwd = COL_FWD, rc = COL_RC),
      labels = c(fwd = "Direct (same strand)",
                 rc  = "Inverted (reverse complement)"),
      drop = FALSE, name = NULL
    ) +
    scale_x_continuous(breaks = brksX, labels = scales::comma) +
    scale_y_continuous(breaks = brksY, labels = scales::comma) +
    coord_fixed(ratio = 1, xlim = c(loX, rightEnd), ylim = c(loY, topEnd),
                expand = FALSE, clip = "off") +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title          = element_text(face = "bold", size = 15),
      plot.subtitle       = element_text(colour = "grey35", size = 11,
                                         margin = margin(b = 10)),
      plot.title.position = "plot",
      axis.title          = element_text(size = 12, colour = "grey20"),
      axis.title.x        = element_text(hjust = (lenX / 2 - loX) / (rightEnd - loX)),
      axis.title.y        = element_text(hjust = (lenY / 2 - loY) / (topEnd - loY)),
      axis.text           = element_text(size = 11, colour = "grey20"),
      axis.ticks          = element_line(colour = "grey40", linewidth = 0.4),
      panel.grid          = element_blank(),
      legend.position     = "bottom",
      legend.text         = element_text(size = 11),
      plot.margin         = margin(12, 16, 8, 12)
    )

  ggsave(fname, plot = p, width = dims$width, height = dims$height,
         units = "in", dpi = DPI, bg = "white")
}

orf_note <- function() {
  sprintf("ORFs > %d aa (%s)", ORF_MIN_AA,
          if (ORF_REQUIRE_ATG) "ATG to stop" else "stop to stop")
}

### Sequence statistics (cheap: reuse the k-mer table)
# Returns a list of label/value pairs shown beside the triangle.
#  - Ends: first and last N_ENDS bases
#  - GC content over A/C/G/T only (Ns ignored)
#  - Distinct k-mers: distinct canonical k-mers / all k-mers (canonical =
#    the smaller of a k-mer and its reverse complement, so inverted copies
#    count as repeats too). 100% = no internal repeat of length >= k.
#  - Repeat-covered bp: share of bases inside a k-mer that occurs more than
#    once (either strand); i.e. how much of the element is internal repeats.
N_ENDS <- 6
seq_stats <- function(s, kmers, k) {
  L <- nchar(s)
  ends <- if (L <= 2 * N_ENDS) s else
    paste0(substr(s, 1, N_ENDS), "...", substr(s, L - N_ENDS + 1, L))
  gc  <- stri_count_regex(s, "[GC]")
  acgt <- stri_count_regex(s, "[ACGT]")

  if (nrow(kmers)) {
    canon <- pmin(kmers$kmer, kmers$rcc)
    distinct <- uniqueN(canon) / length(canon)
    rep_pos <- kmers$pos[duplicated(canon) | duplicated(canon, fromLast = TRUE)]
    # mark rep_pos .. rep_pos + k - 1 as covered (difference array, O(L))
    d <- integer(L + 1)
    if (length(rep_pos)) {
      d <- d + tabulate(rep_pos, nbins = L + 1) -
        tabulate(rep_pos + k, nbins = L + k + 1)[seq_len(L + 1)]
    }
    covered <- sum(cumsum(d)[seq_len(L)] > 0) / L
  } else {
    distinct <- NA; covered <- NA
  }
  pct <- function(x, d = 1) if (is.na(x)) "NA" else sprintf(paste0("%.", d, "f%%"), 100 * x)
  list(
    list(label = "5' ... 3'",                         value = ends, mono = TRUE),
    list(label = "GC content",                   value = pct(if (acgt) gc / acgt else NA)),
    list(label = sprintf("Distinct %d-mers", k), value = pct(distinct)),
    list(label = "Repeat-covered bp",            value = pct(covered))
  )
}

### Triangle self-dotplot
# The upper half of the self-dotplot rotated 45 degrees. A match between
# positions i < j is drawn at x = (i + j) / 2 (midpoint of the two copies) and
# height (j - i) / 2, i.e. a true rotation, so the base of the triangle is the
# sequence itself and lines up with the ORF lanes below it. The y axis is
# labelled with the distance between the two copies (j - i).
# Direct repeats appear as horizontal lines, inverted repeats as vertical ones.
draw_triangle <- function(fwd, inv, L, orfs, title, subtitle, fname, info = NULL) {
  dots <- rbind(fwd[, type := "fwd"], inv[, type := "rc"])
  dots[, `:=`(x = (V1 + V2) / 2, y = (V2 - V1) / 2,
              type = factor(type, levels = c("fwd", "rc")))]
  setorder(dots, -type)            # draw forward matches on top

  # Only lanes (strand + frame) that contain at least one ORF, in the order
  # +1 +2 +3 -1 -2 -3 going down from the triangle
  all_lanes <- lane_names_for(c("+", "-"))
  orfs[, lab := paste0(strand, frame)]
  lanes <- all_lanes[all_lanes %in% orfs$lab]
  orfs[, lane := match(lab, lanes)]

  pad    <- 0.01 * L
  gap    <- 0.015 * L
  lane_h <- LANE_FRAC * L
  apex   <- (L - 1) / 2
  loX <- 1 - pad; hiX <- L + pad
  trackTop <- -gap
  loY <- trackTop - length(lanes) * lane_h - (if (length(lanes)) pad else 0)
  hiY <- apex + pad

  bp_mm    <- PANEL_IN * 25.4 / (hiX - loX)
  dot_size <- min(max(DOT_SIZE_MM, bp_mm), 3)
  dims     <- plot_dims(hiX - loX, hiY - loY)

  # Grid clipped to the triangle: vertical lines at positions, horizontal
  # lines at distances
  brksX <- nice_breaks(L)
  brksD <- nice_breaks(L - 1)
  brksD <- brksD[brksD > 0]
  grid <- rbind(
    data.table(x = brksX, xend = brksX, y = 0,
               yend = pmax(0, pmin(brksX - 1, L - brksX)))[yend > 0],
    data.table(x = 1 + brksD / 2, xend = L - brksD / 2,
               y = brksD / 2, yend = brksD / 2)
  )
  tri <- data.table(x = c(1, (1 + L) / 2, L), y = c(0, apex, 0))

  # Stats panel in the empty corner right of the triangle, top to bottom.
  # Spacing is set in inches and converted to data units, so it looks the
  # same whatever the sequence length.
  info_layers <- list()
  if (length(info)) {
    u_per_in <- (hiX - loX) / PANEL_IN
    step  <- 0.62 * u_per_in             # one label/value block
    gapLV <- 0.24 * u_per_in             # label -> value
    y0 <- hiY - 0.05 * u_per_in
    for (i in seq_along(info)) {
      yl <- y0 - (i - 1) * step
      info_layers <- c(info_layers, list(
        annotate("text", x = hiX, y = yl, label = info[[i]]$label,
                 hjust = 1, vjust = 1, size = 3.4, colour = "grey45"),
        annotate("text", x = hiX, y = yl - gapLV, label = info[[i]]$value,
                 hjust = 1, vjust = 1, size = 4.6, colour = "grey15",
                 fontface = "bold",
                 family = if (isTRUE(info[[i]]$mono)) "mono" else "")
      ))
    }
  }

  p <- ggplot(dots, aes(x = x, y = y, colour = type)) +
    geom_polygon(data = tri, aes(x = x, y = y), inherit.aes = FALSE,
                 fill = "white", colour = NA) +
    geom_segment(data = grid, aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE, colour = "grey92", linewidth = 0.3) +
    geom_point(shape = 16, size = dot_size, stroke = 0) +
    geom_polygon(data = tri, aes(x = x, y = y), inherit.aes = FALSE,
                 fill = NA, colour = "grey40", linewidth = 0.6) +
    track_frame(lanes, trackTop, lane_h, loX, hiX, "bottom", loX - pad) +
    orf_layers(orfs, trackTop, lane_h, L, "bottom") +
    info_layers +
    scale_colour_manual(
      values = c(fwd = COL_FWD, rc = COL_RC),
      labels = c(fwd = "Direct (same strand)",
                 rc  = "Inverted (reverse complement)"),
      drop = FALSE, name = NULL
    ) +
    scale_x_continuous(breaks = brksX, labels = scales::comma) +
    scale_y_continuous(breaks = brksD / 2,
                       labels = scales::comma(brksD)) +
    coord_fixed(ratio = 1, xlim = c(loX, hiX), ylim = c(loY, hiY),
                expand = FALSE, clip = "off") +
    labs(title = title, subtitle = subtitle, x = "Position (bp)",
         y = "Distance between copies (bp)") +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title          = element_text(face = "bold", size = 15),
      plot.subtitle       = element_text(colour = "grey35", size = 11,
                                         margin = margin(b = 10)),
      plot.title.position = "plot",
      axis.title          = element_text(size = 12, colour = "grey20"),
      axis.title.y        = element_text(hjust = (apex / 2 - loY) / (hiY - loY)),
      axis.text           = element_text(size = 11, colour = "grey20"),
      axis.ticks          = element_line(colour = "grey40", linewidth = 0.4),
      panel.grid          = element_blank(),
      legend.position     = "bottom",
      legend.text         = element_text(size = 11),
      plot.margin         = margin(12, 16, 8, 12)
    )

  ggsave(fname, plot = p, width = dims$width, height = dims$height,
         units = "in", dpi = DPI, bg = "white")
}

### Self-dotplot (one sequence)
get_plot <- function(seq) {
  name <- gsub(">", "", seq$name)
  print(name)
  s <- toupper(seq$seq)
  kmers1 <- get_kmers(s, kmer_size)
  kmers1[, rcc := rc(kmer)]

  fwd <- match_pairs(kmers1, "kmer")   # direct repeats
  inv <- match_pairs(kmers1, "rcc")    # inverted repeats
  orfs   <- find_orfs(s)
  seqlen <- nchar(s)

  if (SELF_STYLE == "triangle") {
    n_orf <- nrow(orfs)
    draw_triangle(
      fwd, inv, seqlen, orfs,
      title    = wrap(name),
      subtitle = sprintf("%s bp  |  k = %d  |  %s%s",
                         scales::comma(seqlen), kmer_size, orf_note(),
                         if (n_orf) "" else ": none found"),
      fname = paste0(safe_name(name), ".png"),
      info  = seq_stats(s, kmers1, kmer_size)
    )
    return(invisible(NULL))
  }

  # Square style: mirror across the diagonal so the plot is symmetric
  fwd[, type := "fwd"]
  inv[, type := "rc"]
  dots <- rbind(fwd, fwd[, .(V1 = V2, V2 = V1, type)],
                inv, inv[, .(V1 = V2, V2 = V1, type)])

  draw_dotplot(
    dots, seqlen, seqlen,
    title    = wrap(name),
    subtitle = sprintf("%s bp  |  k = %d  |  %s: + strand on top, - strand on the right",
                       scales::comma(seqlen), kmer_size, orf_note()),
    xlab = "Position (bp)", ylab = "Position (bp)",
    fname = paste0(safe_name(name), ".png"),
    top_orfs   = assign_lanes(orfs, "+"), top_lanes   = lane_names_for("+"),
    right_orfs = assign_lanes(orfs, "-"), right_lanes = lane_names_for("-"),
    diagonal = TRUE
  )
}

### Pairwise dotplot (two sequences)
# The two axes are different sequences, so each track shows all ORFs of the
# sequence on that axis: top = sequence A (x), right = sequence B (y),
# with + strand frames next to the dotplot and - strand frames outside.
get_plot_pair <- function(seqA, seqB) {
  nameA <- gsub(">", "", seqA$name)
  nameB <- gsub(">", "", seqB$name)
  print(paste("Pairwise:", nameA, "vs", nameB))
  sA <- toupper(seqA$seq)
  sB <- toupper(seqB$seq)

  kmersA <- get_kmers(sA, kmer_size)
  kmersB <- get_kmers(sB, kmer_size)
  kmersA[, rcc := rc(kmer)]

  dots <- rbind(
    cross_pairs(kmersA, kmersB, "kmer")[, type := "fwd"],
    cross_pairs(kmersA, kmersB, "rcc")[,  type := "rc"]
  )

  lenX <- nchar(sA)
  lenY <- nchar(sB)
  both <- c("+", "-")
  draw_dotplot(
    dots, lenX, lenY,
    title    = wrap(paste0(nameA, "  vs  ", nameB)),
    subtitle = sprintf("x: %s bp  |  y: %s bp  |  k = %d  |  %s of each sequence along its axis",
                       scales::comma(lenX), scales::comma(lenY), kmer_size, orf_note()),
    xlab = wrap(paste0(nameA, " (bp)"), 60),
    ylab = wrap(paste0(nameB, " (bp)"), 60),
    fname = paste0(safe_name(nameA, 60), "_vs_", safe_name(nameB, 60), ".png"),
    top_orfs   = assign_lanes(find_orfs(sA), both), top_lanes   = lane_names_for(both),
    right_orfs = assign_lanes(find_orfs(sB), both), right_lanes = lane_names_for(both)
  )
}

### Main
args <- commandArgs(trailingOnly = TRUE)

if (!length(args) %in% c(3, 4)) {
  print("usage:  Rscript TEplot.R fasta_file kmer_size  output_folder")
 # print("usage (pairwise dotplot): dtp.R fasta_file1 fasta_file2 kmer_size output_folder")
  quit()
}

two_seqs <- length(args) == 4

if (two_seqs) {
  kmer_size <- as.integer(args[3])
  outputf   <- args[4]
  seqsA     <- ffasta(args[1])
  seqsB     <- ffasta(args[2])
} else {
  kmer_size <- as.integer(args[2])
  outputf   <- args[3]
  seqs      <- ffasta(args[1])
}

dir.create(outputf, showWarnings = FALSE, recursive = TRUE)
setwd(outputf)

if (two_seqs) {
  for (i in seq_len(nrow(seqsA))) {
    for (j in seq_len(nrow(seqsB))) {
      get_plot_pair(seqsA[i], seqsB[j])
    }
  }
} else {
  invisible(lapply(seq_len(nrow(seqs)), function(x) get_plot(seqs[x])))
}

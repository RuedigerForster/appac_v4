
column_names <- c("Sample_Name", "Peak_Name", "Injection_Date",
                  "Air_Pressure", "Raw_Area")

## Maximum per-peak NA fraction a sample may carry: .prepare() drops peaks above
## it, and appac()'s degeneracy guard only flags series at or below it -- both
## MUST use this single value so the guard matches which peaks are actually
## kept.
max_na_fraction <- 0.3

## Okabe-Ito colourblind-safe colours (Okabe & Ito 2008,
## https://jfly.uni-koeln.de/color/).  Keys keep their semantic roles:
##   highlight = corrected / emphasised series, lowlight = raw / background,
##   line = fit & changepoint lines, fill = confidence-band fill.
default_palette <- list(
  highlight_color = "#0072B2",  # blue
  lowlight_color  = "#999999",  # grey
  line_color      = "#D55E00",  # vermillion
  fill_color      = "#56B4E9"   # sky blue
)

## Configurable plot GEOMETRY (point/line sizes, alpha, bins, ...).  Together
## with default_palette (colours) and the .*_plot_theme() functions
## (text/layout)
## this keeps every visual choice in one place; the plot functions read from
## here
## rather than hard-coding values.
default_aes <- list(
  point_size            = 0.6,
  point_alpha           = 0.5,
  reference_linewidth   = 0.4,      # reference / centre / zero line
  fit_linewidth         = 0.6,      # fitted curve / normal overlay / Q-Q line
  changepoint_linewidth = 0.5,
  changepoint_linetype  = "dotted",
  histogram_bins        = 60,
  legend_key_size       = 3,        # legend point size, so the key dots read
  legend_text_rel       = 1.2       # legend text size, relative to base
)

## Okabe-Ito categorical palette, for colouring multiple series (e.g. peaks).
## Ordered so the first colours are high-contrast on a grey panel (the low-
## contrast yellow and black come last), while staying colourblind-safe.
okabe_ito_palette <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2",
                       "#D55E00", "#CC79A7", "#F0E442", "#000000")

## Themes are built on theme_grey(base_size) so all text scales sensibly from a
## single base font size (points); `size` is that base size.
.main_plot_theme <- function(size = 12) {
  ggplot2::theme_grey(base_size = size) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold",
                                              size = ggplot2::rel(1.25)),
      plot.subtitle   = ggplot2::element_text(colour = "grey35"),
      legend.position = "none"
    )
}

.residuals_plot_theme <- function(size = 12) {
  ggplot2::theme_grey(base_size = size) +
    ggplot2::theme(legend.position = "none")
}

.histogram_plot_theme <- function(size = 12) {
  ggplot2::theme_grey(base_size = size) +
    ggplot2::theme(
      axis.title.y    = ggplot2::element_blank(),
      axis.text.y     = ggplot2::element_blank(),
      axis.ticks.y    = ggplot2::element_blank(),
      legend.position = "none"
    )
}

## column names referenced inside aes() are not globals; declare them so
## R CMD check does not flag "no visible binding for global variable".
utils::globalVariables(c(
  "area", "series", "pressure", "date", "residual", "fitted", "theoretical",
  "sample", "deviation", "peak", ".covar"
))

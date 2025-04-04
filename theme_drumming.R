theme_drumming <- function() {
  font <- "Helvetica"
  theme <- list(
    ggplot2::theme(
      ## TITLES
      plot.title = ggplot2::element_text(
        size = 20,
        face = "bold",
        color = "black",
        family = font
      ),
      plot.subtitle = ggplot2::element_text(
        size = 16,
        face = "bold",
        color = "black",
        family = font
      ),
      plot.caption = ggplot2::element_text(
        size = 12,
        face = "italic",
        color = "black",
        family = font
      ),


      ## LEGEND
      legend.title = ggplot2::element_text(family = font, size = 16, color = "#222222"),
      #legend.title.align = 0,
      legend.key = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(
        family = font, size = 16,
        color = "#222222"
      ),
      #legend.position = c(.98, .99),
      #legend.justification = c("right", "top"),
      #legend.box.just = "right",


      ## AXIS
      axis.title = ggplot2::element_text(family = font, size = 16, color = "#222222"),
      axis.text = ggplot2::element_text(family = font, size = 12, color = "#222222"),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),

      ## PANELS
      panel.grid.major.y = ggplot2::element_line(color = "#cbcbcb"),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing.y = ggplot2::unit(1, "lines"),

      ## FACETS
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        family = font, size = 16,
        face = "bold", color = "#222222"
      ),
      strip.placement = "inside"
    ) # end of ggplot2::theme
    ,
    ggplot2::guides(
      colour = guide_legend(label.position = "left"), # second item in list
      fill = guide_legend(label.position = "left")
    ), # third item in list
    scico::scale_colour_scico_d(palette = "roma", begin=0, end=0.7), # fourth item in list
    scico::scale_fill_scico_d(palette = "roma", begin=0, end=0.7) # fifth item in list
  )

  return(theme)
}

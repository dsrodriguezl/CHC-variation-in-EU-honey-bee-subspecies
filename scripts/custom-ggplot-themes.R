# Script to generate custom ggplot settings

library(ggplot2)
library(ggtext)

# Theme for cluster analysis dendrogram plots ----
cluster_theme <- theme(legend.position = "left"
                       , legend.key = element_rect(fill = "white"
                                                   , color = "white")
                       , legend.key.size = unit(3.5, "mm")
                       , legend.spacing.y = unit(0.2, "lines")
                       , legend.title = element_text(size = 10)
                       , legend.text = element_text(size = 8)
                       , panel.background = element_rect(fill = "white")
                       , plot.margin = unit(c(1, 10, 1, 1), "lines")
                       , axis.text.x = element_text(size = 7)
                       , axis.text.y = element_blank()
                       , axis.ticks = element_blank()
                       , axis.title = element_blank())

cluster_theme2 <- theme(legend.position = "left"
                        , legend.key = element_rect(fill = "white"
                                                    , color = "white")
                        , legend.key.size = unit(3.5, "mm")
                        , legend.spacing.y = unit(0.2, "lines")
                        , legend.title = element_text(size = 10)
                        , legend.text = element_text(size = 8)
                        , panel.background = element_rect(fill = "white")
                        , plot.margin = unit(c(1, 5, 1, 1), "lines")
                        , axis.text.x = element_text(size = 7)
                        , axis.text.y = element_blank()
                        , axis.ticks = element_blank()
                        , axis.title = element_blank())

# Theme for NMDS plots ----
nmds_theme <- theme_classic() +
  theme(strip.text = element_markdown(face = "bold"
                                             , size = 10)
        , strip.background = element_blank()
        #, legend.position = "left"
        , legend.key.size = unit(3.5, "mm")
        , legend.spacing.y = unit(0.2, "lines")
        , legend.title = element_text(size = 10)
        , panel.border = element_rect(fill = NA
                                      , color = "black"
                                      , linetype = "dotted"
                                      , size = 0.5))

# Theme for box plots ----
boxplot_theme <- theme_classic() +
  theme(title = element_text(face = "bold"
                             , size = 11)
        , strip.text = element_text(face = "bold"
                                    , size = 10)
        , strip.background = element_blank()
        , legend.position = "none"
        # , legend.justification = "center"
        , legend.key.size = unit(5, "mm")
        , legend.spacing.y = unit(0.2, "lines")
        , legend.title = element_text(face = "plain"
                                      , size = 10)
        , legend.text = element_text(size = 8)
        , axis.text = element_text(face = "plain"
                                   , size = 8)
        , axis.title = element_text(face = "plain"
                                    , size = 10)
        # , axis.title.y = element_blank()
        # , axis.text.y = element_blank()
        # , axis.ticks.y = element_blank()
        )
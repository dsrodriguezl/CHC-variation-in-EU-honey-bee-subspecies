# Load packages
packages <- c("kableExtra"
              , "flextable"
              , "equatags"
              , "here"
              , "ggpubr"
              , "patchwork"
              , "ggtext"
              , "gghalves"
              , "ggh4x"
              , "ragg"
              , "vegan"
              , "dunn.test"
              , "ggdendro"
              , "magrittr"
              , "tidyverse")
pacman::p_load(char = packages)

# Custom functions
source(here("scripts", "functions.R"))
# Ggplot settings
source(here("scripts", "custom-ggplot-themes.R"))

# If climatic data has already been downloaded in a previous session, avoid 
# re-downloading it
if(list.files(here(), pattern = "worldclim_tables.Rdata") |> length() == 0) {
  source(here("scripts", "worldclim_data.R"))
}

# Get Composition data
master_table <- read.csv("composition_data.csv"
                         , row.names = 1)

# Separate data frame with composition of samples
datenc <- master_table |> 
  select(-(Class:RI)) |> 
  t() |> 
  as.data.frame()

# Separate data frame with information of compounds
Comps_data <- master_table |> 
  select(Class:RI)

# Samples group information
grouping_info <- read.csv("samples_list.csv"
                         , row.names = 1
                         , colClasses = "character") |> 
  mutate_all(as.factor)


# NMDS ----
set.seed(92)
# Use samples' composition
sol <- datenc |>
  # Calculate Bray-Curtis dissimilarity between samples
  vegdist(method = "bray") |> 
  # Perform 2D NMDS
  metaMDS(try = 200
          , k = 2
          , trymax = 400
          , autotransform = F)

# Extract NMDS coordinates of samples and merge it in a data frame with their
# group information
NMDS <- scores(sol, display = "sites") |> 
  as.data.frame() |> 
  rownames_to_column("rowname") |> 
  merge(grouping_info |> 
          rownames_to_column("rowname")
        , all = T
        , sort = F) |> 
  column_to_rownames("rowname")

# Plot NMDS
nmds_plot <- NMDS |>  
  ggplot(aes(x = NMDS1, y = NMDS2
             , fill = Subspecies)) +
  geom_point(data = transform(NMDS
                              , Task = NULL)
             , color = "grey85"
               , fill = "grey85"
               , size = 2.5) +
  geom_point(shape = 21
             , color = "grey40"
             , size = 2.5) +
  scale_fill_viridis_d(guide =
                         guide_legend(label.theme =
                                        element_text(face = "italic"
                                                     , size = 8)
                                      , byrow = T
                                      # Force the shape in the legend to
                                      # include filling, so it can show
                                      # the color of the fill
                                      , override.aes = list(shape = 21
                                                            , size = 3.2))) +
  facet_wrap(~Task
             , nrow = 2) +
  coord_equal() +
  nmds_theme

# Export plot
ggsave(here("figs"
            ,"nmds_sbsp_splitted.png")
       , plot = nmds_plot
       , dpi = "retina"
       , width = 4
       , height = 6
       , units = "in"
       , device = agg_png
       , scaling = 1)

# PERMANOVA ----
## Set up permutations
### Number of permutations
PERM <- how(nperm = 1000)
permTsk <- PERM
permSbsp <- PERM

# This defines the permutation blocks as the different subspecies for the 
# comparisons between task performance groups
setBlocks(permTsk) <- with(datenc, grouping_info$Task)
setBlocks(permSbsp) <- with(datenc, grouping_info$Subspecies)

## permutational dispersion test ----

### Dispersion measurement
task_dispersion <- vegdist(datenc
                           , method = "bray") |> 
  betadisper(grouping_info$Task) 

### Dispersion test
set.seed(92)
dispperm_table <- task_dispersion |> 
  permutest(permutations = permSbsp) |> 
  permutest_table(factor = "Nurses - Foragers")

## PERMANOVA test ----
set.seed(92)
PMVSsp <- adonis2(vegdist(datenc
                          , method = "bray") ~ Subspecies
                  , data = grouping_info
                  , permutations = permTsk
                  )

### Getting table to report PERMANIVA results
PMVssp_table <- permanova_table(PMVSsp)

### Pair-wise PERMANOVA
pwpmw_table <- grouping_info |> 
  rownames_to_column("rowname") |> 
  # Merge samples group data with samples compositional data
  merge(datenc |> 
          rownames_to_column("rowname")
        , all = T
        , sort = F) |> 
  column_to_rownames("rowname") |> 
  group_by(Task) |> 
  # Separate forager and nurses
  group_split() |> 
  # Name list items after the group they represent
  set_names(levels(grouping_info$Task)) |> 
  # Iterate through data frames in the list (i.e. Tasks)
  map(function(dfn) {
    # list of all possible subspecies pairs
    pairs <- dfn$Subspecies |>
      levels() |>
      combn(2, simplify = F)
    
    # Collapse pairs in a single string to have subspecies names separated by 
    # "vs" as the names of each entry in the list
    names(pairs) <- pairs |> 
      map(function(x){
        paste(x, collapse = " vs ")
      })
    
    pairs |> 
      # Iterate through pairs
      imap(function(pair, pair_id){
        
        # Data frame with only samples corresponding to the subspecies of 
        # the current pair
        dfn <- dfn |> 
          filter(Subspecies == pair[1] | Subspecies == pair[2]) |> 
          as.data.frame()
        
        # Data frame with group information in dfn
        groups <- dfn |> 
          select(Task:Individual) |> 
          as.data.frame()
        
        # Data frame with compositional data in dfn
        tmp <- dfn |> 
          select(-(Task:Individual)) |> 
          as.data.frame()
        
        # PERMANOVA between current pair
        pmv_table <- withr::with_seed(12345
                                      , adonis2(vegdist(tmp
                                                        , method = 
                                                          "bray") ~ Subspecies
                                                , data = groups
                                                , permutations = PERM
                                      )) |>
          # Reportable PERMANOVA results table
          permanova_table() |> 
          # Get adjusted p value for current pair
          mutate(`adjusted p` = p.adjust(p, method = "fdr", n = length(pairs))
                 , contrast = pair_id
          ) |> 
          rownames_to_column("rowname") |> 
          # Remove redundant row, the total can be directly inferred from
          # the sum of the previous rows
          filter(rowname != "Total")
        
        # Get new name of "F" column
        f <- paste0("F("
                    , paste(pmv_table$Df
                            , collapse = ", ")
                    , ")")
        # Reorder columns
        pmv_table <- pmv_table |> 
          select(contrast, SS:`adjusted p`)
        
        # Replace name of "F" column
        colnames(pmv_table) <- pmv_table |> 
          colnames() |> 
          str_replace("F"
                      , f)
        # PERMANOVE results table for pair of current iteration
        pmv_table
      }) |>  
      # Reduce list of data frames into a single data frame
      reduce(merge
             , all = T
             , sort = F)
  }) |>  
  # Add column indicating Task group to each data frame
  imap(function(df, name) {
    df |> 
      cbind(name)
  }) |> 
  # Reduce list of data frames to a single data frame
  reduce(merge
         , all = T
         , sort = F) |> 
  # Reorder columns
  (function(df){
    df |> 
      set_colnames(c(df |> 
                       select(contrast:`adjusted p`) |> 
                       colnames()
                     , "Task")) |> 
      select(Task, everything())
  })()
  
## Cluster analysis ----


datenc2 <- datenc |> 
  # Merge samples composition with their group information
  merge(grouping_info
        , by = "row.names"
        , all = T) |> 
  # Remove unnecesary columns
  select(-Row.names, -Individual) |> 
  group_by(Task, Subspecies) |>
  # Calculate mean abundance of compounds for each task group of each subspecies
  summarise(across(everything(), mean)) |> 
  ungroup() |> 
  # Column specifying group
  mutate(group = as.factor(paste(Subspecies, Task, sep =  " - "))) |> 
  # Reorder columns
  select(group, everything()) |> 
  as.data.frame()

# Indicate group in row names
row.names(datenc2) <- datenc2$group

# Calculate Bray-Curtis dissimilarity between groups
ddata_all <- vegan::vegdist(datenc2 |> 
                                 select(!group:Subspecies)
                               , method = "bray") |> 
  # Perform hierarchical clustering
  hclust(method  = "ward.D2") |>  
  # Extract data for plotting
  as.dendrogram() |> 
  dendro_data()

# Get labels
dlabs_all <- label(ddata_all)
# Set labels as row names
row.names(dlabs_all) <- dlabs_all$label
# Merge labels with grouping information of each factor
dlabs_all <-  merge(dlabs_all
                       , datenc2 |> 
                      select(group:Subspecies)
                       , by = "row.names"
                       , all = T)

# PLot dendrogram of the hierarchical clustering
all_plot <- ggplot(segment(ddata_all)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = dlabs_all
            , aes(label = Subspecies
                  , x = x
                  , y = y - 0.025
                  , hjust = "left"
                  , fontface = "italic"
                  # , color = Subspecies
            )
            # Keep the text of the same size as legend labels' text
            , size = cluster_theme$legend.text$size * 1.3 / .pt
            , alpha = 0.9) +
  geom_point(data = dlabs_all
             , aes(shape = Task
                   , fill = Subspecies
                   , x = x
                   , y = y - 0.01)
             , size = 2.8) +
  scale_y_reverse() +
  scale_shape_manual(values = c("Foragers" = 21
                                , "Nurses" = 24)) +
  scale_fill_viridis_d(guide = "none") +
  coord_flip(clip = "off") +
  cluster_theme2 + 
  theme(plot.title = element_text(hjust = 0)
        , plot.title.position = "plot"
        , legend.text = element_text(size = cluster_theme$legend.text$size * 1.2)
        , legend.key.size = cluster_theme$legend.key.size * 1.2)

# Export plot
ggsave(here("figs"
            ,"cluster_sbsp.png")
       , plot = all_plot
       , dpi = "print"
       , width = 6
       , height = 4
       , units = "in"
       , device = agg_png
       , scaling = 0.8)

# Analysis by hydrocarbon classes ----
## Data by hydrocarbon classes 
Prop_CompsClass <- cclasses_df(datenc, grouping_info, Comps_data
                               , fuse.methyls = T)

### Hydrocarbons' abundance among groups
Prop_CompsClass_long <- pivot_longer(Prop_CompsClass
                                     , cols = !where(is.factor)
                                     , names_to = "CClass"
                                     , values_to = "abundance")

# Test difference between subspecies for nurses and foragers separately
nurses_Sbsp_cclassestest <- cclasses_test(Prop_CompsClass_long |> 
                                        filter(Task == "Nurses")
                                      , faktor = "Subspecies"
                                      , test = "kw")

foragers_Sbsp_cclassestest <- cclasses_test(Prop_CompsClass_long |> 
                                        filter(Task == "Foragers")
                                      , faktor = "Subspecies"
                                      , test = "kw")

# Correct name for methyl alkanes
Prop_CompsClass_long[
  , "CClass"][Prop_CompsClass_long[
    , "CClass"] == "Methyls"] <- "Methyl alkanes"

# Data frame to set axis scales for nurses plots
scales_nurses <- {list(scale_y_continuous(limits = c(20, 80)
                                         , breaks = seq(20, 80, by = 10))
                      , scale_y_continuous(limits = c(10, 70)
                                           , breaks = seq(10, 70, by = 10))
                      , scale_y_continuous(limits = c(0, 20)
                                           , breaks = seq(0, 20, by = 5))
                      , scale_y_continuous(limits = c(0, 40)
                                           , breaks = seq(0, 40, by = 5)))}

# Data frame to set anotations of significative differences between subspecies
# among nurses
anotation_nurses <- {data.frame(CClass = c(rep("Alkanes"
                                         , times =
                                           Prop_CompsClass_long$Subspecies |> 
                                           levels() |> 
                                           length())
                                     , rep("Alkenes"
                                           , times =
                                             Prop_CompsClass_long$Subspecies |> 
                                             levels() |> 
                                             length())
                                     , rep("Alkadienes"
                                           , times =
                                             Prop_CompsClass_long$Subspecies |> 
                                             levels() |> 
                                             length())
                                     , rep("Methyl alkanes"
                                           , times =
                                             Prop_CompsClass_long$Subspecies |> 
                                             levels() |> 
                                             length()))
                          , Subspecies = rep(Prop_CompsClass_long$Subspecies |> 
                                               levels()
                                             , times = Prop_CompsClass_long$CClass |> 
                                               unique() |> 
                                               length()) |>  as.factor()
                          , label = c(c("A", "AB", "AB", "B", "A", "A")
                                      , c("A", "A", "A", "A", "A", "A")
                                      , c("A", "BC", "BC", "B", "BC", "AC")
                                      , c("A", "AB", "AB", "B", "AB", "AB"))) |> 
  as_tibble()}

# Specifying coordinates of each anotation
anotation_nurses <- {anotation_nurses |> 
  mutate(abundance = 
           c(ifelse(anotation_nurses$CClass == "Alkanes"
                    , "80"
                    , ifelse(anotation_nurses$CClass == "Alkenes"
                             , "70"
                             , ifelse(anotation_nurses$CClass == "Alkadienes"
                                      , "20"
                                      , "40")))) |> 
           as.numeric())}

# Plot differences in abundance of eahc hydrocarbon classes for the nurses
cclasses_nurses <- {ggplot(Prop_CompsClass_long |> 
                             filter(Task == "Nurses")
                           , aes(x = Subspecies 
                                 , y = abundance
                                 , fill = Subspecies)) +
    geom_half_boxplot(color = "grey30"
                      , size = 0.4
                      , alpha = 0.9
                      , side = "r"
                      #, center = T
                      , outlier.shape = NA
                      , errorbar.draw = F
                      , nudge = 0.0) +
    geom_half_point(color = "grey40"
                    , side = "l"
                    , shape = 21
                    , size = 0.75
                    , stroke = 0.2
                    , alpha = 0.9
                    , transformation = position_jitter(seed = 12345
                                                       , width = 0.08
                                                       , height = 0)) +
    scale_fill_viridis_d() +
    scale_color_viridis_d() +
    facet_wrap(vars(factor(CClass
                           , levels = c("Alkanes"
                                        , "Alkenes"
                                        , "Alkadienes"
                                        , "Methyl alkanes")))
               , scales = "free_x") +
    # Set axes scales
    facetted_pos_scales(y = scales_nurses) +
    coord_flip() +
    boxplot_theme +
    theme(panel.border = element_rect(fill = NA
                                      , color = "black"
                                      , linetype = "dotted"
                                      , size = 0.5)
          , legend.text = element_text(size = boxplot_theme$legend.text$size
                                       , face = "italic")
          , axis.text.y = element_text(face = "italic")) +
    labs(y = "Relative abundance (%)"
         , subtitle = "Nurses") +
    # Add significance labels
    geom_text(anotation_nurses
              , mapping = aes(x = Subspecies
                              , y = abundance
                              , label = label)
              , size = boxplot_theme$legend.text$size / .pt
              , fontface = "bold")}

# Set axis scales for foragers plots
scales_foragers <- {list(scale_y_continuous(limits = c(40, 90)
                                           , breaks = seq(40, 90, by = 10))
                        , scale_y_continuous(limits = c(10, 60)
                                             , breaks = seq(10, 60, by = 10))
                        , scale_y_continuous(limits = c(0, 25)
                                             , breaks = seq(0, 25, by = 5))
                        , scale_y_continuous(limits = c(0, 10)
                                             , breaks = seq(0, 10, by = 2)))}

# Data frame to set anotations of significative differences between subspecies
# among foragers
anotation_foragers <- {data.frame(CClass = 
                                    c(rep("Alkanes"
                                          , times =
                                            Prop_CompsClass_long$Subspecies |> 
                                            levels() |> 
                                            length())
                                      , rep("Alkenes"
                                            , times =
                                              Prop_CompsClass_long$Subspecies |> 
                                              levels() |> 
                                              length())
                                      , rep("Alkadienes"
                                            , times =
                                              Prop_CompsClass_long$Subspecies |> 
                                              levels() |> 
                                              length())
                                      , rep("Methyl alkanes"
                                            , times =
                                              Prop_CompsClass_long$Subspecies |> 
                                              levels() |> 
                                              length()))
                                  , Subspecies = 
                                    rep(Prop_CompsClass_long$Subspecies |> 
                                          levels()
                                        , times = Prop_CompsClass_long$CClass |> 
                                          unique() |> 
                                          length()) |>  as.factor()
                                  , label = 
                                    c(c("A", "B", "A", "AB", "AB", "AB")
                                      , c("A", "A", "A", "A", "A", "A")
                                      , c("A", "B", "B", "B", "B", "B")
                                      , c("A", "B", "AC", "C", "BC", "C"))) |> 
    as_tibble()}

# Specifying coordinates of each anotation
anotation_foragers <- {anotation_foragers |> 
    mutate(abundance = 
             c(ifelse(anotation_nurses$CClass == "Alkanes"
                      , "90"
                      , ifelse(anotation_nurses$CClass == "Alkenes"
                               , "60"
                               , ifelse(anotation_nurses$CClass == "Alkadienes"
                                        , "25"
                                        , "10")))) |> 
             as.numeric())}

# Plot differences in abundance of eahc hydrocarbon classes for the nurses
cclasses_foragers <- {ggplot(Prop_CompsClass_long |> 
                            filter(Task == "Foragers")
                          , aes(x = Subspecies 
                                , y = abundance
                                , fill = Subspecies)) +
  geom_half_boxplot(color = "grey30"
                    , size = 0.4
                    , alpha = 0.9
                    , side = "r"
                    , outlier.shape = NA
                    , errorbar.draw = F
                    , nudge = 0.0) +
  geom_half_point(color = "grey40"
                  , side = "l"
                  , shape = 21
                  , size = 0.75
                  , stroke = 0.2
                  , alpha = 0.9
                  , transformation = position_jitter(seed = 12345
                                                     , width = 0.08
                                                     , height = 0)) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  facet_wrap(vars(factor(CClass
                         , levels = c("Alkanes"
                                      , "Alkenes"
                                      , "Alkadienes"
                                      , "Methyl alkanes")))
             , scales = "free_x") +
    # Set axes scales
  facetted_pos_scales(y = scales_foragers) +
  coord_flip() +
  boxplot_theme +
  theme(panel.border = element_rect(fill = NA
                                    , color = "black"
                                    , linetype = "dotted"
                                    , size = 0.5)
        , legend.text = element_text(size = boxplot_theme$legend.text$size
                                     , face = "italic")
        , axis.text.y = element_text(face = "italic")) +
  labs(y = "Relative abundance (%)"
       , subtitle = "Foragers") +
    # Add significance labels
geom_text(anotation_foragers
          , mapping = aes(x = Subspecies
                          , y = abundance
                          , label = label)
          , size = boxplot_theme$legend.text$size / .pt
          , fontface = "bold")}

# Merge plots of nurses and forager into a single compound plot
cclasses_plot <- cclasses_foragers / cclasses_nurses

# Export plot
ggsave(here("figs"
            ,"cclasses_sbsp_splitted.png")
       , plot = cclasses_plot
       , dpi = "retina"
       , width = 5.5
       , height = 7.5
       , units = "in"
       , device = agg_png
       , scaling = 1)

# Chain length analysis ----
# Get Weighted mean chain length of each sample
Prop_chain.length <- cl_df(datenc, Comps_data, grouping_info)

# Test difference in weighted mean chain length between subspecies for nurses
# and foragers separately
nurses_Sbsp_cltest <- cl_test(Prop_chain.length |> 
                            filter(Task == "Nurses")
                          , faktor = "Subspecies", test.type = "kw")

foragers_Sbsp_cltest <- cl_test(Prop_chain.length |> 
                            filter(Task == "Foragers")
                          , faktor = "Subspecies", test.type = "kw")

# Plot difference in mean chain length for nurses
clength_nurses <- ggplot(Prop_chain.length |> 
         filter(Task == "Nurses")
       , aes(x = Subspecies
             , y = cl_w.mean
             , fill = Subspecies)) +
  geom_half_boxplot(color = "grey30"
                    , size = 0.4
                    , alpha = 0.9
                    , side = "r"
                    , outlier.shape = NA
                    , errorbar.draw = F
                    , nudge = 0.0) +
  geom_half_point(color = "grey40"
                  , side = "l"
                  , shape = 21
                  , size = 0.75
                  , stroke = 0.8
                  , alpha = 0.9
                  , range_scale = 0.45
                  , transformation = position_jitter(seed = 12345
                                                     , width = 0.05
                                                     , height = 0)) +
  ylim(26, 32) +
  scale_fill_viridis_d() +
  coord_flip() +
  theme_classic() +
  boxplot_theme +
  theme(legend.text = element_text(size = boxplot_theme$legend.text$size
                                   , face = "italic")
        , axis.text.y = element_text(face = "italic")) +
  labs(y = "Mean CHCs length"
       , subtitle = "Nurses") +
  # Add significance labels
  annotate("text"
           , label = c("AB", "A", "AB", "B", "A", "AB")
           , x = levels(Prop_chain.length$Subspecies)
           , y = 32
           , size = boxplot_theme$legend.text$size / .pt
           , fontface = "bold")

# PLot difference in mean chain length for foragers
clength_foragers <- ggplot(Prop_chain.length |> 
                           filter(Task == "Foragers")
                         , aes(x = Subspecies
                               , y = cl_w.mean
                               , fill = Subspecies)) +
  geom_half_boxplot(color = "grey30"
                    , size = 0.4
                    , alpha = 0.9
                    , side = "r"
                    , outlier.shape = NA
                    , errorbar.draw = F
                    , nudge = 0.0) +
  geom_half_point(color = "grey40"
                  , side = "l"
                  , shape = 21
                  , size = 0.75
                  , stroke = 0.8
                  , alpha = 0.9
                  , range_scale = 0.45
                  , transformation = position_jitter(seed = 12345
                                                     , width = 0.05
                                                     , height = 0)) +
  ylim(25, 30) +
  scale_fill_viridis_d() +
  coord_flip() +
  theme_classic() +
  boxplot_theme +
  theme(legend.text = element_text(size = boxplot_theme$legend.text$size
                                   , face = "italic")
        , axis.text.y = element_text(face = "italic")) +
  labs(y = "Mean CHCs length"
       , subtitle = "Foragers") +
  # Add significance labels
  annotate("text"
           , label = c("AB", "AB", "AB", "A", "B", "AB")
           , x = levels(Prop_chain.length$Subspecies)
           , y = 30
           , size = boxplot_theme$legend.text$size / .pt
           , fontface = "bold")

# Merge plots of nurses and foragers in a single compound plot
clength_plot <- clength_foragers / clength_nurses

# Export plot
ggsave(here("figs"
            ,"clength_sbsp_splitted.png")
       , plot = clength_plot
       , dpi = "retina"
       , width = 6
       , height = 4.5
       , units = "in"
       , device = agg_png
       , scaling = 1)


# climate vs CHC ----
# Get climatic data
load(here("worldclim_tables.Rdata"))

# FUcntion to assign each subspecies to the corresponding country
sbsp_land_match <- function(country) {
  case_when(country == "GERMANY" ~ "A. m. carnica"
            , country == "PORTUGAL" ~ "A. m. iberiensis"
            , country == "ITALY" ~ "A. m. ligustica"
            , country == "BELGIUM" ~ "A. m. mellifera"
            , country == "GREECE" ~ "A. m. macedonica"
            , country == "MALTA" ~ "A. m. ruttneri")
}

# Function to perform Spearmann's correlation tests for climatic variables and 
# CHC composition
spearman_cor_clim_chc <- function(df) {
  df |> 
    group_by(Task) |> 
    # Separate nurses and foragers
    group_split() |> 
    # Iterate through task data frames
    lapply(function(x) {
      # Perform Spearman's correlation test
      x |> 
        pull(4) |>  
        cor.test(x |> 
                   pull(5)
                 , method = "spearman") |>
        # Shape data frame
        (function(l) {
          cbind.data.frame(estimate = l |> 
                             pluck("estimate")
                           , p_value = l |> 
                             pluck("p.value")) |> 
            t() |> 
            as.data.frame()
        })()
    })  |> 
    # Reduce list of data frames to a single data frame
    list_cbind() |> 
    # Specify task group
    set_colnames(unique(df$Task)) |>
    # Set up decimals
    mutate_all(function(x) {
      x |> 
        round(digits = 3) |> 
        format(digits = 3
               , nsmall = 3)
    }) |> 
    # Shape data frame
    rownames_to_column("stat") |> 
    pivot_longer(!stat
                 , names_to = "Task") |> 
    pivot_wider(names_from = "stat"
                , values_from = "value") |> 
    # Format significance reporting
    mutate(signif = case_when(p_value >= 0.05 ~ "n.s"
                              , p_value < 0.05 & p_value > 0.01 ~ "*"
                              , p_value <= 0.01 & p_value > 0.001 ~ "**"
                              , p_value <= 0.001 ~ "***"))
}

# Correlation between olefins to n-alkanes ratio and temperature
olefins_temp <- Prop_CompsClass |> 
  # Calculate olefins abundance by summing abundance of alkenes and alkadienes
  # Calculate olefins to n-alkanes ratio
  mutate(olefins = Alkenes + Alkadienes
         , o_2_alka = olefins / Alkanes) |> 
  # Reorder columns
  select(Task:Individual, o_2_alka) |> 
  # Add temperature and country information
  merge(temperature |> 
          mutate(Subspecies = country |> 
                   sbsp_land_match())
        , all = T
        , sort = F) |>   
  as_tibble() |>
  (function(df) {
    # Calculate correlation
    r2 <- df |> 
      spearman_cor_clim_chc()
    
    # Dispersion plot to report correlation
    df |> 
      ggplot(aes(x = mean_temp_C
                 , y = o_2_alka
      )) +
      geom_point(shape = 21 
                 , aes( fill = Subspecies)
                 , color = "grey40"
                 # , alpha = 0.85
                 , size = 2) +
      scale_fill_viridis_d(guide =
                             guide_legend(label.theme =
                                            element_text(face = "italic"
                                                         , size = 8)
                                          , byrow = T
                                          # Force the shape in the legend to
                                          # include filling, so it can show
                                          # the color of the fill
                                          , override.aes = list(shape = 21
                                                                , size = 3.2))) +
      scale_x_continuous(breaks = seq(7, 19, by = 2), limits = c(7, 19)) +
      scale_y_continuous(breaks = seq(1, 3, by = 1), limits = c(0, 3.1)) +
      geom_text(r2
                , mapping = aes(x = 17
                                , y = 3.1
                                , label = paste("R:"
                                                , estimate
                                                , signif))
                , fill = NA
                , label.color = NA
                , size = (boxplot_theme$legend.text$size / .pt) * 1
                , fontface = "bold") +
      facet_wrap(~Task
                 , nrow = 1) +
      labs(y = "olefins / n-alkanes"
           , x = "Temperature (°C)") +
      nmds_theme
  })()

# Correlation between mean chain length and temperature  
cl_temp <- Prop_chain.length |> 
  select(Task:Subspecies, Individual, cl_w.mean) |> 
  # Add temperature and country information
  merge(temperature |> 
          mutate(Subspecies = country |> 
                   sbsp_land_match())
        , all = T
        , sort = F) |>   
  as_tibble() |> 
  (function(df) {
    # Calculate correlation
    r2 <- df |> 
      spearman_cor_clim_chc()
    
    # Dispersion plot to report correlation
    df |> 
      ggplot(aes(x = mean_temp_C
                 , y = cl_w.mean
      )) +
      geom_point(shape = 21 
                 , aes( fill = Subspecies)
                 , color = "grey40"
                 , size = 2) +
      scale_fill_viridis_d(guide =
                             guide_legend(label.theme =
                                            element_text(face = "italic"
                                                         , size = 8)
                                          , byrow = T
                                          # Force the shape in the legend to
                                          # include filling, so it can show
                                          # the color of the fill
                                          , override.aes = list(shape = 21
                                                                , size = 3.2))) +
      scale_x_continuous(breaks = seq(7, 19, by = 2), limits = c(7, 19)) +
      scale_y_continuous(breaks = seq(24, 32, by = 2), limits = c(24, 32)) +
      geom_text(r2
                , mapping = aes(x = 17
                                , y = 32
                                , label = paste("R:"
                                                , estimate
                                                , signif))
                , fill = NA
                , label.color = NA
                , size = (boxplot_theme$legend.text$size / .pt) * 1
                , fontface = "bold") +
      facet_wrap(~Task
                 , nrow = 1) +
      labs(y = "Mean CHCs length"
           , x = "Temperature (°C)") +
      nmds_theme
  })()

# Merge plots for correlation with temperature
temp_plot <- (olefins_temp / cl_temp) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") 

# Export compound plot
ggsave(here("figs"
            ,"tempVchc.png")
       , plot = temp_plot
       , dpi = "print"
       , width = 6
       , height = 4.5
       , units = "in"
       , device = agg_png
       , scaling = 0.8)

# Correlation between olefins to n-alkanes ratio and precipitation
olefins_precip <- Prop_CompsClass |> 
  # Calculate olefins abundance by summing abundance of alkenes and alkadienes
  # Calculate olefins to n-alkanes ratio
  mutate(olefins = Alkenes + Alkadienes
         , o_2_alka = olefins / Alkanes) |> 
  select(Task:Individual, o_2_alka) |> 
  # Add precipitation and country information
  merge(precipitation |> 
          mutate(Subspecies = country |> 
                   sbsp_land_match())
        , all = T
        , sort = F) |>  
  as_tibble() |> 
  (function(df) {
    # Calculate correlation
    r2 <- df |> 
      spearman_cor_clim_chc()
    
    # Dispersion plot to report correlation
    df |> 
      ggplot(aes(x = precipitation_mm
                 , y = o_2_alka
      )) +
      geom_point(shape = 21 
                 , aes( fill = Subspecies)
                 , color = "grey40"
                 , size = 2) +
      scale_fill_viridis_d(guide =
                             guide_legend(label.theme =
                                            element_text(face = "italic"
                                                         , size = 8)
                                          , byrow = T
                                          # Force the shape in the legend to
                                          # include filling, so it can show
                                          # the color of the fill
                                          , override.aes = 
                                            list(shape = 21
                                                 , size = 3.2))) +
      scale_x_continuous(breaks = seq(40, 75, by = 5), limits = c(40, 75)) +
      scale_y_continuous(breaks = seq(1, 3, by = 1), limits = c(0, 3.1)) +
      geom_text(r2
                , mapping = aes(x = 70
                                , y = 3.1
                                , label = paste("R:"
                                                , estimate
                                                , signif))
                , fill = NA
                , label.color = NA
                , size = (boxplot_theme$legend.text$size / .pt) * 1
                , fontface = "bold") +
      facet_wrap(~Task
                 , nrow = 1) +
      labs(y = "olefins / n-alkanes"
           , x = "Precipitation (mm)") +
      nmds_theme
  })()

# Correlation between mean chain length and precipitation
cl_precip <- Prop_chain.length |> 
  select(Task:Subspecies, Individual, cl_w.mean) |>
  # Add precipitation and country information
  merge(precipitation |> 
          mutate(Subspecies = country |> 
                   sbsp_land_match())
        , all = T
        , sort = F) |>  
  as_tibble() |> 
  (function(df) {
    # Calculate correlation
    r2 <- df |> 
      spearman_cor_clim_chc()
    
    # Dispersion plot to report correlation
    df |> 
      ggplot(aes(x = precipitation_mm
                 , y = cl_w.mean
      )) +
      geom_point(shape = 21 
                 , aes( fill = Subspecies)
                 , color = "grey40"
                 , size = 2) +
      scale_fill_viridis_d(guide =
                             guide_legend(label.theme =
                                            element_text(face = "italic"
                                                         , size = 8)
                                          , byrow = T
                                          # Force the shape in the legend to
                                          # include filling, so it can show
                                          # the color of the fill
                                          , override.aes = 
                                            list(shape = 21
                                                 , size = 3.2))) +
      scale_x_continuous(breaks = seq(40, 75, by = 5), limits = c(40, 75)) +
      scale_y_continuous(breaks = seq(24, 32, by = 2), limits = c(24, 32)) +
      geom_text(r2
                , mapping = aes(x = 70
                                , y = 32
                                , label = paste("R:"
                                                , estimate
                                                , signif))
                , fill = NA
                , label.color = NA
                , size = (boxplot_theme$legend.text$size / .pt) * 1
                , fontface = "bold") +
      facet_wrap(~Task
                 , nrow = 1) +
      labs(y = "Mean CHCs length"
           , x = "Precipitation (mm)") +
      nmds_theme
  })()

# Merge plots for correlation with precipitation
precip_plot <- (olefins_precip / cl_precip) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

# Export compound plot
ggsave(here("figs"
            ,"precipVchc.png")
       , plot = precip_plot
       , dpi = "print"
       , width = 6
       , height = 4.5
       , units = "in"
       , device = agg_png
       , scaling = 0.8)

# Tables with mean and sd per compound of foragers and nurses of each subspecies ----

# Foragers table
foragers_table <- comps_summary(daten = datenc |> 
                         filter(grouping_info$Task == "Foragers")
                       , grouping.info = grouping_info  |> 
                         filter(Task == "Foragers")
                       , Comps.data = Comps_data
                       , grouping.factor = c("Subspecies"
                                             , "Task"))
colnames(foragers_table) <- colnames(foragers_table) |> 
  str_remove_all("_Foragers")

# Nurses table
nurses_table <- comps_summary(daten = datenc |> 
                         filter(grouping_info$Task == "Nurses")
                       , grouping.info = grouping_info  |> 
                         filter(Task == "Nurses")
                       , Comps.data = Comps_data
                       , grouping.factor = c("Subspecies"
                                             , "Task"))
colnames(nurses_table) <- colnames(nurses_table) |> 
  str_remove_all("_Nurses")


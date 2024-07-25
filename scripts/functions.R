# Script to create several custom functions


# PERMANOVA results ----
# Function creating a printable data frame to report the results of a
# PERMANOVA performed with "adonis()" or "adonis2()" from vegan.
permanova_table <- function(PMV){
  require(dplyr)
  require(vegan)
  sum_PMV <- PMV %>% permustats() %>% summary()
  SES <- sum_PMV$z %>% as.data.frame() %>% 
    mutate(faktoren = row.names(sum_PMV$z %>% as.data.frame()))
  colnames(SES) <- c("SES", "faktoren")
  
  PMV <- PMV %>% as.data.frame() %>% mutate(faktoren = row.names(PMV))
  colnames(PMV) <- c("Df", "SS", "R2", "F", "p", "faktoren")
  
  PMV_table <- merge(PMV, SES, by = "faktoren", all.x = T) %>% 
  arrange(factor(faktoren, levels = row.names(PMV))) %>% 
  select(faktoren:F,'SES','p')
  rownames(PMV_table) <- PMV_table$faktoren
  PMV_table %>% select(!faktoren)
}

# permutest results ----
# Function creating a printable data frame to report the results of a
# PERMANOVA performed with "adonis()" or "adonis2()" from vegan.
permutest_table <- function(PMV, factor){
  require(dplyr)
  require(vegan)
  sum_PMV <- PMV %>% permustats() %>% summary()
  SES <- sum_PMV$z %>% as.data.frame() %>% 
    mutate(faktoren = factor)
  colnames(SES) <- c("SES", "faktoren")
  rownames(SES) <- factor
  
  rownames(PMV$tab) <- c(factor, "Residual")
  
  PMV <- PMV$tab %>% 
    mutate(faktoren = row.names(PMV$tab)) %>%
    select(-N.Perm)
  colnames(PMV) <- c("Df", "SS", "MS", "F", "p", "faktoren")
  
  PMV_table <- merge(PMV, SES, by = "faktoren", all.x = T) %>% 
    arrange(factor(faktoren, levels = row.names(PMV))) %>% 
    select(faktoren:F,'SES','p')
  rownames(PMV_table) <- PMV_table$faktoren
  PMV_table %>% select(!faktoren) %>%
    rbind(Total = c(sum(PMV_table$Df)
                    , sum(PMV_table$SS)
                    , sum(PMV_table$MS)
                    , rep(NA, 3)))
}

# Data frame by hydrocarbon classes ----
# Calculates the total relative abundance of each hydrocarbon classes
# for each individual.
cclasses_df <- function(daten, grouping.info, Comps.data, fuse.methyls){
  require(dplyr)
  if(exists("fuse.methyls") == F){
    fuse.methyls = F
  } else{
    fuse.methyls = fuse.methyls
  }
  
  # Data set per hydrocarbon class
  ## Vector to store compund classes` names
  cclasses <- c()
  
  for (c.class in unique(Comps.data$Class)) {
    cclasses <- c(cclasses, c.class) 
    tmp <- as.data.frame(t(daten)) %>% 
      filter(Comps.data$Class == c.class)
    assign(c.class, tmp)
  }
  
  Sampls_index <- data.frame('Sample name' = rownames(daten)
                             , rank = rank(grouping.info$Individual)
                             , stringsAsFactors  =  FALSE)
  
  if(fuse.methyls == T){
    Prop.CompsClass <- data.frame(grouping.info
                                  , Alkanes = colSums(Alkane)
                                  , Alkenes = colSums(Alkene)
                                  , Alkadienes = colSums(Alkadiene)
                                  , Methyls = rbind(Methyl
                                                    , Dimethyl) %>%
                                    colSums()
                                  , row.names  =  Sampls_index$Sample.name)
  } else{
    # Has to be adjusted to let the function work
    Prop.CompsClass <- data.frame(grouping.info
                                  , Alkanes = colSums(Alkane)
                                  , Alkenes = colSums(Alkene)
                                  , Alkadienes = colSums(Alkadiene)
                                  , Monomethyls = colSums(Methyl)
                                  , Dimethyls = colSums(Dimethyl)
                                  , row.names  =  Sampls_index$Sample.name)
  }
  
  Prop.CompsClass %>% as_tibble()
}

# Difference test for hydrocarbon classes abundance between group ----
# Perform test (U or  kw) for the difference in abundance of each hydrocabron class
# between groups of a specified factor (faktor)
cclasses_test <- function(Prop.CompsClass.long, faktor, test){
  require(dplyr)
  
  # Mann-Whitney U test
  if(test == "U"){
    df_list <- c()
    results <- list()
    lap = 1
    for (c.class in unique(Prop.CompsClass.long$CClass)) {
      c.class
      df <- Prop.CompsClass.long %>% filter(CClass == c.class)
      tmp <- wilcox.test(data = df
                             , abundance ~ get(faktor))
      
      results[[lap]] <- tmp
      df_list <- c(df_list, c.class)
      lap = lap+1
    }
  }
  
  # Kruskal-Wallis test
  if(test == "kw"){
    
    require(dunn.test)
    
    df_list <- c()
    results <- list()
    lap = 1
    for (c.class in unique(Prop.CompsClass.long$CClass)) {
      c.class
      df <- Prop.CompsClass.long %>% filter(CClass == c.class)
      tmp <- list(kruskal = kruskal.test(abundance ~ get(faktor), data = df)
                  , dunn = dunn.test(df$abundance
                                     , df %>% pull(faktor)
                                     , kw = FALSE
                                     , table = FALSE
                                     , list = TRUE
                                     , method = 'bh'
                                     , altp = TRUE
                                     ))
      # Post-hoc Dunn's test
      tmp$dunn <- cbind.data.frame(comparisons = tmp$dunn$comparisons
                                   , Z = tmp$dunn$Z
                                   , P = tmp$dunn$altP
                                   , P.adjust = tmp$dunn$altP.adjust
                                   , significant = (tmp$dunn$altP.adjust<=
                                                      0.05)) %>%
        as_tibble
      
      results[[lap]] <- tmp
      df_list <- c(df_list, c.class)
      lap = lap+1
    }
  }
  names(results) <- df_list
  results
}

# Table of results of difference test for hydrocarbon classes among groups ----
# Build printable table to report results of statistical test for the difference
# in abundance of hydrocarbon classes between groups.
cclasses_test_table <- function(test_df_list, test_type){
  names_test_dfs <- test_df_list |> 
    names()
  results_table <- data.frame(statistic = c()
                              , df = c()
                              , p = c())
  for (name in names_test_dfs) {
    # For Kruskal-Wallis test results
    if(test_type == "KW") {
      name_table <- data.frame(Chi2 = test_df_list[[name]][["kruskal"]][["statistic"]] |> 
                                    round(3)
                                  , Df = test_df_list[[name]][["kruskal"]][["parameter"]][["df"]] |> 
                                    as.integer()
                                  , p = test_df_list[[name]][["kruskal"]][["p.value"]] #|> 
                                    # format(digits = 4
                                    #        , nsmall = 3)
                                  )
    }
    
    # For Mann-Whitney U test results
    if(test_type == "U") {
      name_table <- data.frame(W = test_df_list[[name]][["statistic"]][["W"]] |> 
                                 round(3)
                               , p = test_df_list[[name]][["p.value"]] #|> 
                                 # format(digits = 4
                                 #        , nsmall = 3)
      )
    }
    
    # Format decimals for p values
    if(name_table$p < 0.001){
      name_table$p <- "<0.001"
    } else {
      name_table$p <- name_table$p |> 
        round(3)
    }
    
    results_table <- rbind(results_table, name_table)
  }
  names_test_dfs[names_test_dfs == "Methyls"] <- "Methyl alkanes"
  # row.names(results_table) <- names_test_dfs
  results_table <- cbind.data.frame("class" = names_test_dfs
                                    , results_table)
  row.names(results_table) <- NULL
  results_table
}

# Chain length data set ----
## Function to generate a data set with the weighted average chain length per
## individual
cl_df <- function(daten, Comps.data, grouping.info){
  cl.df <- data.frame(row.names = rownames(daten))
  Ncols <- as.character()
  for (i in unique(Comps.data$Chain.length)) {
    #We need to sum the abundance of the compounds#
    Tmp <- as.data.frame(t(daten)) %>%
      filter(Comps.data$Chain.length == i) %>%
      colSums() %>% 
      as.numeric()
    cl.df <- cbind.data.frame(cl.df, Tmp)
    Ncols <- c(Ncols, i)
  }
  rm(i, Tmp)
  colnames(cl.df) <- Ncols
  rm(Ncols)
  
  cl.df <- cbind(grouping.info, cl.df)
  
  merge(grouping.info
        , cl.df %>% 
          pivot_longer(cols = !where(is.factor)
                       , names_to = "Chain.length"
                       , values_to = "Abundance") %>%  
          mutate(Chain.length = as.numeric(Chain.length)) %>%
          group_by(Individual) %>%
          summarise(cl_w.mean = weighted.mean(Chain.length
                                              , Abundance))
        , by = "Individual") %>% as_tibble()
}

# Difference test for weighted mean chain length between groups ----
# Difference test for weighted mean chain length between groups of a specified
# factor (faktor).
cl_test <- function(Prop.chain.length, faktor, test.type){
  # Kruskal-Wallis test
  if (test.type == "kw"){
    require(dplyr)
    require(dunn.test)
    
    test <- list (kw = kruskal.test(data = Prop.chain.length
                                    , cl_w.mean ~ get(faktor))
                  , dunn = dunn.test(x = Prop.chain.length$cl_w.mean
                                     , g = Prop.chain.length %>% pull(faktor)
                                     , kw = FALSE
                                     , table = FALSE
                                     , list = TRUE
                                     , method = 'bh'
                                     , altp = TRUE
                                     ))
    # Post-hoc Dunn's test
    test$dunn <- cbind.data.frame(
      comparisons = test$dunn$comparisons
                                  , Z = test$dunn$Z
                                  , P = #test$dunn$P
                                    test$dunn$altP
                                  , P.adjust = #test$dunn$P.adjusted
                                    test$dunn$altP.adjust
                                  , significant = (test$dunn$altP.adjust<=
                                                     #test$dunn$P.adjusted<=
                                                     #(0.05/2))) %>%
                                                     0.05)) %>%
      as_tibble
  }
  # Mann-Whitney U test
  if (test.type == "U"){
    test <- wilcox.test(data = Prop.chain.length
                        , cl_w.mean ~ get(faktor))
  }
  test
}

# Table of results of groups' hydrocarbon classes abundance difference test ----
# Build printable table to report results of statistical test for the difference
# in weighted mean chain length between groups.
cl_test_table <- function(test_list, test_type){
  if(test_type == "KW") {
    test_table <- data.frame(Chi2 = test_list[["kw"]][["statistic"]][["Kruskal-Wallis chi-squared"]] |> 
                               round(3)
                             , Df = test_list[["kw"]][["parameter"]][["df"]] |> 
                               as.integer()
                             , p = test_list[["kw"]][["p.value"]]
                             )
  }
    
  if(test_type == "U") {
    test_table <- data.frame(W = test_list[["statistic"]][["W"]] |> 
                               round(3)
                             , p = test_list[["p.value"]]
                             )
  }
    
  if(test_table$p < 0.001){
    test_table$p <- test_table$p |> 
      format(digits = 4
             #, nsmall = 3
             , scientific = T)
  } else {
    test_table$p <- test_table$p |> 
      round(3)
    }
  test_table
}

# Summary compounds table ----
## Generates a table with the RI, mean, and sd of each compound for each
## level of the selected factor
comps_summary <- function(daten, grouping.info, Comps.data, grouping.factor){
  require(tidyr)
  require(dplyr)
  
  mean_daten <- cbind.data.frame(grouping.info %>% 
                                   select(all_of(grouping.factor))
                                 , daten) %>%
    group_by(across(where(is.factor))) %>% 
    summarise(across(where(is.numeric)
                     , mean)) %>% 
    pivot_longer(cols = !where(is.factor)
                 , names_to = "Compound"
                 , values_to = "mean")
  mean_daten <- mean_daten %>% mutate(mean = round(mean, digits = 3))
  
  sd_daten <- cbind.data.frame(grouping.info %>% 
                                 select(all_of(grouping.factor))
                               , daten) %>%
    group_by(across(where(is.factor))) %>% 
    summarise(across(where(is.numeric)
                     , sd)) %>% 
    pivot_longer(cols = !where(is.factor)
                 , names_to = "Compound"
                 , values_to = "sd")
  sd_daten <- sd_daten %>% mutate(sd = round(sd, digits = 3))
  
  summarised_daten <- merge(mean_daten
                            , sd_daten
                            , by = c(grouping.factor, "Compound")
                            , sort = F) %>% 
    pivot_wider(names_from = grouping.factor
                , values_from = c(mean, sd)
                , names_glue = paste0(paste0("{"
                                             , grouping.factor
                                             , "}"
                                             , collapse = "_")
                                      , "_{.value}")
                , names_sort = F)
  
  summarised_daten <-summarised_daten %>%
    select(all_of(colnames(summarised_daten) %>%
                    sort()))
  
  summarised_daten <- cbind.data.frame(Compound = summarised_daten %>% 
                                         pull(Compound)
                                       , RI = Comps.data %>% 
                                         pull(RI)
                                       , summarised_daten %>% 
                                         select(-Compound)) |> 
    (function(df) {
      df_comps <- df |> 
        select(Compound:RI)
      
      df_daten <- df |> 
        column_to_rownames("Compound") |> 
        select(-RI) 
      
      df_daten <- df_daten |> 
        filter(rowSums(df_daten) > 0) |> 
        rownames_to_column("Compound")
      
      df_comps <- df_comps |> 
        filter(Compound %in% df_daten$Compound)
      
      df_comps |> 
        merge(df_daten
              , all = T
              , sort = F)
    })()
}

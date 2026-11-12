###########################################################################################
#
# Script dedicated to the study of RNA editing
# Adapted from a script by Matthew Birk
#
###########################################################################################
PATH <- here::here("data", "analyseDiff", "Other", "Soff_head_ES.csv")
# Read Editing Sites tables
readES <- function(PATH){
  a <- read.csv(PATH, header=TRUE)
  colnames(a)[c(2, 4)] = c('transcript', 'ES')
  
  a <- a %>%
    select(-X) %>%
    mutate(es_code = paste(transcript, ES, sep = ':'))
  
  # Calculate the standard deviation of each editing site-sample combination's editing level estimate -------------------------------------------------------
  
  samples = gsub('EL_', '', grep('EL_', colnames(a), value = TRUE))
  
  SDS = as.data.frame(t(apply(a, 1, function(i){
    EL = as.numeric(i[grep('EL_', names(i))])
    cov = as.numeric(i[grep('cov_', names(i))])
    SD = sqrt(EL * (1 - EL) / cov) # based on a binomial distribution of sampling editing events from reads
    names(SD) = paste0('EL_SD_', samples)
    return(SD)
  })))
  
  a = cbind(a, SDS)
  
  return(a)
}

# Tidy the sample table file
tidySample <- function(file){
  
a <-  file %>%
  dplyr::select(!File) %>%
    dplyr::rename(Indiv = SampleName) %>%
    mutate(Sample = paste(Treatment, Indiv, sep="_"))
  
  return(a)
}

# Make correspondance between samples and treatments
SampTreat <- function(file){
  
  ctrl_samples = paste('EL', na.omit(subset(file, Treatment == 'CT8.1')$sample_ID), sep = '_')
  OA_samples = paste('EL', na.omit(subset(file, Treatment == 'CT7.7')$sample_ID), sep = '_')
  Hg_samples = paste('EL', na.omit(subset(file, Treatment == 'Hg8.1')$sample_ID), sep = '_')
  HgOA_samples = paste('EL', na.omit(subset(file, Treatment == 'Hg7.7')$sample_ID), sep = '_')
  
  samples_treatment = rbind(data.frame(treatment = 'ctrl', sample = ctrl_samples),
                            data.frame(treatment = 'OA', sample = OA_samples),
                            data.frame(treatment = 'Hg', sample = Hg_samples),
                            data.frame(treatment = 'HgOA', sample = HgOA_samples))
  
  return(samples_treatment)
}

# Retrieving the known ES from Shoshan2021
retrieveKnownES <- function(PATH){
  ES = read.table(PATH, sep = '\t', header = TRUE, quote = '')
  ES$es_code = paste(ES$Trinity_id, ES$Editing_location, sep = ':')
  
  return(ES)
  
}

# Display coverage histogram
displayCovHist <- function(file){
  
  a <- hist(unlist(as.vector(file[, grep('cov_', colnames(file))])), xlim = c(0, 1000), breaks = seq(0, 1e6, by = 10), xlab = 'coverage') # visualize how many of all ESs across all the samples had how much coverage.
  
  return(a)
}


# Remove low coverage sites
rmLowCov <- function(file, SD, COV){
  
  min_EL_SD_ratio = SD
  min_coverage = COV
  
  tmp2 = apply(file, 1, function(i){ # remove those ES where either EL < 2x its own SD or coverage < 20 reads
    tmp_el = as.numeric(i[grep('EL_SEP2', names(i), value = TRUE)])
    tmp_sd = as.numeric(i[grep('EL_SD', names(i), value = TRUE)])
    tmp_cov = as.numeric(i[grep('cov_', names(i), value = TRUE)])
    
    i[grep('EL_SEP2', names(i), value = TRUE)[which(tmp_el < (tmp_sd * min_EL_SD_ratio))]] = 'fail'
    i[grep('EL_SEP2', names(i), value = TRUE)[which(tmp_cov < min_coverage)]] = 'fail'
    
    return(i)
  })
  
  d_cov = as.data.frame(t(tmp2))
  

  tmp3 = which(!(colnames(d_cov) %in% c('transcript', 'SwissProt_name', 'original_AA', 'target_AA', 'es_code')))
  d_cov[tmp3] = sapply(d_cov[tmp3], as.numeric) # after the numeric columns numeric again
  
  # ----------------------------------------------------------------------------------------------------------
  # THOMAS' ADDITION
  # To properly compare between conditions, we need each site to have enough reliable values in each condition
  # Hence, I propose to filter as follows:
  
  CT_cols <- grep("EL_SEP2.CT8.1", colnames(d_cov), value = TRUE)
  Hg_cols <- grep("EL_SEP2.Hg8.1", colnames(d_cov), value = TRUE)
  OA_cols <- grep("EL_SEP2.CT7.7", colnames(d_cov), value = TRUE)
  HgOA_cols <- grep("EL_SEP2.Hg7.7", colnames(d_cov), value=TRUE)
  
  # Here, adding more filtering to avoid any unconsistent comparison (e.g. NA vs value)
  keep_3in3 <- rowSums(cbind(
    rowSums(!is.na(d_cov[, CT_cols]))   >= 3,
    rowSums(!is.na(d_cov[, Hg_cols]))   >= 3,
    rowSums(!is.na(d_cov[, OA_cols]))   >= 3,
    rowSums(!is.na(d_cov[, HgOA_cols])) >= 3)) >= 3
  
  cat("Option ≥3 dans 3/4   :", sum(keep_3in3), "\n") # 83,503 sites -> reduces a lot but ensure to have only consistent comparisons
  
  data_final <- d_cov[keep_3in3, ] 
  data_final$Type <- ifelse(data_final$original_AA == data_final$target_AA, "Synonymous", "Non-synonymous")
  
  return(data_final)
}

# Build data as matrix
buildMat <- function(file, ST){
  
  EL_cols <- grep("^EL_", names(file), value=TRUE)
  Ref_col <- grep("Editing_Level", names(file), value=TRUE)
  
  EL_mat <- file %>%
    dplyr::select(all_of(Ref_col), all_of(EL_cols)) %>%
    as.matrix()
  
  rownames(EL_mat) <- file$es_code
  colnames(EL_mat) <- gsub("^EL_SEP2.", "", colnames(EL_mat))
  
  EL_mat.long1 <- as.data.frame(EL_mat) %>%
    tibble::rownames_to_column(var="site") %>%
    tidyr::pivot_longer(cols = starts_with(c("CT", "Hg", "Ed")),
                        names_to = "Sample",
                        values_to = "Freq",
                        values_drop_na = FALSE) %>%
    dplyr::select(!starts_with("EL_")) %>%
    tidyr::pivot_wider(names_from = site, values_from = Freq) %>%
    filter(Sample != "Editing_Level")
  
  EL_mat.test <- merge(EL_mat.long1, ST[,2:5], all.x=TRUE)
  
  # Final dataset to be used for statistical analysis
  Final <- EL_mat.test %>%
    dplyr::select(1, 83505:83507, 2:83504)
  
  Final$Treatment <- relevel(factor(Final$Treatment), ref="CT8.1")
  
  return(Final)
}

#file <- tar_read(final.mat)
#file2 <- tar_read(dCov)
#FDR1 <- 0.1
#FDR2 <- 0.05

rundEL <- function(file, file2, FDR1, FDR2){
  
  meta <- file[, c("Sample", "Treatment", "Mercury", "pCO2")] # Metadata table
  
  
  # Fit model + extract pairwise Wald contrasts for every site
  
  run_omnibus <- function(variable, meta) {
    df <- data.frame(myX = variable, meta)
    df <- df[is.finite(df$myX), ]
    if (nrow(df) < 3) return(list(p.value = NA, model = NULL))
    
    model <- glm(myX ~ Treatment, data = df, na.action = na.omit, family=quasibinomial())
    
    # Omnibus Wald test for the Treatment factor (equivalent to LRT)
    p_omnibus <- drop1(model, test = "F")["Treatment", "Pr(>F)"]
    
    return(list(p.value = p_omnibus, model = model))
  }
  
  # Run model
  raw_results <- apply(file[, 5:ncol(file)], MARGIN = 2, FUN = run_omnibus, meta = meta, simplify = FALSE)
  
  # Retrieve p.values and p.adjusted
  pvals <- sapply(raw_results, `[[`, "p.value")
  pvals_fdr <- p.adjust(pvals, method = "fdr")
  
  # Filter significant sites: FDR < 0.1
  sig_sites <- names(pvals_fdr[!is.na(pvals_fdr) & 
                                 pvals_fdr < FDR1])
  
  
  posthoc <- do.call(rbind, lapply(sig_sites, function(site) {
    coefs        <- as.data.frame(summary(raw_results[[site]]$model)$coefficients)
    colnames(coefs)[4] <- "p.value"
    coefs$contrast     <- rownames(coefs)
    coefs$site         <- site
    return(coefs[-1, ])  # drop intercept
  }))
  
  # Computing p.adjusted
  posthoc$p.fdr <- p.adjust(posthoc$p.value, method = "fdr")
  
  # Retrieving final significant sites (FDR < 0.05)
  posthoc.sig <- posthoc %>%
    dplyr::filter(p.fdr < FDR2)

  # Reference level for treatment:
  ref_level <- levels(factor(meta$Treatment))[1]

  # Mean editing par site x treatment for significant sites only
  group_means <- do.call(rbind, lapply(sig_sites, function(site) {
    df <- data.frame(myX = file[[site]], meta)
    df <- df[is.finite(df$myX), ]
    agg <- aggregate(myX ~ Treatment, data = df, FUN = mean)
    data.frame(site = site, Treatment = agg$Treatment, mean_editing = agg$myX)
  }))
  
  posthoc.sig <- posthoc.sig %>%
  dplyr::mutate(treatment = sub("^Treatment", "", contrast)) %>%
  dplyr::left_join(group_means, by = c("site", "treatment" = "Treatment")) %>%
  dplyr::rename(mean_treatment = mean_editing) %>%
  dplyr::left_join(group_means %>% dplyr::filter(Treatment == ref_level) %>%
              dplyr::select(site, mean_ref = mean_editing),
            by = "site") %>%
  dplyr::mutate(mean_diff = mean_treatment - mean_ref,
         odds_ratio = exp(Estimate))
  
  # Now, build the final result table:
  Final.res <- posthoc.sig %>%
    dplyr::select(site, contrast, Estimate, odds_ratio, mean_ref, mean_treatment, mean_diff, p.fdr) %>%
    dplyr::left_join(
      file2[, c("es_code", "transcript", "SwissProt_name", "ES", "original_AA", "target_AA")],
      by = c("site" = "es_code")) %>%
    dplyr::mutate(trend = ifelse(Estimate > 0, "more", "less"),
           type = ifelse(original_AA == target_AA, "synonymous", "non-synonymous"),
           comparison = ifelse(contrast == "TreatmentCT7.7", "OA_vs_CT",
                               ifelse(contrast == "TreatmentHg8.1", "Hg_vs_CT",
                                      ifelse(contrast == "TreatmentHg7.7", "HgOA_vs_CT", NA))),
            diff_pct = mean_diff * 100) %>%
    dplyr::select(comparison, SwissProt_name, transcript, ES, Estimate, odds_ratio, mean_ref, mean_treatment, mean_diff, diff_pct, trend, original_AA, target_AA, type)
  
  polar <- c("R","N","D","C","Q","E","H","K","S","T","Y")
  charg <- c("K", "R", "H", "D", "E")
  Final.res$polarity <- ifelse(Final.res$original_AA %in% polar & Final.res$target_AA %in% polar,
                          "polar_to_polar",
                          ifelse(Final.res$original_AA %in% polar & !(Final.res$target_AA %in% polar),
                                 "polar_to_nonpolar",
                                 ifelse(!(Final.res$original_AA %in% polar) & Final.res$target_AA %in% polar,
                                        "nonpolar_to_polar",
                                        "nonpolar_to_nonpolar")))
  Final.res$charge <- ifelse(Final.res$original_AA %in% polar & !(Final.res$original_AA %in% charg) & Final.res$target_AA %in% polar & Final.res$target_AA %in% charg,
                          "polarNotCharged_to_charged",
                        ifelse(Final.res$original_AA %in% polar & !(Final.res$original_AA %in% charg) & Final.res$target_AA %in% polar & !(Final.res$target_AA %in% charg),
                          "polarNotCharged_to_notcharged",
                        ifelse(Final.res$original_AA %in% polar & !(Final.res$original_AA %in% charg) & !(Final.res$target_AA %in% polar),
                          "polarNotCharged_to_nonpolar",
                        ifelse(!(Final.res$original_AA %in% polar) & Final.res$target_AA %in% polar & !(Final.res$target_AA %in% charg),
                          "nonpolar_to_polarNotCharged",
                        ifelse(!(Final.res$original_AA %in% polar) & Final.res$target_AA %in% polar & Final.res$target_AA %in% charg,
                          "nonpolar_to_polarCharged",
                        ifelse(Final.res$original_AA %in% polar & Final.res$original_AA %in% charg & Final.res$target_AA %in% polar & Final.res$target_AA %in% charg,
                          "charged_to_charged",
                        ifelse(Final.res$original_AA %in% polar & Final.res$original_AA %in% charg & Final.res$target_AA %in% polar & !(Final.res$target_AA %in% charg),
                        "charged_to_polarNotCharged",
                        ifelse(!(Final.res$original_AA %in% polar) & !(Final.res$target_AA %in% polar),
                        "nonpolar_to_nonpolar",
                        "charged_to_nonpolar"))))))))
  
  return(Final.res)
}


#file <- tar_read(dEL)
displayTornado <- function(file){
  
  sub_counts <- file %>%
    dplyr::mutate(substitution = paste(original_AA, target_AA, sep = "→")) %>%
    #dplyr::filter(type != 'synonymous') %>%
    dplyr::group_by(trend) %>%
    dplyr::count(substitution, sort = TRUE) %>%
    dplyr::mutate(N = ifelse(trend == "more", n, -n))
  
  
  # Type tornado plot
  p <- ggplot(sub_counts, aes(x = reorder(substitution, N), y = N, fill=trend)) +
    geom_col(width=0.1) +
    geom_point(size=3) +
    scale_y_continuous(labels=abs) +
    scale_fill_manual(values = c("#046C9A", "#FD6467")) +
    coord_flip() +
    geom_hline(yintercept = 0, color = "grey") +
    labs(x = "Amino acid substitution",
         y = "Number of editing sites",
         fill="diff. editing trend") +
    theme_bw() +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=12))
  
  return(p)
}


displayCircos <- function(file, CONTRAST, PATH){
  
  aa <- c("A","C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T","V", "W", "Y")
  
  full_pairs <- expand.grid(
    original_AA = aa,
    target_AA   = aa,
    stringsAsFactors = FALSE
  )
  
  sub_counts2 <- file %>%
    dplyr::mutate(substitution = paste(original_AA, target_AA, sep = "→")) %>%
    dplyr::filter(#type != 'synonymous',
           comparison == CONTRAST) %>%
    dplyr::group_by(trend) %>%
    dplyr::count(substitution, sort = TRUE) %>%
    dplyr::mutate(N = ifelse(trend == "more", n, -n))
  
  # Matrix
  mat <- sub_counts2 %>%
    tidyr::separate(substitution, into = c("original_AA", "target_AA")) %>%
    dplyr::group_by(original_AA, target_AA) %>%
    dplyr::summarise(n = sum(abs(N)), .groups = "drop") %>%
    dplyr::right_join(full_pairs, by = c("original_AA", "target_AA")) %>%
    dplyr::mutate(n = tidyr::replace_na(n, 0)) %>%
    tidyr::pivot_wider(names_from = target_AA, values_from = n) %>%
    dplyr::select(original_AA, A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y) %>%
    as.data.frame()
  
  rownames(mat) <- mat$original_AA
  mat <- mat[,-1]
  mat <- as.matrix(mat)
  
  # colors
  AA_colors <- structure(c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666",
                   "#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00",
                   "#CAB2D6","#6A3D9A","#FFFF99","#B15928"),
                   names=colnames(mat))
  
  grid_colors <- c(AA_colors, AA_colors)
  
  link.col = data.frame(c("Q","K","E","I","M","R","S","E","Q","K","K","K","T","K","N","I","E","K","Y","L","K","P","K","K","A","S","E","K","L","E","Q","R","Q"),
                    c("R","E","G","V","V","G","G","G","R","R","R","R","A","R","D","M","G","E","C","L","K","P","K","K","A","S","E","K","L","E","Q","R","Q"), 
                    c("#E31A1C","#A6CEE3","#E7298A","#666666","#B2DF8A","#FDBF6F","#FF7F00","#E7298A","#E31A1C","#A6CEE3","#A6CEE3","#A6CEE3","#CAB2D6","#A6CEE3","#33A02C","#666666","#E7298A","#A6CEE3","#B15928",
                      "#1F78B433","#A6CEE333","#FB9A9933","#A6CEE333","#A6CEE333","#1B9E7733","#FF7F0033","#E7298A33","#A6CEE333","#1F78B433","#E7298A33","#E31A1C33","#FDBF6F33","#E31A1C33"))
  
  aa <- as.vector(colnames(mat))
  
  
  circos.clear()
  circos.par(gap.degree = 7)
  
  # plot
  png(PATH, width=3, height=3, res = 300, units = "in")
  chordDiagram(mat,
               diffHeight = -mm_h(2),
               directional = 1,
               direction.type = "arrows",
               link.arr.type = "big.arrow",
               link.arr.length = 0.05,
               order = aa,
               col = link.col,
               grid.col = grid_colors,
               keep.diagonal=FALSE,
               annotationTrack = c("grid", "name"))
  dev.off()
  
  
}

#file1 <- tar_read(d_all)
#file2 <- tar_read(d_sampleOK)
# Calculating the Editing Index
editingIndex <- function(file1, file2){
  
  samples = gsub('EL_', '', grep('EL_SEP2', colnames(file1), value = TRUE))
  samples = grep('EL_SEP2', colnames(file1), value = TRUE)
  
  ctrl_samples = paste('EL_SEP2', file2[startsWith(file2$Sample, "CT8.1"),]$Sample, sep = '_')
  OA_samples = paste('EL_SEP2', file2[startsWith(file2$Sample, "CT7.7"),]$Sample, sep = '_')
  Hg_samples = paste('EL_SEP2', file2[startsWith(file2$Sample, "Hg8.1"),]$Sample, sep = '_')
  HgOA_samples = paste('EL_SEP2', file2[startsWith(file2$Sample, "Hg7.7"),]$Sample, sep = '_')
  
  EIs = base::colSums(file1[, paste0('EL_', samples)] * file1[, paste0('cov_', samples)], na.rm = TRUE) / colSums(file1[, paste0('cov_', samples)], na.rm = TRUE) # calculate the editing index for each sample
  
  EIs = data.frame(sample = samples, EI = EIs, row.names = NULL)
  
  EIs[which(EIs$sample %in% gsub('EL_SEP2_', 'SEP2.', ctrl_samples)), 'treatment'] = 'Control'
  EIs[which(EIs$sample %in% gsub('EL_SEP2_', 'SEP2.', OA_samples)), 'treatment'] = 'OA'
  EIs[which(EIs$sample %in% gsub('EL_SEP2_', 'SEP2.', Hg_samples)), 'treatment'] = 'Hg'
  EIs[which(EIs$sample %in% gsub('EL_SEP2_', 'SEP2.', HgOA_samples)), 'treatment'] = 'Hg + OA'
  
  EIs$treatment = ordered(EIs$treatment, levels = c('Control', 'OA', 'Hg', 'Hg + OA'))
  
  return(EIs)
}

# Testing difference between treatment with Anova:
anovaEI <- function(file){
  
  anova(lm(EI ~ treatment, data = file))
  
}

# Plotting EI
plotEI <- function(file){
  
  ggplot(file, aes(treatment, EI, fill = treatment)) +
    geom_boxplot() +
    geom_point() +
    scale_fill_manual(values=c("#50aa3a", "#ece87c","#30bcd7", '#e97132')) + # modify to match the color scheme used in the paper.
    scale_x_discrete(labels=c("Ctrl", expression(CO[2]), expression(Hg[Low]), expression(CO[2]+Hg[Low]))) +
    theme_bw() +
    theme(legend.position = 'none') +
    labs(x = '', y = 'Editing index')
  
}

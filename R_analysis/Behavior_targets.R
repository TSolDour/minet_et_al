# ==============================================================================
# Minet et al. - Unmasked in a changing ocean: Acidification and mercury impair cuttlefish defense mechanisms
# Main script of the targets pipeline managing Behavioral analyses.
# Created by Antoine Minet and Thomas Sol Dourdin
# All figures used for the publication are exported with the prefix "Figure"
# ==============================================================================



rm(list = ls())
# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(visNetwork)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c("dplyr", "tidyr", "purrr", "readxl", "broom", "tibble", "ggplot2", "ggpattern", "patchwork", "wesanderson",
               "permuco", "RVAideMemoire", "rcompanion", "FactoMineR", "missMDA", "factoextra"))

# Run the R scripts in the R/ folder with your custom functions:
tar_source(here::here("R", "Behavior"))


list(
  # ===========================================================================
  # 1. IMPORT
  # ===========================================================================
  tar_target(file1, here::here("data", "behavior", "Minet_et_al.xlsx"), format = "file"),
  tar_target(data_2020, read_data_2020(file1)),
  tar_target(data_2021, read_data_2021(file1)),
  
  # ===========================================================================
  # 2. ANALYSES : BRIGHTNESS (comparisons to LI=1 + ANOVA + perm. ANOVA)
  # ===========================================================================
  tar_target(res_bright_2020_d15, analyse_variable(data_2020, "Brightness_score15", mu = 1, post_hoc = FALSE)), 
  tar_target(res_bright_2020_d30, analyse_variable(data_2020, "Brightness_score30", mu = 1, post_hoc = FALSE)), 
  tar_target(res_bright_2021_d35, analyse_variable(data_2021, "Brightness_score", mu = 1, post_hoc = FALSE)), 
  
  # 2-way ANOVA (Condition x Day), data 2020
  tar_target(data_2020_long_bright, pivot_two_way(data_2020, cols = c("Brightness_score15","Brightness_score30"), values_to = "Brightness_score")),
  tar_target(res_bright_twoway, twoway_aov(data_2020_long_bright, "Brightness_score", id = "Ind")), #Condition : F=0.868 et p=0.458 ; days F=0.019  p=0.889 ; interaction F=0.323 ;  p=0.809)
  
  
  # ===========================================================================
  # 3. FIGURES : BRIGHTNESS
  # ===========================================================================
  tar_target(panel_bright_d15, make_brightness_panel(data_2020, "Brightness_score15", year = "2020", day = "Day 15", y_label = "Brightness score")),
  tar_target(panel_bright_d30, make_brightness_panel(data_2020, "Brightness_score30", year = "2020", day = "Day 30")),
  tar_target(panel_bright_d35, make_brightness_panel(data_2021, "Brightness_score", year = "2021", day = "Day 35", use_pattern = TRUE)),
  tar_target(fig_brightness, combine_panels(panel_bright_d15, panel_bright_d30, panel_bright_d35)),
  tar_target(fig_brightness_file, save_figure(fig_brightness, here::here("results", "Behavior", "figures", "Figure2a_Brightness.png"))),
  
  # ===========================================================================
  # 4. ANALYSES : DISRUPTIVE (ANOVA + perm. anova + post-hoc)
  # ===========================================================================
  tar_target(res_disr_2020_d15, analyse_variable(data_2020, "Disr15")),
  tar_target(res_disr_2020_d30, analyse_variable(data_2020, "Disr30")),
  tar_target(res_disr_2021_d35, analyse_variable(data_2021, "Disr35")),
  
  # 2-way ANOVA (Condition x Day), data 2020
  tar_target(data_2020_long_disr, pivot_two_way(data_2020, cols = c("Disr15", "Disr30"), values_to = "Disr")),
  tar_target(res_disr_twoway, twoway_perm(data_2020_long_disr, "Disr", id = "Ind")),
  
  # ===========================================================================
  # 5. FIGURES : DISRUPTIVE
  # ===========================================================================
  tar_target(panel_disr_d15, make_disruptive_panel(data_2020, "Disr15", labels = labels_disruptive$d15, year = "2020", day = "Day 15", y_label = "Disruptive score")),
  tar_target(panel_disr_d30, make_disruptive_panel(data_2020, "Disr30", labels = labels_disruptive$d30, year = "2020", day = "Day 30")),
  tar_target(panel_disr_d35, make_disruptive_panel(data_2021, "Disr35", labels = labels_disruptive$d35, year = "2021", day = "Day 35", use_pattern = TRUE)),
  tar_target(fig_disruptive, combine_panels(panel_disr_d15, panel_disr_d30, panel_disr_d35)),
  tar_target(fig_disruptive_file, save_figure(fig_disruptive, here::here("results", "Behavior", "figures", "Figure2b_disruptive.png"))),
  
  # ===========================================================================
  # 6. LATERALITY INDEX (LI) : analysis (comparaisison to 0)
  #    post_hoc = TRUE only when Condition effect is significant
  # ===========================================================================
  tar_target(res_li_2020_d5, analyse_variable(data_2020, "LI5",  mu = 0, post_hoc = FALSE)),
  tar_target(res_li_2020_d14, analyse_variable(data_2020, "LI14", mu = 0, post_hoc = FALSE)),
  tar_target(res_li_2020_d25, analyse_variable(data_2020, "LI25", mu = 0, post_hoc = TRUE)),
  tar_target(res_li_2021_d14, analyse_variable(data_2021, "LI14", mu = 0, post_hoc = FALSE)),
  tar_target(res_li_2021_d31, analyse_variable(data_2021, "LI31", mu = 0, post_hoc = TRUE)),
  
  # 2-way Anova + perm. Anova (Condition x Jour), 2020 : days 5, 14, 25
  tar_target(data_2020_long_li, pivot_two_way(data_2020, c("LI5", "LI14", "LI25"), values_to = "LI")),
  tar_target(res_li_twoway_perm_2020, twoway_perm(data_2020_long_li, "LI", id = "Ind")),
  tar_target(res_li_twoway_aov_2020, twoway_aov(data_2020_long_li, "LI", id = "Ind")),
  
  # 2-way Anova + perm. Anova, 2021 : jours 14, 31
  tar_target(data_2021_long_li, pivot_two_way(data_2021, c("LI14", "LI31"), values_to = "LI")),
  tar_target(res_li_twoway_perm_2021, twoway_perm(data_2021_long_li, "LI", id = "Individu")),
  tar_target(res_li_twoway_aov_2021, twoway_aov(data_2021_long_li, "LI", id = "Individu")),
  
  # ===========================================================================
  # 7. FIGURES : LATERALITY INDEX (5 panneaux, 3 colonnes)
  # ===========================================================================
  tar_target(panel_li_2020_d5, make_li_panel(data_2020, "LI5", year = "2020", day = "Day 5", y_label = "Laterality index")),
  tar_target(panel_li_2020_d14, make_li_panel(data_2020, "LI14", year = "2020", day = "Day 14")),
  tar_target(panel_li_2020_d25, make_li_panel(data_2020, "LI25", year = "2020", day = "Day 25", y_label = "Laterality index", stars = stars_li$d25)),
  tar_target(panel_li_2021_d14, make_li_panel(data_2021, "LI14", year = "2021", day = "Day 14", y_label = "Laterality index", use_pattern = TRUE)),
  tar_target(panel_li_2021_d31, make_li_panel(data_2021, "LI31", year = "2021", day = "Day 31", y_label = "Laterality index", use_pattern = TRUE, stars = stars_li$d31, group_letters = letters_laterality_d31)),
  tar_target(fig_li, combine_panels(panel_li_2020_d5, panel_li_2020_d14, panel_li_2020_d25, panel_li_2021_d14, panel_li_2021_d31, ncol = 3)),
  tar_target(fig_li_file, save_figure(fig_li, here::here("results", "Behavior", "figures", "Figure3b_laterality_index.png"), width = 12, height = 8)),
  
  
  # ===========================================================================
  # 8. ABSOLUTE LATERALITY INDEX (AbsLI) : analyses
  # ===========================================================================
  tar_target(res_absli_2020_d5, analyse_variable(data_2020, "AbsLI5",  post_hoc = FALSE)),
  tar_target(res_absli_2020_d14, analyse_variable(data_2020, "AbsLI14", post_hoc = FALSE)),
  tar_target(res_absli_2020_d25, analyse_variable(data_2020, "AbsLI25", post_hoc = FALSE)),
  tar_target(res_absli_2021_d14, analyse_variable(data_2021, "AbsLI14", post_hoc = FALSE)),
  tar_target(res_absli_2021_d31, analyse_variable(data_2021, "AbsLI31", post_hoc = TRUE)),
  
  # 2-way Anova + perm. Anova - 2020
  tar_target(data_2020_long_absli, pivot_two_way(data_2020, c("AbsLI5", "AbsLI14", "AbsLI25"), values_to = "AbsLI")),
  tar_target(res_absli_twoway_perm_2020, twoway_perm(data_2020_long_absli, "AbsLI", id = "Ind")),
  tar_target(res_absli_twoway_aov_2020, twoway_aov(data_2020_long_absli, "AbsLI", id = "Ind")),
  
  # Pairwise comparison between days 
  tar_target(res_absli_days_pairwise_2020, run_pairwise_perm(data_2020_long_absli, "AbsLI", group_var = "Days")),
  
  # 2-way Anova + perm. Anova - 2021
  tar_target(data_2021_long_absli, pivot_two_way(data_2021, c("AbsLI14", "AbsLI31"), values_to = "AbsLI")),
  tar_target(res_absli_twoway_perm_2021, twoway_perm(data_2021_long_absli, "AbsLI", id = "Individu")),
  tar_target(res_absli_twoway_aov_2021, twoway_aov(data_2021_long_absli, "AbsLI", id = "Individu")),
  
  # ===========================================================================
  # 9. FIGURES : ABSOLUTE LATERALITY INDEX
  # ===========================================================================
  tar_target(panel_absli_2020_d5, make_abs_panel(data_2020, "AbsLI5", year = "2020", day = "Day 5", y_label = "Absolute laterality index")),
  tar_target(panel_absli_2020_d14, make_abs_panel(data_2020, "AbsLI14", year = "2020", day = "Day 14")),
  tar_target(panel_absli_2020_d25, make_abs_panel(data_2020, "AbsLI25", year = "2020", day = "Day 25", y_label = "Absolute laterality index")),
  tar_target(panel_absli_2021_d14, make_abs_panel(data_2021, "AbsLI14", year = "2021", day = "Day 14", y_label = "Absolute laterality index", use_pattern = TRUE)),
  tar_target(panel_absli_2021_d31, make_abs_panel(data_2021, "AbsLI31", year = "2021", day = "Day 31", use_pattern = TRUE, group_letters = letters_laterality_d31)),
  tar_target(fig_absli, combine_panels(panel_absli_2020_d5, panel_absli_2020_d14, panel_absli_2020_d25, panel_absli_2021_d14, panel_absli_2021_d31, ncol = 3)),
  tar_target(fig_absli_file, save_figure(fig_absli, here::here("results", "Behavior", "figures", "Figure3a_absolute_laterality_index.png"), width = 12, height = 8)),
  
  # ===========================================================================
  # 10. VIGILANT SCANNING (2020)
  # ===========================================================================
  tar_target(res_vigi_2020_d5, analyse_variable(data_2020, "vigilance_mono5",  post_hoc = FALSE)),
  tar_target(res_vigi_2020_d14, analyse_variable(data_2020, "vigilance_mono14", post_hoc = FALSE)),
  tar_target(res_vigi_2020_d25, analyse_variable(data_2020, "vigilance_mono25", post_hoc = FALSE)),
  
  # 2-way Anova + perm. Anova
  tar_target(data_2020_long_vigi, pivot_two_way(data_2020, c("vigilance_mono5", "vigilance_mono14", "vigilance_mono25"), values_to = "vigilance_mono")),
  tar_target(res_vigi_twoway_perm, twoway_perm(data_2020_long_vigi, "vigilance_mono", id = "Ind")),
  
  # ===========================================================================
  # 11. FIGURES : VIGILANT SCANNING
  # ===========================================================================
  tar_target(panel_vigi_d5, make_abs_panel(data_2020, "vigilance_mono5", year = "2020", day = "Day 5", y_label = "Monocular vigilance")),
  tar_target(panel_vigi_d14, make_abs_panel(data_2020, "vigilance_mono14", year = "2020", day = "Day 14")),
  tar_target(panel_vigi_d25, make_abs_panel(data_2020, "vigilance_mono25", year = "2020", day = "Day 25")),
  tar_target(fig_vigi, combine_panels(panel_vigi_d5, panel_vigi_d14, panel_vigi_d25)),
  tar_target(fig_vigi_file, save_figure(fig_vigi, here::here("results", "Behavior", "figures", "monocular_vigilance.png"))),
  
  # ===========================================================================
  # 12. EYE USE PREFERENCE : tables et tests
  # ===========================================================================
  tar_target(table_cat_2020, make_categories_table(categories_counts_2020)),
  tar_target(table_cat_2021, make_categories_table(categories_counts_2021)),
  tar_target(res_cat_2020, test_categories(table_cat_2020)),
  tar_target(res_cat_2021, test_categories(table_cat_2021)),
  
  # ===========================================================================
  # 13. FIGURES : EYE USE PREFERENCE
  # ===========================================================================
  tar_target(panel_cat_2020, make_category_plot(categories_to_percent(table_cat_2020), year = "2020", star_x = stars_categories$y2020, y_label = "Eye use preference (%)")),
  tar_target(panel_cat_2021, make_category_plot(categories_to_percent(table_cat_2021), year = "2021", star_x = stars_categories$y2021, show_legend = TRUE)),
  tar_target(fig_categories, combine_panels(panel_cat_2020, panel_cat_2021)),
  tar_target(fig_categories_file, save_figure(fig_categories, here::here("results", "Behavior", "figures", "Figure3c_eye_use_categories.png"), width = 9, height = 4)),
  
  # ===========================================================================
  # 14. LOCOMOTION (2020)
  #     NB : les figures montrent Mobility_percent et thigmotaxis, mais les
  #     tests portent sur Mobility_time et nb_frame_ext (comme dans l'original).
  # ===========================================================================
  tar_target(res_loco_distance, analyse_variable(data_2020, "distance")),
  tar_target(res_loco_mobility, analyse_variable(data_2020, "Mobility_time")),
  tar_target(res_loco_thigmo, analyse_variable(data_2020, "nb_frame_ext", nperm = 50000)),
  
  # ===========================================================================
  # 15. FIGURES: LOCOMOTION
  # ===========================================================================
  tar_target(panel_loco_distance, make_locomotion_panel(data_2020, "distance", ylim = c(0, 1200), y_label = "Distance (cm)", letter_y = 1110)),
  tar_target(panel_loco_mobility, make_locomotion_panel(data_2020, "Mobility_percent", ylim = c(0, 32), y_label = "Time spent mobile (%)", letter_y = 30)),
  tar_target(panel_loco_thigmo, make_thigmotaxis_panel(data_2020)),
  tar_target(fig_locomotion, combine_panels(panel_loco_distance, panel_loco_mobility, panel_loco_thigmo)),
  tar_target(fig_locomotion_file, save_figure(fig_locomotion, here::here("results", "Behavior", "figures", "Figure1_Locomotion.png"))),
  
  # ===========================================================================
  # 16. PREDATION
  # ===========================================================================
  tar_target(table_preda_2020, make_predation_long(data_2020, id_col = "Ind", days = predation_days_2020)),
  tar_target(table_preda_2021, make_predation_long(data_2021, id_col = "Individu", days = predation_days_2021)),
  
  # Day comparison for each condition (Cochran + McNemar)
  tar_target(res_preda_days_2020, test_predation_by_condition(table_preda_2020)),
  tar_target(res_preda_days_2021, test_predation_by_condition(table_preda_2021, post_hoc = predation_posthoc_2021)),
  
  # Binomial (p = 0.5) per condition and day
  tar_target(res_preda_binom_2020, binom_by_day(table_preda_2020)),
  tar_target(res_preda_binom_2021, binom_by_day(table_preda_2021)),
  
  # Condition comparison for each day (chi-square)
  tar_target(res_preda_chisq_2020, chisq_by_day(table_preda_2020)),
  tar_target(res_preda_chisq_2021, chisq_by_day(table_preda_2021)),
  
  # Crabe (2020)
  tar_target(res_crab, test_crab(data_2020)),
  
  # ===========================================================================
  # 17. FIGURES : PREDATION
  # ===========================================================================
  tar_target(fig_preda_2020, make_predation_figure(table_preda_2020, predation_annotations$y2020)),
  tar_target(fig_preda_2020_file, save_figure(fig_preda_2020, here::here("results", "Behavior", "figures", "predation_2020.png"),  width = 9, height = 7)),
  tar_target(fig_preda_2021, make_predation_figure(table_preda_2021, predation_annotations$y2021)),
  tar_target(fig_preda_2021_file, save_figure(fig_preda_2021, here::here("results", "Behavior", "figures", "predation_2021.png"), width = 9, height = 7)),
  
  # ===========================================================================
  # 18. ACP 2020 : ACP, HCPC, clusters tests, figures
  # ===========================================================================
  tar_target(res_pca_2020, run_pca(data_2020, vars = pca_vars_2020)),
  tar_target(res_hcpc_2020, run_hcpc(res_pca_2020$pca)),
  tar_target(tab_clusters_2020, cluster_table(res_hcpc_2020, res_pca_2020$table)),
  tar_target(res_clusters_2020, test_cluster_table(tab_clusters_2020, within_condition = TRUE)),
  tar_target(fig_pca_2020_files, save_pca_figures(res_pca_2020, here::here("results", "Behavior", "figures", "pca_2020"))),
  tar_target(fig_clusters_2020, make_cluster_plot(tab_clusters_2020, annotations = cluster_annotations_2020,
                                                  condition_labels = cluster_condition_labels_2020, palette = "Zissou1")),
  tar_target(fig_clusters_2020_file, save_figure(fig_clusters_2020, here::here("results", "Behavior", "figures", "clusters_2020.png"), width = 6, height = 4)),
  
  # ===========================================================================
  # 19. ACP 2021 : ACP, HCPC, clusters tests, figures
  # ===========================================================================
  tar_target(res_pca_2021, run_pca(data_2021, vars = pca_vars_2021)),
  tar_target(res_hcpc_2021, run_hcpc(res_pca_2021$pca)),
  tar_target(tab_clusters_2021, cluster_table(res_hcpc_2021, res_pca_2021$table)),
  tar_target(res_clusters_2021, test_cluster_table(tab_clusters_2021)),
  tar_target(fig_pca_2021_files, save_pca_figures(res_pca_2021, here::here("results", "Behavior", "figures", "pca_2021"))),
  
  # ===========================================================================
  # 20. ACP 2020-2021 (sheet "acp2020_2021") : same + cluster figures
  # ===========================================================================
  tar_target(data_acp, read_data_acp(file1)),
  tar_target(res_pca_acp, run_pca(data_acp, vars = pca_vars_acp)),
  tar_target(res_hcpc_acp, run_hcpc(res_pca_acp$pca)),
  tar_target(tab_clusters_acp, cluster_table(res_hcpc_acp, res_pca_acp$table)),
  tar_target(res_clusters_acp, test_cluster_table(tab_clusters_acp, within_condition = TRUE)),
  tar_target(fig_pca_acp_files, save_pca_figures(res_pca_acp, here::here("results", "Behavior", "figures", "Figure4_pca_2020_2021"))),
  tar_target(fig_clusters_acp, make_cluster_plot(tab_clusters_acp, annotations = cluster_annotations_acp,
                                                 condition_labels = cluster_condition_labels_acp, palette = "Cavalcanti1")),
  tar_target(fig_clusters_acp_file, save_figure(fig_clusters_acp, here::here("results", "Behavior", "figures", "Figure4c_clusters_2020_2021.png"), width = 8, height = 4),)
)

# Sys.setenv(TAR_PROJECT="Behavior")

# tar_manifest(fields=command)

# tar_visnetwork(physics=TRUE)

# tar_make(ExportTan_BP_Plot)

# tar_read(Tan_BP_Plot)
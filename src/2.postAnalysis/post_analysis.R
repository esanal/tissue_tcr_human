# Main dir
setwd("../../")

# library imports
library(iNEXT)
library(ggplot2)
library(patchwork)
library(dplyr)
library(poweRlaw)
library(GGally)
library(ggprism)
library(GGally)
library(stringr)
library(future.apply)
library(reshape)
library(ComplexHeatmap)
library(vecsets)
library(foreach)
library(doParallel)
require(iNEXT)
require(tictoc)
library(ggExtra)
library(lsa)
library(purrr)
library(igraph)
library(tibble)
library(ggraph)
library(ggpubr)
library(ggsci)
library(scales)
library(stringdist)
library(RColorBrewer)

# custom functions
source("./src/2.postAnalysis/functions.R")
options(show.error.locations = TRUE)

# Clone file directories and where to find them is dictated here!
patterns <- c("*.1.clones.txt")
thresholds <- c("read_th_1")
results_dir <- c("./results/results_combinedGSIDs")

names(patterns) <- thresholds
clone_file_dirs <-
  lapply(patterns, function(pat) {
    list.files(
      path = c(results_dir),
        pattern = pat,
        full.names = TRUE)
  })

# read metadata
metadata <-
  read.csv("./data/metadata/metadata_gsid_cellNumber.csv")
colnames(metadata)[1:4] <- c("individual", "organ", "subset", "cell.count")
metadata$Cell <- metadata$subset

organ_groups <- list(
  non_lymphoid = c("Fat", "Skin", "Blood"),
  lymphoid = c("BM", "Blood"),
  full = c("Blood", "BM", "Skin", "Fat")
)


###############################
#####Analyze per Threshold#####
###############################
for (threshold in c("read_th_1")) {
  # read current threshold's data for all samples
  cat("\n--------------------------------------")
  cat("\n####Reading txt: ", threshold, "####")
  cat("\n--------------------------------------")
  
  clone_dfs <-
    lapply(
      clone_file_dirs[[threshold]],
      function(file) 
        data.frame(data.table::fread(file, sep = "\t"))
    )

  cat("\n1. Reading txt: Done.")
  
  # add IDs from file names
  sample.ids <-
    str_extract(clone_file_dirs[[threshold]], "(?<=/)[^/]+(?=\\.1)")
  
  names(clone_dfs) <- paste(sample.ids, sep = "_")
  # split clones into tra and trb
  cat("\n2. Splitting clones into TRA and TRB.")
  clone_dfs_tra_trb <- lapply(clone_dfs, split_clones)
  clone_dfs_tra <- lapply(clone_dfs_tra_trb, `[[`, 1)
  clone_dfs_trb <- lapply(clone_dfs_tra_trb, `[[`, 2)
  clone_abundance_tra <- lapply(clone_dfs_tra_trb, `[[`, 3)
  clone_abundance_trb <- lapply(clone_dfs_tra_trb, `[[`, 4)
  rm(clone_dfs, clone_dfs_tra_trb)
  cat("\nSplitting clones into TRA and TRB. Done.")
  
  
  
  #####################
  ####  Summarize  ####
  #####################
  cat("\n3. Summarizing clones of TRA and TRB.")
  # summarize each clone df of tra and trb
  summarylist_tra <- sapply(clone_dfs_tra, summary_clones)
  summarylist_trb <- sapply(clone_dfs_trb, summary_clones)
  
  summary_df_tra <- data.frame(t(summarylist_tra), chain = "TRA")
  summary_df_trb <- data.frame(t(summarylist_trb), chain = "TRB")
  
  summary_df_tra$run.id_gs.id <- rownames(summary_df_tra)
  summary_df_trb$run.id_gs.id <- rownames(summary_df_trb)
  
  # merge samples to plot
  ## TRA
  summary_df_tra <-
    data.frame(
      clone_count = unlist(summary_df_tra[["clone_count"]]),
      UMI_count = unlist(summary_df_tra[["UMI_count"]]),
      read_count = unlist(summary_df_tra[["read_count"]]),
      chain = "TRA",
      gs.id = unlist(str_split_i(summary_df_tra[["run.id_gs.id"]], "_", 1)),
      run.id = unlist(str_split_i(summary_df_tra[["run.id_gs.id"]], "_", 1))
    )
  ##TRB
  summary_df_trb <-
    data.frame(
      clone_count = unlist(summary_df_trb[["clone_count"]]),
      UMI_count = unlist(summary_df_trb[["UMI_count"]]),
      read_count = unlist(summary_df_trb[["read_count"]]),
      chain = "TRB",
      gs.id = unlist(str_split_i(summary_df_trb[["run.id_gs.id"]], "_", 1)),
      run.id = unlist(str_split_i(summary_df_trb[["run.id_gs.id"]], "_", 1))
    )
  

  # add metadata columns
  summary_df_tra <-
    left_join(
      summary_df_tra,
      metadata[c("individual", "organ", "subset", "cell.count", "gs.id", "Cell")],
      by = "gs.id")
  
  summary_df_trb <-
    left_join(
      summary_df_trb,
      metadata[c("individual", "organ", "subset", "cell.count", "gs.id", "Cell")],
      by = "gs.id"
    )

  ## add per_cell_metrics
  summary_df_tra <- per_cell_metrics(summary_df_tra)
  summary_df_trb <- per_cell_metrics(summary_df_trb)

  # Export summary metrics of samples
  ## create dir
  summary.path <- file.path("./data/summary", threshold)
  if (!dir.exists(summary.path)) {
    dir.create(summary.path, recursive = TRUE)
  }
  ## save
  write.csv(
    x = rbind(summary_df_tra, summary_df_trb),
    file = paste0(summary.path, "/sample_metrics.csv"),
    row.names = FALSE)

  cat("\n4. Plotting histograms for TRA and TRB.")
  # Plot histograms and scatter
  histograms_tra <-
    suppressWarnings(
      plot_histograms(summary_df_tra, folder = threshold, file_ext = "TRA")
    )
  histograms_trb <-
    suppressWarnings(
      plot_histograms(summary_df_trb, folder = threshold, file_ext = "TRB")
    )
  remove(histograms_tra, histograms_trb)
  cat("\nPlotting histograms for TRA and TRB. Done.")
  
  ####################################
  #### Richness per cell Barplots ####
  ####################################
  summary_df_tra$chain <- "TRA"; summary_df_trb$chain <- "TRB"
  combined_summary <- bind_rows(summary_df_tra, summary_df_trb, .id = "CCHAIN")

  combined_summary$cell_main <- combined_summary$Cell
  # order subsets
  combined_summary$organ <- factor(combined_summary$organ,
    levels = c("Blood", "BM", "Skin", "Fat")
  )
  
  lymphoid_richness_per_cell_tra <-
    plot_richness_bars(chain_to_plot = "TRA", lymphoid_group = "lymphoid")
  lymphoid_richness_per_cell_trb <- 
    plot_richness_bars(chain_to_plot = "TRB", lymphoid_group = "lymphoid")
  non_lymphoid_richness_per_cell_tra <- 
    plot_richness_bars(chain_to_plot = "TRA", lymphoid_group = "non_lymphoid")
  non_lymphoid_richness_per_cell_trb <- 
    plot_richness_bars(chain_to_plot = "TRB", lymphoid_group = "non_lymphoid")
  
  
  richness_lymphoid_tra_path <-
    file.path("./figures/clones_per_cell",
    threshold, "lymphoid_clones_per_cell_tra.pdf"
  )
  richness_lymphoid_trb_path <-
    file.path("./figures/clones_per_cell", threshold, "lymphoid_clones_per_cell_trb.pdf")
  richness_nonlymphoid_tra_path <-
    file.path("./figures/clones_per_cell", threshold, "non_lymphoid_clones_per_cell_tra.pdf")
  richness_nonlymphoid_trb_path <-
    file.path("./figures/clones_per_cell", threshold, "non_lymphoid_clones_per_cell_trb.pdf")
  
  ggsave(richness_lymphoid_tra_path, lymphoid_richness_per_cell_tra,
    create.dir = TRUE, width = 30, height = 16, unit = "cm"
  )
  ggsave(richness_lymphoid_trb_path, lymphoid_richness_per_cell_trb,
    create.dir = TRUE, width = 30, height = 16, unit = "cm"
  )
  ggsave(richness_nonlymphoid_tra_path, non_lymphoid_richness_per_cell_tra,
    create.dir = TRUE, width = 30, height = 16, unit = "cm"
  )
  ggsave(richness_nonlymphoid_trb_path, non_lymphoid_richness_per_cell_trb,
    create.dir = TRUE, width = 30, height = 16, unit = "cm"
  )

  #####################################
  #### Clone Distribution Barplots ####
  #####################################
  # Plot grouped abundances in bar plots
  grouped_abundances_df <-
    bind_rows(
      lapply(names(clone_dfs_tra), function(x) {
        abundance_dataframe <- group_abundances(clone_abundance_tra[[x]])
        abundance_dataframe <-
          bind_cols(
            abundance_dataframe,
            metadata %>%
              filter(gs.id == x) %>% select(individual, organ, subset, cell.count, gs.id),
            chain = "TRA"
          )
        return(abundance_dataframe)
      }),
      lapply(names(clone_dfs_tra), function(x) {
        abundance_dataframe <- group_abundances(clone_abundance_trb[[x]])
        abundance_dataframe <-
          bind_cols(
            abundance_dataframe,
            metadata %>%
              filter(gs.id == x) %>% select(individual, organ, subset, cell.count, gs.id),
            chain = "TRB"
          )
        return(abundance_dataframe)
      })
    )
  
  # fill in cell for downstream
  grouped_abundances_df$Cell <- grouped_abundances_df$subset
  grouped_abundances_df$Group <- factor(grouped_abundances_df$Group, levels = c("1-10", "11-100", "101-1000", ">1000"))
  
  grouped_abundance_tra_plot <-
    ggplot(
      grouped_abundances_df %>% filter(chain == "TRA"),
      aes(fill = Group, y = Abundance_percentage_group, x = Cell)
    ) +
    geom_bar(position = "stack", stat = "identity") +
    facet_grid(rows = vars(individual), cols = vars(organ), scales = "free", space = "free") +
    ylab("% pool size") +
    scale_fill_grey(start = 0.1, end = 0.9) +
    my_theme() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  
  
  grouped_abundance_trb_plot <- 
    ggplot(grouped_abundances_df%>%filter(chain == "TRB"), aes(fill=Group, y=Abundance_percentage_group, x=Cell)) + 
    geom_bar(position="stack", stat="identity") + 
    facet_grid(rows = vars(individual),cols = vars(organ), scales='free', space = "free") + 
    ylab("% pool size") +
    scale_fill_grey(start = 0.1, end = 0.9) +
    my_theme() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  
  # save plots
  clone_dist_path <- file.path("./figures/clone_dist", threshold)
  ggsave(
    filename = paste0(clone_dist_path, "/binned_abundance_tra.pdf"), plot = grouped_abundance_tra_plot,
    width = 13, height = 8
  )
  ggsave(
    filename = paste0(clone_dist_path, "/binned_abundance_trb.pdf"), plot = grouped_abundance_trb_plot,
    width = 13, height = 8
  )
  
  ## Combine all clonetypes and plot individual clone frequencies
  abundances_df <-
    bind_rows(
      lapply(names(clone_dfs_tra), function(x) {
        abundance_dataframe <- clone_abundance_tra[[x]]
        abundance_dataframe_f <- 100 * abundance_dataframe / sum(abundance_dataframe)
        abundance_dataframe <-
          bind_cols(
            Abundance = abundance_dataframe_f,
            metadata %>%
              filter(gs.id == x) %>% select(individual, organ, subset, cell.count),
            chain = "TRA",
            UMI_count = clone_abundance_tra[[x]],
            gs.id = x,
            nSeqCDR3 = clone_dfs_tra[[x]]$nSeqCDR3
          )
        return(abundance_dataframe)
      }),
      lapply(names(clone_dfs_trb), function(x) {
        abundance_dataframe <- clone_abundance_trb[[x]]
        abundance_dataframe_f <- sort(100 * abundance_dataframe / sum(abundance_dataframe), decreasing = TRUE)
        abundance_dataframe <-
          bind_cols(
            Abundance = abundance_dataframe_f,
            metadata %>%
              filter(gs.id == x) %>% select(individual, organ, subset, cell.count),
            chain = "TRB",
            UMI_count = clone_abundance_trb[[x]],
            gs.id = x,
            nSeqCDR3 = clone_dfs_trb[[x]]$nSeqCDR3
          )
        return(abundance_dataframe)
      })
    )

  remove(clone_dist_plots_tra)
  remove(clone_dist_plots_trb)
  # Create IDs for clones
  # abundances_df$clone_id <- factor(seq_len(nrow(abundances_df)))
  
  # Fix cell names to print better at plot
  abundances_df$Cell <- abundances_df$subset

  # clone ids
  abundances_df <- 
    abundances_df %>% group_by(individual, organ, Cell, chain) %>% 
    mutate(clone_id_pg = row_number()) %>%
    ungroup() %>% 
    mutate(nSeqCDR3Encoded = as.factor(nSeqCDR3)) %>%
    mutate(highlight_10 = clone_id_pg <= 10)

  # samples to be filtered
  low_umi_samples <- 
    combined_summary %>%
    filter(UMI_count<500) %>%
    select(individual, organ, chain, gs.id)
  
  # filter low quality samples and save clonotype abundances of all individuals
  write.csv(
    x = abundances_df %>%
      anti_join(low_umi_samples, by = c("individual", "organ", "chain", "gs.id")) %>%
      select(individual, organ, Cell, chain, Abundance, nSeqCDR3),
    paste0("./data/human.clonetypes.", threshold, ".filtered.csv"),
    row.names = FALSE
  )
  # HERE!!!
  individual_clone_abundance_tra_plot <-
    ggplot(
      abundances_df %>% filter(chain == "TRA"),
      aes(
        y = Abundance, x = Cell,
        color = highlight_10, group = -clone_id_pg
      )
    ) +
    geom_col(
      position = "stack",
      fill = "white",
      linewidth = 0.05
    ) +
    scale_color_manual(
      values = c("FALSE" = "black", "TRUE" = "darkred"),
      guide = "none"
    ) +
    facet_grid(
      rows = vars(individual), cols = vars(organ),
      scales = "free", space = "free_x"
    ) +
    ylab("% pool size") +
    my_theme(font_size = 6) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  individual_clone_abundance_trb_plot <-
    ggplot(
      abundances_df %>% filter(chain == "TRB"),
      aes(
        y = Abundance, x = Cell,
        color = highlight_10, group = -clone_id_pg
      )
    ) +
    geom_col(
      position = "stack",
      fill = "white",
      linewidth = 0.05
    ) +
    scale_color_manual(
      values = c("FALSE" = "black", "TRUE" = "darkred"),
      guide = "none"
    ) +
    facet_grid(rows = vars(individual), cols = vars(organ),
      scales = "free", space = "free_x"
    ) +
    ylab("% pool size") +
    my_theme(font_size = 6) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  # save plots
  ggsave(
    filename = paste0(clone_dist_path, "/binless_abundance_tra.pdf"),
    plot = individual_clone_abundance_tra_plot,
    width = 18.4, height = 8, unit = "cm"
  )
  ggsave(
    filename = paste0(clone_dist_path, "/binless_abundance_trb.pdf"),
    plot = individual_clone_abundance_trb_plot,
    width = 18.4, height = 8, unit = "cm"
  )
  
  
  # Select the top 10 clone abundance
  abundances_df_summed <-
    abundances_df %>%
    filter(clone_id_pg <= 10) %>%
    group_by(individual, organ, Cell, chain) %>%
    mutate(Abundance_top_10 = sum(Abundance)) %>%
    filter(clone_id_pg == 10)

  colnames(combined_summary)

  abundances_df_summed <- 
  left_join(
      combined_summary,
      abundances_df_summed[c("Abundance_top_10", "chain", "gs.id")],
      by = c("gs.id", "chain")
    )

  # Export combined metrics
  summary.path <- file.path("./data/summary", threshold)
  if (!dir.exists(summary.path)) {
    dir.create(summary.path, recursive = TRUE)
  }
  write.csv(
    x = abundances_df_summed,
    file = paste0(summary.path, "/individual_tissue_tcr_sample_metrics.csv"),
    row.names = FALSE)
  shape_pool <- c(21, 22, 23, 24, 25,16, 17, 18, 15, 0, 1, 2, 3, 4, 8, 11)
  lymphoid_clone_abundance_summed_plot <-
    ggplot(
      abundances_df_summed %>%
        filter(
          organ %in% organ_groups$lymphoid,
          # cell_count != 20000,
          UMI_count > 500
        ),
      aes(
        y = Abundance_top_10, x = Cell,
        fill = chain
      )
    ) +
    stat_summary(
      fun = "median", geom = "col",
      width = 0.9, linewidth = 0.4,
      #position = position_dodge(width = 0.9)
      fill = NA,
      colour = rgb(0, 0, 0, 0.7),
      alpha = 0.7
    ) +
    geom_line(aes(group = individual),
      linewidth = 0.3
    ) +
    geom_point(
      aes(
        shape = individual,
        fill = CellOrUMI
      ),
      alpha = 1,
      size = 2
    ) +
    ggh4x::facet_nested(
      cols = vars(organ, cell_main), rows = vars(chain),
      scales = "free", space = "free",
      labeller = labeller(
        cell.main = function(x) rep("", length(x))
      )
    ) +
    # theme(strip.text.x = ggtext::element_markdown(size = c(7, 0))) +  # size 0 hides 2nd)
    ylab("% pool size of top 10 clones") +
    coord_cartesian(ylim = c(0, 50), clip = "off") +
    scale_y_continuous(expand = c(0, 0)) +
    xlab("") +
    my_theme(font_size = 7) +
    clone_abundance_summed_plot_theme +
    scale_fill_manual(values = c(umi = "white", "cell" = "gray")) +
    scale_shape_manual(
      values = shape_pool[1:length(unique(abundances_df_summed$individual))]
    ) +
    guides(
      shape = guide_legend(
      override.aes = list(fill = "white", color = "black")
      ),
      fill = "none",
      color = "none"
    ) +
    theme(panel.spacing.y = unit(0.5, "lines"),
          axis.title.x = element_blank(),
          legend.position = "bottom",
          legend.key.size = unit(0.2, "cm"), # shrink the boxes
          legend.spacing.x = unit(0.2, "cm"), # horizontal spacing between items
          legend.text = element_text(size = 6), # shrink text
          axis.ticks.x = element_blank()
        ) #+
    #labs(shape = subject_label)

  ggsave(
    filename = paste0(clone_dist_path, "/lymphoid_top10_pool_size.pdf"),
    plot = lymphoid_clone_abundance_summed_plot,
    width = 18.4, height = 8, units = "cm"
  )

  non_lymphoid_clone_abundance_summed_plot <-
    ggplot(
      abundances_df_summed %>%
        filter(organ %in% organ_groups$non_lymphoid,
        # cell_count != 20000,
        UMI_count > 500
        ),
      aes(
        y = Abundance_top_10, x = Cell,
        fill = chain
      )
    ) +
    stat_summary(
      fun = "median", geom = "col",
      width = 0.9, linewidth = 0.4,
      #position = position_dodge(width = 0.9)
      fill = NA,
      colour = rgb(0, 0, 0, 0.7),
      alpha = 0.7
    ) +
    geom_line(aes(group = individual),
      linewidth = 0.3
    ) +
    geom_point(
      aes(
        shape = individual,
        fill = CellOrUMI
      ),
      alpha = 1,
      size = 2
    ) +
    ggh4x::facet_nested(
      cols = vars(organ, cell_main), rows = vars(chain),
      scales = "free", space = "free",
      labeller = labeller(
        cell.main = function(x) rep("", length(x))
      )
    ) +
    # theme(strip.text.x = ggtext::element_markdown(size = c(7, 0))) +  # size 0 hides 2nd)
    ylab("% pool size of top 10 clones") +
    coord_cartesian(ylim = c(0, 80), clip = "off") +
    scale_y_continuous(expand = c(0, 0)) +
    xlab("") +
    my_theme(font_size = 7) +
    clone_abundance_summed_plot_theme +
    scale_fill_manual(values = c(umi = "white", "cell" = "gray")) +
    scale_shape_manual(
      values = shape_pool[1:length(unique(abundances_df_summed$individual))]
    ) +
    guides(
      fill = "none",
      color = "none"
    ) +
    theme(panel.spacing.y = unit(0.5, "lines"),
          axis.title.x = element_blank(),
          legend.position = "bottom",
          legend.key.size = unit(0.2, "cm"), # shrink the boxes
          legend.spacing.x = unit(0.2, "cm"), # horizontal spacing between items
          legend.text = element_text(size = 6), # shrink text
          axis.ticks.x = element_blank()
        ) #+
    #labs(shape = subject_label)
    
  ggsave(
    filename = paste0(clone_dist_path, "/nonLymphoid_top10_pool_size.pdf"),
    plot = non_lymphoid_clone_abundance_summed_plot,
    width = 13.8, height = 8, units = "cm"
  )
  
} # threshold loop end

####################
########Save########
####################
cat("\nSaving Rdat.")
save.image("./data/post_analysis.Rdat")
cat("Completed!")


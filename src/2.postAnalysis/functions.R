# post analysis functions
# subsamples UMIs
subsample_abundance <- function(data,
                                target,
                                target_type = "umi",
                                name_col = "nSeqCDR3",
                                abundance_col = "uniqueMoleculeCount",
                                keep_all = FALSE,
                                seed = 123) {
  if (!is.null(seed))
    set.seed(seed)
  
  # If data is a data.frame, extract abundances
  if (is.data.frame(data)) {
    abundances <- data[[abundance_col]]
    names(abundances) <- data[[name_col]]
  } else {
    abundances <- data
  }
  
  # calcualte basic metrics
  total_abundance <- sum(abundances)
  total_richness <- length(abundances)
  
  # sample based on umis
  if (target_type == "umi") {
    if (total_abundance < target) {
      stop("Target is larger than total abundance. Can not sample.")
    }
    
    # keeps all clones
    if (keep_all) {
      if (target < length(abundances)) {
        stop(
          "Target is smaller than number of species, cannot retain all species with at least one count."
        )
      }
      
      min_vec <- rep(1, length(abundances))
      names(min_vec) <- names(abundances)
      
      to_sample <- abundances - min_vec
      to_sample[to_sample < 0] <- 0
      
      remaining <- target - length(abundances)
      
      if (remaining > 0) {
        sampled <- table(sample(rep(names(abundances), to_sample), remaining))
        out_counts <- min_vec
        out_counts[names(sampled)] <- out_counts[names(sampled)] + as.numeric(sampled)
      } else {
        out_counts <- min_vec
      }
    } else {
      out_counts <- table(sample(rep(names(abundances), abundances), target))
    }
    return(out_counts)
  }
  
  # target is clones
  if (target_type == "clone") {
    if (target >= total_richness){
      return(table(rep(names(abundances), abundances)))
    }
    # shuffle the clones
    shuffled <- sample(rep(names(abundances), abundances))
    # get the position where each unique clone appears first
    idx_first_appearance <- match(unique(shuffled), shuffled)
    # sort those positions, the tarhet unique will be reached at position idx_first_appearance[target]
    stop_idx <- sort(idx_first_appearance)[target]
    return(table(shuffled[1:stop_idx]))
  }
}

subsample_sequence <- function(sequences, target, seed = 123) {
  if (!is.null(seed))
    set.seed(seed)
  
  richness <- length(sequences)
  if (richness < target) {
    return(sequences)
  }
  else {
    sampled <- sample(sequences, target, replace = FALSE)
    return(sampled)
  }
}

# split clones file into TRA and TRB
## split TRA and TRB
split_clones <- function(clones_df) {
  if (sum(is.na(clones_df$allJHitsWithScore)) > 0)
    stop("Not all rows are valid")
  # clones_df <-  clone_dfs[[49]]
  clones_df$chain <- substr(clones_df$allJHitsWithScore,
                            start = 0,
                            stop = 3)
  clones_tra <- clones_df[clones_df$chain == "TRA", ]
  clones_trb <- clones_df[clones_df$chain == "TRB", ]
  #combine possible nSeqCDR3s!
  clones_tra <-
    clones_tra %>%
    group_by(nSeqCDR3) %>%
    summarise(
      readCount = sum(readCount),
      readFraction = sum(readFraction),
      uniqueMoleculeCount = sum(uniqueMoleculeCount),
      uniqueMoleculeFraction = sum(uniqueMoleculeFraction),
      chain = first(chain)
    ) %>%
    ungroup() %>%
    arrange(-uniqueMoleculeCount) %>%
    mutate(cloneId = row_number())
  
  clones_trb <-
    clones_trb %>%
    group_by(nSeqCDR3) %>%
    summarise(
      readCount = sum(readCount),
      readFraction = sum(readFraction),
      uniqueMoleculeCount = sum(uniqueMoleculeCount),
      uniqueMoleculeFraction = sum(uniqueMoleculeFraction),
      chain = first(chain)
    ) %>%
    ungroup() %>%
    arrange(-uniqueMoleculeCount) %>%
    mutate(cloneId = row_number())

  clone_abundance_tra <- clones_tra$uniqueMoleculeCount
  clone_abundance_trb <- clones_trb$uniqueMoleculeCount
  names(clone_abundance_tra) <- clones_tra$cloneId
  names(clone_abundance_trb) <- clones_trb$cloneId
  ## number of clones
  number_of_clones_tra <- length(clone_abundance_tra)
  number_of_clones_trb <- length(clone_abundance_trb)
  return(
    list(
      clones_tra = clones_tra,
      clones_trb = clones_trb,
      clone_abundance_tra = sort(clone_abundance_tra, decreasing = TRUE),
      clone_abundance_trb = sort(clone_abundance_trb, decreasing = TRUE)
    )
  )
}

summary_clones <- function(clones_df) {
  sum(clones_df$readCount)
  sum(clones_df$uniqueMoleculeCount)
  nrow(clones_df)
  return(list(
    clone_count = nrow(clones_df),
    UMI_count = sum(clones_df$uniqueMoleculeCount),
    read_count = sum(clones_df$readCount)
  ))
}

# Calculate (umi/cel etc.) of combined samples df
# calculate summary stats
per_cell_metrics <- function(df) {
  df$ReadperCell <- df$read_count / df$cell.count
  df$ReadperUMI <- df$read_count / df$UMI_count
  df$UMIperClone <- df$UMI_count / df$clone_count
  df$ReadperClone <- df$read_count / df$clone_count
  df$UMIperCell <- df$UMI_count / df$cell.count
  df$clonePerCell <- df$clone_count / df$cell.count
  df$CellOrUMI <- 
    ifelse(
      df$UMI_count < df$cell.count,
      "umi",
      "cell")
  df$clonePerCellUMI <- ifelse(df$UMI_count < df$cell.count,
    df$clone_count / df$UMI_count,
    df$clone_count / df$cell.count
  )
  return(df)
}

# Summarizes samples into dataframes for TRA and TRB
# Adds metadata information
# Calculates further metrics: UMIperCell etc.
summarize_samples <-
  function(clone_dfs_2_tra, clone_dfs_2_trb) {
    summarylist_tra <- sapply(clone_dfs_2_tra, summary_clones)
    summarylist_trb <- sapply(clone_dfs_2_trb, summary_clones)
    
    summary_df_tra <- data.frame(t(summarylist_tra), chain = "TRA")
    summary_df_trb <- data.frame(t(summarylist_trb), chain = "TRB")
    
    summary_df_tra[, c("cell_numbers", "organ")] <-
      metadata[, c("cell.numbers", "organ")]
    summary_df_trb[, c("cell_numbers", "organ")] <-
      metadata[, c("cell.numbers", "organ")]
    
    # re-merge as a df to plot
    ## TRA
    summary_df_tra <- data.frame(
      clone_count = unlist(summary_df_tra[["clone_count"]]),
      UMI_count = unlist(summary_df_tra[["UMI_count"]]),
      read_count = unlist(summary_df_tra[["read_count"]], ),
      cell_numbers = unlist(summary_df_tra[["cell_numbers"]]),
      organ = unlist(summary_df_tra[["organ"]])
    )
    ##TRB
    summary_df_trb <- data.frame(
      clone_count = unlist(summary_df_trb[["clone_count"]]),
      UMI_count = unlist(summary_df_trb[["UMI_count"]]),
      read_count = unlist(summary_df_trb[["read_count"]], ),
      cell_numbers = unlist(summary_df_trb[["cell_numbers"]]),
      organ = unlist(summary_df_trb[["organ"]])
    )
    ## add per_cell_metrics
    summary_df_tra <- per_cell_metrics(summary_df_tra)
    summary_df_trb <- per_cell_metrics(summary_df_trb)
    
    # add cell ids etc. from metadata
    summary_df_tra[c("mouse", "organ", "cell", "cell count")] <-
      metadata[c("mouse", "organ", "cell.type", "cell.numbers")]
    
  }

# Plot histograms to explore samples
plot_histograms <-
  function(summary_df, folder, file_ext, font_size_on_bar = 2, font_size_on_x = 6) {

    cells_unique <- unique(summary_df$subset)

    n.cells <- length(cells_unique)
    # Define colors
    if (n.cells <= 9) {
      cell_colors <- brewer.pal(n.cells, "Set1")
    } else {
      cell_colors <- colorRampPalette(brewer.pal(9, "Set1"))(n.cells)
    }
    names(cell_colors) <- cells_unique
    summary_df$cell <- factor(summary_df$subset, levels=cells_unique)
    
    # Plot histograms
    # test summary_df <- summary_df_tra
    cell.count_hist <-
      ggplot(data = summary_df, aes(x = subset, y = cell.count, fill = cell)) +
      geom_bar(position = "dodge", stat = "identity") +
      facet_grid(rows = vars(individual), cols = vars(organ), scales = "free_x", space = "free_x") +
      ylab("# of cells") + xlab("") +
      scale_y_continuous() +
      my_theme(font_size = font_size_on_x) +
      scale_fill_manual(values = cell_colors) +
      geom_text(
        size = font_size_on_bar,
        aes(label = round(cell.count, 2)),
        vjust = 0,
        position = position_dodge(width = 0.9)
      ) +
      guides(fill=guide_legend(nrow=2, byrow=FALSE)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
    
    
    read_count_hist <-
      ggplot(data = summary_df, aes(x = subset, y = read_count, fill = cell)) +
      geom_bar(position = "dodge", stat = "identity") +
      facet_grid(rows = vars(individual), cols = vars(organ), scales = "free_x", space = "free_x") +
      ylab("# of reads") + xlab("") +
      scale_y_continuous() +
      my_theme(font_size = font_size_on_x) +
      scale_fill_manual(values = cell_colors) +
      geom_text(
        size = font_size_on_bar,
        aes(label = round(read_count, 2)),
        vjust = 0,
        position = position_dodge(width = 0.9)
      ) + 
      guides(fill=guide_legend(nrow=2, byrow=FALSE)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
    
    read_count_mouse_hist <-
      ggplot(data = summary_df, aes(x = individual, y = read_count, fill = cell)) +
      geom_bar(position = "stack", stat = "identity") +
      ylab("# of reads") + xlab("") +
      scale_y_continuous() +
      my_theme(font_size = font_size_on_x) +
      scale_fill_manual(values = cell_colors) +
      guides(fill=guide_legend(nrow=2, byrow=FALSE)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
    
    umi_count_hist <-
      ggplot(data = summary_df, aes(x = subset, y = UMI_count, fill = cell)) +
      geom_bar(position = "dodge", stat = "identity", na.rm = TRUE) +
      facet_grid(rows = vars(individual), cols = vars(organ), scales = "free_x", space = "free_x") +
      ylab("# of UMIs") + xlab("") +
      my_theme(font_size = font_size_on_x) +
      scale_fill_manual(values = cell_colors) +
      geom_text(
        size = font_size_on_bar,
        aes(label = round(UMI_count, 2)),
        vjust = 0,
        position = position_dodge(width = 0.9),
        na.rm = TRUE
      ) + 
      guides(fill=guide_legend(nrow=2, byrow=FALSE)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
    
    UMI_count_mouse_hist <-
      ggplot(data = summary_df, aes(x = individual, y = UMI_count, fill = subset)) +
      geom_bar(position = "stack", stat = "identity") +
      ylab("# of UMIs") + xlab("") +
      my_theme(font_size = font_size_on_x) +
      scale_fill_manual(values = cell_colors) + 
      guides(fill=guide_legend(nrow=2, byrow=FALSE)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
    
    cell.count_mouse_hist <-
      ggplot(data = summary_df, aes(x = individual, y = cell.count, fill = subset)) +
      geom_bar(position = "stack", stat = "identity") +
      ylab("# of cells") + xlab("") +
      my_theme(font_size = font_size_on_x) +
      scale_fill_manual(values = cell_colors) + 
      guides(fill=guide_legend(nrow=2, byrow=FALSE)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
    
    clones_hist <-
      ggplot(data = summary_df, aes(x = subset, y = clone_count, fill = subset)) +
      geom_bar(position = "dodge", stat = "identity", na.rm = TRUE) +
      facet_grid(rows = vars(individual), cols = vars(organ), scales = "free_x", space = "free_x") +
      geom_text(
        size = font_size_on_bar,
        aes(label = clone_count),
        vjust = 0,
        position = position_dodge(width = 0.9),
        na.rm = TRUE
      ) +
      ylab("# of clones") + xlab("") +
      my_theme(font_size = font_size_on_x) +
      scale_fill_manual(values = cell_colors) + 
      guides(fill=guide_legend(nrow=2, byrow=FALSE)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
    
    readPerCell_hist <-
      ggplot(data = summary_df, aes(x = cell, y = ReadperCell, fill = cell)) +
      geom_bar(position = "dodge", stat = "identity") +
      facet_grid(rows = vars(individual), cols = vars(organ), scales = "free_x", space = "free_x") +
      ylab("Read/Cell") + xlab("") +
      my_theme(font_size = font_size_on_x) +
      scale_fill_manual(values = cell_colors) +
      geom_text(
        size = font_size_on_bar,
        aes(label = round(ReadperCell, 2)),
        vjust = 0,
        position = position_dodge(width = 0.9)
      ) + 
      guides(fill=guide_legend(nrow=2, byrow=FALSE)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
    
    UMIperCell_hist <-
      ggplot(data = summary_df, aes(x = subset, y = UMIperCell, fill = subset)) +
      # geom_hline(yintercept = 1,
      #            linetype = "dashed",
      #            color = "red") +
      geom_bar(position = "dodge", stat = "identity") +
      facet_grid(rows = vars(individual), cols = vars(organ), scales = "free_x", space = "free_x") +
      geom_text(
        size = font_size_on_bar,
        aes(label = round(UMIperCell, 2)),
        vjust = 0,
        position = position_dodge(width = 0.9)
      ) +
      ylab("UMI/Cell") + xlab("") +
      my_theme(font_size = font_size_on_x) +
      scale_fill_manual(values = cell_colors) +
      # scale_y_continuous(breaks = seq(0, 7)) + 
      guides(fill=guide_legend(nrow=2, byrow=FALSE)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
    
    clonePerCell_hist <-
      ggplot(data = summary_df, aes(x = subset, y = clonePerCell, fill = subset)) +
      # geom_hline(yintercept = 1,
      #            linetype = "dashed",
      #            # color = "red"
      #            ) +
      geom_bar(position = "dodge", stat = "identity") +
      facet_grid(rows = vars(individual), cols = vars(organ), scales = "free_x", space = "free_x") +
      geom_text(
        size = font_size_on_bar,
        aes(label = round(clonePerCell, 2)),
        vjust = 0,
        position = position_dodge(width = 0.9)
      ) +
      ylab("Clone/Cell") + xlab("") +
      my_theme(font_size = font_size_on_x) +
      scale_fill_manual(values = cell_colors) + 
    guides(fill=guide_legend(nrow=2, byrow=FALSE)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
    # scale_y_continuous(breaks = seq(0, 1))
    
    clonePerCellorUMI_hist <-
      ggplot(data = summary_df, aes(x = cell, y = clonePerCellUMI, fill = cell)) +
      # geom_hline(yintercept = 1,
      #            linetype = "dashed",
      #            # color = "red"
      #            ) +
      geom_bar(position = "dodge", stat = "identity") +
      facet_grid(rows = vars(individual), cols = vars(organ), scales = "free_x", space = "free_x") +
      geom_text(
        size = font_size_on_bar,
        aes(label = round(clonePerCellUMI, 2)),
        vjust = 0,
        position = position_dodge(width = 0.9)
      ) +
      ylab("Clone/(Cell | UMI)") + xlab("") +
      my_theme(font_size = font_size_on_x) +
      scale_fill_manual(values = cell_colors) + 
      guides(fill=guide_legend(nrow=2, byrow=FALSE)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
    
    # Scatterplots
    UMIperCell_CellCount <-
      ggplot(data = summary_df, aes(x = cell.count, y = UMIperCell, color = subset)) +
      geom_hline(yintercept = 1,
                linetype = "dashed",
                color = "red") +
      geom_point(size = 3) +
      facet_grid(cols = vars(individual)) +
      ylab("UMI/Cell") + xlab("# of cells") +
      my_theme(font_size = font_size_on_x) +
      scale_color_manual(values = cell_colors) +
      scale_y_continuous(breaks = seq(0, 5)) + 
      guides(fill=guide_legend(nrow=2, byrow=FALSE)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
    
    clones_CellCount <-
      ggplot(data = summary_df, aes(x = cell.count, y = clone_count, color = subset)) +
      geom_abline(intercept = 0, slope = 1, color="red", linetype = "dashed") +
      geom_point(size = 2, alpha = 0.7) +
      facet_grid(rows = vars(individual)) +
      ylab("Clones") + xlab("# of cells") +
      my_theme(font_size = font_size_on_x) +
      scale_color_manual(values = cell_colors) + 
      stat_regline_equation(label.x = 2e05, label.y = 20000) + 
      guides(fill=guide_legend(nrow=2, byrow=FALSE))
    
    clones_CellCountBw <-
      ggplot(data = summary_df, aes(x = cell.count, y = clone_count)) +
      geom_abline(intercept = 0, slope = 1, color="red", linetype = "dashed") +
      geom_point(size = 2, alpha = 0.7) +
      facet_grid(rows = vars(individual)) +
      ylab("Clones") + xlab("# of cells") +
      my_theme(font_size = font_size_on_x) +
      scale_color_manual(values = cell_colors) + 
      stat_regline_equation(label.x = 2e05, label.y = 20000)
    
    clones_CellCountBw <- ggscatter(
      summary_df, x = "cell.count", y = "clone_count",
      add = "reg.line", alpha = 0.7
    ) +
      geom_abline(intercept = 0, slope = 1, color="red", linetype = "dashed") +
      facet_wrap(~individual) +
      stat_cor(label.y = 20000, label.x = 1.4e05) +
      stat_regline_equation(label.y = 23000, label.x = 1.4e05)
    
    
    clones_CellCountColor <- ggscatter(
      summary_df, x = "cell.count", y = "clone_count", alpha = 0.7, color = "cell", shape = "individual"
    ) +
      geom_abline(intercept = 0, slope = 1, color="red", linetype = "dashed") +
      stat_cor(label.y = 20000, label.x = 1.4e05) +
      stat_regline_equation(label.y = 23000, label.x = 1.4e05) + 
      guides(fill=guide_legend(nrow=1, byrow=FALSE))
    
    unique_tissues <- c("blood", "bm", "liver", "lung", "pp", "skin", "spleen")
    # Assign distinct shapes (0 to 6 example shapes)
    my_shapes <- c(0, 1, 2, 3, 4, 5, 6)
    clones_CellCountColorCell <- 
      ggscatter(summary_df, x = "cell.count", y = "clone_count", alpha = 0.7,
                color = "cell", shape = "organ") +
      geom_abline(intercept = 0, slope = 1,
                  color = "red", linetype = "dashed") +
      facet_wrap(~individual,scales = "free") +
      stat_cor(label.x.npc = "right", label.y.npc = "bottom", hjust = 1, vjust = 1, size = 2.5) +
      stat_regline_equation(size = 2.5, label.x.npc = "right", label.y.npc = "bottom", hjust = 1, vjust = 0) +
      geom_smooth(method = "lm", se = FALSE, color = "gray", alpha = 0.7) +  # adds regression line per facet
      scale_shape_manual(values = setNames(my_shapes, unique_tissues))
    
    # list of plots
    plots <-
      list(
        cell.count_histA = cell.count_hist,
        read_count_hist = read_count_hist,
        umi_count_hist = umi_count_hist,
        clones_hist = clones_hist,
        cell.count_mouse_hist = cell.count_mouse_hist,
        read_count_mouse_hist = read_count_mouse_hist,
        UMI_count_mouse_hist = UMI_count_mouse_hist,
        readPerCell_hist = readPerCell_hist,
        UMIperCell_hist = UMIperCell_hist,
        clonePerCell_hist = clonePerCell_hist,
        clonePerCellorUMI_hist = clonePerCellorUMI_hist,
        UMIperCell_CellCount = UMIperCell_CellCount,
        clones_CellCount = clones_CellCount,
        clones_CellCountBw = clones_CellCountBw,
        clones_CellCountColor = clones_CellCountColor,
        clones_CellCountColorCell = clones_CellCountColorCell
      )
    
    # save plots
    dir_path <- file.path("./figures/summary/histogram", folder, file_ext)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    
    for (name in names(plots)) {
      ggsave(
        filename = file.path(dir_path, paste0(name, "_", file_ext, ".pdf")),
        plot = plots[[name]],
        units = c("cm"),
        width = 24,
        height = 16
      )
    }
    ggsave(
      filename = file.path(dir_path, paste0("clones_CellCountBw", "_", file_ext, ".pdf")),
      plot = plots[["clones_CellCountBw"]],
      width = 24,
      height = 10
    )
    
    ggsave(
      filename = file.path(dir_path, paste0("clones_CellCountColor", "_", file_ext, ".pdf")),
      plot = plots[["clones_CellCountColor"]],
      width = 8,
      height = 8
    )
    
    ggsave(
      filename = file.path(dir_path, paste0("clones_CellCountColorCell", "_", file_ext, ".pdf")),
      plot = plots[["clones_CellCountColorCell"]],
      width = 12,
      height = 4
    )
    return(plots)
    
  }

plot_richness_bars <-
  function(chain_to_plot, lymphoid_group) {
    shape_pool <- c(21, 22, 23, 24, 25,16, 17, 18, 15, 0, 1, 2, 3, 4, 8, 11)
    ggplot(
      combined_summary %>%
        filter(
          organ %in% organ_groups[[lymphoid_group]],
          chain == chain_to_plot,
          # cell_count != 20000,
          UMI_count > 500
        ),
      aes(y = clonePerCellUMI, x = subset)
    ) +
      stat_summary(
        fun = "median", geom = "col",
        width = 0.9, linewidth = 0.4,
        fill = NA,
        colour = rgb(0, 0, 0, 0.7),
        alpha = 0.7
      ) +
      #geom_line(aes(group = individual),
      #  linewidth = 0.3
      #) +
      geom_point(
        aes(
          shape = individual,
          fill = CellOrUMI
        ),
        alpha = 1, size = 2
      ) +
      ggh4x::facet_nested(
        cols = vars(organ, cell_main),
        scales = "free", space = "free_x",
        labeller = labeller(cell_main = function(x) rep("", length(x)))
      ) +
      ylab("Clone per Cell (or UMI)") +
      my_theme(font_size = 6) +
      clone_abundance_summed_plot_theme +
      scale_fill_manual(values = c(cell = "gray", umi = "white")) +
      scale_color_manual(values = c("white")) +
      scale_shape_manual(
        values = shape_pool[1:length(unique(combined_summary$individual))]
      ) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_cartesian(
        ylim = c(0, 2),
        clip = "off"
      ) +
      guides(
        fill = "none",
        color = "none"
      ) +
      theme(
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.2, "cm"), # shrink the boxes
        legend.spacing.x = unit(0.2, "cm"), # horizontal spacing between items
        legend.text = element_text(size = 6), # shrink text
        axis.ticks.x = element_blank()
      ) +
      labs(shape = "individual") + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  }

# Plot clonotype distribution
clone_distribution <-
  function(clones_txt,
           title = "",
    cell_number = cell.count,
    clone_size_th = 2,
    folder = threshold,
    file_ext = "TRB",
    fraction = FALSE) {
    # order
    #cell_number <- 1000
    # clones_txt <- clone_dfs_tra$`107293-001-013`
    #umi_abundance <- sort(rep, decreasing = TRUE)
    umi_abundance <-
      clones_txt[order(-clones_txt$uniqueMoleculeCount), "uniqueMoleculeCount"]
    # count umis
    umi_count <- sum(umi_abundance)
    # read count
    read_count <- sum(clones_txt$readCount)
    # clone
    clone_count <- nrow(clones_txt)
    # calculate fraction
    umi_abundance_fraction <-
      umi_abundance$uniqueMoleculeCount / umi_count
    
    # Sort fractions descending
    fractions_sorted <- sort(umi_abundance_fraction, decreasing = TRUE)
    
    # Calculate CDF values
    cdf_y <- seq(1, length(fractions_sorted)) / length(fractions_sorted)
    cdf_y <- rev(cumsum(rev(fractions_sorted)))
    
    #ci50
    index_d50 <- which(cumsum(fractions_sorted) >= 0.5)[1]
    # CI50 is the percentage of clones contributing to 50% of the population abundance
    ci50_value <- index_d50 / length(fractions_sorted) * 100
    # DI50 is the percentage of clones contributing to 50% of the population abundance
    di50_value <- index_d50
    
    df_abundance_rank <- data.frame(
      Abundance = umi_abundance$uniqueMoleculeCount,
      #Abundance_fraction = umi_abundance_fraction,
      Rank = seq(1:length(umi_abundance))
    )
    
    
    # Create data frame and group by fraction to plot CDF
    df <- data.frame(fraction = fractions_sorted)
    
    # Count occurrences for each unique fraction
    df_summary <- df %>%
      group_by(fraction) %>%
      summarise(count = n()) %>%
      arrange(desc(fraction)) %>%
      # Calculate CDF as fraction of clones with size >= fraction
      mutate(cdf = cumsum(count) / sum(count),
             cdf2 = (sum(count) - cumsum(count) + 1))
    
    # linear regression on logAbundance-logRank
    model_log <- lm(log10(Abundance) ~ log10(Rank),
                    data = df_abundance_rank %>% filter(Abundance>=clone_size_th))
    # linear regression on CDF
    model_log_cdf <- lm(log10(cdf) ~ log10(fraction), 
                        data = df_summary)
    
    # regression estimates saved in df to plot etc.
    regression_df <- data.frame(
      "intercept" = coef(model_log)[[1]],
      "slope" = coef(model_log)[[2]],
      "CI50" = ci50_value,
      "D50" = di50_value
    )
    regression_df_cdf <- data.frame(
      "intercept" = coef(model_log_cdf)[[1]],
      "slope" = coef(model_log_cdf)[[2]]
      #"Cell number" = cell_number
    )
    
    # plot abundance and slope
    y_text <- max(df_abundance_rank$Abundance) * 0.9
    
    plot_log_log <-
      ggplot(data = df_abundance_rank, aes(x = .data[["Rank"]], y = .data[["Abundance"]])) +
      geom_point() +
      scale_y_log10() +
      scale_x_log10() +
      xlab("Rank") +
      ylab("Clone size") +
      ggtitle(
        paste(
          title,
          "\ncells = ",
          cell_number,
          " reads = ",
          read_count,
          "\numis = ",
          umi_count,
          " clones = ",
          clone_count
        )
      ) +
      geom_abline(aes(
        intercept = model_log$coefficients[1],
        slope = model_log$coefficients[2]
      )) +
      geom_text(
        data = regression_df,
        aes(
          label = paste("y = ", round(slope, 3), "x + ", round(intercept, 3), sep =
                          " "),
          x = Inf,
          y = Inf
        ),
        hjust = "right",
        vjust = "top"
      ) +
      theme_minimal()
    
    # Plot CDF
    legend_df <- data.frame(
      label = c("Data line", "Model fit"),
      x = c(NA, NA),
      y = c(NA, NA)
    )
    plot_cdf <-
      ggplot(df_summary, aes(x = fraction, y = cdf)) +
      geom_line(aes(color = "Data line")) +
      scale_x_log10(
        labels = scales::trans_format("log10", scales::math_format(10^.x)),
        expand = c(0, 0),
        limits = c(10e-6, 10e-1)
      ) +
      scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      geom_abline(
        aes(
          intercept = model_log_cdf$coefficients[1],
          slope = model_log_cdf$coefficients[2],
          color = "Model fit"
        )
      ) +
      scale_color_manual(name = "Legend",
                         values = c("Data line" = "black", "Model fit" = "red")) +
      geom_text(
        data = regression_df_cdf,
        aes(
          label = paste("y = ", round(slope, 3), "x + ", round(intercept, 3), sep =
                          " "),
          x = Inf,
          y = Inf
        ),
        color = "red",
        hjust = "right",
        vjust = "top"
      ) +
      labs(x = "Relative abundance of T cell clone", y = "Probability (≥ clone frequency) \nFraction of clones with at least this abundance") +
      ggtitle(
        paste(
          title,
          "\ncells = ",
          cell_number,
          " reads = ",
          read_count,
          "\numis = ",
          umi_count,
          " clones = ",
          clone_count
        )
      ) +
      theme_classic() +
      annotation_logticks(sides = "bl", outside = TRUE) +
      coord_cartesian(clip = "off") +
      theme(
        plot.margin = margin(5, 5, 20, 5),
        axis.text.x = element_text(margin = margin(t = 10)),
        axis.text.y = element_text(margin = margin(r = 10))
      )
    
    # list plots
    plots <- list(plot_log_log = plot_log_log, plot_cdf = plot_cdf)
    slopes <- list(regression = regression_df, regression_cdf = regression_df_cdf)
    
    # save plots
    dir_path <- file.path("./figures/clone_dist", folder, file_ext)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    
    for (name in names(plots)) {
      ggsave(
        filename = file.path(dir_path, paste0(title, name, "_", file_ext, ".pdf")),
        plot = plots[[name]],
        width = 6,
        height = 4
      )
    }
    return(list(plots = plots, slopes = slopes))
  }

# Plot clonotype distribution for 80k and 20k samples
clone_distribution20vs80Cdf <-
  function(clone_dfs = clone_dfs_tra,
           file_ext = "TRA",
           mouse = "M20",
           cell = "CD4",
           title = "",
           clone_size_th = 10,
           folder = threshold,
           fraction = FALSE) {
    
    # grab 20k and 80k samples
    ids <- metadata %>%
      filter(
        grepl("k$", sample.id), mouse == .env$mouse,
        grepl(cell, cell.type)
      ) %>%
      pull(gs.id, cell.numbers)

    cur_dfs <- clone_dfs[ids]
    names(cur_dfs) <- names(ids)
  
    # umi abundances
    umi_abundance_20k <- sort(cur_dfs[["20000"]]$uniqueMoleculeCount, decreasing = TRUE)
    umi_abundance_80k <- sort(cur_dfs[["80000"]]$uniqueMoleculeCount, decreasing = TRUE)
    # umi_abundance_80k <- cur_dfs[["80000"]][order(-cur_dfs[["80000"]]$uniqueMoleculeCount), "uniqueMoleculeCount"]

    # count umis
    umi_count_20k <- sum(umi_abundance_20k)
    umi_count_80k <- sum(umi_abundance_80k)
    # read count
    read_count_20k <- sum(cur_dfs[["20000"]]$readCount)
    read_count_80k <- sum(cur_dfs[["80000"]]$readCount)
    # clone
    clone_count_20k <- nrow(cur_dfs[["20000"]])
    clone_count_80k <- nrow(cur_dfs[["80000"]])
    # calculate fraction
    umi_abundance_fraction_20k <-
      umi_abundance_20k / umi_count_20k
    umi_abundance_fraction_80k <-
      umi_abundance_80k / umi_count_80k
    
    # Sort fractions descending
    fractions_sorted_20k <- sort(umi_abundance_fraction_20k, decreasing = TRUE)
    fractions_sorted_80k <- sort(umi_abundance_fraction_80k, decreasing = TRUE)
    
    # Calculate CDF values
    cdf_y_20k <- seq(1, length(fractions_sorted_20k)) / length(fractions_sorted_20k)
    cdf_y_20k <- rev(cumsum(rev(fractions_sorted_20k)))
    cdf_y_80k <- seq(1, length(fractions_sorted_80k)) / length(fractions_sorted_80k)
    cdf_y_80k <- rev(cumsum(rev(fractions_sorted_80k)))
    
    df_abundance_rank_20k <- data.frame(
      Abundance = umi_abundance_20k,
      Rank = seq(1:length(umi_abundance_20k))
    )
    df_abundance_rank_80k <- data.frame(
      Abundance = umi_abundance_80k,
      Rank = seq(1:length(umi_abundance_80k))
    )
    
    # Create data frame and group by fraction to plot CDF
    df_20k <- data.frame(fraction = fractions_sorted_20k)
    df_80k <- data.frame(fraction = fractions_sorted_80k)
    
    # Count occurrences for each unique fraction
    df_summary_20k <- df_20k %>%
      group_by(fraction) %>%
      summarise(count = n()) %>%
      arrange(desc(fraction)) %>%
      # Calculate CDF as fraction of clones with size >= fraction
      mutate(cdf = cumsum(count) / sum(count),
            cdf2 = (sum(count) - cumsum(count) + 1))
    df_summary_80k <- df_80k %>%
      group_by(fraction) %>%
      summarise(count = n()) %>%
      arrange(desc(fraction)) %>%
      # Calculate CDF as fraction of clones with size >= fraction
      mutate(cdf = cumsum(count) / sum(count),
            cdf2 = (sum(count) - cumsum(count) + 1))
    
    # linear regression on logAbundance-logRank
    model_log_20k <- lm(log10(Abundance) ~ log10(Rank),
                    data = data.frame(df_abundance_rank_20k %>% filter(Abundance>=clone_size_th)))
    model_log_80k <- lm(log10(Abundance) ~ log10(Rank),
                        data = df_abundance_rank_80k %>% filter(Abundance>=clone_size_th))
    # linear regression on CDF
    model_log_cdf_20k <- lm(log10(cdf) ~ log10(fraction), 
                            data = df_summary_20k)
    model_log_cdf_80k <- lm(log10(cdf) ~ log10(fraction), 
                            data = df_summary_80k)
    
    # regression estimates saved in df to plot etc.
    regression_df_20k <- data.frame(
      "intercept" = coef(model_log_20k)[[1]],
      "slope" = coef(model_log_20k)[[2]]#,
      # "CI50" = ci50_value
    )
    regression_df_80k <- data.frame(
      "intercept" = coef(model_log_80k)[[1]],
      "slope" = coef(model_log_80k)[[2]]#,
      # "CI50" = ci50_value
    )
    regression_df_cdf_20k <- data.frame(
      "intercept" = coef(model_log_cdf_20k)[[1]],
      "slope" = coef(model_log_cdf_20k)[[2]]
      #"Cell number" = cell_number
    )
    regression_df_cdf_80k <- data.frame(
      "intercept" = coef(model_log_cdf_80k)[[1]],
      "slope" = coef(model_log_cdf_80k)[[2]]
      #"Cell number" = cell_number
    )
    
    # plot abundance and slope
    y_text <- max(df_abundance_rank_80k$Abundance) * 0.9
    
    plot_log_log <-
      ggplot(data = df_abundance_rank_20k, aes(x = .data[["Rank"]], y = .data[["Abundance"]])) +
      geom_point(color = "darkblue") +
      geom_point(data = df_abundance_rank_80k, aes(x = .data[["Rank"]], y = .data[["Abundance"]]),
                 color = "darkred", alpha = 0.6) +
      scale_y_log10() +
      scale_x_log10() +
      xlab("Rank") +
      ylab("Clone size") +
      ggtitle(
        paste(
          title,
          "\ncells = ",
          "20k & 80k" ,
          " reads = ",
          paste(read_count_20k, read_count_80k, sep = " - "),
          "\numis = ",
          paste(umi_count_20k,umi_count_80k, sep = " - "),
          " clones = ",
          paste(clone_count_20k,clone_count_80k, sep = " - ")
        )
      ) +
      geom_abline(aes(
        intercept = model_log_20k$coefficients[1],
        slope = model_log_20k$coefficients[2]
      ), color = "darkblue") +
      geom_abline(aes(
        intercept = model_log_80k$coefficients[1],
        slope = model_log_80k$coefficients[2]
      ), color = "darkred") +
      geom_text(
        data = regression_df_20k,
        aes(
          label = paste("y = ", round(slope, 3), "x + ", round(intercept, 3), sep =
                          " "),
          x = 1,
          y = 1
        ),
        hjust = "left",
        vjust = "bottom",
        color = "darkblue"
      ) +
      geom_text(
        data = regression_df_80k,
        aes(
          label = paste("y = ", round(slope, 3), "x + ", round(intercept, 3), sep =
                          " "),
          x = Inf,
          y = Inf
        ),
        hjust = "right",
        vjust = "top",
        color = "darkred"
      ) +
      theme_minimal()
    
    
    # Plot CDF
    legend_df <- data.frame(
      label = c("Data line", "Model fit"),
      x = c(NA, NA),
      y = c(NA, NA)
    )
    plot_cdf <-
      ggplot(df_summary_20k, aes(x = fraction, y = cdf)) +
      geom_line(color = "darkblue") +
      geom_line(data = df_summary_80k, aes(x = fraction, y = cdf), color = "darkred") +
      scale_x_log10(
        labels = scales::trans_format("log10", scales::math_format(10^.x)),
        expand = c(0, 0),
        limits = c(10e-7, 10e-1)
      ) +
      scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      geom_abline(
        aes(
          intercept = model_log_cdf_20k$coefficients[1],
          slope = model_log_cdf_20k$coefficients[2]
        ),
        color = "lightblue"
      ) +
      geom_abline(
        aes(
          intercept = model_log_cdf_80k$coefficients[1],
          slope = model_log_cdf_80k$coefficients[2]
        ),
        color = "brown"
      ) +
      geom_text(
        data = regression_df_cdf_20k,
        aes(label = paste("y = ", round(slope, 3), "x + ", round(intercept, 3), sep = " "),
          x = Inf,
          y = Inf
        ),
        color = "darkblue",
        hjust = "right",
        vjust = "top",
        size = 2.5
      ) +
      geom_text(
        data = regression_df_cdf_80k,
        aes(
          label = paste("y = ", round(slope, 3), "x + ", round(intercept, 3), sep = " "),
          x = 1e-5,
          y = 1e-5
        ),
        color = "darkred",
        hjust = "right",
        vjust = "top",
        size = 2.5
      ) +
      labs(x = "Relative abundance of T cell clone", 
           y = "Probability (≥ clone frequency) \nFraction of clones with at least this abundance") +
      ggtitle(
        paste(
          title,
          "\ncells = ",
          "20k & 80k" ,
          " reads = ",
          paste(read_count_20k, read_count_80k, sep = " - "),
          "\numis = ",
          paste(umi_count_20k,umi_count_80k, sep = " - "),
          " clones = ",
          paste(clone_count_20k,clone_count_80k, sep = " - ")
        )
      ) +
      theme_classic() +
      annotation_logticks(sides = "bl", outside = TRUE) +
      coord_cartesian(clip = "off") +
      theme(
        plot.margin = margin(5, 5, 20, 5),
        axis.text.x = element_text(margin = margin(t = 10)),
        axis.text.y = element_text(margin = margin(r = 10))
      )
    
    # list plots
    plots <- list(plot_log_log = plot_log_log, plot_cdf = plot_cdf)
    #slopes <- list(regression = regression_df, regression_cdf = regression_df_cdf)
    
    # save plots
    dir_path <- file.path("./figures/clone_dist", folder, file_ext)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    
    for (name in names(plots)) {
      ggsave(
        filename = file.path(dir_path, paste0(title, name, "_", file_ext, ".pdf")),
        plot = plots[[name]],
        width = 8,
        height = 6
      )
    }
    return(list(plots = plots))
  }

calc_power_exponent <-
  function(clone_sizes) {
    # Suppose your clone size data vector is named `clone_sizes`
    m <- displ$new(clone_sizes)  # assumes discrete power-law data
    
    # Estimate xmin (minimum value to consider in the fit)
    est <- estimate_xmin(m)
    m$setXmin(est)
    
    # Estimate the power-law exponent alpha
    est <- estimate_pars(m)
    m$setPars(est$pars)
    
    # Extract exponent
    power_law_exponent <- est$pars
    
    # Visualize the fit
    plot(m)
    lines(m, col = "red")
    return(m)
  }

infer_powerlaw_threshold <- function(x) {
  # threshold the min
  m <- displ$new(x)  # assumes discrete power-law data
  # Estimate xmin (minimum value to consider in the fit)
  est <- estimate_xmin(m)
  est <- est$xmin
  if (is.na(est)) {
    est <- 0
  }
  return(est)
}

# Group abundances per defined bins
group_abundances <- function(abundances) {
  sorted_abundances <- sort(abundances, decreasing = TRUE)
  abundance_dataframe <- data.frame(Abundance = sorted_abundances)
  abundance_dataframe <-
    abundance_dataframe %>% mutate(Group = case_when(
      row_number() >= 1 & row_number() <= 10 ~ "1-10",
      row_number() > 10 & row_number() <= 100 ~ "11-100",
      row_number() > 100 & row_number() <= 1000 ~ "101-1000",
      row_number() > 1000 ~ ">1000"
    )) %>% 
    group_by(Group) %>%
    summarise(
      Abundance_group = sum(Abundance, na.rm = TRUE)
    ) %>%
    mutate(
      Abundance_percentage_group = Abundance_group / sum(Abundance_group) * 100
    )
  return(abundance_dataframe)
}




plot_clone_overlap <-
  function(clone_dfs,
           chain,
           collapse_clones = TRUE,
           sample = "umi",
           sample_to = "20k",
           mouse = mouse,
           cell = "CD4",
           cluster_columns = FALSE) {
    # clone_dfs <- data
    # collapse clones or not
    if (collapse_clones == TRUE) {
      # calculate richness
      clone_counts <- lapply(clone_dfs, function(df) {
        # count number of clones
        length(unique(df$nSeqCDR3))
      })
      # calculate UMIs
      df_grouped_seqs <- lapply(clone_dfs, function(df) {
        # count number of clones
        df_grouped <-
          df %>%
          group_by(nSeqCDR3) %>%
          summarise(uniqueMoleculeCount = sum(uniqueMoleculeCount))
      })
      umi_counts <- lapply(df_grouped_seqs, function(df) {
        # count number of clones
        umi_counts <- df$uniqueMoleculeCount
      })
      clone_dfs <- df_grouped_seqs
    } else {
      clone_counts <- lapply(clone_dfs, nrow)
      umi_counts <- lapply(clone_dfs, function(x)
        x$uniqueMoleculeCount)
    }
    
    # convert dfs to list of lists
    # clone_dfs_mouse <- lapply(clone_dfs, function(x) {
    #   as.list(x)
    # })
    clone_dfs_mouse <- clone_dfs
    
    # number of clones row and clones will differ if sampled
    clone_counts_row <- clone_counts
    clone_counts_columns <- clone_counts
    umi_counts_row <- lapply(umi_counts, sum)
    umi_counts_columns <- umi_counts_row
    # calculate to what value to sample to
    if (sample %in% c("tcr", "umi")) {
      if (sample_to == "min") {
        subsample_id <- search_subsample_level(clone_counts)
        subsample_id_umi <- search_subsample_level_umi(umi_counts)
        subsample_to <- clone_counts[[subsample_id]]
        subsample_to_umi <- sum(umi_counts[[subsample_id_umi]])
        reference_row <- grep(subsample_id, names(clone_counts))
        reference_row_umi <- grep(subsample_id_umi, names(umi_counts))
        clone_counts_row[clone_counts_row > subsample_to] <- subsample_to
        umi_counts_row[umi_counts_row > subsample_to_umi] <- subsample_to_umi
      } else if (sample_to == "20k") {
        id_row <- metadata[grepl(paste(cell, "CD69- 20k"), metadata$sample.id) &
                             metadata$mouse == mouse, ]$gs.id
        subsample_to <- clone_counts[[id_row]]
        subsample_to_umi <- sum(umi_counts[[id_row]])
        reference_row <- grep(id_row, names(clone_counts))
        umi_counts_row[umi_counts_row > subsample_to_umi] <- subsample_to_umi
      }
    } else {
      reference_row <- 0
    }
    
    # get names to label
    sample_names <- metadata[metadata$gs.id %in% names(clone_dfs_mouse), ]$sample.id
    
    
    
    # count rows & create matrix with it
    dfs_count <- length(clone_dfs_mouse)
    overlap_matrix <-
      matrix(nrow = dfs_count, ncol = dfs_count)
    overlap_matrix_jaccard <-
      matrix(nrow = dfs_count, ncol = dfs_count)
    
    # fill in the matrix
    # TO DO: collapse_clone == FALSE can not be subsampled yet!!!
    ## Iterate over rows
    for (i in seq(1, dfs_count)) {
      # get the the row in vector_1
      if (collapse_clones == TRUE) {
        if (sample == "tcr") {
          vector_1 <- clone_dfs_mouse[[i]]$nSeqCDR3
          n_i <- length(vector_1)
          if (n_i > subsample_to) {
            vector_1 <- sample(vector_1, subsample_to)
            n_i <- subsample_to
          }
        } else if (sample == FALSE) {
          vector_1 <- clone_dfs_mouse[[i]]$nSeqCDR3
          n_i <- length(vector_1)
        } else if (sample == "umi") {
          vector_1 <- clone_dfs_mouse[[i]]$nSeqCDR3
          n_i <- length(vector_1)
          n_umi <- sum(clone_dfs_mouse[[i]]$uniqueMoleculeCount)
          if (n_umi > subsample_to_umi) {
            vector_1_umi <- subsample_abundance(clone_dfs_mouse[[i]], subsample_to_umi)
            vector_1 <- names(vector_1_umi)
            n_i <- length(vector_1)
            umi_counts_row[i] <- sum(vector_1_umi)
            clone_counts_row[i] <- n_i
          } else {
            vector_1 <- clone_dfs_mouse[[i]]$nSeqCDR3
            n_i <- length(vector_1)
          }
          
        }
      }
      
      # move to the columns while keeping row constant
      for (j in seq(1, dfs_count)) {
        # Find overlaps
        if (collapse_clones == FALSE) {
          # make TCRs unique in each list
          vector_1 <- make.unique(clone_dfs_mouse[[i]]$nSeqCDR3)
          vector_2 <- make.unique(clone_dfs_mouse[[j]]$nSeqCDR3)
          overlap <- length(intersect(vector_1, vector_2))
          n_i <- nrow(clone_dfs_mouse[[i]])
          n_j <- nrow(clone_dfs_mouse[[j]])
        }
        else if (collapse_clones == TRUE) {
          # collapse clones
          # collapsed vector_1 comes from outer loop
          vector_2 <- clone_dfs_mouse[[j]]$nSeqCDR3
          n_j <- length(vector_2)
        }
        # calculate overlap
        overlap <- length(intersect(vector_1, vector_2))
        overlap_matrix[i, j] <- (overlap / n_i) * 100
        union <- length(union(vector_1, vector_2))
        overlap_matrix_jaccard[i, j] <- overlap / union
      }
    }
    
    # Fix names of the rows and columns with updated values
    if (sample == "umi") {
      sample_names_row <- paste0(sample_names,
                                 "\n(",
                                 umi_counts_row,
                                 " - ",
                                 clone_counts_row,
                                 ")")
      sample_names_column <- paste0(sample_names,
                                    "\n(",
                                    umi_counts_columns,
                                    " - ",
                                    clone_counts_columns,
                                    ")")
    } else {
      sample_names_row <- paste0(sample_names, "\n(", clone_counts_row, ")")
      sample_names_column <- paste0(sample_names, "\n(", clone_counts_columns, ")")
    }
    rownames(overlap_matrix) <- str_to_upper(str_split_fixed(sample_names_row, " ", 2)[, 2])
    colnames(overlap_matrix) <- str_to_upper(str_split_fixed(sample_names_column, " ", 2)[, 2])
    rownames(overlap_matrix_jaccard) <- str_to_upper(str_split_fixed(sample_names_row, " ", 2)[, 2])
    colnames(overlap_matrix_jaccard) <- str_to_upper(str_split_fixed(sample_names_column, " ", 2)[, 2])
    
    # split blood
    sample_name_blood <- grep("BLOOD ", colnames(overlap_matrix), value = TRUE)
    
    overlap_matrix_blood <- overlap_matrix[, sample_name_blood]
    overlap_matrix_blood_jaccard <- overlap_matrix_jaccard[, sample_name_blood]
    
    # Plot matrixes
    # Raw overlap
    heatmap_plot <- plot_overlap_matrix(overlap_matrix, colnames(overlap_matrix), highlight_row = reference_row)
    heatmap_plot_jaccard <- plot_overlap_matrix(overlap_matrix_jaccard, sample_names, highlight_row = reference_row)
    
    # save plots
    collapsed <-
      ifelse(collapse_clones,
             "clones_collapsed",
             "clones_not_collapsed")
    if (sample %in% c("tcr", "umi")) {
      sampled <- paste0(sample, "sampled_to_", sample_to)
    } else {
      sampled <- "clones_not_sampled"
    }
    dir_path <- file.path("./figures/clone_overlap", threshold, collapsed, sampled)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    
    
    save_plot <-
      function(plot, file, width, height = 16) {
        path <-
          file.path(
            "./figures/clone_overlap",
            threshold,
            collapsed,
            sampled,
            paste(collapsed, mouse, chain, file, sep = ".")
          )
        pdf(file = path,
            width = width,
            height = height)
        draw(plot)
        dev.off()
      }
    
    # save raws
    plots <- c(heatmap_plot)
    plot_names <- c(paste0(cell, ".pdf"))
    plot_width <- c(9)
    plot_height <- c(9)
    mapply(save_plot, plots, plot_names, plot_width, plot_height)
    # # save jaccard
    # plots <- c(heatmap_CD4_jaccard, heatmap_CD8_jaccard,
    #            heatmap_CD4_blood_jaccard, heatmap_CD8_blood_jaccard)
    # plot_names <- c("CD4.jaccard.pdf", "CD8.jaccard.pdf",
    #                 "CD4.blood.jaccard.pdf", "CD8.blood.jaccard.pdf")
    # mapply(save_plot, plots, plot_names, plot_width, plot_height)
    
    # plots <- c(heatmap_CD4, heatmap_CD8,
    #            heatmap_CD4_blood, heatmap_CD8_blood)
    # plot_names <- c("CD4.pdf", "CD8.pdf", "CD4.blood.pdf", "CD8.blood.pdf")
    # plot_width <- c(9, 9, 3, 3)
    # plot_height <- c(9, 9, 9, 9)
    # mapply(save_plot, plots, plot_names, plot_width, plot_height)
    # # save jaccard
    # plots <- c(heatmap_CD4_jaccard, heatmap_CD8_jaccard,
    #            heatmap_CD4_blood_jaccard, heatmap_CD8_blood_jaccard)
    # plot_names <- c("CD4.jaccard.pdf", "CD8.jaccard.pdf",
    #                 "CD4.blood.jaccard.pdf", "CD8.blood.jaccard.pdf")
    # mapply(save_plot, plots, plot_names, plot_width, plot_height)
  }


# finds the subsample level for the overlap analysis
# used in plot_clone_overlap()
search_subsample_level <- function(clone_counts) {
  max_remaining <- 0
  best_clone_id <- ""
  for (clone_id in names(clone_counts)) {
    clone_count <- clone_counts[[clone_id]]
    # calculate the remaining TCR count
    remaining <- sum(clone_counts >= clone_count) * clone_count
    if (remaining > max_remaining) {
      max_remaining <- remaining
      best_clone_id <- clone_id
    }
  }
  return(best_clone_id)
}



# finds the subsample UMI level for the overlap analysis
# used in plot_clone_overlap()
search_subsample_level_umi <- function(umi_counts) {
  max_remaining <- 0
  best_clone_id <- ""
  umi_counts_sum <- lapply(umi_counts, sum)
  for (umi_id in names(umi_counts_sum)) {
    umi_count <- umi_counts_sum[[umi_id]]
    # calculate the remaining UMI count
    remaining <- sum(umi_counts_sum >= umi_count) * umi_count
    if (remaining > max_remaining) {
      max_remaining <- remaining
      best_clone_id <- umi_id
    }
  }
  return(best_clone_id)
}




# function to plot
# matrix <- overlap_matrix_cd4_blood
# column_labels <-  sample_name_cd4_blood
plot_overlap_matrix <- function(matrix,
                                column_labels,
                                highlight_row = reference_row,
                                highlight_column = 0,
                                heatmap_title = "",
                                cluster = FALSE,
                                col_rot = 90,
                                font_size = 9,
                                matrix_min = NULL,
                                matrix_max = NULL,
                                order = TRUE
) {
  # Cluster columns
  if (length(column_labels) > 1) {
    if (cluster) {
      matrix[matrix == 100] <- NA
      col_dend <- hclust(as.dist(
        1 - cor(matrix, method = "spearman", use = "pairwise.complete.obs")
      ))
      matrix[is.na(matrix)] <- 100
      ordered_columns <- col_dend$order
      ordered_rows <- col_dend$order
    } else if (order == TRUE) {
      col_dend <- FALSE
      ordered_columns <- order(colnames(matrix))
      ordered_rows <- order(rownames(matrix))
    } else if (order == FALSE) {
      col_dend <- FALSE
      ordered_columns <- seq(1,length(colnames(matrix)))
      ordered_rows <- order(rownames(matrix))
    }
  } 
  
  # highlight row and column
  highlights_row <- rep("black", nrow(matrix))
  highlights_row[highlight_row] <- ifelse(highlight_row == 0, "black", "red")
  highlights_column <- rep("black", nrow(matrix))
  highlights_column[highlight_column] <- ifelse(highlight_column == 0, "black", "red")
  
  
  heatmap_plot <- Heatmap(
    matrix,
    cluster_columns = col_dend,
    # cluster_columns = cluster_columns,
    column_labels = column_labels,
    column_order = ordered_columns,
    column_names_side = "top",
    column_dend_side = "bottom",
    column_names_rot = col_rot,
    column_dend_height = unit(2, "cm"),
    row_order = ordered_rows,
    row_names_side = "left",
    row_names_gp = gpar(col = highlights_row, fontsize = 8),
    column_names_gp = gpar(col = highlights_column, fontsize = 8),
    show_heatmap_legend = FALSE,
    col = circlize::colorRamp2(c(0,50,100), c("blue", "white", "red")),
    column_title = heatmap_title,
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (is.matrix(matrix)) {
        val <- matrix[i, j]
        val_min <- matrix_min[i, j]
        val_max <- matrix_max[i, j]
        grid.text(paste0(sprintf("%.1f\n(%.1f - %.1f)", val, val_min, val_max)), x, y, gp = gpar(fontsize = font_size))
      }
      else {
        val <- matrix[i]
        val_min <- matrix_min[i]
        val_max <- matrix_max[i]
        grid.text(paste0(sprintf("%.2f\n(%.2f - %.2f)", val, val_min, val_max)), x, y, gp = gpar(fontsize = font_size))
      }
    }
  )
  return(heatmap_plot)
}



# Calculate predefined groups' overlap
plot_overlap_sampled_and_reference <-
  function(group, group_name, n_bootstrap_samples = 1000) {
    # get the ids of samples to be compared
    id_sam <- metadata %>% filter(
      mouse == .env$mouse,
      grepl(paste0(" ", .env$cell, " "), sample.id, ignore.case = TRUE),
      organ %in% .env$group
    ) %>% pull(gs.id)
    
    # get the ids of reference samples (+2 columns)
    id_ref <- metadata %>% filter(mouse == .env$mouse,
                                  grepl(
                                    paste0("spleen ", cell, " cd69(\\+|\\- 80k)$"),
                                    sample.id,
                                    ignore.case = TRUE
                                  )) %>% pull(gs.id)
    if (chain == "TRA") {
      reference_dfs <- clone_dfs_tra[id_sam]
      reference_columns_dfs <- clone_dfs_tra[id_ref]
    } else {
      reference_dfs <- clone_dfs_trb[id_sam]
      reference_columns_dfs <- clone_dfs_trb[id_ref]
    }
    
    # group by nCDR3 seq
    summarize_clones <- function(df) {
      df %>%
        group_by(nSeqCDR3) %>%
        summarise(uniqueMoleculeCount = sum(uniqueMoleculeCount),
                  .groups = "drop")
    }
    
    reference_dfs <- lapply(reference_dfs, summarize_clones)
    reference_columns_dfs <- lapply(reference_columns_dfs, summarize_clones)
    
    # count # of clonotypes
    tcr_counts <- lapply(reference_dfs, nrow)
    tcr_counts_reference_columns <- lapply(reference_columns_dfs, nrow)
    
    # find which sample has the lowest # of clonotypes
    sampled_to_column <-  which.min(unlist(tcr_counts[id_sam]))
    sampled_to_reference_column <-  which.min(unlist(tcr_counts_reference_columns[id_ref]))
    # find the sample size
    sample_to <- unlist(tcr_counts[id_sam])[[sampled_to_column]]
    # find which REFERENCE sample has the lowest # of clonotypes
    # in the 2nd heatmap
    sample_to_reference_column <- unlist(tcr_counts_reference_columns[id_ref])[[sampled_to_reference_column]]
    # sample data
    options(future.globals.maxSize= 1500*1024^2)
    Sys.setenv(OPENBLAS_NUM_THREADS = "1")
    n_of_cores <- parallel::detectCores()
    plan(multisession, workers = n_of_cores - 4)
    cat("\n\t Sampling:", cell, chain)
    cat("\t Processing samples...")
    # Future for samples
    future_samples <- future({
      lapply(1:n_bootstrap_samples, function(i) {
        mapply(
          subsample_abundance,
          reference_dfs[id_sam],
          sample_to,
          target_type = "clone",
          SIMPLIFY = FALSE,
          seed = i
        )
      })
    })
    cat("\t Processing references...")
    # Future for references
    future_references <- future({
      lapply(1:n_bootstrap_samples, function(i) {
        mapply(
          subsample_abundance,
          reference_columns_dfs[id_ref],
          sample_to_reference_column,
          target_type = "clone",
          SIMPLIFY = FALSE,
          seed = i
        )
      })
    })
    cat("\t Done.")
    
    # Resolve both futures (wait for completion)
    reference_dfs_sampled <- value(future_samples)
    reference_columns_dfs_sampled <- value(future_references)
    plan(sequential)
    # count umis
    umi_counts <- lapply(reference_dfs, function(x) {
      sum(x$uniqueMoleculeCount)
    })
    umi_counts_reference_columns <- lapply(reference_columns_dfs, function(x) {
      sum(x$uniqueMoleculeCount)
    })
    # count umis and tcr after sampling
    umi_counts_sampled <- sapply(1:n_bootstrap_samples, function(i) {
      sapply(reference_dfs_sampled[[i]], sum)
    })
    umi_counts_sampled_median <- apply(umi_counts_sampled, 1, median)
    umi_counts_sampled_min <- apply(umi_counts_sampled, 1, quantile, probs = 0.025)
    umi_counts_sampled_max <- apply(umi_counts_sampled, 1, quantile, probs = 0.975)
    # count umis and tcr after sampling for reference columns
    umi_counts_sampled_reference_columns <- sapply(1:n_bootstrap_samples, function(i) {
      sapply(reference_columns_dfs_sampled[[i]], sum)
    })
    umi_counts_sampled_median_reference_columns <-
      apply(umi_counts_sampled_reference_columns, 1, median)
    umi_counts_sampled_min_reference_columns <-
      apply(umi_counts_sampled_reference_columns, 1, quantile, probs = 0.025)
    umi_counts_sampled_max_reference_columns <-
      apply(umi_counts_sampled_reference_columns, 1, quantile, probs = 0.975)
    
    tcr_counts_sampled <- lapply(1:n_bootstrap_samples, function(i) {
      lapply(reference_dfs_sampled[[i]], length)
    })
    tcr_counts_sampled_reference_columns <- lapply(1:n_bootstrap_samples, function(i) {
      lapply(reference_columns_dfs_sampled[[i]], length)
    })
    
    # count rows & create matrix with it
    reference_dfs_len <- length(reference_dfs)
    overlap_matrix_ref <-
      replicate(
        n_bootstrap_samples,
        matrix(nrow = reference_dfs_len, ncol = reference_dfs_len),
        simplify = FALSE
      )
    
    # count rows & create matrix with it for the reference cols
    reference_columns_dfs_len <- length(reference_columns_dfs)
    overlap_matrix_ref_columns <-
      replicate(
        n_bootstrap_samples,
        matrix(nrow = reference_dfs_len, ncol = reference_columns_dfs_len + 2), # overlap between references: +2
        simplify = FALSE
      )
    
    
    # fill in the matrix
    ## iterate rows
    for (i in 1:reference_dfs_len) {
      vectors_1 <- lapply(1:n_bootstrap_samples, function(x) {
        names(reference_dfs_sampled[[x]][[i]])
      })
      ns_i <- lapply(vectors_1, length)
      
      ## iterate columns for main heatmap
      for (j in 1:reference_dfs_len) {
        # Find overlaps
        vectors_2 <- lapply(1:n_bootstrap_samples, function(x) {
          names(reference_dfs_sampled[[x]][[j]])
        })
        # ns_j <- lapply(vectors_2, length)
        # calculate overlap
        for (x in 1:n_bootstrap_samples) {
          overlap_cur <- length(intersect(vectors_1[[x]], vectors_2[[x]]))
          overlap_matrix_ref[[x]][i, j] <-
            (overlap_cur / ns_i[[x]]) * 100
        }
      }
      
      ## iterate columns for reference column heatmap
      for (jcol in 1:reference_columns_dfs_len) {
        # Find overlaps
        vectors_2_cols <- lapply(1:n_bootstrap_samples, function(x) {
          names(reference_columns_dfs_sampled[[x]][[jcol]])
        })
        # calculate overlap
        for (x in 1:n_bootstrap_samples) {
          overlap_cur <- length(intersect(vectors_1[[x]], vectors_2_cols[[x]]))
          overlap_matrix_ref_columns[[x]][i, jcol] <-
            (overlap_cur / ns_i[[x]]) * 100
        }
      }
      
      ## calculate sample vs reference1 and reference 2 and reference1 vs reference2 overlaps
      # Find overlaps
      vectors_2_cols <- lapply(1:n_bootstrap_samples, function(x) {
        intersect(
          names(reference_columns_dfs_sampled[[x]][[1]]),
          names(reference_columns_dfs_sampled[[x]][[2]]))
      })
      # calculate overlap
      for (x in 1:n_bootstrap_samples) {
        overlap_cur <- length(intersect(vectors_1[[x]], vectors_2_cols[[x]]))
        overlap_matrix_ref_columns[[x]][i, 3] <-
          (overlap_cur / ns_i[[x]]) * 100
        # ref overlaps
        overlap_cur <- length(vectors_2_cols[[x]])
        overlap_matrix_ref_columns[[x]][i, 4] <-
          (overlap_cur / sample_to_reference_column) * 100
        # reference multiplication
        # multip_cur <- overlap_matrix_ref_columns[[x]][i, 1] * overlap_matrix_ref_columns[[x]][i, 2] / 100
        # overlap_matrix_ref_columns[[x]][i, 4] <- multip_cur
      }
    }
    
    
    # calculate overlap matrix metrics
    overlap_matrix_simple <- matrix(nrow = reference_dfs_len, ncol = reference_dfs_len)
    overlap_matrix_simple_min <- matrix(nrow = reference_dfs_len, ncol = reference_dfs_len)
    overlap_matrix_simple_max <- matrix(nrow = reference_dfs_len, ncol = reference_dfs_len)
    for (j in 1:reference_dfs_len) {
      # collect columns
      for (i in 1:reference_dfs_len) {
        vals <- sapply(overlap_matrix_ref, function(m)
          m[i, j])
        overlap_matrix_simple[i, j] <- median(vals, na.rm = TRUE)
        overlap_matrix_simple_min[i, j] <- quantile(vals, 0.025, na.rm = TRUE)
        overlap_matrix_simple_max[i, j] <- quantile(vals, 0.975, na.rm = TRUE)
      }
    }
    # calculate reference columns overlap matrix metrics
    overlap_matrix_simple_reference_columns <- matrix(nrow = reference_dfs_len, ncol = reference_columns_dfs_len+2)
    overlap_matrix_simple_min_reference_columns <- overlap_matrix_simple_reference_columns
    overlap_matrix_simple_max_reference_columns <- overlap_matrix_simple_reference_columns
    for (j in 1:(reference_columns_dfs_len+2)) {
      # collect columns
      for (i in 1:reference_dfs_len) {
        vals <- sapply(overlap_matrix_ref_columns, function(m)
          m[i, j])
        overlap_matrix_simple_reference_columns[i, j] <- median(vals, na.rm = TRUE)
        overlap_matrix_simple_min_reference_columns[i, j] <- quantile(vals, 0.025, na.rm = TRUE)
        overlap_matrix_simple_max_reference_columns[i, j] <- quantile(vals, 0.975, na.rm = TRUE)
      }
    }
    
    overlap_matrix_simple <- as.matrix(overlap_matrix_simple)
    col_names <- metadata %>% filter(gs.id %in% .env$id_sam) %>% pull(sample.id)
    col_names_reference_columns <-
      metadata %>% filter(gs.id %in% .env$id_ref) %>% pull(sample.id)
    row_names <- metadata %>% filter(gs.id %in% .env$id_sam) %>%
      pull(sample.id)
    colnames(overlap_matrix_simple) <- paste(
      str_to_upper(substring(col_names, 5)),
      "\numi = ",
      paste0(
        umi_counts,
        " -> ",
        umi_counts_sampled_median,
        " (",
        umi_counts_sampled_min,
        "-",
        umi_counts_sampled_max,
        ") "
      ),
      "\ntcr = ",
      tcr_counts,
      " -> ",
      sample_to
    )
    rownames(overlap_matrix_simple) <- colnames(overlap_matrix_simple)
    rownames(overlap_matrix_simple_reference_columns) <- rownames(overlap_matrix_simple)
    colnames(overlap_matrix_simple_reference_columns) <- 
      c(
        paste(
          str_to_upper(substring(col_names_reference_columns, 5)),
          "\numi = ",
          paste0(
            umi_counts_reference_columns,
            " -> ",
            umi_counts_sampled_median_reference_columns,
            " (",
            umi_counts_sampled_min_reference_columns,
            "-",
            umi_counts_sampled_max_reference_columns,
            ") "
          ),
          "\ntcr = ",
          tcr_counts_reference_columns,
          " -> ",
          sample_to_reference_column
        ),
        "Overlap with both reference repertoires",
        "Spleen CD69+ and CD69- Clonotype Overlap"
      )
    
    
    
    title <- paste(mouse, cell, chain)
    p_ref <-
      plot_overlap_matrix(
        overlap_matrix_simple,
        column_labels = colnames(overlap_matrix_simple),
        highlight_row = sampled_to_column,
        highlight_column = sampled_to_column,
        heatmap_title = title,
        cluster = FALSE,
        col_rot = 90,
        font_size = 11,
        matrix_min = overlap_matrix_simple_min,
        matrix_max = overlap_matrix_simple_max
      )
    print(p_ref)
    # reference columns
    p_ref_columns <-
      plot_overlap_matrix(
        overlap_matrix_simple_reference_columns,
        column_labels = colnames(overlap_matrix_simple_reference_columns),
        highlight_row = 0,
        highlight_column = sampled_to_reference_column,
        heatmap_title = NULL,
        cluster = FALSE,
        col_rot = 90,
        font_size = 11,
        matrix_min = overlap_matrix_simple_min_reference_columns,
        matrix_max = overlap_matrix_simple_max_reference_columns,
        order = FALSE
      )
    print(p_ref + p_ref_columns)
    # save to file
    dir_path <- file.path("./figures/clone_overlap/sample_to_min/", threshold, group_name)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    # save plot to made files
    pdf(
      file = file.path(dir_path, paste(mouse, cell, chain, "pdf", sep = ".")),
      width = 13,
      height = 6.5
    )
    print(p_ref + p_ref_columns)
    dev.off()
  }




# Calculate CD69- and CD69+ overlap
plot_tissue_cd69_pairs <-
  function(mouse, tissue, cell, 
           chain,
           n_bootstrap_samples = 1,
           sample_type = "umi") { # umi or clone!!!
    # get the ids of samples to be compared
    meta_x_y <- metadata %>% filter(mouse == .env$mouse,
                                    organ == .env$organ,
                                    grepl(.env$cell, Cell, ignore.case = TRUE),
                                    cell.numbers != 20000) #%>% pull(gs.id)
    if ((nrow(meta_x_y) == 0) || (nrow(meta_x_y) == 1)) return()
    # x vs y: get IDs
    x_dataframe_id <- meta_x_y %>% 
      filter(grepl("CD69\\-", Cell, ignore.case = TRUE)) %>% 
      pull(gs.id)# cd69-
    
    y_dataframe_id <- meta_x_y %>% 
      filter(grepl("CD69\\+", Cell, ignore.case = TRUE)) %>% 
      pull(gs.id)# cd69+
    
    
    if (chain == "TRA") {
      x_y_dataframes <- list(x = clone_dfs_tra[[x_dataframe_id]],
                             y = clone_dfs_tra[[y_dataframe_id]])
    } else {
      x_y_dataframes <- list(x = clone_dfs_trb[[x_dataframe_id]],
                             y = clone_dfs_trb[[y_dataframe_id]])
    }
    
    # group by nCDR3 seq
    summarize_clones <- function(df) {
      df %>%
        group_by(nSeqCDR3) %>%
        summarise(uniqueMoleculeCount = sum(uniqueMoleculeCount),
                  .groups = "drop")
    }
    
    x_y_data <- lapply(x_y_dataframes, summarize_clones)
    
    # count # of clonotypes
    tcr_counts <- lapply(x_y_data, nrow)
    umi_counts <- lapply(x_y_data, function(x) sum(x$uniqueMoleculeCount))
    # find which sample has the lowest # of clonotypes
    sampled_to <-  ifelse(sample_type == "clone", 
                          which.min(unlist(tcr_counts)),
                          which.min(unlist(umi_counts)))
    # find the sample size
    sample_to <-  ifelse(sample_type == "clone", 
                         unlist(tcr_counts)[[sampled_to]],
                         unlist(umi_counts)[[sampled_to]])
    # sample data
    #options(future.globals.maxSize= 1500*1024^2)
    #Sys.setenv(OPENBLAS_NUM_THREADS = "1")
    #n_of_cores <- parallel::detectCores()
    #plan(multisession, workers = n_of_cores - 4)
    cat("\n\t Sampling:", cell, chain)
    # Future for samples
    future_samples <- future({
      lapply(1:n_bootstrap_samples, function(i) {
        mapply(
          subsample_abundance,
          x_y_data,
          sample_to,
          target_type = sample_type,
          SIMPLIFY = FALSE,
          seed = i
        )
      })
    })
    
    # Resolve both futures (wait for completion)
    x_y_data_sampled <- value(future_samples)
    #plan(sequential)
    
    # Sort by UMIs after sampling (keep structure as list of lists)
    x_y_data_sampled <- lapply(1:n_bootstrap_samples, function(i) {
      lapply(x_y_data_sampled[[i]], sort, decreasing = TRUE)
    })
    
    # Count UMIs and TCR sums after sampling (list of lists -> numeric vectors)
    umi_counts_sampled <- lapply(1:n_bootstrap_samples, function(i) {
      sapply(x_y_data_sampled[[i]], sum)
    })
    
    # Combine list of named numeric vectors into a matrix (each column a bootstrap sample)
    umi_counts_mat <- do.call(cbind, umi_counts_sampled)
    
    # Calculate statistics across bootstrap samples for each element:
    umi_counts_sampled_median <- apply(umi_counts_mat, 1, median)
    umi_counts_sampled_min <- apply(umi_counts_mat, 1, quantile, probs = 0.025)
    umi_counts_sampled_max <- apply(umi_counts_mat, 1, quantile, probs = 0.975)
    
    
    tcr_counts_sampled <- lapply(1:n_bootstrap_samples, function(i) {
      lapply(x_y_data_sampled[[i]], length)
    })
    
    # Make plot data
    get_values <- function(lst, keys) {
      vals <- sapply(keys, function(k) if (k %in% names(lst)) lst[[k]] else NA)
      return(vals)
    }
    
    sequences <- union(names(x_y_data_sampled[[1]]$x),names(x_y_data_sampled[[1]]$y))
    sequences_overlap <- intersect(names(x_y_data_sampled[[1]]$x),names(x_y_data_sampled[[1]]$y))
    count_sequences_overlap <- length(sequences_overlap)
    
    plot_data <-
      data.frame(
        Sequence = sequences,
        X = get_values(x_y_data_sampled[[1]][["x"]], sequences),
        Y = get_values(x_y_data_sampled[[1]][["y"]], sequences),
        row.names = NULL
      ) %>%
      mutate(
        x_p = X / sum(X, na.rm = TRUE),
        y_p = Y / sum(Y, na.rm = TRUE),
        Group_scatter = case_when((is.na(X) & !is.na(Y) ~ "unique"),
                                  (!is.na(X) & is.na(Y) ~ "unique"),
                                  (!is.na(X) & !is.na(Y) ~ "shared")
        ),
        Group_histogram = case_when((!is.na(X) & !is.na(Y) ~ "overlap"),
                                    (is.na(X) & !is.na(Y) ~ "y"),
                                    (!is.na(X) & is.na(Y) ~ "x"))
      )
    
    plot_data <- plot_data %>% replace(is.na(.), 0)
    
    # plot_data[c("x_p", "y_p")][plot_data[c("x_p", "y_p")]==0] = 0.000001
    # plot the main scatterplot
    colors_grop <- c("all" = "black", "overlap" = "darkorange", "x" = "#619CFF", "y" = "#619CFF",
                     "unique" = "#619CFF", "shared" = "darkorange")
    x_text <- plot_data%>%top_n(10,x_p)%>%slice_head(n = 10)
    y_text <- plot_data%>%top_n(10,y_p)%>%slice_head(n = 10)
    n_x_text <- nrow(x_text)
    n_y_text <- nrow(y_text)
    p1 <- ggplot(plot_data, aes(x = x_p, y = y_p, colour = Group_scatter)) +
      geom_point(size = 2, shape = 21, alpha = 0.75) + 
      geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") +
      scale_x_log10(breaks = 10^seq(-6, 0, by = 1),
                    limits = c(1e-6, 1),
                    labels = scales::trans_format("log10", scales::math_format(10^.x)),
                    # oob = scales::squish_infinite
      ) +
      scale_y_log10(breaks = 10^seq(-6, 0, by = 1),
                    limits = c(1e-6, 1),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      theme(legend.position = "none") +
      labs(x = paste0(cell, " CD69-", "\n(",organ,")"),
           y =  paste0(cell, " CD69+", "\n(",organ,")")) +
      my_theme_p1 + 
      coord_cartesian(clip = 'off') + 
      scale_color_manual(values=colors_grop) + 
      geom_text(data = y_text, 
                aes(x = 1e-6, y=10^seq(0, -1.5, length.out = n_y_text), color = Group_scatter,label = Sequence), 
                hjust = 0, vjust = 1, check_overlap = 0, size = 3) + 
      geom_point(data = y_text, 
                aes(x = x_p, y=y_p, color = Group_scatter)) + 
      geom_text(data = x_text, 
                aes(x = 1e-2, y=10^seq(-4, -5.5, length.out = n_x_text), color = Group_scatter,label = Sequence), 
                hjust = 0, vjust = 1, check_overlap = 0, size = 3) + 
      geom_point(data = x_text, 
                 aes(x = x_p, y=y_p, color = Group_scatter))
    
    g <- ggplotGrob(p1)
    legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
    
    p1 <- p1 + theme(legend.position = "none")
    
    # x distribution
    breaks <- seq(-6, 0, by = 1)
    colors_grop <- c("all" = "black", "overlap" = "orange", "x" = "#619CFF", "y" = "#619CFF")
    
    plot_data_x <- plot_data%>%filter(Group_histogram!="y")
    px <- ggplot() + 
      geom_freqpoly(data = plot_data_x,
                   aes(x = log10(x_p)), alpha=1, linewidth = 1) +
      geom_freqpoly(data = plot_data_x,
                   aes(x = log10(x_p),color=Group_histogram), alpha=0.8, linewidth = 1,
                   ) +
      scale_y_log10(limits = c(1e0, NA),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      my_theme() +
      scale_x_continuous(breaks = seq(-6,0,by = 1),
                         limits = c(-6,0),
                         labels = scales::math_format()) +
      theme(axis.text=element_text(size=6)) + 
      ylab("Count") +
      scale_fill_manual(values=colors_grop)
    
    # y distribution
    plot_data_y <- plot_data%>%filter(Group_histogram!="x")
    py <- ggplot() + 
      geom_freqpoly(data = plot_data_y,
                        aes(x = log10(y_p)), alpha=1, linewidth = 1) +
      geom_freqpoly(data = plot_data_y,
                    aes(x = log10(y_p),color=Group_histogram), alpha=0.8, linewidth = 1,
      ) +
      scale_y_log10(limits = c(1, NA), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      scale_x_continuous(breaks = seq(-6,0,by = 1),
                         limits = c(-6,0),
                         labels = scales::math_format()) + coord_flip() +
      my_theme_y + 
      theme(axis.text=element_text(size=6)) +
      ylab("Count") +
      scale_fill_manual(values=colors_grop)
    
    # Remove texts of some axes: x from x and y from y plots
    px <- px + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) + theme(legend.position = "none")
    py <- py + theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) + theme(legend.position = "none")
    
    # stats
    cell_x_y <- meta_x_y$cell.numbers
    cell_text <- paste("Cell number\n", paste(c("CD69-", "CD69+"), cell_x_y, sep = " = ", collapse = "\n"))
    tcr_count_text <- paste("TCR\n", paste(c("CD69-", "CD69+"), tcr_counts, sep = " = ", collapse = "\n"))
    umi_count_text <- paste("UMI\n", paste(c("CD69-", "CD69+"), umi_counts, sep = " = ", collapse = "\n"))
    
    sampled_umi_count_text <- paste("Sampled UMIs\n", paste(c("CD69-", "CD69+"), umi_counts_sampled[[1]], sep = " = ", collapse = "\n"))
    sampled_tcr_count_text <- paste("Sampled TCRs\n", paste(c("CD69-", "CD69+"), tcr_counts_sampled[[1]], sep = " = ", collapse = "\n"))
    
    
    overlap_text <- paste(c("Overlap"), length(sequences_overlap), sep = " = ", collapse = "\n")
    df_text <- plot_data %>% group_by(Group_histogram) %>% summarise(n = n()) %>%
      mutate(Group_histogram = factor(Group_histogram, levels = c("x", "y", "overlap"))) %>%
      arrange(Group_histogram)
    overlapDetail_text <- paste("Overlap info:\n", paste(c("CD69-", "CD69+", "overlap"), df_text$n, sep = " = ", collapse = "\n"))
    
    allText_left <- paste(cell_text,tcr_count_text, umi_count_text,tcr_count_text, sep="\n")
    allText_right <- paste(sampled_umi_count_text,sampled_tcr_count_text, overlap_text, overlapDetail_text, sep="\n")
    plot_text <-
      ggplot() +
      annotate(
        "text",
        x = 0.1,
        y = 0.9,
        size = 4,
        label = allText_left,
        hjust = 0,
        vjust = 1
      ) +
      annotate(
        "text",
        x = 0.5,
        y = 0.9,
        size = 4,
        label = allText_right,
        hjust = 0,
        vjust = 1
      ) +
      theme_void() + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
    
    design <-
      "AD
       BC"
    final_plot <- (px) + (p1) + (py) + plot_text + plot_layout(design = design, widths = c(1, 0.3),
                                                               heights = c(0.3, 1))
    # p1+inset_element(legend, left = 0.6, bottom = 0.6, right = 1, top = 1)
    # save to fil# save to fil# save to file
    dir_path <- file.path("./figures/clone_overlap_CD69/", threshold, sample_type)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    
    # save plot to made files
    ggsave(
      filename = file.path(dir_path, paste(mouse, organ, cell, chain, "pdf", sep = ".")),
      plot = final_plot,
      width = 15,
      height = 15, 
      device = "pdf"
    )
  }


# Calculate overlap of top 100s 
plot_overlap_100_sequence_overlap <-
  function(group, group_name, n_bootstrap_samples = 1) {
    
    # get the ids of samples to be compared
    id_sam <- metadata %>% filter(
      mouse == .env$mouse,
      grepl(paste0(" ", .env$cell, " "), sample.id, ignore.case = TRUE),
      organ %in% .env$group
    ) %>% pull(gs.id)
    
    # get the ids of reference samples (+2 columns)
    id_ref <- metadata %>% filter(mouse == .env$mouse,
                                  grepl(
                                    paste0("spleen ", cell, " cd69(\\+|\\- 80k)$"),
                                    sample.id,
                                    ignore.case = TRUE
                                  )) %>% pull(gs.id)
    
    if (chain == "TRA") {
      reference_dfs <- clone_dfs_tra[id_sam]
      reference_columns_dfs <- clone_dfs_tra[id_ref]
    } else {
      reference_dfs <- clone_dfs_trb[id_sam]
      reference_columns_dfs <- clone_dfs_trb[id_ref]
    }
    
    # group by nCDR3 seq
    summarize_clones <- function(df) {
      df %>%
        group_by(nSeqCDR3) %>%
        summarise(uniqueMoleculeCount = sum(uniqueMoleculeCount),
                  .groups = "drop") %>%
        slice_max(uniqueMoleculeCount,n=100, with_ties = FALSE)
    }
    
    reference_dfs <- lapply(reference_dfs, summarize_clones)
    reference_columns_dfs <- lapply(reference_columns_dfs, summarize_clones)
    
    # count # of clonotypes
    tcr_counts <- lapply(reference_dfs, nrow)
    tcr_counts_reference_columns <- lapply(reference_columns_dfs, nrow)
    
    # count umis
    umi_counts <- lapply(reference_dfs, function(x) {
      sum(x$uniqueMoleculeCount)
    })
    umi_counts_reference_columns <- lapply(reference_columns_dfs, function(x) {
      sum(x$uniqueMoleculeCount)
    })
    
    # count rows & create main matrix with it
    reference_dfs_len <- length(reference_dfs)
    overlap_matrix_ref <-
      replicate(
        n_bootstrap_samples,
        matrix(nrow = reference_dfs_len, ncol = reference_dfs_len),
        simplify = FALSE
      )
    
    # Create top 10 sequence storage matrix
    overlap_matrix_sequence <- matrix(list(), nrow = reference_dfs_len, ncol = reference_dfs_len)
 
    # count rows & create matrix with it for the reference cols
    reference_columns_dfs_len <- length(reference_columns_dfs)
    overlap_matrix_ref_columns <-
      replicate(
        n_bootstrap_samples,
        matrix(nrow = reference_dfs_len, ncol = reference_columns_dfs_len + 2), # overlap between references: +2
        simplify = FALSE
      )
    
    
    
    # fill in the matrix
    ## iterate rows
    for (i in 1:reference_dfs_len) {
      
      sequences_1 <- reference_dfs[[i]][["nSeqCDR3"]]
      ns_i <- length(sequences_1)
      
      ## iterate columns for main heatmap
      for (j in 1:reference_dfs_len) {
        # Find overlaps
        sequences_2 <- reference_dfs[[j]][["nSeqCDR3"]]
        # combine sequence and abundances
        combined_df <- 
          reference_dfs[[i]] %>% inner_join(reference_dfs[[j]], by = "nSeqCDR3") %>%
          mutate(rank = rank(uniqueMoleculeCount.x)+rank(uniqueMoleculeCount.y))
        
        cosine_current <- 
          cosine(combined_df$uniqueMoleculeCount.x, combined_df$uniqueMoleculeCount.y)
        
        ns_j <- length(sequences_2)
        
        # calculate overlap and store
        overlap_cur <- length(intersect(sequences_1, sequences_2))
        
        overlap_count_cur <- length(overlap_cur)
        overlap_matrix_ref[[1]][i, j] <-
          (overlap_cur / ns_i[[1]]) * 100
        
        # get top 10 sequences and store 
        overlap_matrix_sequence[[i,j]] <- combined_df %>% slice_max(rank,n=10,with_ties = FALSE) %>% pull(nSeqCDR3)
        }
      
      ## iterate columns for reference column heatmap
      for (jcol in 1:reference_columns_dfs_len) {
        # Find overlaps
        combined_df <- 
          reference_dfs[[i]] %>% inner_join(reference_columns_dfs[[jcol]], by = "nSeqCDR3") %>%
          mutate(rank = rank(uniqueMoleculeCount.x)+rank(uniqueMoleculeCount.y))

        cosine_current <- 
          cosine(combined_df$uniqueMoleculeCount.x, combined_df$uniqueMoleculeCount.y)
        
        # calculate overlap of 100 top clones
        overlap_cur <- nrow(combined_df)
        overlap_matrix_ref_columns[[1]][i, jcol] <-
          (overlap_cur / ns_i[[1]]) * 100
      }
      
      ## calculate sample vs reference1 and reference 2 and reference1 vs reference2 overlaps
      # Find overlaps
      combined_df <- 
        reference_dfs[[i]] %>% inner_join(reference_columns_dfs[[1]], by = "nSeqCDR3") %>%
        inner_join(reference_columns_dfs[[2]], by = "nSeqCDR3") %>%
        mutate(rank = rank(uniqueMoleculeCount.x)+rank(uniqueMoleculeCount.y)+rank(uniqueMoleculeCount))
      
      combined_df_refs <- 
        reference_columns_dfs[[1]] %>%
        inner_join(reference_columns_dfs[[2]], by = "nSeqCDR3") %>%
        mutate(rank = rank(uniqueMoleculeCount.x)+rank(uniqueMoleculeCount.y))
      
      
      # calculate overlap
      overlap_matrix_ref_columns[[1]][i, 3] <-
          (nrow(combined_df) / ns_i[[1]]) * 100
      # ref overlaps
      overlap_matrix_ref_columns[[1]][i, 4] <-
        (nrow(combined_df_refs) / 100) * 100
      }
    
    
    # calculate overlap matrix metrics
    overlap_matrix_simple <- matrix(nrow = reference_dfs_len, ncol = reference_dfs_len)
    overlap_matrix_simple_min <- matrix(nrow = reference_dfs_len, ncol = reference_dfs_len)
    overlap_matrix_simple_max <- matrix(nrow = reference_dfs_len, ncol = reference_dfs_len)
    for (j in 1:reference_dfs_len) {
      # collect columns
      for (i in 1:reference_dfs_len) {
        vals <- sapply(overlap_matrix_ref, function(m)
          m[i, j])
        overlap_matrix_simple[i, j] <- median(vals, na.rm = TRUE)
        overlap_matrix_simple_min[i, j] <- quantile(vals, 0.025, na.rm = TRUE)
        overlap_matrix_simple_max[i, j] <- quantile(vals, 0.975, na.rm = TRUE)
      }
    }
    
    # calculate reference columns overlap matrix metrics
    overlap_matrix_simple_reference_columns <- matrix(nrow = reference_dfs_len, ncol = reference_columns_dfs_len+2)
    overlap_matrix_simple_min_reference_columns <- overlap_matrix_simple_reference_columns
    overlap_matrix_simple_max_reference_columns <- overlap_matrix_simple_reference_columns
    
    for (j in 1:(reference_columns_dfs_len+2)) {
      # collect columns
      for (i in 1:reference_dfs_len) {
        vals <- sapply(overlap_matrix_ref_columns, function(m)
          m[i, j])
        overlap_matrix_simple_reference_columns[i, j] <- median(vals, na.rm = TRUE)
        overlap_matrix_simple_min_reference_columns[i, j] <- quantile(vals, 0.025, na.rm = TRUE)
        overlap_matrix_simple_max_reference_columns[i, j] <- quantile(vals, 0.975, na.rm = TRUE)
      }
    }
    
    overlap_matrix_simple <- as.matrix(overlap_matrix_simple)
    col_names <- metadata %>% filter(gs.id %in% .env$id_sam) %>% pull(sample.id)
    col_names_reference_columns <-
      metadata %>% filter(gs.id %in% .env$id_ref) %>% pull(sample.id)
    row_names <- metadata %>% filter(gs.id %in% .env$id_sam) %>%
      pull(sample.id)
    colnames(overlap_matrix_simple) <- paste(
      str_to_upper(substring(col_names, 5)),
      "\numi = ",
      paste0(
        umi_counts
      ),
      "\ntcr = ",
      tcr_counts
    )
    
    rownames(overlap_matrix_simple) <- colnames(overlap_matrix_simple)
    rownames(overlap_matrix_simple_reference_columns) <- rownames(overlap_matrix_simple)
    colnames(overlap_matrix_simple_reference_columns) <- 
      c(
        paste(
          str_to_upper(substring(col_names_reference_columns, 5)),
          "\numi = ",
          paste0(
            umi_counts_reference_columns
          ),
          "\ntcr = ",
          tcr_counts_reference_columns
        ),
        "Overlap with both reference repertoires",
        "Spleen CD69+ and CD69- Clonotype Overlap"
      )
    
    
    
    title <- paste(mouse, cell, chain)
    p_ref <-
      plot_overlap_matrix(
        overlap_matrix_simple,
        column_labels = colnames(overlap_matrix_simple),
        highlight_row = 0,
        highlight_column = 0,
        heatmap_title = title,
        cluster = FALSE,
        col_rot = 90,
        font_size = 11,
        matrix_min = overlap_matrix_simple_min,
        matrix_max = overlap_matrix_simple_max
      )
    print(p_ref)
    # reference columns
    p_ref_columns <-
      plot_overlap_matrix(
        overlap_matrix_simple_reference_columns,
        column_labels = colnames(overlap_matrix_simple_reference_columns),
        highlight_row = 0,
        highlight_column = 0,
        heatmap_title = NULL,
        cluster = FALSE,
        col_rot = 90,
        font_size = 11,
        matrix_min = overlap_matrix_simple_min_reference_columns,
        matrix_max = overlap_matrix_simple_max_reference_columns,
        order = FALSE
      )
    print(p_ref + p_ref_columns)

    # save to file
    dir_path <- file.path("./figures/clone_overlap/top100/", threshold, group_name)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    # save plot to made files
    pdf(
      file = file.path(dir_path, paste(mouse, cell, chain, "pdf", sep = ".")),
      width = 13,
      height = 6.5
    )
    print(p_ref + p_ref_columns)
    dev.off()
  }


# Calculate overlap of top 100 to (min of samples) sampled 
plot_overlap_100_and_sampled <-
  function(group, group_name, n_bootstrap_samples = 1000) {
    
    # get the ids of samples to be compared
    id_sam <- metadata %>% filter(
      mouse == .env$mouse,
      grepl(paste0(" ", .env$cell, " "), sample.id, ignore.case = TRUE),
      organ %in% .env$group
    ) %>% pull(gs.id)
    
    # get the ids of reference samples (+2 columns)
    id_ref <- metadata %>% filter(mouse == .env$mouse,
                                  grepl(
                                    paste0("spleen ", cell, " cd69(\\+|\\- 80k)$"),
                                    sample.id,
                                    ignore.case = TRUE
                                  )) %>% pull(gs.id)
    if (chain == "TRA") {
      reference_dfs <- clone_dfs_tra[id_sam]
      reference_columns_dfs <- clone_dfs_tra[id_ref]
    } else {
      reference_dfs <- clone_dfs_trb[id_sam]
      reference_columns_dfs <- clone_dfs_trb[id_ref]
    }
    
    # group by nCDR3 seq
    summarize_clones <- function(df) {
      df %>%
        group_by(nSeqCDR3) %>%
        summarise(uniqueMoleculeCount = sum(uniqueMoleculeCount),
                  .groups = "drop")
    }
    
    reference_dfs <- lapply(reference_dfs, summarize_clones)
    reference_columns_dfs <- lapply(reference_columns_dfs, summarize_clones)
    
    # count # of clonotypes
    tcr_counts <- lapply(reference_dfs, nrow)
    tcr_counts_reference_columns <- lapply(reference_columns_dfs, nrow)
    
    # find which sample has the lowest # of clonotypes
    sampled_to_column <-  which.min(unlist(tcr_counts[id_sam]))
    sampled_to_reference_column <-  which.min(unlist(tcr_counts_reference_columns[id_ref]))
    # find the sample size
    sample_to <- unlist(tcr_counts[id_sam])[[sampled_to_column]]
    # find which REFERENCE sample has the lowest # of clonotypes
    # in the 2nd heatmap
    sample_to_reference_column <- unlist(tcr_counts_reference_columns[id_ref])[[sampled_to_reference_column]]
    # sample data
    options(future.globals.maxSize= 1500*1024^2)
    Sys.setenv(OPENBLAS_NUM_THREADS = "1")
    n_of_cores <- parallel::detectCores()
    plan(multisession, workers = n_of_cores - 4)
    cat("\n\t Sampling:", cell, chain)
    cat("\t Processing samples...")
    # Future for samples
    future_samples <- future({
      lapply(1:n_bootstrap_samples, function(i) {
        mapply(
          subsample_abundance,
          reference_dfs[id_sam],
          sample_to,
          target_type = "clone",
          SIMPLIFY = FALSE,
          seed = i
        )
      })
    })
    cat("\t Processing references...")
    # Future for references
    future_references <- future({
      lapply(1:n_bootstrap_samples, function(i) {
        mapply(
          subsample_abundance,
          reference_columns_dfs[id_ref],
          sample_to_reference_column,
          target_type = "clone",
          SIMPLIFY = FALSE,
          seed = i
        )
      })
    })
    cat("\t Done.")
    
    # Resolve both futures (wait for completion)
    reference_dfs_sampled <- value(future_samples)
    reference_columns_dfs_sampled <- value(future_references)
    plan(sequential)
    
    # count umis
    umi_counts <- lapply(reference_dfs, function(x) {
      sum(x$uniqueMoleculeCount)
    })
    
    umi_counts_row <- lapply(reference_dfs, function (x) {
      x %>% 
      slice_max(uniqueMoleculeCount,n = 100, with_ties = FALSE) %>%
      summarise(count = sum(uniqueMoleculeCount)) %>%
      pull(count)})
    
    umi_counts_reference_columns <- lapply(reference_columns_dfs, function(x) {
      sum(x$uniqueMoleculeCount)
    })
    
    # count umis and tcr after sampling
    umi_counts_sampled <- sapply(1:n_bootstrap_samples, function(i) {
      sapply(reference_dfs_sampled[[i]], sum)
    })
    umi_counts_sampled_median <- apply(umi_counts_sampled, 1, median)
    umi_counts_sampled_min <- apply(umi_counts_sampled, 1, quantile, probs = 0.025)
    umi_counts_sampled_max <- apply(umi_counts_sampled, 1, quantile, probs = 0.975)
    
    # count umis and tcr after sampling for reference columns
    umi_counts_sampled_reference_columns <- sapply(1:n_bootstrap_samples, function(i) {
      sapply(reference_columns_dfs_sampled[[i]], sum)
    })
    umi_counts_sampled_median_reference_columns <-
      apply(umi_counts_sampled_reference_columns, 1, median)
    umi_counts_sampled_min_reference_columns <-
      apply(umi_counts_sampled_reference_columns, 1, quantile, probs = 0.025)
    umi_counts_sampled_max_reference_columns <-
      apply(umi_counts_sampled_reference_columns, 1, quantile, probs = 0.975)
    
    tcr_counts_sampled <- lapply(1:n_bootstrap_samples, function(i) {
      lapply(reference_dfs_sampled[[i]], length)
    })
    tcr_counts_sampled_reference_columns <- lapply(1:n_bootstrap_samples, function(i) {
      lapply(reference_columns_dfs_sampled[[i]], length)
    })
    
    # count rows & create matrix with it
    reference_dfs_len <- length(reference_dfs)
    overlap_matrix_ref <-
      replicate(
        n_bootstrap_samples,
        matrix(nrow = reference_dfs_len, ncol = reference_dfs_len),
        simplify = FALSE
      )
    
    # count rows & create matrix with it for the reference cols
    reference_columns_dfs_len <- length(reference_columns_dfs)
    overlap_matrix_ref_columns <-
      replicate(
        n_bootstrap_samples,
        matrix(nrow = reference_dfs_len, ncol = reference_columns_dfs_len + 2), # overlap between references: +2
        simplify = FALSE
      )
    
    
    # fill in the matrix
    ## iterate rows
    for (i in 1:reference_dfs_len) {
      vectors_1 <- reference_dfs[[i]] %>% 
        slice_max(uniqueMoleculeCount,n = 100, with_ties = FALSE) %>%
        pull(nSeqCDR3)

      ns_i <- length(vectors_1)
      
      ## iterate columns for main heatmap
      for (j in 1:reference_dfs_len) {
        # Find overlaps
        vectors_2 <- lapply(1:n_bootstrap_samples, function(x) {
          names(reference_dfs_sampled[[x]][[j]])
        })
        ns_j <- lapply(vectors_2, length)
        # calculate overlap
        for (x in 1:n_bootstrap_samples) {
          overlap_cur <- length(intersect(vectors_1, vectors_2[[x]]))
          overlap_matrix_ref[[x]][i, j] <-
            (overlap_cur / ns_i) * 100
        }
      }
      
      ## iterate columns for reference column heatmap
      for (jcol in 1:reference_columns_dfs_len) {
        # Find overlaps
        vectors_2_cols <- lapply(1:n_bootstrap_samples, function(x) {
          names(reference_columns_dfs_sampled[[x]][[jcol]])
        })
        # calculate overlap
        for (x in 1:n_bootstrap_samples) {
          overlap_cur <- length(intersect(vectors_1, vectors_2_cols[[x]]))
          overlap_matrix_ref_columns[[x]][i, jcol] <-
            (overlap_cur / ns_i) * 100
        }
      }
      
      ## calculate sample vs reference1 and reference 2 and reference1 vs reference2 overlaps
      # Find overlaps
      vectors_2_cols <- lapply(1:n_bootstrap_samples, function(x) {
        intersect(
          names(reference_columns_dfs_sampled[[x]][[1]]),
          names(reference_columns_dfs_sampled[[x]][[2]]))
      })
      # calculate overlap
      for (x in 1:n_bootstrap_samples) {
        overlap_cur <- length(intersect(vectors_1, vectors_2_cols[[x]]))
        overlap_matrix_ref_columns[[x]][i, 3] <-
          (overlap_cur / ns_i) * 100
        # ref overlaps
        overlap_cur <- length(vectors_2_cols[[x]])
        overlap_matrix_ref_columns[[x]][i, 4] <-
          (overlap_cur / sample_to_reference_column) * 100
      }
    }
    
    
    # calculate overlap matrix metrics
    overlap_matrix_simple <- matrix(nrow = reference_dfs_len, ncol = reference_dfs_len)
    overlap_matrix_simple_min <- matrix(nrow = reference_dfs_len, ncol = reference_dfs_len)
    overlap_matrix_simple_max <- matrix(nrow = reference_dfs_len, ncol = reference_dfs_len)
    for (j in 1:reference_dfs_len) {
      # collect columns
      for (i in 1:reference_dfs_len) {
        vals <- sapply(overlap_matrix_ref, function(m)
          m[i, j])
        overlap_matrix_simple[i, j] <- median(vals, na.rm = TRUE)
        overlap_matrix_simple_min[i, j] <- quantile(vals, 0.025, na.rm = TRUE)
        overlap_matrix_simple_max[i, j] <- quantile(vals, 0.975, na.rm = TRUE)
      }
    }
    # calculate reference columns overlap matrix metrics
    overlap_matrix_simple_reference_columns <- matrix(nrow = reference_dfs_len, ncol = reference_columns_dfs_len+2)
    overlap_matrix_simple_min_reference_columns <- overlap_matrix_simple_reference_columns
    overlap_matrix_simple_max_reference_columns <- overlap_matrix_simple_reference_columns
    for (j in 1:(reference_columns_dfs_len+2)) {
      # collect columns
      for (i in 1:reference_dfs_len) {
        vals <- sapply(overlap_matrix_ref_columns, function(m)
          m[i, j])
        overlap_matrix_simple_reference_columns[i, j] <- median(vals, na.rm = TRUE)
        overlap_matrix_simple_min_reference_columns[i, j] <- quantile(vals, 0.025, na.rm = TRUE)
        overlap_matrix_simple_max_reference_columns[i, j] <- quantile(vals, 0.975, na.rm = TRUE)
      }
    }
    
    overlap_matrix_simple <- as.matrix(overlap_matrix_simple)
    col_names <- metadata %>% filter(gs.id %in% .env$id_sam) %>% pull(sample.id)
    col_names_reference_columns <-
      metadata %>% filter(gs.id %in% .env$id_ref) %>% pull(sample.id)
    row_names <- metadata %>% filter(gs.id %in% .env$id_sam) %>%
      pull(sample.id)
    
    # Colnames
    colnames(overlap_matrix_simple) <- paste(
      str_to_upper(substring(col_names, 5)),
      "\numi = ",
      paste0(
        umi_counts,
        " -> ",
        umi_counts_sampled_median,
        " (",
        umi_counts_sampled_min,
        "-",
        umi_counts_sampled_max,
        ") "
      ),
      "\ntcr = ",
      tcr_counts,
      " -> ",
      sample_to
    )
    
    # rownames
    rownames(overlap_matrix_simple) <- paste(
      str_to_upper(substring(col_names, 5)),
      "\numi = ",
      paste0(
        umi_counts,
        " -> ",
        umi_counts_row
      ),
      "\ntcr = top",
      100
    )
    # rownames(overlap_matrix_simple) <- colnames(overlap_matrix_simple)
    rownames(overlap_matrix_simple_reference_columns) <- rownames(overlap_matrix_simple)
    colnames(overlap_matrix_simple_reference_columns) <- 
      c(
        paste(
          str_to_upper(substring(col_names_reference_columns, 5)),
          "\numi = ",
          paste0(
            umi_counts_reference_columns,
            " -> ",
            umi_counts_sampled_median_reference_columns,
            " (",
            umi_counts_sampled_min_reference_columns,
            "-",
            umi_counts_sampled_max_reference_columns,
            ") "
          ),
          "\ntcr = ",
          tcr_counts_reference_columns,
          " -> ",
          sample_to_reference_column
        ),
        "Overlap with both reference repertoires",
        "Spleen CD69+ and CD69- Clonotype Overlap"
      )
    
    
    
    title <- paste(mouse, cell, chain)
    p_ref <-
      plot_overlap_matrix(
        overlap_matrix_simple,
        column_labels = colnames(overlap_matrix_simple),
        highlight_row = 0,
        highlight_column = sampled_to_column,
        heatmap_title = title,
        cluster = FALSE,
        col_rot = 90,
        font_size = 11,
        matrix_min = overlap_matrix_simple_min,
        matrix_max = overlap_matrix_simple_max
      )
    print(p_ref)
    # reference columns
    p_ref_columns <-
      plot_overlap_matrix(
        overlap_matrix_simple_reference_columns,
        column_labels = colnames(overlap_matrix_simple_reference_columns),
        highlight_row = 0,
        highlight_column = sampled_to_reference_column,
        heatmap_title = NULL,
        cluster = FALSE,
        col_rot = 90,
        font_size = 11,
        matrix_min = overlap_matrix_simple_min_reference_columns,
        matrix_max = overlap_matrix_simple_max_reference_columns,
        order = FALSE
      )
    print(p_ref + p_ref_columns)
    # save to file
    dir_path <- file.path("./figures/clone_overlap/top100_vs_sample_to_min/", threshold, group_name)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    # save plot to made files
    pdf(
      file = file.path(dir_path, paste(mouse, cell, chain, "pdf", sep = ".")),
      width = 13,
      height = 6.5
    )
    print(p_ref + p_ref_columns)
    dev.off()
  }

# Calculate overlap of top 100 to top "n" (min of samples) clones of samples
plot_overlap_100_and_top_min <-
  function(group, group_name, 
           n_bootstrap_samples = 1000, 
           query_n = 100,
           reference_n = 500) {
    # group <- c("bm", "spleen", "pp", "blood")
    # group_name <- "lymphoid"
    # get the ids of samples to be compared
    id_sam <- metadata %>% filter(
      mouse == .env$mouse,
      grepl(paste0(" ", .env$cell, " "), sample.id, ignore.case = TRUE),
      organ %in% .env$group
    ) %>% pull(gs.id)
    
    # get the ids of reference samples (+2 columns)
    id_ref <- metadata %>% filter(mouse == .env$mouse,
                                  grepl(
                                    paste0("spleen ", cell, " cd69(\\+|\\- 80k)$"),
                                    sample.id,
                                    ignore.case = TRUE
                                  )) %>% pull(gs.id)
    if (chain == "TRA") {
      reference_dfs <- clone_dfs_tra[id_sam]
      reference_columns_dfs <- clone_dfs_tra[id_ref]
    } else {
      reference_dfs <- clone_dfs_trb[id_sam]
      reference_columns_dfs <- clone_dfs_trb[id_ref]
    }
    
    # group by nCDR3 seq
    summarize_clones <- function(df) {
      df %>%
        group_by(nSeqCDR3) %>%
        summarise(uniqueMoleculeCount = sum(uniqueMoleculeCount),
                  .groups = "drop")
    }
    
    reference_dfs <- lapply(reference_dfs, summarize_clones)
    reference_columns_dfs <- lapply(reference_columns_dfs, summarize_clones)
    
    # count # of clonotypes
    tcr_counts <- lapply(reference_dfs, nrow)
    tcr_counts_reference_columns <- lapply(reference_columns_dfs, nrow)
    
    # find which sample has the lowest # of clonotypes
    sampled_to_column <-  which.min(unlist(tcr_counts[id_sam]))
    sampled_to_reference_column <-  which.min(unlist(tcr_counts_reference_columns[id_ref]))
    
    # find the sample size
    # sample_to <- unlist(tcr_counts[id_sam])[[sampled_to_column]]
    sample_to <- 
      ifelse(reference_n,
             reference_n,
             unlist(tcr_counts[id_sam])[[sampled_to_column]])
    
    # sampled_to_column is not allowed to be less than 100!
    if(sample_to < 100) sample_to <- 100
    
    
    # find which REFERENCE sample has the lowest # of clonotypes
    # in the 2nd heatmap
    sample_to_reference_column <- unlist(tcr_counts_reference_columns[id_ref])[[sampled_to_reference_column]]
    
    # pick rows' top 100 clones data
    reference_dfs_top100_row <-
      lapply(reference_dfs, function(x) {
        x %>% slice_max(uniqueMoleculeCount, n = query_n, with_ties = FALSE)
      })
    
    # pick columns' top n (sample_to) clones data
    reference_dfs_top_n_columns <-
      lapply(reference_dfs, function(x) {
        x %>% slice_max(uniqueMoleculeCount, n = sample_to, with_ties = FALSE)
      })
    
    # sample data
    options(future.globals.maxSize= 1500*1024^2)
    Sys.setenv(OPENBLAS_NUM_THREADS = "1")
    n_of_cores <- parallel::detectCores()
    plan(multisession, workers = n_of_cores - 4)
    cat("\n\t Sampling:", cell, chain)
    cat("\t Processing references...")
    # Future for references
    future_references <- future({
      lapply(1:n_bootstrap_samples, function(i) {
        mapply(
          subsample_abundance,
          reference_columns_dfs[id_ref],
          sample_to_reference_column,
          target_type = "clone",
          SIMPLIFY = FALSE,
          seed = i
        )
      })
    })
    cat("\t Done.")
    
    # Resolve both futures (wait for completion)
    reference_columns_dfs_sampled <- value(future_references)
    plan(sequential)
    
    # count umis
    umi_counts <- lapply(reference_dfs, function(x) {
      sum(x$uniqueMoleculeCount)
    })
    
    umi_counts_row <- lapply(reference_dfs_top100_row, function (x) {
      x %>% 
        summarise(count = sum(uniqueMoleculeCount)) %>%
        pull(count)})
    

    umi_counts_col <- lapply(reference_dfs_top_n_columns, function (x) {
      x %>% 
        summarise(count = sum(uniqueMoleculeCount)) %>%
        pull(count)})
    
    umi_counts_reference_columns <- lapply(reference_columns_dfs, function(x) {
      sum(x$uniqueMoleculeCount)
    })

    
    # count umis and tcr after sampling for reference columns
    umi_counts_sampled_reference_columns <- sapply(1:n_bootstrap_samples, function(i) {
      sapply(reference_columns_dfs_sampled[[i]], sum)
    })
    umi_counts_sampled_median_reference_columns <-
      apply(umi_counts_sampled_reference_columns, 1, median)
    umi_counts_sampled_min_reference_columns <-
      apply(umi_counts_sampled_reference_columns, 1, quantile, probs = 0.025)
    umi_counts_sampled_max_reference_columns <-
      apply(umi_counts_sampled_reference_columns, 1, quantile, probs = 0.975)
    
    tcr_counts_col <- lapply(reference_dfs_top_n_columns, function(i) {
      nrow(i)
    })
    tcr_counts_sampled_reference_columns <- lapply(1:n_bootstrap_samples, function(i) {
      lapply(reference_columns_dfs_sampled[[i]], length)
    })
    
    # count rows & create matrix with it
    reference_dfs_len <- length(reference_dfs_top100_row)
    overlap_matrix_ref <- matrix(nrow = reference_dfs_len, ncol = reference_dfs_len)

    # count rows & create matrix with it for the reference cols
    reference_columns_dfs_len <- length(reference_columns_dfs)
    overlap_matrix_ref_columns <-
      replicate(
        n_bootstrap_samples,
        matrix(nrow = reference_dfs_len, ncol = reference_columns_dfs_len + 2), # overlap between references: +2
        simplify = FALSE
      )
    
    
    # fill in the matrix
    ## iterate rows
    for (i in 1:reference_dfs_len) {
      vectors_1 <- reference_dfs_top100_row[[i]] %>%
        pull(nSeqCDR3)
      
      ns_i <- length(vectors_1)
      
      ## iterate columns for main heatmap
      for (j in 1:reference_dfs_len) {
        # Find overlaps
        vectors_2 <- reference_dfs_top_n_columns[[j]] %>% pull(nSeqCDR3)
        
        ns_j <- length(vectors_2)
        
        # calculate overlap
        overlap_cur <- length(intersect(vectors_1, vectors_2))
          overlap_matrix_ref[i, j] <-
            (overlap_cur / ns_i) * 100
      }
      
      ## iterate columns for reference column heatmap
      for (jcol in 1:reference_columns_dfs_len) {
        # Find overlaps
        vectors_2_cols <- lapply(1:n_bootstrap_samples, function(x) {
          names(reference_columns_dfs_sampled[[x]][[jcol]])
        })
        # calculate overlap
        for (x in 1:n_bootstrap_samples) {
          overlap_cur <- length(intersect(vectors_1, vectors_2_cols[[x]]))
          overlap_matrix_ref_columns[[x]][i, jcol] <-
            (overlap_cur / ns_i) * 100
        }
      }
      
      ## calculate sample vs reference1 and reference2 and reference1 vs reference2 overlaps
      # Find overlaps
      vectors_2_cols <- lapply(1:n_bootstrap_samples, function(x) {
        intersect(
          names(reference_columns_dfs_sampled[[x]][[1]]),
          names(reference_columns_dfs_sampled[[x]][[2]]))
      })
      # calculate overlap
      for (x in 1:n_bootstrap_samples) {
        overlap_cur <- length(intersect(vectors_1, vectors_2_cols[[x]]))
        overlap_matrix_ref_columns[[x]][i, 3] <-
          (overlap_cur / ns_i) * 100
        # ref overlaps
        overlap_cur <- length(vectors_2_cols[[x]])
        overlap_matrix_ref_columns[[x]][i, 4] <-
          (overlap_cur / sample_to_reference_column) * 100
      }
    }
    
    # calculate overlap matrix metrics
    overlap_matrix_simple <- matrix(nrow = reference_dfs_len, ncol = reference_dfs_len)
    overlap_matrix_simple_min <- matrix(nrow = reference_dfs_len, ncol = reference_dfs_len)
    overlap_matrix_simple_max <- matrix(nrow = reference_dfs_len, ncol = reference_dfs_len)
    for (j in 1:reference_dfs_len) {
      # collect columns
      for (i in 1:reference_dfs_len) {
        vals <- overlap_matrix_ref[i, j]
        overlap_matrix_simple[i, j] <- median(vals, na.rm = TRUE)
        overlap_matrix_simple_min[i, j] <- quantile(vals, 0.025, na.rm = TRUE)
        overlap_matrix_simple_max[i, j] <- quantile(vals, 0.975, na.rm = TRUE)
      }
    }
    # calculate reference columns overlap matrix metrics
    overlap_matrix_simple_reference_columns <- matrix(nrow = reference_dfs_len, ncol = reference_columns_dfs_len+2)
    overlap_matrix_simple_min_reference_columns <- overlap_matrix_simple_reference_columns
    overlap_matrix_simple_max_reference_columns <- overlap_matrix_simple_reference_columns
    for (j in 1:(reference_columns_dfs_len+2)) {
      # collect columns
      for (i in 1:reference_dfs_len) {
        vals <- sapply(overlap_matrix_ref_columns, function(m)
          m[i, j])
        overlap_matrix_simple_reference_columns[i, j] <- median(vals, na.rm = TRUE)
        overlap_matrix_simple_min_reference_columns[i, j] <- quantile(vals, 0.025, na.rm = TRUE)
        overlap_matrix_simple_max_reference_columns[i, j] <- quantile(vals, 0.975, na.rm = TRUE)
      }
    }
    
    overlap_matrix_simple <- as.matrix(overlap_matrix_simple)
    col_names <- metadata %>% filter(gs.id %in% .env$id_sam) %>% pull(sample.id)
    col_names_reference_columns <-
      metadata %>% filter(gs.id %in% .env$id_ref) %>% pull(sample.id)
    row_names <- metadata %>% filter(gs.id %in% .env$id_sam) %>%
      pull(sample.id)
    
    # Colnames
    colnames(overlap_matrix_simple) <- paste(
      str_to_upper(substring(col_names, 5)),
      "\numi = ",
      paste0(
        umi_counts,
        " -> ",
        umi_counts_col
      ),
      "\ntcr = ",
      tcr_counts,
      " ->  top",
      sample_to
    )
    
    # rownames
    rownames(overlap_matrix_simple) <- paste(
      str_to_upper(substring(col_names, 5)),
      "\numi = ",
      paste0(
        umi_counts,
        " -> ",
        umi_counts_row
      ),
      "\ntcr = ",
      tcr_counts,
      " ->  top",
      100
    )
    
    # rownames(overlap_matrix_simple) <- colnames(overlap_matrix_simple)
    rownames(overlap_matrix_simple_reference_columns) <- rownames(overlap_matrix_simple)
    colnames(overlap_matrix_simple_reference_columns) <- 
      c(
        paste(
          str_to_upper(substring(col_names_reference_columns, 5)),
          "\numi = ",
          paste0(
            umi_counts_reference_columns,
            " -> ",
            umi_counts_sampled_median_reference_columns,
            " (",
            umi_counts_sampled_min_reference_columns,
            "-",
            umi_counts_sampled_max_reference_columns,
            ") "
          ),
          "\ntcr = ",
          tcr_counts_reference_columns,
          " -> ",
          sample_to_reference_column
        ),
        "Overlap with both reference repertoires",
        "Spleen CD69+ and CD69- Clonotype Overlap"
      )
    
    
    
    title <- paste(mouse, cell, chain)
    p_ref <-
      plot_overlap_matrix(
        overlap_matrix_simple,
        column_labels = colnames(overlap_matrix_simple),
        highlight_row = 0,
        highlight_column = sampled_to_column,
        heatmap_title = paste(title, ""),
        cluster = FALSE,
        col_rot = 90,
        font_size = 11,
        matrix_min = overlap_matrix_simple_min,
        matrix_max = overlap_matrix_simple_max
      )
    print(p_ref)
    # reference columns
    p_ref_columns <-
      plot_overlap_matrix(
        overlap_matrix_simple_reference_columns,
        column_labels = colnames(overlap_matrix_simple_reference_columns),
        highlight_row = 0,
        highlight_column = sampled_to_reference_column,
        heatmap_title = NULL,
        cluster = FALSE,
        col_rot = 90,
        font_size = 11,
        # matrix_min = overlap_matrix_simple_min_reference_columns,
        # matrix_max = overlap_matrix_simple_max_reference_columns,
        order = FALSE
      )
    print(p_ref + p_ref_columns)
    # save to file
    dir_path <- file.path("./figures/clone_overlap/top100_vs_top_n/", threshold, group_name)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    dir_path_csv <- file.path("./figures/clone_overlap/top100_vs_top_n/", threshold, group_name, "overlap_matrix_data")
    if (!dir.exists(dir_path_csv)) {
      dir.create(dir_path_csv, recursive = TRUE)
    }
    # save plot to made files
    pdf(
      file = file.path(dir_path, paste(mouse, cell, chain, "pdf", sep = ".")),
      width = 13,
      height = 6.5
    )
    print(p_ref + p_ref_columns)
    dev.off()
    # export overlap matrix
    first_parts <- sapply(str_split(colnames(overlap_matrix_simple), "\numi"), function(x) x[1])
    
    colnames(overlap_matrix_ref) <- first_parts
    rownames(overlap_matrix_ref) <- first_parts
    write.csv(x = overlap_matrix_ref,
              file.path(dir_path_csv, paste(mouse, cell, chain, "csv", sep = ".")))
    
    # return overlap matrix
    column_sampled_to <- 
      ifelse(reference_n, reference_n, sampled_to_column)
    
    return(list(overlap_matix = overlap_matrix_ref, column_sampled_to = column_sampled_to))
  }


    
# Calculate overlap of top 100 to top "n" (min of samples) clones of samples
plot_overlap_top_100_and_top_200 <-
  function(group, group_name, 
           query_n = 100,
           reference_n = 200) {
    # group <- c("bm", "spleen", "pp", "blood")
    # group_name <- "lymphoid"
    # get the ids of samples to be compared
    id_sam <- metadata %>% filter(
      mouse == .env$mouse,
      grepl(paste0(" ", .env$cell, " "), sample.id, ignore.case = TRUE),
      organ %in% .env$group
    ) %>% pull(gs.id)
    
    # Get the chain data
    if (chain == "TRA") {
      reference_dfs <- clone_dfs_tra[id_sam]
    } else {
      reference_dfs <- clone_dfs_trb[id_sam]
    }
    
    # group by nCDR3 seq
    summarize_clones <- function(df) {
      df %>%
        group_by(nSeqCDR3) %>%
        summarise(uniqueMoleculeCount = sum(uniqueMoleculeCount),
                  .groups = "drop")
    }
    
    reference_dfs <- lapply(reference_dfs, summarize_clones)

    # count # of clonotypes
    tcr_counts <- lapply(reference_dfs, nrow)

    # find which sample has the lowest # of clonotypes
    sampled_to_column <-  which.min(unlist(tcr_counts[id_sam]))

    # find the sample size
    # sample_to <- unlist(tcr_counts[id_sam])[[sampled_to_column]]
    selected_top_n_clones_in_column <- 
      ifelse(reference_n,
             reference_n,
             unlist(tcr_counts[id_sam])[[sampled_to_column]])
    
    # selected_top_n_clones_in_column is not allowed to be less than 100!
    if(selected_top_n_clones_in_column < 100) selected_top_n_clones_in_column <- 100
    
    # pick rows' top 100 clones data
    reference_dfs_top100_row <-
      lapply(reference_dfs, function(x) {
        if(nrow(x) < query_n) return(x %>% slice_max(uniqueMoleculeCount, n = 0, with_ties = FALSE))
        x %>% slice_max(uniqueMoleculeCount, n = query_n, with_ties = FALSE)
      })
    
    # pick columns' top n clones data
    reference_dfs_top_n_columns <-
      lapply(reference_dfs, function(x) {
        if(nrow(x) < selected_top_n_clones_in_column) return(x %>% slice_max(uniqueMoleculeCount, n = 0, with_ties = FALSE))
        x %>% slice_max(uniqueMoleculeCount, n = selected_top_n_clones_in_column, with_ties = FALSE)
      })
    
    # count umis
    umi_counts <- lapply(reference_dfs, function(x) {
      sum(x$uniqueMoleculeCount)
    })
    
    umi_counts_row <- lapply(reference_dfs_top100_row, function (x) {
      x %>% 
        summarise(count = sum(uniqueMoleculeCount)) %>%
        pull(count)})
    
    
    umi_counts_col <- lapply(reference_dfs_top_n_columns, function (x) {
      x %>% 
        summarise(count = sum(uniqueMoleculeCount)) %>%
        pull(count)})

    # count tcrs
    tcr_counts_col <- lapply(reference_dfs_top_n_columns, function(i) {
      nrow(i)
    })
    
    # count rows & create matrix with it
    reference_dfs_len <- length(reference_dfs_top100_row)
    overlap_matrix_ref <- matrix(nrow = reference_dfs_len, ncol = reference_dfs_len)
    
    
    # fill in the matrix
    ## iterate rows
    for (i in 1:reference_dfs_len) {
      
      vectors_1 <- reference_dfs_top100_row[[i]] %>%
        pull(nSeqCDR3)
      
      ns_i <- length(vectors_1)
      
      if (ns_i == 0) {
        overlap_matrix_ref[i, ] <- NA
        next
      }
      
      ## iterate columns for main heatmap
      for (j in 1:reference_dfs_len) {
        
        # Find overlaps
        vectors_2 <- reference_dfs_top_n_columns[[j]] %>% pull(nSeqCDR3)
        
        ns_j <- length(vectors_2)
        
        # Set overlap to NA if column is dropped out
        if (ns_j == 0) {
          overlap_matrix_ref[i, j] <- NA
          next
        }
        
        # calculate overlap
        overlap_cur <- length(intersect(vectors_1, vectors_2))
        overlap_matrix_ref[i, j] <-
          (overlap_cur / ns_i) * 100
        
      }
      
    }
    
    # calculate overlap matrix metrics
    overlap_matrix_simple <- matrix(nrow = reference_dfs_len, ncol = reference_dfs_len)
    overlap_matrix_simple_min <- matrix(nrow = reference_dfs_len, ncol = reference_dfs_len)
    overlap_matrix_simple_max <- matrix(nrow = reference_dfs_len, ncol = reference_dfs_len)
    for (j in 1:reference_dfs_len) {
      # collect columns
      for (i in 1:reference_dfs_len) {
        vals <- overlap_matrix_ref[i, j]
        overlap_matrix_simple[i, j] <- median(vals, na.rm = TRUE)
        overlap_matrix_simple_min[i, j] <- quantile(vals, 0.025, na.rm = TRUE)
        overlap_matrix_simple_max[i, j] <- quantile(vals, 0.975, na.rm = TRUE)
      }
    }
    
    overlap_matrix_simple <- as.matrix(overlap_matrix_simple)
    
    col_names <- metadata %>% filter(gs.id %in% .env$id_sam) %>% pull(sample.id)
    
    row_names <- metadata %>% filter(gs.id %in% .env$id_sam) %>%
      pull(sample.id)
    
    # Colnames
    colnames(overlap_matrix_simple) <- paste(
      str_to_upper(substring(col_names, 5)),
      "\numi = ",
      paste0(
        umi_counts,
        " -> ",
        umi_counts_col
      ),
      "\ntcr = ",
      tcr_counts,
      " ->  top",
      selected_top_n_clones_in_column
    )
    
    # rownames
    rownames(overlap_matrix_simple) <- paste(
      str_to_upper(substring(col_names, 5)),
      "\numi = ",
      paste0(
        umi_counts,
        " -> ",
        umi_counts_row
      ),
      "\ntcr = ",
      tcr_counts,
      " ->  top",
      query_n
    )
    

    title <- paste(mouse, cell, chain)
    p_ref <-
      plot_overlap_matrix(
        overlap_matrix_simple,
        column_labels = colnames(overlap_matrix_simple),
        highlight_row = 0,
        highlight_column = sampled_to_column,
        heatmap_title = paste(title, ""),
        cluster = FALSE,
        col_rot = 90,
        font_size = 11,
        matrix_min = overlap_matrix_simple_min,
        matrix_max = overlap_matrix_simple_max
      )
    print(p_ref)

    # save to file
    base_path_analysis <- paste0("./figures/clone_overlap/top100_vs_top_", selected_top_n_clones_in_column)
    dir_path <- file.path(base_path_analysis, threshold, group_name)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    dir_path_csv <- file.path(base_path_analysis, threshold, group_name, "overlap_matrix_data")
    if (!dir.exists(dir_path_csv)) {
      dir.create(dir_path_csv, recursive = TRUE)
    }
    # save plot to made files
    pdf(
      file = file.path(dir_path, paste(mouse, cell, chain, "pdf", sep = ".")),
      width = 13,
      height = 6.5
    )
    print(p_ref)
    dev.off()
    
    # export overlap matrix
    first_parts <- sapply(str_split(colnames(overlap_matrix_simple), "\numi"), function(x) x[1])
    
    colnames(overlap_matrix_ref) <- first_parts
    rownames(overlap_matrix_ref) <- first_parts
    write.csv(x = overlap_matrix_ref,
              file.path(dir_path_csv, paste(mouse, cell, chain, "csv", sep = ".")))
    
    # return overlap matrix
    column_sampled_to <- 
      ifelse(reference_n, reference_n, sampled_to_column)
    
    return(list(overlap_matix = overlap_matrix_ref, column_sampled_to = column_sampled_to))
  }









plot_overlap_network <-
  function(clone_dfs = clone_dfs_tra, top_n_clones = 10, chain = chain, group){
  # Makes the sequence newtork of top shared clones
  sample_ids <- metadata %>% filter(gs.id %in% names(clone_dfs)) %>% pull(sample.id)
  
  names(clone_dfs) <- sample_ids
  
  selected_dfs <- map(clone_dfs, ~ select(.x, nSeqCDR3, uniqueMoleculeCount))
  renamed_dfs <- map2(selected_dfs, sample_ids, ~ {
    df <- .x
    df <- df %>% slice_max(uniqueMoleculeCount, n = top_n_clones, with_ties = FALSE)
    colnames(df)[colnames(df) == "uniqueMoleculeCount"] <- paste0(str_trim(word(.y, start = 2, end = -1)))
    df
  })
  
  renamed_dfs_700 <- map2(selected_dfs, sample_ids, ~ {
    df <- .x
    df <- df %>% slice_max(uniqueMoleculeCount, n = 700, with_ties = FALSE)
    colnames(df)[colnames(df) == "uniqueMoleculeCount"] <- paste0(str_trim(word(.y, start = 2, end = -1)))
    df
  })

  
  merged_df <- reduce(renamed_dfs, full_join, by = "nSeqCDR3")
  merged_df_700 <- reduce(renamed_dfs_700, full_join, by = "nSeqCDR3")
  
  subs <- merged_df_700 %>%
    filter(nSeqCDR3 %in% merged_df$nSeqCDR3) %>%
    slice(match(merged_df$nSeqCDR3, nSeqCDR3)) %>%
    filter(!is.na(nSeqCDR3))
  # subs[merged_df$nSeqCDR3,]
  # fdf <- merged_df %>% left_join(merged_df_700, by = "nSeqCDR3", suffix = c("s", "l"))
  
  # Convert to adjacency matrix - drop nSeqCDR3 column and convert to matrix
  adj_matrix <- merged_df %>% 
    column_to_rownames(var = "nSeqCDR3") %>%
    as.matrix()
  
  adj_matrix_f <- subs %>% 
    column_to_rownames(var = "nSeqCDR3") %>%
    as.matrix()
  # Replace NA with 0
  adj_matrix[is.na(adj_matrix)] <- 0
  adj_matrix_f[is.na(adj_matrix_f)] <- 0
  
  # If values are counts, consider binarizing edges (optional)
  adj_matrix_binary <- (adj_matrix > 0) * 1
  adj_matrix_binary_f <- (adj_matrix_f > 0) * 2
  colnames(adj_matrix_binary_f) <- colnames(adj_matrix_binary)
  
  final_df <- adj_matrix_binary + adj_matrix_binary_f
  # write.csv(igraph::as_data_frame(network, what = "edges"), file = "~/Desktop/test.csv")
  # Create bipartite graph from biadjacency matrix
  # network <- graph_from_biadjacency_matrix(adj_matrix_binary,directed = FALSE)
  network <- graph_from_biadjacency_matrix(final_df,directed = FALSE,weighted = TRUE)
  
  V(network)$size <- ifelse(V(network)$type, 2, 1)  
  E(network)$edge_type <- E(network)$weight
  
  dir_path <- file.path("./figures/clone_overlap/network/", threshold)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  # save plot to made files
  # pn <- ggraph(network, layout = "kk") +
  #   geom_edge_link(alpha = 0.4) +
  #   geom_node_point(aes(color = factor(type))) +
  #   geom_node_text(aes(label = ifelse(V(network)$type, V(network)$name, NA)), repel = TRUE, size = 3) +
  #   theme_void() +
  #   theme(legend.position = "none")
  
  pn <- 
    ggraph(network, layout = "auto") +
    geom_edge_link(aes(color = factor(edge_type)),alpha = 0.4) +
    geom_node_point(aes(color = factor(type), size = 2)) +  # map size aesthetic
    geom_node_text(aes(label = ifelse(type, name, NA)), repel = TRUE, size = 2) +
    theme_void() +
    theme(legend.position = "none") +
    scale_size_identity() +  # interpret size values literally
    scale_edge_color_manual(values = c("2" = "blue", "3" = "red"),
                            name = "Edge Type",
                            labels = c("2" = "Type 2", "3" = "Type 3"))
  # pn
  ggsave(
    filename = file.path(dir_path, paste(mouse, chain, group, ".color.pdf", sep = ".")),
    plot = pn,
    width = 6,
    height = 6,
    device = "pdf"
  )
  
  write.csv(igraph::as_data_frame(network, what = "edges"), 
            file = file.path(dir_path, paste(mouse, chain, group, ".color.csv", sep = ".")))
  
  }



plotLargeOverlaps <-
  function(mouse = "M23", gs.ids.current = cd4_gs.id,
  chain.current = "TRA", top.n = 10, ySclae = "log"){

    # Grab top 10 seqs of given ids and chain
    all.seqs <- abundances_df

    # define colors per seq
    all_unique_seqs <- unique(all.seqs$nSeqCDR3)
    n_seqs <- length(all_unique_seqs)

    # Generate a palette with high variance
    set.seed(42) # Ensure the random shuffle is the same every time you run the script
    large_palette <- scales::hue_pal()(n_seqs) 
    large_palette <- sample(large_palette) # Shuffle so similar colors aren't next to each other

    # Map them to the sequences
    color_map <- setNames(large_palette, all_unique_seqs)

    # Generate abetter colo mapping
    dist_matrix <- stringdistmatrix(all_unique_seqs, all_unique_seqs, method = "hamming")

    # 2. Perform Hierarchical Clustering to find the best order
    # This groups similar sequences together
    hc <- hclust(as.dist(dist_matrix), method = "average")
    ordered_indices <- hc$order
    ordered_seqs <- all_unique_seqs[ordered_indices]

    # 3. Generate a sequential palette (not shuffled)
    # hue_pal()(n) provides a rainbow-like progression
    sequential_palette <- hue_pal()(length(all_unique_seqs))

    # 4. Map the ordered sequences to the sequential colors
    # This ensures sequence A and sequence B (which are similar) get colors 
    # that are next to each other in the color wheel.
    color_map <- setNames(sequential_palette, ordered_seqs)



    all.seqs$organ <- 
      factor(
        all.seqs$organ,
        levels = c("blood", "spleen", "bm", "pp", "liver", "lung", "skin")
        )

    all.seqs <- 
      all.seqs %>%
        filter(
          gs.id %in% gs.ids.current,
          chain == chain.current,
          cell.numbers != 20000
          ) %>%
        arrange(organ) %>%
        mutate(organ.cell = interaction(organ, Cell, sep = "\n", lex.order = TRUE, drop = TRUE)) %>%        # 3. Filter for your top sequences
        group_by(gs.id, mouse, chain, organ, cell.type, Cell, organ.cell, nSeqCDR3) %>%
        summarise(Abundance = sum(Abundance), UMI_count = sum(UMI_count)) %>%
        mutate(
          clone.rank = rank(-Abundance, ties.method = "first"),
          in.top.n = clone.rank <= top.n
        ) %>%
        arrange(organ, clone.rank) 

    # Plot per organ.cell group
    current.organ.cell.groups <-
      unique(all.seqs %>% pull(organ.cell))
    
    # list of plots to be saved
    plots.current <- list()
    plots.current[chain.current] <- list()

    for (current.organ.cell.group in current.organ.cell.groups) {
      # current.organ.cell.group <- current.organ.cell.groups[[2]]
      # Pick top 10 (n) of the group
      top.10.seqs.current <- 
        all.seqs %>%
        filter(organ.cell == current.organ.cell.group) %>%
        slice_max(Abundance,n=top.n, with_ties = FALSE) %>%
      ungroup()

      organ.current <- top.10.seqs.current$organ[[1]]
      cell.current <- str_replace_all(top.10.seqs.current$Cell[[1]], "\n", " ")
      
      # fill sequences outside top 10 abundances from the main abundances_df
      top.10.seqs.list <- top.10.seqs.current %>% pull(nSeqCDR3)
      bottom.seqs.found <-
      all.seqs %>%
        filter(organ.cell != current.organ.cell.group) %>%
        filter(nSeqCDR3 %in% top.10.seqs.list) %>%
        ungroup()
      

      # Combine sequences with measured abundances and fill zeros if not measured
      # in other groups
      seqs.to.plot <- 
        top.10.seqs.current %>%
        bind_rows(bottom.seqs.found) %>%
        tidyr::complete(nSeqCDR3, organ.cell, fill = list(Abundance = 0, in.top.n = FALSE)) %>%
        mutate()




      # 1. Identify your unique x-axis groups
      x_groups <- unique(seqs.to.plot$organ.cell)
      # 2. Find the index of every other group
      even_groups <- which(seq(along = x_groups) %% 2 == 0)

      # plot the current organ.cell group
      c.plot <- 
        ggplot(seqs.to.plot,
          aes(
            x = organ.cell, y = Abundance,
            group = nSeqCDR3, color = nSeqCDR3,
            fill = nSeqCDR3
            )) +
        geom_rect(
          data = data.frame(x = even_groups),
          aes(xmin = x - 0.5, xmax = x + 0.5, ymin = -Inf, ymax = Inf),
          fill = "grey90", alpha = 0.5, inherit.aes = FALSE
        ) +
        geom_rect(
          data = data.frame(x = current.organ.cell.group),
          aes(xmin = x - 0.5, xmax = x + 0.5, ymin = -Inf, ymax = Inf),
          fill = "darkred", alpha = 0.5, inherit.aes = FALSE
        ) +
        geom_line(alpha = 0.7, linewidth = 1) +
        scale_color_manual(values = color_map) +
        # --- LOG OPTION FOR ZEROS --
        scale_y_continuous(
          # sigma = 0.0005 tells ggplot to start 'logging' very early
          trans = pseudo_log_trans(sigma = 0.0001, base = 10),
          breaks = c(0, 0.001, 0.01, 0.1, 1, 10, 100),
          limits = c(0, 100), 
          labels = function(x) {
            # 1. Handle 0 separately for a clean look
            if (length(x) == 0) return(character(0))
            labs <- ifelse(x == 0, "0", scales::label_log()(x))
            
            # 2. Remove the leading "+" from the start of the label
            # and remove "+" from inside the exponent (the math format part)
            labs <- gsub("^\\+", "", labs)      # Removes leading +
            labs <- gsub("\\(\\+", "\\(", labs) # Removes + inside power notation if present
            
            return(parse(text = labs))
          },
          expand = expansion(mult = c(0.02, 0))
        ) +
        # ----------------------------
        # Solid points for in top
        geom_point(
          data = subset(seqs.to.plot, in.top.n == TRUE),
          aes(shape = in.top.n, fill = nSeqCDR3),
          #fill = "white",
          alpha = 1,
          size = 2
          ) +
        # Hollow points for exactly 0
        geom_point(
          data = subset(seqs.to.plot, in.top.n == FALSE),
          aes(shape = in.top.n),
          fill = "white",
          alpha = 1,
          size = 2
          ) +
        scale_fill_manual(values = color_map) +
        scale_shape_manual(values = c("TRUE" = 21, "FALSE" = 21)) +
        theme_classic(base_size = 12) +
        theme(
          legend.position = "none",
          panel.grid.major.y = element_line(color = "grey95"),
          axis.title.x = element_blank()
        ) +
        # ylab(str_replace_all(current.organ.cell.group, "\n", " ")) +
        annotate("text", label = str_replace_all(current.organ.cell.group, "\n", " "),
        x = -Inf, y = Inf, 
        hjust = -0.2, vjust = 1.5, size = 5)

        plot_path <- 
          file.path(
            "~", "Dropbox", "Research", "tissue_tcr", "figures", "clone_overlap", "top10_per_organ", threshold,
            paste(mouse, organ.current, cell.current, chain.current, "pdf", sep = ".")
            )
        
        ggsave(plot_path, c.plot, width = 30, units = "cm", create.dir = TRUE)

        # save plots to a list
        clean_name <- stringr::str_replace_all(current.organ.cell.group, "\n", ".")
        plots.current[[chain.current]][[clean_name]] <- c.plot
        write.csv(x = abundances_df, "~/Sync/test.abundance.csv", row.names = FALSE, sep = "\t")
    }
    return(plots.current)
  }

remove_x_axis <- function(plot_list) {
  n <- length(plot_list)
  if (n > 1) {
    for (i in 1:(n-1)) {
      plot_list[[i]] <- plot_list[[i]] + 
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
    }
  }
  return(plot_list)
}


# scientific notation in plots
fancy_scientific <- function(l) {
  l <- format(l, scientific = TRUE)            # Convert to scientific notation
  l <- gsub("^(.*)e", "'\\1'e", l)            # Keep digits before exponent
  l <- gsub("e", "%*%10^", l)                 # Replace e with x10^
  parse(text = l)                               # Parse as plotmath expression
}

my_theme <-
  function(font_size = 8){
    list(
      theme_classic(base_size = font_size, base_rect_size = 0) +
        theme(
          panel.grid.major.x = element_blank(),
          # remove vertical major grid lines
          panel.grid.major.y = element_line(color = "gray80", linewidth = 0.2),
          # faint gray horizontal lines
          panel.grid.minor = element_blank(),
          coord_cartesian(clip = "off"),
          axis.text=element_text(size=font_size),
          legend.position="bottom"
        )
    )
  }

my_theme_p1 <-
  list(
    theme_classic(base_size = 9, base_rect_size = 0) +
      theme(
        panel.grid.major.x = element_line(color = "gray80", linewidth = 0.2),
        # remove vertical major grid lines
        panel.grid.major.y = element_line(color = "gray80", linewidth = 0.2),
        # faint gray horizontal lines
        panel.grid.minor = element_blank(),
        coord_cartesian(clip = "off"),
        axis.text=element_text(size=5),
        legend.position="bottom"
      )
  )

my_theme_y <-
  list(
    theme_classic(base_size = 9, base_rect_size = 0) +
      theme(
        panel.grid.major.y = element_blank(),
        # remove vertical major grid lines
        panel.grid.major.x = element_line(color = "gray80", linewidth = 0.2),
        # faint gray horizontal lines
        panel.grid.minor = element_blank(),
        coord_cartesian(clip = "off"),
        axis.text=element_text(size=5)
      )
  )


# Plot themes individually
clone_abundance_summed_plot_theme <- 
  clone_abundance_summed_plot_theme <- theme(
    strip.text = element_text(size = 7),
    axis.text = element_text(size = 7),
    strip.text.y.left = element_text(angle = 0),
    legend.position = "bottom",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.margin = margin(0, 0, 0, 0, "pt"),
    legend.box.margin = margin(0, 0, 0, 0, "pt"),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white"),
    strip.background = element_rect(fill = "transparent", color = "white")
  )
 
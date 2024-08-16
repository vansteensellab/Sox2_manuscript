# Christine: 
# The function below takes a bigwig as genomic ranges (imported with rtracklayer::import.bw) and bins a given region of interest into fixed size bins. 
# Then returns a new genomic ranges with the mean score for those bins. 
BinAv_on_ROI <- function(data_grange, region_oi_gr, binsize){
  #get the chromosome(s) to tile
  chrOI = as.character(unique(seqnames(region_oi_gr))) 
  
  #tile a fake genome with only these chromosomes
  tib_seqlengths = as_tibble(region_oi_gr) %>%
    group_by(seqnames) %>%
    summarize(endpoint = max(end))
  
  seqlength_OI = tib_seqlengths$endpoint
  names(seqlength_OI) = tib_seqlengths$seqnames
  
  #create bins
  bins <- tileGenome(seqlengths = seqlength_OI,
                     tilewidth = binsize,
                     cut.last.tile.in.chrom=TRUE)
  
  #calculate 'coverage' from the BigWig granges
  cov.bw = GenomicRanges::coverage(data_grange, weight = "score")
  cov.bw_chrOI = cov.bw[chrOI]
  
  #calculate binned average
  BinAv_all = binnedAverage(bins, cov.bw_chrOI, "score_binned")
  BinAv_ROI = subsetByOverlaps(BinAv_all, region_oi_gr)
  BinAv_ROI
}


# From Tom: (slight adaptation by Christine)
# Takes a tibble of binned data (with fixed bin size), returns a plot (which might be smoothened)
# adaptation 20240301: add an option to give a fixed y axis 
#  (e.g. when you want to manually plot multiple ranges on the same scale)
# adaptation 20240621: option to remove the line at y = 0 (for ChIP/ATAC tracks, which always start at 0 and the thick line gets in the way)
PlotDataTracksLight <- function(input_tib, exp_names, color_groups, 
                                centromeres = NULL, plot_chr = "chr1", 
                                plot_start = 1, plot_end = 500e6,
                                smooth = 1, color_list = NULL,
                                fix_yaxis = F, aspect_ratio = NULL,
                                ylimits_fixed = NULL,
                                lighten_negative = F, raster = F,
                                squish = F, plot_y0_line = T) {
  
  # Get the scores for the samples
  tib_plot <- input_tib %>%
    dplyr::select(seqnames, start, end, 
                  all_of(exp_names))
  
  if (smooth > 1) {
    tib_plot <- tib_plot %>%
      mutate_at(all_of(exp_names), 
                runmean, k = smooth)
  }
  
  # Filter for plotting window
  tib_plot <- tib_plot %>%
    filter(seqnames == plot_chr,
           start >= plot_start,
           end <= plot_end)
  
  # Gather
  tib_gather <- tib_plot %>%
    gather(key, value, -seqnames, -start, -end) %>%
    mutate(key = factor(key, levels = exp_names),
           fill_column = color_groups[match(key, names(color_groups))],
           fill_column = factor(fill_column, levels = unique(color_groups)))
  
  # If desired, make negative values a lighter shade to improve visualization
  if (lighten_negative) {
    tib_gather <- tib_gather %>%
      mutate(fill_column = interaction(fill_column,
                                       value < 0))
  }
  
  # Plot
  if( is.null(ylimits_fixed)) {
    ylimits <- quantile(tib_gather$value, c(0.001, 0.999), na.rm = T)
  } else {
    ylimits = ylimits_fixed
  }
  
  fill_column <- NULL
  
  plt <- tib_gather %>% 
    ggplot(aes(fill = fill_column))
  
  
  # Set all counts tracks to the same limits
  if (raster) {
    plt <- plt + 
      rasterize(geom_rect(aes(xmin = start / 1e6, xmax = end / 1e6, 
                              ymin = 0, ymax = value)),
                dpi = 300)
  } else {
    plt <- plt + 
      geom_rect(aes(xmin = start / 1e6, xmax = end / 1e6, 
                    ymin = 0, ymax = value))
  }
  
  if(plot_y0_line){
    plt <- plt + 
      geom_hline(yintercept = 0, linewidth = 0.5)
  }
  
  plt <- plt + 
    facet_grid(key ~ ., scales = "free_y") +
    xlab(paste0(plot_chr, " (Mb)")) +
    ylab("Score") +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0.025, 0.025)) +
    theme_classic()
  
  #Christine: if you want to squish panels closer in the final plot_grid they shouldn't have big margins
  if(squish){
    plt <- plt +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  }
  
  # Centromeres
  if (! is.null(centromeres)) {
    centromeres.plt <- as_tibble(centromeres) %>%
      filter(seqnames == plot_chr) %>%
      gather(key_centromeres, value, start, end)
    
    plt <- plt +
      geom_line(data = centromeres.plt, 
                aes(x = value / 1e6, y = 0),
                color = "black", linewidth = 2)
  }
  
  if (! is.null(color_list)) {
    
    if (lighten_negative) {
      color_list <- c(color_list,
                      lighten(color_list, amount = 0.5))
    }
    
    colors <- levels(tib_gather$fill_column)
    
    color_list <- color_list[1:length(colors)]
    names(color_list) <- colors
    #colors <- recode(colors, !!!color_list)
    
    plt <- plt +
      scale_fill_manual(values = color_list, guide = "none")
  } else {
    plt <- plt +
      scale_fill_brewer(palette = "Set1", guide = "none")
  }
  
  if (fix_yaxis) {
    plt <- plt + 
      coord_cartesian(xlim = c(max(c(plot_start,
                                     min(tib_gather$start))) / 1e6,
                               min(c(plot_end,
                                     max(tib_gather$end))) / 1e6),
                      ylim = ylimits)
  } else {
    plt <- plt + 
      coord_cartesian(xlim = c(max(c(plot_start,
                                     min(tib_gather$start))) / 1e6,
                               min(c(plot_end,
                                     max(tib_gather$end))) / 1e6))
  }
  
  if (! is.null(aspect_ratio)) {
    plt <- plt +
      theme(aspect.ratio = aspect_ratio)
  }
  
  plot(plt)
  
}
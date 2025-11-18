suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(plotgardener) # For locus plots
  library(TxDb.Mmusculus.UCSC.mm39.knownGene) 
  library(org.Mm.eg.db) # For gene ID mapping
  #library(igraph)         
  library(GenomicRanges)
  library(RColorBrewer)
})

# Get inputs: regions, scans, thresholds, and gene info
sig_regions <- read.csv("../cc_gwas/data/processed/joinLoci/trait_qtl/miQTL/all_significant_regions_summary.csv")
all_scans <- readRDS("../cc_gwas/data/processed/trait_qtl/all_scans.rds")
all_thresholds <- readRDS("../cc_gwas/data/processed/trait_qtl/all_thresholds.rds")
merged_gene_info <- read.csv("../cc_gwas/data/processed/joinLoci/geneTables/multTrait_cis-eQTL_nrvmExp.csv")

# format and filter the regions table
sig_regions <- sig_regions |>
  mutate(
    across(c(upper_pos_lod_drop, peak_pos, lower_pos_lod_drop, max_lod), as.numeric),
    chr = as.character(chr)
  ) %>%
  drop_na(upper_pos_lod_drop, lower_pos_lod_drop, chr, trait, drug)

chr3_example <- sig_regions %>%
    filter(chr == "3" & trait == "EF.21") |>
    mutate( across(c(upper_pos_lod_drop, peak_pos, lower_pos_lod_drop, max_lod), as.numeric),
        chr = as.character(chr),
        start_bp = as.integer(floor(upper_pos_lod_drop * 1e6)),
        end_bp   = as.integer(floor(lower_pos_lod_drop * 1e6))) |>
  drop_na(upper_pos_lod_drop, lower_pos_lod_drop, chr, trait, drug) |>
  dplyr::select(-upper_pos_lod_drop, -lower_pos_lod_drop)


# Subset the scan and threshold data for this locus
scan_key <- paste0(chr3_example$trait, "_", chr3_example$drug)
alt_key <- paste0(chr3_example$trait, "_Ctrl")
chr3_scan <- all_scans[[scan_key]]
chr3_threshold <- all_thresholds[[paste0(scan_key, "_threshold")]]

alt_scan <- all_scans[[alt_key]]
alt_threshold <- all_thresholds[[paste0(alt_key, "_threshold")]]

#### Set plotting parameters ####
mm39 <- assembly(
  Genome = "mm39_GRCm39",
  TxDb   = "TxDb.Mmusculus.UCSC.mm39.knownGene",
  OrgDb  = "org.Mm.eg.db"
)
# Boundaries for the plot
plot_start_bp <- chr3_example$start_bp - 5*1e5 # Add padding
plot_end_bp   <- chr3_example$end_bp + 5*1e5 # Add padding
bounds_bp <- c(plot_start_bp, plot_end_bp)

params_genome <- pgParams(
    assembly   = mm39, 
    chrom      = paste0("chr", chr3_example$chr),
    chromstart = plot_start_bp,
    chromend   = plot_end_bp
)

PLOT_DIMS <- list(page_width = 9, page_height = 6.5, res = 300)
PLOT_PARAMS <- list(x = 4.25, plot_width = 8, plot_height = 1, plot_y = 0.5)
GENE_DIMS <- list(height = 2, y_offset = 0.5, label_offset = 0.1)
MIN_YLIM <- 5
strain_colors <- rep(c("#ff0", "#888", "#f88", "#11f", "#0cf", "#0a0", "#f00", "#90e"), length.out = 8)
founders <- c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ", "NZO/H1LtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ")

dir_path <- file.path("figures")
plot_file_name <- file.path(dir_path, "fig2a.pdf")
if(!dir.exists(dir_path)){dir.create(dir_path)}

# Prepare the scan data for plotting
chr3_pg <- tibble(
    marker = names(chr3_scan$LOD),
    chr    = as.character(chr3_scan$chr),
    pos    = chr3_scan$pos$Mb * 1e6, 
    lod    = chr3_scan$LOD
    ) %>%
filter(chr == chr3_example$chr, !is.na(pos), !is.na(lod)) %>% # Filter for current chromosome
transmute(chrom = paste0("chr", chr), pos, p = 10^(-lod))

# Prepare the alternative scan data for plotting
alt_pg <- tibble(
    marker = names(alt_scan$LOD),
    chr    = as.character(alt_scan$chr),
    pos    = alt_scan$pos$Mb * 1e6, 
    lod    = alt_scan$LOD
    ) %>%
filter(chr == chr3_example$chr, !is.na(pos), !is.na(lod)) %>% # Filter for current chromosome
transmute(chrom = paste0("chr", chr), pos, p = 10^(-lod))

# Determine y-axis limits for miQTL plot
miqtl_ylim <- c(0, max(c(-log10(chr3_pg$p), alt_threshold, chr3_threshold, 5), na.rm = TRUE) + 1)

# Get list of top genes in locus 
gene_table_key <-paste0(chr3_example$trait, ":", chr3_example$drug)
genes_in_locus <- merged_gene_info[grepl(gene_table_key, merged_gene_info$trait_drug), ]
top_genes_in_locus <- genes_in_locus |>
  filter(avgenes_cpm > 5) |> 
  transmute(
    gene = mouse_gene_symbol,
    color = "#a03b60"
  )|> # make sure we label the genes mentioned in the text
  arrange(match(gene, c("Pdlim5", "Cisd2", setdiff(gene, c("Pdlim5", "Cisd2")))))

#top_genes_in_locus <- data.frame(gene= c("Pdlim5", "Cisd2"),
#                                color = "#a03b60")

############################
#### plotGardener setup ####
############################
pdf(plot_file_name, width = PLOT_DIMS$page_width, height = PLOT_DIMS$page_height)
pageCreate(width = PLOT_DIMS$page_width, height = PLOT_DIMS$page_height, default.units = "inches", showGuides = FALSE)

  ## Plot miQTL ##
  alt_miqtl_plot <- plotManhattan(
    data = alt_pg, params = params_genome,
    range = miqtl_ylim, trans = "-log10", sigVal = 10^(-alt_threshold),
    x              = PLOT_PARAMS$x,
    y              = PLOT_PARAMS$plot_y ,
    width          = PLOT_PARAMS$plot_width,
    height         = PLOT_PARAMS$plot_height,
    just = c("center", "top"), xfield = "", yfield = "",
    fill = "#c68866", sigCol = "#c76428", sigLine = FALSE, baseline = FALSE,
    default.units = "inches"
  )
  
miqtl_plot <- plotManhattan(
    data = chr3_pg, params = params_genome,
    range = miqtl_ylim, trans = "-log10", sigVal = 10^(-chr3_threshold),
    x              = PLOT_PARAMS$x,
    y              = PLOT_PARAMS$plot_y ,
    width          = PLOT_PARAMS$plot_width,
    height         = PLOT_PARAMS$plot_height,
    just = c("center", "top"), xfield = "pos", yfield = "p",
    fill = "#a1c6d1", sigCol = "#53add2", sigLine = TRUE, baseline = TRUE,
    default.units = "inches"
  )
annoYaxis(plot = miqtl_plot, at = pretty(miqtl_ylim),
            axisLine     = TRUE,
            fontsize     = 8,
            main         = FALSE)

plotText(label = "LOD",  
           x       = 2 * PLOT_PARAMS$x + 0.1,
           y       = PLOT_PARAMS$plot_y + PLOT_PARAMS$plot_height / 2,
           rot     = 270,
           fontsize= 8,
           just    = "center",
           default.units = "in")

pos_in_range_logical <- chr3_scan$pos$Mb * 10^6 >= plot_start_bp & chr3_scan$pos$Mb * 10^6 <= plot_end_bp
chr_logical <- chr3_scan$chr == "3"
marker_positions_bp <- chr3_scan$pos$Mb[pos_in_range_logical & chr_logical] * 10^6

markers <- data.frame(marker = chr3_scan$loci[pos_in_range_logical & chr_logical],
                      start  = marker_positions_bp) |> 
  arrange(start)

allele_effects_matrix <- chr3_scan$allele.effects[,markers$marker, drop = FALSE] # drop=FALSE to handle single marker case
allele_effects_transposed <- t(chr3_scan$allele.effects) |> as.data.frame()
allele_effects_transposed$marker <- rownames(allele_effects_transposed)


num_strains <- ncol(allele_effects_transposed) -1 
chromosome_name <- paste0("chr", chr3_example$chr)
founder_strains <- rownames(chr3_scan$allele.effects)

plot_data_list <- lapply(1:num_strains, function(strain){
  curr_strain <- colnames(allele_effects_transposed)[[strain]]
  temp_df <- allele_effects_transposed |> 
    dplyr::select(marker, founder_strains[strain]) %>% 
    filter(marker %in% markers$marker)
  
  temp_df <- left_join(temp_df, markers, by = "marker")
  
  temp_df <- temp_df |> 
    mutate(chrom = chromosome_name) |> 
    arrange(start)
  colnames(temp_df)[2] <- "score"
  
  # Need to change data points into ranges, so make a p value continue downstream until the next 
  temp_df$end <- c(temp_df$start[2:nrow(temp_df)] - 1L, 
                    temp_df$start[  nrow(temp_df)] + 1)   
  temp_df <- temp_df |> 
    dplyr::select(chrom, start, end, score)
  return(temp_df)
  }
)
names(plot_data_list) <- founder_strains

# Get the max abs haplotype effects to set y axis 
max_allele_effect <- plot_data_list |> rbindlist() |> pull(score) |> abs() |> max() 
max_allele_effect <- round(max_allele_effect * 1.1, digits = 2)

#strain_colors <- rep(c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"), length.out = num_strains)
signalPlots <- c()

#------- Plot the strain effects ------#
signalPlots[[1]] <- plotgardener::plotSignal( #:length(plot_data_list)
  data = plot_data_list[[1]],
  params = params_genome,
  range = c(-max_allele_effect,max_allele_effect), 
  linecolor = strain_colors[1],
  fill = NA, 
  x              = PLOT_PARAMS$x,
  y              = PLOT_PARAMS$plot_y + PLOT_PARAMS$plot_height + 0.2,
  width          = PLOT_PARAMS$plot_width,
  height         = PLOT_PARAMS$plot_height,
  just = c("center", "top"), default.units = "inches",
  baseline = TRUE, 
  baseline.color = "grey")

annoYaxis(plot = signalPlots[[1]], 
          at = c(-max_allele_effect, 0, max_allele_effect),
          axisLine     = TRUE,
          fontsize     = 8,
          main         = FALSE)
lapply(2:8, function(strain_data){
  signalPlots[[strain_data]] <- plotgardener::plotSignal( #:length(plot_data_list)
    data = plot_data_list[[strain_data]],
    params = params_genome,
    range = c(-max_allele_effect,max_allele_effect), 
    linecolor = strain_colors[strain_data],
    fill = NA, 
    x              = PLOT_PARAMS$x,
    y              = PLOT_PARAMS$plot_y + PLOT_PARAMS$plot_height + 0.2,
    width          = PLOT_PARAMS$plot_width,
    height         = PLOT_PARAMS$plot_height,
    just = c("center", "top"), default.units = "inches",
    baseline = TRUE, 
    baseline.color = "grey")
}
)

plotText(label = "Founder Effects",  
          x       = 2 * PLOT_PARAMS$x + 0.1,
          y       = PLOT_PARAMS$plot_y + 1.5*PLOT_PARAMS$plot_height + 0.2,
          rot     = 270,
          fontsize= 8,
          just    = c("center","center"),
          default.units = "in")
# --- Plot Genes ---
# This requires `genes_in_locus` to be prepared
genes_y_pos <- PLOT_PARAMS$plot_y + 2*PLOT_PARAMS$plot_height + 2*0.2 # Position below pyLMM
legendPlot <- plotLegend(
  legend = founders,
  fill = strain_colors,
  border = FALSE, 
  x = 1, y = genes_y_pos+2, width = 2, height = 1,
  fontsize = 8,
  just = c("left", "center"),
  orientation = "v",
  default.units = "inches"
)

legendPlot <- plotLegend(
  legend = c("Meet basic criteria"),
  fill = "#a03b60",
  border = FALSE, 
  x = 3, y = genes_y_pos+1.5, width = 5, height = 0.15,
  fontsize = 8,
  just = c("left", "top"),
  orientation = "v",
  default.units = "inches"
)

gene_plot <- plotGenes(
  params = params_genome,
  x = 4.25, y = genes_y_pos, width = 8, height = 1, # Adjust height as needed
  just = c("center", "top"), default.units = "inches",
  geneOrder = top_genes_in_locus$gene,
  fontsize = 8,
  geneHighlights = top_genes_in_locus,
  geneBackground = "darkgrey",)

# --- Add Genome Label ---
annoGenomeLabel(
  plot = gene_plot, # Attach to one of the plots, e.g., gene_plot
  params = params_genome,
  x = 4.25, y = genes_y_pos + 1, # Below last track
  scale = "Mb", fontsize = 10,
  just = c("center", "top"), default.units = "inches")

dev.off()  


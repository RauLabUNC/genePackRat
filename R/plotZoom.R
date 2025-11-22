#' Generate LocusZoom Plot from locusPackRat Project
#'
#' Creates a LocusZoom-style visualization plot from project data
#'
#' @param region_id Character: ID of the region to plot (NULL for first region or gene-based region)
#' @param project_dir Character: Path to locusPackRat project directory
#' @param scan_table Character: Name of supplementary table containing scan/QTL data
#' @param scan_object List: A raw scan object (e.g. scan_ctrl) to use instead of loading a table file
#' @param threshold Numeric: A specific significance threshold value (e.g. threshold_ctrl)
#' @param highlight_genes Character vector: Gene symbols to highlight in the plot
#' @param highlight_selection Character: Expression to select genes to highlight
#' @param priority_table Character: Name of supplementary table containing gene priorities/scores
#' @param signal_tables Character vector: Names of supplementary tables to plot as signal tracks
#' @param show_founders Logical: Whether to plot founder strain effects/haplotypes
#' @param founder_table Character: Name of supplementary table containing founder data
#' @param output_file Character: Path for output PDF (auto-generated if NULL)
#' @param assembly Assembly object for plotgardener (auto-detected from project if NULL)
#' @param layout Character: Page layout "landscape" (default) or "portrait"
#' @param plot_params List of plot parameters (uses defaults based on layout if NULL)
#'
#' @return Invisible TRUE on success
#'
#' @importFrom data.table fread fwrite setDT setnames as.data.table
#' @importFrom jsonlite read_json
#' @importFrom plotgardener pageCreate plotManhattan plotGenes plotRanges plotText plotSignal annoYaxis annoHighlight plotGenomeLabel
#'
#' @export
generateLocusZoomPlot <- function(
  region_id = NULL,
  project_dir = ".",
  scan_table = NULL,
  scan_object = NULL,
  threshold = NULL,
  highlight_genes = NULL,
  highlight_selection = NULL,
  priority_table = NULL,
  signal_tables = NULL,
  show_founders = FALSE,
  founder_table = "founder_haplotypes",
  output_file = NULL,
  assembly = NULL,
  layout = "landscape",
  plot_params = NULL
) {

  # Check required packages
  if (!requireNamespace("plotgardener", quietly = TRUE)) {
    stop("Package 'plotgardener' is required for plotting.")
  }

  # Define color palettes
  STRAIN_COLORS <- c(
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D", "#666666"
  )
  GENE_HIGHLIGHT_COLOR <- "#e34a33"
  GENE_BACKGROUND_COLOR <- "#fdbb84"
  LOCI_COLORS <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5")

  # Set default plot parameters
  if (is.null(plot_params)) {
    if (layout == "portrait") {
      plot_params <- list(page_width = 8.5, page_height = 11, x = 4.25, plot_width = 7.5, plot_height = 1.5, track_height = 0.8, plot_y = 1.0)
    } else {
      plot_params <- list(page_width = 10.5, page_height = 5.5, x = 4.25, plot_width = 8, plot_height = 1, track_height = 0.5, plot_y = 0.5)
    }
  }

  # --- Load project configuration ---
  packrat_dir <- file.path(project_dir, ".locusPackRat")
  if (!dir.exists(packrat_dir)) stop("No .locusPackRat directory found. Run initPackRat() first.")
  config <- jsonlite::read_json(file.path(packrat_dir, "config.json"))

  message(sprintf("Generating LocusZoom plot (%s) from %s %s %s project...", layout, config$species, config$genome, config$mode))

  # --- Determine region to plot ---
  if (config$mode == "region") {
    regions <- fread(file.path(packrat_dir, "input/regions.csv"))
    if (!is.null(region_id)) {
      locus_info <- regions[region_id == region_id]
      if (nrow(locus_info) == 0) stop("Region ID '", region_id, "' not found")
    } else {
      locus_info <- regions[1]
      region_id <- locus_info$region_id
    }
  } else {
    # Gene mode - create synthetic region
    all_genes <- fread(file.path(packrat_dir, "input/genes.csv"))
    all_genes <- all_genes[!is.na(chr) & !is.na(start) & !is.na(end)]
    if (nrow(all_genes) == 0) stop("No genes with coordinates found")
    
    chr_counts <- table(all_genes$chr)
    main_chr <- names(chr_counts)[which.max(chr_counts)]
    chr_genes <- all_genes[chr == main_chr]
    
    locus_info <- data.table(
      region_id = paste0("chr", main_chr, "_all_genes"),
      chr = main_chr,
      start = min(chr_genes$start, na.rm = TRUE),
      end = max(chr_genes$end, na.rm = TRUE)
    )
  }

  # =========================================================================
  # --- Load / Parse Scan Data ---
  # =========================================================================
  scan_data <- NULL
  threshold_val <- if(!is.null(threshold)) threshold else 5

  # 1. Adapter for Raw Object (List structure provided by user)
  if (!is.null(scan_object)) {
    message("  Using provided scan object...")
    
    # Check if it matches the List of 18 structure
    if (is.list(scan_object) && all(c("LOD", "pos", "chr") %in% names(scan_object))) {
      
      # Extract relevant vectors
      # Handle Position: List usually has Mb, we need bp
      pos_vec <- if("Mb" %in% names(scan_object$pos)) {
        scan_object$pos$Mb * 1e6 # Convert Mb to bp
      } else {
        scan_object$pos # Assume generic vector or bp
      }
      
      # Create standardized data.table
      scan_data <- data.table(
        chr = as.character(scan_object$chr),
        pos = as.numeric(pos_vec),
        lod = as.numeric(scan_object$LOD)
      )
      
      # Add p-value if present
      if ("p.value" %in% names(scan_object)) {
        scan_data[, p_value := as.numeric(scan_object$p.value)]
      } else {
        scan_data[, p_value := 10^(-lod)] # Fallback
      }
      
      # Note: allele.effects are available in this object (scan_object$allele.effects)
      
    } else {
      warning("  Provided scan_object does not match expected structure (List with $LOD, $pos$Mb, $chr). Attempting generic coercion.")
      scan_data <- as.data.table(scan_object)
    }
  } 
  # 2. Fallback to File Loading
  else if (!is.null(scan_table)) {
    scan_file <- file.path(packrat_dir, "supplementary", paste0(scan_table, ".csv"))
    if (file.exists(scan_file)) {
      scan_data <- fread(scan_file)
      message("  Loaded scan data from table: ", scan_table)
      
      # Look for threshold in file if not provided in args
      if (is.null(threshold)) {
        if ("threshold" %in% names(scan_data)) threshold_val <- unique(scan_data$threshold)[1]
        else if ("significance_threshold" %in% names(scan_data)) threshold_val <- unique(scan_data$significance_threshold)[1]
      }
    } else {
      warning("Scan table '", scan_table, "' not found")
    }
  } else {
    # 3. Auto-detect (original logic)
    available_tables <- listPackRatTables(project_dir)
    scan_tables <- available_tables[grepl("scan|qtl|lod", table_name, ignore.case = TRUE)]
    if (nrow(scan_tables) > 0) {
      scan_table <- scan_tables$table_name[1]
      message("  Auto-detected scan table: ", scan_table)
      scan_data <- fread(file.path(packrat_dir, "supplementary", paste0(scan_table, ".csv")))
      if (is.null(threshold) && "threshold" %in% names(scan_data)) threshold_val <- unique(scan_data$threshold)[1]
    }
  }
  
  # =========================================================================

  # --- Load genes in the region ---
  genes_file <- file.path(packrat_dir, "input/genes.csv")
  if (file.exists(genes_file)) {
    all_genes <- fread(genes_file)
    genes_in_locus <- all_genes[chr == locus_info$chr & end >= locus_info$start & start <= locus_info$end]
    message("  Found ", nrow(genes_in_locus), " genes in region")
  } else {
    genes_in_locus <- data.table()
  }

  # --- Determine genes to highlight (Priority / Selection / Manual) ---
  top_genes_in_locus <- data.table()
  
  if (!is.null(highlight_selection) && nrow(genes_in_locus) > 0) {
     tryCatch({
      top_genes_in_locus <- genes_in_locus[eval(parse(text = highlight_selection))]
      message(sprintf("  Highlighting %d genes based on selection: %s", nrow(top_genes_in_locus), highlight_selection))
    }, error = function(e) warning("Failed to evaluate highlight_selection: ", e$message))
  } else if (!is.null(priority_table)) {
    priority_file <- file.path(packrat_dir, "supplementary", paste0(priority_table, ".csv"))
    if (file.exists(priority_file)) {
      priority_data <- fread(priority_file)
      if ("gene_symbol" %in% names(priority_data) && nrow(genes_in_locus) > 0) {
        genes_in_locus <- merge(genes_in_locus, priority_data, by = "gene_symbol", all.x = TRUE)
        score_cols <- names(priority_data)[grepl("score|priority|rank", names(priority_data), ignore.case = TRUE)]
        if (length(score_cols) > 0) {
          setorderv(genes_in_locus, score_cols[1], order = -1)
          top_genes_in_locus <- genes_in_locus[1:min(10, .N)]
        }
      }
    }
  } else if (!is.null(highlight_genes)) {
    top_genes_in_locus <- genes_in_locus[gene_symbol %in% highlight_genes]
  } else if (nrow(genes_in_locus) > 0) {
    top_genes_in_locus <- genes_in_locus[1:min(10, .N)]
  }

  gene_highlights <- if (nrow(top_genes_in_locus) > 0) data.table(gene = top_genes_in_locus$gene_symbol, color = GENE_HIGHLIGHT_COLOR) else NULL

  # --- Overlapping regions ---
  overlapping_loci <- data.table()
  if (config$mode == "region" && exists("regions")) {
    overlapping_loci <- regions[chr == locus_info$chr & start <= locus_info$end & end >= locus_info$start & region_id != locus_info$region_id]
  }

# --- Set up genome assembly ---
  if (is.null(assembly)) {
    
    # 1. Define the requirements for each genome
    if (config$genome == "mm39") {
      req_tx <- "TxDb.Mmusculus.UCSC.mm39.knownGene"
      req_org <- "org.Mm.eg.db"
      genome_build <- "mm39"
    } else if (config$genome == "hg38") {
      req_tx <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
      req_org <- "org.Hs.eg.db"
      genome_build <- "hg38"
    } else if (config$genome == "mm10") {
      req_tx <- "TxDb.Mmusculus.UCSC.mm10.knownGene"
      req_org <- "org.Mm.eg.db"
      genome_build <- "mm10"
    } else if (config$genome == "hg19") {
      req_tx <- "TxDb.Hsapiens.UCSC.hg19.knownGene"
      req_org <- "org.Hs.eg.db"
      genome_build <- "hg19"
    } else {
      stop("Unsupported genome in config: ", config$genome)
    }

    # check required packages are installed to see gene tracts
    missing_pkgs <- character()
    if (!requireNamespace(req_tx, quietly = TRUE)) missing_pkgs <- c(missing_pkgs, req_tx)
    if (!requireNamespace(req_org, quietly = TRUE)) missing_pkgs <- c(missing_pkgs, req_org)

    if (length(missing_pkgs) > 0) {
      install_cmd <- paste0("BiocManager::install(c('", paste(missing_pkgs, collapse = "', '"), "'))")
      stop("To plot genes for ", config$genome, ", you must install the following annotation package(s): ", 
           paste(missing_pkgs, collapse = ", "),
           "\nRun this command: ", install_cmd)
    }
    # 3. Create the assembly (Now we know packages are safe)
    assembly <- plotgardener::assembly(
      Genome = genome_build,
      TxDb = req_tx,
      OrgDb = req_org
    )
  }

  # --- Output File Setup ---
  if (is.null(output_file)) {
    output_dir <- file.path(packrat_dir, "output/plots")
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    suffix <- if(layout == "portrait") "_portrait" else ""
    plot_name <- if (!is.null(region_id)) paste0("locuszoom_", region_id, suffix, ".pdf") else paste0("locuszoom_chr", locus_info$chr, suffix, ".pdf")
    output_file <- file.path(output_dir, plot_name)
  }

  # ================= PLOTTING =================
  message("Creating plot: ", output_file)
  plot_start_bp <- max(0, locus_info$start - 5e5)
  plot_end_bp <- locus_info$end + 5e5
  bounds_bp <- c(locus_info$start, locus_info$end)

  pdf(output_file, width = plot_params$page_width, height = plot_params$page_height)
  plotgardener::pageCreate(width = plot_params$page_width, height = plot_params$page_height, default.units = "inches", showGuides = FALSE)
  params_genome <- plotgardener::pgParams(assembly = assembly, chrom = paste0("chr", locus_info$chr), chromstart = plot_start_bp, chromend = plot_end_bp)
  current_y <- plot_params$plot_y

  # 1. Manhattan Plot (using prepared scan_data)
  if (!is.null(scan_data) && nrow(scan_data) > 0) {
    plot_data <- copy(scan_data)
    if ("chromosome" %in% names(plot_data)) setnames(plot_data, "chromosome", "chr")
    if ("position" %in% names(plot_data)) setnames(plot_data, "position", "pos")
    
    # Ensure chromosome format matches params_genome
    if (!grepl("^chr", as.character(plot_data$chr[1]))) plot_data[, chrom := paste0("chr", plot_data$chr)] else plot_data[, chrom := chr]
    plot_data <- plot_data[chrom == paste0("chr", locus_info$chr)]

    if ("lod" %in% names(plot_data)) plot_data[, p := 10^(-lod)] else if ("p_value" %in% names(plot_data)) setnames(plot_data, "p_value", "p")
    
    # Plotting
    ylim <- c(0, max(-log10(plot_data$p), -log10(10^(-threshold_val)), 5, na.rm = TRUE) + 1)
    miqtl_plot <- plotgardener::plotManhattan(
      data = plot_data, params = params_genome, range = ylim, trans = "-log10",
      sigVal = 10^(-threshold_val), x = plot_params$x, y = current_y,
      width = plot_params$plot_width, height = plot_params$plot_height, just = c("center", "top"),
      fill = "#a6cee3", sigCol = "#1f78b4", sigLine = TRUE, baseline = TRUE, default.units = "inches"
    )
    plotgardener::annoYaxis(plot = miqtl_plot, at = pretty(ylim), axisLine = TRUE, fontsize = 8)
    plotgardener::plotText(label = "LOD Score", x = plot_params$x - plot_params$plot_width/2 - 0.3, y = current_y + plot_params$plot_height / 2, rot = 90, fontsize = 8, just = "center", default.units = "inches")
    
    # Highlight region
    plotgardener::annoHighlight(
      plot = miqtl_plot, chrom = paste0("chr", locus_info$chr),
      chromstart = floor(min(bounds_bp)), chromend = ceiling(max(bounds_bp)),
      fill = "#fb9a99", y = current_y, height = plot_params$plot_height,
      just = c("left", "top"), default.units = "inches", alpha = 0.2, params = params_genome
    )
    current_y <- current_y + plot_params$plot_height + 0.2
  }

  # 2. Signal Tracks
  if (!is.null(signal_tables)) {
    for (sig_name in signal_tables) {
      sig_file <- file.path(packrat_dir, "supplementary", paste0(sig_name, ".csv"))
      if (file.exists(sig_file)) {
        sig_data <- fread(sig_file)
        if (!"chrom" %in% names(sig_data) && "chr" %in% names(sig_data)) setnames(sig_data, "chr", "chrom")
        if (!grepl("^chr", as.character(sig_data$chrom[1]))) sig_data[, chrom := paste0("chr", chrom)]
        score_col <- intersect(names(sig_data), c("score", "value", "lod", "intensity"))[1]
        if (!is.na(score_col)) setnames(sig_data, score_col, "score")
        
        if ("score" %in% names(sig_data)) {
          message("  Adding signal track: ", sig_name)
          plotgardener::plotSignal(
            data = sig_data, params = params_genome, x = plot_params$x, y = current_y,
            width = plot_params$plot_width, height = plot_params$track_height, just = c("center", "top"),
            linecolor = "#377eb8", fill = "#377eb8", default.units = "inches"
          )
          plotgardener::plotText(label = sig_name, x = plot_params$x - plot_params$plot_width/2 - 0.3, y = current_y + plot_params$track_height / 2, rot = 90, fontsize = 8, just = "center", default.units = "inches")
          current_y <- current_y + plot_params$track_height + 0.1
        }
      }
    }
  }

  # 3. Founder Effects
  if (show_founders) {
    founder_file <- file.path(packrat_dir, "supplementary", paste0(founder_table, ".csv"))
    if (file.exists(founder_file)) {
      founder_data <- fread(founder_file)
      strain_col <- intersect(names(founder_data), c("strain", "founder", "genotype"))[1]
      if (!is.na(strain_col)) {
        message("  Adding founder track: ", founder_table)
        setnames(founder_data, strain_col, "strain_id")
        if ("chr" %in% names(founder_data)) setnames(founder_data, "chr", "chrom")
        if (!grepl("^chr", as.character(founder_data$chrom[1]))) founder_data[, chrom := paste0("chr", chrom)]
        
        strain_names <- c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ", "NZO/HILtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ")
        mapped_colors <- STRAIN_COLORS; names(mapped_colors) <- strain_names
        founder_data[, color := "#999999"]
        for (sn in strain_names) founder_data[grepl(sn, strain_id, fixed = TRUE), color := mapped_colors[sn]]
        
        plotgardener::plotRanges(
          data = founder_data, params = params_genome, fill = founder_data$color, linecolor = NA,
          x = plot_params$x, y = current_y, width = plot_params$plot_width, height = plot_params$track_height, just = c("center", "top"), default.units = "inches"
        )
        plotgardener::plotText(label = "Founders", x = plot_params$x - plot_params$plot_width/2 - 0.3, y = current_y + plot_params$track_height / 2, rot = 90, fontsize = 8, just = "center", default.units = "inches")
        current_y <- current_y + plot_params$track_height + 0.2
      }
    }
  }

  # 4. Overlapping Regions
  if (nrow(overlapping_loci) > 1) {
    loci_palette <- function(n) if (n <= length(LOCI_COLORS)) LOCI_COLORS[1:n] else grDevices::colorRampPalette(LOCI_COLORS)(n)
    plotgardener::plotRanges(
      data = overlapping_loci, params = params_genome, fill = plotgardener::colorby("region_id", palette = loci_palette),
      x = plot_params$x, y = current_y, width = plot_params$plot_width, height = 0.5, just = c("center", "top"), default.units = "inches"
    )
    plotgardener::plotText(label = "Overlapping Regions", x = plot_params$x - plot_params$plot_width/2 - 0.3, y = current_y + 0.25, rot = 90, fontsize = 8, just = "center", default.units = "inches")
    current_y <- current_y + 0.7
  }

  # 5. Genes
  gene_order <- if(nrow(top_genes_in_locus) > 0) top_genes_in_locus$gene_symbol else NULL
  tryCatch({
    plotgardener::plotGenes(
      params = params_genome, x = plot_params$x, y = current_y, width = plot_params$plot_width, height = 1,
      just = c("center", "top"), default.units = "inches", geneOrder = gene_order, fontsize = 6,
      geneHighlights = gene_highlights, geneBackground = GENE_BACKGROUND_COLOR
    )
  }, error = function(e) {
    message("Note: Gene highlighting failed, plotting without highlights")
    plotgardener::plotGenes(params = params_genome, x = plot_params$x, y = current_y, width = plot_params$plot_width, height = 1, just = c("center", "top"), default.units = "inches", fontsize = 6)
  })

  plotgardener::plotGenomeLabel(params = params_genome, x = plot_params$x, y = current_y + 1.1, length = plot_params$plot_width, just = c("center", "top"), default.units = "inches")
  
  title_text <- if (!is.null(region_id)) paste0("LocusZoom Plot: ", region_id) else paste0("LocusZoom Plot: Chr", locus_info$chr, ":", round(locus_info$start/1e6, 1), "-", round(locus_info$end/1e6, 1), " Mb")
  plotgardener::plotText(label = title_text, x = plot_params$page_width / 2, y = 0.3, fontsize = 12, fontface = "bold", just = "center", default.units = "inches")
  plotgardener::plotText(label = paste0(config$species, " ", config$genome, " | ", config$mode, " mode"), x = plot_params$page_width / 2, y = 0.5, fontsize = 10, just = "center", default.units = "inches")

  dev.off()
  message("Plot saved to: ", output_file)
  invisible(TRUE)
}


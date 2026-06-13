main <- function() {
  tryCatch({
    cat("\014")
    rm(list = ls())
    while (!is.null(dev.list())) dev.off()
    print("Starting script execution.")
    
    library(DESeq2)
    library(tcltk)
    
    intro_message <- function() {
      tt <- tktoplevel()
      tkwm.title(tt, "RNA-Seq Analysis Software")
      tcl("wm", "attributes", tt, topmost = TRUE)
      msg <- "DESeq2 Analysis:
        To use this software, you need to have the RNA-Seq count data (post-alignment) ready in a CSV file where the rows are genes and columns are samples.
      
      The first column should be titled 'Gene Symbol' and contain gene names. The next columns should represent samples of different groups. For example, if you have 3 groups (A, B, and C) where:
        - Group A has 10 samples
      - Group B has 15 samples
      - Group C has 10 samples
      
      The columns should be organized and named as follows:
        A (1), A (2), ..., A (10), B (1), ..., B (15), C (1), ..., C (10).
      
      Note: The software will recognize samples within a group by the unique name before parentheses (e.g., 'A', 'B', 'C'). The numbers inside parentheses indicate the sample count for each group.
      
      Please make sure your CSV file is formatted accordingly."
      tkgrid(tklabel(tt, text = msg, wraplength = 600, justify = "left"), padx = 20, pady = 20)
      
      user_choice <- tclVar("continue")
      button_frame <- tkframe(tt)
      tkgrid(tkbutton(button_frame, text = "Continue", command = function() {
        tclvalue(user_choice) <- "continue"
        tkdestroy(tt)
      }), padx = 10, pady = 10)
      tkgrid(tkbutton(button_frame, text = "Exit and Fix", command = function() {
        tclvalue(user_choice) <- "exit"
        tkdestroy(tt)
      }), padx = 10, pady = 10)
      tkgrid(button_frame)
      tkwait.window(tt)
      if (tclvalue(user_choice) == "exit") return(FALSE)
      return(TRUE)
    }
    
    if (!intro_message()) stop("Execution stopped by user.")
    
    load_data <- function(filepath) {
      data <- read.csv(filepath, header = TRUE, check.names = FALSE)
      gene_symbols <- data[, 1]
      counts <- data[, -1]
      
      # Order by row sums to keep the most expressed duplicates
      row_sums <- rowSums(counts)
      data_ordered <- data[order(row_sums, decreasing = TRUE), ]
      
      # Remove duplicate gene symbols (keep the first occurrence)
      data_unique <- data_ordered[!duplicated(data_ordered[, 1]), ]
      
      # Extract updated counts and gene symbols
      counts <- data_unique[, -1]
      gene_symbols <- data_unique[, 1]
      
      rownames(counts) <- gene_symbols
      return(counts)
    }
    
    
    detect_groups <- function(column_names) {
      groups <- sub(" \\(.+\\)", "", column_names)
      return(groups)
    }
    
    confirm_groups <- function(column_names) {
      groups <- sub(" \\(.+\\)", "", column_names)
      group_counts <- table(groups)
      
      tt <- tktoplevel()
      tkwm.title(tt, "Confirm Groups")
      tkgrid(tklabel(tt, text = "Detected Groups and Sample Counts:"))
      
      for (g in names(group_counts)) {
        label_text <- sprintf("%s: %d samples", g, group_counts[[g]])
        tkgrid(tklabel(tt, text = label_text, justify = "left"), padx = 10, sticky = "w")
      }
      
      decision <- tclVar("no")
      button_frame <- tkframe(tt)
      
      tkgrid(tkbutton(button_frame, text = "Confirm", command = function() {
        tclvalue(decision) <- "yes"
        tkdestroy(tt)
      }), padx = 10)
      
      tkgrid(tkbutton(button_frame, text = "Decline", command = function() {
        tclvalue(decision) <- "no"
        tkdestroy(tt)
      }), padx = 10)
      
      tkgrid(button_frame, pady = 10)
      tkwait.window(tt)
      
      return(tclvalue(decision) == "yes")
    }
    
    
    select_groups <- function(groups) {
      unique_groups <- unique(groups)
      tt <- tktoplevel()
      tkwm.title(tt, "Select Baseline and Treatment")
      
      baseline <- tclVar(unique_groups[1])
      treatment <- tclVar(unique_groups[2])
      tkgrid(tklabel(tt, text = "Select Baseline:"))
      tkgrid(ttkcombobox(tt, values = unique_groups, textvariable = baseline))
      tkgrid(tklabel(tt, text = "Select Treatment:"))
      tkgrid(ttkcombobox(tt, values = unique_groups, textvariable = treatment))
      tkgrid(tkbutton(tt, text = "Confirm", command = function() {
        if (tclvalue(baseline) == tclvalue(treatment)) {
          tkmessageBox(message = "Baseline and Treatment must differ.")
        } else {
          tkdestroy(tt)
        }
      }))
      tkwait.window(tt)
      return(list(baseline = tclvalue(baseline), treated = tclvalue(treatment)))
    }
    
    run_deseq2 <- function(counts, groups, baseline, treated) {
      group_factor <- factor(groups)
      design_df <- data.frame(row.names = colnames(counts), group = group_factor)
      
      # Create DESeqDataSet
      dds <- DESeqDataSetFromMatrix(countData = counts, colData = design_df, design = ~group)
      
      # Relevel the group factor to set baseline as the reference level
      dds$group <- relevel(dds$group, ref = baseline)
      
      # Run DESeq2
      dds <- DESeq(dds)
      
      # Get the results for the specified contrast (treated vs baseline)
      res <- results(dds, contrast = c("group", treated, baseline))
      
      # Order results by p-value
      resOrdered <- res[order(res$pvalue), ]
      
      # Convert results to a dataframe
      resDF <- as.data.frame(resOrdered)
      
      # Add GeneSymbol column
      resDF$GeneSymbol <- rownames(resDF)
      
      # Rename columns for compatibility
      colnames(resDF)[colnames(resDF) == "log2FoldChange"] <- "logFC"
      colnames(resDF)[colnames(resDF) == "pvalue"] <- "PValue"
      colnames(resDF)[colnames(resDF) == "padj"] <- "FDR"
      
      # Apply Bonferroni correction to p-values
      resDF$Bonferroni_pvalue <- p.adjust(resDF$PValue, method = "bonferroni")
      
      # Move GeneSymbol to the first column
      resDF <- resDF[, c("GeneSymbol", setdiff(names(resDF), "GeneSymbol"))]
      
      # Save the full table with both FDR and Bonferroni p-values
      write.csv(resDF, file = paste0("DESeq2_", baseline, "_vs_", treated, "_with_FDR_and_Bonferroni.csv"), row.names = FALSE)
    }
    
    
    ask_comparisons <- function() {
      tt <- tktoplevel()
      tkwm.title(tt, "Number of Comparisons")
      tkgrid(tklabel(tt, text = "How many comparisons?"))
      num_comparisons <- tclVar("1")
      tkgrid(tkentry(tt, textvariable = num_comparisons))
      tkgrid(tkbutton(tt, text = "Confirm", command = function() tkdestroy(tt)))
      tkwait.window(tt)
      return(as.integer(tclvalue(num_comparisons)))
    }
    
    filepath <- file.choose()
    counts <- load_data(filepath)
    column_names <- colnames(counts)
    groups <- detect_groups(column_names)
    if (!confirm_groups(column_names)) stop("Group confirmation declined.")
    
    
    num_comp <- ask_comparisons()
    for (i in 1:num_comp) {
      selected <- select_groups(groups)
      run_deseq2(counts, groups, selected$baseline, selected$treated)
    }
    
    print("Script completed successfully.")
  }, error = function(e) {
    print(sprintf("An error occurred: %s", e$message))
  })
}

main()






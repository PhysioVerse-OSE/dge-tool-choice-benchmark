main <- function() {
  tryCatch({
    print("Starting script execution.")




    
    cat("\014")
    # Clear the environment
    rm(list = ls())
    while (!is.null(dev.list())) dev.off()
    print("Script has started running.")
    
    library(edgeR)
    library(tcltk)
    library(rstudioapi)
    
    # Set working directory to the script's directory
   # setwd(dirname(getActiveDocumentContext()$path))
    
    # Verify the current working directory
   # print(getwd())
    
    
    print("Guideline")
    
    ###################### Step 0: Guideline #################
    # Function to display an introductory message
    intro_message <- function() {
      require(tcltk)
      
      # Create the main window
      tt <- tktoplevel()
      tkwm.title(tt, "RNA-Seq Analysis Software")
      tcl("wm", "attributes", tt, topmost = TRUE)
      # Add the message to the window
      msg <- "GLMQL-MAS:
      
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
      
      # Variable to track the user's choice
      user_choice <- tclVar("continue")
      
      # Create the buttons
      button_frame <- tkframe(tt)
      
      # Continue Button
      continue_btn <- tkbutton(button_frame, text = "Continue", command = function() {
        tclvalue(user_choice) <- "continue"  # Set choice to continue
        tkdestroy(tt)  # Close the window
      })
      
      # Exit and Fix Button
      exit_btn <- tkbutton(button_frame, text = "Exit and Fix", command = function() {
        tclvalue(user_choice) <- "exit"  # Set choice to exit
        tkdestroy(tt)  # Close the window
      })
      
      # Arrange buttons in the window
      tkgrid(continue_btn, padx = 10, pady = 10)
      tkgrid(exit_btn, padx = 10, pady = 10)
      tkgrid(button_frame)
      
      # Wait for the user to click a button
      tkwait.window(tt)
      
      # Check user choice
      if (tclvalue(user_choice) == "exit") {
        cat("Please fix the input CSV file and rerun the script.\n")
        return(FALSE)  # Exit the function and signal to stop further execution
      }
      
      return(TRUE)  # Signal to continue with the script
    }
    
    # Call the intro_message function at the start of your script
    if (!intro_message()) {
      stop("Execution stopped by user. Fix the input CSV file and rerun the script.")
    }
    
    # Add further code execution here (only if the user clicks "Continue")
    cat("Script continues...\n")
    
    
    
    
    
    
    ############################ Step 1: Uploading the csv file and creating DGEList ############
    
    # Function to load data and remove duplicate genes
    load_data <- function(filepath) {
      # Load the CSV file
      data <- read.csv(filepath, header = TRUE, check.names = FALSE)
      
      # Extract gene symbols and counts
      gene_symbols <- data[, 1]  # First column as gene symbols
      counts <- data[, -1]       # All other columns as counts
      
      # Create a DGEList object
      dge <- DGEList(counts = as.matrix(counts), genes = data.frame(GeneSymbol = gene_symbols))
      
      
      # Order by total count (row sums) and remove duplicates
      o <- order(rowSums(dge$counts), decreasing = TRUE)
      dge <- dge[o, ]
      d <- duplicated(dge$genes$GeneSymbol)
      dge <- dge[!d, ]
      
      # Set row names to unique gene symbols
      rownames(dge$counts) <- dge$genes$GeneSymbol
      
      return(dge)
    }
    
    # Function to automatically detect groups from column names
    detect_groups <- function(column_names) {
      # Extract the group names (before parentheses)
      groups <- sub(" \\(.*\\)", "", column_names)
      
      # Count the number of replicates for each group
      group_info <- table(groups)
      return(as.list(group_info))
    }
    
    # Function to confirm detected groups
    confirm_groups <- function(group_info) {
      # Create a GUI window
      tt <- tktoplevel()
      tkwm.title(tt, "Confirm Groups")
      
      # Display group information for user confirmation
      tkgrid(tklabel(tt, text = "Detected Groups and Sample Counts:"))
      for (group_name in names(group_info)) {
        tkgrid(tklabel(tt, text = sprintf("%s: %d samples", group_name, group_info[[group_name]])))
      }
      
      # Variable to capture user decision
      user_decision <- tclVar("no")
      
      # Buttons for confirmation or decline
      tkgrid(tkbutton(tt, text = "Confirm", command = function() {
        tclvalue(user_decision) <- "yes"
        tkdestroy(tt)
      }))
      tkgrid(tkbutton(tt, text = "Decline", command = function() {
        tclvalue(user_decision) <- "no"
        tkdestroy(tt)
      }))
      
      # Wait for the user's response
      tkwait.window(tt)
      
      return(tclvalue(user_decision) == "yes")
    }
    
    # Main function to run the analysis
    run_analysis <- function() {
      # Step 1: Upload CSV file
      filepath <- file.choose()
      dge <- load_data(filepath)
      
      # Step 2: Detect groups automatically
      column_names <- colnames(dge$counts)
      group_info <- detect_groups(column_names)
      
      # Step 3: Confirm detected groups with the user
      if (confirm_groups(group_info)) {
        # Assign group labels to the columns
        group_labels <- rep(NA, ncol(dge$counts))
        for (group_name in names(group_info)) {
          # Get the indices of columns belonging to this group
          group_indices <- which(sub(" \\(.*\\)", "", column_names) == group_name)
          group_labels[group_indices] <- group_name
        }
        
        # Step 4: Assign the group labels to the DGEList
        dge <- DGEList(counts = dge$counts, group = factor(group_labels))
        
        
        
        
        # Apply filtering using filterByExpr
#        keep <- filterByExpr(dge$counts, group = dge$samples$group)
 #       dge <- dge[keep, , keep.lib.sizes = FALSE]
#        
        # Print summary of the DGEList
        print(dge)
        
        return(dge)
      } else {
        tkmessageBox(title = "Error", message = "Please fix the CSV file and try again.")
        return(NULL)
      }
    }
    
    # Execute the main function
    dge <- run_analysis()
    
    table(dge$samples$group)
    
    
    ##################### Step 2: TMM Normalization ########
    
    # Perform TMM normalization
    dge <- calcNormFactors(dge)

    # Calculate normalized counts
    norm_factors <- dge$samples$norm.factors  # Retrieve normalization factors
    lib_sizes <- dge$samples$lib.size         # Retrieve library sizes
    
    # Calculate effective library sizes
    effective_lib_sizes <- lib_sizes * norm_factors
    
    # Compute TMM normalized counts (counts per million - CPM)
    tmm_normalized_counts <- cpm(dge, normalized.lib.sizes = TRUE, prior.count = 0, log = FALSE)
    
    # Save the TMM normalized counts as a CSV file
    write.csv(tmm_normalized_counts, file = "TMM_normalized_counts.csv", row.names = TRUE)
    
    # Optional: Print a message or preview the data
    cat("TMM-normalized counts saved to 'TMM_normalized_counts.csv'.")
    

    
    ################ Step 3: Baseline and Treated Groups should be defined  #######################
    # Function to select baseline and treated groups and print their sizes
    select_groups <- function(dge) {
      require(tcltk)  # Ensure the tcltk library is loaded
      group_names <- levels(dge$samples$group)
      
      # GUI setup
      tt <- tktoplevel()
      tkwm.title(tt, "Select Baseline and Treatment Groups")
      
      # Set default selection directly in tclVar initialization
      baseline <- tclVar(group_names[1])
      treatment <- tclVar(if (length(group_names) > 1) group_names[2] else group_names[1])
      
      # Adding widgets to the window
      tkgrid(tklabel(tt, text = "Select Baseline Group:"))
      baseline_menu <- ttkcombobox(tt, values = group_names, textvariable = baseline, state = "readonly")
      tkgrid(baseline_menu)
      
      tkgrid(tklabel(tt, text = "Select Treatment Group:"))
      treatment_menu <- ttkcombobox(tt, values = group_names, textvariable = treatment, state = "readonly")
      tkgrid(treatment_menu)
      
      # Confirmation button
      confirm_btn <- tkbutton(tt, text = "Confirm", command = function() {
        if (tclvalue(baseline) == tclvalue(treatment)) {
          tkmessageBox(title = "Error", message = "Baseline and Treatment groups must be different.")
        } else {
          tkdestroy(tt)
        }
      })
      tkgrid(confirm_btn)
      
      tkwait.window(tt)
      
      # Extracting the sizes of the selected groups
      baseline_group_size <- sum(dge$samples$group == tclvalue(baseline))
      treatment_group_size <- sum(dge$samples$group == tclvalue(treatment))
      
      print(paste("Size of Baseline group (", tclvalue(baseline), "): ", baseline_group_size, sep = ""))
      print(paste("Size of Treated group (", tclvalue(treatment), "): ", treatment_group_size, sep = ""))
      
      return(list(baseline = tclvalue(baseline), treated = tclvalue(treatment)))
    }
    
    ################ Step 4: GLM Quasi-Likelihood F-test #######################
    
    # Function to perform GLM with Quasi-Likelihood F-test
    perform_glm_qlf <- function(dge, selected_groups) {
      # Subset the DGEList object for the selected groups
      dge_subset <- dge[, dge$samples$group %in% c(selected_groups$baseline, selected_groups$treated)]
      
      # Relevel the group factor to set the baseline as the reference level
      dge_subset$samples$group <- relevel(dge_subset$samples$group, ref = selected_groups$baseline)
      
      # Design matrix for GLM
      design <- model.matrix(~ group, data = dge_subset$samples)
      
      # Estimate dispersions
      dge_subset <- estimateDisp(dge_subset, design)
      
      # Fit the GLM model with Quasi-Likelihood
      fit <- glmQLFit(dge_subset, design)
      
      # Perform the Quasi-Likelihood F-test
      qlf <- glmQLFTest(fit, coef = 2)  # coef = 2 compares the treated group against the baseline group
      
      # View the top differentially expressed genes
      top_results <- topTags(qlf, n = Inf)
      
      # Print the results to the console
      print(head(top_results$table))
      top_results$table$Bonferroni_pvalue <- p.adjust(top_results$table$PValue, method = "bonferroni")
      
      
      # Save the results to a CSV file
      write.csv(top_results$table, file = "differential_expression_results.csv", row.names = TRUE)
      
      # Optional: Plot an MA plot to visualize differential expression
      # Correct MA Plot
      plotMD(qlf, main = "MA Plot", xlab = "Average Log CPM", ylab = "Log Fold Change")
      abline(h = c(-1, 1), col = "blue")
      return(top_results)
    }
    
    # Function to ask the user how many comparisons they want to perform
    ask_comparisons <- function() {
      require(tcltk)
      
      # GUI to ask for the number of comparisons
      tt <- tktoplevel()
      tkwm.title(tt, "Number of Comparisons")
      
      tkgrid(tklabel(tt, text = "How many comparisons would you like to perform?"))
      
      num_comparisons <- tclVar("1")  # Default value
      
      entry <- tkentry(tt, textvariable = num_comparisons, width = 10)
      tkgrid(entry)
      
      confirm_btn <- tkbutton(tt, text = "Confirm", command = function() {
        tkdestroy(tt)
      })
      tkgrid(confirm_btn)
      
      tkwait.window(tt)
      
      return(as.integer(tclvalue(num_comparisons)))
    }
    
    ############################## Main Execution #################################
    
    # Ask user for the number of comparisons
    num_comparisons <- ask_comparisons()
    for (i in 1:num_comparisons) {
      cat("\nPerforming comparison", i, "of", num_comparisons, "\n")
      
      # Call the `select_groups` function to allow the user to select baseline and treated groups
      selected_groups <- select_groups(dge)  # Select baseline and treated groups
      results <- perform_glm_qlf(dge, selected_groups)
      
      #### Save Results ######
      
      # Generate a meaningful filename based on the selected groups
      filename <- paste0("EdgeR_", selected_groups$baseline, "_vs_", selected_groups$treated, ".csv")
      
      # Extract the results table
      results_table <- results$table
      
      # Add the row names (Gene Symbol) as a column and include a header
      results_table <- data.frame(GeneSymbol = rownames(results_table), results_table)
      
      # Write the results to a CSV file
      write.csv(results_table, file = filename, row.names = FALSE)
      
      # Print a message indicating the file was saved
      cat("Differential expression results saved to:", filename, "\n")
    }
    
    
    
    
    print("Script completed successfully.")
    
  }, error = function(e) {
    print(sprintf("An error occurred: %s", e$message))
  })
}

# Call the main function to execute
main()

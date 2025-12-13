# Author: Sorro Sun
# Date: 2025-11-29
# Description: Shiny App for Protein Secondary Structure Prediction
# Dependencies: Chou-Fasman, Improved Chou-Fasman, GOR
# Modes: Single Prediction (CF, ICF, GOR), Comparison, Benchmark (Multi-Data)

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# =========================================================================
# 1. Parsing Functions
# =========================================================================

# --- Parse interval output from Chou-Fasman ---
parse_cf_output <- function(output_lines, seq_len) {
  alpha_lines <- grep("Alpha.*:", output_lines, value = TRUE)
  beta_lines <- grep("Beta.*:", output_lines, value = TRUE)
  
  extract_indices <- function(line) {
    if (length(line) == 0) return(numeric(0))
    matches <- str_match_all(line, "\\{(\\d+)\\s+(\\d+)\\}")[[1]]
    if (nrow(matches) == 0) return(numeric(0))
    
    idx <- c()
    starts <- as.numeric(matches[, 2])
    ends <- as.numeric(matches[, 3])
    
    for (i in seq_along(starts)) {
      # Adapt Go (0-based) to R (1-based)
      s <- max(1, starts[i] + 1)
      e <- min(seq_len, ends[i] + 1)
      if (s <= e) idx <- c(idx, s:e)
    }
    return(idx)
  }
  
  h_indices <- if(length(alpha_lines) > 0) extract_indices(tail(alpha_lines, 1)) else c()
  e_indices <- if(length(beta_lines) > 0) extract_indices(tail(beta_lines, 1)) else c()
  
  res_vec <- rep("C", seq_len)
  res_vec[h_indices] <- "H"
  res_vec[e_indices] <- "E" 
  
  return(paste0(res_vec, collapse = ""))
}

# --- Parse string output (Improved CF & GOR) ---
parse_string_output <- function(output_lines, seq_len) {
  marker_idx <- grep("Prediction Result:", output_lines, fixed = TRUE)
  candidate <- NULL
  
  if (length(marker_idx) > 0) {
    possible_line_idx <- marker_idx[length(marker_idx)] + 1 
    if (possible_line_idx <= length(output_lines)) {
      candidate <- trimws(output_lines[possible_line_idx])
    }
  }
  
  if (is.null(candidate) || nchar(gsub("-", "C", candidate)) != seq_len) {
    valid_lines <- grep("^[HECX-]+$", trimws(output_lines), value = TRUE)
    if (length(valid_lines) > 0) {
      lens <- nchar(valid_lines)
      best_idx <- which.min(abs(lens - seq_len))
      candidate <- valid_lines[best_idx]
    }
  }
  
  if (!is.null(candidate)) {
    candidate <- gsub("-", "C", candidate)
    candidate <- gsub(" ", "", candidate)
    return(candidate)
  }
  return(NULL)
}

# Helper to calculate accuracy (reused in multiple places)
# Note: This function specifically excludes 'X' from accuracy calculation
calc_accuracy_score <- function(p, e) {
  if (is.null(e)) return(0) # No experiment data
  len <- min(nchar(p), nchar(e))
  if(len == 0) return(0)
  
  p_vec <- strsplit(substr(p, 1, len), "")[[1]]
  e_vec <- strsplit(substr(e, 1, len), "")[[1]]
  
  # Only compare positions where experimental structure is NOT 'X'
  idx <- which(e_vec != "X")
  
  if (length(idx)==0) return(0)
  mean(p_vec[idx] == e_vec[idx])
}

# =========================================================================
# 2. UI Definition
# =========================================================================
ui <- fluidPage(
  titlePanel("Protein Secondary Structure Prediction"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Upload .txt file as input", accept = c("text/plain", ".txt")),
      
      # --- UPDATED HELP TEXT SECTION ---
      helpText("Input types:"),
      
      helpText("1. Amino Acid Sequence (1 line)"),
      helpText("-> Modes 1-3 available"),
      
      helpText("2. Amino Acid Sequence + Structure from Experiment (2 lines)"),
      helpText("-> All Modes available"),
      
      helpText("3. Multiple pairs of Name + Sequence + Structure "),
      helpText("-> Mode 5 available"),
      # ---------------------------------
      
      radioButtons("mode", "Select Mode:",
                   choices = list(
                     "1. Chou-Fasman" = "cf",
                     "2. Improved CF" = "icf",
                     "3. GOR" = "gor",
                     "4. Comparison" = "comp",
                     "5. Benchmark" = "bench" 
                   )),
      hr(),
      
      # Legend: Only show for single prediction modes (1-3)
      conditionalPanel(
        condition = "input.mode != 'comp' && input.mode != 'bench'",
        h4("Legend"),
        tags$ul(
          tags$li(span("Helix (H)", style = "color: #E69191; font-weight: bold;")), 
          # Removed background-color: #555 and other box styles
          tags$li(span("Sheet (E)", style = "color: #F9DF91; font-weight: bold;")), 
          tags$li(span("Loop/Coil (C)", style = "color: #A8CFE8; font-weight: bold;")), 
          tags$li(span("Disordered (X)", style = "color: #808080; font-weight: bold;")) 
        )
      ),
      
      # Note for Benchmark
      conditionalPanel(
        condition = "input.mode == 'bench'",
        helpText("Note: This module is only effective with multi-data input.")
      ),
      
      hr(),
      h4("Environment Check"),
      verbatimTextOutput("env_check"),
      h4("Execution Log"),
      verbatimTextOutput("debug_log")
    ),
    
    mainPanel(
      uiOutput("structurePlotUI"),
      h4("Data Info"),
      verbatimTextOutput("seq_info")
    )
  )
)

# =========================================================================
# 3. Server Logic
# =========================================================================
server <- function(input, output, session) {
  
  logs <- reactiveVal("System Ready.")
  
  add_log <- function(msg) {
    timestamp <- format(Sys.time(), "%H:%M:%S")
    full_msg <- paste0("[", timestamp, "] ", msg)
    print(full_msg) 
    current_logs <- isolate(logs())
    logs(paste(current_logs, full_msg, sep = "\n"))
  }
  
  find_executable <- function(base_folder, sub_folders = c(""), candidates) {
    wd <- normalizePath(getwd(), winslash = "/")
    if (.Platform$OS.type == "windows") candidates <- paste0(candidates, ".exe")
    for (sub in sub_folders) {
      for (name in candidates) {
        path <- if (sub == "") file.path(wd, base_folder, name) else file.path(wd, base_folder, sub, name)
        if (file.exists(path)) return(path)
      }
    }
    return(NULL)
  }
  
  output$env_check <- renderText({
    cf_path <- find_executable("chou_fasman", c(""), c("chou_fasman", "main"))
    icf_path <- find_executable("improved_chou_fasman", c(""), c("improved_chou_fasman", "main"))
    gor_path <- find_executable("GOR", c("", "predict"), c("gor_predict", "predict", "GOR", "main"))
    gor_dir <- file.path(getwd(), "GOR")
    gor_model_path <- file.path(gor_dir, "gor_model.json")
    
    status_msg <- paste(
      paste("CF:", if(!is.null(cf_path)) "FOUND" else "MISSING"),
      paste("ICF:", if(!is.null(icf_path)) "FOUND" else "MISSING"),
      paste("GOR Exe:", if(!is.null(gor_path)) "FOUND" else "MISSING"),
      paste("GOR Model:", if(file.exists(gor_model_path)) "FOUND" else "MISSING"),
      sep = "\n"
    )
    return(status_msg)
  })
  
  # --- 1. Read Data (Smart Parsing for 3 formats) ---
  dataInput <- reactive({
    req(input$file1)
    add_log(paste("Reading file:", input$file1$name))
    
    tryCatch({
      lines <- readLines(input$file1$datapath, warn = FALSE)
      # Clean up lines
      lines <- lines[trimws(lines) != ""]
      
      N <- length(lines)
      if (N == 0) {
        add_log("Error: File content is empty.")
        return(NULL)
      }
      
      pairs <- list()
      
      # --- Smart Format Detection ---
      is_3_line_format <- FALSE
      
      # Heuristic: Check if divisible by 3. 
      # Also check if line 1 looks like a name (short) and line 2 looks like sequence (long, amino acids)
      if (N %% 3 == 0) {
        line1 <- trimws(lines[1]) # Potential ID
        line2 <- trimws(lines[2]) # Potential Seq
        # If ID is shorter than Seq OR ID contains numbers (common in PDB IDs like 5V4TB)
        if (nchar(line1) < nchar(line2) || grepl("[0-9]", line1)) {
          is_3_line_format <- TRUE
        }
      }
      
      if (is_3_line_format) {
        # Format: Name / Sequence / Structure
        num_sets <- N / 3
        add_log(paste("Detected 3-line format (Name+Seq+Str). Found", num_sets, "entries."))
        for (i in 1:num_sets) {
          pairs[[i]] <- list(
            id = trimws(lines[3*i - 2]),  # Line 1: Name
            seq = trimws(lines[3*i - 1]), # Line 2: Sequence
            exp = trimws(lines[3*i])      # Line 3: Structure
          )
        }
        
      } else if (N %% 2 == 0) {
        # Format: Sequence / Structure
        num_sets <- N / 2
        add_log(paste("Detected 2-line format (Seq+Str). Found", num_sets, "entries."))
        for (i in 1:num_sets) {
          pairs[[i]] <- list(
            id = paste("Protein", i), # Default ID
            seq = trimws(lines[2*i - 1]),
            exp = trimws(lines[2*i])
          )
        }
        
      } else {
        # Format: Single Sequence (1 line)
        add_log("Detected Single Sequence (No Experimental Structure).")
        pairs[[1]] <- list(id = "Protein 1", seq = trimws(lines[1]), exp = NULL)
      }
      
      return(pairs)
      
    }, error = function(e) {
      add_log(paste("Error reading file:", e$message))
      return(NULL)
    })
  })
  
  # --- AUTO-UPDATE UI: Enable/Disable Modes based on Data ---
  observe({
    data_list <- dataInput()
    if (is.null(data_list)) return()
    
    # Logic: 
    # If > 1 dataset -> Force Benchmark (5), disable others.
    
    if (length(data_list) > 1) {
      # Multi-data: Only Benchmark allowed
      updateRadioButtons(session, "mode", 
                         choices = list("5. Benchmark (Multi-Data)" = "bench"),
                         selected = "bench")
      add_log("Multi-data detected: Switched to Benchmark mode.")
      
    } else {
      # Single data
      if (is.null(data_list[[1]]$exp)) {
        # No structure: Prediction only (1-3)
        updateRadioButtons(session, "mode", 
                           choices = list(
                             "1. Chou-Fasman" = "cf",
                             "2. Improved CF" = "icf",
                             "3. GOR" = "gor"
                           ),
                           selected = "cf")
        add_log("Single sequence detected: Switched to Prediction mode (1-3).")
      } else {
        # Full data: All modes enabled
        updateRadioButtons(session, "mode", 
                           choices = list(
                             "1. Chou-Fasman" = "cf",
                             "2. Improved CF" = "icf",
                             "3. GOR" = "gor",
                             "4. Comparison" = "comp",
                             "5. Benchmark" = "bench" 
                           ),
                           selected = "cf")
        add_log("Full dataset detected: All modes enabled.")
      }
    }
  })
  
  # --- 2. Run Predictions ---
  predictions <- eventReactive(c(dataInput(), input$mode), {
    data_list <- dataInput()
    if (is.null(data_list)) return(NULL)
    
    # Executable paths
    exe_cf <- find_executable("chou_fasman", c(""), c("chou_fasman", "main"))
    exe_icf <- find_executable("improved_chou_fasman", c(""), c("improved_chou_fasman", "main"))
    exe_gor <- find_executable("GOR", c("", "predict"), c("gor_predict", "predict", "GOR", "main"))
    model_path <- file.path(normalizePath(getwd(), winslash = "/"), "GOR", "gor_model.json")
    
    # Internal function to run one prediction set
    run_one_set <- function(seq) {
      n <- nchar(seq)
      res_cf <- paste(rep("C", n), collapse="")
      res_icf <- paste(rep("C", n), collapse="")
      res_gor <- paste(rep("C", n), collapse="")
      
      # Temp files
      input_txt <- tempfile(fileext = ".txt")
      cat(seq, file = input_txt) 
      input_txt_abs <- normalizePath(input_txt, winslash = "/", mustWork = FALSE)
      
      input_fasta <- tempfile(fileext = ".fasta")
      cat(paste0(">Query\n", seq, "\n"), file = input_fasta)
      input_fasta_abs <- normalizePath(input_fasta, winslash = "/", mustWork = FALSE)
      
      # CF
      if (!is.null(exe_cf)) {
        out <- system2(exe_cf, args = c(shQuote(input_txt_abs)), stdout = TRUE, stderr = TRUE)
        res_cf <- parse_cf_output(out, n)
      }
      
      # ICF
      if (!is.null(exe_icf)) {
        out <- system2(exe_icf, args = c(shQuote(input_txt_abs)), stdout = TRUE, stderr = TRUE)
        parsed <- parse_string_output(out, n)
        if (!is.null(parsed) && nchar(parsed) == n) res_icf <- parsed
      }
      
      # GOR
      if (!is.null(exe_gor) && file.exists(model_path)) {
        args <- c("-model", shQuote(model_path), "-fasta", shQuote(input_fasta_abs))
        out <- system2(exe_gor, args = args, stdout = TRUE, stderr = TRUE)
        parsed <- parse_string_output(out, n)
        if (!is.null(parsed) && nchar(parsed) == n) res_gor <- parsed
      }
      
      unlink(c(input_txt, input_fasta))
      return(list(cf=res_cf, icf=res_icf, gor=res_gor))
    }
    
    if (input$mode == "bench") {
      # --- Benchmark Mode ---
      if (is.null(data_list[[1]]$exp)) {
        return(list(type = "error", msg = "Benchmark requires experimental data."))
      }
      
      results_df <- data.frame()
      withProgress(message = 'Running Benchmark...', value = 0, {
        for (i in seq_along(data_list)) {
          incProgress(1/length(data_list), detail = paste("Protein", i))
          d <- data_list[[i]]
          p <- run_one_set(d$seq)
          
          # Accuracies
          acc_cf <- calc_accuracy_score(p$cf, d$exp)
          acc_icf <- calc_accuracy_score(p$icf, d$exp)
          acc_gor <- calc_accuracy_score(p$gor, d$exp)
          
          # Use parsed ID for the table
          tmp <- data.frame(
            ID = d$id, 
            Method = c("Chou-Fasman", "Improved Chou-Fasman", "GOR"),
            Acc = c(acc_cf, acc_icf, acc_gor)
          )
          results_df <- rbind(results_df, tmp)
        }
      })
      add_log("Benchmark Completed.")
      return(list(type = "bench", data = results_df))
      
    } else {
      # --- Single Mode ---
      add_log("Running Single Prediction...")
      d <- data_list[[1]]
      p <- run_one_set(d$seq)
      p$type <- "single"
      return(p)
    }
  })
  
  output$debug_log <- renderText({ logs() })
  
  output$seq_info <- renderText({
    req(dataInput())
    pairs <- dataInput()
    
    info <- paste("Total Datasets Loaded:", length(pairs), 
                  "\nFirst Sequence Length:", nchar(pairs[[1]]$seq),
                  "\nFormat ID Sample:", pairs[[1]]$id)
    
    if (is.null(pairs[[1]]$exp)) {
      info <- paste(info, "\nExperimental Structure: NOT FOUND (Sequence only mode)")
    } else {
      info <- paste(info, "\nExperimental Structure: PRESENT")
    }
    return(info)
  })
  
  output$structurePlotUI <- renderUI({
    h <- "300px" 
    if (!is.null(input$mode) && (input$mode == "comp" || input$mode == "bench")) {
      h <- "500px" 
    }
    plotOutput("structurePlot", height = h)
  })
  
  output$structurePlot <- renderPlot({
    req(predictions(), dataInput())
    preds <- predictions()
    pairs <- dataInput()
    
    # Error handling
    if (!is.null(preds$type) && preds$type == "error") {
      ggplot() + 
        annotate("text", x = 0.5, y = 0.5, label = preds$msg, size = 6, color = "red") + 
        theme_void()
      return()
    }
    
    # Colors
    cols <- c("H" = "#E69191", "E" = "#F9DF91", "C" = "#A8CFE8", "X" = "#808080")
    method_cols <- c("Chou-Fasman" = "#E69191", 
                     "Improved Chou-Fasman" = "#F9DF91", 
                     "GOR" = "#A8CFE8")
    
    # Helper to draw strip
    draw_strip <- function(str, title) {
      if (is.null(str)) return(NULL)
      
      full_len <- nchar(str)
      str_sub <- str
      chars <- strsplit(str_sub, "")[[1]]
      df <- data.frame(Pos = 1:length(chars), SS = factor(chars, levels = c("H","E","C","X")))
      
      ggplot(df, aes(x=Pos, y=1, fill=SS)) + geom_tile() +
        scale_fill_manual(values = cols, drop=FALSE) +
        labs(title=paste(title, "(Full Length:", full_len, "residues)")) + theme_minimal() + 
        theme(axis.text.y=element_blank(), legend.position="bottom")
    }
    
    if (preds$type == "bench") {
      # --- Benchmark Plot ---
      df <- preds$data
      df$Method <- factor(df$Method, levels = c("Chou-Fasman", "Improved Chou-Fasman", "GOR"))
      
      # Ensure IDs are ordered by input order, not alphabetical
      df$ID <- factor(df$ID, levels = unique(df$ID))
      
      ggplot(df, aes(x=ID, y=Acc, fill=Method)) + 
        geom_col(position = position_dodge(width = 0.8), width = 0.7, color="black") +
        scale_y_continuous(labels=scales::percent, limits=c(0,1.1)) +
        scale_fill_manual(values = method_cols) +
        labs(title = "Benchmark: Accuracy Comparison Across Datasets", x = "Protein Name", y = "Accuracy") +
        theme_classic() +
        theme(text = element_text(size = 14))
      
    } else {
      # --- Single Mode Plots ---
      exp <- pairs[[1]]$exp # Might be NULL
      
      if (!is.null(input$mode)) {
        if (input$mode == "cf") draw_strip(preds$cf, "Chou-Fasman")
        else if (input$mode == "icf") draw_strip(preds$icf, "Improved CF")
        else if (input$mode == "gor") draw_strip(preds$gor, "GOR")
        else if (input$mode == "comp") {
          # Comparison Mode
          if (is.null(exp)) {
            ggplot() + 
              annotate("text", x = 0.5, y = 0.5, label = "Comparison Mode requires Experimental Structure.", size = 6) + 
              theme_void()
          } else {
            method_order <- c("Chou-Fasman", "Improved Chou-Fasman", "GOR")
            accs <- data.frame(
              Method = factor(method_order, levels = method_order),
              Acc = c(calc_accuracy_score(preds$cf, exp), 
                      calc_accuracy_score(preds$icf, exp), 
                      calc_accuracy_score(preds$gor, exp))
            )
            ggplot(accs, aes(x=Method, y=Acc, fill=Method)) + geom_col() +
              geom_text(aes(label=scales::percent(Acc, 0.1)), vjust=-0.5, size=5) +
              scale_y_continuous(labels=scales::percent, limits=c(0,1.1)) +
              scale_fill_manual(values = method_cols) + 
              labs(y = "Accuracy") + # Changed y-axis label to "Accuracy"
              theme_classic() +
              theme(text = element_text(size = 14), legend.position = "none")
          }
        }
      }
    }
  })
}

shinyApp(ui, server)
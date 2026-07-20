

#path

project_dir <- "/Users/bohe/Desktop/Metabolomics"

log2_dir <- file.path(
  project_dir,
  "Data_processing",
  "Final_feature_tables"
)

differential_dir <- file.path(
  project_dir,
  "Results",
  "Differential_analysis"
)

dir.create(
  file.path(differential_dir, "positive"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(differential_dir, "negative"),
  recursive = TRUE,
  showWarnings = FALSE
)




run_differential_analysis <- function(
    input_file,
    sample_prefix,
    output_dir,
    output_prefix
) {
  

  feature_data <- read.delim(
    input_file,
    check.names = FALSE
  )
  

  sample_columns <- grep(
    paste0("^", sample_prefix),
    colnames(feature_data),
    value = TRUE
  )
  
  

  feature_information <- feature_data[
    ,
    !colnames(feature_data) %in% sample_columns,
    drop = FALSE
  ]
  
  

  intensity_matrix <- as.matrix(
    feature_data[
      ,
      sample_columns,
      drop = FALSE
    ]
  )
  
  storage.mode(intensity_matrix) <- "numeric"
  

  sample_group <- ifelse(
    grepl("_D_", sample_columns),
    "D",
    "L"
  )
  

  d_columns <- sample_columns[
    sample_group == "D"
  ]
  
  l_columns <- sample_columns[
    sample_group == "L"
  ]
  

  mean_d <- rowMeans(
    feature_data[
      ,
      d_columns,
      drop = FALSE
    ],
    na.rm = TRUE
  )
  

  mean_l <- rowMeans(
    feature_data[
      ,
      l_columns,
      drop = FALSE
    ],
    na.rm = TRUE
  )
  
  

  log2_fc <- mean_d - mean_l
  
  

  p_value <- apply(
    intensity_matrix,
    1,
    function(feature_values) {
      
      d_values <- feature_values[
        sample_group == "D"
      ]
      
      l_values <- feature_values[
        sample_group == "L"
      ]
      
      d_values <- d_values[
        is.finite(d_values)
      ]
      
      l_values <- l_values[
        is.finite(l_values)
      ]
      
      

      if (
        length(d_values) < 2 ||
        length(l_values) < 2
      ) {
        return(NA_real_)
      }
      
      
      #
      if (
        length(unique(c(d_values, l_values))) == 1
      ) {
        return(1)
      }
      
      
      # Welch t
      test_result <- tryCatch(
        t.test(
          d_values,
          l_values,
          var.equal = FALSE
        ),
        error = function(e) NULL
      )
      
      
      if (is.null(test_result)) {
        return(NA_real_)
      }
      
      
      test_result$p.value
    }
  )
  
  
  # BH
  fdr <- p.adjust(
    p_value,
    method = "BH"
  )
  
  

  status <- ifelse(
    !is.na(fdr) &
      fdr < 0.05 &
      log2_fc >= 1,
    "Up",
    ifelse(
      !is.na(fdr) &
        fdr < 0.05 &
        log2_fc <= -1,
      "Down",
      "Not_significant"
    )
  )
  
  

  differential_results <- cbind(
    feature_information,
    Mean_D = mean_d,
    Mean_L = mean_l,
    log2FC = log2_fc,
    P_value = p_value,
    FDR = fdr,
    Status = status
  )
  
  
  # FDR
  differential_results <- differential_results[
    order(
      differential_results$FDR,
      differential_results$P_value,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  
  # select significant feature
  significant_results <- differential_results[
    differential_results$Status %in% c(
      "Up",
      "Down"
    ),
    ,
    drop = FALSE
  ]
  
  

  write.table(
    differential_results,
    file = file.path(
      output_dir,
      paste0(
        output_prefix,
        "_all_features.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  
  # save significant feature
  write.table(
    significant_results,
    file = file.path(
      output_dir,
      paste0(
        output_prefix,
        "_significant_features.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  

  result_summary <- data.frame(
    Category = c(
      "Total",
      "Up",
      "Down",
      "Not_significant"
    ),
    Number = c(
      nrow(differential_results),
      sum(
        differential_results$Status == "Up",
        na.rm = TRUE
      ),
      sum(
        differential_results$Status == "Down",
        na.rm = TRUE
      ),
      sum(
        differential_results$Status == "Not_significant",
        na.rm = TRUE
      )
    )
  )
  
  
  write.table(
    result_summary,
    file = file.path(
      output_dir,
      paste0(
        output_prefix,
        "_differential_summary.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  
  return(differential_results)
}



# positive analysis

positive_results <- run_differential_analysis(
  input_file = file.path(
    log2_dir,
    "MTBLS12932_positive_feature_table_log2.tsv"
  ),
  sample_prefix = "POS_",
  output_dir = file.path(
    differential_dir,
    "positive"
  ),
  output_prefix = "MTBLS12932_positive"
)

# negative analysis

negative_results <- run_differential_analysis(
  input_file = file.path(
    log2_dir,
    "MTBLS12932_negative_feature_table_log2.tsv"
  ),
  sample_prefix = "NEG_",
  output_dir = file.path(
    differential_dir,
    "negative"
  ),
  output_prefix = "MTBLS12932_negative"
)



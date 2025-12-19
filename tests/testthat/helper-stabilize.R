#' Stabilize a character vector for snapshot testing
#'
#' Finds all number-like substrings in a character vector, converts them
#' to numeric, formats them to a consistent precision, and substitutes them
#' back into the strings.
#'
#' @param text_vector The character vector to stabilize.
#' @param sig_figs The number of significant figures to keep.
#'
#' @return The stabilized character vector.
#' Stabilize a character vector for snapshot testing (Robust Version)
#'
#' Deconstructs each string into number and non-number parts, processes the
#' numbers, and reconstructs the string to perfectly preserve structure.
#'
#' @param text_vector The character vector to stabilize.
#' @param sig_figs The number of significant figures to keep.
#'
#' @return The stabilized character vector.
stabilize_text_snapshot <- function(text_vector, sig_figs = 8) {

  # A robust regex to find numbers (including integers, decimals, and sci-notation)
  num_regex <- "[+-]?\\d*\\.?\\d+(?:[Ee][+-]?\\d+)?"

  sapply(text_vector, FUN = function(line) {
    if (!is.character(line) || is.na(line)) {
      return(line)
    }

    # 1. Find the positions of all numbers
    matches <- gregexpr(num_regex, line)
    if (matches[[1]][1] == -1) {
      return(line) # No numbers found, return line as-is
    }

    # 2. Extract the number strings and the non-number "scaffolding" separately
    number_strings <- regmatches(line, matches)[[1]]
    scaffolding <- regmatches(line, matches, invert = TRUE)[[1]]

    # 3. Process only the number strings
    processed_numbers <- sapply(number_strings, function(num_str) {
      num <- as.numeric(num_str)
      # Use sprintf for reliable, non-padded formatting
      sprintf(paste0("%.", sig_figs, "g"), num)
    })

    # 4. Weave the scaffolding and processed numbers back together
    # The result is an alternating sequence:
    # scaffold[1], processed_num[1], scaffold[2], processed_num[2], ...
    result <- character(length(scaffolding) + length(processed_numbers))
    result[seq(1, by = 2, length.out = length(scaffolding))] <- scaffolding
    result[seq(2, by = 2, length.out = length(processed_numbers))] <- processed_numbers

    paste(result, collapse = "")

  }, USE.NAMES = FALSE)
}



#' Smart wrapper to stabilize any object for snapshot testing
#'
#' Automatically selects the correct stabilization method based on the input type.
#' - If `x` is a character vector, it stabilizes numbers within the text.
#' - If `x` is any other object (list, data.frame, numeric), it recursively
#'   rounds all numeric elements.
#'
#' @param x The object to stabilize.
#' @param digits The number of decimal places for numeric objects.
#' @param sig_figs The number of significant figures for text-based numbers.
#'
#' @return The stabilized object.
#' Recursively stabilize any object for snapshot testing
#'
#' This function intelligently handles any R object by applying the correct
#' stabilization rule at every level of a nested structure.
#' - Numeric vectors/matrices are rounded.
#' - Character vectors have numbers within them rounded and formatted.
#' - Lists and data.frames are recursively processed element by element or column by column.
#'
#' @param x The object to stabilize.
#' @param digits The number of decimal places for numeric objects.
#' @param sig_figs The number of significant figures for text-based numbers.
#'
#' @return The stabilized object.
stabilize <- function(x, digits = 8, sig_figs = 8) {
  # Base case 1: If it's a numeric vector or matrix, round it.
  if (is.numeric(x)) {
    return(round(x, digits))
  }

  # Base case 2: If it's a character vector, stabilize the text.
  if (is.character(x)) {
    return(stabilize_text_snapshot(x, sig_figs = sig_figs))
  }

  # Recursive step for data.frames: apply stabilize to each column.
  if (is.data.frame(x)) {
    # The x[] is important to preserve the data.frame class
    x[] <- lapply(x, stabilize, digits = digits, sig_figs = sig_figs)
    return(x)
  }

  # Recursive step for lists: apply stabilize to each element.
  if (is.list(x)) {
    return(lapply(x, stabilize, digits = digits, sig_figs = sig_figs))
  }

  # Fallback for any other type (factors, dates, etc.)
  return(x)
}

# Add this to your helper-stabilize.R file
# You can delete the old clean_condition_message function.

#' Extracts the core message(s) from a formatted testthat condition snapshot
#' @param text The formatted character string from the snapshot.
#' @return The core message text, stripped of version-specific headers.
extract_core_messages <- function(text) {
  # Split the full text into individual lines
  lines <- strsplit(text, "\n")[[1]]

  # Keep only the lines that are the actual message, filtering out headers
  # from both old and new testthat/rlang versions.
  core_lines <- grep("^(<|Error|Warning|!|\\s*$)", lines, value = TRUE, invert = TRUE)

  # Trim whitespace and paste back together
  paste(trimws(core_lines), collapse = "\n")
}

#' Standardize ggplot_build() data for snapshotting
#'
#' Handles changes in ggplot2's internal data structures across versions,
#' like the renaming of the 'linewidth' aesthetic to 'size'.
#'
#' @param plot_data_list The list of data frames from `ggplot_build(p)$data`.
#' @return A list of data frames with standardized column names.
standardize_plot_data <- function(plot_data_list) {
  lapply(plot_data_list, function(df) {
    # In ggplot2 v3.4.0, `size` aesthetic for lines was renamed to `linewidth`.
    # We standardize it back to `size` for snapshot consistency.
    if ("linewidth" %in% names(df)) {
      names(df)[names(df) == "linewidth"] <- "size"
    }
    df
  })
}

# Add this to helper-stabilize.R

#' Stabilize random temp file paths in model code snapshots
#'
#' Finds the $DATA line in a model file and replaces the random filepath
#' with a consistent placeholder. It also stabilizes the rest of the object.
#'
#' @param result_list A list returned from updateFREMmodel, containing $data and $model.
#' @return A stabilized list suitable for snapshotting.
stabilize_model_paths <- function(result_list) {
  # First, stabilize the data components recursively as usual
  stable_list <- stabilize(result_list)

  # Then, specifically fix the random path in the model code
  model_code <- stable_list$model
  if (!is.null(model_code)) {
    data_line_index <- grep("^\\$DATA", model_code)
    if (length(data_line_index) > 0) {
      # Use gsub to replace the random path with a placeholder
      stable_list$model[data_line_index] <- gsub(
        pattern = "(\\$DATA\\s+).*( IGNORE=@)",
        replacement = "\\1[placeholder_path]\\2",
        x = model_code[data_line_index]
      )
    }
  }

  return(stable_list)
}

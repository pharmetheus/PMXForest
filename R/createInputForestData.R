#' Create Input Data Frame for a Forest Plot
#'
#' @description This function transforms a specially formatted list of covariates
#'   and their values into a data frame suitable for generating a forest plot.
#'   It is designed to work with the output of the `getCovStats` function.
#'
#' @param listCovs A list where each element represents a covariate or a group
#'   of covariates.
#'   \itemize{
#'     \item For simple (univariate) effects, elements should be a named vector
#'       (e.g., `list(WT = c(70, 110))`).
#'     \item For grouped (multivariate) effects, elements should be a nested,
#'       named list of vectors (e.g., for a categorical variable `GENO`,
#'       `list(GENO = list(GENO_2 = c(0,1,0), GENO_3 = c(0,0,1)))`).
#'   }
#' @param iMiss The numeric value that should be used to fill in columns that are
#'   not applicable for a given row. Defaults to -99.
#'
#' @return A data frame where each row represents a specific covariate value
#'   (or combination of values for grouped effects). The data frame contains a
#'   `COVARIATEGROUPS` column derived from the names in the input list, which is
#'   used to group effects in a forest plot.
#' @export
#'
#' @examples
#' # --- Example 1: Simple (Univariate) Covariates ---
#' # Each element in the list becomes its own group in the output data frame.
#' simple_covs <- list(
#'   SEX = c(0, 1),
#'   WT = c(70, 110)
#' )
#' df_simple <- createInputForestData(simple_covs)
#' print(df_simple)
#'
#'
#' # --- Example 2: Grouped (Multivariate) Covariates ---
#' # A nested list is used to group related effects. The outer list name ("GENO")
#' # becomes the value in the COVARIATEGROUPS column for all derived rows.
#' # This is the standard format for multi-level categorical variables.
#' grouped_covs <- list(
#'   GENO = list(
#'     GENO_2 = c(0, 1, 0), # Represents 3 levels, with level 1 as reference (0,0)
#'     GENO_3 = c(0, 0, 1)
#'   )
#' )
#' df_grouped <- createInputForestData(grouped_covs)
#' print(df_grouped)
#'
#'
#' # --- Example 3: Full Workflow with getCovStats ---
#' # The primary use case is to transform the output of getCovStats into a
#' # data frame ready for plotting.
#'
#' # 1. Create sample data
#' sample_data <- data.frame(
#'   ID = c(1, 1, 2, 3, 4, 5),
#'   WT = c(60.5, 61.0, 70.2, 80.8, 65.1, 90.3),
#'   SEX = c(0, 0, 1, 0, 1, 0),
#'   RACE = c(1, 1, 2, 3, 1, 2)
#' )
#'
#' # 2. Get covariate statistics using getCovStats. This creates a list
#' # containing both simple (WT, SEX) and nested (RACE) elements.
#' stats_list <- getCovStats(
#'   data = sample_data,
#'   covariates = c("WT", "SEX", "RACE"),
#'   minLevels = 4 # Ensure WT is continuous and RACE is categorical
#' )
#'
#' # 3. Use the output from getCovStats to create the final data frame.
#' forest_plot_data <- createInputForestData(stats_list)
#' print(forest_plot_data)
createInputForestData <- function(listCovs, iMiss = -99) {
  fixedname <- "COVARIATEGROUPS"
  cstrUniqueCols <- c()

  # Use seq_along for safe iteration over lists, including empty ones.
  for (i in seq_along(listCovs)) {
    if (is.list(listCovs[[i]])) {
      # It's safer to also use seq_along on the inner loop
      for (j in seq_along(listCovs[[i]])) {
        cstrUniqueCols <- c(cstrUniqueCols, names(listCovs[[i]][j]))
      }
    } else {
      cstrUniqueCols <- c(cstrUniqueCols, names(listCovs[i]))
    }
  }

  cstrUniqueCols <- unique(cstrUniqueCols)
  ncols <- length(cstrUniqueCols)

  # Create a template row
  df <- data.frame(matrix(NA, nrow = 1, ncol = ncols + 1))
  names(df) <- c(cstrUniqueCols, fixedname)

  # Use seq_along for the main processing loop as well
  for (i in seq_along(listCovs)) {
    if (is.list(listCovs[[i]])) {
      strNameBig <- names(listCovs[i])
      strVal1 <- listCovs[[i]][[1]]
      if (length(strVal1) > 0) {
        for (k in 1:length(strVal1)) {
          dfr <- df[1, ]
          dfr[, ] <- iMiss
          for (j in seq_along(listCovs[[i]])) {
            strName <- names(listCovs[[i]][j])
            # Check for out-of-bounds access in case of unequal vector lengths
            if (k <= length(listCovs[[i]][[j]])) {
              if (!is.na(listCovs[[i]][[j]][k])) dfr[[strName]] <- listCovs[[i]][[j]][k]
            }
          }
          dfr[[fixedname]] <- strNameBig
          df <- rbind(df, dfr)
        }
      }
    } else {
      vecr <- listCovs[[i]]
      strName <- names(listCovs[i])
      if (length(vecr) > 0) {
        for (k in seq_along(vecr)) {
          dfr <- df[1, ]
          dfr[, ] <- iMiss
          dfr[[strName]] <- vecr[k]
          dfr[[fixedname]] <- strName
          df <- rbind(df, dfr)
        }
      }
    }
  }

  # Remove the initial template row
  if (nrow(df) > 1) {
    df <- df[-1, , drop = FALSE]
  } else {
    df <- df[0, , drop = FALSE]
  }

  # FIX: Reset row names for clean, consistent output
  row.names(df) <- NULL

  return(df)
}

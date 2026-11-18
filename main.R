# Internal helper: encode DNA 5-mers as categorical factor columns.
#' Encode DNA 5-mers as Categorical Factor Columns
#'
#' Splits each DNA 5-mer into five nucleotide position columns and casts each
#' column to a factor with levels \code{A}, \code{T}, \code{C} and \code{G},
#' which is the encoding the bundled random forest model was trained on.
#'
#' This function is an internal helper used by \code{\link{prediction_multiple}}
#' and is not exported.
#'
#' @param dna_strings Character vector of five-base DNA sequences, e.g.
#'   \code{"GGACA"}. Every element must have the same number of bases.
#' @return A data frame with one factor column per nucleotide position, named
#'   \code{nt_pos1} to \code{nt_pos5}, and one row per input sequence. Each
#'   column is a factor with levels \code{A}, \code{T}, \code{C}, \code{G}.
#' @keywords internal
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict m6A Status for a Batch of Feature Rows
#'
#' Uses the bundled random forest model to predict the m6A probability and
#' status for every row of a feature table.
#'
#' @param ml_fit Fitted random forest model, typically loaded from the
#'   bundled file with \code{readRDS(system.file("extdata", "rf_fit.rds",
#'   package = "m6APrediction"))}.
#' @param feature_df Data frame with one row per candidate site and the
#'   required columns \code{gc_content}, \code{RNA_type}, \code{RNA_region},
#'   \code{exon_length}, \code{distance_to_junction},
#'   \code{evolutionary_conservation} and \code{DNA_5mer}. An error is raised
#'   if any of these columns is missing.
#' @param positive_threshold Numeric score threshold between 0 and 1, used to
#'   label sites; defaults to 0.5. A probability strictly greater than the
#'   threshold is labelled \code{"Positive"}; equality is \code{"Negative"}.
#' @return The input \code{feature_df} with two appended columns:
#'   \code{predicted_m6A_prob}, the numeric probability of the Positive class,
#'   and \code{predicted_m6A_status}, a character label that is either
#'   \code{"Positive"} or \code{"Negative"}.
#' @import randomForest
#' @importFrom stats predict
#' @export
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'   package = "m6APrediction", mustWork = TRUE))
#' example_df <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'   package = "m6APrediction", mustWork = TRUE))
#' prediction_multiple(ml_fit, example_df[1:5, ], positive_threshold = 0.6)
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df)))
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[, "Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df)
}

#' Predict m6A Status for a Single Candidate Site
#'
#' Convenience wrapper around \code{\link{prediction_multiple}} for one site,
#' returning the predicted probability and status as a named character vector.
#'
#' @inheritParams prediction_multiple
#' @param gc_content Numeric GC content of the region, e.g. \code{0.6}.
#' @param RNA_type Character RNA biotype; one of \code{"mRNA"},
#'   \code{"lincRNA"}, \code{"lncRNA"} or \code{"pseudogene"}.
#' @param RNA_region Character region of the transcript; one of \code{"CDS"},
#'   \code{"intron"}, \code{"3'UTR"} or \code{"5'UTR"}.
#' @param exon_length Numeric exon length. The course data are already
#'   log2-transformed, so pass the transformed value used during training.
#' @param distance_to_junction Numeric distance to the nearest splice junction,
#'   also already log2-transformed.
#' @param evolutionary_conservation Numeric evolutionary conservation score.
#' @param DNA_5mer Character five-base sequence context, e.g. \code{"GGACA"}.
#'   The model was trained on contexts with A in the central position.
#' @return A named character vector of length two with elements
#'   \code{predicted_m6A_prob} (the probability rendered as a character) and
#'   \code{predicted_m6A_status} (\code{"Positive"} or \code{"Negative"}).
#' @export
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'   package = "m6APrediction", mustWork = TRUE))
#' prediction_single(ml_fit, gc_content = 0.6, RNA_type = "mRNA",
#'   RNA_region = "CDS", exon_length = 12, distance_to_junction = 5,
#'   evolutionary_conservation = 0.8, DNA_5mer = "GGACA",
#'   positive_threshold = 0.5)
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  pred_outcome <- prediction_multiple(ml_fit = ml_fit,
                                      data.frame(gc_content = gc_content,
                                                 RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                                                 RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                                                 exon_length = exon_length,
                                                 distance_to_junction = distance_to_junction,
                                                 evolutionary_conservation = evolutionary_conservation,
                                                 DNA_5mer = DNA_5mer),
                                      positive_threshold = positive_threshold)
  returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector)
}

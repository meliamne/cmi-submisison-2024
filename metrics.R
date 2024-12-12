library(tidymodels)
library(tidyverse)
library(parsnip)
library(ggplot2)
library(recipes)
library(yardstick)

# Define the Spearman correlation metric
spearman_cor_abs_impl <- function(truth, estimate, case_weights) {
  abs(cor(truth, estimate, method = "spearman"))
}

spearman_cor_abs_vec=function(truth, estimate, na_rm=TRUE,case_weights=NULL){
	check_numeric_metric(truth, estimate, case_weights)
	if (na_rm) {
		result <- yardstick_remove_missing(truth, estimate, case_weights)
		truth <- result$truth
		estimate <- result$estimate
		case_weights <- result$case_weights
	} else if (yardstick_any_missing(truth, estimate, case_weights)) {
		return(NA_real_)
	}
	spearman_cor_abs_impl(truth, estimate, case_weights = case_weights)
}

spearman_cor_abs <- function(data, ...) {
  UseMethod("spearman_cor_abs")
}

spearman_cor_abs <- new_numeric_metric(spearman_cor_abs, direction = "maximize")
	spearman_cor_abs.data.frame <- function(data, truth, estimate, na_rm = TRUE, case_weights = NULL, ...) {
	numeric_metric_summarizer(
		name = "spearman_cor_abs",
		fn = spearman_cor_abs_vec,
		data = data,
		truth = !!enquo(truth),
		estimate = !!enquo(estimate),
		na_rm = na_rm,
		case_weights = !!enquo(case_weights)
	)
}

#' bnec_newdata.drc
#' 
#' Create a dataset for predictions
#'
#' @inheritParams bnec_newdata
#' @inherit bnec_newdata description return examples
#' 
#' @importFrom stats model.frame
#' @noRd
#' @export
bnec_newdata.drc <- function(x, precision = 100, x_range = NA) {
  check_args_newdata(precision, x_range)
  data <- model.frame(formula(x), data = x$data) |> 
    data.frame()
  x_var <- colnames(data)[ncol(data)] # assume x is last column of model.frame
  fit <- x$fit
  x_vec <- data[, x_var]
  if (any(is.na(x_range))) {
    x_seq <- seq(min(x_vec), max(x_vec), length = precision)
  } else {
    x_seq <- seq(min(x_range), max(x_range), length = precision)
  }
  newdata <- data.frame(x_seq)
  names(newdata) <- x_var
  #fam_tag <- fit$family$family
  #custom_name <- check_custom_name(fit$family)
  #if (fam_tag == "binomial" || custom_name == "beta_binomial2") {
  #  trials_var <- attr(data, "bnec_pop")[["trials_var"]]
  #  newdata[[trials_var]] <- 1
  #}
  newdata
}

#' @noRd
newdata_eval <- function(object, precision, x_range) {
  # Just need one model to extract and generate data
  # since all models are considered to have the exact same raw data.
  # if (inherits(object, "bayesmanecfit")) {
  #   model_set <- names(object$mod_fits)
  #   object <- suppressMessages(pull_out(object, model = model_set[1]))
  # }
  data <- model.frame(formula(object), object$data)
  #bnec_pop_vars <- attr(data, "bnec_pop")
  newdata <- bnec_newdata(object, precision = precision, x_range = x_range)
  x_var <- colnames(data)[ncol(data)] # assume x is last column of model.frame  
  x_vec <- newdata[[x_var]]
  list(newdata = newdata, x_vec = x_vec)
}

#' @noRd
newdata_eval_fitted <- function(object, precision, x_range, make_newdata,
                                fct_eval, ...) {
  # Just need one model to extract and generate data
  # since all models are considered to have the exact same raw data.
  if (inherits(object, "bayesmanecfit")) {
    model_set <- names(object$mod_fits)
    object <- suppressMessages(pull_out(object, model = model_set[1]))
  }
  data <- model.frame(object$bayesnecformula, object$fit$data)
  bnec_pop_vars <- attr(data, "bnec_pop")
  dot_list <- list(...)
  if ("newdata" %in% names(dot_list) && make_newdata) {
    stop("You cannot provide a \"newdata\" argument and set",
         " make_newdata = TRUE at the same time. Please use one or another.",
         " See details in help file ?", fct_eval)
  }
  if (!("newdata" %in% names(dot_list))) {
    if (make_newdata) {
      newdata <- bnec_newdata(object, precision = precision, x_range = x_range)
      x_vec <- newdata[[bnec_pop_vars[["x_var"]]]]
      if ("re_formula" %in% names(dot_list)) {
        message("Argument \"re_formula\" ignored and set to NA because",
                " function bnec_newdata cannot guess random effect structure.")
      }
      re_formula <- NA
    } else {
      newdata <- NULL
      x_vec <- pull_brmsfit(object)$data[[bnec_pop_vars[["x_var"]]]]
      precision <- "from raw data"
      if (!("re_formula" %in% names(dot_list))) {
        re_formula <- NULL
      } else {
        re_formula <- dot_list$re_formula
      }
    }
  } else {
    newdata <- dot_list$newdata
    x_vec <- newdata[[bnec_pop_vars[["x_var"]]]]
    precision <- "from user-specified newdata"
    if (!("re_formula" %in% names(dot_list))) {
      re_formula <- NULL
    } else {
      re_formula <- dot_list$re_formula
    }
  }
  list(newdata = newdata, x_vec = x_vec, precision = precision,
       re_formula = re_formula)
}

#' min_abs
#' @param x A \code{\link[base]{numeric}} vector.
#' @return A \code{\link[base]{numeric}} vector.
#' @noRd
min_abs <- function(x) {
  which.min(abs(x))
}

#' @noRd
clean_names <- function(x) {
  paste0("Q", gsub("%", "", names(x), fixed = TRUE))
}

#' @noRd
#' @importFrom chk chk_numeric
check_args_newdata <- function(precision, x_range) {
  chk_numeric(precision)
  if (!is.na(x_range[1])) {
    chk_numeric(x_range)
  }  
}




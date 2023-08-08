#' Extracts the predicted NSEC value from supported classes.
#'
#' @param object An object of class currently supported. Includes 
#' \code{\link{bayesnecfit}} or
#' \code{\link{bayesmanecfit}} returned by \code{\link{bnec}} or  
#' \code{\link{drc}} returned by \code{\link{drc}}.
#' @param sig_val Probability value to use as the lower quantile to test
#' significance of the predicted posterior values.
#' against the lowest observed concentration (assumed to be the control), to
#' estimate NEC as an interpolated NOEC value from smooth ECx curves.
#' @param precision The number of unique x values over which to find NSEC -
#' large values will make the NSEC estimate more precise.
#' @param posterior A \code{\link[base]{logical}} value indicating if the full
#' posterior sample of calculated NSEC values should be returned instead of
#' just the median and 95 credible intervals. Not relevant to \code{\link{drc}}
#' @param hormesis_def A \code{\link[base]{character}} vector, taking values
#' of "max" or "control". See Details.
#' @param xform A function to apply to the returned estimated concentration
#' values.
#' @param x_range A range of x values over which to consider extracting NSEC.
#' @param prob_vals A vector indicating the probability values over which to
#' return the estimated NSEC value. Defaults to 0.5 (median) and 0.025 and
#' 0.975 (95 percent credible intervals).
#'
#' @details For \code{hormesis_def}, if "max", then NSEC values are calculated
#' as a decline from the maximum estimates (i.e. the peak at NEC);
#' if "control", then ECx values are calculated relative to the control, which
#' is assumed to be the lowest observed concentration. Not current supported for 
#'  \code{\link{drc}} fits.
#' 
#' The argument \code{precision} controls how precisely the
#' \code{\link{ecx}} or \code{\link{nsec}} value is estimated, with 
#' argument \code{x_range} allowing estimation beyond the existing range of
#' the observed data (otherwise the default range) which can be useful in a
#' small number of cases. There is also no reasonable case where estimating
#' these from the raw data would be of value, because both functions would
#' simply return one of the treatment concentrations, making NOEC a better
#' metric in that case.
#'
#' @return A vector containing the estimated NSEC value, including upper and
#' lower 95% credible interval bounds.
#'
#' @examples
#' \donttest{
#' library(bayesnec)
#'
#' data(manec_example)
#' nsec(manec_example)
#' }
#'
#' @export
#' 
nsec <- function(object, sig_val = 0.01, precision = 1000,
                 posterior = FALSE, x_range = NA, hormesis_def = "control",
                 xform = identity, prob_vals = c(0.5, 0.025, 0.975)) {
  UseMethod("nsec")
}

#' @inheritParams nsec
#'
#' @param object An object of class \code{\link{drc}} returned by
#' \code{\link{drc}}.
#'
#' @inherit nsec
#' 
#' @importFrom chk chk_logical chk_numeric
#' 
#' @noRd
#'
#' @export
nsec.drc <- function(object, sig_val = 0.01, precision = 1000,
                             x_range = NA,
                             hormesis_def = "control", xform = identity,
                             prob_vals = c(0.5, 0.025, 0.975)) {
  chk_numeric(sig_val)
  chk_numeric(precision)

  if (length(sig_val)>1) {
    stop("You may only pass one sig_val")  
  }
  if ((hormesis_def %in% c("max", "control")) == FALSE) {
    stop("type must be one of \"max\" or \"control\" (the default). ",
         "Please see ?ecx for more details.")
  }
  if(!inherits(xform, "function")) { 
    stop("xform must be a function.")}  
  if (length(prob_vals) < 3 | prob_vals[1] < prob_vals[2] |
      prob_vals[1] > prob_vals[3] | prob_vals[2] > prob_vals[3]) {
    stop("prob_vals must include central, lower and upper quantiles,",
         " in that order.")
  }

  # if (length(grep("ecx", object$model)) > 0) {
  #   mod_class <- "ecx"
  # } else {
  #   mod_class <- "nec"
  # }
  newdata_list <- newdata_eval(
    object, precision = precision, x_range = x_range
  )
  p_samples <- predict(object, newdata = newdata_list$newdata,
                       interval = "confidence", level = prob_vals[3]-prob_vals[2])
  x_vec <- newdata_list$x_vec
  
  # calculate the reference level
  ref_dat <- data.frame(min(newdata_list$newdata))
  colnames(ref_dat) <- colnames(newdata_list$newdata)
  reference <- predict(object, newdata = ref_dat,
                       interval = "confidence" , 
                       level = 1-(sig_val*2))["Lower"]

  #   if (grepl("horme", object$model)) {
  #   n <- seq_len(nrow(p_samples))
  #   p_samples <- do_wrapper(n, modify_posterior, object, x_vec,
  #                           p_samples, hormesis_def, fct = "rbind")
  #   nec_posterior <- as_draws_df(object$fit)[["b_nec_Intercept"]]
  #   if (hormesis_def == "max") {
  #     reference <- quantile(apply(p_samples, 2, max), probs = sig_val)
  #   }
  # }
  nsec_out <- apply(p_samples, 2, nsec_fct, reference, x_vec)
  # formula <- formula(object)
  # x_call <- x$call
  # if (inherits(x_call, "call")) {
  #   x_call[[2]] <- str2lang("nsec_out")
  #   nsec_out <- eval(x_call)
  # }
  if (inherits(xform, "function")) {
    nsec_out <- xform(nsec_out)
  }
  nsec_estimate <-nsec_out
  names(nsec_estimate) <- clean_names(nsec_estimate)
  attr(nsec_estimate, "precision") <- precision
  attr(nsec_out, "precision") <- precision
  attr(nsec_estimate, "sig_val") <- sig_val
  attr(nsec_out, "sig_val") <- sig_val
  attr(nsec_estimate, "toxicity_estimate") <- "nsec"
  attr(nsec_out, "toxicity_estimate") <-  "nsec"
  nsec_estimate

}


#' @noRd
nsec_fct <- function(y, reference, x_vec) {
  x_vec[min_abs(y - reference)]
}

#' @inheritParams nsec
#'
#' @param object An object of class \code{\link{bayesnecfit}} returned by
#' \code{\link{bnec}}.
#'
#' @inherit nsec
#' 
#' @importFrom stats quantile
#' @importFrom brms as_draws_df posterior_epred
#' @importFrom chk chk_logical chk_numeric
#' 
#' @noRd
#'
#' @export
nsec.bayesnecfit <- function(object, sig_val = 0.01, precision = 1000,
                             posterior = FALSE, x_range = NA,
                             hormesis_def = "control", xform = identity,
                             prob_vals = c(0.5, 0.025, 0.975)) {
  chk_numeric(sig_val)
  chk_numeric(precision)
  chk_logical(posterior)
  if (length(sig_val)>1) {
    stop("You may only pass one sig_val")  
  }
  if ((hormesis_def %in% c("max", "control")) == FALSE) {
    stop("type must be one of \"max\" or \"control\" (the default). ",
         "Please see ?ecx for more details.")
  }
  if(!inherits(xform, "function")) { 
    stop("xform must be a function.")}  
  if (length(prob_vals) < 3 | prob_vals[1] < prob_vals[2] |
      prob_vals[1] > prob_vals[3] | prob_vals[2] > prob_vals[3]) {
    stop("prob_vals must include central, lower and upper quantiles,",
         " in that order.")
  }
  if (length(grep("ecx", object$model)) > 0) {
    mod_class <- "ecx"
  } else {
    mod_class <- "nec"
  }
  newdata_list <- newdata_eval(
    object, precision = precision, x_range = x_range
  )
  p_samples <- posterior_epred(object, newdata = newdata_list$newdata,
                               re_formula = NA)
  x_vec <- newdata_list$x_vec
  reference <- quantile(p_samples[, 1], sig_val)
  if (grepl("horme", object$model)) {
    n <- seq_len(nrow(p_samples))
    p_samples <- do_wrapper(n, modify_posterior, object, x_vec,
                            p_samples, hormesis_def, fct = "rbind")
    nec_posterior <- as_draws_df(object$fit)[["b_nec_Intercept"]]
    if (hormesis_def == "max") {
      reference <- quantile(apply(p_samples, 2, max), probs = sig_val)
    }
  }
  nsec_out <- apply(p_samples, 1, nsec_fct, reference, x_vec)
  formula <- object$bayesnecformula
  x_str <- grep("crf(", labels(terms(formula)), fixed = TRUE, value = TRUE)
  x_call <- str2lang(eval(parse(text = x_str)))
  if (inherits(x_call, "call")) {
    x_call[[2]] <- str2lang("nsec_out")
    nsec_out <- eval(x_call)
  }
  if (inherits(xform, "function")) {
    nsec_out <- xform(nsec_out)
  }
  nsec_estimate <- quantile(unlist(nsec_out), probs = prob_vals)
  names(nsec_estimate) <- clean_names(nsec_estimate)
  attr(nsec_estimate, "precision") <- precision
  attr(nsec_out, "precision") <- precision
  attr(nsec_estimate, "sig_val") <- sig_val
  attr(nsec_out, "sig_val") <- sig_val
  attr(nsec_estimate, "toxicity_estimate") <- "nsec"
  attr(nsec_out, "toxicity_estimate") <-  "nsec"
  if (!posterior) {
    nsec_estimate
  } else {
    nsec_out
  }
}

#' @inheritParams nsec
#'
#' @param object An object of class \code{\link{bayesmanecfit}} returned by
#' \code{\link{bnec}}.
#'
#' @inherit nsec
#' 
#' @importFrom stats quantile
#'
#' @noRd
#'
#' @export
nsec.bayesmanecfit <- function(object, sig_val = 0.01, precision = 1000,
                               posterior = FALSE, x_range = NA,
                               hormesis_def = "control", xform = identity,
                               prob_vals = c(0.5, 0.025, 0.975)) {
  if (length(sig_val)>1) {
    stop("You may only pass one sig_val")  
  }
  sample_nsec <- function(x, object, sig_val, precision,
                          posterior, hormesis_def,
                          x_range, xform, prob_vals, sample_size) {
    mod <- names(object$mod_fits)[x]
    target <- suppressMessages(pull_out(object, model = mod))
    out <- nsec(target, sig_val = sig_val, precision = precision,
                posterior = posterior, hormesis_def = hormesis_def,
                x_range = x_range, xform = xform, prob_vals = prob_vals)
    n_s <- as.integer(round(sample_size * object$mod_stats[x, "wi"]))
    sample(out, n_s)
  }
  sample_size <- object$sample_size
  to_iter <- seq_len(length(object$success_models))
  nsec_out <- sapply(to_iter, sample_nsec, object, sig_val, precision,
                     posterior = TRUE, hormesis_def, x_range,
                     xform, prob_vals, sample_size)
  nsec_out <- unlist(nsec_out)
  nsec_estimate <- quantile(nsec_out, probs = prob_vals)
  names(nsec_estimate) <- clean_names(nsec_estimate)
  attr(nsec_estimate, "precision") <- precision
  attr(nsec_out, "precision") <- precision
  attr(nsec_estimate, "sig_val") <- sig_val
  attr(nsec_out, "sig_val") <- sig_val
  attr(nsec_estimate, "toxicity_estimate") <- "nsec"
  attr(nsec_out, "toxicity_estimate") <-  "nsec"
  if (!posterior) {
    nsec_estimate
  } else {
    nsec_out
  }
}

#' @noRd
nsec_fct <- function(y, reference, x_vec) {
  x_vec[min_abs(y - reference)]
}

#' @inheritParams nsec
#'
#' @param object An object of class \code{\link{brmsfit}} returned by
#' \code{\link{brms}}.
#' @param x_var A character indicating the name of the predictor (x) data in object
#' @param group_var A character indicating the name of the grouping variable in object
#' @param by_group A logical indicating if nsec values should be returned for 
#' each level in group_var, or marginalised across all groups.
#'
#' @inherit nsec
#' 
#' @importFrom stats quantile
#' @importFrom brms as_draws_df posterior_epred
#' @importFrom chk chk_logical chk_numeric
#' 
#' @noRd
#'
#' @export
nsec.brmsfit <- function(object, 
                        x_var, 
                        group_var, 
                        by_group = TRUE,
                        probs = c(0.025, 0.5, 0.975),
                        precision = 1000,
                        sig_val = 0.01,
                        posterior = FALSE, 
                        x_range = NA,
                        horme = FALSE,
                        hormesis_def = "control", 
                        xform = identity
){
  chk_numeric(sig_val)
  chk_numeric(precision)
  chk_logical(posterior)
  if (length(sig_val)>1) {
    stop("You may only pass one sig_val")  
  }
  if ((hormesis_def %in% c("max", "control")) == FALSE) {
    stop("type must be one of \"max\" or \"control\" (the default). ",
         "Please see ?ecx for more details.")
  }
  if(!inherits(xform, "function")) { 
    stop("xform must be a function.")}  
  if (length(prob_vals) < 3 | prob_vals[1] < prob_vals[2] |
      prob_vals[1] > prob_vals[3] | prob_vals[2] > prob_vals[3]) {
    stop("prob_vals must include central, lower and upper quantiles,",
         " in that order.")
  }
  if (missing(group_var)) {
    stop("group_var must be supplied.")    
  }
  if(is.na(x_range)){
    x_range = range(object$data[x_var])
  }
  x_vec <- seq(min(x_range), max(x_range), length=precision)
  
  groups <-  unlist(unique(object$data[group_var]))
  out_vals <- lapply(groups, FUN = function(g){
    dat_list <- list(x_vec,
                     g) 
    names(dat_list) <- c(x_var, group_var)
    pred_dat <- expand.grid(dat_list)
    
    p_samples <- posterior_epred(object, newdata = pred_dat, re_formula = NA)
    reference <- quantile(p_samples[, 1], sig_val)
    
    if (horme) {
      n <- seq_len(nrow(p_samples))
      p_samples <- bayesnec:::do_wrapper(n, bayesnec:::modify_posterior, object, x_vec,
                                         p_samples, hormesis_def, fct = "rbind")
      nec_posterior <- as_draws_df(object$fit)[["b_nec_Intercept"]]
      if (hormesis_def == "max") {
        reference <- quantile(apply(p_samples, 2, max), probs = sig_val)
      }
    }    
    
    nsec_out <- apply(p_samples, 1, nsec_fct, reference, x_vec)
    unlist(nsec_out)
    
  })
  
  if(by_group & posterior){   
    names(out_vals) <- groups
    out_vals <- out_vals |> bind_cols() |> 
      pivot_longer(everything(), names_to = group_var, values_to = "NSEC")
  }
  
  if(!by_group & posterior){   
    out_vals <- as.numeric(unlist(out_vals))
  }
  
  if(by_group & !posterior){   
    names(out_vals) <- groups
    out_vals <- lapply(out_vals, quantile, probs = probs) |> 
      bind_rows(.id = group_var)
  }
  
  if(!by_group & !posterior){   
    out_vals <- quantile(unlist(out_vals), probs = probs)
  }
  
  attr(out_vals, "precision") <- precision
  attr(out_vals, "sig_val") <- sig_val
  attr(out_vals, "toxicity_estimate") <- "nsec"
  return(out_vals)
}



#' @inheritParams nsec
#'
#' @param object An object of class \code{\link{drc}} returned by
#' \code{\link{drc}}.
#'
#' @inherit nsec details seealso return examples
#' 
#' @importFrom bayesnec nsec
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

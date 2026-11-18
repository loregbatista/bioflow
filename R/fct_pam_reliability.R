#' Reliability matrix used to weight the PAM selection index
#'
#' Builds a designation x trait matrix of reliabilities in [0, 1]. Per trait:
#' (1) use the reported `reliability` when at least one value is usable;
#' (2) otherwise (e.g. designation fitted as fixed effect, so reliability is NA or
#' the 1e-6 placeholder) estimate a proxy from the standard errors:
#' Vg = var(predictedValue) - mean(stdError^2), r2 = 1 - stdError^2 / Vg;
#' (3) otherwise return 1 (no penalty) so the index is not zeroed.
#'
#' @param preds long data frame with designation, trait, predictedValue and
#'   optionally reliability and stdError.
#' @param designations character vector, row order of the result.
#' @param traits character vector, column order of the result.
#' @param placeholder reliabilities at or below this value count as "not estimated".
#' @return numeric matrix with attribute "source" (reliability / proxy / none per trait).
#' @noRd
pam_reliability_matrix <- function(preds, designations, traits, placeholder = 1e-5) {
  m <- matrix(1, nrow = length(designations), ncol = length(traits),
              dimnames = list(designations, traits))
  src <- stats::setNames(rep("none", length(traits)), traits)
  has_rel <- "reliability" %in% colnames(preds)
  has_se <- all(c("stdError", "predictedValue") %in% colnames(preds))
  for (t in traits) {
    d <- preds[preds$trait == t, , drop = FALSE]
    d <- d[!duplicated(d$designation), , drop = FALSE]
    idx <- match(designations, d$designation)
    rel <- if (has_rel) suppressWarnings(as.numeric(d$reliability[idx])) else rep(NA_real_, length(designations))
    if (any(is.finite(rel) & rel > placeholder)) {
      rel[!is.finite(rel)] <- 0
      m[, t] <- pmax(0, pmin(1, rel))
      src[t] <- "reliability"
      next
    }
    if (has_se) {
      se <- suppressWarnings(as.numeric(d$stdError[idx]))
      pv <- suppressWarnings(as.numeric(d$predictedValue[idx]))
      ok <- is.finite(se) & is.finite(pv)
      if (sum(ok) >= 3) {
        vg <- stats::var(pv[ok]) - mean(se[ok]^2)
        if (is.finite(vg) && vg > 0) {
          r <- 1 - se^2 / vg
          r[!is.finite(r)] <- 0
          m[, t] <- pmax(0, pmin(1, r))
          src[t] <- "proxy"
        }
      }
    }
  }
  attr(m, "source") <- src
  m
}

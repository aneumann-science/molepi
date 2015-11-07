# R functions to calculate environemtnal correlation
#
# Formula from Trzaskowski, M., Yang, J., Visscher, P. M., & Plomin, R. (2013). DNA evidence for strong genetic stability and increasing heritability of intelligence from age 7 to 12. Molecular psychiatry, 19(3), 380-384.
#
# Author: Alexander Neumann (a.neumann@erasmusmc.nl)
# Date: 07.11.2015

#' Calculate the residual correlation from a fitted bivariate GREML model
#'
#' \code{write_gpf_syntax} Bivariate GREML as implemented in GCTA outputs
#' genetic correlation but not residual correlation. This function calculates
#' the residual correlation based on the variances and covariances of both
#' traits.
#'
#' Formula from Trzaskowski, M., Yang, J., Visscher, P. M., & Plomin, R. (2013).
#' DNA evidence for strong genetic stability and increasing heritability
#' of intelligence from age 7 to 12. Molecular psychiatry, 19(3), 380-384.
#'
#' @param Ve1 A numeric value. The residual variance of trait 1.
#' @param Ve2 A numeric value. The residual variance of trait 2.
#' @param Ce A numeric value. The residual covariance between
#'        trait 1 and 2.
#'
#' @return Residual correlation between trait 1 and 2
res_cor <- function(Ve1, Ve2, Ce) {
  return(Ce/(sqrt(Ve1)*sqrt(Ve2)))
}

#' Calculate the standard error for the
#' residual correlation from a fitted bivariate GREML model
#'
#' \code{write_gpf_syntax} Bivariate GREML as implemented in GCTA outputs
#' genetic correlation but not residual correlation. This function calculates
#' the standard error of the residual correlation based on the variances and
#' covariances of both traits. See covariance matrix of GCTA output.
#'
#' Formula from Trzaskowski, M., Yang, J., Visscher, P. M., & Plomin, R. (2013).
#' DNA evidence for strong genetic stability and increasing heritability
#' of intelligence from age 7 to 12. Molecular psychiatry, 19(3), 380-384.
#'
#' @param Ve1 A numeric value. The residual variance of trait 1.
#' @param VarVe1 A numeric value. Sampling variance of trait 1.
#'        Can be found in 1,1 of covariance matrix.
#' @param Ve2 A numeric value. The residual variance of trait 2.
#' @param VarVe2 A numeric value.
#'        The sampling variance of trait 2.
#'        Can be found in 2,2 of covariance matrix.
#' @param Ce A numeric value. The residual covariance between
#'        trait 1 and 2.
#' @param VarCe A numeric value. Sampling Variance of Ce [square(SE of C(e)_tr2)].
#'        Can be found in 3,3 of covariance matrix.
#' @param CovVe1Ve2 A numeric value. Sampling covariance between Ve1 and Ve2.
#'        Can be found in 5,4 in Covariance Matrix
#' @param CovVe1Ce A numeric value. Sampling covariance between Ve1 and Ce.
#'        Can be found in 6,4 in Covariance Matrix.
#' @param CovVe2Ce A numeric value. Sampling covariance between Ve2 and Ce.
#'        Can be found in 6,5 in Covariance Matrix
#'
#' @return Standard error of residual correlation between trait 1 and 2
res_se <- function(Ve1, VarVe1, Ve2, VarVe2, Ce, VarCe, CovVe1Ve2, CovVe1Ce, CovVe2Ce) {
  Re <- res_cor(Ve1, Ve2, Ce)
  VarRe <- Re*Re*((VarVe1)/(4*Ve1*Ve1)+(VarVe2)/(4*Ve2*Ve2) + (VarCe)/(Ce*Ce) + CovVe1Ve2/(2*Ve1*Ve2)-CovVe1Ce/(Ve1*Ce)
                  -CovVe2Ce/(Ve2*Ce))
  return(sqrt(VarRe))
}

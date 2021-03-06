% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lss_calc.R
\name{lss_calc}
\alias{lss_calc}
\title{Calculations for SN}
\usage{
lss_calc(x)
}
\arguments{
\item{x}{Object of class lssVarReg (output from \code{lssVarReg}).}
}
\value{
dataframe containing:
\itemize{
\item\code{y}: y variable
\item\code{x}: x variable
\item\code{eta}: \eqn{\eta} or fitted location estimates
\item\code{omega}: \eqn{\omega} or fitted scale estimates
\item\code{shape}: \eqn{\alpha} or fitted shape estimates
 \item\code{predicted mean}: fitted mean estimates
 \item\code{predicted variance}: fitted variance estimates
 \item\code{Predicted skewness}: fitted skewness estimates
 \item\code{stand.res2}: Squared standardised residuals
}
}
\description{
\code{lss_calc} performs calculations for transforming SN data (location, scale and shape) to mean, variance and skew. This function is utilised by other, higher level functions.
}

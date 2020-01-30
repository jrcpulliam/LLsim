# simHelpers.R

#-------------#
#-- Helpers --#
#-------------#

#' Generate IDs for all population members
#'
#' \code{generateIDs} creates random IDs for each individual in the simulated
#' population. This function is used internally by \code{simpleSim}.
#'
#' @return Returns a random 5-letter string.
GenerateIDs <- function() paste(sample(letters,5,T),collapse = "") 

#' Determine incubation period
#'
#' \code{get_inc} selects \code{number} random values for the incubation period
#' from a Gamma distribution with the specified shape and rate parameters.
#'
#' Currently, \code{get_inc} is a wrapper for the \code{rgamma} function from
#' the \code{stats} package. The intent is to implement this more generically in
#' the future to allow use of distributions other than Gamma.
#'
#' @param number The number of random number draws to be made (integer).
#' @param shape The shape parameter of the Gamma distribution to draw from
#'   (numeric).
#' @param rate The rate parameter of the Gamma distribution to draw from
#'   (numeric).
get_inc <- function(number, shape, rate){rgamma(number,shape,rate)} #@ Allow non-gamma options

#' Determine infectious period
#'
#' \code{get_inf} selects \code{number} random values for the infectious period
#' from a Gamma distribution with the specified shape and rate parameters.
#'
#' Currently, \code{get_inf} is a wrapper for the \code{rgamma} function from
#' the \code{stats} package. The intent is to implement this more generically in
#' the future to allow use of distributions other than Gamma.
#'
#' @param number The number of random number draws to be made (integer).
#' @param shape The shape parameter of the Gamma distribution to draw from
#'   (numeric).
#' @param rate The rate parameter of the Gamma distribution to draw from
#'   (numeric).
#'   
get_inf <- function(number, shape, rate){rgamma(number,shape,rate)} #@ Allow non-gamma options

#' Calculate total transmission rate
#'
#' \code{calcTransRate} calculates the total transmission rate at a given point
#' in time
#'
#' @param beta The transmission coefficient (can be thought of as the per capita
#'   contact rate times the probability of trasnmission given contact between a
#'   susceptible and an infectious individual) (numeric, positive).
#' @param S The number of susceptible individuals in the population at a given
#'   point in time (numeric, positive).
#' @param I The number of infectious individuals in the population at a given
#'   point in time (numeric, positive).
#' @param N The total population size to be used in the calculation of the
#'   transmission rate (numeric, positive). Can be either the population size at
#'   the beginning of the epidemic (as in \code{simpleSim}) or the population
#'   size at the time point of interest (i.e., the initial population size minus
#'   the number of individuals who have died).
#'   
calcTransRate <- function(beta, S, I, N){beta * S * I / N}
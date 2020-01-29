# simHelpers.R

#-------------#
#-- Helpers --#
#-------------#

## Generate IDs for all population members
GenerateIDs <- function() paste(sample(letters,5,T),collapse = "") 

## Determine incubation period
get_inc <- function(number, shape, rate){rgamma(number,shape,rate)} #@ Allow non-gamma options

## Determine infectious period
get_inf <- function(number, shape, rate){rgamma(number,shape,rate)} #@ Allow non-gamma options

## Calculate total transmission rate
calcTransRate <- function(beta, S, I, N){beta * S * I / N}
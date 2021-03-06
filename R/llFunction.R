# llFunctions.R

#' Convert simulation output to line list
#'
#' \code{convertLineList} takes simulation output on cases, converts numeric
#' time to dates, and adds an observation layer to output observed cases on a
#' specified date.
#'
#' TODO: Add more details later
#'
#' @param cases A data frame generated by \code{simpleSim} with rows
#'   representing all cases that occur in the outbreak simulation, having at
#'   minimum the following columns: \code{id} (unique individual ID), 
#' @param obs_prob The individual-level probabilty that a case is reported
#'   (numeric)
#' @param mean_report_delay The average number of days between disease onset and
#'   reporting to a health facility (numeric)
#' @param mean_confirm_delay The average number of days between reporting to a
#'   health facility and confirmation of a case (numeric)
#' @param zero_date Character string representing the date on which the initial
#'   case was infected (YYYY-MM-DD)
#' @param obs_date Character string representing the date on which observation
#'   of the state of the line list takes place (YYYY-MM-DD; default = system
#'   date)
#'
#' @return Returns a data frame the provides the line list, as observed on
#'   \code{obs_date}
#'
#'   TODO: Add more detail
#'
#' @export
#' 
createLineList <- function(cases, obs_prob, mean_report_delay, mean_confirm_delay, zero_date, obs_date = as.character(Sys.Date())){

    zero_date <- as.Date(zero_date, format = '%Y-%m-%d')
    obs_date <- as.Date(obs_date, format = '%Y-%m-%d')
    
    cases$observed <- rbinom(nrow(cases), 1, obs_prob)
    
    cases$reported <- ifelse(cases$observed, cases$doo + rpois(nrow(cases), lambda = mean_report_delay), NA)
    cases$doc <- ifelse(cases$observed, cases$reported + rpois(nrow(cases), lambda = mean_confirm_delay), NA)

    cases$reported <- as.Date(floor(cases$reported), origin = zero_date)
    
    cases$doi <- as.Date(floor(cases$doi), origin = zero_date)
    cases$doo <- as.Date(floor(cases$doo), origin = zero_date)
    
    cases$dod <- ifelse(cases$died == TRUE, cases$dor, NA)
    cases$dod <- as.Date(floor(cases$dod), origin = zero_date)
    
    cases$dor <- as.Date(floor(cases$dor), origin = zero_date) # Date of removal
    
    cases$doc <- as.Date(floor(cases$doc), origin = zero_date)
    cases$status <- ifelse(cases$doc <= obs_date, 'confirmed', 'suspected')
    cases$status[cases$reportDate > cases$deathDate] <- 'probable'
    
    obs <- subset(cases, observed == TRUE, select = c(id, doo, dod, reported, doc, status, source))
    names(obs) <- c('caseID', 'onsetDate', 'deathDate', 'reportDate', 'confirmedDate', 'status', 'source')
    dat <- subset(obs, reportDate <= obs_date)
    dat$source[!dat$source %in% dat$caseID] <- NA
    
    dat$deathDate[dat$deathDate > obs_date] <- NA
    
    return(dat)
}

# llFunctions.R

createLineList <- function(cases, obsProb, meanReportDelay, meanConfirmDelay, zeroDate, obsDate = as.character(Sys.Date()), browse = FALSE){

    if(browse) browser()
    zeroDate <- as.Date(zeroDate, format = '%Y-%m-%d')
    obsDate <- as.Date(obsDate, format = '%Y-%m-%d')
    
    cases$observed <- rbinom(nrow(cases), 1, obsProb)
    
    cases$reported <- ifelse(cases$observed, cases$doo + rpois(nrow(cases), lambda = meanReportDelay), NA)
    cases$reported <- as.Date(floor(cases$reported), origin = zeroDate)
    
    cases$doi <- as.Date(floor(cases$doi), origin = zeroDate)
    cases$doo <- as.Date(floor(cases$doo), origin = zeroDate)
    
    cases$dod <- ifelse(cases$died == TRUE, cases$dor, NA)
    cases$dod <- as.Date(floor(cases$dod), origin = zeroDate)
    
    cases$dor <- as.Date(floor(cases$dor), origin = zeroDate)
    
    cases$doc <- ifelse(cases$observed, cases$reported + rpois(nrow(cases), lambda = meanConfirmDelay), NA)
    cases$doc <- as.Date(floor(cases$doc), origin = zeroDate)
    cases$status <- ifelse(cases$doc <= obsDate, 'confirmed', 'suspected')
    cases$status[cases$reportDate > cases$deathDate] <- 'probable'
    
    obs <- subset(cases, observed == TRUE, select = c(id, doo, dod, reported, status))
    names(obs) <- c('caseID', 'onsetDate', 'deathDate', 'reportDate', 'status')
    dat <- subset(obs, reportDate <= obsDate)
    
    dat$deathDate[dat$deathDate > obsDate] <- NA
    
    return(list(cases = cases, linelist = dat))
}

# llFunctions.R

createLL <- function(cases, full, paramVals, zeroDate, obsDate = '2019-12-15', browse = FALSE){
  with(paramVals,{
    
    if(browse) browser()
    zeroDate <- as.Date(zeroDate, format = '%Y-%m-%d')
    obsDate <- as.Date(obsDate, format = '%Y-%m-%d')
    
    if(full == FALSE){
      cases$reported <- ifelse(cases$observed, cases$doo + rpois(nrow(cases),lambda = meanReportDelay_initial), NA)
      cases$reported <- as.Date(floor(cases$reported), origin = zeroDate)
      
      cases$doi <- as.Date(floor(cases$doi), origin = zeroDate)
      cases$doo <- as.Date(floor(cases$doo), origin = zeroDate)
      
      cases$dod <- ifelse(cases$died == TRUE, cases$dor, NA)
      cases$dod <- as.Date(floor(cases$dod), origin = zeroDate)
      
      cases$dor <- as.Date(floor(cases$dor), origin = zeroDate)
    }
    obs <- subset(cases, observed == TRUE, select = c(id, hh, doo, dod, reported))
    names(obs) <- c('caseID','householdID', 'onsetDate', 'deathDate', 'reportDate')
    dat <- subset(obs, reportDate <= obsDate)
    
    dat$status <- NA
    dat$status <- ifelse((obsDate - dat$reportDate) >= rpois(nrow(dat), lambda = meanConfDelay), 'confirmed', 'suspected') #@ will be messy if second report date close to first or long meanConfDelay!!!
    dat$status[dat$reportDate > dat$deathDate] <- 'probable'
    dat$deathDate[dat$deathDate > obsDate] <- NA
    
    return(list(cases = cases, linelist = dat))
  })
}

## runSim.R
##
## Characteristics: stochastic, continuous time, includes heterogenous offspring
## distribution for secondary cases
##
## Output: line listing (date of onset, date of death if applicable), known 
## exposure events (id, type, time), simulation inputs (including seed and parameters),
## acutal infection tree

runSim <- function(parms, cont, simFun = simpleSim){
  out <- simFun(parms, cont)
  return(out)
}

# simpleSim <- function(parms, cont){
#   # PARMS: inc_mean_param, inf_mean_param, R_0, case_fatality, pop_size,
#   # inc_shape_param = 1, inf_shape_param = 1
#   # CONT: initial_cases = 1, max_time = 365, seed = 20200121, initDate
#   with(c(parms, cont),{

simpleSim <- function(inc_mean_param, inf_mean_param, R_0, case_fatality, pop_size, inc_shape_param = 1, inf_shape_param = 1, initial_cases = 1, max_time = 365, seed = 20200129, browse = FALSE){  

  set.seed(seed)
  
  # Define derived parameters
  
  inc_rate_param <- inc_shape_param / inc_mean_param
  inf_rate_param <- inf_shape_param / inf_mean_param
  
  beta <- R_0 / inf_mean_param
  
  ###------------###
  ### Initialize ###
  ###------------###
  
  ## Create population
  pop <- data.frame(
    id = replicate(n = pop_size, GenerateIDs())
    , doi = NA
    , doo = NA
    , dor = NA
    , ill = FALSE
    , died = 0
    , source = NA
  )
  
  ## Ensure indvidual identifiers are unique
  while(sum(duplicated(pop$id)) > 0){
    new.ids <- replicate(n = sum(duplicated(pop$id)), GenerateIDs())
    levels(pop$id) <- c(levels(pop$id),new.ids) # Why is this a factor?
    pop$id[duplicated(pop$id)] <- new.ids
  } #@ test_that IDs are unique
  
  # Initialize cases
  tmp <- sample(1:pop_size, initial_cases, replace = FALSE)
  pop$doi[tmp] <- 0
  pop$doo[tmp] <- get_inc(initial_cases, inc_shape_param, inc_rate_param)
  pop$dor[tmp] <- pop$doo[tmp] + get_inf(initial_cases, inf_shape_param, inf_rate_param)
  pop$ill[tmp] <- TRUE
  
  ## Initialize simulation
  time <- pop$doo[tmp]
  if(browse) browser()
  
  ## Simulation
  while(time < max_time & any(is.na(pop$dor))){
    tot_trans_rate <- calcTransRate(beta, S = sum(is.na(pop$doi)), I = sum(pop$ill), N = pop_size)
    time_next_trans_event <- time + ifelse(tot_trans_rate == 0, NA, rexp(1,tot_trans_rate))
    next_transmitter <- ifelse(tot_trans_rate == 0, NA, 
                               ifelse(sum(pop$ill) == 1, which(pop$ill == TRUE),
                                      sample(which(pop$ill == TRUE), 1, FALSE)))
    time_next_symp_onset <- ifelse(all(is.na(pop$doo)|pop$doo <= time), Inf, min(subset(pop, doo>time)$doo, na.rm = T))
    time_next_removal <- ifelse(all(is.na(pop$dor)|pop$dor <= time), Inf, min(subset(pop, dor>time)$dor, na.rm = T))
    
    ## Determining time and type of next event
    event_times <- c(time_next_trans_event,
                     time_next_symp_onset,
                     time_next_removal)
    event_type <- c("transmission",
                    "symptom_onset",
                    "removal")[which.min(event_times)]
    
    ## Implement next event
    switch(event_type,
           transmission = {
             inf_whom <- ifelse(sum(is.na(pop$doi)) == 1, which(is.na(pop$doi)), sample(which(is.na(pop$doi)), 1, FALSE))
             time <- time_next_trans_event
             
             pop$doi[inf_whom] <- time
             pop$doo[inf_whom] <- time + get_inc(1, inc_shape_param, inc_rate_param)
             pop$dor[inf_whom] <- pop$doo[inf_whom] + get_inf(1, inf_shape_param, inf_rate_param)
             
             ## Track infection source
             pop$source[inf_whom] <- as.character(pop$id[next_transmitter])
           },
           symptom_onset = {
             time <- time_next_symp_onset
             
             pop$ill[which(pop$doo == time_next_symp_onset)] <- TRUE
           },
           removal = {
             who <- which(pop$dor == time_next_removal)
             time <- time_next_removal
             
             pop$ill[who] <- FALSE
             if(runif(1) < case_fatality){
               pop$died[who] <- TRUE
             }
           },
           warning("Unknown event type.")
    )
  }
  
  cases <- subset(pop, !is.na(doi))
  
  return(list(cases = cases, pop = pop))
  # })
}
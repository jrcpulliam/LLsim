## runSim.R
##
## Characteristics: stochastic, continuous time, includes heterogenous offspring
## distribution for secondary cases
##
## Output: line listing (date of onset, date of death if applicable), known 
## exposure events (id, type, time), simulation inputs (including seed and parameters),
## acutal infection tree

runSim <- function(parms, cont, simFun = simpleSim, browse = FALSE){
  out <- simFun(parms, cont)
  return(out)
}

stoopidSim <- function(parms, cont){
  with(c(parms, cont),{
    
  })
}

simpleSim <- function(parms, cont){
  # PARMS: inc_mean_param, inf_mean_param, R_0, case_fatality, pop_size,
  # inc_shape_param = 1, inf_shape_param = 1
  # CONT: initial_cases = 1, max_time = 365, seed = 20200121, initDate
  with(c(parms, cont),{
    
    # Define derived parameters
    
    inc_rate_param <- inc_shape_param / inc_mean_param
    inf_rate_param <- inf_shape_param / inf_mean_param
    
    beta <- inf_mean_param * R_0 / pop_size
    
    ###------------###
    ### Initialize ###
    ###------------###
    
    ## Create population
    pop <- data.frame(
      id = replicate(n = pop_size, GenerateIDs())
      , doi = NA
      , doo = NA
      , dor = NA
      , ill = F
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
    pop$ill[tmp] <- 1
    
    ## Initialize simulation
    time <- 0

    ## Simulation
    while(time < MAXTIME){
      tot_trans_rate <- calcTransRate(beta)
      time_next_trans_event <- ifelse(tot_trans_rate == 0, NA, rexp(1,tot_trans_rate))
      next_transmitter <- ifelse(tot_trans_rate == 0, NA, sample(which(pop$ill), 1, FALSE))
      time_next_symp_onset <- ifelse(all(is.na(pop$doo)|pop$doo <= time), Inf, min(subset(pop, doo>time)$doo, na.rm = T))
      time_next_removal <- ifelse(all(is.na(pop$dor)|pop$dor <= time), Inf, min(subset(pop,dor>time)$dor, na.rm = T))
      
      ## Determining time and type of next event
      event_times <- c(time_next_trans_event,
                       time_next_symp_onset,
                       time_next_removal)
      event_type <- c("secondary_transmission",
                      "symptom_onset",
                      "removal")[which.min(event_times)]
      
      ## Implement next event
      switch(event_type,
             secondary_transmission = {
               inf_whom <- sample(which(is.na(pop$doi)), 1, FALSE)
               time <- time_next_trans_event
               
               pop$doi[inf_whom] <- time
               pop$doo[inf_whom] <- time + get_inc(1, inc.shape, inc.rate)
               pop$dor[inf_whom] <- pop$doo[inf_whom] + get_inf(1, inf.shape, inf.rate)

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
               
               population$ill[who] <- FALSE
               if(runif(1) < case_fatality){
                 population$died[who] <- TRUE  #@ output
               }
             },
             warning("Unknown event type.")
      )
    }
    
    cases <- subset(pop, !is.na(doi))
    
    return(list(cases = cases, pop = pop))
  })
}

# unclearSim <- function(inc_mean_param, inf_mean_param, R_0, pop_size = 500, inc_shape_param = 1, inf_shape_param = 1, initial_cases = 1, max_time = 365, seed = 20200121, browse = FALSE){
#   
#   # Define derived parameters
#   
#   inc_rate_param <- inc_shape_param / inc_mean_param
#   inf_rate_param <- inf_shape_param / inf_mean_param
#   
#   beta <- inf_mean_param * R_0 / pop_size
#   
#   ###------------###
#   ### Initialize ###
#   ###------------###
#   
#   ## Create population
#   pop <- data.frame(
#     id = replicate(n = pop_size, GenerateIDs())
#     , doi = NA
#     , doo = NA
#     , dor = NA
#     , ill = F
#     , died = 0
#     , source = NA
#   )
#   
#   ## Ensure indvidual identifiers are unique
#   while(sum(duplicated(pop$id)) > 0){
#     new.ids <- replicate(n = sum(duplicated(pop$id)), GenerateIDs())
#     levels(pop$id) <- c(levels(pop$id),new.ids) # Why is this a factor?
#     pop$id[duplicated(pop$id)] <- new.ids
#   }
#   
#   # Initialize cases
#   tmp <- sample(1:pop_size, initial_cases, replace = FALSE)
#   pop$doi[tmp] <- 0
#   pop$doo[tmp] <- get_inc(initial_cases, inc_shape_param, inc_rate_param)
#   pop$dor[tmp] <- pop$doo[tmp] + get_inf(initial_cases, inf_shape_param, inf_rate_param)
#   pop$ill[tmp] <- 1
#   ################## STOPPED HERE 2020-01-21 16:56 SAST
#   ## First event is zoonotic risk event (may or may not produce infection)
#   ## Calculate zoonotic hazards, time of next zoonotic risk event, and type
#   ## of next zoonotic risk event
#   z.hazards <- UpdateZoonoticHazards(nn = num, population = pop)
#   time.next.risk.event <- 0 # Always first event and occurs at time 0
#   type.next.risk.event <- sample(LETTERS[1:NUM.ZOO.RISKS], 1, FALSE, z.hazards)
#   exposure.event <- NULL
#   
#   ## Calculate h-to-h hazards (0), time of next h-to-h event (Inf),  and source of next h-to-h event (NA)
#   c.hazards <- 0 # UpdateContactHazards(pop)
#   time.next.trans.event <- Inf # time.to.trans.event <- ifelse(sum(c.hazards) == 0,Inf,rexp(1,sum(c.hazards))) 
#   next.transmitter <- NA # ifelse(sum(c.hazards) == 0,NA,sample(which(pop$ill),1,F,apply(c.hazards,2,sum)))
#   
#   ## Intialize time
#   time <- 0
#   ## Initialize counter for zoonotic risk events
#   num.events <- 0
#   
#   ###------------###
#   ### Simulation ###
#   ###------------###
#   
#   # if(browse) browser()
#   
#   while(time < MAXTIME){
#     ## Determining time and type of next event
#     event.times <- c(time.next.risk.event,
#                      time.next.trans.event,
#                      ifelse(all(is.na(pop$doo)|pop$doo <= time), Inf, min(subset(pop,doo>time)$doo,na.rm = T)),
#                      ifelse(all(is.na(pop$dor)|pop$dor <= time), Inf, min(subset(pop,dor>time)$dor,na.rm = T)))
#     event.type <- c("zoonotic.exposure",
#                     "secondary.transmission",
#                     "symptom.onset",
#                     "recovery.or.death")[which.min(event.times)]
#     
#     ## Implement next event
#     switch(event.type,
#            zoonotic.exposure = {
#              num.events <- num.events + 1
#              out <- ZoonoticRiskEvent(type.next.risk.event, num.events, num, prob, pop, time, inc.shape.param, inc.rate.param, inf.shape.param, inf.rate.param, z.hazards)
#              pop <- out$pop
#              time <- time.next.risk.event
#              time.next.risk.event <- out$time.nre
#              type.next.risk.event <- out$type.nre
#              exposure.event <- rbind(exposure.event, out$event)
#            },
#            secondary.transmission = {
#              out <- SecondaryInfection(pop, time, inc.shape.param, inc.rate.param, inf.shape.param, inf.rate.param, contact.hh, contact.neighbor, contact.other, contacts, next.transmitter, c.hazards)
#              pop <- out$pop
#              time <- time.next.trans.event
#              time.next.trans.event <- out$time.nte
#              next.transmitter <- out$next.transmitter
#            },
#            symptom.onset = {
#              out <- BecomeInfectious(pop, event.times, contact.hh, contact.neighbor, contact.other, contacts, next.transmitter)
#              pop <- out$pop
#              time <- out$time
#              time.next.trans.event <- out$time.nte
#              next.transmitter <- out$next.transmitter
#              c.hazards <- out$c.hazards
#            },
#            recovery.or.death = {
#              out <- LeaveInfectious(num, pop, event.times, prob.mort, time.next.risk.event, type.next.risk.event, contact.hh, contact.neighbor, contact.other, contacts, z.hazards, c.hazards)
#              pop <- out$pop
#              time <- out$time
#              time.next.risk.event <- out$time.nre
#              type.next.risk.event <- out$type.nre
#              z.hazards <- out$z.hazards
#              time.next.trans.event <- out$time.nte
#              next.transmitter <- out$next.transmitter
#              c.hazards <- out$c.hazards
#              
#            },
#            warning("Unknown event type.")
#     )
#     # print(time)      
#   }
#   cases <- subset(pop,!is.na(source))
#   cases$observed <- rbinom(nrow(cases), 1, obs.prob)
#   
#   return(list(cases = cases, pop = pop, contacts = contacts, numZooEvents = num.events))
# }

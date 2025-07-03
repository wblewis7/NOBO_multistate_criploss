# R and NIMBLE code for performing a Bayesian multistate survival model estimating
#   survival and crippling loss for populations of northern bobwhite (Colinus 
#   virginianus) in the Southeastern United States.
# From: Cramer, A. S., J. A. Rectenwald, W. B. Lewis, M. Hazelbaker, A. Jackson,
#   H. Lott, G. Beane, D. C. Sisson, and J. A. Martin. Estimating Crippling Loss
#   from Harvest with Multi-State Models.
# Contact: Will Lewis, University of Georgia, wblewis7@gmail.com


# Bobwhite were tracked via radiotelemetry over 120-day periods (November - February) 
#   in three consecutive winters (2021-2022, 2022-2023, 2023-2024) and four study 
#   sites in the Southeastern United States. Birds were tracked to assess survival
#   every few days. Harvest occurred periodically at each site, with records kept
#   of birds killed and which birds were shot at on a given day.
# We analyze capture histories using a multistate survival model (Schaub & Pradel. 2004.
#   Ecology). We classify birds into three latent ecological states:
#   1) Alive
#   2) Dead from harvest
#   3) Dead from natural causes not due to harvest (e.g., predation, exposure)
# Birds transition from 1 to 2 based on the harvest mortality rate (SH) and from
#   1 to 3 based on (1-SH) and the daily natural survival rate (phi).
# We model phi based on additive effects of study site and year. Data were collected
#   at four study sites: Central Florida (CF), Tall Timbers (TT), Livingston Place (LP),
#   and the Albany Quail Project (AQP). Data were also collected over three winters
#   starting in 2021, 2022, and 2023. We modeled effects of TT, LP, AQP, 2022, and 2023
#   on phi, with CF 2021 corresponding to the intercept.
# We assume that SH will decline exponentially with time since hunt (dsh), with 
#   SH being near 0 by 2-week post-hunt. This is indicated with the parameter I,
#   which is given a value of 1 if the day is within 14 days of a hunt and 0 otherwise.
# Observation depends on the ecological state and is classified into 5 categories:
#   1) Observed alive via radiotelemetry.
#   2) Observed dead via natural causes with radiotelemetry.
#   3) Observed harvested on hunt days (dsh=0).
#   4) Observed dead after hunt day via crippling loss from radiotelemetry.
#   5) Not observed (i.e., was not tracked or harvested on a given day)
# Alive birds can either be in observation states 1 or 5 based on theta, the probability
#   of tracking a bird on a given day.
# Birds dead from harvest can be in observation states 2, 3, 4, or 5. They can be in
#   state 5 based on 1-theta. They can be in state 3 based on theta and whether or not
#   it is a hunt day (H=1). They can be in state 4 based on theta, not being a hunt day
#   (H=0), and the probability of finding a crippled bird before it is scavenged (FC,
#   i,e, the probability that correctly assigned crippled). They can be in state 2 based
#   on theta, H=0, and 1-FC (i.e., bird was crippled but was scavenged before found).
# Birds dead from natural causes can either be in observation states 2 or 5 based on theta.



require(nimble)

load("NOBO.multistate.criploss.data.gzip")

str(NOBO.multistate.criploss.data)

# y: a 1707x120 multistate capture history. Rows correspond to individual bobwhite
#   while rows correspond to days of the study period. Values of y correspond to
#   the 5 observation states described above.
# nid: the total number of bobwhite tracked in the data (1707).
# TT: a binary vector denoting if each bird in nind was located at TT.
# LP: a binary vector denoting if each bird in nind was located at LP.
# AQP: a binary vector denoting if each bird in nind was located at AQP.
# Y2022: a binary vector denoting if each bird in nind was tracked during the winter
#   starting in 2022.
# Y2023: a binary vector denoting if each bird in nind was tracked during the winter
#   starting in 2023.
# first: a vector denoting the tracking period (1-120) in which each bird in nind
#   was first released.
# last: a vector denoting the tracking period (1-120) in which each bird in nind
#   was last tracked.
# dsh: a matrix of same dimensions as y giving the number of days since the most
#   recent hunt for each bird (rows) and each day of the season (columns). Values
#   of 0 represent hunt days. Days before any hunt in a given year were given a
#   value of 108. Mortality risk was only calculated when dsh<=14.
# I: a binary matrix of same dimensions as y denoting whether or not each day of
#   the season (columns) for each bird (rows) was within 14 days of the hunt day
#   (i.e., dsh<=14).
# H: a binary matrix of same dimensions as y denoting whether or not each day of
#   the season (columns) for each bird (rows) was a hunt day (i.e., dsh=0).

mod.data <- list(y = NOBO.multistate.criploss.data$y)

mod.constants <- list(nind = NOBO.multistate.criploss.data$nind,
                      TT = NOBO.multistate.criploss.data$TT,
                      LP = NOBO.multistate.criploss.data$LP,
                      AQP = NOBO.multistate.criploss.data$AQP,
                      Y2022 = NOBO.multistate.criploss.data$Y2022,
                      Y2023 = NOBO.multistate.criploss.data$Y2023,
                      first = NOBO.multistate.criploss.data$first,
                      last = NOBO.multistate.criploss.data$last,
                      dsh = NOBO.multistate.criploss.data$dsh,
                      I = NOBO.multistate.criploss.data$I,
                      H = NOBO.multistate.criploss.data$H)

NOBO_multistate_criploss <- nimbleCode({
  
  # Priors
  FC ~ dunif(0,1) # Probability of finding a crippled bird and assigning to the correct observation state
  theta ~ dunif(0,1) # Probability of tracking a bird on a given day
  g0 ~ dunif(0,1) # Baseline value of harvest rate on hunt day
  sigma ~ dunif(0,100) # Informs scale of exponential relationship between dsh and SH
  b.0 ~ dnorm(0, sd=100) # Phi intercept (CF, 2021)
  b.1 ~ dnorm(0, sd=100) # Effect of TT on phi
  b.2 ~ dnorm(0, sd=100) # Effect of LP on phi
  b.3 ~ dnorm(0, sd=100) # Effect of AQP on phi
  b.4 ~ dnorm(0, sd=100) # Effect of the 2022-2023 winter on phi
  b.5 ~ dnorm(0, sd=100) # Effect of the 2023-2024 winter on phi

  
  #   Ecological states are:
  #     1) Alive 
  #     2) Dead from harvest
  #     3) Dead from natural causes unrelated to harvest
  
  #   Observation states are:
  #     1) Detected alive via radiotelemetry
  #     2) Observed dead during telemetry from natural causes
  #     3) Observed harvested on the hunt day
  #     4) Observed dead after the hunt with mortality attributed to harvest
  #     5) Not observed
  
  
  for(i in 1:nind){
    
    # Know alive at first capture period
    z[i, first[i]] <- 1
    
    # Linear model for natural survival
    logit(phi[i]) <- b.0 + b.1 * TT[i] + b.2 * LP[i] + b.3 * AQP[i] + b.4 * Y2022[i] + b.5 * Y2023[i]
    
    # Looping through periods
    for(t in (first[i]+1):last[i]){
      
      # Calculating risk of harvest mortality based on exponential decline since hunt day (14 days max)
      S.H[i,t] <- g0 * exp(-dsh[i,t]/sigma) * I[i,t]
      
      
      # State transition matrix. Transition from state in t (rows) to state in t+1 (columns)
      ps[1,1:3,i,t] <- c((1 - S.H[i,t]) * phi[i],
                         S.H[i,t],
                         (1 - S.H[i,t]) * (1 - phi[i]))
      ps[2,1:3,i,t] <- c(0, 
                         1,
                         0)
      ps[3,1:3,i,t] <- c(0,
                         0,
                         1)
      
      # First order state process
      z[i,t] ~ dcat(ps[z[i,t-1], 1:3, i, t])
      
      
      
      
      # Observation matrix. Observation state (columns) depends on ecological state (rows)
      po[1,1:5,i,t] <- c(theta, 
                         0, 
                         0, 
                         0, 
                         1 - theta)
      po[2,1:5,i,t] <- c(0, 
                         theta * (1 - FC) * (1 - H[i,t]), 
                         theta * H[i,t], 
                         theta * FC * (1 - H[i,t]), 
                         1 - theta)
      po[3,1:5,i,t] <- c(0,
                         theta,
                         0,
                         0,
                         1 - theta) 
      
      # Observation process
      y[i,t] ~ dcat(po[z[i,t], 1:5, i, t])
    }
  }
})


# Setting initial values. Need to ensure that birds area given ecological state of
#   2 (dead harvest) only if found dead and there was recent hunting

znits <- matrix(NA, nrow = nrow(mod.data$y),
                ncol = ncol(mod.data$y))

for(i in 1:nrow(znits)){
  if(mod.data$y[i,mod.constants$last[i]]==1){
    # Censured birds
    znits[i,mod.constants$first[i]:mod.constants$last[i]] <- 1
  } else{
    # For dead birds, filling in up to last known alive with 1s, randomly select
    #   true fate date for any unknown fate dates (between tracking events)
    maxalive <- max(which(mod.data$y[i,]==1))
    znits[i,mod.constants$first[i]:maxalive] <- 1
    if(mod.data$y[i,mod.constants$last[i]]==3){
      # Killed and recovered on hunt
      znits[i,(maxalive+1):mod.constants$last[i]] <- 1
      znits[i,mod.constants$last[i]] <- 2 # Formulation allows for birds which were harvested the data after tracking
    } else if(mod.data$y[i,mod.constants$last[i]]==4){
      # Found crippled on telemetry, can't die until at least day 1 post harvest
      daysest <- mod.constants$I[i,(maxalive+1):mod.constants$last[i]] * mod.constants$dsh[i,(maxalive+1):mod.constants$last[i]]
      daysest_prob <- daysest/sum(daysest)
      whichfate <- sample(1:length(daysest), 1, prob=daysest_prob)
      zvals <- rep(1, times=length(daysest))
      zvals[whichfate:length(zvals)] <- 2
      znits[i,(maxalive+1):mod.constants$last[i]] <- zvals
    } else{
      # Birds found dead on telemetry, could be dead natually or crippled (not on day of harvest)
      daysest <- mod.constants$I[i,(maxalive+1):mod.constants$last[i]] * mod.constants$dsh[i,(maxalive+1):mod.constants$last[i]]
      # Randomly selecting dead either naturally (probability each day in daysest)
      #   or via harvest (based on values of 1 in daysest)
      fatevals <- c(2,3)
      fatetype <- sample(fatevals, 1, prob=c(sum(daysest),length(daysest))/(sum(daysest)+length(daysest)))
      if(fatetype == 2){
        # If assigned as harvested
        daysest_prob <- daysest/sum(daysest)
        whichfate <- sample(1:length(daysest), 1, prob=daysest_prob)
      } else{
        # If assigned as dead via natural mortality
        whichfate <- sample(1:length(daysest), 1)
      }
      zvals <- rep(1, times=length(daysest))
      zvals[whichfate:length(zvals)] <- fatetype
      znits[i,(maxalive+1):mod.constants$last[i]] <- zvals
    }
  }
}

inits <- list(z = znits,
              g0 = 0.2,
              sigma = 1,
              theta = 0.5,
              b.0 = runif(1, 5, 7),
              b.1 = runif(1, -1, 1),
              b.2 = runif(1, -1, 1),
              b.3 = runif(1, -1, 1),
              b.4 = runif(1, -1, 1),
              b.5 = runif(1, -1, 1),
              FC = 0.15)

monitors <- c("FC", "theta", "g0", "sigma", "b.0", "b.1", "b.2", "b.3", "b.4", "b.5", "z")


NOBO.multistate.cirploss.model <- nimbleModel(code = NOBO_multistate_criploss,
                                              data = mod.data,
                                              constants = mod.constants,
                                              inits = inits)

NOBO.multistate.cirploss.mcmc <- nimbleMCMC(model = NOBO.multistate.cirploss.model,
                                            niter = 100000, nchains = 3, nburnin = 4000, thin=10,
                                            samplesAsCodaMCMC=TRUE, monitor=monitors)

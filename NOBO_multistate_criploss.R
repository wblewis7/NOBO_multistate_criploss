# R and NIMBLE code for performing a Bayesian multistate survival model estimating
#   survival and crippling loss for populations of northern bobwhite (Colinus 
#   virginianus) in the Southeastern United States based on radiotelemetry and
#   harvest data.
# From: Cramer, A. S., J. A. Rectenwald, W. B. Lewis, M. A. Hazelbaker, A. L. Jackson,
#   H. Lott, G. Beane, D. C. Sisson, M. R. Kunkel, N. M. Nemeth, and J. A. Martin.
#   Estimating crippling loss from harvest with multistate models: A case study on
#   northern bobwhites.
# Contact: Will Lewis, University of Georgia, wblewis7@gmail.com


# Bobwhite were tracked via radiotelemetry over 120-day periods (November 1  - February 28) 
#   in three consecutive winters (2021-2022, 2022-2023, 2023-2024) and four study 
#   sites in the Southeastern United States. Birds were tracked to assess survival
#   every few days. Harvest occurred periodically at each site, with records kept
#   of birds killed and which birds were shot at on a given day. Birds were tracked 
#   from tagging until either mortality, tag failure, or the end of the study season.
#   Tracking did not occur on hunt days for a given bird, so birds could only be 
#   detected on hunt days if shot and recovered. We assume 100% detectability with transmitters. 
# We analyze daily capture histories using a multistate survival model (Schaub & Pradel.
#   2004. Ecology). We classify birds into three latent ecological states (z):
#   1) Alive
#   2) Dead from harvest
#   3) Dead from natural causes not due to harvest (e.g., predation, exposure)
# Birds transition from 1 to 2 based on the daily harvest mortality rate (theta) and from
#   1 to 3 based on (1-theta) and the daily natural survival rate (phi).
# We model phi based on additive effects of study site and year. Data were collected
#   at four study sites: Central Florida (CFL), Tall Timbers (TT), Livingston Place (LP),
#   and the Albany Quail Project (AQP). Data were also collected over three winters
#   starting in 2021, 2022, and 2023. We modeled site-specific intercepts and 
#   effects of the winters starting in 2022 and 2023 on phi, with 2021 corresponding
#   to the intercept.
# We assume that theta will decline exponentially with time since hunt (D), with 
#   theta being near 0 by 2-weeks post-hunt. This is indicated with the parameter I,
#   which is given a value of 1 if the day is within 14 days of a hunt and 0 otherwise.
# Observation depends on the ecological state and is classified into 4 categories:
#   1) Observed alive via radiotelemetry.
#   2) Observed dead via natural causes with radiotelemetry.
#   3) Observed harvested on hunt days (D=0).
#   4) Observed dead after hunt day via crippling loss from radiotelemetry.
# Alive birds can only be observed in observation state 1.
# Birds dead from harvest can be observed in observation states 2, 3, or 4. 
#   On hunt days (denoted with binary variable H, with H=1 representing hunt days)
#   birds can only be observed in observation state 3. On non-hunt days (H=0),
#   birds can be in either observation states 2 or 4 based on the probability (psi) of
#   a carcass being correctly identified as dead from crippling rather than from
#   natural causes (i.e., not scavenged before found by observer).
# Birds dead from natural causes can only be observed in observation state 2.
# The latent ecological state is estimated each day between bird tagging/release
#   and either found dead or censoring; however, birds were only observed every
#   few days during this period. We use an indexing variable A to only incorporate
#   data likelihoods for days with observations.



require(nimble)

load("NOBO.multistate.criploss.data.gzip")

str(NOBO.multistate.criploss.data)

# y: A 1707x120 multistate capture history matrix. Rows correspond to individual bobwhite
#   while columns correspond to days of the study period. Values of y correspond to
#   the 4 observation states described above. Columns refer to the 120 days within the
#   study season (November 1 - February 28) in each of the three years. Note that capture
#   histories were generated for each study season (winters of 2021, 2022, 2023) and stacked,
#   so that the first column represents the observation state on the first day of each study
#   season (November 1).
# N: The total number of bobwhite tracked in the dataset (1707).
# site: A vector giving the study site for each bird in N; 1 = CFL, 2 = TT, 3 = LP, 4 = AQP.
# N.site: The number of unique study sites (4).
# year: A vector giving the study year for each bird in N. Integers correspond to the 
#   winters starting in 2021 (1), 2022 (2), and 2023 (3).
# N.year: The number of unique years in the study (3).
# first: A vector denoting the day of the study season (1-120) in which each bird in N
#   was first capture, tagged, and released.
# last: A vector denoting the day of the study season (1-120) in which each bird in N
#   was last tracked (i.e., period in which the bird was found dead, or after
#   which the bird was censored due to radio failure or end of the study season).
# A: A 1707x55 matrix giving the days of the study season in which each bird in N (rows)
#   was observed, starting with first and ending with last. For instance, the tenth
#   bird was tagged on day 1, observed alive again on day 3, observed alive again on
#   day 5, etc., until the 38th observation on day 85 when it was shot and recovered
#   during a hunt.
# N.A: A vector giving the total number of observations for each bird in N. Corresponds
#   to the length of non-NA values for each row in A.
# D: A matrix of same dimensions as y giving, for each day of the season (columns) and
#   each bird in N (rows), the number of days since the most recent hunting event that the
#   bird was exposed to. Values of 0 represent hunt days. Days before the first hunting event
#   for a given bird were denoted with 108; however, these values do not affect the likelihood
#   because mortality risk is only calculated when D<=14 via the indexing variable I, below.
# I: A binary matrix of same dimensions as y denoting whether or not each day of
#   the season (columns) for each bird in N (rows) was within 14 days of the most recent
#   hunt day (i.e., D<=14).
# H: A binary matrix of same dimensions as y denoting whether or not each day of
#   the season (columns) for each bird in N (rows) was a hunt day (i.e., D=0).





# Multistate Survival Model ----------------------------------------------------
# ------------------------------------------------------------------------------

mod.data <- list(y = NOBO.multistate.criploss.data$y)

mod.constants <- within(NOBO.multistate.criploss.data, rm(y))


# Model
NOBO_multistate_criploss <- nimbleCode({
  
  # Priors
  psi ~ dunif(0,1) # Probability of finding a crippled bird and assigning to the correct observation state (i.e., not scavenged)
  g0 ~ dunif(0,1) # Informs baseline probability of harvest mortality (theta)
  for(q in 1:N.site){
    b.0[q] ~ dnorm(0, sd=100) # Site-specific intercepts of natural survival (logit-scale)
  }
  b.1[1] <- 0
  for(v in 2:N.year){
    b.1[v] ~ dnorm(0, sd=100) # Year effects on natural survival (logit-scale), setting first year (2021) to 0
  }
  sigma ~ dunif(0,100) # Informs scale of relationship between D and theta
  
  
  #   Ecological states are:
  #     1) Alive 
  #     2) Dead from harvest
  #     3) Dead from natural causes unrelated to harvest
  
  #   Observation states are:
  #     1) Detected alive via radiotelemetry
  #     2) Observed dead during telemetry from natural causes
  #     3) Observed harvested on the hunt day
  #     4) Observed dead after the hunt with mortality attributed to harvest
  
  # Daily survival (phi) from natural causes is a factor of site and year
  for(s in 1:N.site){
    for(u in 1:N.year){
      logit(phi[s,u]) <- b.0[s] + b.1[u]
    }
  }
  
  # Ecological State (z). Calculating for each day from first to last.
  for(i in 1:N){
    
    # Known alive at first capture period
    z[i, first[i]] <- 1
    
    # Looping through days
    for(t in (first[i]+1):last[i]){
      
      # Calculating risk of harvest mortality (theta) based on exponential decline since hunt day (D).
      # Enforcing that risk of harvest mortality only occurs within 14 days of the hunt day via the
      #   indicator variable I.
      theta[i,t] <- g0 * exp(-D[i,t] / sigma) * I[i,t]
      
      
      # State transition matrix. Transition from state in t (rows) to state in t+1 (columns).
      # Birds in dead states (2,3) stay in the same state. Alive birds can stay alive based on
      #   not suffering harvest mortality (1-theta) and the daily natural survival rate (phi), can
      #   suffer harvest mortality (theta), or can suffer natural mortality (1-phi) if not 
      #   harvested (1-theta).
      ps[1,1:3,i,t] <- c((1 - theta[i,t]) * phi[site[i],year[i]],
                         theta[i,t],
                         (1 - theta[i,t]) * (1 - phi[site[i],year[i]]))
      ps[2,1:3,i,t] <- c(0, 
                         1,
                         0)
      ps[3,1:3,i,t] <- c(0,
                         0,
                         1)
      
      # First order state process
      z[i,t] ~ dcat(ps[z[i,t-1], 1:3, i, t])
    }
    
    
    # Have ecological states for each day first:last, though birds generally observed every
    #   few days until mortality or censor. Indexing with A to only include data likelihoods
    #   on days when a bird was observed.
    for(a in 2:N.A[i]){
      
      # Observation matrix. Observation state (columns) depends on ecological state (rows).
      # Alive birds can only be observed via telemetry, so can only be observed alive.
      # Birds dead from natural causes can only be observed dead on telemetry from natural causes.
      # On hunt days, only birds dead from harvest can be observed, and only in the dead on hunt
      #   day observation state.
      # On non-hunt days, birds dead from harvest will be found via radiotelemetry but can
      #   be classified correctly as dead from harvest or incorrectly as dead from natural causes
      #   based on psi. This parameter allows for incorrect classification if crippled carcasses 
      #   are scavenged before finding via telemetry.
      po[1,1:4,i,a] <- c(1, 
                         0, 
                         0, 
                         0)
      po[2,1:4,i,a] <- c(0, 
                         (1 - psi) * (1 - H[i,A[i,a]]), 
                         H[i,A[i,a]], 
                         psi * (1 - H[i,A[i,a]]))
      po[3,1:4,i,a] <- c(0,
                         1,
                         0,
                         0) 
      
      # Observation process
      y[i,A[i,a]] ~ dcat(po[z[i,A[i,a]], 1:4, i, a])
      
    }
  }
})


# Initial Values

znits <- matrix(NA, nrow = nrow(mod.data$y),
                ncol = ncol(mod.data$y))

for(i in 1:nrow(znits)){
  if(mod.data$y[i,mod.constants$last[i]]==1){
    # Censored birds
    znits[i,mod.constants$first[i]:mod.constants$last[i]] <- 1
  } else{
    # For dead birds, filling in up to last known alive with 1s, randomly select
    #   true fate date for any days between last known alive date and date found
    #   dead or censored
    maxalive <- max(which(mod.data$y[i,]==1))
    znits[i,mod.constants$first[i]:maxalive] <- 1
    if(mod.data$y[i,mod.constants$last[i]]==3){
      # Killed and recovered on hunt
      znits[i,(maxalive+1):mod.constants$last[i]] <- 1
      znits[i,mod.constants$last[i]] <- 2 # Formulation allows for birds which were harvested the data after tracking
    } else if(mod.data$y[i,mod.constants$last[i]]==4){
      # Found crippled on telemetry, setting mortality at least 1-day post harvest
      daysest <- mod.constants$I[i,(maxalive+1):mod.constants$last[i]] * mod.constants$D[i,(maxalive+1):mod.constants$last[i]]
      daysest_prob <- daysest/sum(daysest)
      whichfate <- sample(1:length(daysest), 1, prob=daysest_prob)
      zvals <- rep(1, times=length(daysest))
      zvals[whichfate:length(zvals)] <- 2
      znits[i,(maxalive+1):mod.constants$last[i]] <- zvals
    } else{
      # Birds found dead on telemetry, could be dead naturally or crippled (not on day of harvest)
      daysest <- mod.constants$I[i,(maxalive+1):mod.constants$last[i]] * mod.constants$D[i,(maxalive+1):mod.constants$last[i]]
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
              b.0 = runif(mod.constants$N.site, 5, 7),
              b.1 = c(NA, runif(mod.constants$N.year-1, -1, 1)),
              psi = 0.15)


monitors <- c("psi", "g0", "sigma", "b.0", "b.1")


NOBO.multistate.criploss.model <- nimbleModel(code = NOBO_multistate_criploss,
                                              data = mod.data,
                                              constants = mod.constants,
                                              inits = inits)

NOBO.multistate.criploss.mcmc <- nimbleMCMC(model = NOBO.multistate.criploss.model,
                                            niter = 100000, nchains = 3, nburnin = 4000, thin=10,
                                            samplesAsCodaMCMC=TRUE, monitor=monitors)
#-------------------------------------------------------------------------------






# Post-Processing --------------------------------------------------------------
#-------------------------------------------------------------------------------
NOBO.mcmc <- rbind(NOBO.multistate.criploss.mcmc[[1]],
                   NOBO.multistate.criploss.mcmc[[2]],
                   NOBO.multistate.criploss.mcmc[[3]])

# Predicting harvest mortality risk from days 0 - 14 post-hunt (theta.pred)
D.range <- 0:14
mcmc.theta <- NOBO.mcmc[,c("g0","sigma")]
theta.pred <- matrix(NA, nrow=nrow(mcmc.theta), ncol=length(D.range))
for(i in 1:nrow(theta.pred)){
  theta.pred[i,] <- mcmc.theta[i,1] * exp(-D.range / mcmc.theta[i,2])
}

# Predicted daily risk of harvest mortality (theta.pred)
apply(theta.pred, 2, quantile, probs=c(0.025,0.5,0.975))

# Cumulative harvest mortality risk (delta)
delta <- 1 - apply(1-theta.pred, 1, prod)
quantile(delta, probs=c(0.025,0.5,0.975))

# Cumulative crippling loss (tau)
tau <- delta - theta.pred[,1]
quantile(tau, probs=c(0.025,0.5,0.975))

# Crippling loss expressed as proportion of total harvest (eta.M.k)
eta.M.k <- tau/delta
quantile(eta.M.k, probs=c(0.025,0.5,0.975))

# Crippling loss expressed as proportion of birds killed and recovered (eta.M.r)
eta.M.r <- tau/theta.pred[,1]
quantile(eta.M.r, probs=c(0.025,0.5,0.975))

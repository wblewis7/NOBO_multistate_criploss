# R and NIMBLE code for performing a Bayesian Binomial regression estimating
#   crippling loss for populations of northern bobwhite (Colinus 
#   virginianus) in the Southeastern United States based on field observations
#   during hunts.
# From: Cramer, A. S., J. A. Rectenwald, W. B. Lewis, M. A. Hazelbaker, A. L. Jackson,
#   H. Lott, G. Beane, D. C. Sisson, M. R. Kunkel, N. M. Nemeth, and J. A. Martin.
#   Estimating crippling loss from harvest with multistate models: A case study on
#   northern bobwhites.
# Contact: Will Lewis, University of Georgia, wblewis7@gmail.com


# Bobwhite were hunted in three consecutive winters (2021-2022, 2022-2023, 2023-2024) 
#   and four study sites in the Southeastern United States: Central Florida (CFL),
#   Tall Timbers (TT), Livingston Place (LP), and the Albany Quail Project (AQP).
# Hunting occurred from mid-November - February. Each site was divided into multiple
#   hunting courses and, on a given day, generally only 1 course was hunted. Courses
#   were not hunted more than 6 times/season, generally with 10-14 days between 
#   hunts at a course.
# On hunts, coveys were located by trained pointing dogs and shot at by hunters
#   with small-bore shotguns. Hunting parties consisted of a dog handler and
#   assistants, no more than two shooters at a time, a wagon driver (when 
#   applicable), and an independent researcher.
# When a covey was located, the researcher and 2-3 members of the hunting party
#   observed the hunt and recorded the number of birds shot and recovered by 
#   hunters (r) and the number of birds shot but unrecovered (c).
# We model the number of unrecovered crippled birds (c) as a Binomial process
#   based on the total number of shot birds n (sum of c and r) and the probability
#   of crippling loss (p).
# We model p based on effects of study site and year.



require(nimble)

load("NOBO.fieldobserved.criploss.data.gzip")

str(NOBO.fieldobserved.criploss.data)

# c: A vector denoting the number of crippled, unrecovered birds on each hunt.
# N.F: The total number of hunts (146).
# r: A vector denoting the number of shot and recovered birds on each hunt.
# site.F: A vector giving the study site for each hunt; 1 = CFL, 2 = TT, 3 = LP, 4 = AQP.
# N.site: The number of unique study sites (4).
# year.F: A vector giving the study year for each hunt. Integers correspond to the 
#   winters starting in 2021 (1), 2022 (2), and 2023 (3).
# N.year: The number of unique years in the study (3).





# Binomial Crippling Loss Model ------------------------------------------------
# ------------------------------------------------------------------------------

mod.data <- list(c = NOBO.fieldobserved.criploss.data$c)

mod.constants <- list(N.F = NOBO.fieldobserved.criploss.data$N.F,
                      n = NOBO.fieldobserved.criploss.data$c + NOBO.fieldobserved.criploss.data$r,
                      site.F = NOBO.fieldobserved.criploss.data$site.F,
                      N.site = NOBO.fieldobserved.criploss.data$N.site,
                      year.F = NOBO.fieldobserved.criploss.data$year.F,
                      N.year = NOBO.fieldobserved.criploss.data$N.year)


# Model
NOBO_fieldobserved_criploss <- nimbleCode({
  
  # Priors
  for(s in 1:N.site){
    a.0[s] ~ dnorm(0, sd=100)  # Site-specific intercepts of crippling loss probability (logit-scale)
  }
  a.1[1] <- 0
  for(y in 2:N.year){
    a.1[y] ~ dnorm(0, sd=100) # Year effects on crippling loss probability (logit-scale), setting first year (2021) to 0
  }
  
  # Likelihood
  for (f in 1:N.F) {
    logit(p[f]) <- a.0[site.F[f]] + a.1[year.F[f]]
    c[f] ~ dbin(p[f], n[f])
  }
  
})


# Initial Values

inits <- function() list(
  a.0 = runif(mod.constants$N.site, -1, 1),
  a.1 = c(NA, runif(mod.constants$N.year-1, -1, 1))
)


monitors <- c("a.0", "a.1")


NOBO.fieldobserved.criploss.model <- nimbleModel(code = NOBO_fieldobserved_criploss,
                                              data = mod.data,
                                              constants = mod.constants,
                                              inits = inits())

NOBO.fieldobserved.criploss.mcmc <- nimbleMCMC(model = NOBO.fieldobserved.criploss.model,
                                            niter = 100000, nchains = 3, nburnin = 4000, thin=10,
                                            samplesAsCodaMCMC=TRUE, monitor=monitors)
#-------------------------------------------------------------------------------






# Post-Processing --------------------------------------------------------------
#-------------------------------------------------------------------------------
NOBO.mcmc.FO <- rbind(NOBO.fieldobserved.criploss.mcmc[[1]],
                      NOBO.fieldobserved.criploss.mcmc[[2]],
                      NOBO.fieldobserved.criploss.mcmc[[3]])

# Predicting crippling loss at each site and year
mat <- as.matrix(data.frame(a0.1 = rep(c(1,0,0,0), each=3),
                            a0.2 = rep(c(0,1,0,0), each=3),
                            a0.3 = rep(c(0,0,1,0), each=3),
                            a0.4 = rep(c(0,0,0,1), each=3),
                            a1.1 = rep(c(1,0,0), times=4),
                            a1.2 = rep(c(0,1,0), times=4),
                            a1.3 = rep(c(0,0,1), times=4)))
p.pred <- plogis(mat %*% t(NOBO.mcmc.FO)) 

# Average crippling loss across sites and years (eta.F.k)
eta.F.k <- colMeans(p.pred)
quantile(eta.F.k, probs=c(0.025,0.5,0.975))

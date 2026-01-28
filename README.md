# NOBO_multistate_criploss
R and NIMBLE code for: Estimating Crippling Loss from Harvest with Multi-State Models.
----
Authors: Amanda S. Cramer, Justin A. Rectenwald, William B. Lewis*, Michael A. Hazelbaker, Alexander L. Jackson, Holly Lott, Geoff Beane, D. Clay Sisson, Melanie R. Kunkel, Nicole M. Nemeth, and James A. Martin
* corresponding author: William B. Lewis, University of Georgia, wblewis7@gmail.com
---

---

# Metadata

Data to estimate crippling loss in northern bobwhite (Colinus virginianus) were collected over 3 consecutive winters (2021-2022, 2022-2023, and 2023-2024) at four sites in the Southeastern United States: Central Florida (CF), Tall Timbers (TT), Livingston Place (LP), and the Albany Quail Project (AQP). Crippling loss was assessed through two methods: 1) a multistate survival model incorporating radiotelemetry and harvest data to estimate survival and crippling loss; and 2) a binomial regression model estimating crippling loss from field observations of th number of birds crippled vs. recovered.

<br>

----

# Multistate Survival Model

We use a multistate surival incorporating radiotelemetry and harvest data to estimate survival, harvest mortality, and crippling loss. Birds were captured in the fall (October - November) and deployed with radiotransmitters, then tracked every few days until either the bird was found dead or was censored due to radio failure or being alive at the end of the 120-day study period on February 28th. Harvest occurred periodically at each site, with records kept of birds killed and which birds were shot on a given day. Birds were not tracked on days in which hunting occurred in their immediate area, so birds could only be detected on hunt days if they were shot and recovered. We assume 100% detectability with transmitters.
For each bird, we estimate the latent ecological state for each day between the first capture date and the last observation, which occurred when a bird was found dead or after which an alive bird was censored. We define three ecological states: 1) alive, 2) dead from harvest, and 3) dead from natural causes.
We model birds as transitioning from state 1 to 2 based on the daily harvest mortality rate (theta), and from 1 to 3 based on surviving harvest (1-theta) and succumbing to natural sources of mortality (1-phi). We model phi based on site-specific intercepts and year effects. We model theta as declining exponentially based on the number of days since the bird was last shot at by hunters (D). We assume that the risk of harvest mortality will decline to near 0 by 2-weeks post-harvest, which we denote with an indicator variable (I) taking the value of 1 if the day is within 2-weeks of the hunt day.
Birds were only observed every few days; therefore, we use an indexing variable A to only include data likelihoods from days with observations. For each day of observations in A, we classify birds into 4 observation states based on the ecological state: 1) observed alive via radiotelemetry, 2) observed dead via natural causes with radiotelemetry, 3) observed harvested on hunt days, and 4) observed dead via radiotelemetry after hunt days from crippling loss. Alive birds can only be observed in observation state 1, while birds dead from natural causes can only be observed in state 2. Birds dead from harvest can be observed in states 2, 3, or 4. On hunt days (denoted with binary variable H, with H=1 representing hunt days), birds can only be observed in observation state 3. On non-hunt days (H=0), birds can be in either observation states 2 or 4 based on the probability (psi) of a carcass being correctly identified as dead from crippling rather than from natural causes (i.e., not scavenged before found by observer).

## NOBO.multistate.criploss.data.gzip
Data for estimating survival, harvest mortality, and crippling loss with the multistate survival model are stored in the file "NOBO.multistate.criploss.data.gzip". 
## y
A 1707x120 multistate capture history matrix. Rows correspond to individual bobwhite while columns correspond to days of the study period. Values of y correspond to the 5 observation states described below.
## nid
The total number of bobwhite tracked in the data (1707).
## TT
A binary vector denoting if each bird in nind was located at TT.
## LP
A binary vector denoting if each bird in nind was located at LP.
## AQP
A binary vector denoting if each bird in nind was located at AQP.
## Y2022
A binary vector denoting if each bird in nind was tracked during the winter starting in 2022.
## Y2023
A binary vector denoting if each bird in nind was tracked during the winter starting in 2023.
## first
A vector denoting the tracking period (1-120) in which each bird in nind was first released.
## last
A vector denoting the tracking period (1-120) in which each bird in nind was last tracked.
## dsh
A matrix of same dimensions as y giving the number of days since the most recent hunt for each bird (rows) and each day of the season (columns). Values of 0 represent hunt days. Days before any hunt in a given year were given a value of 108. Mortality risk was only calculated when dsh<=14.
## I
A binary matrix of same dimensions as y denoting whether or not each day of the season (columns) for each bird (rows) was within 14 days of the hunt day (i.e., dsh<=14).
## H
A binary matrix of same dimensions as y denoting whether or not each day ofthe season (columns) for each bird (rows) was a hunt day (i.e., dsh=0).

<br />
<br />

# NOBO_multistate_criploss.R
R and NIMBLE code for running the Bayesian mulistate model. The model incorporates three latent ecological states: 1) alive, 2) dead from harvest, 3) dead from natural causes not due to harvest (e.g., predation, exposure). Birds transition from ecological state 1 to 2 based on the harvest mortality rate (SH) and from state 1 to 3 based on (1-SH) and the daily natural survival rate (phi). We model phi based on effects of TT, LP, AQP, 2022, and 2023, with CF 2021 corresponding to the intercept. We assume that SH will decline exponentially with time since hunt (dsh), with SH being near 0 by 2-weeks post-hunt. This is indicated with the parameter I, which is given a value of 1 if the day is within 14 days of a hunt and 0 otherwise.
The observation state in y depends on the ecological state and the detection process. We classify observation state into 5 categories: 1) observed alive, 2) observed dead via natural causes, 3) observed harvested on hunt days (dsh=0), 4) observed dead after hunt day attributable to crippling loss, and 5) not observed (i.e., was not tracked or harvested on a given day). Alive birds could either be in observation states 1 or 5 based on theta, the probability of tracking a bird on a given day. Birds dead from harvest could be in observation states 2, 3, 4, or 5. They could be in state 5 based on 1-theta. They could be in state 3 based on theta and whether or not it was a hunt day (H=1). They could be in state 4 based on theta, not being a hunt day (H=0), and the probability of finding a crippled bird before it was scavenged (FC, i,e, the probability that correctly assigned crippled). They could be in state 2 based on theta, H=0, and 1-FC (i.e., bird was crippled but was scavenged before found). Birds dead from natural causes could either be in observation states 2 or 5 based on theta.

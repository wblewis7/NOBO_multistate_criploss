# NOBO_multistate_criploss
R and NIMBLE code for: Estimating Crippling Loss from Harvest with Multi-State Models.
----
Authors: Amanda S. Cramer, Justin A. Rectenwald, William B. Lewis*, Michael Hazelbaker, Alex Jackson, Holly Lott, Geoffry Beane, D. Clay Sisson, and James A. Martin
* corresponding author: William B. Lewis, University of Georgia, wblewis7@gmail.com
---

---

# Metadata

# NOBO.multistate.criploss.data.gzip
Data for estimating survival, harvest mortality, and crippling loss with multistate survival models for populations of northern bobwhite (Colinus virginianus) in the Southeastern United States are stored in the file "NOBO.multistate.criploss.data.gzip". Data were collected at four study sites: Central Florida (CF), Tall Timbers (TT), Livingston Place (LP), and the Albany Quail Project (AQP). Data were also collected over 3 consecutive winters (2021-2022, 2022-2023, and 2023-2024). Birds were captured in the fall (October-Novemeber) and deployed with radiotransmitters. Birds were tracked every few days over 120-day periods (November-February) in each year. Harvest occurred periodically at each site and within each year, with records kept of birds killed and which birds were shot at on a given day.
## y
A 1707x120 multistate capture history. Rows correspond to individual bobwhite while rows correspond to days of the study period. Values of y correspond to the 5 observation states described below.
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
A binary matrix of same dimensions as y denoting whether or not each day ofthe season (columns) for each bird (rows) was a hunt day (i.e., dsh=14).

<br />
<br />

# NOBO_multistate_criploss.R
R and NIMBLE code for running the Bayesian mulistate model. The model incorporates three latent ecological states: 1) alive, 2) dead from harvest, 3) dead from natural causes not due to harvest (e.g., predation, exposure). Birds transition from ecological state 1 to 2 based on the harvest mortality rate (SH) and from state 1 to 3 based on (1-SH) and the daily natural survival rate (phi). We model phi based on effects of TT, LP, AQP, 2022, and 2023, with CF 2021 corresponding to the intercept. We assume that SH will decline exponentially with time since hunt (dsh), with SH being near 0 by 2-week post-hunt. This is indicated with the parameter I, which is given a value of 1 if the day is within 14 days of a hunt and 0 otherwise.
The observation state in y depends on the ecological state and the detection process. We classify observation state into 5 categories: 1) observed alive, 2) observed dead via natural causes, 3) observed harvested on hunt days (dsh=0), 4) observed dead after hunt day attributable to crippling loss, and 5) not observed (i.e., was not tracked or harvested on a given day). Alive birds can either be in observation states 1 or 5 based on theta, the probability of tracking a bird on a given day. Birds dead from harvest can be in observation states 2, 3, 4, or 5. They can be in state 5 based on 1-theta. They can be in state 3 based on theta and whether or not it is a hunt day (H=1). They can be in state 4 based on theta, not being a hunt day (H=0), and the probability of finding a crippled bird before it is scavenged (FC, i,e, the probability that correctly assigned crippled). They can be in state 2 based on theta, H=0, and 1-FC (i.e., bird was crippled but was scavenged before found). Birds dead from natural causes can either be in observation states 2 or 5 based on theta.

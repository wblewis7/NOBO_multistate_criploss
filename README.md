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
<br>

----

# Multistate Survival Model

We use a multistate surival incorporating radiotelemetry and harvest data to estimate survival, harvest mortality, and crippling loss. Birds were captured in the fall (October - November) and deployed with radiotransmitters, then tracked every few days until either the bird was found dead or was censored due to radio failure or being alive at the end of the 120-day study period on February 28th. Harvest occurred periodically at each site, with records kept of birds killed and which birds were shot on a given day. Birds were not tracked on days in which hunting occurred in their immediate area, so birds could only be detected on hunt days if they were shot and recovered. We assume 100% detectability with transmitters. <br />
For each bird, we estimate the latent ecological state for each day between the first capture date and the last observation, which occurred when a bird was found dead or after which an alive bird was censored. <br />
We define three ecological states: 
1. Alive
2. Dead from harvest
3. Dead from natural causes  
We model birds as transitioning from state 1 to 2 based on the daily harvest mortality rate (theta), and from 1 to 3 based on surviving harvest (1-theta) and succumbing to natural sources of mortality (1-phi). We model phi based on site-specific intercepts and year effects. We model theta as declining exponentially based on the number of days since the bird was last shot at by hunters (D). We assume that the risk of harvest mortality will decline to near 0 by 2-weeks post-harvest, which we denote with an indicator variable (I) taking the value of 1 if the day is within 2-weeks of the hunt day. <br />
Birds were only observed every few days; therefore, we use an indexing variable A to only include data likelihoods from days with observations. For each day of observations in A, we classify birds into 4 observation states based on the ecological state: <br />
1. Observed alive via radiotelemetry
2. Observed dead via natural causes with radiotelemetry
3. Observed harvested on hunt days
4. Observed dead via radiotelemetry after hunt days from crippling loss <br />
Alive birds can only be observed in observation state 1, while birds dead from natural causes can only be observed in state 2. Birds dead from harvest can be observed in states 2, 3, or 4. On hunt days (denoted with binary variable H, with H=1 representing hunt days), birds can only be observed in observation state 3. On non-hunt days (H=0), birds can be in either observation states 2 or 4 based on the probability (psi) of a carcass being correctly identified as dead from crippling rather than from natural causes (i.e., not scavenged before found by observer).

## NOBO.multistate.criploss.data.gzip
Data for estimating survival, harvest mortality, and crippling loss with the multistate survival model are stored in the file "NOBO.multistate.criploss.data.gzip". 
### y
A 1707x120 multistate capture history matrix. Rows correspond to individual bobwhite while columns correspond to days of the study period. Values of y correspond to the 4 observation states described above. Columns refer to the 120 days within the study season (November 1 - February 28) in each of the three years. Note that capture histories were generated for each study season (winters of 2021, 2022, 2023) and stacked, so that the first column represents the observation state on the first day of each study season (November 1).
### N
The total number of bobwhite tracked in the dataset (1707).
### site
A vector giving the study site for each bird in N; 1 = CFL, 2 = TT, 3 = LP, 4 = AQP.
### N.site
The number of unique study sites (4).
### year
A vector giving the study year for each bird in N. Integers correspond to the winters starting in 2021 (1), 2022 (2), and 2023 (3).
### N.year
The number of unique years in the study (3).
### first
A vector denoting the day of the study season (1-120) in which each bird in N was first capture, tagged, and released.
### last
A vector denoting the day of the study season (1-120) in which each bird in N was last tracked (i.e., period in which the bird was found dead, or after which the bird was censored due to radio failure or end of the study season).
### A
A 1707x55 matrix giving the days of the study season in which each bird in N (rows) was observed, starting with first and ending with last. For instance, the tenth bird was tagged on day 1, observed alive again on day 3, observed alive again on day 5, etc., until the 38th observation on day 85 when it was shot and recovered during a hunt.
### N.A
A vector giving the total number of observations for each bird in N. Corresponds to the length of non-NA values for each row in A.
### D
A matrix of same dimensions as y giving, for each day of the season (columns) and each bird in N (rows), the number of days since the most recent hunting event that the bird was exposed to. Values of 0 represent hunt days. Days before the first hunting event for a given bird were denoted with 108; however, these values do not affect the likelihood because mortality risk is only calculated when D<=14 via the indexing variable I.
### I
A binary matrix of same dimensions as y denoting whether or not each day of the season (columns) for each bird in N (rows) was within 14 days of the most recent hunt day (i.e., D<=14).
### H
A binary matrix of same dimensions as y denoting whether or not each day of the season (columns) for each bird in N (rows) was a hunt day (i.e., D=0).

<br />

## NOBO_multistate_criploss.R
R and NIMBLE code for running the Bayesian mulistate model and for post-processing to calculate cumulative mortality and rates of crippling loss.


<br>
<br>

----

# Field-Observed Binomial Crippling Loss Model

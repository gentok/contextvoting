% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abm_functions_latest.R
\name{oneset}
\alias{oneset}
\title{All One Stage Simulation in Action}
\usage{
oneset(alike_threshold = 0.2, segregation_radius = 1,
  signalbias_mean = 0.2, signalbias_sd = 0.2, signalbias_size = NULL,
  signalbias_ingroup = 0.5, signaldifficulty = 3,
  updatecommonfactor_sd = 0.1, updatecommonfactor_mean = 0,
  updategroupfactor_sd = 0, updategroupfactor_mean = 0,
  updatepersonalfactor_sd = 0, updatepersonalfactor_mean = 0,
  partyloc_noise = 0, partycapacity_noise = 0.5,
  partyblue_loc_direction = 0, partyred_loc_direction = 0,
  partycapacity_direction_mean = 0, partyupdate_type = "oneshot",
  context_type = "region", contextradius = 5, contexttimetoagg = 2,
  b_idedist = 0.198, b_capacity = 0.716, b_incumbent = 1.526,
  b_pvi_nation = -9.827, b_pvi_local = 5.156, b_pvi_sublocal = 1.663,
  b_pid = 1.215, maxt_segregation = 50, n_elections = 30,
  n_drop_elections = 5, n_voters = 2300, spacelim = c(54, 54),
  sizegroup = c(0.5, 0.5), regionlevel = c(3, 3),
  subregionlevel = c(9, 9), startcommonfactor = NULL,
  startgroupfactor = 0.8, startpersonalfactor = 1,
  groupfactor_type = "constant", knowmean1 = 0.65, knowmean2 = 0.65,
  knowspread = 0.85, partyblue_startloc = -0.8,
  partyred_startloc = 0.8, partystartloc_noise = 0,
  partyblue_startcapacity = 0, partyred_startcapacity = 0,
  partystartcapacity_noise = 0.5, burnin_redvoteprob_sd = 0.05,
  burnin_redvoteprob_cor = 0.19, seedval = NULL, drawplot = FALSE,
  run_type = "all", segregated_soc = NULL)
}
\arguments{
\item{alike_threshold}{1 - tolerance in the segregation game. The default = 0.2;}

\item{segregation_radius}{Search radius to look for neightbors The default = 1; 8 people around}

\item{signalbias_mean}{The average common directional bias (i.e., positive value advantages red 
party and vice versa) in signal. The default is 0.5.}

\item{signalbias_sd}{The SD common directional bias. The bias is drawn from normal distribution each time.
The default is 0.5.}

\item{signalbias_size}{The extent of fixed common directional bias. 
If assigned; mean and sd arguments above will be ignored.}

\item{signalbias_ingroup}{The extent of fixed common ingroup bias (i.e., Evaluate the ingroup candidate better).
The default is 0.}

\item{signaldifficulty}{The extent of noise level in signal. Given the default values of knowledge, 
ideology, and party position, the default difficulty of 3 produces idedist with SD 1.9; Consistent with ANES.}

\item{updatecommonfactor_sd}{The movement in common factor. SD. The default is 0.1.}

\item{updatecommonfactor_mean}{The movement in common factor. Mean. The default is 0.}

\item{updategroupfactor_sd}{The movement in group influence on ideology. SD. The default is 0.}

\item{updategroupfactor_mean}{The movement in group influence on ideology. Mean. The default is 0.}

\item{updatepersonalfactor_sd}{The SD for the movemebnt in personal influence on ideology. SD. The default is 0.}

\item{updatepersonalfactor_mean}{The SD for the movemebnt in personal influence on ideology. Mean. The default is 0.}

\item{partyloc_noise}{Random noise in updating locations of party ideologies. The default is 0.}

\item{partycapacity_noise}{Random noise in updating party capacity. The default is 0.5. 
The average score in ANES is 0.7 but it produces too much movement in voting.}

\item{partyblue_loc_direction}{Direction and expected quantity of the movement in blue party ideology.}

\item{partyred_loc_direction}{Direction and expected quantity of the movement in red party ideology.}

\item{partycapacity_direction_mean}{The average direction and expected quantity of the movement in party capacity. 
(Positive vaue advantages red party and negative value advantages blue party).}

\item{partyupdate_type}{The type of updating procedure. 
If \code{"oneshot"}, the update is always given over \code{startloc} and \code{startcapacity} defined in 
\code{\link{setparty}} at each election. If \code{"online"}, the update is given over the ideology and 
capacity realized in the previous election. If \code{"fromstart"} the update is given over the 
ideology and capacity realized in the first election of the simulation run.}

\item{context_type}{"region" (default) or "network" (currently not used).}

\item{contextradius}{Radius to seach neighbors. Only applied when \code{context_type=="network"} (currently not used).}

\item{contexttimetoagg}{How many previous elections to average. Either 1 or 2 (default).}

\item{b_idedist}{The influence of ideological distance. The default = 0.198, from ANES.}

\item{b_capacity}{The influence of capacity The default = 0.716, from ANES.}

\item{b_incumbent}{The influence if incumebent. The default = 1.526, from ANES (EDV result).}

\item{b_pvi_nation}{The influence of national PVI. The default = -9.827, from ANES (EDV result).}

\item{b_pvi_local}{The influence of regional PVI. The default = 5.156, from CCES pvi_state.}

\item{b_pvi_sublocal}{The influence of sub-regional PVI. The default = 1.663, from CCES pvi_state.}

\item{b_pid}{The influence of group membership. The default = 1.215, from ANES.}

\item{maxt_segregation}{Max time of movement in the segregation game. The default is 50.}

\item{n_elections}{Number of elections in each simulation run. The default is 30.}

\item{n_drop_elections}{Number of first N elections dropped from the quality assessment. The default is 5.}

\item{n_voters}{Interger. Number of agents. The default is 2300.}

\item{spacelim}{Integer vector of length 2. Grid space size 
(number of rows and number of columns, respectively). The default is c(54,54).}

\item{sizegroup}{Numeric Vector. Proportion of agents in each group. 
Length defines the number of groups. Must adds up to 1. The default is c(0.5,0.5).}

\item{regionlevel}{Integer vector of length 2. Level to Split Region (Horizontaly & Vertically).
The default is c(3,3).}

\item{subregionlevel}{Ingeger vector of length 2. Level to Split Sub-Region (Horizontaly & Vertically).
The default is c(9,9).}

\item{startcommonfactor}{The starting value of The common factor in ideology (intercept).
The default is NULL (draw randomly).}

\item{startgroupfactor}{The starting value of the group influence on ideology.
The default is 0.8 (from ANES, skewed slightly towards right in the reality).}

\item{startpersonalfactor}{The starting value of the random personal factor in ideology 
(The SD of the random error to the model). The defauls is 1 (to get ideology distribution with SD 1, from ANES).}

\item{groupfactor_type}{If \code{"constant"}, \code{groupfactor} is the constant coefficient for \code{group}.
If \code{"diverging"}, \code{groupfactor} is the SD for folded normal distribution that randomly 
assigns coeffcient for \code{group}.}

\item{knowmean1}{The lowest average level of knowledge for Group Blue/-1. 0.65 is the overall mean in ANES (stable).}

\item{knowmean2}{The highest average level of knowledge for Group Red/1. 0.65 is the overall mean in ANES (stable).}

\item{knowspread}{The spread (SD for normal distribution) of knowledg within a group.
The default is 0.85 (To get approx 0.25 SD after conversion, referred to ANES).}

\item{partyblue_startloc}{Starting location of blue party ideology. The default is -0.8 (from ANES).}

\item{partyred_startloc}{Starting location of red party ideology. The default is 0.8 (from ANES).}

\item{partystartloc_noise}{Random noise given to the starting locations of party ideologies. The default is 0.}

\item{partyblue_startcapacity}{Starting level of blue party capacity. The default is 0.}

\item{partyred_startcapacity}{Starting level of red party capacity. The default is 0.}

\item{partystartcapacity_noise}{Random noise given to the starting level of party capacities. The default is 0.5.}

\item{burnin_redvoteprob_sd}{The spread for the "burn-in" period (i.e., first two elections) voting behaviors. 
The default is 0.05, referred to SD in ANES National PVI. See details.}

\item{burnin_redvoteprob_cor}{The correlation between voting and group in "burn-in" period voting behaviors. 
the default is 0.19. referred to ANES. See details.}

\item{seedval}{The random number seed.}

\item{drawplot}{If True; draw plots of results}

\item{run_type}{The default is "all" (to run both segregation and election models) 
Can also be "segregation" or "election" (which runs only one of them).}

\item{segregated_soc}{The pre-defined society object. If not NULL and run_type == "election", predefined society is used in simulation.}
}
\value{
List object with society characteristics, including following elements:
\itemize{
 \item \code{grid} The grid of voters (with group membership at each cell)
 \item \code{grid_region} The grid of region identifier at each cell.
 \item \code{grid_subregion} The grid of sub-region identifier at each cell. 
 \item \code{grid_coords} All grid coordinates. 
 \item \code{group} The numeric vector of group memberships 
 \item \code{region_stat} Aggregated group statistics by region. 
 \item \code{subregion_stat} Aggregated group statistics by region.
 \item \code{grid_votes2} The grid of voting decisions in t-2 election.
 \item \code{grid_votes1} The grid of voting decisions in t-1 election.
 \item \code{ideology} Created by \code{\link{setideology}}, may be updated by \code{\link{updateideology}}
 \item \code{knowledge} Created by \code{\link{setknowledge}}
 \item \code{context} Created by \code{\link{getregion}}
 \item \code{signal} Created by \code{\link{getsignal}}
 \item \code{incumbent} Character of incumbent party. Either \code{"blue"} or \code{"red"}.
 \item \code{votes} Created by \code{\link{letsvote}}
 \item \code{quality} Aggregated quality of decisions.
}
}
\description{
Conducting one stage of election.
}
\details{
In the first two elections, the model does not have the sufficient 
number of previous elections to aggregate. If the society does not cotain the information of 
two previous elections, the hypothetical votes in first two election is 
radomly provided by the following procedure.
\itemize{
  \item Draw the reference value of voting for red for each group of voters by 
  multi-variate normal distribution with mean 0 and the covariance with 
  group membership (i.e., -1, or 1) set by \code{burnin_redvoteprob_cor}.
  \item Multiply the reference values by \code{burnin_redvoteprob_sd} and 
  add 0.5 to generate the probability of red vote for each group.
  \item Draw vote of voters in each group by the Bernouli distribution with probability 
  drawin in the previous step. 
}
}

#'
#' Initialize the society
#'
#' @description Initializing the Society for Segregation Model
#'
#' @param number Interger. Number of agents.
#' @param sizegroup Numeric Vector. Proportion of agents in each group. 
#' Length defines the number of groups. Must adds up to 1.
#' @param spacelim Integer vector of length 2. Grid space size 
#' (number of rows and number of columns, respectively). 
#' @param regionlevel Integer vector of length 2. Level to Split Region (Horizontaly & Vertically).
#' @param subregionlevel Ingeger vector of length 2. Level to Split Sub-Region (Horizontaly & Vertically).
#' @param seedval Random seed value for the simulation.
#'
#' @return List object with society characteristics.
#' \itemize{
#'   \item \code{grid} The grid of voters (with group membership at each cell)
#'   \item \code{grid_region} The grid of region identifier at each cell.
#'   \item \code{grid_subregion} The grid of sub-region identifier at each cell. 
#'   \item \code{grid_coords} All grid coordinates. 
#'   \item \code{group} The numeric vector of group memberships 
#'   \item \code{region_stat} Aggregated group statistics by region. 
#'   \item \code{subregion_stat} Aggregated group statistics by region.
#' }
#' 
#' @importFrom graphics abline 
#' @importFrom graphics image
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics title
#' @importFrom stats pnorm qnorm rbinom rnorm runif
#' @importFrom utils tail

#' 
#' @export

initialsoc <- function(number = 2300, # number of agents
                       sizegroup = c(0.5, 0.5), # proportion of agents in each group 
                       spacelim = c(54,54), # Space Size. Need 
                       regionlevel = c(3,3), # Level to Split Region (Horizontaly & Vertically)
                       subregionlevel = c(9,9), # Level to Split Sub-Region (Horizontaly & Vertically)
                       seedval = 356) 
{
  
  if (sum(sizegroup)!=1) stop("sizegroup must adds up to 1")

  # Define Groups
  groupid <- seq(1,length(sizegroup),1)
  groupnumber <- floor(number*sizegroup)
  if (sum(groupnumber) < number) {
    set.seed(seedval-1)
    addid <- sample(groupid, number-sum(groupnumber))
    for (i in groupid) {
      groupnumber[which(groupid==i)] <- 
        groupnumber[which(groupid==i)] + length(addid[addid==i])
    }
    if (number != sum(groupnumber)) stop("Ooups, something is wrong.")
    warning(paste("number cannot be divided exactly proportional to sizegroup. 
                   thus number of voters in each group is approximated"))
  }
  groupset <- cbind(groupid, groupnumber)
  group <- c(rep(0,(spacelim[1]*spacelim[2])-number), 
             unlist(lapply(split(groupset,seq(nrow(groupset))), function(k) rep(k[1],k[2]))))
  # Permute the order
  set.seed(seedval)
  group <- sample(group, spacelim[1]*spacelim[2], replace=F)
  
  # Draw A Distribution of Group in Society
  grid <- matrix(group, nrow=spacelim[1], ncol=spacelim[2])
  
  # All coordinates in the space
  grid_coords <- do.call("rbind",
                         lapply(seq(1,ncol(grid),1), 
                                function(k) cbind(seq(1,nrow(grid),1), rep(k,nrow(grid)))))
  
  # Set Region in the space
  rowdiv <-  floor(spacelim[1]/regionlevel[1])
  set.seed(seedval)
  rowremain <- sample(regionlevel[1],spacelim[1]-regionlevel[1]*rowdiv) 
  rowregion <- sort(c(rep(seq(1,regionlevel[1],1),each=rowdiv),rowremain))
  coldiv <-  floor(spacelim[2]/regionlevel[2])
  set.seed(seedval) # Using the same seed value
  colremain <- sample(regionlevel[2],spacelim[2]-regionlevel[2]*coldiv) 
  colregion <- sort(c(rep(seq(1,regionlevel[2],1),each=coldiv),colremain))
  tmp <- do.call("rbind",lapply(colregion, function(k) cbind(rowregion, rep(k,spacelim[1]))))
  region <- as.numeric(as.factor(paste(tmp[,1],tmp[,2])))
  # Draw A Distribution of Group in Society
  grid_region <- matrix(region, nrow=spacelim[1], ncol=spacelim[2])
  regioncharacter <- function(k) {
    total <- length(region[region==k])
    groupshare <- as.numeric(table(factor(group[region==k],levels=c(0,groupid))))
    res <- c(groupshare,total)
    names(res) <- c("g0",paste0("g",groupid),"total")
    return(res)
  }
  region_stat <- as.data.frame(t(sapply(seq(1,max(region),1), regioncharacter)))
  
  # Set Sub-Region in the space
  rowdiv <-  floor(spacelim[1]/subregionlevel[1])
  set.seed(seedval)
  rowremain <- sample(subregionlevel[1],spacelim[1]-subregionlevel[1]*rowdiv) 
  rowsubregion <- sort(c(rep(seq(1,subregionlevel[1],1),each=rowdiv),rowremain))
  coldiv <-  floor(spacelim[2]/subregionlevel[2])
  set.seed(seedval) # Using the same seed value
  colremain <- sample(subregionlevel[2],spacelim[2]-subregionlevel[2]*coldiv) 
  colsubregion <- sort(c(rep(seq(1,subregionlevel[2],1),each=coldiv),colremain))
  tmp <- do.call("rbind",lapply(colsubregion, function(k) cbind(rowsubregion, rep(k,spacelim[1]))))
  subregion <- as.numeric(as.factor(paste(tmp[,1],tmp[,2])))
  # Draw A Distribution of Group in Society
  grid_subregion <- matrix(subregion, nrow=spacelim[1], ncol=spacelim[2])
  subregioncharacter <- function(k) {
    total <- length(subregion[subregion==k])
    groupshare <- as.numeric(table(factor(group[subregion==k],levels=c(0,groupid))))
    res <- c(groupshare,total)
    names(res) <- c("g0",paste0("g",groupid),"total")
    return(res)
  }
  subregion_stat <- as.data.frame(t(sapply(seq(1,max(subregion),1), subregioncharacter)))

  # Return a list
  return(list(grid=grid, 
              grid_region=grid_region,  
              grid_subregion=grid_subregion,  
              coords=grid_coords, happy_tracker=c(), 
              group=group[group>0], 
              region_stat=region_stat,
              subregion_stat=subregion_stat))
  
}

#' Collecting Neighbor Locations
#'
#' @description Collecting the location of neighbors.
#'
#' @param coords Coordinates of the agent. 
#' @param radius Search radius of neighbor.
#' @param spacelim Integer vector of length 2. Grid space size 
#' (number of rows and number of columns, respectively). 
#' @param includeme If TRUE, collect the agent's own value.
#'
#' @return Coordinates of neightbors
#' 
#' @importFrom arrangements permutations
#'
#' @export

get_neighbors <- function(coords, radius, spacelim, includeme=FALSE) {
  
  # Throws error if radius is too large
  if (radius > min(spacelim)) stop("radius larger than the space")
  
  # Location of Neightbors relative to current position
  searchradius <- seq(-radius,radius,1)
  locperm <- permutations(length(searchradius),2, replace=TRUE)
  nbrs <- t(apply(locperm, 1, function(k) searchradius[k]))
  
  # Not to include the agent-self
  if (includeme==FALSE) nbrs <- nbrs[-which(nbrs[,1]==0 & nbrs[,2]==0),]
  
  # Get neighbors' coordinates
  nbrs <- cbind(nbrs[,1] + coords[1],nbrs[,2] + coords[2])
  
  # Jump to the other side if over the limit
  nbrs[,1] <- ifelse(nbrs[,1] < 1, spacelim[1] + nbrs[,1], 
                     ifelse(nbrs[,1] > spacelim[1], nbrs[,1] - spacelim[1],
                            nbrs[,1]))
  nbrs[,2] <- ifelse(nbrs[,2] < 1, spacelim[2] + nbrs[,2], 
                     ifelse(nbrs[,2] > spacelim[2], nbrs[,2] - spacelim[2],
                            nbrs[,2]))
  
  # Return Neighbors' coordinates
  nbrs
}

#' Collecting Region Ingroup Share
#'
#' @description Recoriding the share of ingroup members in region for each agent.
#'
#' @param group Integer. Group membership of the agent. 
#' @param region Integer. Region the agent resides. 
#' @param rs \code{region_stat} exported in the list generated by \code{\link{initialsoc}} function.
#'
#' @return Numeric. The share of ingroup members in the region.
#'
#' @export

collect_region <- function(group, region, rs) {
  
  # Identify the location of person
  g <- group
  r <- region
  # g <- allvalues[allvalues > 0]
  # r <- region[allvalues > 0]
  
  # Give Ingroup stats for individuals
  gr <- data.frame(g=g, r=r)
  gr <- split(gr, seq(nrow(gr)))
  getingroupstatus <- function(k) {
    num <- rs[,paste0("g",as.numeric(k[1]))][as.numeric(k[2])]
    num <- ifelse(is.na(num),0,num)
    dem <- rs[,"total"][as.numeric(k[2])]
    num/dem
  }
  ingroupshare <- sapply(gr, getingroupstatus)
  ingroupshare <- unlist(ingroupshare)
  
  # Return ingroup share
  ingroupshare
  
}

#' Collecting Sub-Region Ingroup Share
#'
#' @description Recoriding the share of ingroup members in sub-region for each agent.
#'
#' @param group Integer. Group membership of the agent. 
#' @param subregion Integer. Sub-region the agent resides. 
#' @param rs \code{subregion_stat} exported in the list generated by \code{\link{initialsoc}} function.
#'
#' @return Numeric. The share of ingroup members in the region.
#'
#' @export

collect_subregion <- function(group, subregion, rs) {
  
  # Identify the location of person
  g <- group
  r <- subregion
  # g <- allvalues[allvalues > 0]
  # r <- subregion[allvalues > 0]
  
  # Give Ingroup stats for individuals
  gr <- data.frame(g=g, r=r)
  gr <- split(gr, seq(nrow(gr)))
  getingroupstatus <- function(k) {
    num <- rs[,paste0("g",as.numeric(k[1]))][as.numeric(k[2])]
    num <- ifelse(is.na(num),0,num)
    dem <- rs[,"total"][as.numeric(k[2])]
    num/dem
  }
  ingroupshare <- sapply(gr, getingroupstatus)
  ingroupshare <- unlist(ingroupshare)
  
  # Return ingroup share
  ingroupshare
  
}

#' Stage Assessment to Aggregate Happiness by Neighborhood
#'
#' @description Aggregating the happiness status of agents in the society.
#'
#' @param soc The list of society characteristics initiated by \code{\link{initialsoc}}.
#' @param radius Integer. Search radius to look for neighbors.  
#' @param alike_threshold 1 - tolerance in the segregation model.
#' @param initiate Default is \code{FALSE}. If \code{TRUE}, do not produce \code{happy_tracker}.
#'
#' @return The list of society characteristics with the happiness measures. 
#' \code{happy_indicators} records the proportion of happy agents and
#' \code{happy_tracker} indicates the history of \code{happy_indicators}.
#'
#' @export

happysoc <- function(soc, 
                   radius, 
                   alike_threshold,
                   initiate=FALSE) {
  
  # Collect values from all cells
  allvalues <- as.numeric(soc$grid)
  # Coordinates with a person (value > 0)
  active_cells <- soc$coords[allvalues > 0,]
  
  # Collecting My and Neighbors' values 
  collect_values <- function(coords) {
    nbrs <- get_neighbors(coords=coords, radius=radius, spacelim=dim(soc$grid))
    nbrs_values <- apply(nbrs, 1, function(k) soc$grid[k[1],k[2]])
    list(me=soc$grid[coords[1],coords[2]], nbrs=nbrs_values)
  }
  cxtvalues <- apply(active_cells, 1, collect_values)

  # Assess if Happy or not
  happy <- function(cxt) {
    amhappy <- 1
    # If somebody is in radius
    if (sum(cxt$nbrs)>0) {
      # Happy if proportion of "alike" neighbors go above threshold
      amhappy <- ifelse(sum(cxt$me==cxt$nbrs)/sum(cxt$nbrs>0) >= alike_threshold, 1, 0)
    } 
    return(amhappy)
  }
  
  # Proportion of happy individuals
  soc$happy_indicators <- sapply(cxtvalues, happy)
  if (initiate==FALSE) {
    soc$happy_tracker <- append(soc$happy_tracker,sum(soc$happy_indicators)/length(soc$happy_indicators))
  }

  return(soc)
}

#' Stage Assessment to Aggregate Happiness by Region
#'
#' @description Aggregating in-group share status of agents by region 
#'
#' @param soc The list of society characteristics initiated by \code{\link{initialsoc}}. 
#'
#' @return The list of society characteristics with the regional in-group share measures.
#' Update \code{region_stat} and generate \code{ingroupshare} and \code{grid_ingroupshare}. 
#'
#' @export

# Stage Assessment to Aggregate Happiness by Region
happyregion <- function(soc) {
  
  # Collect values from all cells
  allvalues <- as.numeric(soc$grid)
  region <- as.numeric(soc$grid_region)
  # Subsets
  r <- region[allvalues > 0]
  g <- allvalues[allvalues > 0]
  # Coordinates with a person (value > 0)
  active_cells <- soc$coords[allvalues > 0,]
  
  # Collect Region Stats
  regioncharacter <- function(k) {
    total <- length(r[r==k])
    groupn <- sapply(seq(1,max(g),1), function(j) length(g[g==j & r==k]))
    res <- as.numeric(c(k, groupn, total))
    return(res)
  }
  rs <- as.data.frame(t(sapply(seq(1,max(r),1), regioncharacter)))
  colnames(rs) <- c("region",paste0("g",seq(1,max(g),1)),"total")
  
  # Collecting My and Neighbors' values
  cxtshare <- collect_region(g,r,rs)
  
  # Ingroup Share
  soc$region_stat <- rs
  soc$ingroupshare <- cxtshare
  allcxts <- rep(NA, length(allvalues))
  allcxts[allvalues>0] <- cxtshare
  soc$grid_ingroupshare <- matrix(allcxts, nrow=dim(soc$grid)[1], ncol=dim(soc$grid)[2])
  
  return(soc)
}

#' Stage Assessment to Aggregate Happiness by Sub-Region
#'
#' @description Aggregating in-group share status of agents by sub-region 
#'
#' @param soc The list of society characteristics initiated by \code{\link{initialsoc}} 
#'
#' @return The list of society characteristics with the regional in-group share measures.
#' Update \code{subregion_stat} and generate \code{ingroupshare_sub} and \code{grid_ingroupshare_sub}. 
#'
#' @export

happysubregion <- function(soc) {
  
  # Collect values from all cells
  allvalues <- as.numeric(soc$grid)
  subregion <- as.numeric(soc$grid_subregion)
  # Subsets
  r <- subregion[allvalues > 0]
  g <- allvalues[allvalues > 0]
  # Coordinates with a person (value > 0)
  active_cells <- soc$coords[allvalues > 0,]
  
  # Collect subregion Stats
  subregioncharacter <- function(k) {
    total <- length(r[r==k])
    groupn <- sapply(seq(1,max(g),1), function(j) length(g[g==j & r==k]))
    res <- as.numeric(c(k, groupn, total))
    return(res)
  }
  rs <- as.data.frame(t(sapply(seq(1,max(r),1), subregioncharacter)))
  colnames(rs) <- c("subregion",paste0("g",seq(1,max(g),1)),"total")
  
  # Collecting My and Neighbors' values
  cxtshare <- collect_subregion(g,r,rs)
  
  # Ingroup Share
  soc$subregion_stat <- rs
  soc$ingroupshare_sub <- cxtshare
  allcxts <- rep(NA, length(allvalues))
  allcxts[allvalues>0] <- cxtshare
  soc$grid_ingroupshare_sub <- matrix(allcxts, nrow=dim(soc$grid)[1], ncol=dim(soc$grid)[2])
  
  return(soc)
}

#' Move unhappy agents
#'
#' @description Move unhappy agents from the current cell.
#'
#' @param soc The list of society characteristics initiated by \code{\link{initialsoc}}.
#' @param seedval The seeding value to control randomization. 
#'
#' @return The list of society characteristics with the updated locations of agents. 
#'
#' @export

# Move unhappy agents 
moveunhappy <- function(soc, seedval=NULL) {

  # Collect values from all cells
  allvalues <- as.numeric(soc$grid)
  
  # Location of unhappy cells
  unhappy_loc <- which(allvalues>0)[soc$happy_indicators==0]

  if (length(unhappy_loc)>0) {
    
    # Identify Available Locations
    inactive_loc <- which(allvalues==0)
    available_loc <- c(inactive_loc, unhappy_loc)
    
    # Draw new locations to move
    set.seed(seedval)
    moveto_loc <- sample(length(available_loc))[1:length(unhappy_loc)]
    moveto_loc <- available_loc[moveto_loc]
    
    # Move unhappy agents
    newvalues <- allvalues
    newvalues[unhappy_loc] <- 0 # Clean up original locations
    newvalues[moveto_loc] <- allvalues[unhappy_loc] # Add new values

    # Return new society
    soc$grid <-matrix(newvalues, 
                      nrow=nrow(soc$grid), ncol=ncol(soc$grid))
    soc$group <- newvalues[newvalues>0]

  } 

  return(soc)
}

#' Set Continuous Ideology
#'
#' @description Initialize ideology of agents with the linear model.
#'
#' @param group The numeric vector of the group memberships of the agents.
#' @param commonfactor The common factor in ideology (intercept).
#' @param groupfactor The group influence on ideology.
#' @param personalfactor The random personal factor in  ideology (The SD of the random error to the model).
#' @param groupfactor_type If \code{"constant"}, \code{groupfactor} is the constant coefficient for \code{group}.
#' If \code{"diverging"}, \code{groupfactor} is the SD for folded normal distribution that randomly 
#' assigns coeffcient for \code{group}.
#' @param b_seed The random number seed to draw group influence used when \code{groupfactor_type=="constant"}.
#' @param e_seed The random number seed to draw personal influence. 
#'
#' @return The list of following elements: 
#' \itemize{
#'   \item \code{y} The agents' ideologies/
#'   \item \code{x} The group values rescaled in -1 to 1 range.
#'   \item \code{a} The common influence = \code{commmonfactor}
#'   \item \code{b} The vector of the group influcence
#'   \item \code{e} The vector of the personal influence
#'   \item \code{tracker} The tracker table of ideology.
#' }
#' 
#' @details The ideology for agent $i$ $(y_i)$ is calculated by the linear function: 
#' $y = a + b_i * x_i + e_i$.
#' 
#' @examples 
#' 
#' set.seed(30)
#' x <- sample(rep(c(-1,1),each=1000))
#' table(x)
#' 
#' # Ideology Highly Sorted by Group
#' ide <- setideology(x, 0, 1.5, 1)
#' par(mfrow=c(1,3),mar=c(2,2,2,2))
#' hist(ide$y, xlim=c(-6,6)); mean(ide$y); sd(ide$y)
#' hist(ide$y[x==-1],xlim=c(-6,6)); sd(ide$y[x==-1])
#' hist(ide$y[x==1],xlim=c(-6,6)); sd(ide$y[x==1])
#' 
#' # Ideology Moderately Sorted by Group
#' ide <- setideology(x, 0, 1, 1)
#' par(mfrow=c(1,3),mar=c(2,2,2,2))
#' hist(ide$y, xlim=c(-6,6)); mean(ide$y); sd(ide$y)
#' hist(ide$y[x==-1],xlim=c(-6,6)); sd(ide$y[x==-1])
#' hist(ide$y[x==1],xlim=c(-6,6)); sd(ide$y[x==1])
#' 
#' # Ideology Weakly Sorted by Group
#' ide <- setideology(x, 0, 0.8, 1)
#' par(mfrow=c(1,3),mar=c(2,2,2,2))
#' hist(ide$y, xlim=c(-6,6)); mean(ide$y); sd(ide$y)
#' hist(ide$y[x==-1],xlim=c(-6,6)); sd(ide$y[x==-1])
#' hist(ide$y[x==1],xlim=c(-6,6)); sd(ide$y[x==1])
#' 
#' @export

setideology <- function(group, commonfactor=0, groupfactor=1, personalfactor=1, 
                        groupfactor_type="constant", # or diverging
                        b_seed=NULL, e_seed=NULL) {
  
  # Rescale group to -1 to 1 scale
  x <- (((group - min(group))/(max(group)-min(group))) - 0.5)*2
  
  # Set and store random seed values (if not set)
  if (is.null(b_seed)) b_seed <- sample(-2^15:2^15, 1)
  if (is.null(e_seed)) e_seed <- sample(-2^15:2^15, 1)
  
  # Common Factor
  a <- commonfactor
  # Group Factor
  if (groupfactor_type=="constant") {
    b_seed <- NA
    b <- groupfactor
  } else if (groupfactor_type=="diverging") {
    set.seed(b_seed)
    b <- abs(rnorm(length(x),sd=abs(groupfactor)))
  }
  # Personal Factor
  set.seed(e_seed)
  e <- rnorm(length(x),sd=personalfactor)
  
  # Generate Ideology
  y <- a + b*x + e
  
  # Store Parameters
  param <- cbind(0,commonfactor,groupfactor,personalfactor,
                 b_seed,e_seed)
  colnames(param) <- c("t","commonfactor",
                       "groupfactor",
                       "personalfactor",
                       "b_seed","e_seed")
  param <- as.data.frame(param)
  
  return(list(y=y, x=x, a=a, b=b, e=e, tracker=param))
  
}

#' Update Continuous Ideology
#'
#' @description Update ideology of agents with the linear model.
#'
#' @param ideology The list of ideology parameters initiated by \code{\link{setideology}}.
#' @param commonfactor The movement in common factor.
#' @param groupfactor The movement in group influence on ideology.
#' @param personalfactor The SD for the movemebnt in personal influence on ideology.
#' @param b_seed The random number seed to draw group influence used when 
#' \code{groupfactor_type=="constant"} when the ideolog is initiated.
#' @param e_seed The random number seed to draw personal influence. 
#'
#' @return The updated list of ideology parameters 
#' \itemize{
#'   \item \code{y} The agents' ideologies. 
#'   \item \code{x} The group values rescaled in -1 to 1 range.
#'   \item \code{a} $a_{t-1} + $\code{commonfactor}
#'   \item \code{b} $b_{t-1} + $\code{groupfactor} or $b_{i,t-1} + abs(N(0,abs(groupfactor)))$
#'   \item \code{e} $e_{i, t-1} + N(0, personalfactor)$ 
#'   \item \code{tracker} The tracker table of ideology.
#' }
#' 
#' @details The ideology for agent $i$ $(y_i)$ is recalculated by the linear function: 
#' $y = a + b_i * x_i + e_i$.
#' 
#' @examples
#' 
#' # initiate ideology
#' set.seed(30)
#' x <- sample(rep(c(-1,1),each=1000))
#' ide <- setideology(x, 0, 0.8, 1)
#' 
#' # update
#' idex <- updateideology(ide, 0)
#' par(mfrow=c(1,3),mar=c(2,2,2,2))
#' hist(ide$y); mean(ide$y); sd(ide$y)
#' hist(idex$y); mean(idex$y); sd(idex$y)
#' plot(ide$y, idex$y);
#' 
#' @export

updateideology <- function(ideology, 
                           commonfactor=0, 
                           groupfactor=0, 
                           personalfactor=0.2,
                           b_seed=NULL, e_seed=NULL) {
  
  x <- ideology$x
  
  # Set and store random seed values (if not set)
  if (is.null(b_seed)) b_seed <- sample(-2^15:2^15, 1)
  if (is.null(e_seed)) e_seed <- sample(-2^15:2^15, 1)
  
  # Common Factor Update
  a <- commonfactor <- ideology$a + commonfactor
  # Group Factor Update
  if (is.na(ideology$tracker$b_seed[1])) {
    b_seed <- NA
    b <- groupfactor <- ideology$b + groupfactor
  } else {
    set.seed(b_seed)
    b <- ideology$b + abs(rnorm(length(x),sd=abs(groupfactor)))
  }
  # Personal Factor Update
  set.seed(e_seed)
  e <- ideology$e + rnorm(length(x),mean=0,sd=personalfactor)
  
  # Generate Ideology
  y <- a + b*x + e
  
  # Store Parameters
  param <- cbind(nrow(ideology$tracker),
                 commonfactor, groupfactor,personalfactor,
                 b_seed,e_seed)
  colnames(param) <- c("t","commonfactor",
                       "groupfactor",
                       "personalfactor",
                       "b_seed","e_seed")
  param <- as.data.frame(param)
  
  ideology$y <- y
  ideology$a <- a
  ideology$b <- b
  ideology$e <- e
  ideology$tracker <- rbind(ideology$tracker,param)
  
  return(ideology)
  
}

#' Set Knowledge
#'
#' @description Set knowledge levels of agents (Now only applicable with two groups).
#'
#' @param group The vector of group membership of agents.
#' @param mean1 The lowest average level of knowledge for a group.
#' @param mean2 The highest average level of knowledge for a group.
#' @param spread The spread (SD for normal distribution) of knowledg within a group.
#' @param kn_seed Random seed value to draw knowledge
#' 
#' @return The updated list of ideology parameters 
#' \itemize{
#'   \item \code{knp} Knowledge level converted to 0-1 range by inverse-standard-normal function. 
#'   \item \code{kn} The raw knowledge score
#'   \item \code{kn_seed} random seed value used.
#' }
#' 
#' @details If two groups, the raw knowledge score is generated by the normal distribution with 
#' mean \code{mean1} and standard deviation \code{spread} for the first group; 
#' the normal distribution with mean \code{mean1} and standard deviation \code{spread} 
#' for the second group. If more than two groups, additional group means are assigned by the 
#' points that equally divides the interval between \code{mean1} and \code{mean2}.
#' 
#' @examples
#' 
#' # initiate groups
#' set.seed(30)
#' x <- sample(rep(c(-1,1),each=1000))
#' 
#' # Set Knowledge
#' 
#' # sd of 0.85 produce approx. 0.25 SD after conversion.
#' kn <- setknowledge(x, 0.65, 0.65, 0.85)
#' par(mfrow=c(2,2),mar=c(2,2,2,2))
#' hist(kn$p); mean(kn$p); sd(kn$p)
#' hist(kn$y); mean(kn$y); sd(kn$y)
#' hist(kn$p[x==-1]); mean(kn$p[x==-1]); sd(kn$p[x==-1])
#' hist(kn$p[x==1]); mean(kn$p[x==1]); sd(kn$p[x==1])
#' 
#' # sd of 0.8 produce approx. 0.25 SD after conversion.
#' kn <- setknowledge(x, 0.60, 0.70, 0.85)
#' par(mfrow=c(2,2),mar=c(2,2,2,2))
#' hist(kn$p); mean(kn$p); sd(kn$p)
#' hist(kn$y); mean(kn$y); sd(kn$y)
#' hist(kn$p[x==-1]); mean(kn$p[x==-1]); sd(kn$p[x==-1])
#' hist(kn$p[x==1]); mean(kn$p[x==1]); sd(kn$p[x==1])
#' 
#' @export

setknowledge <- function(group, mean1=0.60, mean2=0.70, spread=0.8, kn_seed=NULL) {
  
  if (is.null(kn_seed)) kn_seed <- sample(-2^15:2^15, 1)
  
  group <- (((group - min(group))/(max(group)-min(group))) - 0.5)*2
  groupN <- length(unique(group))
  kn <- rep(NA, length(group))
  
  if (groupN==2) {
    set.seed(kn_seed)
    seedagain <- sample(-2^15:2^15,2)
    set.seed(seedagain[1])
    kn[group==-1] <- rnorm(length(group[group==-1]),qnorm(mean1),spread)
    set.seed(seedagain[2])
    kn[group==1] <- rnorm(length(group[group==1]),qnorm(mean2),spread)
  } else if (groupN>=2) {
    set.seed(kn_seed)
    seedagain <- sample(-2^15:2^15,groupN)
    means <- seq(mean1, mean2, length=groupN)
    for (i in 1:groupN) {
      kn[group==sort(unique(group))[i]] <- 
        rnorm(length(group[group==sort(unique(group))[i]]),qnorm(means[i]),spread)
    }  
  } else {
    stop("Number of unique values in group must be 2 or more.")
  }

  knp <- pnorm(kn)
  
  return(list(p=knp, y=kn, seed=kn_seed))
  
}

#' Set Characteristics of Two Parties
#'
#' @description Set ideologies and capacity of parties (i.e., blue and red) running in election.
#'
#' @param blue_startloc Starting location of blue party ideology. 
#' @param red_startloc Starting location of red party ideology.
#' @param startloc_noise Random noise given to the starting locations of party ideologies.
#' @param blue_startcapacity Starting level of blue party capacity.
#' @param red_startcapacity Starting level of red party capacity.
#' @param startcapacity_noise Random noise given to the starting level of party capacities.
#' @param startincumbent Incumbent in the first election. If \code{NULL}, randomly chosen.
#' @param blue_loc_seed Random number seed to draw blue party ideology.
#' @param red_loc_seed Random number seed to draw red party ideology.
#' @param blue_capacity_seed Random number seed to draw blue party capacity.
#' @param red_capacity_seed Random number seed to draw red party capacity.
#' 
#' @return The updated list of ideology parameters 
#' \itemize{
#'   \item \code{loc} ideology locations of two parties. 
#'   \item \code{capacity} capacity levels of two parties.
#'   \item \code{tracker} Tracker of the history of ideogoloy and capacity.
#' }
#' 
#' @export

setparty <- function(blue_startloc=-0.8,
                     red_startloc=0.8,
                     startloc_noise=0.1,
                     blue_startcapacity=0,
                     red_startcapacity=0,
                     startcapacity_noise=0.1,
                     startincumbent = NULL,
                     blue_loc_seed=NULL,
                     red_loc_seed=NULL,
                     blue_capacity_seed=NULL,
                     red_capacity_seed=NULL) {
  
  if (is.null(blue_loc_seed)) blue_loc_seed <- sample(-2^15:2^15, 1)
  if (is.null(red_loc_seed)) red_loc_seed <- sample(-2^15:2^15, 1)
  if (is.null(blue_capacity_seed)) blue_capacity_seed <- sample(-2^15:2^15, 1)
  if (is.null(red_capacity_seed)) red_capacity_seed <- sample(-2^15:2^15, 1)
  if (is.null(startincumbent)) startincumbent <- ifelse(rbinom(1,1,0.5)==1,"red","blue")
  
  # Ideological Location of Parties
  set.seed(blue_loc_seed)
  blue_loc <- rnorm(1,mean=blue_startloc,sd=startloc_noise)
  set.seed(red_loc_seed)
  red_loc <- rnorm(1,mean=red_startloc,sd=startloc_noise)
  
  # capacity Level of Parties (Randomly Set only if Incumbent)
  if (startincumbent=="blue") {
    set.seed(blue_capacity_seed)
    blue_capacity <- rnorm(1,mean=blue_startcapacity,sd=startcapacity_noise)
    red_capacity_seed <- NA
    red_capacity <- red_startcapacity
  } else if (startincumbent=="red") {
    blue_capacity_seed <- NA
    blue_capacity <- blue_startcapacity
    set.seed(red_capacity_seed)
    red_capacity <- rnorm(1,mean=red_startcapacity,sd=startcapacity_noise)
  } else {
    stop("Invalid startincumbent! Must be 'blue' or 'red'.")
  }
  
  trackdt <- as.data.frame(cbind(1,NA,
                                 blue_loc,red_loc,blue_capacity,red_capacity,
                                 blue_startloc,red_startloc,startloc_noise,
                                 blue_startcapacity,red_startcapacity,startcapacity_noise,
                                 blue_loc_seed,red_loc_seed,
                                 blue_capacity_seed,red_capacity_seed))
  colnames(trackdt) <- c("t","incumbent",
                         "blue_loc","red_loc",
                         "blue_capacity","red_capacity",
                         "blue_loc_direction","red_loc_direction","loc_noise",
                         "blue_capacity_direction","red_capacity_direction","capacity_noise",
                         "blue_loc_seed","red_loc_seed",
                         "blue_capacity_seed","red_capacity_seed")
  trackdt$incumbent <- factor(startincumbent, levels=c("red","blue"))
  
  return(list(loc=list(blue=blue_loc, red=red_loc), 
              capacity=list(blue=blue_capacity,red=red_capacity),
              tracker=trackdt))

}

#' Update Characteristics of Two Parties
#'
#' @description Set ideologies and capacity of parties (i.e., blue and red) running in election.
#' 
#' @param party The lists generated initially by \code{\link{setparty}}.
#' @param incumbent Incumbent party in the current election. Either \code{"blue"} or \code{"red"}.
#' @param loc_noise Random noise in updating locations of party ideologies.
#' @param capacity_noise Random noise in updating party capacity.
#' @param blue_loc_direction Direction and expected quantity of the movement in blue party ideology. 
#' @param red_loc_direction Direction and expected quantity of the movement in red party ideology.
#' @param blue_capacity_direction Direction and expected quantity of the movement in blue party capacity.
#' @param red_capacity_direction Direction and expected quantity of the movement in red party capacity.
#' @param update_type The type of updating procedure. 
#' If \code{"oneshot"}, the update is always given over \code{startloc} and \code{startcapacity} defined in 
#' \code{\link{setparty}} at each election. If \code{"online"}, the update is given over the ideology and 
#' capacity realized in the previous election. If \code{"fromstart"} the update is given over the 
#' ideology and capacity realized in the first election of the simulation run.  
#' @param blue_loc_seed Random number seed to draw blue party ideology.
#' @param red_loc_seed Random number seed to draw red party ideology.
#' @param blue_capacity_seed Random number seed to draw blue party capacity.
#' @param red_capacity_seed Random number seed to draw red party capacity.
#' 
#' @return The updated list of ideology parameters: 
#' \itemize{
#'   \item \code{loc} ideology locations of two parties. 
#'   \item \code{capacity} capacity levels of two parties.
#'   \item \code{tracker} Tracker of the history of ideogoloy and capacity.
#' }
#' 
#' @examples
#' 
#' require(ggplot2)
#' 
#' par <- setparty()
#' for (t in 1:14) {
#'   incumbent <- sample(c("red","blue"),1)
#'   par <- updateparty(par,incumbent, 0, 0.7)
#' }
#' # party
#' ggplot(par$tracker,aes(x=t)) +
#'   geom_line(aes(y=blue_loc,color="2.blue")) +
#'   geom_line(aes(y=red_loc,color="1.red")) +
#'   theme_bw()
#' # capacity
#' ggplot(par$tracker,aes(x=t)) +
#'   geom_line(aes(y=blue_capacity,color="2.blue")) +
#'   geom_line(aes(y=red_capacity,color="1.red")) +
#'   theme_bw()
#' 
#' @export

updateparty <- function(party, incumbent,
                        loc_noise = 0, capacity_noise=0.7,
                        blue_loc_direction=0, red_loc_direction=0, 
                        blue_capacity_direction=0, red_capacity_direction=0,
                        update_type = "oneshot",
                        blue_loc_seed=NULL,red_loc_seed=NULL,
                        blue_capacity_seed=NULL,red_capacity_seed=NULL) {
  
  if (is.null(blue_loc_seed)) blue_loc_seed <- sample(-2^15:2^15, 1)
  if (is.null(red_loc_seed)) red_loc_seed <- sample(-2^15:2^15, 1)
  if (is.null(blue_capacity_seed)) blue_capacity_seed <- sample(-2^15:2^15, 1)
  if (is.null(red_capacity_seed)) red_capacity_seed <- sample(-2^15:2^15, 1)
  
  if (update_type=="online") {
    blue_loc_base <- party$loc$blue
    red_loc_base <- party$loc$red
    blue_capacity <- blue_capacity_base <- party$capacity$blue
    red_capacity <- red_capacity_base <- party$capacity$red
  } else if (update_type=="oneshot") {
    blue_loc_base <- party$tracker$blue_loc_direction[1]
    red_loc_base <- party$tracker$red_loc_direction[1]
    blue_capacity <- blue_capacity_base <- party$tracker$blue_capacity_direction[1]
    red_capacity <- red_capacity_base <- party$tracker$red_capacity_direction[1]
  } else if (update_type=="fromstart") {
    lue_loc_base <- party$tracker$blue_loc[1]
    red_loc_base <- party$tracker$red_loc[1]
    blue_capacity <- blue_capacity_base <- party$tracker$blue_capacity[1]
    red_capacity <- red_capacity_base <- party$tracker$red_capacity[1]
  } else {
    stop("invalid update_type! It must be 'online', 'oneshot', or 'fromstart'.")
  }
  
  # Policy Location
  if (loc_noise>0) {
    set.seed(blue_loc_seed)
    blue_loc_plus <- rnorm(1,blue_loc_direction,loc_noise)
    set.seed(red_loc_seed)
    red_loc_plus <- rnorm(1,red_loc_direction,loc_noise)
  } else {
    blue_loc_seed <- red_loc_seed <- NA
    blue_loc_plus <- blue_loc_direction
    red_loc_plus <- red_loc_direction
  }
  
  blue_loc <- blue_loc_base + blue_loc_plus
  red_loc <- red_loc_base + red_loc_plus

  # capacity Level (Movement only when incumbent)
  if (incumbent=="blue") {
    red_capacity_seed <- NA
    set.seed(blue_capacity_seed)
    blue_capacity_plus <- rnorm(1,blue_capacity_direction,capacity_noise)
    blue_capacity <- blue_capacity_base + blue_capacity_plus
  } else if (incumbent=="red") {
    blue_capacity_seed <- NA
    set.seed(red_capacity_seed)
    red_capacity_plus <- rnorm(1,red_capacity_direction,capacity_noise)
    red_capacity <- red_capacity_base + red_capacity_plus
  } else {
    stop("invalid incumbent value!")
  }
  
  # Tracking Data
  trackline <- as.data.frame(cbind(t=nrow(party$track)+1, incumbent=NA,
                                   blue_loc, red_loc, 
                                   blue_capacity, red_capacity,
                                   blue_loc_direction,red_loc_direction,loc_noise,
                                   blue_capacity_direction,red_capacity_direction,capacity_noise,
                                   blue_loc_seed,red_loc_seed,
                                   blue_capacity_seed,red_capacity_seed))
  trackline$incumbent <- incumbent
  
  trackdt <- rbind(party$track, trackline)
  
  return(list(loc=list(blue=blue_loc, red=red_loc), 
              capacity=list(blue=blue_capacity, red=red_capacity), 
              tracker=trackdt))

}

#' Get Preference Signals
#'
#' @description Obtain preference signals for voters.
#' 
#' @param party The lists generated initially by \code{\link{setparty}} and may be updated by \code{\link{updateparty}}.
#' @param ideology The lists generated initially by \code{\link{setideology}} and may be updated by \code{\link{updateideology}}.
#' @param knowledge The lists generated by \code{\link{setknowledge}}.
#' @param difficulty The extent of noise level in signal.
#' @param bias_common The extent of common directional bias (i.e., positive value advantages red party and vice versa) in signal. 
#' @param bias_ingroup The extent of common ingroup bias (i.e., Evaluate the ingroup candidate better).
#' @param pastbias Currently not used.  
#' @param idedist_blue_signal_seed The random number seed to draw the signal regarding ideological distance from blue party.
#' @param idedist_red_signal_seed The random number seed to draw the signal regarding ideological distance from red party.
#' @param capacity_signal_seed The random number seed to draw the signal regarding the capacity of incumbent party.
#' 
#' @return The list of signal parameters: 
#' \itemize{
#'   \item \code{signal} The list of signals regarding ideological distance from parties and incumbent capacity.
#'   \item \code{truth} The list of truth regarding ideological distance from parties and incumbent capacity.
#'   \item \code{seed} Random number seeds used.
#' }
#' 
#' @examples
#' 
#' # initiate ideology, knowledge and party
#' set.seed(30)
#' x <- sample(rep(c(-1,1),each=1000))
#' ide <- setideology(x, 0, 0.8, 1)
#' kn <- setknowledge(x, 0.60, 0.70, 0.85)
#' par <- setparty()
#' 
#' # Get signals
#' sig <- getsignal(par,ide,kn,3) # Difficulty of 3.3 creates approx. 1.9 SD
#' sd(sig$signal$idedist_blue) # around 1.9
#' sd(sig$signal$idedist_red) # around 1.9
#' par(mfrow=c(2,2),mar=c(2,2,2,2))
#' plot(sig$truth$idedist_blue,sig$signal$idedist_blue)
#' plot(sig$truth$idedist_red,sig$signal$idedist_red)
#' plot(sig$truth$idedist_red[kn$p>0.6],
#'      sig$signal$idedist_red[kn$p>0.6])
#' plot(sig$truth$idedist_red[kn$p<0.6],
#'      sig$signal$idedist_red[kn$p<0.6])
#' 
#' @export

getsignal <- function(party, 
                      ideology, 
                      knowledge, 
                      difficulty = 3, 
                      bias_common = 0, 
                      bias_ingroup = 0, 
                      pastbias=NULL,
                      idedist_blue_signal_seed=NULL,
                      idedist_red_signal_seed=NULL,
                      capacity_signal_seed=NULL) {
  
  if (is.null(idedist_blue_signal_seed)) idedist_blue_signal_seed <- sample(-2^15:2^15, 1)
  if (is.null(idedist_red_signal_seed)) idedist_red_signal_seed <- sample(-2^15:2^15, 1)
  if (is.null(capacity_signal_seed)) capacity_signal_seed <- sample(-2^15:2^15, 1)
  
  # Ingroup Bias
  bias_ingroup_blue <- bias_ingroup*(ideology$x==-1)
  bias_ingroup_red <- bias_ingroup*(ideology$x==1)
  
  # True Ideology
  idedist_blue <- ideology$y - party$loc$blue 
  idedist_red <- ideology$y - party$loc$red
  # Blue Ideology Signal
  set.seed(idedist_blue_signal_seed)
  idedist_blue_signal <- idedist_blue + 
    sapply(knowledge$p, function(k) rnorm(1, bias_common*(1-k), sd=(1-k)*difficulty))
  idesist_blue_signal <- idedist_blue_signal*(1-bias_ingroup_blue*(1-knowledge$p))
  # Red Ideology Signal  
  set.seed(idedist_red_signal_seed)
  idedist_red_signal <- idedist_red + 
    sapply(knowledge$p, function(k) rnorm(1, bias_common*(1-k), sd=(1-k)*difficulty))
  idesist_red_signal <- idedist_red_signal*(1-bias_ingroup_red*(1-knowledge$p))
  
  # True Capacity
  capacity_redadv <- party$capacity$red - party$capacity$blue
  # Capacity Signal
  set.seed(capacity_signal_seed)
  capacity_signal <- capacity_redadv + 
    sapply(knowledge$p, function(k) rnorm(1, bias_common*(1-k), sd=(1-k)*difficulty))
    capacity_signal[capacity_signal<0] <- capacity_signal[capacity_signal<0] * 
    (1-bias_ingroup_red[capacity_signal<0]*(1-knowledge$p[capacity_signal<0]))
  capacity_signal[capacity_signal>0] <- capacity_signal[capacity_signal>0] * 
    (1-bias_ingroup_blue[capacity_signal>0]*(1-knowledge$p[capacity_signal>0]))
  
  return(list(signal = list(idedist_blue=idedist_blue_signal,
                            idedist_red=idedist_red_signal,
                            capacity_redadv=capacity_signal,
                            bias_common=bias_common,
                            bias_ingroup=bias_ingroup),
              truth = list(idedist_blue=idedist_blue,
                           idedist_red=idedist_red,
                           capacity_redadv=capacity_redadv),
              seed = list(idedist_blue=idedist_blue_signal_seed,
                          idedist_red=idedist_red_signal_seed,
                          capacity_redadv=capacity_signal_seed)))
  
} 

#' Get Neighborhood Context
#'
#' @description Get neighborhood context from society data
#' 
#' @param soc The lists generated initially by \code{\link{initialsoc}}. 
#' Values of previous elections votes (i.e., \code{"grid_votes1"} and \code{"grid_votes2"}) must be included.
#' @param personloc The coordinates with agents.
#' @param radius Search radius to collect neighborhood preferences.
#' @param timetoagg How many previous elections to average. Either 1 or 2 (default). 
#' 
#' @return The list with:
#' \itemize{
#'   \item \code{local} Local vote share in previous elections for voters 
#'   (i.e., the voting patterns of neightbors in search radius)
#'   \item \code{national} National vote share in previous elections. 
#'   \item \code{detail} Detailed results.
#' }
#' 
#' @export

getcontext <- function(soc, personloc, radius, timetoagg=2) {
  
  if (!timetoagg%in%c(1,2)) stop("Invalid timetoagg value. It must be 1 or 2")
  
  # Collect Votes in Local Context
  collect_votes <- function(coords, election) {
    nbrs <- get_neighbors(coords=coords, radius=radius, spacelim=dim(soc$grid))
    nbrs_values <- apply(nbrs, 1, function(k) soc[[election]][k[1],k[2]])
    list(me=soc[[election]][coords[1],coords[2]], nbrs=nbrs_values)
  }
  cxtvotes <- apply(soc$coords[personloc,], 1, collect_votes, election="grid_votes1")
  if (timetoagg==2) {
    cxtvotes2 <- apply(soc$coords[personloc,], 1, collect_votes, election="grid_votes2")
  }

  # Calcualte Vote Share (of paty red) in Local Context
  localshare <- function(cxt) {
    if (all(is.na(cxt$nbrs))) {
      voteshare <- 0.5
    } else {
      voteshare <- sum(cxt$nbrs, na.rm=TRUE)/length(cxt$nbrs[!is.na(cxt$nbrs)]) 
    }
    return(voteshare)
  }
  # Store in Society Object
  votes_localshare <- sapply(cxtvotes, localshare) 
  if (timetoagg==2) {
    votes_localshare2 <- sapply(cxtvotes2, localshare)
  }
  
  # National Vote Share (of party red)
  votes_nationalshare <- 
    sum(as.numeric(soc$grid_votes1[personloc]))/length(personloc)
  if (timetoagg==2) {
    votes_nationalshare2 <- sum(as.numeric(soc$grid_votes2[personloc]))/length(personloc)
  }
  
  # Return
  if (timetoagg==2) {
    votes_localshare <- (votes_localshare + votes_localshare2)/2
    votes_nationalshare <- (votes_nationalshare + votes_nationalshare2)/2
    cxtvotes <- list(cxtvotes, cxtvotes2)
  }
  return(list(local = votes_localshare,
              national = votes_nationalshare,
              detail = cxtvotes))
}

#' Get Regional Context for Voting
#'
#' @description Get regional context from society data
#' 
#' @param soc The lists generated initially by \code{\link{initialsoc}}. 
#' Values of previous elections votes (i.e., \code{"grid_votes1"} and \code{"grid_votes2"}) must be included.
#' @param personloc The coordinates with agents.
#' @param timetoagg How many previous elections to average. Either 1 or 2 (default). 
#' 
#' @return The list with:
#' \itemize{
#'   \item \code{sublocal} Sub-regional vote share in previous elections (Standard contextual voting Used).
#'   \item \code{local} Regional vote share in previous elections (Standard contextual voting used).
#'   \item \code{national} National vote share in previous elections (Standard contextual voting used). 
#'   \item \code{sublocal1} Sub-regional vote share in previous elections (Local only contextual voting Used).
#'   \item \code{local1} Regional vote share in previous elections (Local only contextual voting used).
#'   \item \code{national1} National vote share in previous elections (Local only contextual voting used). 
#'   \item \code{sublocal2} Sub-regional vote share in previous elections (National only contextual voting Used).
#'   \item \code{local2} Regional vote share in previous elections (National only contextual voting used).
#'   \item \code{national2} National vote share in previous elections (National only contextual voting used). 
#' }
#' 
#' @export

getregion <- function(soc, personloc, timetoagg=2) {
  
  if (!timetoagg%in%c(1,2)) stop("Invalid timetoagg value. It must be 1 or 2")
  
  sr <- as.numeric(soc$grid_subregion)[personloc]
  r <- as.numeric(soc$grid_region)[personloc]
  v <- as.numeric(soc$grid_votes1)[personloc]
  if (!is.null(soc[["votes"]])) {
    v1x <- as.numeric(soc$votes$vote$actual1[,ncol(soc$votes$vote$actual1)])
    v1y <- as.numeric(soc$votes$vote$actual2[,ncol(soc$votes$vote$actual2)])
  } else {
    v1x <- v1y <- v
  }

  # Collect Region Stats
  regioncharacter <- function(k,v,r) {
    total <- length(r[r==k])
    voteshare <- mean(v[r==k]) # red party share
    res <- as.numeric(c(k, voteshare, total))
    return(res)
  }
  rs <- as.data.frame(t(sapply(seq(1,max(r),1), regioncharacter, v=v, r=r)))
  rsx <- as.data.frame(t(sapply(seq(1,max(r),1), regioncharacter, v=v1x, r=r)))
  rsy <- as.data.frame(t(sapply(seq(1,max(r),1), regioncharacter, v=v1y, r=r)))
  colnames(rs) <- colnames(rsx) <- colnames(rsy) <- 
    c("region","voteshare","totalvotes")
  srs <- as.data.frame(t(sapply(seq(1,max(sr),1), regioncharacter, v=v, r=sr)))
  srsx <- as.data.frame(t(sapply(seq(1,max(sr),1), regioncharacter, v=v1x, r=sr)))
  srsy <- as.data.frame(t(sapply(seq(1,max(sr),1), regioncharacter, v=v1y, r=sr)))
  colnames(srs) <- colnames(srsx) <- colnames(srsy) <- 
    c("subregion","voteshare","totalvotes")
  
  # Match SubRegion and Vote Share for Individuals
  votes_sublocalshare <- srs$voteshare[match(sr, srs$subregion)]
  votes_sublocalsharex <- srs$voteshare[match(sr, srsx$subregion)]
  votes_sublocalsharey <- srs$voteshare[match(sr, srsy$subregion)]
  
  # Match Region and Vote Share for Individuals
  votes_localshare <- rs$voteshare[match(r, rs$region)]
  votes_localsharex <- rs$voteshare[match(r, rsx$region)]
  votes_localsharey <- rs$voteshare[match(r, rsy$region)]
  
  # National Vote Share (of party red)
  votes_nationalshare <- sum(v)/length(personloc)
  votes_nationalsharex <- sum(v1x)/length(personloc)
  votes_nationalsharey <- sum(v1y)/length(personloc)
  
  # Aggregate and averge t-2
  if (timetoagg==2) {
    v2 <- as.numeric(soc$grid_votes2)[personloc]
    if (!is.null(soc[["votes"]])) {
      if (ncol(soc$votes$vote$actual)>1) {
        v2x <- as.numeric(soc$votes$vote$actual1[,(ncol(soc$votes$vote$actual1)-1)])
        v2y <- as.numeric(soc$votes$vote$actual2[,(ncol(soc$votes$vote$actual2)-1)])
      } else {
        v2x <- v2y <- v2
      }
    } else {
      v2x <- v2y <- v2
    }
    rs2 <- as.data.frame(t(sapply(seq(1,max(r),1), regioncharacter, v=v2, r=r)))
    rs2x <- as.data.frame(t(sapply(seq(1,max(r),1), regioncharacter, v=v2x, r=r)))
    rs2y <- as.data.frame(t(sapply(seq(1,max(r),1), regioncharacter, v=v2y, r=r)))
    colnames(rs2) <- colnames(rs2x) <- colnames(rs2y) <- 
      c("region","voteshare","totalvotes")
    srs2 <- as.data.frame(t(sapply(seq(1,max(sr),1), regioncharacter, v=v2, r=sr)))
    srs2x <- as.data.frame(t(sapply(seq(1,max(sr),1), regioncharacter, v=v2x, r=sr)))
    srs2y <- as.data.frame(t(sapply(seq(1,max(sr),1), regioncharacter, v=v2y, r=sr)))
    colnames(srs2) <- colnames(srs2x) <- colnames(srs2y) <- 
      c("subregion","voteshare","totalvotes")
    votes_sublocalshare2 <- srs2$voteshare[match(sr, srs2$subregion)]
    votes_sublocalshare2x <- srs2$voteshare[match(sr, srs2x$subregion)]
    votes_sublocalshare2y <- srs2$voteshare[match(sr, srs2y$subregion)]
    votes_localshare2 <- rs2$voteshare[match(r, rs2$region)]
    votes_localshare2x <- rs2$voteshare[match(r, rs2x$region)]
    votes_localshare2y <- rs2$voteshare[match(r, rs2y$region)]
    votes_nationalshare2 <- sum(v2)/length(personloc)
    votes_nationalshare2x <- sum(v2x)/length(personloc)
    votes_nationalshare2y <- sum(v2y)/length(personloc)
    votes_sublocalshare <- (votes_sublocalshare + votes_sublocalshare2)/2
    votes_sublocalsharex <- (votes_sublocalsharex + votes_sublocalshare2x)/2
    votes_sublocalsharey <- (votes_sublocalsharey + votes_sublocalshare2y)/2
    votes_localshare <- (votes_localshare + votes_localshare2)/2
    votes_localsharex <- (votes_localsharex + votes_localshare2x)/2
    votes_localsharey <- (votes_localsharey + votes_localshare2y)/2
    votes_nationalshare <- (votes_nationalshare + votes_nationalshare2)/2
    votes_nationalsharex <- (votes_nationalsharex + votes_nationalshare2x)/2
    votes_nationalsharey <- (votes_nationalsharey + votes_nationalshare2y)/2
  }
  
  # Return
  return(list(sublocal = votes_sublocalshare,
              local = votes_localshare,
              national = votes_nationalshare,
              sublocal1 = votes_sublocalsharex,
              sublocal2 = votes_sublocalsharey,
              local1 = votes_localsharex,
              local2 = votes_localsharey,
              national1 = votes_nationalsharex,
              national2 = votes_nationalsharey))
}

#' Vote in Election
#'
#' @description Voters make voting decisions.
#' 
#' @param soc The lists generated initially by \code{\link{initialsoc}}. 
#' Also following elements must be contained within \code{soc} object:
#' \itemize{
#'  \item \code{knowledge} Created by \code{\link{setknowledge}}
#'  \item \code{ideology} Created by \code{\link{setideology}}, may be updated by \code{\link{updateideology}}
#'  \item \code{context} Created by \code{\link{getregion}}
#'  \item \code{signal} Created by \code{\link{getsignal}}
#'  \item \code{incumbent} Character of incumbent party. Either \code{"blue"} or \code{"red"}.
#' }
#' @param b_idedist The influence of ideological distance. The default = 0.198, # From ANES
#' @param b_capacity The influence of capacity The default = 0.716, # From ANES
#' @param b_incumbent The influence if incumebent. The default = 1.526, # From ANES (EDV result)
#' @param b_pvi_nation The influence of national PVI. The default = -9.827, # From ANES (EDV result)
#' @param b_pvi_local The influence of regional PVI. The default = 5.156, # From CCES pvi_state
#' @param b_pvi_sublocal The influence of sub-regional PVI. The default = 1.663, # From CCES pvi_state
#' @param b_pid The influence of group membership. The default = 1.215, # From ANES
#' @param vote_seed The random number seed.
#' 
#' @return The list with:
#' \itemize{
#'   \item \code{vote} The list of vote decisions by different voting functions.
#'   \item \code{prob} The list of vote probabilities by different voting functions.
#'   \item \code{pref} The list of raw vote preferences by different voting functions. 
#'   \item \code{noidprob} The list of vote probabilities (if group membership has no effect).
#'   \item \code{b} The list of parameters.
#'   \item \code{tracker} Tracker of voting history.
#'   \item \code{vote_seed} The random number seed used. 
#' }
#' 
#' @importFrom faraway ilogit
#' 
#' @export

letsvote <- function(soc,
                     b_idedist = 0.198, # From ANES
                     b_capacity = 0.716, # From ANES
                     b_incumbent = 1.526, # From ANES (EDV result)
                     b_pvi_nation = -9.827, # From ANES (EDV result)
                     b_pvi_local = 5.156, # From CCES pvi_state
                     b_pvi_sublocal = 1.663, # From CCES pvi_state
                     b_pid = 1.215, # From ANES
                     vote_seed=NULL
                     ) {
  
  knowledge <- soc$knowledge
  signal <- soc$signal
  context <- soc$context
  ideology <- soc$ideology
  incumbent <- soc$incumbent
  
  if(is.null(vote_seed)) vote_seed <- sample(-2^15:2^15,1)
  
  # Context Based Voting (Probability)
  if (incumbent == "blue") b_incumbent <- -b_incumbent
  cxtpref <- b_incumbent * 0.5  + b_pvi_nation * ((context$national) - 0.5) +
    b_pvi_local * (context$local-context$national) + 
    b_pvi_sublocal * (context$sublocal-context$local)
    
  # no balancing
  cxtpref1 <- b_incumbent + 
    b_pvi_local * (context$local1 - context$national1) + 
    b_pvi_sublocal * (context$sublocal1 - context$local1)
  # no local context
  cxtpref2 <- b_incumbent + b_pvi_nation * (context$national2 - 0.5)
  
  # Identity Based Voting
  pidpref <- b_pid * ideology$x
  
  # Preference Based Voting 
  # Signaled Preference
  signalpref <- b_idedist*(signal$signal$idedist_blue^2 - signal$signal$idedist_red^2) + 
    b_capacity*signal$signal$capacity_redadv
  # Signaled Preference
  wsignalpref <- knowledge$p * signalpref
  # Informed preference
  informedpref <- b_idedist*(signal$truth$idedist_blue^2 - signal$truth$idedist_red^2) + 
    b_capacity*signal$truth$capacity_redadv
  # Actual Voting Preference
  actualpref <- wsignalpref + (1-knowledge$p) * cxtpref + pidpref
  actualpref1 <- wsignalpref + (1-knowledge$p) * cxtpref1 + pidpref
  actualpref2 <- wsignalpref + (1-knowledge$p) * cxtpref2 + pidpref
  
  # Probabilities of Voting
  cxtprob <- ilogit(cxtpref + pidpref)
  cxtprob1 <- ilogit(cxtpref1 + pidpref)
  cxtprob2 <- ilogit(cxtpref2 + pidpref)
  signalprob <- ilogit(signalpref + pidpref)
  wsignalprob <- ilogit(wsignalpref + pidpref)
  informedprob <- ilogit(informedpref + pidpref)
  noidcxtprob <- ilogit(cxtpref)
  noidcxtprob1 <- ilogit(cxtpref1)
  noidcxtprob2 <- ilogit(cxtpref2)
  noidsignalprob <- ilogit(signalpref)
  noidwsignalprob <- ilogit(wsignalpref)
  noidinformedprob <- ilogit(informedpref)
  
  # Actual Voting Probability
  actualprob <- ilogit(actualpref)
  actualprob1 <- ilogit(actualpref1)
  actualprob2 <- ilogit(actualpref2)
  
  # Voting
  set.seed(vote_seed)
  cxtvote <- sapply(cxtprob, function(k) rbinom(1,1,prob=k))
  set.seed(vote_seed)
  cxtvote1 <- sapply(cxtprob1, function(k) rbinom(1,1,prob=k))
  set.seed(vote_seed)
  cxtvote2 <- sapply(cxtprob2, function(k) rbinom(1,1,prob=k))
  set.seed(vote_seed)
  signalvote <- sapply(signalprob, function(k) rbinom(1,1,prob=k))
  set.seed(vote_seed)
  wsignalvote <- sapply(wsignalprob, function(k) rbinom(1,1,prob=k))
  set.seed(vote_seed)
  actualvote <- sapply(actualprob, function(k) rbinom(1,1,prob=k))
  set.seed(vote_seed)
  actualvote1 <- sapply(actualprob1, function(k) rbinom(1,1,prob=k))
  set.seed(vote_seed)
  actualvote2 <- sapply(actualprob2, function(k) rbinom(1,1,prob=k))
  set.seed(vote_seed)
  informedvote <- sapply(informedprob, function(k) rbinom(1,1,prob=k))
  set.seed(vote_seed)
  noidinformedvote <- sapply(noidinformedprob, function(k) rbinom(1,1,prob=k))
  
  
  # Return Results
  if (is.null(soc[["votes"]])) {
    
    track_vote <- cbind(1,mean(actualvote),mean(actualvote1),mean(actualvote2),
                        mean(signalvote),mean(wsignalvote),
                        mean(cxtvote),mean(cxtvote1),mean(cxtvote2),
                        mean(informedvote),mean(noidinformedvote))
    colnames(track_vote) <- c("t","actual","actual1","actual2",
                              "signal","wsignal","context",
                              "context1","context2","informed","noidinformed")
    track_vote <- as.data.frame(track_vote)

    return(list(vote=list(actual=data.frame(t1=actualvote),
                          actual1=data.frame(t1=actualvote1),
                          actual2=data.frame(t1=actualvote2),
                          signal=data.frame(t1=signalvote),
                          wsignal=data.frame(t1=wsignalvote),
                          context=data.frame(t1=cxtvote),
                          context1=data.frame(t1=cxtvote1),
                          context2=data.frame(t1=cxtvote2),
                          informed=data.frame(t1=informedvote),
                          noidinformed=data.frame(t1=noidinformedvote)), 
                prob=list(actual=data.frame(t1=actualprob),
                          actual1=data.frame(t1=actualprob1),
                          actual2=data.frame(t1=actualprob2),
                          signal=data.frame(t1=signalprob),
                          wsignal=data.frame(t1=wsignalprob),
                          context=data.frame(t1=cxtprob),
                          context1=data.frame(t1=cxtprob1),
                          context2=data.frame(t1=cxtprob2),
                          informed=data.frame(t1=informedprob)), 
                pref=list(actual=data.frame(t1=actualpref),
                          actual1=data.frame(t1=actualpref1),
                          actual2=data.frame(t1=actualpref2),
                          signal=data.frame(t1=signalpref),
                          wsignal=data.frame(t1=wsignalpref),
                          context=data.frame(t1=cxtpref),
                          context1=data.frame(t1=cxtpref1),
                          context2=data.frame(t1=cxtpref2),
                          informed=data.frame(t1=informedpref),
                          pid=data.frame(t1=pidpref)),
                noidprob=list(signal=data.frame(t1=noidsignalprob),
                              wsignal=data.frame(t1=noidwsignalprob),
                              context=data.frame(t1=noidcxtprob),
                              context1=data.frame(t1=noidcxtprob1),
                              context2=data.frame(t1=noidcxtprob2),
                              informed=data.frame(t1=noidinformedprob)),
                b = list(idedist=b_idedist,
                         capacity=b_capacity,
                         incumbent=b_incumbent,
                         pvi_nation=b_pvi_nation,
                         pvi_local=b_pvi_local,
                         pid=b_pid),
                tracker = track_vote,
                vote_seed=vote_seed))
  } else {
    
    telec <- ncol(soc$votes$vote$actual) + 1
    
    setdf <- function(pastname,newname) {
      mydf <- eval(parse(text = paste0(
        "data.frame(", pastname,", t", telec, "=", newname, ")")))
      mydf
    }
    
    track_vote <- cbind(telec,mean(actualvote),mean(actualvote1),mean(actualvote2),
                        mean(signalvote),mean(wsignalvote),
                        mean(cxtvote),mean(cxtvote1),mean(cxtvote2),
                        mean(informedvote),mean(noidinformedvote))
    colnames(track_vote) <- c("t","actual","actual1","actual2",
                              "signal","wsignal","context",
                              "context1","context2","informed","noidinformed")
    track_vote <- as.data.frame(track_vote)
    track_vote <- rbind(soc$votes$tracker,track_vote)

    return(list(vote=list(actual=setdf("soc$votes$vote$actual","actualvote"),
                          actual1=setdf("soc$votes$vote$actual1","actualvote1"),
                          actual2=setdf("soc$votes$vote$actual2","actualvote2"),
                          signal=setdf("soc$votes$vote$signal","signalvote"),
                          wsignal=setdf("soc$votes$vote$wsignal","wsignalvote"),
                          context=setdf("soc$votes$vote$context","cxtvote"),
                          context1=setdf("soc$votes$vote$context1","cxtvote1"),
                          context2=setdf("soc$votes$vote$context2","cxtvote2"),
                          informed=setdf("soc$votes$vote$informed","informedvote"),
                          noidinformed=setdf("soc$votes$vote$noidinformed","noidinformedvote")), 
                prob=list(actual=setdf("soc$votes$prob$actual","actualprob"),
                          actual1=setdf("soc$votes$prob$actual1","actualprob1"),
                          actual2=setdf("soc$votes$prob$actual2","actualprob2"),
                          signal=setdf("soc$votes$prob$signal","signalprob"),
                          wsignal=setdf("soc$votes$prob$wsignal","wsignalprob"),
                          context=setdf("soc$votes$prob$context","cxtprob"),
                          context1=setdf("soc$votes$prob$context1","cxtprob1"),
                          context2=setdf("soc$votes$prob$context2","cxtprob2"),
                          informed=setdf("soc$votes$prob$informed","informedprob")), 
                pref=list(actual=setdf("soc$votes$pref$actual","actualpref"),
                          actual1=setdf("soc$votes$pref$actual1","actualpref1"),
                          actual2=setdf("soc$votes$pref$actual2","actualpref2"),
                          signal=setdf("soc$votes$pref$signal","signalpref"),
                          wsignal=setdf("soc$votes$pref$wsignal","wsignalpref"),
                          context=setdf("soc$votes$pref$context","cxtpref"),
                          context1=setdf("soc$votes$pref$context1","cxtpref1"),
                          context2=setdf("soc$votes$pref$context2","cxtpref2"),
                          informed=setdf("soc$votes$pref$informed","informedpref"),
                          pid=setdf("soc$votes$pref$pid","pidpref")), 
                noidprob=list(signal=setdf("soc$votes$noidprob$signal","noidsignalprob"),
                          wsignal=setdf("soc$votes$noidprob$wsignal","noidwsignalprob"),
                          context=setdf("soc$votes$noidprob$context","noidcxtprob"),
                          context1=setdf("soc$votes$noidprob$context1","noidcxtprob1"),
                          context2=setdf("soc$votes$noidprob$context2","noidcxtprob2"),
                          informed=setdf("soc$votes$noidprob$informed","noidinformedprob")), 
                b = list(idedist=append(soc$votes$b$idedist,b_idedist),
                         capacity=append(soc$votes$b$capacity,b_capacity),
                         incumbent=append(soc$votes$b$incumbent,b_incumbent),
                         pvi_nation=append(soc$votes$b$pvi_nation,b_pvi_nation),
                         pvi_local=append(soc$votes$b$pvi_local,b_pvi_local),
                         pid=append(soc$votes$b$pid,b_pid)),
                tracker = track_vote,
                vote_seed = append(soc$votes$vote_seed,vote_seed)))
  }

}

#' All One Stage Voting in Action 
#'
#' @description Conducting one stage of election.
#' 
#' @param soc The lists generated initially by \code{\link{initialsoc}}. 
#' @param commonfactor The common factor in ideology (intercept).
#' @param groupfactor The group influence on ideology.
#' @param personalfactor The random personal factor in  ideology (The SD of the random error to the model).
#' @param updatecommonfactor The movement in common factor.
#' @param updategroupfactor The movement in group influence on ideology.
#' @param updatepersonalfactor The SD for the movemebnt in personal influence on ideology.
#' @param groupfactor_type If \code{"constant"}, \code{groupfactor} is the constant coefficient for \code{group}.
#' If \code{"diverging"}, \code{groupfactor} is the SD for folded normal distribution that randomly 
#' assigns coeffcient for \code{group}.
#' @param knowmean1 The lowest average level of knowledge for a group.
#' @param knowmean2 The highest average level of knowledge for a group.
#' @param knowspread The spread (SD for normal distribution) of knowledg within a group.
#' @param partyblue_startloc Starting location of blue party ideology. 
#' @param partyred_startloc Starting location of red party ideology.
#' @param partystartloc_noise Random noise given to the starting locations of party ideologies.
#' @param partyblue_startcapacity Starting level of blue party capacity.
#' @param partyred_startcapacity Starting level of red party capacity.
#' @param partystartcapacity_noise Random noise given to the starting level of party capacities.
#' @param partyloc_noise Random noise in updating locations of party ideologies.
#' @param partycapacity_noise Random noise in updating party capacity.
#' @param partyblue_loc_direction Direction and expected quantity of the movement in blue party ideology. 
#' @param partyred_loc_direction Direction and expected quantity of the movement in red party ideology.
#' @param partyblue_capacity_direction Direction and expected quantity of the movement in blue party capacity.
#' @param partyred_capacity_direction Direction and expected quantity of the movement in red party capacity.
#' @param partyupdate_type The type of updating procedure. 
#' If \code{"oneshot"}, the update is always given over \code{startloc} and \code{startcapacity} defined in 
#' \code{\link{setparty}} at each election. If \code{"online"}, the update is given over the ideology and 
#' capacity realized in the previous election. If \code{"fromstart"} the update is given over the 
#' ideology and capacity realized in the first election of the simulation run.  
#' @param signaldifficulty The extent of noise level in signal.
#' @param signalbias_common The extent of common directional bias (i.e., positive value advantages red party and vice versa) in signal. 
#' @param signalbias_ingroup The extent of common ingroup bias (i.e., Evaluate the ingroup candidate better).
#' @param context_type "region" (default) or "network" (currently not used).
#' @param contextradius Radius to seach neighbors. Only applied when \code{context_type=="network"} (currently not used). 
#' @param contexttimetoagg How many previous elections to average. Either 1 or 2 (default). 
#' @param b_idedist The influence of ideological distance. The default = 0.198, # From ANES
#' @param b_capacity The influence of capacity The default = 0.716, # From ANES
#' @param b_incumbent The influence if incumebent. The default = 1.526, # From ANES (EDV result)
#' @param b_pvi_nation The influence of national PVI. The default = -9.827, # From ANES (EDV result)
#' @param b_pvi_local The influence of regional PVI. The default = 5.156, # From CCES pvi_state
#' @param b_pvi_sublocal The influence of sub-regional PVI. The default = 1.663, # From CCES pvi_state
#' @param b_pid The influence of group membership. The default = 1.215, # From ANES
#' @param burnin_redvoteprob_sd The spread for the "burn-in" period (i.e., first two elections) voting behaviors. 
#' See details.
#' @param burnin_redvoteprob_cor The correlation between voting and group in "burn-in" period voting behaviors.
#' See details.
#' @param seedval The random number seed.
#' 
#' @return List object with society characteristics, including following elements:
#' \itemize{
#'  \item \code{grid} The grid of voters (with group membership at each cell)
#'  \item \code{grid_region} The grid of region identifier at each cell.
#'  \item \code{grid_subregion} The grid of sub-region identifier at each cell. 
#'  \item \code{grid_coords} All grid coordinates. 
#'  \item \code{group} The numeric vector of group memberships 
#'  \item \code{region_stat} Aggregated group statistics by region. 
#'  \item \code{subregion_stat} Aggregated group statistics by region.
#'  \item \code{grid_votes2} The grid of voting decisions in t-2 election.
#'  \item \code{grid_votes1} The grid of voting decisions in t-1 election.
#'  \item \code{ideology} Created by \code{\link{setideology}}, may be updated by \code{\link{updateideology}}
#'  \item \code{knowledge} Created by \code{\link{setknowledge}}
#'  \item \code{context} Created by \code{\link{getregion}}
#'  \item \code{signal} Created by \code{\link{getsignal}}
#'  \item \code{incumbent} Character of incumbent party. Either \code{"blue"} or \code{"red"}.
#'  \item \code{votes} Created by \code{\link{letsvote}}
#' }
#' 
#' @details In the first two elections, the model does not have the sufficient 
#' number of previous elections to aggregate. If the society does not cotain the information of 
#' two previous elections, the hypothetical votes in first two election is 
#' radomly provided by the following procedure.
#' \itemize{
#'   \item Draw the reference value of voting for red for each group of voters by 
#'   multi-variate normal distribution with mean 0 and the covariance with 
#'   group membership (i.e., -1, or 1) set by \code{burnin_redvoteprob_cor}.
#'   \item Multiply the reference values by \code{burnin_redvoteprob_sd} and 
#'   add 0.5 to generate the probability of red vote for each group.
#'   \item Draw vote of voters in each group by the Bernouli distribution with probability 
#'   drawin in the previous step. 
#' } 
#' 
#' @importFrom MASS mvrnorm
#' 
#' @export

prefsoc <- function(soc,              
                    commonfactor=0.1,
                    groupfactor = 0.8, 
                    personalfactor = 1,
                    updatecommonfactor=0.1, 
                    updategroupfactor=0, 
                    updatepersonalfactor=0,
                    groupfactor_type="constant", # or diverging
                    knowmean1 = 0.65, # Also 0.6
                    knowmean2 = 0.65, # Also 0.7
                    knowspread = 0.85,
                    partyblue_startloc = -0.8, 
                    partyred_startloc = 0.8, 
                    partystartloc_noise=0.5,
                    partyblue_startcapacity=0,
                    partyred_startcapacity=0,
                    partystartcapacity_noise=0.5,
                    partyloc_noise = 0, # 0.5? 
                    partycapacity_noise = 0, # 1?
                    partyblue_loc_direction = 0, 
                    partyred_loc_direction = 0, 
                    partyblue_capacity_direction = 0, 
                    partyred_capacity_direction = 0,
                    partyupdate_type = "oneshot",
                    signaldifficulty = 3, 
                    signalbias_common = 0, 
                    signalbias_ingroup = 0, # 0 to 1
                    context_type = "region", # or "network"
                    contextradius = 5,
                    contexttimetoagg = 2,
                    b_idedist = 0.198, # From ANES
                    b_capacity = 0.716, # From ANES
                    b_incumbent = 1.526, # From ANES (EDV result)
                    b_pvi_nation = -9.827, # From ANES (EDV result)
                    b_pvi_local = 5.156, # From CCES pvi_state
                    b_pvi_sublocal = 1.663, # From CCES pvi_state
                    b_pid = 1.215, # From ANES
                    burnin_redvoteprob_sd = 0.045,
                    burnin_redvoteprob_cor = 0.19,
                    seedval = NULL
) {
  
  if (is.null(seedval)) {
    seedset <- NULL #sample(-2^15:2^15, 14)
  } else {
    if (length(seedval)==1) {
      set.seed(seedval)
      seedset <- sample(-2^15:2^15, 14)
    } else if (length(seedval)==14) {
      seedset <- seedval
    } else {
      stop("seedval must be NULL, length of 1, or length of 14!")
    }
  } 

  # Extract values
  allvalues <- as.numeric(soc$grid)
    ## voting Data
  personloc <- which(allvalues > 0)

  # Create data for first two elections
  if (length(soc[["grid_votes1"]])==0) {
    
    set.seed(seedset[14])
    tmp = mvrnorm(n=1, mu=c(0, 0), 
                  Sigma=matrix(c(1, 
                                 burnin_redvoteprob_cor, 
                                 burnin_redvoteprob_cor, 
                                 1), nrow=2), 
                  empirical=FALSE)
    probred1 <- (tmp[1]*burnin_redvoteprob_sd) + 0.5
    probred2 <- (tmp[2]*burnin_redvoteprob_sd) + 0.5
    
    # t-2 election 
    set.seed(seedset[1])
    randomvotes2 <- rbinom(length(allvalues[personloc]), 1, 
                           prob=probred1)
    votes2 <- rep(NA, length(allvalues))
    votes2[personloc] <- randomvotes2
    soc$grid_votes2 <- matrix(votes2, nrow=dim(soc$grid)[1], ncol=dim(soc$grid)[2])
    
    # t-1 election
    set.seed(seedset[2])
    randomvotes1 <- rbinom(length(allvalues[personloc]), 1, 
                           prob=probred2)
    votes1 <- rep(NA, length(allvalues))
    votes1[personloc] <- randomvotes1
    soc$grid_votes1 <- matrix(votes1, nrow=dim(soc$grid)[1], ncol=dim(soc$grid)[2])
    soc$incumbent <- ifelse(mean(randomvotes1)>0.5,"red","blue")
    
    # Set (True) Ideology
    soc$ideology <- setideology(soc$group, 
                                commonfactor, 
                                groupfactor, 
                                personalfactor, 
                                groupfactor_type,
                                b_seed=seedset[3], e_seed=seedset[4])
    ideology <- rep(NA, length(allvalues))
    ideology[personloc] <- soc$ideology$y
    soc$grid_ideology <- matrix(ideology, nrow=dim(soc$grid)[1], ncol=dim(soc$grid)[2])
    
    # Set Knowledge
    soc$knowledge <- setknowledge(soc$ideology$x, knowmean1, knowmean2, knowspread, kn_seed = seedset[5])
    knowledge <- rep(NA, length(allvalues))
    knowledge[personloc] <- soc$knowledge$p
    soc$grid_knowledge <-matrix(knowledge, nrow=dim(soc$grid)[1], ncol=dim(soc$grid)[2])
    
    # Set Party
    soc$party <- setparty(partyblue_startloc,
                          partyred_startloc,
                          partystartloc_noise,
                          partyblue_startcapacity,
                          partyred_startcapacity,
                          partystartcapacity_noise,
                          startincumbent = soc$incumbent,
                          blue_loc_seed=seedset[6],
                          red_loc_seed=seedset[7],
                          blue_capacity_seed=seedset[8],
                          red_capacity_seed=seedset[9])

    # Get Preference Signal
    soc$signal <- getsignal(soc$party, soc$ideology, soc$knowledge,
                            signaldifficulty, 
                            signalbias_common, 
                            signalbias_ingroup, #pastbias=NULL,
                            idedist_blue_signal_seed=seedset[10],
                            idedist_red_signal_seed=seedset[11],
                            capacity_signal_seed=seedset[12])
    
  } else {
    
    # Update Ideology
    soc$ideology <- updateideology(soc$ideology, 
                                   updatecommonfactor, 
                                   updategroupfactor, 
                                   updatepersonalfactor,
                                   b_seed=seedset[3], e_seed=seedset[4])
    
    # Update party
    soc$party <- updateparty(soc$party, soc$incumbent,
                             partyloc_noise, 
                             partycapacity_noise,
                             partyblue_loc_direction, 
                             partyred_loc_direction, 
                             partyblue_capacity_direction, 
                             partyred_capacity_direction,
                             partyupdate_type,
                             blue_loc_seed=seedset[6],red_loc_seed=seedset[7],
                             blue_capacity_seed=seedset[8],red_capacity_seed=seedset[9])
    
    # Get Preference Signal
    soc$signal <- getsignal(soc$party, soc$ideology, soc$knowledge,
                            signaldifficulty, 
                            signalbias_common, 
                            signalbias_ingroup, #pastbias=soc$signal$signal$bias,
                            idedist_blue_signal_seed=seedset[10],
                            idedist_red_signal_seed=seedset[11],
                            capacity_signal_seed=seedset[12])
    
  }

  
  # Get Context
  if (context_type=="network") {
    soc$context <- getcontext(soc, personloc, contextradius, contexttimetoagg)
  } else if (context_type=="region") {
    soc$context <- getregion(soc, personloc, contexttimetoagg)
  } else {stop("invalid context_type!")}
  
  # Vote in Election
  soc$votes <- letsvote(soc,
                        b_idedist = b_idedist, 
                        b_capacity = b_capacity,
                        b_incumbent = b_incumbent,
                        b_pvi_nation = b_pvi_nation,
                        b_pvi_local = b_pvi_local,
                        b_pvi_sublocal = b_pvi_sublocal,
                        b_pid = b_pid,
                        vote_seed=seedset[13])
  
  # Update grids
  soc$grid_votes2 <- soc$grid_votes1
  # Actual Votes
  votes <- rep(NA, length(allvalues))
  votes[personloc] <- as.numeric(soc$votes$vote$actual[,ncol(soc$votes$vote$actual)])
  soc$grid_votes1 <- matrix(votes, nrow=dim(soc$grid)[1], ncol=dim(soc$grid)[2])
  # Hypothetical Vote Probabilities (Signal)
  signal <- rep(NA, length(allvalues))
  signal[personloc] <- as.numeric(soc$votes$noidprob$signal[,ncol(soc$votes$vote$actual)])
  soc$grid_signal <-matrix(signal, nrow=dim(soc$grid)[1], ncol=dim(soc$grid)[2])
  # Hypothetical Vote Probabilities (Context)
  context <- rep(NA, length(allvalues))
  context[personloc] <- as.numeric(soc$votes$noidprob$context[,ncol(soc$votes$vote$actual)])
  soc$grid_context <-matrix(context, nrow=dim(soc$grid)[1], ncol=dim(soc$grid)[2])
  # Hypothetical Vote Probabilities (Informed)
  informed <- rep(NA, length(allvalues))
  informed[personloc] <- as.numeric(soc$votes$prob$informed[,ncol(soc$votes$vote$actual)])
  soc$grid_informed <-matrix(informed, nrow=dim(soc$grid)[1], ncol=dim(soc$grid)[2])
  
  # New Incumbent
  soc$incumbent <- ifelse(mean(soc$votes$vote$actual[,ncol(soc$votes$vote$actual)])>0.5,"red","blue")

  return(soc)
  
}

#' All One Stage Simulation in Action 
#'
#' @description Conducting one stage of election.
#' 
#' @param alike_threshold 1 - tolerance in the segregation game. The default = 0.2;
#' @param segregation_radius Search radius to look for neightbors The default = 1; 8 people around
#' @param signalbias_mean The average common directional bias (i.e., positive value advantages red 
#' party and vice versa) in signal. The default is 0.5.
#' @param signalbias_sd The SD common directional bias. The bias is drawn from normal distribution each time.
#' The default is 0.5.
#' @param signalbias_size The extent of fixed common directional bias. 
#' If assigned; mean and sd arguments above will be ignored.
#' @param signalbias_ingroup The extent of fixed common ingroup bias (i.e., Evaluate the ingroup candidate better).
#' The default is 0.
#' @param signaldifficulty The extent of noise level in signal. Given the default values of knowledge, 
#' ideology, and party position, the default difficulty of 3 produces idedist with SD 1.9; Consistent with ANES.
#' @param updatecommonfactor_sd The movement in common factor. SD. The default is 0.1.
#' @param updatecommonfactor_mean The movement in common factor. Mean. The default is 0.
#' @param updategroupfactor_sd The movement in group influence on ideology. SD. The default is 0.
#' @param updategroupfactor_mean The movement in group influence on ideology. Mean. The default is 0.
#' @param updatepersonalfactor_sd The SD for the movemebnt in personal influence on ideology. SD. The default is 0.
#' @param updatepersonalfactor_mean The SD for the movemebnt in personal influence on ideology. Mean. The default is 0.
#' @param partyloc_noise Random noise in updating locations of party ideologies. The default is 0.
#' @param partycapacity_noise Random noise in updating party capacity. The default is 0.5. 
#' The average score in ANES is 0.7 but it produces too much movement in voting.
#' @param partyblue_loc_direction Direction and expected quantity of the movement in blue party ideology. 
#' @param partyred_loc_direction Direction and expected quantity of the movement in red party ideology.
#' @param partycapacity_direction_mean The average direction and expected quantity of the movement in party capacity. 
#' (Positive vaue advantages red party and negative value advantages blue party).
#' @param partyupdate_type The type of updating procedure. 
#' If \code{"oneshot"}, the update is always given over \code{startloc} and \code{startcapacity} defined in 
#' \code{\link{setparty}} at each election. If \code{"online"}, the update is given over the ideology and 
#' capacity realized in the previous election. If \code{"fromstart"} the update is given over the 
#' ideology and capacity realized in the first election of the simulation run.  
#' @param context_type "region" (default) or "network" (currently not used).
#' @param contextradius Radius to seach neighbors. Only applied when \code{context_type=="network"} (currently not used). 
#' @param contexttimetoagg How many previous elections to average. Either 1 or 2 (default). 
#' @param b_idedist The influence of ideological distance. The default = 0.198, from ANES.
#' @param b_capacity The influence of capacity The default = 0.716, from ANES.
#' @param b_incumbent The influence if incumebent. The default = 1.526, from ANES (EDV result).
#' @param b_pvi_nation The influence of national PVI. The default = -9.827, from ANES (EDV result).
#' @param b_pvi_local The influence of regional PVI. The default = 5.156, from CCES pvi_state.
#' @param b_pvi_sublocal The influence of sub-regional PVI. The default = 1.663, from CCES pvi_state.
#' @param b_pid The influence of group membership. The default = 1.215, from ANES.
#' @param maxt_segregation Max time of movement in the segregation game. The default is 50.
#' @param n_elections Number of elections in each simulation run. The default is 30.
#' @param n_drop_elections Number of first N elections dropped from the quality assessment. The default is 5.
#' @param n_voters Interger. Number of agents. The default is 2300.
#' @param sizegroup Numeric Vector. Proportion of agents in each group. 
#' Length defines the number of groups. Must adds up to 1. The default is c(0.5,0.5).
#' @param spacelim Integer vector of length 2. Grid space size 
#' (number of rows and number of columns, respectively). The default is c(54,54).
#' @param regionlevel Integer vector of length 2. Level to Split Region (Horizontaly & Vertically).
#' The default is c(3,3).
#' @param subregionlevel Ingeger vector of length 2. Level to Split Sub-Region (Horizontaly & Vertically).
#' The default is c(9,9).
#' @param startcommonfactor The starting value of The common factor in ideology (intercept).
#' The default is NULL (draw randomly).
#' @param startgroupfactor The starting value of the group influence on ideology.
#' The default is 0.8 (from ANES, skewed slightly towards right in the reality). 
#' @param startpersonalfactor The starting value of the random personal factor in ideology 
#' (The SD of the random error to the model). The defauls is 1 (to get ideology distribution with SD 1, from ANES).
#' @param groupfactor_type If \code{"constant"}, \code{groupfactor} is the constant coefficient for \code{group}.
#' If \code{"diverging"}, \code{groupfactor} is the SD for folded normal distribution that randomly 
#' assigns coeffcient for \code{group}.
#' @param knowmean1 The lowest average level of knowledge for Group Blue/-1. 0.65 is the overall mean in ANES (stable).
#' @param knowmean2 The highest average level of knowledge for Group Red/1. 0.65 is the overall mean in ANES (stable).
#' @param knowspread The spread (SD for normal distribution) of knowledg within a group.
#' The default is 0.85 (To get approx 0.25 SD after conversion, referred to ANES).
#' @param partyblue_startloc Starting location of blue party ideology. The default is -0.8 (from ANES).
#' @param partyred_startloc Starting location of red party ideology. The default is 0.8 (from ANES).
#' @param partystartloc_noise Random noise given to the starting locations of party ideologies. The default is 0.
#' @param partyblue_startcapacity Starting level of blue party capacity. The default is 0.
#' @param partyred_startcapacity Starting level of red party capacity. The default is 0. 
#' @param partystartcapacity_noise Random noise given to the starting level of party capacities. The default is 0.5.
#' @param burnin_redvoteprob_sd The spread for the "burn-in" period (i.e., first two elections) voting behaviors. 
#' The default is 0.05, referred to SD in ANES National PVI. See details.
#' @param burnin_redvoteprob_cor The correlation between voting and group in "burn-in" period voting behaviors. 
#' the default is 0.19. referred to ANES. See details.
#' @param seedval The random number seed.
#' @param drawplot If True; draw plots of results
#' @param run_type The default is "all" (to run both segregation and election models) 
#' Can also be "segregation" or "election" (which runs only one of them).
#' @param segregated_soc The pre-defined society object. If not NULL and run_type == "election", predefined society is used in simulation.
#' 
#' @return List object with society characteristics, including following elements:
#' \itemize{
#'  \item \code{grid} The grid of voters (with group membership at each cell)
#'  \item \code{grid_region} The grid of region identifier at each cell.
#'  \item \code{grid_subregion} The grid of sub-region identifier at each cell. 
#'  \item \code{grid_coords} All grid coordinates. 
#'  \item \code{group} The numeric vector of group memberships 
#'  \item \code{region_stat} Aggregated group statistics by region. 
#'  \item \code{subregion_stat} Aggregated group statistics by region.
#'  \item \code{grid_votes2} The grid of voting decisions in t-2 election.
#'  \item \code{grid_votes1} The grid of voting decisions in t-1 election.
#'  \item \code{ideology} Created by \code{\link{setideology}}, may be updated by \code{\link{updateideology}}
#'  \item \code{knowledge} Created by \code{\link{setknowledge}}
#'  \item \code{context} Created by \code{\link{getregion}}
#'  \item \code{signal} Created by \code{\link{getsignal}}
#'  \item \code{incumbent} Character of incumbent party. Either \code{"blue"} or \code{"red"}.
#'  \item \code{votes} Created by \code{\link{letsvote}}
#'  \item \code{quality} Aggregated quality of decisions.
#' }
#' 
#' @details In the first two elections, the model does not have the sufficient 
#' number of previous elections to aggregate. If the society does not cotain the information of 
#' two previous elections, the hypothetical votes in first two election is 
#' radomly provided by the following procedure.
#' \itemize{
#'   \item Draw the reference value of voting for red for each group of voters by 
#'   multi-variate normal distribution with mean 0 and the covariance with 
#'   group membership (i.e., -1, or 1) set by \code{burnin_redvoteprob_cor}.
#'   \item Multiply the reference values by \code{burnin_redvoteprob_sd} and 
#'   add 0.5 to generate the probability of red vote for each group.
#'   \item Draw vote of voters in each group by the Bernouli distribution with probability 
#'   drawin in the previous step. 
#' } 
#' 
#' @importFrom RColorBrewer brewer.pal
#' 
#' @export

oneset <- function(# Segregation Parameters
                   alike_threshold = 0.2,
                   segregation_radius = 1, # 9*9 - 1 = 80 people around
                   # Voting Parameters
                   signalbias_mean = 0.2,
                   signalbias_sd = 0.2,
                   signalbias_size = NULL, # 1? If assigned, mean and sd arguments will be ignored.
                   signalbias_ingroup = 0.5, # Constant Bias to discount in-group disadvantage
                   signaldifficulty = 3, ## Given Set knowledge, ideology, party position, this difficulty produces idedist with SD 1.9, Consistent with ANES.
                   updatecommonfactor_sd = 0.1,
                   updatecommonfactor_mean = 0,
                   updategroupfactor_sd = 0,
                   updategroupfactor_mean = 0,
                   updatepersonalfactor_sd = 0,
                   updatepersonalfactor_mean = 0,
                   partyloc_noise = 0, ## No Movement
                   partycapacity_noise = 0.5, # 0.7 in ANES but it produces too much movements.
                   partyblue_loc_direction = 0, 
                   partyred_loc_direction = 0, 
                   partycapacity_direction_mean = 0, 
                   partyupdate_type = "oneshot", #
                   context_type = "region", # or "network"
                   contextradius = 5, # only used if context_type is "network"
                   contexttimetoagg = 2, # 1 or 2
                   # Coefficients for Voting Calculus
                   b_idedist = 0.198, # From ANES
                   b_capacity = 0.716, # From ANES
                   b_incumbent = 1.526, # From ANES (EDV result)
                   b_pvi_nation = -9.827, # From ANES (EDV result)
                   b_pvi_local = 5.156, # From CCES pvi_state
                   b_pvi_sublocal = 1.663, # From CCES pvi_state
                   b_pid = 1.215, # From ANES
                   # Setting of Simulation Runs
                   maxt_segregation=50,
                   n_elections=30,
                   n_drop_elections=5,
                   # Setting Up Society 
                   n_voters = 2300,
                   spacelim=c(54,54), 
                   sizegroup=c(0.5,0.5), # 50 % each
                   regionlevel=c(3,3), # 3 by 3 (9) regions
                   subregionlevel=c(9,9), # 9 by 9 (81) regions
                   # Starting Values of Preference
                   startcommonfactor=NULL, ## If NULL, use the value from Draw
                   startgroupfactor = 0.8, ### Reference to ANES, but reality is skewed toward right 
                   startpersonalfactor = 1, ### To get ideology distribution with SD=1.
                   groupfactor_type = "constant", # or diverging
                   knowmean1 = 0.65, # OR 0.6 Average Knowledge mean for Group Blue/-1 Referred from ANES (stable)
                   knowmean2 = 0.65, # OR 0.7 Average Knowledge mean for Group Red/1 Refeerred from ANES (stable)
                   knowspread = 0.85, ### To get approx 0.25 SD after conversion
                   partyblue_startloc = -0.8, ### Referred to Average Mean of ANES Pary Loc Over the Years  
                   partyred_startloc = 0.8, ### Referred to Average Mean of ANES Pary Loc Over the Years
                   partystartloc_noise = 0, # No Movement
                   partyblue_startcapacity = 0,
                   partyred_startcapacity = 0,
                   partystartcapacity_noise = 0.5, # 0.7 in ANES but it produces too much movements.
                   burnin_redvoteprob_sd = 0.05, # SD for the probability of Red vote two "burn-in" elections occured before the simulation. Referred to SD in ANES National PVI.
                   burnin_redvoteprob_cor = 0.19, # Correlation between Red vote share in two "burn-in" elections. 
                   # Other Options
                   seedval = NULL, # Set Seed Values Manually (Need 1 + maxt_segregation + 19*n_elections in length)
                   drawplot = FALSE, # If True, draw plots of results
                   run_type = "all", # can also be "segregation" or "election"
                   segregated_soc = NULL # If not NULL and run_type == "election", predefined society is used in simulation
                   )
{
  
  nseedforelection <- 19
  
  if (run_type!="election" | is.null(segregated_soc)) {

    seedquant <- 1 + maxt_segregation + nseedforelection*n_elections
    
    if (is.null(seedval)) {
      seedset <- NULL #sample(-2^15:2^15,seedquant)
    } else if (length(seedval)==seedquant) {
      seedset <- seedval
    } else if (length(seedval)==1) {
      set.seed(seedval)
      seedset <- sample(-2^15:2^15,seedquant)
    } else {
      stop("invalid seedval argument!")
    }

    # Initialize Society
    soc <- initialsoc(number=n_voters, 
                      sizegroup=sizegroup, 
                      spacelim=spacelim, 
                      regionlevel=regionlevel, 
                      subregionlevel=subregionlevel, 
                      seedval=seedset[1])
    soc <- happyregion(soc)
    soc <- happysubregion(soc)
    soc <- happysoc(soc, radius=segregation_radius, 
                    alike_threshold=alike_threshold, 
                    initiate=TRUE)
    soc$seed_all <- seedset
    soc$maxt_segregation <- maxt_segregation
    
  } else {
    
    
    soc <- segregated_soc
    if (!is.null(soc[["grid_votes1"]])) stop("Election Data Already Exist!")
    maxt_segregation <- soc$maxt_segregation
    seedquant <- 1 + maxt_segregation + nseedforelection*n_elections
    if (length(soc$seed_all) != seedquant) {
      if (is.null(seedval)) {
        seedset <- NULL #sample(-2^15:2^15,seedquant - (1+maxt_segregation))
      } else if (length(seedval)==(seedquant - (1+maxt_segregation))) {
        seedset <- seedval
      } else if (length(seedval)==1) {
        set.seed(seedval)
        seedset <- sample(-2^15:2^15, seedquant - (1+maxt_segregation))
      } else {
        stop("invalid seedval argument!")
      }
      seedset <- c(soc$seed_all[seq(1,(1+maxt_segregation),1)],seedset)
      soc$seed_all <- seedset
    } else {
      seedset <- soc$seed_all
    }
    
  }
  
  # Set colors for plotting
  setcolors <- c("#ffffb3","#bebada","#fb8072") #,"#8dd3c7"
  setcolors2 <- brewer.pal(n = 8, name = "PiYG")
  setcolors3 <- rev(brewer.pal(n = 8, name = "RdBu"))
  
  if (run_type %in% c("all","segregation")) {
    
    # Run Segregation Model
    t <- happy_tracker <- 0
    
    while(happy_tracker < 1 & t < maxt_segregation) {
      
      if (length(soc$happy_tracker)==0) {
        happy_tracker = 0
      } else {
        happy_tracker <- tail(soc$happy_tracker,1) 
      }
      t <- t+1
      
      cat(paste("Segregation Iteration",t, "\n"))
      
      print(soc$region_stat)
      
      cat("\n\n")
      
      # Move unhappy agents
      soc <- moveunhappy(soc, seedset[t+1]) 
      # Examine happiness of agents
      soc <- happyregion(soc)
      soc <- happysubregion(soc)
      soc <- happysoc(soc, radius=segregation_radius, 
                      alike_threshold=alike_threshold)
      
      # Record Number of Iteration
      soc$timetohappysoc <- t
      
      # plot Segregation Pattern
      if (drawplot==TRUE) {
        par(mfrow=c(2,2),mar=c(2,2,2,2))
        image(soc$grid,col=setcolors,axes=F)
        title(main = "Micro Segregation")
        image(soc$grid_ingroupshare,axes=F,col=setcolors2, zlim=c(0,1))
        title(main = "Regional Segregation")
        image(soc$grid_ingroupshare_sub,axes=F,col=setcolors2, zlim=c(0,1))
        title(main = "Sub-Regional Segregation")
        plot(runif(maxt_segregation,0,1),ylab="percent happy",xlab="time",col="white",ylim=c(0,1))
        lines(soc$happy_tracker,oma = c(0, 0, 2, 0),col="red")
        title(main = "Happiness Tracker")
      }
      
    }

  }
  if (run_type=="segregation") return(soc)
  
  elecseeds <- seedset[(maxt_segregation+1):length(seedset)]
  seedstart <- -(nseedforelection-1)
  
  for (t in 1:n_elections) {
    
    seedstart <- seedstart + nseedforelection
    seedend <- seedstart + (nseedforelection-1)
    useseed <- elecseeds[seedstart:seedend]
    
    cat(paste("Election",t,"\n"))
    
    # Common Factor Update
    if (updategroupfactor_mean!=0 | updatecommonfactor_sd!=0) {
      set.seed(useseed[16])
      commonfactor <- rnorm(1,mean=updatecommonfactor_mean,sd=updatecommonfactor_sd)
      cat(paste("Common Factor:",round(commonfactor,3),"\n"))
    } else {
      commonfactor <- 0
    }
    if (is.null(startcommonfactor)) startcommonfactor <- commonfactor

    # Group Factor Update
    if (updategroupfactor_mean!=0 | updategroupfactor_sd!=0) {
      set.seed(useseed[17])
      groupfactor <- rnorm(1,mean=updategroupfactor_mean,sd=updategroupfactor_sd)
      cat(paste("Group Factor:",round(groupfactor,3),"\n"))
    } else {
      groupfactor <- 0
    }
    if (is.null(startgroupfactor)) startgroupfactor <- groupfactor
    
    # Personal Factor Update
    if (updatepersonalfactor_mean!=0 | updatepersonalfactor_sd!=0) {
      set.seed(useseed[18])
      personalfactor <- rnorm(1,mean=updatepersonalfactor_mean,sd=updatepersonalfactor_sd)
      cat(paste("Personal Factor:",round(personalfactor,3),"\n"))
    } else {
      personalfactor <- 0
    }
    if (is.null(startpersonalfactor)) startpersonalfactor <- personalfactor
    
    # Party Capacity Direction
    if (partycapacity_direction_mean!=0 | partycapacity_noise!=0) {
      set.seed(useseed[19])
      partycapacity_direction <- rnorm(1,mean=partycapacity_direction_mean,
                                       sd=partycapacity_noise)
      partyblue_capacity_direction <- - partycapacity_direction
      partyred_capacity_direction <- partycapacity_direction
      cat(paste("Party Capacity Change (Red Advantage):", round(partycapacity_direction,3),"\n"))
    } else {
      partycapacity_direction <- 0
      partyblue_capacity_direction <- 0
      partyred_capacity_direction <- 0
    }

    # Draw Signal Bias
    if (!is.null(signalbias_size)) {
      set.seed(useseed[15])
      signalbias <- sample(c(-signalbias_size,signalbias_size), 1) 
      cat(paste("Signal Bias:",round(signalbias,3),"\n"))
    } else if (signalbias_mean=="balancecommon") {
      set.seed(useseed[15])
      signalbias <- rnorm(1,mean= -commonfactor, sd=signalbias_sd)
      cat(paste("Signal Bias:",round(signalbias,3),"\n"))
    } else if (signalbias_mean=="balancecapacity") {
      set.seed(useseed[15])
      signalbias <- rnorm(1,mean= -partycapacity_direction, sd=signalbias_sd)
      cat(paste("Signal Bias:",round(signalbias,3),"\n"))
    } else if (signalbias_mean!=0 | signalbias_sd!=0) {
      set.seed(useseed[15])
      signalbias <- rnorm(1,mean=signalbias_mean,sd=signalbias_sd)
      cat(paste("Signal Bias:",round(signalbias,3),"\n"))
    } else {
      signalbias <- 0
    }
    
    # Stage Simulation
    soc <- prefsoc(soc,                      
                   commonfactor = startcommonfactor,
                   groupfactor = startgroupfactor, 
                   personalfactor = startpersonalfactor,
                   updatecommonfactor=commonfactor, 
                   updategroupfactor=groupfactor, 
                   updatepersonalfactor=personalfactor,
                   groupfactor_type = groupfactor_type,
                   knowmean1 = knowmean1, 
                   knowmean2 = knowmean2, 
                   knowspread = knowspread, 
                   partyblue_startloc = partyblue_startloc,  
                   partyred_startloc = partyred_startloc,
                   partystartloc_noise = partystartloc_noise,
                   partyblue_startcapacity = partyblue_startcapacity,
                   partyred_startcapacity = partyred_startcapacity,
                   partystartcapacity_noise = partystartcapacity_noise, 
                   partyloc_noise = partyloc_noise, 
                   partycapacity_noise = partycapacity_noise,
                   partyblue_loc_direction = partyblue_loc_direction, 
                   partyred_loc_direction = partyred_loc_direction, 
                   partyblue_capacity_direction = partyblue_capacity_direction, 
                   partyred_capacity_direction = partyred_capacity_direction,
                   partyupdate_type = partyupdate_type,
                   signaldifficulty = signaldifficulty, 
                   signalbias_common = signalbias, 
                   signalbias_ingroup = signalbias_ingroup, 
                   context_type = context_type,
                   contextradius = contextradius,
                   contexttimetoagg = contexttimetoagg,
                   b_idedist = b_idedist, 
                   b_capacity = b_capacity,
                   b_incumbent = b_incumbent,
                   b_pvi_nation = b_pvi_nation,
                   b_pvi_local = b_pvi_local,
                   b_pvi_sublocal = b_pvi_sublocal,
                   b_pid = b_pid,
                   burnin_redvoteprob_sd = burnin_redvoteprob_sd,
                   burnin_redvoteprob_cor = burnin_redvoteprob_cor,
                   seedval = useseed[1:14]
    ) 

    if (drawplot==TRUE) {
      par(mfrow=c(1,3),mar=c(2,2,2,2))
      image(soc$grid_signal,col=setcolors3,axes=F,zlim=c(0,1))
      title(main = "Signal")
      image(soc$grid_context,col=setcolors3,axes=F,zlim=c(0,1))
      title(main = "Context")
      plot(runif(n_elections,0,1),ylab="percent red",xlab="time",col="white",ylim=c(0.35,0.65))
      abline(h=0.5, lty=2)
      lines(soc$votes$tracker$informed,oma = c(0, 0, 2, 0),col="black")
      lines(soc$votes$tracker$signal,oma = c(0, 0, 2, 0),col="green")
      lines(soc$votes$tracker$wsignal,oma = c(0, 0, 2, 0),col="orange")
      lines(soc$votes$tracker$context,oma = c(0, 0, 2, 0),col="blue")
      lines(soc$votes$tracker$actual,oma = c(0, 0, 2, 0),col="red")
      title(main = "Voting Tracker")
    }

  }
  
  # Quality of Decisions
  soc$quality <- list()
  
  # Vote Tracker
  tmp <- soc$votes$tracker
  tmp <- tmp[-c(1:n_drop_elections),] # Dropping first few elections
  # Deviation from Informed Vote Share
  deviation_informed <- colMeans(abs(tmp[,2:9]-tmp[,10]))
  deviation_noidinformed <- colMeans(abs(tmp[,2:9]-tmp[,11]))
  # Same winner as Informed Voter
  tmp <- t(apply(tmp,1, function(k) ifelse(k>=0.5,1,0)))
  samewinner_informed <- colMeans(abs(tmp[,2:9]==tmp[,10]))
  samewinner_noidinformed <- colMeans(abs(tmp[,2:9]==tmp[,11]))
  # Aggregate Preferences
  tmp <- (tmp-0.5)*2
  aggpref_informed <- colMeans(soc$votes$pref$informed + soc$votes$pref$pid)[-c(1:n_drop_elections)]
  aggpref_noidinformed <- colMeans(soc$votes$pref$informed)[-c(1:n_drop_elections)]
  # Aggregate Loss In Preferences
  aggloss_informed <- colMeans(aggpref_informed*cbind(tmp[,2:9])) - mean(aggpref_informed*tmp[,10])
  aggloss_noidinformed <- colMeans(aggpref_noidinformed*cbind(tmp[,2:9])) - mean(aggpref_informed*tmp[,11])
  # Aggregate Preferences & Loss by Group
  aggpref_informed_blue <- colMeans(soc$votes$pref$informed[soc$ideology$x==-1,] + 
                                      soc$votes$pref$pid[soc$ideology$x==-1,])[-c(1:n_drop_elections)]
  aggpref_noidinformed_blue <- colMeans(soc$votes$pref$informed[soc$ideology$x==-1,])[-c(1:n_drop_elections)]
  aggloss_informed_blue <- colMeans(aggpref_informed_blue*cbind(tmp[,2:9])) - mean(aggpref_informed_blue*tmp[,10])
  aggloss_noidinformed_blue <- colMeans(aggpref_noidinformed_blue*cbind(tmp[,2:9])) - mean(aggpref_informed_blue*tmp[,11])
  aggpref_informed_red <- colMeans(soc$votes$pref$informed[soc$ideology$x==1,] + 
                                     soc$votes$pref$pid[soc$ideology$x==1,])[-c(1:n_drop_elections)]
  aggpref_noidinformed_red <- colMeans(soc$votes$pref$informed[soc$ideology$x==1,])[-c(1:n_drop_elections)]
  aggloss_informed_red <- colMeans(aggpref_informed_red*cbind(tmp[,2:9])) - mean(aggpref_informed_red*tmp[,10])
  aggloss_noidinformed_red <- colMeans(aggpref_noidinformed_red*cbind(tmp[,2:9])) - mean(aggpref_informed_red*tmp[,11])
  
  soc$quality$informed <- data.frame(deviation=deviation_informed,
                                     samewinner=samewinner_informed,
                                     aggloss = aggloss_informed,
                                     aggloss_blue = aggloss_informed_blue,
                                     aggloss_red = aggloss_informed_red)
  soc$quality$noidinformed <- data.frame(deviation=deviation_noidinformed,
                                         samewinner=samewinner_noidinformed,
                                         aggloss = aggloss_noidinformed,
                                         aggloss_blue = aggloss_noidinformed_blue,
                                         aggloss_red = aggloss_noidinformed_red)
  
  cat("\nComparison with Informed Votes: \n")
  
  print(round(soc$quality$informed,3))
  
  cat("\n")
  
  cat("Comparison with (No PID) Informed Votes: \n")
  
  print(round(soc$quality$noidinformed,3))
  
  cat("\n\n")
  
  return(soc)

}


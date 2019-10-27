## Research Project Final Implementation
## Author: Shubham S. Diwe
## Genetic Algorithm for Hierarchical Optimization of Time Series Forecasts

#######----------Importing required packages----------#######
library(hts)
library(forecast)
library(TSPred)
library(rlist)
library(dplyr)
library(ggplot2)
set.seed(12345)
#######----------Reading Data and Preprocessing----------#######
data = infantgts

df <- readRDS('mtsSLS_MM_SAMPLE.RData') # as.matrix(allts(data))
dim(df)
# year <- seq(1933, 2003)
# 
# df <- cbind(year, df)

df[is.na(df)] <- 0
# df <- df[, c(1, 13:28)]
# df <- na.omit(df)

R <- dim(df)[1]
C <- dim(df)[2]

R.tr <- round(0.85*R, 0)                                                                               # Rows in training data

train <- df[c(1:R.tr), ]                                                                               # Train data
test <- df[c((R.tr+1):R), ]                                                                            # Test data

fcp <- nrow(test)                                                                                      # forecast period

#######----------Bottom Level Forecasting with simple ETS----------#######
bts.fc <- matrix(0, nrow = fcp, ncol = C)                                                              # Bottom Level Time series Forecast

for (num in 1:C) {
  
  bts.fc[, num] <- forecast(ets(train[, num]), h = fcp)$mean
  
}

# bts.fc

# ets.err <- sapply((2:C), function(y) {sMAPE(test[, y], bts.fc[, y])})
# mean.ets.err <- mean(ets.err)
# mean.ets.err

ets.err <- sapply((2:C), function(y) {sMAPE(test[, y], bts.fc[, y])})
mean.ets.err <- mean(ets.err)
mean.ets.err

#######----------Prerequisites for the Hierarchy----------#######
num.node <- 4                                                                                       # Number of nodes in the 2nd level

mat.node <- matrix(0, nrow = num.node, ncol = (C-1))                                                # Matrix for hierarchy structure

vec.node <- vector(mode = 'numeric', length = num.node)                                             # Number of time series in each node

#######----------Generation of Hierarchies----------#######
hierarchy <- function(nr = num.node, nc = (C-1)) {                                                  # Create Matrix for Hierarchy Structure
  
  nodes <- matrix(0, nc, nr)
  
  nodes[cbind(1:nc, sample(nr, nc, TRUE))] <- 1
  
  nodes <- t(nodes)
  
  return(nodes)
  
} # hierarchy

#######----------Sorting Data as per the Hierarchies----------#######
sort.data <- function(nm = mat.node, data = train) {
  
  hi <- matrix(data[, 1])
  
  for (rn in 1:num.node) {                                                                          # For every row of the matrix 
    
    for (cn in 1:(C-1)) {                                                                           # For every included time series
      
      if (nm[rn, cn] == 1) {
        
        hi <- cbind(hi, data[, (cn+1)])
        
      } # if
      
    } # for cn
    
    vec.node[rn] <- sum(nm[rn, ])
    
  } # for rn
  
  return(list(nm, hi, vec.node))
  
} # sort.data

#######----------Creating Set of Hierarchies for Initiation----------#######
num.hi <- 30

hierarchies <- list()

for (nh in 1:num.hi) {
  
  nih <- hierarchy()
  
  if (all(rowSums(nih) >= 1)) {
    hierarchies[[nh]] <- sort.data(nm = nih)
  }
  
  hierarchies <- Filter(Negate(is.null), hierarchies)
  
}

# str(hierarchies)
# hierarchies[1]

#######----------Fitness or Selection Function -> Evaluation----------#######
selection <- function(n.ms, h.ts, v.ns, node = num.node) {
  
  sel.hts <- hts(h.ts[, 2:C], nodes = list(node, v.ns))
  
  agg.hts <- aggts(sel.hts, levels = c(0, 1))
  
  agg.fc <- matrix(0, nrow = fcp, ncol = ncol(agg.hts))
  
  for (jcol in 1:ncol(agg.hts)) {
    
    agg.fc[, jcol] <- forecast(ets(agg.hts[, jcol]), h = fcp)$mean
    
  } # for forecasting each aggregated series
  
  temp.mat <- sort.data(nm = n.ms, data = bts.fc)[[2]]
  
  merg.mat <- cbind(agg.fc, temp.mat[, 2:C])
  
  fc.mat <- combinef(fcasts = merg.mat, nodes = sel.hts$nodes, keep = 'bottom')
  
  re.test <- sort.data(nm = n.ms, data = test)[[2]]                                     # Rearrange the test time series as per new hierarchy
  
  s.mape <- mean(sapply(1:ncol(n.ms), function(nj) {sMAPE(re.test[, (nj+1)], fc.mat[, nj])}))
  
  return(s.mape)
  
} # Selection Function

#######----------Crossover Function----------#######
crossover <- function(h1, h2) {
  
  new.h <- matrix(0, nrow = num.node, ncol = ncol(h1))
  
  s.num <- round(runif(1, 0, 1), digits = 0)
  
  if (s.num == 1) {
    
    for (colj in 1:ncol(new.h)) {
      
      r.num <- round(runif(1, 0, 1), digits = 0)
      
      if (r.num == 1) {
        
        new.h[, colj] <- h1[, colj]
          
      } # Replace with column of H1
      else {
        
        new.h[, colj] <- h2[, colj]
        
      } # Replace with column of H2
      
    } # for every column
    
  } # Columnwise Crossover
  else {
    
    hf.c = floor(ncol(h1)/2)
    
    new.h[, c(1:hf.c)] <- h1[, c(1:hf.c)]
    new.h[, c((hf.c+1):ncol(h1))] <- h2[, c((hf.c+1):ncol(h1))]
    
  } # Merge halves of two hierarchies
  
  if ((all(colSums(new.h) == 1)) & (all(rowSums(new.h) >= 1))) {
    
    return(new.h)
  
  } # return crossovered hierarchy
  else {
    
    return(h1)
    
  } # return best of the original hierarchy
  
} # CrossOver Function

#######----------Mutation Function----------#######
mutation <- function(hi1) {
  
  new.hi = matrix(0, nrow = num.node, ncol = ncol(hi1))
  
  for (cj in 1:ncol(hi1)) {
    
    m.num <- round(runif(1, 0, 1), digits = 0)
    
    n.num <- sample(1:nrow(hi1), 1)
  
    if (m.num == 1) {
      
      new.hi[n.num, cj] <- 1
      
    } # Mutate the column 
    else {
      
      new.hi[, cj] <- hi1[, cj]
      
    } # return original column
    
  } # For each column/time series
  
  if ((all(colSums(new.hi) == 1)) & (all(rowSums(new.hi) >= 1))) {
    
    return(new.hi)
    
  }
  else {
    
    return(hi1)
    
  }
  
} # Mutation FUnction

#######----------Sort the  Hierarcies based on the Error----------#######
sort.err <- function(h.list) {
  
  #merged.list <- lapply(1:length(h.list), function(x) {list.append(h.list[x], err.list[x])})
  
  sorted.list <- h.list[order(sapply(h.list, '[[', 2))]
  
  return(sorted.list)
  
}

#######----------Genetic Algorithm for Hierarchies----------#######
gahts <- function(tsh = hierarchies, iters = 100) {
  
  obj.co <- list()
  
  min.err <- list()
  
  ret.obj <- list()
  
  for (xi in 1:length(tsh)) {
    
    nx1 <- tsh[[xi]][[1]]
    hx1 <- tsh[[xi]][[2]]
    vx1 <- tsh[[xi]][[3]]
    
    sel.x1 <- selection(n.ms = nx1, h.ts = hx1, v.ns = vx1)
    
    obj.co <- list.append(obj.co, list(nx1, sel.x1))
    
  }
  
  for (ri in 1:iters) {
    # print(ri)
    yum <- (ri %% 100)
    
    for (ci in 1:((length(tsh) - 1))) {
      # print(ri, ci)
      one <- ci
      two <- (ci + 1)
      
      # if (ci < length(tsh)) {
      #   
        h.one <- tsh[[ci]][[1]]
        h.two <- tsh[[(ci + 1)]][[1]]
      #   
      # } # Normal sequential crossover until the last hierarchy
      # else {
      #   
      #   h.one <- tsh[[ci]][[1]]
      #   h.two <- tsh[[1]][[1]]
      #   
      # } # Last hierarchy and the first hierarchy crossover
      
      h.co <- crossover(h1 = h.one, h2 = h.two)                                # Crossver of hierarchies
      
      h.sd <- sort.data(nm = h.co)                                             # Sort data as per the crossovered hierarchy
      
      h.mat <- h.sd[[2]]
      n.vec <- h.sd[[3]]
      
      h.sl <- selection(n.ms = h.co, h.ts = h.mat, v.ns = n.vec)               # Compute error on the crossover hierarchy
      
      co.app <- list(h.co, h.sl)
      
      obj.co <- list.append(obj.co, co.app)                                    # Store the crossovered hierarchy
      # obj.co[[ci]][[2]] <- h.sl                                              # Store the corresponding error
      
      m.prob <- runif(1, 0, 1)                                                 # Create a random number to decide for mutation
      
      if (m.prob <= 0.1) {                                                    # Execute only if decided for mutation
        
        obj.mu <- list()
        
        h.mu <- mutation(h.one)                                                # Mutation of hierarchy
        
        sd.h <- sort.data(nm = h.co)                                           # Sort data as per the mutated hierarchy
        
        mat.h <- sd.h[[2]]
        vec.n <- sd.h[[3]]
        
        sl.h <- selection(n.ms = h.mu, h.ts = mat.h, v.ns = vec.n)             # Compute error on mutated hierarchy
        # print(ri)
        mu.app <- list(h.mu, sl.h)
        
        obj.mu <- list.append(obj.mu, mu.app)                                  # Store mutated hierarchy
        # obj.mu[[1]][[2]] <- sl.h                                             # Store the corresponding error
        obj.mu <- unlist(obj.mu, recursive = F)
        
        obj.co <- list.append(obj.co, obj.mu)                                  # Add to the storage
        
      } # 5% chance of mutation
      
    } # for each of the created hierarchies
    
    h.se <- sort.err(h.list = obj.co)                                          # Sort list based on the errors
    
    min.err <- list.append(min.err, h.se)                                      # Store all sorted results                              
    
    ret.obj[[ri]] <- h.se[[1]]                                                 # Store the best hierarchy for this iteration
    
    if (yum == 0){                                                             # After every 100 iterations
      
      sort.min <- min.err[order(sapply(min.err, '[[', 2))][c(1:20)]            # Collect best 20 hierarchies
      
      list.remove(min.err, c(1:100))                                           # Empty the storage
      
      h.se <- list.append(h.se, sort.min)                                      # Reintroduce the top 20 previous hierarchies
      
      ret.obj <- sort.min                                                      # Keep the best 20 hierarchies so for
      
    } # Reintroduce 20 of last 100 minimum error hierarchies
    print(ri)
    tsh = h.se[c(1:20)]                                                        # Reassign the hierarchies object to new hierarcchies        
    
  } # Number of iterations or till convergence
  
  return(ret.obj)
  
} # genetic algorithm for hierarchies

htsga <- gahts(tsh = hierarchies, iters = 5)

length(htsga)

htsga.er <- htsga[order(sapply(htsga, '[[', 2))]
htsga.er[[1]][[2]]

#######----------References----------#######
## 1. Generate binary matrix [https://stackoverflow.com/questions/41945748/generate-a-random-matrix-in-r-with-m-columns-and-n-rows-where-rows-sum-to-1]
## 2. Remove Nulls from a list [https://stackoverflow.com/questions/16896376/extract-non-null-elements-from-a-list-in-r]
## 3. Sort list of lists [https://stackoverflow.com/questions/53265425/sort-a-list-of-lists-based-on-an-internal-list-variable-in-r]
## 4. MASE in R [https://www.rdocumentation.org/packages/Metrics/versions/0.1.4/topics/mase]
## 5. HTS in R [https://cran.r-project.org/web/packages/hts/hts.pdf]
## 6. Forecast.gts from HTS [https://www.rdocumentation.org/packages/hts/versions/5.1.5/topics/forecast.gts]
## 7. Append to list of lists [https://stackoverflow.com/questions/14848172/appending-a-list-to-a-list-of-lists-in-r]
## 8. Extract values from uniform list [https://stackoverflow.com/questions/6907102/how-do-i-extract-values-from-uniform-list-in-r]
## 9. Extract forecasts only [https://stackoverflow.com/questions/50435939/extract-only-the-forecasted-values-from-forecast]



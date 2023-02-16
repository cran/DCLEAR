## ----knitr_options, echo=FALSE, results=FALSE---------------------------------
library(knitr)
opts_chunk$set(fig.width = 12)

## ----loading, include=FALSE---------------------------------------------------
library(DCLEAR)
library(phangorn)
library(dplyr)
library(purrr)
library(ape)


## ----simulate1----------------------------------------------------------------
m = 30  # number of targets
acell = as.integer(rep(1,m)) ## mother cell with an unmutated state
mu_d = 0.1  ## mutation rate (ex 0.01 to 0.3)
d = 3  ## number of cell division
n_s = 5 ##  the number of possible mutational outcomes (1:n_s)
p_d = 0.05 ## dropout probability

nmstrings = c( '0', '-', LETTERS[1:n_s] ) 
sim_tree = list()
sim_tree[[1]] = acell

k = 2
while(k < 2^d) {
    ## codes for point mutation
    mother_cell = sim_tree[[k%/%2]]
    mu_loc = runif(m) < mu_d
    mutation_cites = (mother_cell == 1) &  mu_loc
    n_mut = sum(mutation_cites)
    if (n_mut > 0) {
        mother_cell[mutation_cites] = as.integer(sample(n_s, n_mut, replace = T)+2)
    }
    
    ## codes for dropout
    dropout_cites = runif(m) < p_d
    if (sum(dropout_cites) > 2 ) {
        dropout_between = sample(which(dropout_cites), 2 )
        mother_cell[dropout_between[1]:dropout_between[2]] = as.integer(2)
    }
    
    sim_tree[[k]] = mother_cell
    k = k+1
}

## ----simulate2----------------------------------------------------------------
1:7 %>% map(~paste(nmstrings[sim_tree[[.]]], collapse=""))

## ----parameters---------------------------------------------------------------
set.seed(1)
mu_d1 = c( 30, 20, 10, 5, 5, 1, 0.01, 0.001)
mu_d1 = mu_d1/sum(mu_d1)
simn = 100 # number of cell samples
m = 200  ## number of targets
mu_d = 0.03 # mutation rate
d = 12 # number of cell division
p_d = 0.005 # dropout probability

## ----simulatedata-------------------------------------------------------------
sD = sim_seqdata(sim_n = simn, m = m, mu_d = mu_d, d = d, n_s = length(mu_d1), outcome_prob = mu_d1, p_d = p_d )

## ----simulatedata_seq---------------------------------------------------------
sD$seqs

## ----simulatedata_tree--------------------------------------------------------
class(sD$tree)

## ----hamming------------------------------------------------------------------
distH = dist.hamming(sD$seqs)

## ----treeconst----------------------------------------------------------------
TreeNJ = NJ(distH)
TreeFM = fastme.ols(distH)

## ----eval---------------------------------------------------------------------
print( RF.dist(TreeNJ, sD$tree, normalize = TRUE) )
print( RF.dist(TreeFM, sD$tree, normalize = TRUE) )

## ----initial_wh1--------------------------------------------------------------
InfoW = -log(mu_d1)  ## entropy as their weights
InfoW[1:2] = 1   ## weight for initial state is 1 , weight for '-' is set to 1. 
InfoW[3:7] = 4.5 ## Constant weight for outcome states
dist_wh1 = WH(sD$seqs, InfoW)
#dist_wh1 = dist_weighted_hamming(sD$seqs, InfoW, dropout = FALSE)


## ----treeconst_wh1------------------------------------------------------------
TreeNJ_wh1 = NJ(dist_wh1)
TreeFM_wh1 = fastme.ols(dist_wh1)

## ----eval_wh1-----------------------------------------------------------------
print( RF.dist(TreeNJ_wh1, sD$tree, normalize = TRUE) )
print( RF.dist(TreeFM_wh1, sD$tree, normalize = TRUE) )

## ----initial_wh2--------------------------------------------------------------
InfoW = -log(mu_d1)  ## entropy as their weights
InfoW[1] = 1  # weight for initial state is 1
InfoW[2] = 12  # weight for interval dropout, '----', is set to 12. 
InfoW[3:7] = 3 # Constant weight for outcome states
dist_wh2 = WH(sD$seqs, InfoW, dropout=TRUE)

## ----treeconst_wh2------------------------------------------------------------
TreeNJ_wh2 = NJ(dist_wh2)
TreeFM_wh2 = fastme.ols(dist_wh2)

## ----eval_wh2-----------------------------------------------------------------
print( RF.dist(TreeNJ_wh2, sD$tree, normalize = TRUE) )
print( RF.dist(TreeFM_wh2, sD$tree, normalize = TRUE) )


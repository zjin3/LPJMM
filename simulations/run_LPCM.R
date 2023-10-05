
folder = "output_LPCM"
if (dir.exists(folder) == FALSE) { dir.create(folder) }

n.GoF = 1000  # number of simulations used in Goodness-of-fit test

load("data/sim1_data.rda")

library(latentnet)
# transfer to edge attributes:
x.edge = matrix(NA, N, N)
for (i in 1:N) {
  for (j in 1:N) {
    x.edge[i,j] = abs(x[i] - x[j])
  }
}
library(network)
my_seed = 246
my_net <- as.network(Adj)
my_fit <- ergmm(my_net ~ edgecov(x.edge) + euclidean(d = 2, G = 5), 
                tofit = c("mcmc"),
                seed = set.seed(my_seed))


library("mclust")
rand_index_chain = apply(my_fit$sample$Z.K, MARGIN = c(1), 
                         FUN = mclust::adjustedRandIndex,
                         y = true.g)
mean_rand <- mean(rand_index_chain)
mean_rand

cl_summary_est <- function(cl.mat = my_fit$sample$Z.K, N=100){
  clustering_summary = array(dim = c(4, N))
  est.g1 = GreedyEPL::MinimiseEPL(cl.mat)$decision
  # calculate Posterior Similarity Matrix:
  my.psm <- mcclust::comp.psm(cl.mat)
  est.g2 = mcclust::minbinder(my.psm)$cl
  est.g3 = mcclust::maxpear(my.psm)$cl
  est.g4 = mcclust::medv(my.psm)
  # collect all:
  all.g = rbind(est.g1, est.g2, est.g3, est.g4)
  method = c("GreedyEPL", "minbinder", "maxpear", "medv")
  rownames(all.g) = method
  
  clustering_summary = all.g
  return(clustering_summary)
}

cluster_sum <- cl_summary_est()

GreedyEPL_rand <- mclust::adjustedRandIndex(cluster_sum["GreedyEPL",], true.g)
MinBinder_rand <- mclust::adjustedRandIndex(cluster_sum["minbinder",], true.g)
MaxPEAR_rand   <- mclust::adjustedRandIndex(cluster_sum["maxpear",], true.g)
Medv_rand      <- mclust::adjustedRandIndex(cluster_sum["medv",], true.g)

rand_index_table <- as.matrix(apply(cluster_sum, MARGIN = 1, FUN = mclust::adjustedRandIndex, true.g))
colnames(rand_index_table) <- "adjusted Rand index"
rand_index_table


library(tidyverse)
n_categ = apply(cluster_sum, MARGIN = 1, FUN = n_distinct)
rand_index_table = cbind(rand_index_table, n_categ)
rand_index_table

z.chain = my_fit$sample$Z
g.chain = my_fit$sample$Z.K
beta.chain = my_fit$sample$beta
save(z.chain,
     g.chain,
     beta.chain,
     cluster_sum,
     file = glue::glue("{folder}/Markov_chain_and_cl_summary.rda"))


# Plots 
mycol = c("deepskyblue", "orange", "lavenderblush3", "forestgreen",
          "coral", "cyan3", "darkblue", "darksalmon",
          "darkkhaki", "violetred", "lemonchiffon3",
          "darkmagenta")


z.chain = my_fit$sample$Z
n.sim   = dim(z.chain)[1]
Z.ref   = Z
z.transform = array(NA, dim = dim(z.chain))  ### [n.sim, N, K]
for(j in 1:n.sim){  
  z.transform[j,,] <- MCMCpack::procrustes(z.chain[j,,], 
                                           Z.ref, 
                                           translation = T,
                                           dilation = T)$X.new
}

Z.mean = apply(z.transform, MARGIN = c(2,3), FUN=mean)

# par(mfrow=c(2,3))
# true.g:
plot(Z.mean, xlab = "Z_1", ylab = "Z_2", pch = 19, col = mycol[true.g],
     main = glue::glue("estimated z (true g)"))
text(Z.mean, as.character(1:N), cex=0.6, pos=4, col="gray")

# GreedyEPL:
plot(Z.mean, xlab = "Z_1", ylab = "Z_2", pch = 19,
     col = mycol[cluster_sum["GreedyEPL",]],
     main = glue::glue("estimated z (GreedyEPL)"))
text(Z.mean, as.character(1:N), cex=0.6, pos=4, col="gray")

# minbinder:
plot(Z.mean, xlab = "Z_1", ylab = "Z_2", pch = 19,
     col = mycol[cluster_sum["minbinder",]],
     main = glue::glue("estimated z (minbinder)"))
text(Z.mean, as.character(1:N), cex=0.6, pos=4, col="gray")

# maxpear:
plot(Z.mean, xlab = "Z_1", ylab = "Z_2", pch = 19,
     col = mycol[cluster_sum["maxpear",]],
     main = glue::glue("estimated z (maxpear)"))
text(Z.mean, as.character(1:N), cex=0.6, pos=4, col="gray")




# goodness-of-fit test:
library(igraph)
density = c()
transit = c()
assorta = c()

set.seed(my_seed)

for (i in 1:n.GoF){
  sim_net <- simulate(my_fit)
  adj_mat = igraph::graph_from_adjacency_matrix(as.matrix.network(sim_net))
  density = c(density, network::network.density(sim_net))
  transit = c(transit, igraph::transitivity(adj_mat))
  assorta = c(assorta, igraph::assortativity_nominal(adj_mat, types = true.g))
}
dens_mean = round(mean(density), 4)
# 0.15313,  true = 0.15
tran_mean = round(mean(transit), 4)
# 0.54656, true = 0.50
asso_mean = round(mean(assorta), 4)
# 0.5469478,  true = 0.55

GoF_table = matrix(c(dens_mean, tran_mean, asso_mean,
                     0.1531, 0.5049, 0.5512), 
                   3, 2, byrow = FALSE)
row.names(GoF_table) <- c("density", "transitivity", "assortativity")
colnames(GoF_table) <- c("mean", "truth")
GoF_table

GoF_chain = rbind(density, transit, assorta)
row.names(GoF_chain) = c("density", "transitivity", "assortativity")



z_dist_fun <- function(Z.est, Z.true = Z, N = 100) {
  sum_e_dist = 0
  for (i in 1:N){
    dist_e = norm(Z.est[i,] - Z.true[i,], type = "2")
    sum_e_dist = sum_e_dist + dist_e
  }
  return(sum_e_dist)
}
sum_euc_dist = z_dist_fun(Z.mean, Z)
sum_euc_dist

save(rand_index_table,
     GoF_table,
     GoF_chain,
     mean_rand,
     Z.mean,
     sum_euc_dist,
     file = glue::glue("{folder}/GoF_rand_and_Zmean.rda"))

# Adjusted rand index:
print(rand_index_table)

# Goodness-of-fit results:
print(GoF_table)

# sum of Euclidean distance of estimated latent locations:
print(sum_euc_dist)





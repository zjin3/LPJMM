
folder = "output_CSBM"
if (dir.exists(folder) == FALSE) { dir.create(folder) }

library(sbm)
load("data/sim1_data.rda")
# convert nodal to edge attributes
cov.mat = matrix(NA, 100 ,100)
for (i in 1:100) {
  for (j in 1:100) {
    cov.mat[i,j] = abs(x[i]-x[j])
  }
}

my.sbm2 = estimateSimpleSBM(Adj, model = "bernoulli",
                            covariates = list(cov.mat),
                            estimOptions = list(plot = FALSE,
                                                nbBlocksRange = c(5,5)))
sbm.g2 = my.sbm2$memberships

# adjusted rand index:
adj_rand_index <- mclust::adjustedRandIndex(true.g, sbm.g2)
adj_rand_index
# 0.707269

col_temp = c("orange", "coral", "lavenderblush3", 
             "deepskyblue", "forestgreen", "darkblue")

plot(Z, xlab = "Z_1", ylab = "Z_2", pch = 19, 
     col = col_temp[sbm.g2], 
     main = glue::glue("based on SBM and true Z"))
text(Z, as.character(1:N), cex=0.6, pos=4, col="gray")

cl_sbm = my.sbm2$memberships

save(my.sbm2, cl_sbm, 
     file = glue::glue("{folder}/sbm_res.rda"))




# goodness-of-fit test
set.seed(123)
n.sim = 1000

library(network)
library(igraph)

P_sbm = my.sbm2$predict()

density.chain = c()
tran.chain = c()
asot.chain = c()
for (n in 1:n.sim) {
  Adj.t = matrix(NA, N, N)
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      Adj.t[i, j] = sample(0:1, 1, prob = c(1 - P_sbm[i, j], P_sbm[i, j]))
      Adj.t[j, i] = Adj.t[i, j]
    }
  }
  diag(Adj.t) = rep(0, N)
  
  G.t = as.network(Adj.t)
  density.t = network.density(G.t)
  density.chain = c(density.chain, density.t)
  
  G.t = igraph::graph_from_adjacency_matrix(Adj.t, diag = F)
  tran.t = transitivity(G.t)
  tran.chain = c(tran.chain, tran.t)
  
  asot.t = assortativity_nominal(G.t, types = true.g)
  asot.chain = c(asot.chain, asot.t)
}

density.sbm      = mean(density.chain)
transitivity.sbm = mean(tran.chain)
asot.sbm         = mean(asot.chain)

g = as.network(Adj)
density.true = network::network.density(g)

g = igraph::graph_from_adjacency_matrix(Adj, diag = F)
transitivity.true = igraph::transitivity(g)
asot.true = assortativity_nominal(g, types = true.g)


GoF_table = matrix(c(density.true, transitivity.true, asot.true,
                     density.sbm, transitivity.sbm, asot.sbm), 
                   nrow = 3, ncol = 2)
rownames(GoF_table) = c("density", "transitivity", "assortativity")
colnames(GoF_table) = c("truth", "mean")
GoF_table

save(density.chain = density.chain,
     transit.chain = tran.chain,
     assorta.chain = asot.chain,
     GoF_table,
     adj_rand_index,
     file = glue::glue("{folder}/GoF_sbm.rda"))

# adjusted rand index:
print(adj_rand_index)

# goodness-of-test result
print(GoF_table)

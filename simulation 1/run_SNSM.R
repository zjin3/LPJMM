
model_name = "SNSM"

n.adapt  = 20000
n.update = 20000
n.sim    = 10000
n.chains = 1
iter.keep= 1:n.sim
thin     = 10

H_test = "na"

# save results:
folder = glue::glue("output_{model_name}")
if (dir.exists(folder) == FALSE) { dir.create(folder) }

# load data
load("data/sim1_data.rda")
Adj.layer = array(Adj, dim = c(N, N, 1))
L = 1                       # single layer
isdirected = "undirected"   # symmetric network

# packages, output txt file and jags setup:
source("source/specific_setup.R")

# generate Markov chain:
my_seed = 246
my_jags <- jags.model(model_path, data_prior, 
                      inits = list(".RNG.name" = "base::Marsaglia-Multicarry",
                                   ".RNG.seed" = my_seed, a = 5, b = -2),
                      n.chains = n.chains,
                      n.adapt = n.adapt)
update(my_jags, n.iter = n.update)
jags_res <- jags.samples(my_jags, variable.names = record_para, n.iter = n.sim)
# save(jags_res, file = glue("{folder}/jags_res.rda"))


# CI table and estimated latent position:
source("source/analyze_jags_single_layer.R")


# post-process group membership:
if (exists("H") == T & exists("true.g") == T & "categ" %in% record_para) {
  source("source/cluster_analysis.R")
}


# goodness-of-fit test:
if ("padj"  %in% record_para) {
  source("source/goodness_of_fit.R")
}


### calculate log likelihood
source("source/loglikelihood_SNSM.R")

print(sessionInfo())




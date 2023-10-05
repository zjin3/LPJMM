model_name = "LPJMM"

n.adapt  = 20000
n.update = 20000
n.sim    = 10000
n.chains = 1
iter.keep= 1:n.sim
thin     = 10

H_test = 5

folder = glue::glue("output_{model_name}")
if (dir.exists(folder) == FALSE) {
  dir.create(folder)
}

# load data
load("data/sim4_data.rda")
L = dim(Adj)[3]
Adj.layer = Adj
isdirected = "undirected"   

# jags setup:
source("source/specific_setup.R")

# generate Markov chain:
my_seed = 110
my_jags <- jags.model(model_path, data_prior, 
                      inits = list(".RNG.name" = "base::Marsaglia-Multicarry",
                                   ".RNG.seed" = my_seed,
                                   a = c(5,  3),
                                   b = c(-2, 1)),
                      n.chains = n.chains,
                      n.adapt = n.adapt)
update(my_jags, n.iter = n.update)
jags_res <- jags.samples(my_jags, 
                         variable.names = record_para, 
                         n.iter = n.sim)


compare = matrix(c(1,2), 1, 2)  # compare theta_i/theta_j
n.compare = nrow(compare)
source("source/analyze_jags_multi_layer.R")



# post-process group membership:
if (exists("H") == T & exists("true.g") == T & "categ" %in% record_para) {
  source("source/cluster_analysis.R")
}


# goodness-of-fit test:
if ("padj"  %in% record_para) {
  source("source/goodness_of_fit.R")
}


### calculate log likelihood
source("source/loglikelihood_LPJMM.R")

print(sessionInfo())





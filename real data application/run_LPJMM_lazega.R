model_name = "LPJMM"

n.adapt  = 20000
n.update = 20000
n.sim    = 10000
n.chains = 1
iter.keep= 1:n.sim
thin     = 10

H_test = 8

folder = glue::glue("output_cl={H_test}")
if (dir.exists(folder) == FALSE) {
  dir.create(folder)
}

# load data
load("data/lazega_post_data.rda")
isdirected = "directed"
x = x.pca

new.g  = split(1:N, list(all.group$practice, all.group$office))
true.g = rep(NA, N)
for (g in 1:length(new.g)){
  true.g[new.g[[g]]] = g
}
rm(new.g)

cl = mycol[true.g]
H = dplyr::n_distinct(true.g)

# jags setup:
source("source/specific_setup.R")

# generate Markov chain:
my_seed = 246
my_jags <- jags.model(model_path, data_prior, 
                      inits = list(".RNG.name" = "base::Marsaglia-Multicarry",
                                   ".RNG.seed" = my_seed),
                      n.chains = n.chains,
                      n.adapt = n.adapt)
update(my_jags, n.iter = n.update)
jags_res <- jags.samples(my_jags, 
                         variable.names = record_para, 
                         n.iter = n.sim)


compare = matrix(c(1, 2,
                   1, 3, 
                   2, 3), 
                 nrow = 3, ncol = 2, byrow = TRUE)  # compare theta_1/theta_2
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






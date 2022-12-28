library(tidyverse)
library(glue)
library(boot)
library(coda)
library(R2WinBUGS)
library(rjags)
library(magrittr)
library(glue)
library(MCMCpack)    ### to use procrustes


var_Z = 1

if (model_name == "LPJMM"){
  model_path = "jags_model/jags_Jin_cluster.txt"
  record_para = c("beta", "sigmasq", "tausq", "phi", 
                  "a", "b", "theta", "z", 
                  "padj", "categ")
  
  
  if (exists("H_test") == FALSE){
    H_test = 5
  }
  data_prior <- list(adj = Adj.layer,         ### adjacent matrix N by N
                     x   = as.vector(x),      ### the observed N attributes s(z_i)
                     sigmasq_prior1 = 2,      ### IG distribution.
                     sigmasq_prior2 = 1,     ### mean
                     tausq_prior1   = 2,      ### IG distribution
                     tausq_prior2   = 1,     ### mean
                     phi_lower      = 0.001,  ### uniform distribution
                     phi_upper      = 1,
                     amean_prior    = rep(0, L), ### length L
                     avar_prior     = 9,   
                     bmean_prior    = rep(0, L), ### length L
                     bvar_prior     = 9,   
                     theta_prior1   = 1,     ### mean = shape/rate.
                     theta_prior2   = 1,      ### Gamma (shape,rate).
                     alpha          = rep(1, H_test),
                     mu_prior1      = rep(0, K),
                     mu_prior2      = 2*var_Z/3,
                     kappa_prior1   = 3,      ### Inverse Gamma
                     kappa_prior2   = 2*var_Z/3,     ### mean of true kappasq
                     N = N,
                     K = K,         ### latent space dim
                     H = H_test,    ### assume N > H.
                     L = L,
                     I = diag(rep(1, N)) )
  tag_cluster  = H_test
  illustration = glue("The following results are based on
                       {model_name}'s model with {tag_cluster} clustering")
}



if (model_name == "SNSM"){
  model_path  = "jags_model/jags_Ciminelli.txt"
  
  record_para = c("beta", "sigmasq", "tausq", "phi", "a", "b", "z", 
                  "padj")
  data_prior <- list(adj            = Adj.layer,           ### adjacent matrix N by N
                     x              = as.vector(x),  ### the observed N attributes s(z_i)
                     sigmasq_prior1 = 2,
                     sigmasq_prior2 = 1,
                     tausq_prior1   = 2,
                     tausq_prior2   = 1,
                     phi_lower      = 0.001,
                     phi_upper      = 0.5,
                     amean_prior    = rep(0, L),      # length L
                     avar_prior     = 9,
                     bmean_prior    = rep(0, L),     # length L
                     bvar_prior     = 9,
                     mu_z           = rep(0, K),
                     var_z          = var_Z,
                     N = N,
                     K = K,         ### latent space dim, and H = 1
                     L = L,
                     I = diag(rep(1, N)) )
  illustration = glue("The following results are based on
                       {model_name}'s model")
}



print(glue("number of adaptation: 
            => {n.adapt} "));  

print(glue("number of burn-in: 
            => {n.update} "));  

print(glue("number of saved iterations: 
            => {n.sim} "));  

print(glue("jags model is stored in the file: 
            => {model_path} "));  


res_path  = glue("{folder}/{model_name}_cl={H_test}_var_Z={var_Z}.txt")





# calculate log likelihood
# SNSM

log.likeli.Ciminelli <- function(x, y, a, b, z){
  library("pracma")
  log_likeli = 0
  for(i in 1:(N-1)){
    for (j in (i+1):N){
      arg = a + b * abs(x[i] - x[j]) - pracma::Norm(z[i,]-z[j,])
      pij = ifelse(y[i,j] != 0, pnorm(arg), 1-pnorm(arg))
      log_pij = log(pij)
      log_likeli = log_likeli + log_pij
    }
  }
  return(log_likeli)
}


par(mfrow = c(n.chains, 1))
for (chain in 1:n.chains){
  log.likeli.chain = 0
  for (l in 1:L){
    a.chain = jags_res$a[l,,chain]        ### length(a.chain) = n.sim
    b.chain = jags_res$b[l,,chain]        ### length(b.chain) = n.sim
    z.chain = jags_res$z[,,,chain]       ### dim(z.chain) = [N, K, n.sim]
    
    log.likeli.layer = sapply(seq(1, n.sim, thin), 
                              function(n){log.likeli.Ciminelli(x, Adj.layer[,,l], a.chain[n], b.chain[n], 
                                                               z.chain[,,n])})
  }
  log.likeli.chain = log.likeli.chain + log.likeli.layer
  
  plot(log.likeli.chain, type = "l",
       ylab = "log likelihood",
       xlab = "iteration",
       main = glue('log likelihood for chain {chain}'))
}


save(log.likeli.chain,
     file = glue::glue("{folder}/loglikeli_chain_thin={thin}.rda"))

### calculate log likelihood


library(pracma)

log.likeli <- function(x = x, 
                       y = Adj.layer, 
                       a = jags_res$a[, iter, chain], 
                       b = jags_res$b[, iter, chain], 
                       theta = jags_res$theta[, iter, chain],  
                       z = jags_res$z[,, iter, chain]){
  log_likeli = 0
  
  for (l in 1:L) {
    for (i in 1:N) {
      for (j in setdiff(1:N, i)) {
        arg = a[l] + b[l] * abs(x[i] - x[j]) - theta[l] * norm(z[i,]-z[j,], type = "2")
        # length(arg) = 1
        pij = ifelse(Adj.layer[i,j,l]!= 0, pnorm(arg), 1-pnorm(arg))
        log_pij = log(pij)
        log_likeli = log_likeli + log_pij
      }
    }
  }
  
  return(log_likeli)
}



log.likeli.chain = c()

for (chain in 1:n.chains){
  ll.chain.i = sapply(seq(1, n.sim, thin), 
                      function(iter){log.likeli(x = x, 
                                                y = Adj.layer,
                                                a = jags_res$a[, iter, chain], 
                                                b = jags_res$b[, iter, chain], 
                                                theta = jags_res$theta[, iter, chain],  
                                                z = jags_res$z[,, iter, chain])})
  log.likeli.chain = rbind(log.likeli.chain, ll.chain.i)
  plot(ll.chain.i, type = "l",
       ylab = "log likelihood",
       xlab = "iteration",
       main = glue::glue('log likelihood for chain {chain} thin = {thin}'))
}

row.names(log.likeli.chain) = paste("chain ", 1:n.chains)

save(log.likeli.chain,
     file = glue::glue("{folder}/loglikeli_chain_thin={thin}.rda"))




capture.output(glue("seed = {my_seed}"),
               cat("----------------------------------------------------------","\n"),
               glue('n.update = {n.update}'),
               glue('n.sim = {n.sim}'),
               cat("----------------------------------------------------------","\n"),
               file = res_path, append = FALSE, split = TRUE)


record_para = names(jags_res)   
draw_para = c("beta", "phi", "sigmasq","tausq")

# CI table:
sum_para <- NULL

if (rlang::is_empty(draw_para) == FALSE){
  for (i in seq_along(draw_para)) {
    k = which(record_para == draw_para[i])
    para.chain = jags_res[[k]]   ### [ 1, n.sim, n.chains]
    
    for (j in 1:n.chains){
      single_chain <- c(para_name = noquote(draw_para[i]),
                        mean      = round(mean(para.chain[,,j]), 4),
                        round(quantile(para.chain[,,j], probs = 0.025), 4),
                        round(quantile(para.chain[,,j], probs = 0.975), 4))
      sum_para <- rbind(sum_para, chain = single_chain)
    }
  }
}

sum_a <- NULL
for(l in 1:L) {
  for (j in 1:n.chains){
    a_chain <- jags_res$a[l,,j]
    sum_a_layer <- c(para_name = glue::glue("a_{l}"),
                     mean = round(mean(a_chain),4),
                     round(quantile(a_chain, probs = 0.025),4),
                     round(quantile(a_chain, probs = 0.975),4))
    sum_a <- rbind(sum_a, sum_a_layer)
  }
}

sum_b <- NULL
for(l in 1:L) {
  for (j in 1:n.chains){
    b_chain <- jags_res$b[l,,j]
    sum_b_layer <- c(para_name = glue::glue("b_{l}"),
                     mean = round(mean(b_chain), 4),
                     round(quantile(b_chain, probs = 0.025), 4),
                     round(quantile(b_chain, probs = 0.975), 4))
    sum_b <- rbind(sum_b, sum_b_layer)
  }
}

sum_theta <- NULL
for(l in 1:L) {
  for (j in 1:n.chains){
    theta_chain <- jags_res$theta[l,,j]
    sum_theta_layer <- c(para_name = glue::glue("theta_{l}"),
                         mean = round(mean(theta_chain),4),
                         round(quantile(theta_chain, probs = 0.025),4),
                         round(quantile(theta_chain, probs = 0.975),4))
    sum_theta <- rbind(sum_theta, sum_theta_layer)
  }
}




# ratio of theta_i / theta_j
sum_ratio <- NULL
for (k in 1:n.compare) {
  for (j in 1:n.chains){
    theta_1 = jags_res$theta[compare[k,1],,j]
    theta_2 = jags_res$theta[compare[k,2],,j]
    ratio = theta_1/theta_2
    sum_ratio_layer <- c(para_name = glue::glue("ratio_{k}"),
                         mean = round(mean(ratio),4),
                         round(quantile(ratio, probs = 0.025),4),
                         round(quantile(ratio, probs = 0.975),4))
    sum_ratio = rbind(sum_ratio, sum_ratio_layer)
  }
}


# proportion of sigma^2 / sigma^2 + tau^2 
sum_propo <- NULL
for (j in 1:n.chains) {
  tausq = jags_res$tausq[1,,j]
  sigmasq = jags_res$sigmasq[1,,j]
  propo = sigmasq / (sigmasq + tausq)
  sum_propo_layer <- c(para_name = glue::glue("propo"),
                       mean = round(mean(propo),4),
                       round(quantile(propo, probs = 0.025),4),
                       round(quantile(propo, probs = 0.975),4))
  sum_propo = rbind(sum_propo, sum_propo_layer)
}



# summary of CI table:
sum_para = rbind(sum_para, sum_a, sum_b, sum_theta, sum_ratio, sum_propo)
para_table <- as.table(sum_para)
rownames(para_table) <- rep(paste0("chain ", 1:n.chains), times = nrow(sum_para)/n.chains)
para_table


capture.output(glue('Analyze posterior samples:'),
               para_table,
               cat("----------------------------------------------------------","\n"),
               file = res_path, append = TRUE, split = TRUE)






# draw and calculate latent positions Z:
mycol = c("deepskyblue", "orange", "lavenderblush3", "forestgreen",
          "coral", "cyan3", "darkblue", "darksalmon",
          "darkkhaki", "violetred", "lemonchiffon3",
          "darkmagenta")

par(mfrow = c(1,1))

z.mcmc = jags_res$z   ### [N, K, n.sim, n.chains]
z.post.mean.chain = array(NA, dim = c(N, K, n.chains))

for (i in 1:n.chains){
  if (!exists("Z")){
    Z.ref = z.mcmc[,,n.sim,i]
    note = "use last iteration as reference"
  } else {
    Z.ref = Z
    note = "use true Z as reference"
  }
  z.chain_i   <- z.mcmc[,,,i]
  z.transform <- array(NA, dim = dim(z.chain_i))  ### [N, K, n.sim]
  for(j in 1:n.sim){   
    z.trans <- MCMCpack::procrustes(z.chain_i[,,j], 
                                    Z.ref, # not Z because Z is unavailable.
                                    translation = T,
                                    dilation = T)$X.new
    z.transform[, , j] <- z.trans
  }
  z.post.mean = apply(z.transform, 1:2, mean)
  z.post.mean.chain[,,i] = z.post.mean
  
  if ( exists("true.g") == T){
    cl = mycol[true.g]
    plot(z.post.mean, xlab = "Z_1", ylab = "Z_2", pch = 19, col = cl, 
         main = glue("estimated latent positions (chain {i})
                      color based on true grouping"))
    text(z.post.mean, as.character(1:N), cex=0.6, pos=4, col="gray")
  } else if ("categ" %in% record_para) {
    cl.res = jags_res$categ[, iter.keep, ]  # dim: [N, n.iter, n.chains]
    cl.sample = cl.res[,,i]                 # i-th chain
    cl.sample = t(cl.sample)   
    my.psm <- mcclust::comp.psm(cl.sample)
    est.g = mcclust::maxpear(my.psm)$cl
    method = "maxpear"
    plot(z.post.mean, xlab = "Z_1", ylab = "Z_2", pch = 19, col = mycol[est.g], 
         main = glue("estimated latent positions (chain {i}) 
                      color based on {method} using MCMC"))
    if (exists("vertex.names") == TRUE){
      text(z.post.mean, vertex.names, cex=0.6, pos=4, col="gray")
    } else {
      text(z.post.mean, as.character(1:N), cex=0.6, pos=4, col="gray")
    }
  } else {
    plot(z.post.mean, xlab = "Z_1", ylab = "Z_2", pch = 19, col = "orange", 
         main = glue("estimated latent positions (chain {i})"))
    text(z.post.mean, as.character(1:N), cex=0.6, pos=4, col="gray")
  }
}



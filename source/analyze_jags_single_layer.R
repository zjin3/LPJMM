
capture.output(glue("seed = {my_seed}"),
               cat("----------------------------------------------------------","\n"),
               glue('n.update = {n.update}'),
               glue('n.sim = {n.sim}'),
               cat("----------------------------------------------------------","\n"),
               file = res_path, append = FALSE, split = TRUE)



# find 1-dim parameters:
record_para = names(jags_res)
rmv_para <- NULL
for (i in seq_along(record_para)){
  if (length(jags_res[[i]]) != n.chains * n.sim){
    rmv_para <- c(rmv_para, record_para[i])
  }
}
draw_para <- setdiff(record_para, rmv_para)

# CI table:
if (rlang::is_empty(draw_para) == FALSE){
  sum_para <- NULL
  for (i in seq_along(draw_para)) {
    k = which(record_para == draw_para[i])
    para.chain = jags_res[[k]]   ### [ 1, n.sim, n.chains]
    
    for (j in 1:n.chains){
      single_chain <- c(para_name = noquote(draw_para[i]),
                        mean      = mean(para.chain[,,j]),
                        quantile(para.chain[,,j], probs = 0.05),
                        quantile(para.chain[,,j], probs = 0.95))
      sum_para <- rbind(sum_para, chain = single_chain)
    }
  }
  para_table <- as.table(sum_para)
  rownames(para_table) <- rep(paste0("chain ", 1:n.chains), times = length(draw_para))
  para_table
  
  
  capture.output(glue('Analyze posterior samples:'),
                 para_table,
                 cat("----------------------------------------------------------","\n"),
                 file = res_path, append = TRUE, split = TRUE)
}


# plots of estimated Z:

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







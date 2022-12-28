# requires "H", "true.g" and that "categ" is recorded

g.mcmc.chain = jags_res$categ[,iter.keep,1]


# clustering analysis based on each iteration of g:

library(cluster)    
library(factoextra) 

analyze_g_chain <- function(g.chain){       
  mean_score = c()
  score.chain = apply(g.chain, MARGIN = 2, 
                      FUN = mclust::adjustedRandIndex,
                      y = true.g)
  # plot:
  plot(score.chain, type = "l",
       main = glue("Rand index chain based using MCMC result"))
  mean_score = c(mean_score, mean(score.chain))
  
  mean_score = t(as.matrix(mean_score))
  colnames(mean_score) = "Adj Rand index"
  row.names(mean_score)= "mean"
  return(mean_score)
}

avg_score_g_mcmc <- analyze_g_chain(g.mcmc.chain) 

capture.output(glue('Analyze of clustering based on MCMC:'),
               round(avg_score_g_mcmc, 4),
               cat("----------------------------------------------------------","\n"),
               file = res_path, append = TRUE, split = TRUE)




# clustering analysis based on clustering estimates

library(GreedyEPL)
library(mcclust)

cl_summary <- function(g.chain, seed = 1234){
  set.seed(seed)
  cl.mat = t(g.chain) # required each *row* to be a clustering
  
  est.g1 = GreedyEPL::MinimiseEPL(cl.mat)$decision
  
  my.psm <- mcclust::comp.psm(cl.mat)
  est.g2 = mcclust::minbinder(my.psm)$cl
  est.g3 = mcclust::maxpear(my.psm)$cl
  est.g4 = mcclust::medv(my.psm)
  
  clustering_summary = rbind(est.g1, est.g2, est.g3, est.g4)
  rownames(clustering_summary) = c("GreedyEPL", "minbinder", "maxpear", "medv")
  
  return(clustering_summary)
}


clustering_summary_mcmc = cl_summary(g.mcmc.chain)

save(clustering_summary_mcmc,
     z.post.mean.chain,
     file = glue::glue("{folder}/clustering_summary.rda"))


# plot estimated Z based on g estimates:
Z.temp = z.post.mean.chain[,,1]
par(mfrow = c(2,2))
for (name in row.names(clustering_summary_mcmc)){
  plot(Z.temp, xlab = "", ylab = "", pch = 19, 
       col = mycol[clustering_summary_mcmc[name,]], 
       main = glue::glue("{name}"))
  text(Z.temp, as.character(1:N), cex=0.6, pos=4, col="gray")
}
par(mfrow=c(1,1))


# calculate rand index:
analyze_g_estimators <- function(cluster.summary){ 
  res = c()
  library("mclust")
  scores = apply(cluster.summary, MARGIN = 1, 
                 FUN = mclust::adjustedRandIndex,
                 y = true.g)
  res = rbind(res, scores)
  rownames(res) = "Adjusted Rand index"
  return(res)
}

sum_score_g_mcmc = analyze_g_estimators(clustering_summary_mcmc)

capture.output(glue('Analyze of summary of clustering based on MCMC:'),
               round(sum_score_g_mcmc, 4),
               cat("----------------------------------------------------------","\n"),
               file = res_path, append = TRUE, split = TRUE)


capture.output(glue('number of clusters:'),
               apply(clustering_summary_mcmc, MARGIN = 1, FUN = n_distinct),
               cat("----------------------------------------------------------","\n"),
               file = res_path, append = TRUE, split = TRUE)





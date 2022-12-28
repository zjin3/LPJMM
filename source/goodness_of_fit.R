
library(network)
library(igraph)

stopifnot("padj"  %in% record_para)

# truth:
density.true = c()
for (l in 1:L){
  g = as.network(Adj.layer[,,l])
  density.true = c(density.true, network.density(g))
}
density.true

transit.true = c()
for (l in 1:L){
  g = graph_from_adjacency_matrix(Adj.layer[,,l], mode = "directed", diag = F)
  transit.true = c(transit.true, transitivity(g))
}
transit.true

if ( exists("true.g") == T ){
  assorta.true = c()
  for (l in 1:L){
    g = graph_from_adjacency_matrix(Adj.layer[,,l], mode = "directed", diag = F)
    assorta.true = c(assorta.true, assortativity_nominal(g, types = true.g))
  }
  assorta.true
}

# Goodness of fit:
set.seed(123)

density.array = c()  # along L
density.mean  = c()  # along L 

transit.array = c()
transit.mean  = c()

assorta.array = c()
assorta.mean  = c()


for (l in 1:L){
  Padj.keep = jags_res$padj[, , l, iter.keep, 1]
  density.chain = c()
  transit.chain = c()
  assorta.chain = c()
  
  for (t in seq_along(iter.keep)){
    Padj = Padj.keep[,,t]
    
    if (isdirected == "directed"){
      Adj.t = apply(Padj, c(1,2), function(x){sample(0:1, 1, prob = c(1-x, x))})
      diag(Adj.t) = rep(0, N)
    } else {
      Adj.t = matrix(NA, N, N)
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          Adj.t[i, j] = sample(0:1, 1, prob = c(1 - Padj[i, j], Padj[i, j]))
          Adj.t[j, i] = Adj.t[i, j]
        }
      }
      diag(Adj.t) = rep(0, N)
    }
    
    G.t = as.network(Adj.t)
    density.t = network.density(G.t)
    density.chain = c(density.chain, density.t)
    
    G2.t = graph_from_adjacency_matrix(Adj.t, mode = "directed", diag = F)
    tran.t = transitivity(G2.t)
    transit.chain = c(transit.chain, tran.t)
    
    if ( exists("true.g") == T ){
      asot.t = assortativity_nominal(G2.t, types = true.g)
      assorta.chain = c(assorta.chain, asot.t)
    }
  }
  
  density.array = rbind(density.array, density.chain)
  transit.array = rbind(transit.array, transit.chain)
  if ( exists("true.g") == T ){
    assorta.array = rbind(assorta.array, assorta.chain)
  }
  
  density.mean = c(density.mean, mean(density.chain, na.rm = T)) # collect for all layers
  transit.mean = c(transit.mean, mean(transit.chain, na.rm = T))
  if ( exists("true.g") == T ){
    assorta.mean = c(assorta.mean, mean(assorta.chain, na.rm = T)) 
  }

  par(mfrow = c(3,1))
  # plots
  plot(density(density.chain, na.rm = T), 
       main = glue("density chain for layer {l}
                    true density = {format(density.true[l], digits = 3)}"))
  abline(v = density.true[l])
  
  plot(density(transit.chain, na.rm = T), 
       main = glue("transitivity chain for layer {l}
                   true transitivity = {format(transit.true[l], digits = 3)}"))
  abline(v = transit.true[l])
  
  if ( exists("true.g") == T ){
    plot(density(assorta.chain, na.rm = T), 
         main = glue("assortativity chain for layer {l}
                   true assortativity = {format(assorta.true[l], digits = 3)}"))
    abline(v = assorta.true[l])
  }
}
par(mfrow = c(1,1))

density.mean = cbind(density.true, density.mean)
row.names(density.mean)  = paste0(glue("layer "), 1:L)
row.names(density.array) = paste0(glue("layer "), 1:L)


transit.mean = cbind(transit.true, transit.mean)
row.names(transit.mean)  = paste0(glue("layer "), 1:L)
row.names(transit.array) = paste0(glue("layer "), 1:L)


if ( exists("true.g") == T ){
  assorta.mean = cbind(assorta.true, assorta.mean)
  row.names(assorta.mean)  = paste0(glue("layer "), 1:L)
  row.names(assorta.array) = paste0(glue("layer "), 1:L)
}


capture.output(glue("Goodness of fit (mean density)"),
               density.mean,
               cat("----------------------------------------------------------","\n"),
               glue("Goodness of fit (mean transitivity)"),
               transit.mean,
               cat("----------------------------------------------------------","\n"),
               file = res_path, append = TRUE, split = TRUE)

if ( exists("true.g") == T ){
  capture.output(glue("Goodness of fit (mean assortativity)"),
                 assorta.mean,
                 cat("----------------------------------------------------------","\n"),
                 file = res_path, append = TRUE, split = TRUE)
}


if (exists("true.g") == T) {
  save(density.array,
       transit.array,
       assorta.array,
       file = glue::glue("{folder}/GoF_chain.rda"))
} else {
  save(density.array,
       transit.array,
       file = glue::glue("{folder}/GoF_chain.rda"))
}

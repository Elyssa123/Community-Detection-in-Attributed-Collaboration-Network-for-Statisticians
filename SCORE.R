#' SCORE 
#' We thank the code provided by
#' @article{ji2016coauthorship, title={Coauthorship and citation networks for statisticians}, author={Ji, Pengsheng and Jin, Jiashun}, journal={The Annals of Applied Statistics}, volume={10}, number={4}, pages={1779--1812}, year={2016}, publisher={Institute of Mathematical Statistics}}

library(igraph)
 
score = function(adj, K)
{
  # input: 
  #    adj - symmetric adjacency matrix of a connected and undirected network/graph
  #    K   - the number of communities
  # output:  the community lables, starting from 1, 2, ...

  
  # check connectivity
  if(!is.connected(graph.adjacency(adj))) stop("The network based on the adjacency matrix is not connected. Use the clusters function in the igraph package to extract the (maximal connected) giant component.")
  # make sure K > 1
  if(K < 2) stop("K must be bigger than 1.")
  
  
  # check symmetry for adj
  if(sum(adj != t(adj))>0) stop("The adjacency matrix is not symmetric.")
  
  x = eigen(adj, symmetric=T)$vectors
  r = x[,2:K]/x[,1]
  
  m = nrow(adj)
  threshold = log(m) 
  r[r > threshold] =  threshold
  r[r < -threshold] = -threshold
  label = kmeans(r, centers=K, nstart=20, iter=100)$cluster
  
  return(label)
}

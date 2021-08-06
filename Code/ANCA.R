#' ANCA
#' We thank the code provided by
#' @INPROCEEDINGS{10.1007/978-3-319-72150-7_20, author="Falih, Issam and Grozavu, Nistor and Kanawati, Rushed and Bennani, Youn{\`e}s", title="ANCA : Attributed Network Clustering Algorithm", booktitle="Complex Networks {\&} Their Applications VI", year="2018", publisher="Springer International Publishing", address="Cham", pages="241--252",isbn="978-3-319-72150-7"}
#' @INPROCEEDINGS{8489183, author={I. {Falih} and N. {Grozavu} and R. {Kanawati} and Y. {Bennani} and B. {Matei}}, booktitle={2018 International Joint Conference on Neural Networks (IJCNN)}, title={Collaborative Multi-View Attributed Networks Mining}, year={2018}, pages={1-8}, doi={10.1109/IJCNN.2018.8489183}, ISSN={2161-4407}, month={July}}

library(rsvd)

select_rank <- function(g, k, num1, num2, discrete_attr, continuous_attr, contextView=TRUE){
  
  # input: 
  #    g                -  graph 
  #    K                -  the number of communities
  #    num1             -  the set of numbers of the large eigenvectors of topology similarity matrix 
  #    num2             -  the set of numbers of the large eigenvectors of attribute similarity matrix 
  #    discrete_attr    -  the names of nodes' discrete attributes
  #    continuous_attr  -  the names of nodes' continuous attributes
  # output: the best num1 and num2
  
  topoMatrix = compute_topoMatrix(g)  # compute topology similarity matrix 
  attrMatrix = compute_attrMatrix(g, discrete_attr = discrete_attr, continuous_attr = continuous_attr)  # compute attribute similarity matrix
  
  
  modu_result <- array(NA, dim = c(length(num1),length(num2)))
  dens_result <- array(NA, dim = c(length(num1),length(num2)))
  cond_result <- array(NA, dim = c(length(num1),length(num2)))
  inte_result <- array(NA, dim = c(length(num1),length(num2)))
  
  # conduct ANCA by using each number in num1 and num2
  for (i in 1:length(num1)) {
    for (j in 1:length(num2)) {
      set.seed(400)  
      result1 <- ANCA_discrete_continuous(g, nbCluster=k, 
                                          nbEigenVectorSeeds=num1[i], nbEigenVectorView=num2[j],
                                          topoMatrix=topoMatrix, attrMatrix=attrMatrix,
                                          discrete_attr=discrete_attr, continuous_attr=continuous_attr,
                                          contextView=contextView) 
      comp_inde <- summary.partition.result(g, result1)  # summary the partition result
      
      modu_result[i,j] <- comp_inde['Modularity']
      dens_result[i,j] <- comp_inde['Density']
      cond_result[i,j] <- comp_inde['Conductance']
      inte_result[i,j] <- internal.density.clusters(g, result1)
    }
  }
  
  all_result <- modu_result+dens_result+cond_result+inte_result
  all_max <- which(all_result==max(all_result), arr.ind=TRUE)   # find out where the maximum is 
  num1_best <- num1[all_max[1,'row']]
  num2_best <- num2[all_max[1,'col']]

  result = c(num1_best,num2_best)
  names(result) = c('num1_best','num2_best')
  return(result)
}


compute_topoMatrix <- function(graph, seedsFunction = DetectSeeds_Centrality, 
                                verbose=FALSE, structuralSimilarityFunction = distances){
  
  # input: 
  #    graph                         -  graph 
  #    seedsFunction                 -  the function to select seeds
  #    structuralSimilarityFunction  -  the function to measure structural similarity 
  # output: topology similarity matrix  
  
  Seeds  <- seedsFunction(graph)                         # select seeds
  if(verbose)
    cat(paste(length(Seeds)," seeds founded ",sep=""),sep = "\n")   
  if(length(Seeds)==0)
    stop("No seeds-nodes found ")
  if(verbose)
    cat("Step 2 : Characterize each node with the set of seeds using measure",sep = "\n")
  
  # compute the topology similarity matrix  (the length of shortest path from each node to the seed node)
  topoMatrix =  structuralSimilarityFunction(graph,to=Seeds)
  
  # If the network is not connected, then replace Inf in the topology similarity matrix  with the maximum value+1 in the topology similarity matrix 
  if(!is.connected(graph)){
    topoMatrix[topoMatrix==Inf] = max(topoMatrix[is.finite(topoMatrix)]) + 1
  }
  if(verbose)
    cat("Step 3 : Matrix factorization on topological distance",sep = "\n")
  return(topoMatrix)
}




compute_attrMatrix <- function(graph, discrete_attr, continuous_attr,contextSimilarityFunction_discrete=similarity.matchingCoefficient, 
                                contextSimilarityFunction_continuous=similarity.euclidean2){
  
  # input: 
  #    graph                                 -  graph 
  #    discrete_attr                         -  the names of nodes' discrete attributes
  #    continuous_attr                       -  the names of nodes' continuous attributes
  #    contextSimilarityFunction_discrete    -  the function to compute the similarity of discrete attributes
  #    contextSimilarityFunction_continuous  -  the function to compute the similarity of continuous attributes
  # output: attribute similarity matrix 
  
  # compute discrete attribute similarity matrix
  attrMatrix1 =  contextSimilarityFunction_discrete(graph,discrete_attr,normalization=FALSE)
  # compute continuous attribute similarity matrix
  attrMatrix2 =  contextSimilarityFunction_continuous(graph,continuous_attr)
  # sum up
  attrMatrix = attrMatrix1 + attrMatrix2
  return(attrMatrix)
}




ANCA_discrete_continuous <- function(graph, nbCluster=100, nbEigenVectorSeeds=2, nbEigenVectorView=2,
                                      topoMatrix, attrMatrix, contextView=TRUE, structuralView=TRUE,
                                      iter.max=100, normalization=TRUE, verbose=FALSE,
                                      discrete_attr, continuous_attr){
  
  # input: 
  #    graph               -  graph 
  #    nbCluster           -  the number of the communities 
  #    nbEigenVectorSeeds  -  the number of the large eigenvectors of topology similarity matrix 
  #    nbEigenVectorView   -  the number of the large eigenvectors of attribute similarity matrix 
  #    topoMatrix          -   the topology similarity matrix
  #    attrMatrix          -   the attribute similarity matrix
  #    discrete_attr       -  the names of nodes' discrete attributes
  #    continuous_attr     -  the names of nodes' continuous attributes
  # output: community detection result  
  
  if(!is.igraph(graph))
    stop("Graph should be an igraph graph")
  if(!is.connected(graph)&& verbose)
    cat("Graph is not connected",sep = "\n")
  if(contextView && length(vertex_attr_names(graph))<=1)   
    stop("Graph must have a vertex attribute")
  if(!("name" %in% vertex_attr_names(graph)))              
    stop("Graph must have a name vertex attribute")
  n = vcount(graph)                                        
  verticesName = V(graph)$name                             
  if(structuralView){
    topoMatrix = topoMatrix
    UTopo <- rsvd(topoMatrix, k=nbEigenVectorSeeds)$u     # extract the top k large eigenvectors of the topology similarity matrix
  }
  
  if(contextView){
    attrMatrix = attrMatrix
    UAttr <- rsvd(attrMatrix,k = nbEigenVectorView)$u     # extract the top k large eigenvectors of the attribute similarity matrix
    U <- cbind(UTopo,UAttr)     
  }
  else
    U = UTopo
  
  # Normalisation
  if( verbose)
    cat("Step 5 : Normalization",sep = "\n")
  if(normalization)
    U <- t(t(U)/sqrt(colSums(U^2)))     
  
  # Cluster each row via K-means
  if( verbose)
    cat("Step 6 : Cluster each row",sep = "\n")
  
  if(length(nbCluster)==1){
    vect <- kmeans(U,nbCluster,iter.max = iter.max)$cluster
    names(vect) <- as.vector(get.vertex.attribute(graph,name = "name"))
    return(vect)
  }
  else{
    partition = list()
    for(k in nbCluster){
      vect <- kmeans(U,k,iter.max = iter.max)$cluster
      names(vect) <- as.vector(get.vertex.attribute(graph,name = "name"))
      partition = c(partition, list(vect))
    }
    return(partition)
  }
  
}


summary.partition.result <- function(graph,partition){
  
  # input: 
  #    graph      -  graph 
  #    partition  -  the result of community detection
  # output: evaluation of community detection: modularity, entropy, density, conductance, size, avgsize
  
  if(!is.igraph(graph))
    stop("Graph should be an igraph graph")
  
  result = vector()
  # modularity
  if(is.list(partition))
    result= c(result,round(modularity(graph,groups2memb(partition)),4))
  else
    result= c(result,round(modularity(graph,partition),4))     
  # entropy
  result= c(result,round(mean(entropy.attribute(graph,partition)),4))
  # density
  result= c(result,round(density.clusters(graph,partition),4))
  # conductance
  result= c(result,round(conductance(graph,partition),4))
  
  if(!is.list(partition)){
    result= c(result,length(memb2groups(partition)))
    result= c(result,round(mean(sapply(memb2groups(partition),length))))
  }
  else{
    result= c(result,length(partition))
    result= c(result,round(mean(sapply(partition,length))))
  }
  
  names(result) = c("Modularity","Entropy","Density","Conductance","size","avgsize")
  return(result)
}




DetectSeeds_Centrality <- function(graph,topPct=0.15, downPct=0.05){
  
  # input: 
  #    graph    -  graph 
  #    topPct   -  the percentage of the high centrality 
  #    downPct  - the percentage of the low centrality 
  # output: seeds
  
  if(!is.igraph(graph))
    stop("not an igraph graph")
  Seeds <- vector()
  # extract nodes with high/low page-rank centrality, eigenvector centrality and degree centrality respectively 
  for(centrality in c("centralization.evcent(graph)$vector","page_rank(graph)$vector", "centralization.degree(graph)$res")){
    nodesCentralityList <-  eval(parse(text=centrality))
    names(nodesCentralityList) <- V(graph)$name
    Seeds <- c(Seeds, names(sort(nodesCentralityList,decreasing = TRUE)[1:round(topPct*vcount(graph))]))
    Seeds <- c(Seeds,names(sort(nodesCentralityList,decreasing = FALSE)[1:round(downPct*vcount(graph))]))
  }
  return(unique(Seeds))
}



similarity.matchingCoefficient <- function(graph, discrete_attr, normalization=TRUE){
  
  # input: 
  #    graph           -  graph 
  #    discrete_attr   -  the names of nodes' discrete attributes
  # output: the attribute similarity matrix of discrete attributes
  
  attr_all = getAll.attribute(graph)  # get all nodes attribute of a given graph
  attr = attr_all[discrete_attr]
  if(normalization)
    return(do.call(rbind,lapply(V(graph)$name, function(i){rowSums(sapply(lapply(seq(ncol(attr)), function(j){ attr[,j]==attr[i,j] }),rbind))}))/ncol(attr))
  else
    return(do.call(rbind,lapply(V(graph)$name, function(i){rowSums(sapply(lapply(seq(ncol(attr)), function(j){ attr[,j]==attr[i,j] }),rbind))})))
}


similarity.euclidean2 <- function(graph, continuous_attr){
  
  # input: 
  #    graph           -  graph 
  #    continuous_attr   -  the names of nodes' continuous attributes
  # output: the attribute similarity matrix of continuous attributes
  
  attr_all = getAll.attribute(graph)  # get all nodes attribute of a given graph
  attr = attr_all[continuous_attr]
  attr = as.data.frame(lapply(attr,normalize))
  sim_matr <- matrix(0,ncol = nrow(attr),nrow = nrow(attr))
  sim_matr[lower.tri(sim_matr)] = dist(attr,method = 'manhattan')
  sim_matr <- sim_matr + t(sim_matr)
  sim_matr <- length(continuous_attr)-sim_matr
  return(sim_matr)
}

# min-max normalization
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}



memb2groups <- function(vect){
  
  # input: 
  #    vect  -  the membership vector of a set of object
  # output: a list of numeric vector each one represent a community
  
  L <- list()
  if(is.null(names(vect))){
    for(i in unique(vect))
      L <- c(L, list(which(vect==i)))
  }
  else{
    for(i in unique(vect))
      L <- c(L,list(names(which(vect==i))))
  }
  return(L)
  
}


#' get all attribute
#'
#' @description return all vertices attribute of a given graph
#' @usage getAll.attribute(graph,binarization, rowNames)
#' @param graph : an igraph graph
#' @param binarization : logical
#' @param rowNames : vector of nodes name
#' @return returns a dataframe containing for each vertex (in rows) all the attribute value ( in columns)
#' @author  Issam Falih <issam.falih@lipn.univ-paris13.fr>
#' @examples
#' graph = DBLP10K
#' getAll.attribute(graph)
#' @export



getAll.attribute <- function(graph, rowNames=NULL){
  
  # input: 
  #    graph     -  graph
  #    rowNames  -  vector of nodes name
  # output: a dataframe containing for each vertex (in rows) all the attribute value (in columns)
  
  if(!is.igraph(graph))
    stop("not an igraph graph")
  if(is.null(vertex_attr(graph)$name))
    stop("graph must have a name vertex attribute ")
  if(is.null(rowNames))
    rowNames = V(graph)$name
  attribute_df <- as.data.frame(sapply(setdiff(vertex_attr_names(graph),c("name" ,"id")),function(attrName){get.vertex.attribute(graph,name = attrName)}))
  row.names(attribute_df) <- rowNames
  return(attribute_df)
}
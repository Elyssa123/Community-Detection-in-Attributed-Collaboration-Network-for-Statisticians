#' Comparison metrics
#' We thank the code provided by
#' @INPROCEEDINGS{10.1007/978-3-319-72150-7_20, author="Falih, Issam and Grozavu, Nistor and Kanawati, Rushed and Bennani, Youn{\`e}s", title="ANCA : Attributed Network Clustering Algorithm", booktitle="Complex Networks {\&} Their Applications VI", year="2018", publisher="Springer International Publishing", address="Cham", pages="241--252",isbn="978-3-319-72150-7"}
#' @INPROCEEDINGS{8489183, author={I. {Falih} and N. {Grozavu} and R. {Kanawati} and Y. {Bennani} and B. {Matei}}, booktitle={2018 International Joint Conference on Neural Networks (IJCNN)}, title={Collaborative Multi-View Attributed Networks Mining}, year={2018}, pages={1-8}, doi={10.1109/IJCNN.2018.8489183}, ISSN={2161-4407}, month={July}}


conductance <- function(graph, partition){

  # input: 
  #    graph            -  graph 
  #    partition        -  the partition of nodes
  # output: the conductance of the partition
  
  if(!is.igraph(graph))
    stop("Should be applied on a igraph graph object")
  cond = vector()
  if(is.list(partition)){
    adj = get.adjacency(graph,type = "both",names = TRUE)
    verticesNames = colnames(adj)
    for(com in partition)
      cond = c(cond, sum(adj[com,verticesNames[!(verticesNames %in% com)]])/min( sum(adj[com,]),sum(adj[verticesNames[!(verticesNames %in% com)],])))
    return(1 - sum(cond,na.rm = TRUE)/length(partition))
  }
  else if(is.vector(partition)){
    conductance(graph,memb2groups(partition))
  }
  else
    stop("argument partition should be a list or a vector")
}




entropy.attribute <- function(graph, partition,categorical=TRUE){
  
  # input: 
  #    graph            -  graph 
  #    partition        -  the partition of nodes
  # output: the entropy of the partition
  
  if(!is.igraph(graph))
    stop("Should be applied on a igraph graph object")
  if(is.list(partition)){
    entrop <- vector()
    if(categorical){
      for(attribute_name in setdiff(vertex_attr_names(graph),c("name","id"))){
        attribute_domaine = unique(get.vertex.attribute(graph,attribute_name))
        entropy = 0
        
        for(cluster in partition){
          subg <- induced_subgraph(graph,cluster)
          pct = vcount(subg)/vcount(graph)
          for(domaine in attribute_domaine){
            p = (length(which(vertex_attr(subg, attribute_name)==domaine))/vcount(subg))
            if(p>0){
              entropy = entropy + pct *  (p * log2(p))
            }
          }
        }
        entrop <- c(entrop, round(- entropy,3))
        #print(paste("Entropy for the attribute : ",attribute_name," = ", round(- entropy,2),sep=""))
      }
      #print(paste("The general entropy is  = ", round(mean(entrop),2),sep=""))
    }
    else{
      for(attribute_name in setdiff(vertex_attr_names(graph),c("name","id"))){
        attribute_domaine = unique(get.vertex.attribute(graph,attribute_name))
        entropy = 0
        
        for(cluster in partition){
          subg <- induced_subgraph(graph,cluster)
          pct = vcount(subg)/vcount(graph)
          p = mean(vertex_attr(subg, attribute_name))
          if(p>0)
            entropy = entropy + pct *  (p * log2(p))
        }
        entrop <- c(entrop, round(- entropy,3))
      }
    }
    names(entrop) <- setdiff(vertex_attr_names(graph),c("id","name"))
    return(entrop)
  }
  else if(is.vector(partition)){
    entropy.attribute(graph, memb2groups(partition))
  }
  else
    stop("argument should be a list or vector ")
}



density.clusters <- function(graph, partition){
  
  # input: 
  #    graph            -  graph 
  #    partition        -  the partition of nodes
  # output: the density of the partition
  
  if(!is.igraph(graph))
    stop("Should be applied on a igraph graph object")
  d = 0
  if(is.list(partition)){
    for(cluster in partition){
      d = d + ecount(induced_subgraph(graph,cluster))
    }
    return(d/ecount(graph))
  }
  else if(is.vector(partition)){
    density.clusters(graph,memb2groups(partition))
  }
  else
    stop("argument partition should be a list or a vector")
}





internal.density.clusters <- function(graph, partition){
  
  # input: 
  #    graph            -  graph 
  #    partition        -  the partition of nodes
  # output: the internal density of the partition
  
  if(!is.igraph(graph))
    stop("Should be applied on a igraph graph object")
  d = 0
  if(is.list(partition)){
    for(cluster in partition){
      if(length(cluster)<=1)
        next
      m =   ecount(induced_subgraph(graph,cluster))
      n = length(cluster)
      d = d + ((2 * m )/(n * (n - 1)))
    }
    return(d/length(partition))
  }
  else if(is.vector(partition)){
    density.clusters(graph,memb2groups(partition))
  }
  else
    stop("argument partition should be a list or a vector")
}


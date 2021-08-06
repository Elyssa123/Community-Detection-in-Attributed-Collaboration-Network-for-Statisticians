###################################################################
#
# This is the main code to produce the results in the paper 
# "Community Detection in Attributed Collaboration Network for Statisticians".
#
####################################################################


library(foreach)
library(igraph)
library(CatEncoders)  # one-hot encoder

source('ANCA.R')
source('ECV.R')
source('comparison_metrics.R')
source('SCORE.R')

load('coauthorship(43journal)_core3.RData')


# 连通分量
component <- clusters(g_core)
component$no  # the number of components

# 提取最大连通分量
g_core3_comp1 <- induced_subgraph(g_core, vids = V(g_core)$name[component$membership==1])
vcount(g_core3_comp1)
W_core3_comp1 <- as_adjacency_matrix(g_core3_comp1)


####################################################
#
#   Estimating the number of communities by ECV
#
####################################################

registerDoSEQ()
result <- foreach(k = 1:200) %dopar% {
  set.seed(300+k)
  system.time(random.est <- ECV.undirected.Rank.weighted(W_core3_comp1,50,B=3,holdout.p=0.1,soft=FALSE,fast=TRUE))
  # tmp <- which.min(random.est$sse)
  tmp <- random.est$sse
}
result
result <- matrix(unlist(result),ncol=50,byrow=TRUE)
result_k <- apply(result, 1, function(x){which.min(x)}) 
mean(result_k)



####################################################
#
#   Community detection by ANCA
#
####################################################

num1 <- c(39,40,41) # best 40
num2 <- c(5,6,7) # best 6

ANCA_result <- select_rank(g_core3_comp1,48, num1, num2,discrete_attr=c('prolific','region_common'),continuous_attr=vertex_attr_names(g_core3_comp1)[3:32])
ANCA_result  # find the best num1 and num2

topoMatrix = compute_topoMatrix(g_core3_comp1)  # compute the topology similarity matrix of the graph
attrMatrix = compute_attrMatrix(g_core3_comp1,discrete_attr=c('prolific','region_common'),continuous_attr=vertex_attr_names(g_core3_comp1)[3:32])  # compute the attribute similarity matrix of the graph
set.seed(400)

# conduct ANCA algorithm
ANCA_com <- ANCA_discrete_continuous(g_core3_comp1, nbCluster=48, topoMatrix=topoMatrix, attrMatrix=attrMatrix,nbEigenVectorSeeds=ANCA_result['num1_best'], 
                                    nbEigenVectorView=ANCA_result['num2_best'],discrete_attr=c('prolific','region_common'),continuous_attr=vertex_attr_names(g_core3_comp1)[3:32])
summary.partition.result(g_core3_comp1, ANCA_com)  # summary the partition result of ANCA



################################################################
#
#   Comparison of ANCA, SCORE, leading eigenvector and k-means
#
################################################################

#---------------------------------------Leading eigenvector----------------------------------
LE_cluster_result <- cluster_leading_eigen(g_core3_comp1)
names(LE_cluster_result$membership) <- as.vector(get.vertex.attribute(g_core3_comp1, name = "name"))
LE_result1 <- summary.partition.result(g_core3_comp1, LE_cluster_result$membership)   # summary the partition result of ANCA
LE_result2 <- internal.density.clusters(g_core3_comp1, LE_cluster_result$membership)  # compute the internal density of the partition

LE_result <- c(LE_result1[c(1,3,4)],LE_result2)
names(LE_result) = c("Modularity","Density","Conductance","Internal Density")
LE_result


LE_comm_num = 37

#-------------------------------------------ANCA---------------------------------------------
num1 <- c(8,9)
num2 <- c(1,2,3)
ANCA_result_10 <- select_rank(g_core3_comp1,10, num1, num2, discrete_attr=c('prolific','region_common'),continuous_attr=vertex_attr_names(g_core3_comp1)[3:32])

num1 <- c(9,13,20)
num2 <- c(1,3,4)
ANCA_result_20 <- select_rank(g_core3_comp1,20, num1, num2,discrete_attr=c('prolific','region_common'),continuous_attr=vertex_attr_names(g_core3_comp1)[3:32])

num1 <- c(28,29,30,31,32,33)
num2 <- c(1,2)
ANCA_result_30 <- select_rank(g_core3_comp1,30, num1, num2,discrete_attr=c('prolific','region_common'),continuous_attr=vertex_attr_names(g_core3_comp1)[3:32])

num1 <- c(38,39)
num2 <- c(14,15,16)
ANCA_result_40 <- select_rank(g_core3_comp1,40, num1, num2,discrete_attr=c('prolific','region_common'),continuous_attr=vertex_attr_names(g_core3_comp1)[3:32])

num1 <- c(41,42)
num2 <- c(19,20,21)
ANCA_result_50 <- select_rank(g_core3_comp1,50, num1, num2,discrete_attr=c('prolific','region_common'),continuous_attr=vertex_attr_names(g_core3_comp1)[3:32])

num1 <- c(46,48,49)
num2 <- c(18,19,20)
ANCA_result_LE_comm_num <- select_rank(g_core3_comp1,LE_comm_num, num1, num2,discrete_attr=c('prolific','region_common'),continuous_attr=vertex_attr_names(g_core3_comp1)[3:32])

ANCA_result_df <- data.frame(community_num=c(10,20,30,40,50,LE_comm_num),Modularity=NA,Density=NA,Conductance=NA,Internal_Density=NA)

ANCA_result_df[1,c(2:5)] <- ANCA_result_10[3:6]
ANCA_result_df[2,c(2:5)] <- ANCA_result_20[3:6]
ANCA_result_df[3,c(2:5)] <- ANCA_result_30[3:6]
ANCA_result_df[4,c(2:5)] <- ANCA_result_40[3:6]
ANCA_result_df[5,c(2:5)] <- ANCA_result_50[3:6]
ANCA_result_df[6,c(2:5)] <- ANCA_result_LE_comm_num[3:6]
ANCA_result_df


#-------------------------------------------SCORE-------------------------------------------
measure_SCORE_result <- function(g,k){
  W <- as_adjacency_matrix(g)
  cluster_result <- score(W,k)
  names(cluster_result) <- as.vector(get.vertex.attribute(g, name = "name"))
  SCORE_result1 <- summary.partition.result(g, cluster_result)
  SCORE_result2 <- internal.density.clusters(g, cluster_result)
  SCORE_result <- c(SCORE_result1[c(1,3,4)],SCORE_result2)
  names(SCORE_result) = c("Modularity","Density","Conductance","Internal Density")
  return(SCORE_result)
}
set.seed(400)
SCORE_result_10 <- measure_SCORE_result(g_core3_comp1,10)
set.seed(400)
SCORE_result_20 <- measure_SCORE_result(g_core3_comp1,20)
set.seed(400)
SCORE_result_30 <- measure_SCORE_result(g_core3_comp1,30)
set.seed(400)
SCORE_result_40 <- measure_SCORE_result(g_core3_comp1,40)
set.seed(400)
SCORE_result_50 <- measure_SCORE_result(g_core3_comp1,50)
set.seed(400)
SCORE_result_LE_comm_num <- measure_SCORE_result(g_core3_comp1,LE_comm_num)

SCORE_result_df <- data.frame(community_num=c(10,20,30,40,50,LE_comm_num),Modularity=NA,Density=NA,Conductance=NA,Internal_Density=NA)

SCORE_result_df[1,c(2:5)] <- SCORE_result_10
SCORE_result_df[2,c(2:5)] <- SCORE_result_20
SCORE_result_df[3,c(2:5)] <- SCORE_result_30
SCORE_result_df[4,c(2:5)] <- SCORE_result_40
SCORE_result_df[5,c(2:5)] <- SCORE_result_50
SCORE_result_df[6,c(2:5)] <- SCORE_result_LE_comm_num
SCORE_result_df


#-------------------------------------------k-means-------------------------------------------
attr_df <- getAll.attribute(g_core3_comp1)
attr_df[1:31]=as.data.frame(lapply(attr_df[1:31],normalize))
attr_df$region_common <- as.character(attr_df$region_common)
total_onehot <- OneHotEncoder.fit(as.matrix(attr_df$region_common))  # one-hot encoder
region_onehot <- as.data.frame(transform(total_onehot,attr_df['region_common'],sparse=FALSE))

attr_df_new <- cbind(attr_df[,-32],region_onehot)
dim(attr_df_new)

measure_kmeans_result <- function(g, attributes_info, k){
  vect <- kmeans(attributes_info,k,iter.max = 100)$cluster
  kmeans_result1 <- summary.partition.result(g, vect)
  kmeans_result2 <-  internal.density.clusters(g, vect)
  kmeans_result <- c(kmeans_result1[c(1,3,4)],kmeans_result2)
  names(kmeans_result) = c("Modularity","Density","Conductance","Internal Density")
  return(kmeans_result)
}


set.seed(400)
kmeans_result_10 <- measure_kmeans_result(g_core3_comp1,attr_df_new,10)
set.seed(400)
kmeans_result_20 <- measure_kmeans_result(g_core3_comp1,attr_df_new,20)
set.seed(400)
kmeans_result_30 <- measure_kmeans_result(g_core3_comp1,attr_df_new,30)
set.seed(400)
kmeans_result_40 <- measure_kmeans_result(g_core3_comp1,attr_df_new,40)
set.seed(400)
kmeans_result_50 <- measure_kmeans_result(g_core3_comp1,attr_df_new,50)
set.seed(400)
kmeans_result_LE_comm_num <- measure_kmeans_result(g_core3_comp1,attr_df_new,LE_comm_num)

kmeans_result_df <- data.frame(community_num=c(10,20,30,40,50,LE_comm_num),Modularity=NA,Density=NA,Conductance=NA,Internal_Density=NA)

kmeans_result_df[1,c(2:5)] <- kmeans_result_10
kmeans_result_df[2,c(2:5)] <- kmeans_result_20
kmeans_result_df[3,c(2:5)] <- kmeans_result_30
kmeans_result_df[4,c(2:5)] <- kmeans_result_40
kmeans_result_df[5,c(2:5)] <- kmeans_result_50
kmeans_result_df[6,c(2:5)] <- kmeans_result_LE_comm_num
kmeans_result_df
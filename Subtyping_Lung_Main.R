#R code for "A topological approach for cancer subtyping from gene expression data".
# Ref:
# Rafique, Omar, and A. H. Mir. "A topological approach for cancer subtyping from gene expression data." 
# Journal of Biomedical Informatics 102 (2020): 103357.



#Libraries
  library(Rtsne)
  library(CancerSubtypes)
  library("plot3D")
  library("robustHD")
  library("hyperSpec")
  par(mar=c(1,1,1,1))
  graphics.off()
  set.seed(235)
  
#Obtain the lung gene expression data as matrix. The data is in the form of samples as cols and genes(featues) as rows
  lung_data <- as.matrix(read.delim("Lung/LUNG_Gene_Expression.txt"))
  survivalData = read.table("Lung/LUNG_Survival.txt",sep="\t",header=TRUE)

  dim(lung_data)
  dim(survivalData)
  d1 <- (dim(lung_data)[1])
  d2 <- (dim(lung_data)[2])
#remove the first column because we dont need it
  dummy1 <- lung_data[1:d1,2:d2]
#SOme preprocessing
  index=which(is.na(dummy1))
  dummy2 <- data.imputation(dummy1,fun="median")
  dummy = data.normalization(dummy2,type="feature_Mean",log2=FALSE)
  dim(dummy)
#transpose the data
  lung_mat <- t(as.matrix(dummy))
  dim(lung_mat)
# remove patients with survival longer than 5 years
  toRem = which(survivalData$Survival>1825)
  if(length(toRem)>0) survivalData =  survivalData[-toRem,]
  if(length(toRem)>0) lung_mat = lung_mat[-toRem,]
  dim(lung_mat)
  dim(survivalData)
  plot(lung_mat)
  
#Pearson distance calculation  
  pd = pearson.dist(lung_mat)
  str(pd)
  
  
#Apply t-SNE dimensionlaity reduction  
  lung_tsne <- Rtsne(pd, dim=1, 
                     perplexity =5,
                     theta = 0.715,
                     #check_duplicates = FALSE,
                     pca = FALSE,
                     partial_pca = FALSE,
                     max_iter = 1000,
                     verbose = getOption("verbose", FALSE), 
                     is_distance = TRUE,
                     Y_init = NULL,
                     pca_center = FALSE, 
                     pca_scale = FALSE,
                     normalize = FALSE)
  #str(lung_tsne)
  
  plot(lung_tsne$Y[,1])
  fil <- as.data.frame(lung_tsne$Y[,1])

  
#Some more libraries  
  library(TDAmapper)
  library(ggplot2)

  
  ########################################################################
#We compute that mapper graph for the data with t-SNE filter as follows:

  
  
  map <- mapper1D(distance_matrix = pd,
                  filter_values = fil,
                  num_intervals = 4,
                  percent_overlap = 10,
                  num_bins_when_clustering = 10)
  
  
  
  map
  
  
  
  
#Plotting the mapper graph
  library(igraph)
  First.Example.graph <- graph.adjacency(map$adjacency, mode="undirected")
  plot(First.Example.graph, layout = layout.auto(First.Example.graph) )
  
#####################################################################
#Now lets make groups. This is a manual step as we have to observe which clusters are formed from the graph.
#G11 stores the data smaples which are in a cluster. The cluster may be a vertex or a group of vertices
#g1 is the same as G11 when the repeated elements are removed
#Group_1 is the same as g1 but with elements sorte din order

  G11 <- c(map$points_in_vertex[[1]])
  g1 <- unique(G11)
  Group_1 <- sort(g1, decreasing = FALSE)
  
  G22 <- c(map$points_in_vertex[[2]])
  g2 <- unique(G22)
  Group_2 <- sort(g2, decreasing = FALSE)
  
  G33 <- c(map$points_in_vertex[[3]])
  g3 <- unique(G33)
  Group_3 <- sort(g3, decreasing = FALSE)
  
  G44 <- c(map$points_in_vertex[[4]])
  g4 <- unique(G44)
  Group_4 <- sort(g4, decreasing = FALSE)
  
  
  
#Now lets work with the clinical data
  str(survivalData)
  sur_dat <- as.matrix(survivalData)
  
#Now extract the TIME and STATUS from the survival data. They are generally the 2nd and 3rd column.
#We ned to extract all rows....number of rows = dim(lung_mat)[1]
  dim_lm1 <- (dim(lung_mat)[1])
  sur_data <- sur_dat[1:dim_lm1,2:3]

#below sur_data[Group_1,] extracts the rows of surv_data, the indicies of whihc are given in Group_1
#therefore Group.1 will contain all the rows of the original survival data that were categorised in the first group above.
  Group.1 <-sur_data[Group_1,]
#Below, we create a vector of all ones(equal to the number of elements in group 1)
  if(is.null(dim(Group.1)[1])) {
    print("dim(Group.1)[1] = NULL")
    array_null <- append(array_null,0)
    next
  }
  ones <- c()
  ones[1:dim(Group.1)[1]] <- 1
  
#below we combine the survival data with its cluster numbering
  G1 <- cbind(Group.1,ones)
  
#Repeat the above for all groups/clusters
  
#GROUP 2
  Group.2 <-sur_data[Group_2,]
#Below, we create a vector of all ones(equal to the number of elements in group 2)
  if(is.null(dim(Group.2)[1])) {
    print("dim(Group.1)[1] = NULL")
    array_null <- append(array_null,0)
    next
  }
  twos <- c()
  twos[1:dim(Group.2)[1]] <- 2
#below we combine the survival data with its cluster numbering
  G2 <- cbind(Group.2,twos)
#Group 3 
  Group.3 <-sur_data[Group_3,]
  #Below, we create a vector of all ones(equal to the number of elements in group 3)
  if(is.null(dim(Group.3)[1])) {
    print("dim(Group.3)[1] = NULL")
    array_null <- append(array_null,0)
    next
  }
  threes <- c()
  threes[1:dim(Group.3)[1]] <- 3
  #below we combine the survival data with its cluster numbering
  G3 <- cbind(Group.3,threes)

#Group 4  
  Group.4 <-sur_data[Group_4,]
  #Below, we create a vector of all ones(equal to the number of elements in group 4)
  if(is.null(dim(Group.4)[1])) {
    print("dim(Group.4)[1] = NULL")
    array_null <- append(array_null,0)
    next
  }
  fours <- c()
  fours[1:dim(Group.4)[1]] <- 4
  #below we combine the survival data with its cluster numbering
  G4 <- cbind(Group.4,fours)
  ##############################################
#below we make our final data object containing Time, status and event
  data_for_kme1 <- rbind(G1,G2,G3,G4) #this may contain more elements than the original dataset as some elements are shared between vertices
  
  
# Write CSV in R:write the survival dat alongwith a column of cluster indicies in t .csv file
  write.csv(data_for_kme1, file = "1_Lung_Clustered_survival_TSNE_TDA.csv")
  str(data_for_kme1)
# Loading
  library("survminer")
# Fit survival curves
  require("survival")
#Read the survivla data from the csv file
  lung2 <- read.csv("1_Lung_Clustered_survival_TSNE_TDA.csv")
  str(lung2)
  colnames(lung2)[4] <- c("C")
  lung2[1:5,1:4]
#fit the survivla data into the function
  fitQ <- survfit(Surv(Survival,Death) ~ C, data = lung2)
#calculate Pvalue
  p <- surv_pvalue(fitQ)
  pval <- c(p[2])
  
# Drawing Survival curves
  ggsurvplot(fitQ,
             pval = TRUE, 
             palette = c("red","blue","green","black"),title="A (Lung)", xlab ="Time(Days)", ylab="Survival Probability",  break.time.by =250,
             surv.plot.height=.6, font.x = c(23,  "black"),font.y = c(23,  "black"),
             font.tickslab = c(18, "plain", "black"),pval.size=7)
  
  # Reduce(intersect, list(Group_1,Group_2))
  # Reduce(intersect, list(Group_1,Group_3))
  # Reduce(intersect, list(Group_2,Group_3))
  # Reduce(intersect, list(Group_1,Group_4))
  # Reduce(intersect, list(Group_2,Group_4))
  # Reduce(intersect, list(Group_3,Group_4))
  
  #########################################################
  
#The code below build on the code by Angy89 available at https://github.com/angy89/RobustSparseCorrelation
  
  
  require(survival)
  require(survRM2)
  
  ## Computes L1 based measures of separation between survival curves
  ## note: higher ==> more separation
  ##
  ## Inputs:
  ##    data = data.frame with samples units on rows and two colums with names
  ##           'Survival' that is  survivial time
  ##           'Death'    that is a 0-1 variable indicating death
  ##
  ## Outputs: a list with multiple separation measures
  ##
  ##
  survival.curves.separation <- function(data, cluster, tau){
    
    ans <- list()
    tmp           <- table(cluster)
    ClusterLabel  <- as.numeric(names(tmp))
    ClusterSize   <- as.numeric(tmp)
    K             <- length(ClusterLabel)
    n             <- sum(ClusterSize)
    
    ## Set a time grid
    time_grid  <- min(data$Survival) : min(max(data$Survival), tau)
    time_scale <- 1/diff(range(time_grid))
    
    ## Create estimated survival curve matrix: [time x clusters]
    H <- matrix(NA, nrow=length(time_grid), ncol=K)
    colnames(H) <- paste('Cluster_', ClusterLabel, sep='')
    
    ## Estimate the KM curve clusterwise on a common support
    for(k in 1:K){
      ## Compute Kaplan-Meier estimator on the kth cluster
      km    <- survfit(Surv(Survival, Death)~1, data = data[cluster==ClusterLabel[k] , ] )
      ## Construct the KM estimator function
      KMfun <- stepfun(x=km$time[-1], y=km$surv)
      H[,k] <- KMfun(time_grid)
    }
    
    ## construct matrix of pairwise L1 distances
    D <- matrix(0, ncol=K, nrow=K)
    for (i in 1:K){
      for(j in 1:K)
        if(i!=j){
          D[i,j] <- D[j,i] <- sum( abs( H[ , i]  -  H[ , j] ))
        }
    }
    ## Some scaling is given so that these numbers are somewhow interpretable
    ## for the same number of clusters independently of the time interval
    D <- D * time_scale
    
    
    ## Metric 1: min pairwise L1 distance
    iut <- which(upper.tri(D, diag=FALSE))
    ans$L1min <- min(D[iut])
    
    ## Metric 2: compute the summed L1 distance of each
    ## cluster to the nearest one
    diag(D)   <- NA
    ans$L1sum <- sum( D[iut] ) / length(iut)
    
    
    return(ans)
  }
  
  
  ## Computes the discrepancy between survival curves in temrs of RMST
  ## (restricted mean survival time).
  ##
  ## For all these measures: higher ==> more separation
  ##
  ## Inputs:
  ##    data = data.frame with samples units on rows and two colums with names
  ##           'Survival' that is  survivial time
  ##           'Death'    that is a 0-1 variable indicating death
  ##           tau        truncation time, default 5 years, but if the min is less
  ##                      this is set to minimum across groups
  ##
  ## Outputs: a list with multiple separation measures
  ##
  ##
  rmst.separation <- function(data, cluster, tau) {
    ans <- list()
    tmp           <- table(cluster)
    ClusterLabel  <- as.numeric(names(tmp))
    ClusterSize   <- as.numeric(tmp)
    K             <- length(ClusterLabel)
    n             <- sum(ClusterSize)
    
    
    ## Compute the minimum of the largest observed time in each of the two groups
    max.time <- rep(0, K)
    for (k in 1:K) {
      max.time[k] <- max(data$Survival[cluster == ClusterLabel[k]])
    }
    TAU  <- min(max.time, tau)
    
    
    ## Names
    ##    * RMST = restricted mean survival time
    ##    * LER  = life expectancy ratio
    ##
    ## LER: Life Expectancy Ratio  Matrix
    ##      LER[i,j] = LER[j,i] = max{RMST[i] / RMST[j], RMST[j] / RMST[i]}
    ##      note that here we don't have a baseline group so we define the ratio
    ##      always using in the denominator the group that have smaller RMST
    ##
    ## LED: Life Expectancy Difference
    ##    LED[i,j] = LED[j,i] = abs(RMST[i] - RMST[j])
    ##    note that here we don't have a baseline group so we define tha abs difference
    ##
    LER <- LED <-  matrix(0, ncol = K, nrow = K)
    for (i in 1:K) {
      for (j in 1:K)
        if (i != j) {
          ## First select data from  the two groups
          idx <- { cluster == ClusterLabel[i] | cluster == ClusterLabel[j]  }
          x   <- data[idx,]
          ##  Create a 0-1 vector, with gr==1 if cluster==ClusterLabel[i]
          gr0  <- cluster[idx]
          gr   <- ifelse(gr0 == ClusterLabel[i], 1, 0)
          u    <- rmst2(time = x$Survival, status = x$Death, arm = gr, tau = TAU)
          
          rmst_i <- u$RMST.arm1$rmst[1]
          rmst_j <- u$RMST.arm0$rmst[1]
          
          LER[i,j]  <- LER[j, i] <- max(rmst_i / rmst_j, rmst_j / rmst_i  )
          LED[i, j] <- LED[j, i] <- abs(rmst_i - rmst_j)
        }
    }
    
    ## index of the upper triangle
    iut <- which(upper.tri(LER, diag = FALSE))
    
    ## metric: min of pairwise LER discrepancy
    ans$LERmin <- min(LER[iut])
    
    ## metric: scaled summed pairwise LER discrepancy
    ans$LERsum <- sum(LER[iut]) / length(iut)
    
    ## metric: min of pairwise LED discrepancy
    ans$LEDmin <- min(LED[iut])
    
    ## metric: scaled summed pairwise LED discrepancy
    ans$LEDsum <- sum(LED[iut]) / length(iut)
    
    
    return(ans)
  }
  
  ###
  
  survPval = function(SurvDat,CLg,nYears=5){
    
    fit <- survfit(Surv(Survival, Death) ~ CLg,data = SurvDat,subset = Survival < (365 * nYears))
    suv <- survminer::ggsurvplot(fit, risk.table = TRUE, risk.table.height = 0.5,
                                 xlim = c(0,5000), break.time.by = 500, pval = TRUE)
    pVal <- survminer:::surv_pvalue(fit,method = "survdiff",data=SurvDat)
    return(list(fit=fit,suv = suv,pVal=pVal))
    
  }
  #########################################################
  ########################################################
  ########################################################
  library(survival)
  survRes = list()
  nYears = 5
  tau = 1825
  dim(lung2)
  
  l2d1 <- dim(lung2)[1]
  ClustGroup = lung2[1:l2d1,4:4]
  SurvDat = lung2[1:l2d1,2:3]
  CLg = ClustGroup
  
  LNormDist = c(unlist(survival.curves.separation(SurvDat, CLg,tau)),
                unlist(rmst.separation(SurvDat, CLg,tau)))
  
  str(LNormDist)
  LNormDist[5]
#Cox Analysis
  fitcox <- coxph(Surv(Survival,Death) ~ as.factor(C), data = lung2)
  summary(fitcox)
  ftest <- cox.zph(fitcox)
  ftest
  ggcoxzph(ftest,font.main = 10,font.x = c(8,  "black"),font.y = c(8,  "black"),
               font.tickslab = c(12, "plain", "black"),caption = "A-Lung")
  LNormDist[5]
  #################################
  # #below sur_data[Group_1,] extracts the rows of surv_data, the indicies of whihc are given in Group_1
  # #therefore Group.1 will contain all the rows of the original survival data that were categorised in the first group above.
  Group.A <-lung_mat[Group_1,]
  #Below, we create a vector of all ones(equal to the number of elements in group 1)
  if(is.null(dim(Group.A)[1])) {
    print("dim(Group.A)[1] = NULL")
    array_null <- append(array_null,0)
    next
  }
  ones.A <- c()
  ones.A[1:dim(Group.A)[1]] <- 1
  #below we combine the survival data with its cluster numbering
  GA <- cbind(Group.A,ones.A)
  str(GA)
  
  
  Group.B <-lung_mat[Group_2,]
  if(is.null(dim(Group.B)[1])) {
    print("dim(Group.B)[1] = NULL")
    array_null <- append(array_null,0)
    next
  }
  ones.B <- c()
  ones.B[1:dim(Group.B)[1]] <- 2
  GB <- cbind(Group.B,ones.B)
  str(GA)
  
  
  Group.C <-lung_mat[Group_3,]
  if(is.null(dim(Group.C)[1])) {
    print("dim(Group.C)[1] = NULL")
    array_null <- append(array_null,0)
    next
  }
  ones.C <- c()
  ones.C[1:dim(Group.C)[1]] <- 3
  GC <- cbind(Group.C,ones.C)
  str(GC)
  
  Group.D <-lung_mat[Group_4,]
  if(is.null(dim(Group.D)[1])) {
    print("dim(Group.D)[1] = NULL")
    array_null <- append(array_null,0)
    next
  }
  ones.D <- c()
  ones.D[1:dim(Group.D)[1]] <- 4
  GD <- cbind(Group.D,ones.D)
  str(GD)
  GD[1:5,1:5]
  Group.D[1:5,1:5]
  str(Group.D)
  ##############################################
#Enrichment Analysis
  library(otrimle) 
  library(survRM2)
  library(survminer)
  library(survival)
  require(org.Hs.eg.db)
  library(limma)
  require(clusterProfiler)
  
  
  source("survival_separation_5y.R") #File from https://github.com/angy89/RobustSparseCorrelation
  source("suvival_analysis_on_clusterings.R") #File from https://github.com/angy89/RobustSparseCorrelation
  source("gene_mapping.R") #File from https://github.com/angy89/RobustSparseCorrelation
  source("pathway_analysis.R") #File from https://github.com/angy89/RobustSparseCorrelation
  
  clusterR <- c(ones.A,ones.B,ones.C,ones.D)
  str(clusterR)
  reddata <- rbind(Group.A,Group.B,Group.C,Group.D)
  str(reddata)
  reddata[1:2,1:3] 
  #detach("dplyr")
  FR = pathway_analysis(cluster=clusterR,reddata)
  View(FR$FR)
  View(FR$FR$ID)
  #Dotplot of relevant pathways
  FR$dp  

 

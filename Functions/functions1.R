## written by: Andrew Hinton (andrew84@email.unc.edu)

calcLogRatio <-
function(df,Weighted=F){
  df1 = finalModelProcess(featureMatrix = df[,-1],labels = df[,1],
                          permuteLabel = F,permuteFeatures = F)
  df2 = df1$allData
  fn = length(colnames(df2[,-1]))
  mat = matrix(rep(0,fn*fn),nrow = fn)
  colnames(mat) = colnames(df2[,-1])
  rownames(mat) =  colnames(df2[,-1])
  mat[lower.tri(mat,diag = T)]=-1000000000000000000999
  mat = subset(melt(mat),value!=-1000000000000000000999)
  mat = unite(data = mat,col = Ratio,sep = "___",c("Var1","Var2"))
  keyRats = separate(data.frame(mat),1,into = c("Num","Denom"),sep = "___",remove = F)
  el_= data.frame(keyRats$Num,keyRats$Denom,keyRats$Ratio)
  df1 = data.frame(Status = df[,1], getLogRatios(z = df2,ratioList  = el_,weighted = Weighted))

  return(df1)
}

catchError <-
function(expr){
  tryCatch(expr,
           error = function(e){
             T
           },
           warning = function(w){
             F
           },
           finally = {
             F
           })

}



dcvScores2 <-
function(raTable.train,nfolds = 5,numKNN = 5,seed_ = 08272008,startCluster = T,cores_ = 10){

  ############ Required Functions
  theMethod2 = function(class1,class2,K,exportToGephi,distrDiff = T,
                        center_=F,plotTF=FALSE,Sim_=FALSE,Normalize=TRUE,freqInterval = .2){
    #===========================================================================
    #============= Class 1 =====================================================
    #===========================================================================
    #calculate LR Means Matrix
    class1_logRatio.means = lrMeans(class1,center = center_) #<----------------
    colnames(class1_logRatio.means) = colnames(class1)
    rownames(class1_logRatio.means) = colnames(class1)
    message('Class1 Log Ratio Means Calculated...')

    #calculate Varation Matrix
    class1_logRatio.varMat = as.matrix(class1)
    class1_logRatio.varMat = acomp(class1_logRatio.varMat)
    class1_logRatio.varMat = variation(class1_logRatio.varMat)#<----------------
    class1_logRatio.varMat[class1_logRatio.varMat<0]=0 #protect againist precision error
    class1_logRatio.varMat.scaled = exp(-sqrt(class1_logRatio.varMat))#<----------------
    message('Class1 Variation Matrix Calculated...')

    #Calculate SNR matrix
    c1.sd = sqrt(class1_logRatio.varMat)
    c1.sd[c1.sd==0]=min(c1.sd[c1.sd>0])
    class1_logRatio.SNR = abs(class1_logRatio.means)/c1.sd#<----------------
    diag(class1_logRatio.SNR) = 0
    diag(c1.sd) = 0

    #===========================================================================
    #============= Class 2 =====================================================
    #===========================================================================
    #calculate LR Means Matrix
    class2.logRatio.means = lrMeans(class2,center = center_) #<----------------
    colnames(class2.logRatio.means) = colnames(class2)
    rownames(class2.logRatio.means) = colnames(class2)
    message('Class2 Log Ratio Means Calculated...')

    #calculate Varation Matrix
    class2.logRatio.varMat = as.matrix(class2)
    class2.logRatio.varMat = acomp(class2.logRatio.varMat)
    class2.logRatio.varMat = variation(class2.logRatio.varMat)#<----------------
    class2.logRatio.varMat[class2.logRatio.varMat<0]=0 #protect againist precision error
    class2.logRatio.varMat.scaled = exp(-sqrt(class2.logRatio.varMat))#<----------------
    message('Class2 Log Ratio Variation Matrix Calculated...')
    #class2.logRatio.totVar = sum(class2.logRatio.varMat[upper.tri(class2.logRatio.varMat)])

    #Calculate SNR matrix
    c2.sd = sqrt(class2.logRatio.varMat)
    c2.sd[c2.sd==0]=min(c2.sd[c2.sd>0])
    class2.logRatio.SNR = abs(class2.logRatio.means)/c2.sd#<----------------
    diag(class2.logRatio.SNR) = 0
    diag(c2.sd) = 0

    ########## F Test ################
    #ftest across all
    comp.df = rbind(class1,class2)
    c12 = acomp(comp.df)
    sst = variation(c12)*(nrow(c12)-1)
    xBar1 = class1_logRatio.means
    xBar2 = class2.logRatio.means
    X_Bar = lrMeans(comp.df)
    SSTR = nrow(class1)*(xBar1-X_Bar)^2+nrow(class2)*(xBar2-X_Bar)^2
    SSE = (nrow(class1)-1)*class1_logRatio.varMat+(nrow(class2)-1)*class2.logRatio.varMat
    MST = sst/(nrow(c12)-1)
    MSTR = SSTR
    MSE = SSE/(nrow(c12)-2)
    MSE[MSE==0]=min(SSE[SSE>0])
    F_ = MSTR/MSE
    diag(F_) =0
    F_ = F_/max(F_)
    message("F-statistic Calculated")
    ###########################################


    #===========================================================================
    #============= Differential Analysis=====================================================
    #===========================================================================

    message('Differential Analysis Started...')
    diffMean = abs(class1_logRatio.means - class2.logRatio.means)
    diffVar = abs(c1.sd - c2.sd)
    diffSNR = abs(class1_logRatio.SNR -class2.logRatio.SNR)

    if(distrDiff==T){
      distrDistance = distComparison(class1,class2,intervalSize_ = freqInterval,Center_  = center_)
    }

    if(Normalize==TRUE){
      #zscore
      diffSNR = (diffSNR-mean(diffSNR))/sd(diffSNR)
      diffVar = (diffVar-mean(diffVar))/sd(diffVar)
      diffMean = (diffMean-mean(diffMean))/sd(diffMean)
      F_ = (F_-mean(F_))/sd(F_)

      if(distrDiff==T){
        distrDistance =(distrDistance-mean(distrDistance))/sd(distrDistance)
      }
    }

    if (distrDiff==T){
      list( diffMean = diffMean,
            diffSNR  = diffSNR,
            Fstat = F_,
            diffVar = diffVar,
            distrDistance = distrDistance)
    }  else{
      list( diffMean = diffMean,
            diffSNR  = diffSNR,
            Fstat = F_,
            diffVar = diffVar
      )
    }


  }

  theMethod3 = function(df,k,Siml = FALSE,freqInterval_ = .2,cent = F,distrDiff_=T){
    class1 = df%>%
      filter(Status==1)%>%
      dplyr::select(-Status)

    class2 = df%>%
      filter(Status==0)%>%
      dplyr::select(-Status)

    p = theMethod2(class1,class2,K = k,center_ = cent,distrDiff=distrDiff_,exportToGephi = FALSE,plotTF = FALSE,Sim_ = Siml,freqInterval = freqInterval_)
    return(p)
  }

  imputeZeros = function(df.cdata2){
    zrs = apply(df.cdata2, 1 , function(x) sum(!x))
    #consistnet imp Factor; fixed across study
    numer = 1/100000000000
    #numer<min(df.cdata1[df.cdata1])#should be true
    min(df.cdata2)
    #rsums = rowSums(df.cdata1)
    for (r in 1:nrow(df.cdata2)){
      impFactor = numer
      for(c in 1:ncol(df.cdata2)){
        #check if count equal 0
        if(df.cdata2[r,c]==0){
          df.cdata2[r,c]=impFactor
        }
        else{
          df.cdata2[r,c]=df.cdata2[r,c]*(1-impFactor*zrs[r])
        }
      }
    }

    #hist(zrs/ncol(df.cdata2),main="Distribution of %Species Abundance=0 \\n per Sample")
    #print(rowSums(df.cdata2)) #should be all 1's
    #print(min(df.cdata2)== numer)#shlould be true

    return(df.cdata2)
  }

  lrMeans = function(df,center = F){
    #cnames = list(colnames(df),colnames(df))
    z = as.matrix(df)
    if(center==T){
      h1 = colMeans(z)
      h1 = outer(h1,h1)
      h1.mat = lapply(1:ncol(z), function(x)(matrix(rep(t(h1[,x]),nrow(z)),nrow = nrow(z),byrow = TRUE)))
      s=lapply(1:ncol(z), function(x) log(z/z[,x])*h1.mat[[x]])
      s=sapply(1:ncol(z),function(x) colMeans(s[[x]]))
    }else{
      s=lapply(1:ncol(z), function(x) z/z[,x])
      s=sapply(1:ncol(z),function(x) colMeans(log(s[[x]])))
    }

    matrix(s,ncol(z),ncol(z))
  }

  knnADJtoSYM = function(knnADJ_MAT){
    adj = knnADJ_MAT
    adjT = knnADJ_MAT
    adjT  =t(adj)
    lt = adj[lower.tri(adj)]
    uT = adjT[lower.tri(adjT)]
    df = data.frame(lt = lt,ut = uT)
    message("dplyr Start")
    df = df%>%
      mutate(s = if_else(lt==ut,lt,pmax(ut,lt)))
    message("dplyr Finished")

    adj[lower.tri(adj)]=df$s
    adj = t(adj)
    adj[lower.tri(adj)]=df$s
    #isSymmetric(adj1)

    return(adj)
  }

  kNN=function(adj_mat,K,plot_TrueFalse){
    k=K
    sigEst = 200
    dist_mat = 1/adj_mat^2
    #dist_mat = exp(-((dist_mat)^2/sigEst))
    diag(dist_mat) = 0
    #dist_mat <- as.matrix(dist(Mat, method = "euclidean", upper = TRUE, diag=TRUE))
    nrst <- lapply(1:nrow(dist_mat), function(i) k.nearest.neighbors(i, 1/dist_mat, k = k))
    w <- matrix(nrow = dim(dist_mat), ncol=dim(dist_mat)) ## all NA right now
    w[is.na(w)] <- 0 ## populate with 0
    for(i in 1:length(nrst)) for(j in nrst[[i]]) w[i,j] = dist_mat[i,j]
    Adj2=w
    Net=graph.adjacency(Adj2,mode='undirected',weighted=TRUE)
    FinalAdj=get.adjacency(Net,type='both',sparse=FALSE,attr = "weight")
    diag(FinalAdj)=0

    FinalAdj = FinalAdj

    w.graph = graph.adjacency(FinalAdj,mode="undirected",weighted = TRUE)
    w.graph.simplified = igraph::simplify(w.graph, remove.loops = TRUE,
                                          edge.attr.comb = igraph_opt("edge.attr.comb"))
    #Plot if TRUE
    if (plot_TrueFalse == TRUE){
      #l = layout.kamada.kawai(w.graph.simplified)
      l = layout.circle(w.graph.simplified)

      plot(w.graph.simplified,
           layout = l,
           main = "Test",
           vertex.size = 1,
           vertex.label = NA,
           edge.width =1.25,
           edge.color =rgb(1,.25,0,.1)
      )
    }
    eigen_analy = eigen_centrality(w.graph.simplified,weights = E(w.graph.simplified)$weight)
    eigenC_react=as.matrix(eigen_analy$vector)
    wgh.lv = cluster_louvain(w.graph.simplified, weights = E(w.graph.simplified)$weight)

    list(Graph = w.graph.simplified,
         AdjMatrix = FinalAdj,
         EigenVC = eigenC_react,
         LouvainClust = wgh.lv
    )
  }

  knn_graph = function(adj_mat,K,plot_TrueFalse=FALSE,sim_=FALSE){

    #Function Variables ###########################
    plt = plot_TrueFalse
    adj = adj_mat

    adj2 = adj_mat
    k = K

    #Similarity of Dissim
    if(sim_ == TRUE){
      diag(adj)=Inf# for smallest
    }else{
      diag(adj)=0
      adj = -adj
    }

    #Choose K
    k_nn = t(seq(1,k,by=1))

    #Replace with Adj Mat with Ranks

    # for (r in 1:nrow(adj)){
    #   tt= as.matrix(rank(adj[r,]))
    #   tt= ifelse(tt%in%k_nn,1,0)
    #   tt=tt*adj2[r,]
    #   adj2[r,]=tt
    # }


    adj2 = sapply(1:nrow(adj), function(x) (rankFunction(adj[x,],k_nn)))
    adj2 = t(-adj2)
    rownames(adj2) = colnames(adj2)

    diag(adj2)=0

    #Make Matric Symetric
    adj =  knnADJtoSYM(adj2)

    #Create Graph
    w.graph = graph.adjacency(adj,mode="undirected",weighted = TRUE)
    w.graph.simplified = igraph::simplify(w.graph, remove.loops = TRUE,
                                          edge.attr.comb = igraph_opt("edge.attr.comb"))
    #Plot if TRUE
    if (plt == TRUE){
      #l = layout.kamada.kawai(w.graph.simplified)
      l = layout.circle(w.graph.simplified)

      plot(w.graph.simplified,
           layout = l,
           main = "Test",
           vertex.size = 1,
           vertex.label = NA,
           edge.width =1.25,
           edge.color =rgb(1,.25,0,.1)
      )
    }

    # eigen_analy = eigen_centrality(w.graph.simplified,weights = E(w.graph.simplified)$weight)
    # eigenC_react=as.matrix(eigen_analy$vector)
    # wgh.lv = cluster_louvain(w.graph.simplified, weights = E(w.graph.simplified)$weight)

    list(Graph = w.graph.simplified,
         AdjMatrix = adj
         #EigenVC = eigenC_react,
         #LouvainClust = wgh.lv
    )
  }
  #####################################



  ################  Within-Fold DCV ##########################
  #<<< DCV PARAMETERS >>>>>
  ddDist = F
  lrs.df = c()
  lrs_median.df = c()
  dfc_df = data.frame()
  dpart = kfoldDataPartition(raTable.train,nfolds,seed = seed_)
  #<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  #asetuo parallel

  if(startCluster==T){
    stopImplicitCluster()
    cores_ = if_else(cores_>nfolds,nfolds,cores_)
    cl <- makeCluster(cores_)
    registerDoParallel(cl)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())
    caret:::requireNamespaceQuietStop
  }

  if(nfolds>1){
    lrs.df=foreach(lala=1:nfolds,.combine = rbind,
                   .packages = c("dplyr","tibble","tidyverse","compositions","data.table","doParallel","ranger")) %dopar% {


                     # #Run Core DCV Function##############################
                     # for (lala  in 1:nfolds){
                     x = dpart[[lala]][[1]]
                     p= theMethod3(x,numKNN,freqInterval_ = .105,cent = F,distrDiff_ = ddDist)
                     message("Differential Graph Complete")
                     message("\\n")
                     gc()
                     #Get Variation Matrices
                     ph = list()
                     ph[[1]]=p$diffMean
                     ph[[2]] = p$diffSNR
                     ph[[3]] = p$Fstat
                     ph[[4]] = p$diffVar

                     if(ddDist==T){
                       ph[[5]] = p$distrDistance

                       names(ph) = c(
                         "vdmean",
                         "vdsnr",
                         "vdfstat",
                         "vdvar",
                         "vddistr")

                       numInc = 5

                     }else{

                       names(ph) = c(
                         "vdmean",
                         "vdsnr",
                         "vdfstat",
                         "vdvar")

                       numInc = 4
                     }


                     #Flatten Matrices
                     ph[[1]][upper.tri(ph[[1]],diag = T)]=-1000000000000000000999
                     ph_ = subset(melt(ph[[1]]),value!=-1000000000000000000999)
                     dfc = unite(data = ph_,col = Ratio,sep = "___",c("Var1","Var2"))
                     for(x in 2:numInc){
                       ph[[x]][upper.tri(ph[[x]],diag = T)]=-1000000000000000000999
                       ph_ = subset(melt(ph[[x]]),value!=-1000000000000000000999)
                       ph_ = unite(data = ph_,col = Ratio,sep = "___",c("Var1","Var2"))
                       dfc = left_join(dfc,ph_,by = "Ratio")
                     }
                     colnames(dfc) = c("Ratio",names(ph))
                     dfc$fold = lala
                     dfc_df = rbind(dfc_df,dfc)



                     #pca on flatten variation
                     #pc = prcomp(dfc[,-1],center = T,scale. = T)
                     #summary(pc)
                     #pcrot = abs(pc$rotation)
                     #Compute linear combination on first pc
                     #comb = sapply(2:ncol(dfc), function(x) pcrot[x-1]*dfc[,x])
                     comb = sapply(2:ncol(dfc), function(x) dfc[,x])

                     #Combine zscores to rank ratios
                     # lrs = data.frame(Ratio = dfc$Ratio,value=rowSums(comb),fold = lala)
                     # lrs.df = rbind(lrs.df,lrs)
                     #
                     data.frame(Ratio = dfc$Ratio,value=rowSums(comb),fold = lala)

                   }

    lrs.df = lrs.df %>%
      spread(key ="fold",value = "value" )
    lrs.df$rowmean = rowMeans(lrs.df[,-1])
  }else{
    lrs.df=foreach(lala=1:nfolds,.combine = rbind,
                   .packages = c("dplyr","tibble","tidyverse","compositions","data.table","doParallel","ranger")) %dopar% {


                     # #Run Core DCV Function##############################
                     # for (lala  in 1:nfolds){
                     x = dpart[[lala]][[2]]

                     p= theMethod3(x,numKNN,freqInterval_ = .105,cent = F,distrDiff_ = ddDist)
                     message("Differential Graph Complete")
                     message("\\n")
                     gc()
                     #Get Variation Matrices
                     ph = list()
                     ph[[1]]=p$diffMean
                     ph[[2]] = p$diffSNR
                     ph[[3]] = p$Fstat
                     ph[[4]] = p$diffVar

                     if(ddDist==T){
                       ph[[5]] = p$distrDistance

                       names(ph) = c(
                         "vdmean",
                         "vdsnr",
                         "vdfstat",
                         "vdvar",
                         "vddistr")

                       numInc = 5

                     }else{

                       names(ph) = c(
                         "vdmean",
                         "vdsnr",
                         "vdfstat",
                         "vdvar")

                       numInc = 4
                     }


                     #Flatten Matrices
                     ph[[1]][upper.tri(ph[[1]],diag = T)]=-1000000000000000000999
                     ph_ = subset(melt(ph[[1]]),value!=-1000000000000000000999)
                     dfc = unite(data = ph_,col = Ratio,sep = "___",c("Var1","Var2"))
                     for(x in 2:numInc){
                       ph[[x]][upper.tri(ph[[x]],diag = T)]=-1000000000000000000999
                       ph_ = subset(melt(ph[[x]]),value!=-1000000000000000000999)
                       ph_ = unite(data = ph_,col = Ratio,sep = "___",c("Var1","Var2"))
                       dfc = left_join(dfc,ph_,by = "Ratio")
                     }
                     colnames(dfc) = c("Ratio",names(ph))
                     dfc$fold = lala
                     dfc_df = rbind(dfc_df,dfc)



                     #pca on flatten variation
                     #pc = prcomp(dfc[,-1],center = T,scale. = T)
                     #summary(pc)
                     #pcrot = abs(pc$rotation)
                     #Compute linear combination on first pc
                     #comb = sapply(2:ncol(dfc), function(x) pcrot[x-1]*dfc[,x])
                     comb = sapply(2:ncol(dfc), function(x) dfc[,x])

                     #Combine zscores to rank ratios
                     # lrs = data.frame(Ratio = dfc$Ratio,value=rowSums(comb),fold = lala)
                     # lrs.df = rbind(lrs.df,lrs)
                     #
                     data.frame(Ratio = dfc$Ratio,value=rowSums(comb),fold = lala)

                   }

    lrs.df = lrs.df %>%
      spread(key ="fold",value = "value" )
    lrs.df$rowmean = lrs.df[,2]
  }


  if(startCluster==T){
    stopCluster(cl)
  }

  return(list(lrs =lrs.df,dfc = dfc_df))
}





energyMaximization <-
function(mst.df,keyFeatures,labels,sz,nreps_energy = 1e5,eps = 0.1,patience = 50, targetFeats = NULL,earlyStop = F){
  tr_df = keyFeatures

  ## Esta
  d.compdat = parallelDist::parDist(as.matrix(tr_df),method = "euclidean")
  library(energy)
  energy_aitch = etest2(x = d.compdat,sizes = sz$Freq,distance = T,R = nreps.Energy)
  estat = energy_aitch$statistic
  ii = energy_aitch$perms
  maxF = as.numeric((estat-mean(ii))/sd(ii))


  ##
  message("Forward Selection ranking of logratio sets..")
  rollmean.df = data.frame()
  maxRLM = 0
  impTime = 0


  cs = colSums(mst.df)
  featOrder = sort(cs,decreasing = T)
  #base
  cn = names(featOrder[1:3])
  ph = subset(tr_df,select =  cn)
  baseSet = 1:3
  astat = data.frame()
  baseSet.list = list()
  astat = rbind(astat,data.frame(Value = maxF,optF = maxF,numFeats = length(baseSet)))
  baseSet.list[[1]] = baseSet
  bs = 2
  ss = 1

  for(x in 4:length(featOrder)){
    message(x, " of ", nrow(mst.df) )
    newSet = c(baseSet,x)
    cn = names(featOrder[newSet])
    ph = subset(tr_df,select =  cn)

    ##--------------------------##
    ### Compute estat
    ##--------------------------##
    d.compdat = parallelDist::parDist(as.matrix(ph),method = "euclidean")
    energy_aitch = etest2(x = d.compdat,sizes = sz$Freq,distance = T,R = nreps_energy)
    estat = energy_aitch$statistic
    ii = energy_aitch$perms
    newF = as.numeric((estat-mean(ii))/sd(ii))


    diff_ = newF-maxF
    if(diff_>eps){
      baseSet = c(baseSet,x)
      maxF = newF
      # maxF.mn = normF
      # maxF.dsp = normDisp
    }

    baseSet.list[[bs]] = baseSet
    bs = bs + 1


    ##--------------------------##
    ### Energy
    ##--------------------------##
    astat = rbind(astat,data.frame(Value = newF,optF = maxF,numFeats = length(baseSet)))

    #Breaking criteria - if target features are reached
    if(!is.null(targetFeats)){
      if(targetFeats<=length(baseSet)){
        break
      }
    }


    message("#Features = ",length(baseSet),", maxF = ",round(maxF,4)," , impTime = ",impTime)
    #stop if perm test and test statistic has been exceeded
    if(earlyStop){
      if(earlyStop_Fthreshold<maxF){
        break
      }
    }

    if(diff_>eps){
      impTime = 0
    }else{
      impTime = impTime+1
    }
    if(impTime>patience){
      break
    }
  }

  ### Post Process and store performance data
  avec = c(astat$Value)
  plot(avec)
  avec.df =astat
  mx = which.max(astat$Value)
  phh = mst.df[mx,mst.df[mx,]==1]
  cn = names(featOrder[baseSet.list[[mx]]])
  feats = subset(tr_df,select =  cn)

  return(list(features = feats,energyMaximization = astat,energy = avec))
}

etest2 <-
function (x, sizes, distance = FALSE, method = c("original",
                                                          "discoB", "discoF"), R)
{
  method <- match.arg(method)
  if (method == "discoB" || method == "discoF") {
    g <- as.factor(rep(1:length(sizes), sizes))
    return(disco(x, factors = g, distance = distance, index = 1,
                 R = R, method = method))
  }
  nsamples <- length(sizes)
  if (nsamples < 2)
    return(NA)
  if (min(sizes) < 1)
    return(NA)
  if (!is.null(attr(x, "Size")))
    distance <- TRUE
  x <- as.matrix(x)
  if (NROW(x) != sum(sizes))
    stop("nrow(x) should equal sum(sizes)")
  if (distance == FALSE && nrow(x) == ncol(x))
    warning("square data matrix with distance==FALSE")
  d <- NCOL(x)
  if (distance == TRUE)
    d <- 0
  str <- "Multivariate "
  if (d == 1)
    str <- "Univariate "
  if (d == 0)
    str <- ""
  e0 <- 0
  repl <- rep(0, R)
  pval <- 1
  b <- .C("ksampleEtest", x = as.double(t(x)), byrow = as.integer(1),
          nsamples = as.integer(nsamples), sizes = as.integer(sizes),
          dim = as.integer(d), R = as.integer(R), e0 = as.double(e0),
          e = as.double(repl), pval = as.double(pval), PACKAGE = "energy")
  names(b$e0) <- "E-statistic"
  sz <- paste(sizes, collapse = " ", sep = "")
  methodname <- paste(str, length(sizes), "-sample E-test of equal distributions",
                      sep = "")
  dataname <- paste("sample sizes ", sz, ", replicates ",
                    R, sep = "")
  e <- list(call = match.call(), method = methodname, statistic = b$e0, perms = b$e,
            p.value = b$pval, data.name = dataname)
  class(e) <- "htest"
  e
}


fastImputeZeroes <-
function(df.cdata2, impFactor = 1e-11){

  require(doParallel)

  cn = colnames(df.cdata2)
  rn = rownames(df.cdata2)
  df.cdata2 = as.matrix(compositions::clo(df.cdata2))

  if(is.null(impFactor)){
    impFactor =min(df.cdata2[df.cdata2>0])/10
  }

  df_imputed = foreach(i  = 1:nrow(df.cdata2),.combine = rbind) %dopar% {
    sampleData = df.cdata2[i,]
    nz = sum(sampleData==0)
    sampleData[sampleData==0] = impFactor
    sampleData[sampleData!=0] =  sampleData[sampleData!=0]*(1-impFactor*nz)
  }

  colnames(df_imputed) = cn
  rownames(df_imputed) = rn
  return(df_imputed)
}


finalDCV <-
  function(logRatioMatrix,includeInfoGain = F,nfolds = 5,numRepeats = 1,seed_ = 08272008){

    cvDCV = data.frame()
    for(r in 1:numRepeats){
      foldData = kfoldDataPartition(df = logRatioMatrix,
                                    kfold = nfolds,
                                    permuteLabel = F,
                                    seed = r)

      for(f in 1:nfolds){
        ####################################################
        ##  select fold data for train and test splilt
        trainData = (foldData[[f]]$xtrain_combinedFolds[,-1])
        ytrain = factor(foldData[[f]]$xtrain_combinedFolds[,1])
        classes = as.character(unique(ytrain))
        #####################################################
        #Group 1
        g1 = trainData[ytrain==classes[1],]
        g1Means = colMeans(g1)
        g1Var = matrixStats::colVars(as.matrix(g1))
        names(g1Var) = names(g1Means)
        n1 = nrow(g1)
        #Group 2
        g2 = trainData[ytrain==classes[2],]
        g2Means = colMeans(g2)
        g2Var = matrixStats::colVars(as.matrix(g2))
        names(g2Var) = names(g2Means)
        n2 = nrow(g2)
        #### Metrics
        ## Tstat
        tstat = abs( (g1Means-g2Means) / sqrt( (g1Var/n1) + (g2Var/n2) ) )
        ## F-ratio
        sm = colMeans(trainData)
        expVar =n1*(g1Means - sm)^2 + n2*(g2Means - sm)^2
        unexpVar = g2Var + g1Var
        fRatio = expVar / unexpVar
        #AIC glm
        y = data.frame(rbind(g1,g2))
        ## KS Statistic
        yy = as.matrix(y)
        levs = unique(ytrain)
        mat1 = yy[ytrain==levs[1],]
        mat2 = yy[ytrain==levs[2],]
        KS.df = data.frame(KS = K_S(mt = mat1,mt2 = mat2))
        ## GLM Deviance
        glmDev = data.frame(glmDev = 1/Rfast::logistic_only(x = as.matrix(yy),y = if_else(ytrain==levs[1],1,0))) ## Deviance sub for IC (BIC/AIC)
        #Combine
        dfc = data.frame(Ratio = names(g1Means),tstat,fRatio,KS.df,glmDev)
        #Information Gain
        if(includeInfoGain){
          ig = FSelectorRcpp::information_gain(x = trainData,y = ytrain,type = "gainratio" ,discIntegers = T)$importance
          ig[is.na(ig)]=0
          dfc$ig = ig
        }
        #Scale
        normDFC = data.frame(Ratio = dfc$Ratio, apply(dfc[,-1],2, scale))
        normDFC$rowmean = rowSums(normDFC[,-1],na.rm = T)
        normDFC$fold = f
        cvDCV = rbind(cvDCV,normDFC)
      }

    }
    ## Aggregate
    dcv = aggregate(rowmean ~ Ratio, data = cvDCV,
                    FUN = function(x) mean(x) )

    return(list(lrs = dcv))

  }


finalDCV1 <-
function(logRatioMatrix,includeInfoGain = F,nfolds = 5,numRepeats = 1,seed_ = 08272008){

  cvDCV = data.frame()
  for(r in 1:numRepeats){
    foldData = kfoldDataPartition(df = logRatioMatrix,
                                  kfold = nfolds,
                                  permuteLabel = F,
                                  seed = r)

    for(f in 1:nfolds){
      ####################################################
      ##  select fold data for train and test splilt
      trainData = (foldData[[f]]$xtrain_combinedFolds[,-1])
      ytrain = factor(foldData[[f]]$xtrain_combinedFolds[,1])
      classes = as.character(unique(ytrain))
      #####################################################
      #Group 1
      g1 = trainData[ytrain==classes[1],]
      g1Means = colMeans(g1)
      g1Var = matrixStats::colVars(as.matrix(g1))
      names(g1Var) = names(g1Means)
      n1 = nrow(g1)
      #Group 2
      g2 = trainData[ytrain==classes[2],]
      g2Means = colMeans(g2)
      g2Var = matrixStats::colVars(as.matrix(g2))
      names(g2Var) = names(g2Means)
      n2 = nrow(g2)
      #### Metrics
      ## Tstat
      tstat = abs( (g1Means-g2Means) / sqrt( (g1Var/n1) + (g2Var/n2) ) )
      ## F-ratio
      sm = colMeans(trainData)
      expVar =n1*(g1Means - sm)^2 + n2*(g2Means - sm)^2
      unexpVar = g2Var + g1Var
      fRatio = expVar / unexpVar
      #AIC glm
      y = data.frame(rbind(g1,g2))
      ## KS Statistic
      yy = as.matrix(y)
      levs = unique(ytrain)
      mat1 = yy[ytrain==levs[1],]
      mat2 = yy[ytrain==levs[2],]
      KS.df = data.frame(KS = K_S(mt = mat1,mt2 = mat2))
      ## GLM Deviance
      glmDev = data.frame(glmDev = 1/Rfast::logistic_only(x = as.matrix(yy),y = if_else(ytrain==levs[1],1,0))) ## Deviance sub for IC (BIC/AIC)
      ## PlS Weights
      useBayes   <- caret::plsda( x = yy,y = ytrain, ncomp = 1)
      w = useBayes$loading.weights[,]
      #Combine
      dfc = data.frame(Ratio = names(g1Means),tstat,fRatio,KS.df,glmDev,w)
      #Information Gain
      if(includeInfoGain){
        ig = FSelectorRcpp::information_gain(x = trainData,y = ytrain,type = "gainratio" ,discIntegers = T)$importance
        ig[is.na(ig)]=0
        dfc$ig = ig
      }
      #Scale
      normDFC = data.frame(Ratio = dfc$Ratio, apply(dfc[,-1],2, scale))
      normDFC$rowmean = rowSums(normDFC[,-1],na.rm = T)
      normDFC$fold = f
      cvDCV = rbind(cvDCV,normDFC)
  }

  }
  ## Aggregate
  dcv = aggregate(rowmean ~ Ratio, data = cvDCV,
                   FUN = function(x) mean(x) )

  return(list(lrs = dcv,scores = cvDCV))

}


finalModelProcess <-
function(featureMatrix,labels,permuteLabel=F,permuteFeatures=F){


  #Data Preprocess - Closure/Normalization ###########################
  xtrain = clo(featureMatrix)

  ytrain = labels
  if(permuteLabel==T){
    ytrain = sample(ytrain)
  }else{
    ytrain = ytrain
  }
  ############################################################################
  ###  Permute Labels to CHeck for statistically significant pattern ###
  ############################################################################
  classes = as.character(unique(ytrain))
  ytrain = factor(as.vector(as.matrix(ytrain)),levels = c(classes[1],classes[2]))
  if(permuteFeatures==T){
    c.names = colnames(xtrain)
    ncol=length(c.names)
    posCases = xtrain[ytrain==classes[1],]
    negCases = xtrain[ytrain==classes[2],]
    posCases = posCases[,sample(1:ncol,ncol,replace = F)]
    negCases = negCases[,sample(1:ncol,ncol,replace = F)]
    colnames(posCases) = c.names
    colnames(negCases)=c.names
    xtrain=rbind(posCases,negCases)
  }
  ytrain1 = as.factor(if_else(ytrain==classes[1],1,0))



  ############################################################################
  ####    Data Preprocessing - Zero Imputation ####
  ############################################################################
  xtrain1_imputed = fastImputeZeroes(clo(xtrain))
  allTrain.method = data.frame(Status = ytrain1,xtrain1_imputed)
  ############# ZERO IMPUTATION (END) ###########################

  return(list(allData = allTrain.method,
              y_train = ytrain)
  )

}


getFeats_bySparsity <-
function(compmat,mxSparsePercent = .5){
  zeroFeats = compmat==0
  zeroFeats = colSums(zeroFeats)
  sparseThreshold = round(mxSparsePercent*nrow(compmat))
  keep = zeroFeats[zeroFeats<=sparseThreshold]
  return(names(keep))
}

getLogratioFromList <-
function(Ratio,raMatrix,Class){
  Ratio = data.frame(Ratio)
  keyRats = separate(Ratio,1,into = c("Num","Denom"),sep = "___",remove = F)
  el_= data.frame(keyRats$Num,keyRats$Denom,keyRats$Ratio)
  ad = finalModelProcess(raMatrix,labels = Class)
  getLogRatios(ad$allData,el_)
}

getLogRatios <-
function(z,ratioList,weighted = F){
  #Check To make sure status variable is removed
  if(sum(colnames(z)=="Status")>0){
    c = which(colnames(z)=="Status")
    z = z[,-c]
  }

  #Create feature to location library
  features = data.frame(loc = 1:ncol(z),f = colnames(z))

  #Map Num and denom feature locations
  num = data.frame(f = ratioList[,1])
  num = left_join(num,features)
  den = data.frame(f = ratioList[,2])
  den = left_join(den,features)

  #num and denom feature table with location
  c = cbind(num,den)

  if(weighted==T){
    #weights
    cm = colMeans(z)
    #compute logratios
    s=lapply(1:nrow(c), function(x) cm[c[x,2]] * cm[c[x,4]] * log( z[,c[x,2]] / z[,c[x,4]] ) )
    s = as.data.table(s)
    colnames(s) = as.character(ratioList[,3])
    s = as.data.frame(s)
  }else{
    #compute logratios
    s=lapply(1:nrow(c), function(x) log( z[,c[x,2]] / z[,c[x,4]] ) )
    s = as.data.table(s)
    colnames(s) = as.character(ratioList[,3])
    s = as.data.frame(s)
  }

return(s)
}


imputeZeros <-
function(df.cdata2){
  zrs = apply(df.cdata2, 1 , function(x) sum(!x))
  #consistnet imp Factor; fixed across study
  numer = 1/100000000000
  #numer<min(df.cdata1[df.cdata1])#should be true
  min(df.cdata2)
  #rsums = rowSums(df.cdata1)
  for (r in 1:nrow(df.cdata2)){
    impFactor = numer
    for(c in 1:ncol(df.cdata2)){
      #check if count equal 0
      if(df.cdata2[r,c]==0){
        df.cdata2[r,c]=impFactor
      }
      else{
        df.cdata2[r,c]=df.cdata2[r,c]*(1-impFactor*zrs[r])
      }
    }
  }

  #hist(zrs/ncol(df.cdata2),main="Distribution of %Species Abundance=0 \\n per Sample")
  #print(rowSums(df.cdata2)) #should be all 1's
  #print(min(df.cdata2)== numer)#shlould be true

  return(df.cdata2)
}


kfoldDataPartition <-
function(df,kfold,seed = 08272008,permuteLabel = F){
  set.seed(seed)

  if(permuteLabel==T){
    df[,1] = sample(df[,1])
  }
  f = createFolds(df[,1],k = kfold,list = FALSE)
  df$fold = f
  df_ = df
  df_ = data.frame(ID = rownames(df_),df_)

  #Partition Data
  foldData = list()
  ph_list = list()

  for (j in 1:kfold){
    xtrain = df%>%
      filter(fold != j)%>%
      dplyr::select(-fold)

    xtrain_ID = df_%>%
      filter(fold != j)%>%
      dplyr::select(ID,1,2,fold)

    xtest = df%>%
      filter(fold == j)%>%
      dplyr::select(-fold)

    xtest_ID = df_%>%
      filter(fold == j)%>%
      dplyr::select(ID,1,2,fold)

    ph_list[[1]] = xtrain
    names(ph_list)[1] = "xtrain_combinedFolds"
    ph_list[[2]] = xtest
    names(ph_list)[2] = "xtest_kthFold"
    ph_list[[3]] = xtrain_ID
    names(ph_list)[3] = "xtrain_IDs"
    ph_list[[4]] = xtest_ID
    names(ph_list)[4] = "xtest_IDs"

    foldData[[j]] = ph_list
    names(foldData)[j] = paste("fold_",j,sep = "")
  }

  return(foldData)
}



mstAll <-
function(featMatrix,dcvRanking){
  dcvRanking = dcvRanking %>%
    filter(Ratio %in% colnames(featMatrix))
  keyRats = separate(dcvRanking,1,into = c("Num","Denom"),sep = "___",remove = F)
  el_= data.frame(keyRats$Num,keyRats$Denom,keyRats$Ratio)
  g = igraph::graph_from_edgelist(as.matrix(el_[,1:2]),directed = T)
  E(g)$weight = dcvRanking$rowmean
  g <-  igraph::minimum.spanning.tree(g,weights = -E(g)$weight)
  Ratios = data.frame(igraph::get.edgelist(g,names = T))
  Ratios_na = data.frame(paste(Ratios[,1],Ratios[,2],sep = "___"))
  el_ = cbind(Ratios,Ratios_na)
  colnames(el_)[1:3] = c("Num","Denom","Ratio")

  return(ratios = subset(featMatrix,select = el_$Ratio))
}

mstTest <-
  function( trainData1,ytrain1,featureSetSize = c(5,10,25,35,50,100),testModels = "lda",useTarSelection = T,mxNumRatios=NULL,
                    cv_method = "repeatedcv",num_Repeats = 10,nf = 10,seed_ = 08272008,dcv_Folds = 5 ,dcv_Repeats = 1,
                    permLabels = F,energyMaxParms){


  set.seed(seed_)
  ###---------------------------------------------------------------------------------###
  ## Impute Zeros  ####
  ###---------------------------------------------------------------------------------###
  trainData1 = fastImputeZeroes( compositions::clo(trainData1) , impFactor  = impFact )
  featureList = list()
  energyResults.list =list()
  modelList = list()
  allFeature_energyStats = list()


  ###---------------------------------------------------------------------------------###
  ## Target Feature Set ####
  ###---------------------------------------------------------------------------------###
  if(permLabels){
    yt =  sample(ytrain1)
    ami_ = aricode::AMI(yt,ytrain1)
  }else{
    yt = ytrain1
    ami_ = NA
  }

  ## targeted feature selection
  if(useTarSelection){
    features = rfeSelection_byPart2(gaData = trainData1,
                                    xtest = trainData1,ytest = yt,
                                    #numParts = tarFeats,
                                    dcvFolds_ = dcv_Folds,dcvRepeats_ = dcv_Repeats,
                                    ytrain = yt,
                                    nreps = 3,useInfoGain = energyMaxParms$useIG)
  }else{
    featureSetSize = ncol(trainData1)
  }


  ## Compute Metrics by Feature Set Size
  perf = data.frame()
  resampledAUC = list()
  cc = 1
  for(fts in featureSetSize){
    ###------------------------------###
    ## Top N
    ###------------------------------###
    tarFeats = fts#ncol(dat[,-1])
    if(useTarSelection){
      nodeStrength1 = features %>%
        top_n(tarFeats,Str)
      raMatrix = trainData1[,nodeStrength1$Node]
    }else{
      raMatrix = trainData1
    }

    ###------------------------------###
    ## MST and Energy Max
    ###------------------------------###
    tbl = data.frame(Status = yt,raMatrix)
    test = testDCV1(raMatrix = tbl,fullTable = trainData1,
                    permLabels = F,
                    targetFeats = mxNumRatios,
                    seed = seed_,scale_Estat = energyMaxParms$scale_Estat,
                    patience = energyMaxParms$patienceInterations,
                    nreps_energy = energyMaxParms$nreps.Energy,
                    infoGain_ = energyMaxParms$useIG,
                    mxFeats = energyMaxParms$maxFeatures,
                    dcvMethod = energyMaxParms$dcvMethod,
                    dcvFolds = dcv_Folds,dcvReps = dcv_Repeats
    )
    trainData2 = test$keyFeatures
    estat =  max(test$aitchVec$optF)
    energyResults.list[[cc]] = test
    normE = normalizedEnergy(lrMat = trainData2,labels = test$Labels)
    allFeature_energyStats[[cc]] = test$allFeatures_Estats



    ###------------------------------###
    ## Train
    ###------------------------------###
    y_train = test$Labels
    mdls = trainML_Models(trainLRs =  data.frame(trainData2),
                          testLRs =  data.frame(trainData2),
                          ytrain = y_train, y_test = y_train,
                          cvMethod = cv_method,mtry_ = 1,numFolds = nf,numRepeats = num_Repeats,
                          testIDs = NULL,
                          bagModels = F,sampleSize = round(percentBag*nrow(trainData1)),
                          models = testModels)


    ###------------------------------###
    ## Results
    ###------------------------------###
    ph = data.frame(mdls$performance)
    resampledAUC[[cc]] = mdls$models
    featureList[[cc]] = data.frame(Status = y_train,trainData2)
    modelList[[cc]] = mdls
    cc = cc+1
    ph$feats = fts
    ph$estat = estat
    ph$raw_estat = normE$Estat
    ph$norm_esat = normE$H
    ph$rawE = test$rawEstat
    ph$ami = ami_
    ph$numRatios = ncol(trainData2)
    parts = str_split(colnames(trainData2),pattern = "___",simplify = T)
    parts = unique(c(parts[,1],parts[,2]))
    ph$numParts = length(parts)
    if(permLabels){
      ph$labels = "Permuted"
    }else{
      ph$labels = "True"
    }
    ph$seed = seed_
    message(fts)
    perf = rbind(perf,ph)
  }


  return(list(features = featureList ,allFeature_Energy = allFeature_energyStats,
              performance = perf,resampledMetricList = resampledAUC,energyResults  = energyResults.list,finalModel = modelList))

  }


normalizedEnergy <-
function(lrMat,labels){
  mat = data.frame(Labels = labels,lrMat) %>%
    arrange(desc(Labels))
  classes = unique(labels)
  cc = combinat::combn2(as.character(classes))
  W = data.frame(Count = (table(labels)),Weight = clo(table(labels)))
  N = sum(W$Count.Freq)

  Estat = c()
  H = c()
  for(c in 1:nrow(cc)){
    l = c(cc[c,1],cc[c,2])
    mat.ph = mat[mat$Labels%in%l,]
    mat.ph$Labels = factor(mat.ph$Labels)
    levs = data.frame(Labels = unique(mat.ph$Labels))
    labels = mat.ph$Labels
    tb = data.frame((table(mat.ph$Labels )))
    colnames(tb)[1] = "Labels"
    sz = left_join(levs,tb)
    ## Compute Distance
    d = parallelDist::parDist(as.matrix(mat.ph[,-1]),method = "minkowski",p=1.99) ## p here is eqv to the alpha in the (rizzo/szekely - energy distanmce advanced review feb 2016 review paper)
    d = as.matrix(d)
    #energy::eqdist.e(x = d,distance = T,sizes = sz$Freq,method = "discoB")

    ## Compute A
    a1  = 1:sz[1,2]
    a2  = (sz[1,2]+1):nrow(mat.ph)
    A = d[a1,a2]
    A = as.vector(A)
    A = mean(A)
    ## Compute B
    B = d[a1,a1]
    B= as.vector(B)
    B = sum(B) / length(a1)^2
    ## Compute C
    C = d[a2,a2]
    C = as.vector(C)
    C = mean(C)
    ## Compute and Norm Energy
    E = 2*A-B-C
    T_ = (length(a1)*length(a2)) / ((length(a1)+length(a2)))
    estat =E*T_
    h = E/(2*A)
    ## Compute weights
    w = sum(sz$Freq)/(2*N)
    ## Weight stats
    estat = w*estat
    h = h*w

    H[c] = h
    Estat[c] = estat

  }

  h = sum(H)
  e = sum(Estat)


  return(list(H = h, Estat = e))

}


rfeSelection_byPart <-
function(gaData = trainData1,
         xtest = testData1,
         ytest,useFeatureRFE = F,
         ytrain = ytrain1,
         numParts = 10,
         nreps_ = 3,numSets = 10,minPerc = 0.05 ,
         ntrees = 5000,
         m_try = 2,useStepwise = F,useInfoGain = T){

  message("Compute DCV Features")
  alldat = finalModelProcess(featureMatrix = gaData,labels = ytrain)
  ## Updated DCV 2
  dcv = finalDCV(logRatioMatrix =  calcLogRatio(alldat$allData),includeInfoGain = useInfoGain )
  dcv = dcv$lrs
  ## Scale Scores
  trainData.md = preProcess(data.frame(dcv$rowmean),method = "range",rangeBounds = c(0,1) )
  scaledScore = predict(trainData.md, data.frame(dcv$rowmean))
  #node strength
  el = data.frame(Ratio = dcv$Ratio,Score = scaledScore[,1])
  el = separate(data = el,col = 1,into = c("num","denom"),sep = "___",remove = F)
  g = graph_from_edgelist(as.matrix(el[,2:3]))
  E(g)$weight = el$Score
  nodeStrength = data.frame(Node = names(strength(g)),Str = strength(g)) %>%
    arrange(desc(Str))
  nodeStrength1 = nodeStrength %>%
    top_n(numParts,Str)
  trainLRs1 = calcLogRatio(df = data.frame(Status = ytrain,gaData[,nodeStrength1$Node]))
  #aucScore = replicate(nreps,expr = rangerFunction(trainLRs1,ntrees,mtry_,trainLRs1[,1]),simplify = T)
  xtrain1 = calcLogRatio(df = data.frame(Status = ytrain,gaData[,nodeStrength1$Node]))[,-1]
  xtest1 = calcLogRatio(df = data.frame(Status = ytest,xtest[,nodeStrength1$Node]))[,-1]
  # Store DCV Dervied Features
  dcvFeatures = list(train = xtrain1,test = xtest1,
                     raMatrix.train = gaData[,nodeStrength1$Node],
                     raMatrix.test = xtest[,nodeStrength1$Node] )  ### *******
  message('Compute RFE on DCV Features')

  # DCV Features with FS
  dcvFeatures_wFS = NULL
  if(useFeatureRFE){
    if(numParts<3){
      dcvFeatures_wFS = NULL
    }else{
      suppressMessages(suppressWarnings({
        fs = rfeSelection.ByMetric(train_ratio = data.frame(xtrain1),mtry_ = m_try,
                                   test_ratio = xtest1,
                                   ytrain,sets = numSets,nreps = nreps_,minPercentFeatReturn = minPerc)
      }))

      xtrain1 = fs$reducedTrainRatios
      xtest1 = fs$reducedTestRatio
      dcvFeatures_wFS = list(train = xtrain1,test = xtest1,raMatrix = gaData[,nodeStrength1$Node])  ### *******
    }
  }

  rfeFeatures = NULL
  rfeFeatures_wFS = NULL


  ## Stepwise Backwards Elimination
  if(useStepwise){
    #Base AUC
    trainLRs1 = calcLogRatio(df = data.frame(Status = ytrain,gaData))
    aucScore = replicate(nreps,expr = rangerFunction(trainLRs1,ntrees,m_try,trainLRs1[,1]),simplify = T)
    mean(aucScore)
    ###
    message('Compute stepwise backwards elimination')
    loss.df = foreach(i = 1:ncol(gaData),.combine = rbind,.export = c("calcLogRatio","finalModelProcess","getLogRatios","fastImputeZeroes"))%dopar%{
      suppressMessages(suppressWarnings({
        trainLRs1 = calcLogRatio(df = data.frame(Status = ytrain,gaData[,-i]))
        aucScore = replicate(nreps,expr = rangerFunctionLoss(trainLRs1,ntrees,m_try,trainLRs1[,1]),simplify = T)
        #loss.df = rbind(loss.df,data.frame(Feature = colnames(gaData)[i],Loss = mean(aucScore)))
      }))
      data.frame(Feature = colnames(gaData)[i],Loss = mean(aucScore))
    }
    i = which.max(loss.df$Loss)
    ph = gaData[,-i]
    nc = ncol(ph)
    while(nc>numParts){
      ph1 = ph
      loss.df =  data.frame()
      loss.df = foreach(i = 1:ncol(ph1),.combine = rbind,.export = c("calcLogRatio","finalModelProcess","getLogRatios","fastImputeZeroes"))%dopar%{
        suppressMessages(suppressWarnings({
          trainLRs1 = calcLogRatio(df = data.frame(Status = ytrain,ph1[,-i]))
          aucScore = replicate(nreps,expr = rangerFunctionLoss(trainLRs1,ntrees,m_try,trainLRs1[,1]),simplify = T)
          #loss.df = rbind(loss.df,data.frame(Feature = colnames(ph1)[i],Loss = mean(aucScore)))
          data.frame(Feature = colnames(ph1)[i],Loss = mean(aucScore))
        }))
      }
      k = which.max(loss.df$Loss)
      ph = ph1[,-k]
      nc = ncol(ph)
      message(nc)
    }
    ## Final Loss
    trainLRs1 = calcLogRatio(df = data.frame(Status = ytrain,ph))
    aucScore = replicate(nreps,expr = rangerFunctionLoss(trainLRs1,ntrees,m_try,trainLRs1[,1]),simplify = T)
    mean(aucScore)
    ## Test
    testPh = xtest[,colnames(ph)]
    xtrain1 = calcLogRatio(df = data.frame(Status = ytrain,ph))[,-1]
    xtest1 = calcLogRatio(df = data.frame(Status = ytest,testPh))[,-1]
    rfeFeatures = list(train = xtrain,test = xtest) ### *******
    ## Additional Feature Selection
    message('Compute stepwise baskwards elimination RFE...')
    suppressMessages(suppressWarnings({
      fs = rfeSelection.ByMetric(train_ratio = xtrain1,
                                 test_ratio = xtest1,
                                 ytrain,sets = 10,nreps = 3,
                                 minPercentFeatReturn = .05)
    }))
    fs$trainPerformance
    rfeFeatures_wFS = list(train = fs$reducedTrainRatios,test = fs$reducedTestRatio) ### *******
  }


  return(list(features = nodeStrength1,DCV = dcv,featureRFE = rfeFeatures,
              featureRFE_wFS = rfeFeatures_wFS,
              dcv_Feats = dcvFeatures,
              dcv_featsFS = dcvFeatures_wFS))

}

rfeSelection_byPart2 <-
function(gaData = trainData1,
          xtest = testData1,
          ytest,useFeatureRFE = F,
         seed = 08272008,dcvFolds_,dcvRepeats_,
          ytrain = ytrain1,
          numParts = 10,
          nreps_ = 3,numSets = 10,minPerc = 0.05 ,
          ntrees = 5000,
          m_try = 2,useStepwise = F,useInfoGain = T){

message("Compute DCV Features")


# ## Updated DCV 2
# dcv = finalDCV(logRatioMatrix =  calcLogRatio(alldat$allData),includeInfoGain = useInfoGain )
# dcv = dcv$lrs
#
#pairwise componnets
labels = ytrain
c = combinat::combn2(as.character(unique(labels)))
#data.frame
dat_ = data.frame(Status = ytrain,fastImputeZeroes(gaData))
lgRatios = calcLogRatio(dat_)

#DCV
comb_lrs = foreach(i = 1:nrow(c),.combine = rbind,.packages = c("igraph","compositions","dplyr","foreach","reshape2"))%do%{
cc = c[i,]
ph = lgRatios %>%
filter(Status %in% cc)
message("Compute DCV Scores...",cc)
dcv = finalDCV(logRatioMatrix = ph,includeInfoGain = useInfoGain,seed_ = seed,nfolds = dcvFolds_ ,numRepeats = dcvRepeats_ )
lrs = dcv$lrs
data.frame(Ratio = lrs$Ratio,Imp = lrs$rowmean,pxComparison = i,lrs[,-1])
}

lrs_true = comb_lrs %>%
dplyr::group_by(Ratio) %>%
dplyr::summarise(rowmean = mean(Imp))
dcv = data.frame(lrs_true)


## Scale Scores
trainData.md = preProcess(data.frame(dcv$rowmean),method = "range",rangeBounds = c(0,1) )
scaledScore = predict(trainData.md, data.frame(dcv$rowmean))
#node strength
el = data.frame(Ratio = dcv$Ratio,Score = scaledScore[,1])
el = separate(data = el,col = 1,into = c("num","denom"),sep = "___",remove = F)
g = graph_from_edgelist(as.matrix(el[,2:3]))
E(g)$weight = el$Score
nodeStrength = data.frame(Node = names(strength(g)),Str = strength(g)) %>%
arrange(desc(Str))

return(nodeStrength)
}




trainML_Models <-
  function(trainLRs,testLRs,ytrain,y_test,testIDs=NULL,
           models = c("ranger","svmRadial","gbm","pls"),
           startCluster = F,numcores = 10,numRepeats = 3,numFolds = 5,
           cvMethod = "repeatedcv",ntrees = 2000,mtry_ = 1,
           bagModels = F,sampleSize,
           seed = 08272008){

    finModels = list()
    train_control <- trainControl(method=cvMethod,
                                  repeats = numRepeats,
                                  number=numFolds,seeds = NULL,
                                  classProbs = TRUE,
                                  savePredictions = T,
                                  allowParallel = TRUE,
                                  summaryFunction = caret::multiClassSummary
    )

    #Start Parallel Cluster
    if(startCluster==T){
      clus <- makeCluster(numcores)
      registerDoParallel(clus)
      clusterCall(clus, function(x) .libPaths(x), .libPaths())
      caret:::requireNamespaceQuietStop
    }

    performance_ = data.frame()
    i  = 1
    nsamps = nrow(testLRs)
    prediction_matrix_ =data.frame()
    set.seed(seed)
    models_ = paste(models,1:length(models),sep = "")

    for(mdl in models){

      if(bagModels==T){
        bm = data.frame(Status = ytrain,trainLRs)
        bm_i = sample(1:nrow(trainLRs),size = sampleSize,replace = F)
        trainLRs1 = bm[bm_i,-1]
        rownames(trainLRs1) = paste("ID_",bm_i)
        ytrain1 = bm[bm_i,1]
      }else{
        ytrain1 = ytrain
        trainLRs1 = trainLRs
        rownames(trainLRs1) = paste0("ID_",1:nrow(trainLRs))
      }

      if(mdl%in%c("ranger")){
        glm.mdl1=train(x = trainLRs1 ,
                       y = ytrain1,
                       metric = "ROC",
                       max.depth = 0,
                       method = "ranger",num.trees=ntrees,
                       importance = "permutation",
                       tuneGrid = expand.grid(min.node.size=1,splitrule = c("gini"), mtry = c(2,round(sqrt(ncol(trainLRs1))))),
                       trControl = train_control
        )
      }else if(mdl%in%c("xgbTree")){
        glm.mdl1=train(x = trainLRs1 ,
                       y = ytrain1,
                       metric = "ROC",
                       method = "xgbTree",tuneGrid = expand.grid(nrounds = c(50,150),
                                                                 max_depth = c(3,6),
                                                                 eta = c(0.01,.3),
                                                                 gamma = c(0),
                                                                 colsample_bytree = c(0.05,.1,.2),
                                                                 min_child_weight = 1,
                                                                 subsample = c(.8)),
                       trControl = train_control)

      } else if(mdl%in%c("rangerE")){
        glm.mdl1=train(x = trainLRs1 ,
                       y = ytrain1,
                       metric = "ROC",
                       max.depth = 0,
                       method = "ranger",num.trees=ntrees,
                       importance = "permutation",
                       tuneGrid = expand.grid(min.node.size=1,
                                              splitrule = c("extratrees"),
                                              mtry = c(1,2,round(sqrt(ncol(trainLRs1))))),
                       trControl = train_control
        )
      } else if(mdl%in%c("knn")){
        glm.mdl1=train(x = trainLRs1 ,
                       y = ytrain1,
                       metric = "ROC",
                       method = "knn",
                       tuneGrid = expand.grid(k = 2:round(sqrt(x = nrow(trainLRs1)))),
                       trControl = train_control
        )}else {
          glm.mdl1=train(x = trainLRs1 ,
                         y = ytrain1,
                         method=mdl,
                         metric = "ROC",
                         trControl = train_control
          )
        }

      getTrainPerf(glm.mdl1)
      preds = predict.train(glm.mdl1,testLRs, type= "prob")
      preds$model = models_[i]
      if(!is.null(testIDs)){
        preds = cbind(testIDs,preds)
      }
      prediction_matrix_ = rbind(prediction_matrix_,preds)
      ph1 = cbind(data.frame((data.frame(getTrainPerf(glm.mdl1)))))
      ph1$method = models_[i]
      #colnames(ph1) = paste(models[i],colnames(ph1),sep = "_")
      performance_ = rbind(performance_,ph1)



      message(models_[i])
      finModels[[models_[i]]] = glm.mdl1


      i = i+1




    }

    if(startCluster==T){
      stopCluster(clus)
    }


    #ph = data.frame(t(performance_))
    #colnames(ph) = models

    return(list(performance = performance_,
                predictionMatrix = prediction_matrix_,
                models = finModels)
    )

  }







testDCV <-
function(raMatrix,permLabels = F, seed = 08272008,infoGain_ = T,
                   nreps_energy = 300000,dcvWghtd = F,earlyStop=F,
                   earlyStop_Fthreshold,pvalSelect = F,eps = 0.1,
                   rollAverage = 20, patience = 5,nr = 20,mxFeats = 250,
                   assumeF = F,dcvMethod = 1,wDCV =F,targetFeats = NULL,perms_f = 1){

  #Parms (need to clean up and replace within function)
  dat_train1 = raMatrix
  permuteLabels = permLabels
  useInfoGain = infoGain_
  dcvWeighted = dcvWghtd
  adjust_wPval = pvalSelect
  nreps.perm = nr
  maxFeatures = mxFeats

  if(permuteLabels){
    set.seed(seed)
    ysamp = sample(dat_train1[,1])
    AMI_ = aricode::AMI(ysamp,dat_train1[,1])
    dat_train1[,1] = ysamp
  }

  dat_train1$i  = 1:nrow(dat_train1)
  dat_train1 = data.frame(Status = dat_train1[,1], dat_train1[,-1])
  dat_train1 = dat_train1 %>%
    arrange((Status))
  ids = dat_train1$i
  dat_train1 = dat_train1 %>%
    dplyr::select(-i)
  dat_train1$Status = factor(dat_train1$Status)


  levs = data.frame(Labels = unique(dat_train1[,1]))
  labels = dat_train1[,1]
  tb = data.frame((table(dat_train1[,1])))
  colnames(tb)[1] = "Labels"
  sz = left_join(levs,tb)

  message("Compute DCV...")
  #################################################-
  ## DCV Rankings
  #################################################-
  #pairwise componnets
  labels = dat_train1[,1]
  c = combinat::combn2(as.character(unique(labels)))
  #data.frame
  dat_ = data.frame(Status = labels,dat_train1[,-1])
  if(dcvMethod==2){
    lgRatios = calcLogRatio(dat_train1)
  }
  #DCV
  comb_lrs = foreach(i = 1:nrow(c),.combine = rbind,.packages = c("igraph","compositions","dplyr","foreach","reshape2"))%do%{
    cc = c[i,]

    if(dcvMethod==2){
      ph = lgRatios %>%
        filter(Status %in% cc)
    }else{
      ph = dat_ %>%
        filter(Status %in% cc)
      #calculate DCV
      ph.train = finalModelProcess(featureMatrix = ph[,-1],labels = ph[,1],
                                   permuteLabel = F,permuteFeatures = F)
    }

    message("Compute DCV Scores...")

    ###################################-
    #Rank with DCV
    ##################################-
    if(dcvMethod==1){
      dcv =  dcvScores2(raTable.train = ph.train$allData,nfolds = 5,
                        seed_ = 08272008,cores_ = 10,startCluster = F)
    }else if (dcvMethod==2){
      dcv = finalDCV(logRatioMatrix = ph,includeInfoGain = useInfoGain )
      # system.time({
      # dcv1 = evolvedDCV( raMatrix = ph.train$allData[,-1],labels =  ph.train$allData[,1],includeInfoGain = useInfoGain,useWeights = wDCV)
      # })
    }
    ###################################-

    #Should Features scores be weighted
    if(dcvWeighted==T){
      message("Apply weight correction to DCV Scores...")
      dcv = weightedDCV(ratable = ph.train$allData[,-1],lrs = dcv$lrs)
    }else{
      dcv = dcv$lrs
    }
    message("DCV Scores Calculated Sucessfully.")
    #lrs = newDCV(df = ph)$lrs
    lrs = dcv
    data.frame(Ratio = lrs$Ratio,Imp = lrs$rowmean,pxComparison = i,lrs[,-1])
  }
  lrs = comb_lrs %>%
    dplyr::group_by(Ratio) %>%
    summarise_all(.funs = mean)
  lrs_true = comb_lrs %>%
    dplyr::group_by(Ratio) %>%
    dplyr::summarise(rowmean = mean(Imp))
  lrs_true = data.frame(lrs_true)
  #################################################-


  #################################################-
  ## NULL DCV Rankings
  #################################################-
  function2Export = c("fastImputeZeroes","finalModelProcess","dcvScores2","evolvedDCV","weightedDCV","createFolds","kfoldDataPartition","finalDCV")
  if(adjust_wPval==T){
    message("Compute DCV signficance..")
    lrs.df = foreach(jj = 1:nreps.perm,.combine = rbind,.export = function2Export,
                     .packages = c("igraph","reshape2","compositions","foreach","tidyverse"))%dopar%{
                       message(jj)
                       set.seed(jj)
                       #pairwise componnets
                       labels = sample(dat_train1[,1])
                       c = combinat::combn2(as.character(unique(labels)))
                       #data.frame
                       dat_ = data.frame(Status = labels,dat_train1[,-1])
                       #DCV
                       comb_lrs = foreach(i = 1:nrow(c),.export = function2Export,.combine = rbind, c("igraph","reshape2","compositions","foreach","tidyverse"))%do%{
                         if(dcvMethod==2){
                           ph = lgRatios %>%
                             filter(Status %in% cc)
                         }else{
                           ph = dat_ %>%
                             filter(Status %in% cc)
                           #calculate DCV
                           ph.train = finalModelProcess(featureMatrix = ph[,-1],labels = ph[,1],
                                                        permuteLabel = F,permuteFeatures = F)
                         }

                         message("Compute DCV Scores...")

                         ###################################-
                         #Rank with DCV
                         ##################################-
                         if(dcvMethod==1){
                           dcv =  dcvScores2(raTable.train = ph.train$allData,nfolds = 5,
                                             seed_ = 08272008,cores_ = 10,startCluster = F)
                         }else if (dcvMethod==2){
                           dcv = finalDCV(logRatioMatrix = ph,includeInfoGain = useInfoGain )
                           # system.time({
                           # dcv1 = evolvedDCV( raMatrix = ph.train$allData[,-1],labels =  ph.train$allData[,1],includeInfoGain = useInfoGain,useWeights = wDCV)
                           # })
                         }
                         ###################################-

                         #Should Features scores be weighted
                         if(dcvWeighted==T){
                           message("Apply weight correction to DCV Scores...")
                           dcv = weightedDCV(ratable = ph.train$allData[,-1],lrs = dcv$lrs)
                         }else{
                           dcv = dcv$lrs
                         }
                         message("DCV Scores Calculated Sucessfully.")
                         #lrs = newDCV(df = ph)$lrs
                         lrs = dcv
                         data.frame(Ratio = lrs$Ratio,Imp = lrs$rowmean,pxComparison = i)
                       }
                       lrs = comb_lrs %>%
                         group_by(Ratio) %>%
                         dplyr::summarise(rowmean = mean(Imp))
                       lrs = data.frame(lrs)
                       lrs$seed = jj

                       lrs
                     }
    #################################################-
    alpha = 1*(1/(nreps.perm+1))
    lrs.df1  = left_join(lrs.df,lrs_true,by = "Ratio")
    lrs.grandMean = lrs.df %>%
      group_by(Ratio) %>%
      dplyr::summarise(gMean = mean(rowmean),gSD = sd(rowmean))
    lrs.df1 = lrs.df1 %>%
      mutate(signf = rowmean.x>rowmean.y) %>%
      group_by(Ratio) %>%
      dplyr::summarise(p = (sum(signf)+1)/(get("nreps.perm")+1))
    keep = lrs.df1 %>%
      filter(p<=alpha)
    keep = left_join(keep,lrs_true)
    keep = keep %>%
      arrange(desc(rowmean))
    keep = left_join(keep,lrs.grandMean)
    keep = keep %>%
      mutate(score = (rowmean-gMean)/gSD)%>%
      arrange(desc(rowmean))
    lrs.df2 = left_join(lrs.df1,lrs_true)
  }else{
    keep = lrs_true %>%
      arrange(desc(rowmean))
  }
  #######################################################-


  if(nrow(keep)>maxFeatures){
    end = maxFeatures
  }else{
    end  = nrow(keep)
  }

  #Compute LRS
  keep = keep[1:end,]
  el = separate(keep,col = 1,into = c("Num","Denom"),sep = "___")
  el = data.frame(el[,1:2],Ratio = paste0(el$Num,"___",el$Denom) )
  dat.train = finalModelProcess(dat_train1[,-1],labels = dat_train1[,1])
  tr_df = getLogRatios(z = dat.train$allData[,-1],ratioList = el)
  #reorder by rank
  tr_df = tr_df %>%
    select(one_of(keep$Ratio))


  #################################################-
  ##  DCV derived maximum spanning tree ####
  #################################################-
  mst.df = data.frame()
  featureSearch = 1#if_else(is.null(targetFeats),1,2)
  message("Compute DCV derived maximum spanning tree of logratios...")
  switch(featureSearch,
         {
           message("Optimize number of logratios..")

           bool = T
           x = 3
           i = 1
           while(bool){
             ph = mstAll(tr_df[,1:x],dcvRanking = keep)
             rd = vegan::rda(data.frame(easyCODA::CLR(dat_train1[,-1],weight = F)$LR)~.,data = ph)
             varExpl.subset = rd$CCA$tot.chi / rd$tot.chi
             i = i+1
             x = x+1
             if(round(varExpl.subset,2)==1 | x>=end){
               bool = F
             }
             mst.df = rbind(mst.df, data.frame(nFeats = x-1,pres = 1,Ratio = colnames(ph)))
             message(r = i-1,"__varExplained = ",round(varExpl.subset,2),"__NumRatios = ",ncol(ph))
           }
         },
         {
           message("Perform targeted-n logratio search")
           nfeats = 1
           x = 3
           while (nfeats<targetFeats) {
             ph = mstAll(tr_df[,1:x],dcvRanking = keep)
             mst.df = data.frame(nFeats = ncol(ph),pres = 1,Ratio = colnames(ph))
             nfeats = ncol(ph)
             message(x)
             x = x + 1
           }
         }
  )
  mst.df = mst.df %>%
    filter(pres==1)
  mst.df = spread(mst.df,key = "Ratio",value = "pres",fill = 0)
  mst.df = distinct(mst.df[,-1])
  lvl = unique(labels)
  t = table(labels)
  n1 = as.numeric(t[lvl[1]])
  n2 = as.numeric(t[lvl[2]])


  #################################################-
  ##  Forward Selection ranking ####
  #################################################-
  message("Forward Selection ranking of logratio sets..")
  rollmean.df = data.frame()
  maxRLM = 0
  impTime = 0


  cs = colSums(mst.df)
  featOrder = sort(cs,decreasing = T)
  #base
  cn = names(featOrder)
  ph = subset(tr_df,select =  cn)




  ### Final Output
  output = list(numFeats = ncol(ph),IDs = ids,
                keyFeatures = ph,Labels = labels,
                dcvScores = keep,dcvAll = lrs,
                #aitchVec = avec.df,
                mstData = mst.df,
                size = sz
                #metrics = data.frame(astat = aitchStat)
                #metrics = data.frame(astat = aitchStat,  meansF = astat$means[mx],astat$disp[mx])
  )

  return(output)
}



testDCV1 <-
  function(raMatrix,fullTable,dcvTable = NULL,permLabels = F, scale_Estat = "norm",seed = 08272008,infoGain_ = F,dcvFolds=5,dcvReps = 1,nreps_energy = 300000,earlyStop=F,earlyStop_Fthreshold,eps = 1e-5, patience = 5,mxFeats = 250,dcvMethod = 1,targetFeats = NULL){
    set.seed(seed)
    #Parms (need to clean up and replace within function)
    dat_train1 = raMatrix
    permuteLabels = permLabels
    useInfoGain = infoGain_
    maxFeatures = mxFeats

    if(permuteLabels){
      set.seed(seed)
      ysamp = sample(dat_train1[,1])
      AMI_ = aricode::AMI(ysamp,dat_train1[,1])
      dat_train1[,1] = ysamp
    }

    dat_train1$i  = 1:nrow(dat_train1)
    dat_train1 = data.frame(Status = dat_train1[,1], dat_train1[,-1])
    dat_train1 = dat_train1 %>%
      arrange((Status))
    ids = dat_train1$i
    dat_train1 = dat_train1 %>%
      dplyr::select(-i)
    dat_train1$Status = factor(dat_train1$Status)


    levs = data.frame(Labels = unique(dat_train1[,1]))
    labels = dat_train1[,1]
    tb = data.frame((table(dat_train1[,1])))
    colnames(tb)[1] = "Labels"
    sz = left_join(levs,tb)

    message("Compute DCV...")


    #################################################
    ## DCV Rankings
    #################################################
    #pairwise componnets
    labels = dat_train1[,1]
    c = combinat::combn2(as.character(unique(labels)))
    #data.frame
    dat_ = data.frame(Status = labels,dat_train1[,-1])
    if(dcvMethod==2){
      lgRatios = calcLogRatio(dat_train1)
    }

    if(is.null(dcvTable)){
      #DCV
      comb_lrs = foreach(i = 1:nrow(c),.combine = rbind,.packages = c("igraph","compositions","dplyr","foreach","reshape2"))%do%{
        cc = c[i,]

        if(dcvMethod==2){
          ph = lgRatios %>%
            filter(Status %in% cc)
        }else{
          ph = dat_ %>%
            filter(Status %in% cc)
          #calculate DCV
          ph.train = finalModelProcess(featureMatrix = ph[,-1],labels = ph[,1],
                                       permuteLabel = F,permuteFeatures = F)
        }

        message("Compute DCV Scores...")

        ###################################
        #Rank with DCV
        ##################################
        if(dcvMethod==1){
          dcv =  dcvScores2(raTable.train = ph.train$allData,nfolds = 5,
                            seed_ = 08272008,cores_ = 10,startCluster = F)
        }else if (dcvMethod==2){
          dcv = finalDCV(logRatioMatrix = ph,includeInfoGain = useInfoGain,numRepeats = dcvReps,nfolds = dcvFolds ,seed_ = seed)
        }
        ###################################

        message("DCV Scores Calculated Sucessfully.")
        #lrs = newDCV(df = ph)$lrs
        lrs = dcv$lrs
        data.frame(Ratio = lrs$Ratio,Imp = lrs$rowmean,pxComparison = i,lrs[,-1])
      }
      lrs = comb_lrs %>%
        dplyr::group_by(Ratio) %>%
        summarise_all(.funs = mean)
      lrs_true = comb_lrs %>%
        dplyr::group_by(Ratio) %>%
        dplyr::summarise(rowmean = mean(Imp))
      lrs_true = data.frame(lrs_true)
      keep = lrs_true %>%
        arrange(desc(rowmean))
    }else{
      lrs = dcvTable
      ratioNames = colnames(lgRatios[,-1])
      keep = subset(lrs,subset = lrs[,1] %in% ratioNames)
    }




    #################################################


    if(nrow(keep)>maxFeatures){
      end = maxFeatures
    }else{
      end  = nrow(keep)
    }

    #Compute LRS
    keep = keep[1:end,]
    el = separate(keep,col = 1,into = c("Num","Denom"),sep = "___")
    el = data.frame(el[,1:2],Ratio = paste0(el$Num,"___",el$Denom) )
    dat.train = finalModelProcess(dat_train1[,-1],labels = dat_train1[,1])
    tr_df = getLogRatios(z = dat.train$allData[,-1],ratioList = el)
    #reorder by rank
    tr_df = tr_df %>%
      select(one_of(keep$Ratio))


    #################################################-
    ##  DCV derived maximum spanning tree ####
    #################################################-
    mst.df = data.frame()
    featureSearch = 1#if_else(is.null(targetFeats),1,2)
    message("Compute DCV derived maximum spanning tree of logratios...")
    switch(featureSearch,
           {
             message("Optimize number of logratios..")

             bool = T
             x = 3
             i = 1
             while(bool){
               ph = mstAll(tr_df[,1:x],dcvRanking = keep)
               rd = vegan::rda(data.frame(easyCODA::CLR(fullTable[,-1],weight = F)$LR)~.,data = ph)
               varExpl.subset = rd$CCA$tot.chi / rd$tot.chi
               i = i+1
               x = x+1
               if(round(varExpl.subset,2)==1 | x>=end){
                 bool = F
               }
               mst.df = rbind(mst.df, data.frame(nFeats = x-1,pres = 1,Ratio = colnames(ph)))
               message(r = i-1,"__varExplained = ",round(varExpl.subset,2),"__NumRatios = ",ncol(ph))
             }
           },
           {
             message("Perform targeted-n logratio search")
             nfeats = 1
             x = 3
             while (nfeats<targetFeats) {
               ph = mstAll(tr_df[,1:x],dcvRanking = keep)
               mst.df = data.frame(nFeats = ncol(ph),pres = 1,Ratio = colnames(ph))
               nfeats = ncol(ph)
               message(x)
               x = x + 1
             }
           }
    )
    mst.df = mst.df %>%
      filter(pres==1)
    mst.df = spread(mst.df,key = "Ratio",value = "pres",fill = 0)
    mst.df = distinct(mst.df[,-1])
    lvl = unique(labels)
    t = table(labels)
    n1 = as.numeric(t[lvl[1]])
    n2 = as.numeric(t[lvl[2]])


    #################################################-
    ##  Forward Selection ranking ####
    #################################################-
    message("Forward Selection ranking of logratio sets..")
    rollmean.df = data.frame()
    maxRLM = 0
    impTime = 0
    cs = colSums(mst.df)
    featOrder = sort(cs,decreasing = T)


    ###-----------------------------------------*
    ## All features
    ##-----------------------------------------*
    cn = names(featOrder)
    ph = subset(tr_df,select =  cn)
    lbs = labels
    d.compdat = parallelDist::parDist(as.matrix(ph),method = "euclidean")
    ## Stats with all Features
      set.seed(seed)
      energy_aitch = etest2(x = d.compdat,sizes = sz$Freq,distance = T,R = nreps_energy)
      estat = energy_aitch$statistic
      allFeatures_rawEstat = estat ##
      ii = energy_aitch$perms
      allFeatures_scaledF = as.numeric((estat-mean(ii))/sd(ii)) ##
      allFeatures_normF = normalizedEnergy(lrMat = ph,labels = lbs)$H ##
      allFeaturesEnergyStats = list(rawEnergy = allFeatures_rawEstat,scaledEnergy = allFeatures_scaledF,normEnergy = allFeatures_normF )




    ###----------------------------*
    ## Stepwise Forward Selection
    ###------------------------------*
    cn = names(featOrder[1:3])
    ph = subset(tr_df,select =  cn)

    ###----------------------------*
    ### Compute estat
    ###----------------------------*
    lbs = labels
    d.compdat = parallelDist::parDist(as.matrix(ph),method = "euclidean")
    if(scale_Estat == "scale"){
      set.seed(seed)
      energy_aitch = etest2(x = d.compdat,sizes = sz$Freq,distance = T,R = nreps_energy)
      estat = energy_aitch$statistic
      e_stat = estat
      ii = energy_aitch$perms
      newF = as.numeric((estat-mean(ii))/sd(ii))
    }else if (scale_Estat == "norm"){
      newF = normalizedEnergy(lrMat = ph,labels = lbs)$H
      e_stat = newF
    }else{
      energy_aitch = etest2(x = d.compdat,sizes = sz$Freq,distance = T,R = 1)
      newF = energy_aitch$statistic
      e_stat = newF
    }

    maxF = newF
    baseSet = 1:3
    astat = data.frame()
    baseSet.list = list()
    #astat = rbind(astat,data.frame(Value = maxF,means = normF,disp = normDisp,optF = maxF,numFeats = length(baseSet)))
    astat = rbind(astat,data.frame(Value = maxF,optF = maxF,numFeats = length(baseSet)))
    baseSet.list[[1]] = baseSet
    bs = 2
    ss = 1

    for(x in 4:length(featOrder)){
      message(x, " of ", nrow(mst.df) )
      newSet = c(baseSet,x)
      cn = names(featOrder[newSet])
      ph = subset(tr_df,select =  cn)

      ###----------------------------*
      ### Compute estat
      ###----------------------------*
      d.compdat = parallelDist::parDist(as.matrix(ph),method = "euclidean")

      if(scale_Estat == "scale"){
        energy_aitch = etest2(x = d.compdat,sizes = sz$Freq,distance = T,R = nreps_energy)
        estat = energy_aitch$statistic
        e_stat = estat
        ii = energy_aitch$perms
        newF = as.numeric((estat-mean(ii))/sd(ii))
      }else if (scale_Estat == "norm"){
        newF = normalizedEnergy(lrMat = ph,labels = lbs)$H
        e_stat = newF
      }else{
        energy_aitch = etest2(x = d.compdat,sizes = sz$Freq,distance = T,R = 1)
        newF = energy_aitch$statistic
        e_stat = newF
      }


      diff_ = newF-maxF
      if(diff_>eps){
        baseSet = c(baseSet,x)
        maxF = newF
        nullDist = ii
        e_stat = estat
        # maxF.mn = normF
        # maxF.dsp = normDisp
      }

      baseSet.list[[bs]] = baseSet
      bs = bs + 1


      ###---------------------------------------------------*
      ### Energy
      ###--------------------------------------------------*
      astat = rbind(astat,data.frame(Value = newF,optF = maxF,numFeats = length(baseSet)))

      #Breaking criteria - if target features are reached
      if(!is.null(targetFeats)){
        if(targetFeats<=length(baseSet)){
          break
        }
      }

      message("#Features = ",length(baseSet),", maxF = ",round(maxF,4)," , impTime = ",impTime)
      #stop if perm test and test statistic has been exceeded
      if(earlyStop){
        if(earlyStop_Fthreshold<maxF){
          break
        }
      }

      if(diff_>eps){
        impTime = 0
      }else{
        impTime = impTime+1
      }
      if(impTime>patience){
        break
      }


    }

    ###----------------------------*
    ### Post Process and store performance data
    ###----------------------------*
    avec = c(astat$Value)
    #plot(avec)
    avec.df =astat
    mx = which.max(astat$Value)
    phh = mst.df[mx,mst.df[mx,]==1]
    cn = names(featOrder[baseSet.list[[mx]]])
    feats = subset(tr_df,select =  cn)




    ###----------------------------*
    ## <<<< Aitch Statistic >>>> ####
    ###----------------------------*
    aitchStat = astat$Value[mx]
    ###----------------------------*



    ### Final Output
    output = list(numFeats = ncol(feats),IDs = ids,nullDistr_estat = ii,rawEstat = e_stat,allFeatures_Estats = allFeaturesEnergyStats,
                  keyFeatures = feats,Labels = labels,
                  dcvScores = keep,dcvAll = lrs,
                  aitchVec = avec.df,mstData = mst.df,
                  metrics = data.frame(astat = aitchStat)
                  #metrics = data.frame(astat = aitchStat,  meansF = astat$means[mx],astat$disp[mx])
    )

    return(output)
  }





processCompData <-
function(tbl,minPrevalence = 0.9,permuteLabel = F){
  print(table(tbl[,1]))
  minorityClass = names(which.min(table(tbl[,1])))
  majorityClass = as.character(unique(tbl[tbl[,1]!=minorityClass,1]))
  keep = getFeats_bySparsity(tbl[,-1],mxSparsePercent = minPrevalence)
  fcol = colnames(tbl)[1]
  tbl = subset(tbl,select = c(fcol,keep))
  tbl = tbl[rowSums(tbl[,-1])>0,]
  factor = 1
  ph = compositions::clo(tbl[,-1])
  impFact = min(ph[ph>0]) / factor
  tbl = data.frame(Status = factor(tbl[,1]), tbl[,-1]   )
  #permuteLabels
  permuteLabel = F
  if(permuteLabel==T){
    tbl$Status = sample(tbl$Status)
  }
  return(list(processedData = tbl,impFactor = impFact,minClss = minorityClass,majClass = majorityClass))
}




rfeSelection_byPart3 <-
  function(gaData = trainData1,
           xtest = testData1,
           ytest,useFeatureRFE = F,
           seed = 08272008,dcvFolds_,dcvRepeats_,
           ytrain = ytrain1,
           numParts = 10,
           nreps_ = 3,numSets = 10,minPerc = 0.05 ,
           ntrees = 5000,
           m_try = 2,useStepwise = F,useInfoGain = T){

    message("Compute DCV Features")


    # ## Updated DCV 2
    # dcv = finalDCV(logRatioMatrix =  calcLogRatio(alldat$allData),includeInfoGain = useInfoGain )
    # dcv = dcv$lrs
    #
    #pairwise componnets
    labels = ytrain
    c = combinat::combn2(as.character(unique(labels)))
    #data.frame
    dat_ = data.frame(Status = ytrain,fastImputeZeroes(gaData))
    lgRatios = calcLogRatio(dat_)

    #DCV
    comb_lrs = foreach(i = 1:nrow(c),.combine = rbind,.packages = c("igraph","compositions","dplyr","foreach","reshape2"))%do%{
      cc = c[i,]
      ph = lgRatios %>%
        filter(Status %in% cc)
      message("Compute DCV Scores...",cc)
      dcv = finalDCV(logRatioMatrix = ph,includeInfoGain = useInfoGain,seed_ = seed,nfolds = dcvFolds_ ,numRepeats = dcvRepeats_ )
      lrs = dcv$lrs
      data.frame(Ratio = lrs$Ratio,Imp = lrs$rowmean,pxComparison = i,lrs[,-1])
    }

    lrs_true = comb_lrs %>%
      dplyr::group_by(Ratio) %>%
      dplyr::summarise(rowmean = mean(Imp))
    dcv = data.frame(lrs_true)


    ## Scale Scores
    trainData.md = preProcess(data.frame(dcv$rowmean),method = "range",rangeBounds = c(0,1) )
    scaledScore = predict(trainData.md, data.frame(dcv$rowmean))
    #node strength
    el = data.frame(Ratio = dcv$Ratio,Score = scaledScore[,1])
    el = separate(data = el,col = 1,into = c("num","denom"),sep = "___",remove = F)
    g = graph_from_edgelist(as.matrix(el[,2:3]))
    E(g)$weight = el$Score
    nodeStrength = data.frame(Node = names(strength(g)),Str = strength(g)) %>%
      arrange(desc(Str))

    return(list(node_Strength = nodeStrength,DCV = dcv))
  }




mstTest1 <-
  function( trainData1,ytrain1,str_table,dcv_Table = NULL,featureSetSize = c(5,10,25,35,50,100),testModels = "lda",useTarSelection = T,mxNumRatios=NULL,
            cv_method = "repeatedcv",num_Repeats = 10,nf = 10,seed_ = 08272008,dcv_Folds = 5 ,dcv_Repeats = 1,
            permLabels = F,energyMaxParms){


    set.seed(seed_)
    ###---------------------------------------------------------------------------------###
    ## Impute Zeros  ####
    ###---------------------------------------------------------------------------------###
    trainData1 = fastImputeZeroes( compositions::clo(trainData1) , impFactor  = impFact )
    featureList = list()
    energyResults.list =list()
    modelList = list()
    allFeature_energyStats = list()


    ###---------------------------------------------------------------------------------###
    ## Target Feature Set ####
    ###---------------------------------------------------------------------------------###
    if(permLabels){
      yt =  sample(ytrain1)
      ami_ = aricode::AMI(yt,ytrain1)
    }else{
      yt = ytrain1
      ami_ = NA
    }

    ## targeted feature selection
    if(useTarSelection){
      features = str_table
    }else{
      featureSetSize = ncol(trainData1)
    }


    ## Compute Metrics by Feature Set Size
    perf = data.frame()
    resampledAUC = list()
    cc = 1
    for(fts in featureSetSize){
      ###------------------------------###
      ## Top N
      ###------------------------------###
      tarFeats = fts#ncol(dat[,-1])
      if(useTarSelection){
        nodeStrength1 = features %>%
          top_n(tarFeats,Str)
        raMatrix = trainData1[,nodeStrength1$Node]
      }else{
        raMatrix = trainData1
      }

      ###------------------------------###
      ## MST and Energy Max
      ###------------------------------###
      tbl = data.frame(Status = yt,raMatrix)
      test = testDCV1(raMatrix = tbl,fullTable = trainData1,dcvTable = dcv_Table,
                      permLabels = F,
                      targetFeats = mxNumRatios,
                      seed = seed_,scale_Estat = energyMaxParms$scale_Estat,
                      patience = energyMaxParms$patienceInterations,
                      nreps_energy = energyMaxParms$nreps.Energy,
                      infoGain_ = energyMaxParms$useIG,
                      mxFeats = energyMaxParms$maxFeatures,
                      dcvMethod = energyMaxParms$dcvMethod,
                      dcvFolds = dcv_Folds,dcvReps = dcv_Repeats
      )
      trainData2 = test$keyFeatures
      estat =  max(test$aitchVec$optF)
      energyResults.list[[cc]] = test
      normE = normalizedEnergy(lrMat = trainData2,labels = test$Labels)
      allFeature_energyStats[[cc]] = test$allFeatures_Estats



      ###------------------------------###
      ## Train
      ###------------------------------###
      y_train = test$Labels
      mdls = trainML_Models(trainLRs =  data.frame(trainData2),
                            testLRs =  data.frame(trainData2),
                            ytrain = y_train, y_test = y_train,
                            cvMethod = cv_method,mtry_ = 1,numFolds = nf,numRepeats = num_Repeats,
                            testIDs = NULL,
                            bagModels = F,sampleSize = round(percentBag*nrow(trainData1)),
                            models = testModels)


      ###------------------------------###
      ## Results
      ###------------------------------###
      ph = data.frame(mdls$performance)
      resampledAUC[[cc]] = mdls$models
      featureList[[cc]] = data.frame(Status = y_train,trainData2)
      modelList[[cc]] = mdls
      cc = cc+1
      ph$feats = fts
      ph$estat = estat
      ph$raw_estat = normE$Estat
      ph$norm_esat = normE$H
      ph$rawE = test$rawEstat
      ph$ami = ami_
      ph$numRatios = ncol(trainData2)
      parts = str_split(colnames(trainData2),pattern = "___",simplify = T)
      parts = unique(c(parts[,1],parts[,2]))
      ph$numParts = length(parts)
      if(permLabels){
        ph$labels = "Permuted"
      }else{
        ph$labels = "True"
      }
      ph$seed = seed_
      message(fts)
      perf = rbind(perf,ph)
    }


    return(list(features = featureList ,allFeature_Energy = allFeature_energyStats,
                performance = perf,resampledMetricList = resampledAUC,energyResults  = energyResults.list,finalModel = modelList))

  }




featureSlectionPerformance = function(tbl,modelName = "pls",Method_Name = "Feat_Selection",scenario = "Scenario_X",cvFolds = 5,cvRepeats = 10){

  ##-------------------------*
  ## Scaled E
  ##-------------------------*
  nreps_energy = 1e5
  mat = data.frame(tbl) %>%
    arrange(desc(Status))
  classes = unique(tbl[,1])
  cc = combinat::combn2(as.character(classes))
  labels = tbl$Status
  W = data.frame(Count = (table(labels)),Weight = clo(table(labels)))
  N = sum(W$Count.Freq)
  d = parallelDist::parDist(as.matrix(mat[,-1])) ## p here is eqv to the alpha in the (rizzo/szekely - energy distanmce advanced review feb 2016 review paper)
  energy_aitch = etest2(x = d,sizes = W$Count.Freq,distance = T,R = nreps_energy)
  estat = energy_aitch$statistic
  ii = energy_aitch$perms
  allFeatures_scaledF = as.numeric((estat-mean(ii))/sd(ii)) ##

  ## PLS Model
  mdl = trainML_Models(trainLRs =  data.frame(tbl[,-1]),
                       testLRs = data.frame(tbl[,-1]),
                       ytrain = tbl[,1], y_test = tbl[,1],
                       cvMethod = "repeatedcv",mtry_ = 1,
                       numFolds = cvFolds,
                       numRepeats = cvRepeats,
                       testIDs = NULL,
                       models = modelName)



  #dispersion test
  d = dist(tbl[,-1])
  labels = tbl[,1]
  mod = vegan::betadisper(d,group = labels)
  bd = anova(mod)
  # permutest(mod)
  # mod.HSD <- TukeyHSD(mod)
  # plot(mod,)
  # boxplot(mod)
  # plot(TukeyHSD(mod))
  #plot(mod, ellipse = F, hull = FALSE,label = F) # 1 sd data ellipse
  #permanova
  a.df = data.frame(Type = labels)
  pmv = vegan::adonis2(d~Type,data = a.df,permutations = 1000)
  ## ANOSIM
  ano = vegan::anosim(x = d,grouping = labels)

  ## Energy
  mat.ph = data.frame(tbl) %>%
    arrange(desc(Status))
  mat.ph[,1] = factor(mat.ph[,1])
  levs = data.frame(Labels = unique(mat.ph[,1]))
  labels = mat.ph[,1]
  tb = data.frame((table(mat.ph[,1] )))
  colnames(tb)[1] = "Labels"
  sz = left_join(levs,tb)
  d = dist(mat.ph[,-1])
  #energy::eqdist.etest(d,sizes = sz$Freq,distance = T,R = 100,method = "discoB")
  enf = energy::disco(x = d,factors = mat.ph$Status,distance = T,R = 100)
  ##-----------------------------------------*
  ## LR Network ####
  ##-----------------------------------------*
  feature.df = tbl
  #### Kruskall Test
  krus.test_df = data.frame()
  # tbl =  mst.empirical$features[[1]]
  lrs_ = feature.df[,-1]
  cnames = colnames(lrs_)
  for(i in 1:ncol(lrs_)){
    ph = kruskal.test(x = lrs_[,i],g  = factor(feature.df[,1]))
    ph = data.frame(Ratio =cnames[i],pval = ph$p.value,Statistic = ph$statistic )
    krus.test_df = rbind(krus.test_df,ph)
  }
  #False Discovery Rate estimation
  krus.test_df$p.adjust = p.adjust(krus.test_df$pval,method = "BH")
  pval_level = 0.05
  fdrLevel =  max(krus.test_df$p.adjust[krus.test_df$pval <= pval_level])
  krus.test_df$signf = if_else(krus.test_df$p.adjust>0.05,F,T)
  #Importance df
  imp.df = data.frame(Ratio = krus.test_df$Ratio,Imp = krus.test_df$Statistic)
  keyRats = separate(imp.df,1,into = c("Num","Denom"),sep = "___",remove = F)
  el_= data.frame(keyRats$Num,keyRats$Denom,keyRats$Ratio)
  g = igraph::graph_from_edgelist(as.matrix(el_[,1:2]),directed = T)
  g = igraph::simplify(g, remove.loops = TRUE,
                       edge.attr.comb = igraph_opt("edge.attr.comb"))
  E(g)$weight = if_else((imp.df$Imp)<0,0,(imp.df$Imp))
  pp=plot(g,layout = layout_with_fr,vertex.size = log(1/clo(strength(g)))+1,
       vertex.label.cex = .75,edge.curved = .2,edge.width = E(g)$weight*.15,
       edge.arrow.size = .25,edge.arrow.width = 1)


  return(list(
    performance =  data.frame(Method = Method_Name,Scenario = scenario,
                              NumRatios = ncol(tbl[,-1]),NumParts = length(strength(g)),
                              PermanovaF = pmv$F[1],
                              betaDispF = bd$`F value`[1],
                              AnosimR = ano$statistic,
                              EnergyF = enf$statistic,
                              PLSDA_AUC = mdl$performance$TrainAUC,
                              ScaledF = allFeatures_scaledF,
                              netDiamByStrength = diameter(g)/length(strength(g))),
    graph = g,
    graphPlot = pp)
  )

}

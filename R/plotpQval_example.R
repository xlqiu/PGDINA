#' @title Mesa plot for polytomous attribute Q-matrix validation
#'
#' @description The mesa plot was first proposed by de la Torre and Ma (2016) for graphically illustrating the best q-vector(s) for each item.
#' The q-vector on the edge of the mesa is likely to be the best q-vector.
#' In this function, the mesa plot for polytomous attribute Q-matrix is developed.
#'
#' @param x model object of class \code{pQval}
#' @param item a vector specifying which item(s) the plots are drawn for
#' @param type types of the plot. It can be \code{"best"} or \code{"all"}. If \code{"best"},
#'     for all q-vectors requiring the same number of attributes, only the one with the largest PVAF
#'     is plotted, which means \eqn{K_j} q-vectors are plotted; If \code{"all"}, all q-vectors
#'     will be plotted.
#' @param eps the cutoff for PVAF. If not \code{NULL}, it must be a value between 0 and 1. A horizontal line will be drawn accordingly.
#' @param no.qvector the number of q vectors that need to be plotted when \code{type="all"}. The default is 10,
#'        which means the 10 q vectors with the largest PVAFs are plotted.
#' @param data.label logical; To show data label or not?
#' @param original.q.label logical; print the label showing the original q-vector or not?
#' @param auto.ylim logical; create y range automatically or not?
#' @param ... additional arguments passed to \code{plot} function
#' @seealso \code{\link{pQval}}
#' @examples
#'\dontrun{
#' dat <- sim18GDINA$simdat
#' Q <- sim18GDINA$simQ
#' Q[18,] <- c(0,1,1)
#' mod1 <- GDINA(dat = dat,Q = Q,model = "GDINA",mono.constraint = TRUE)
#' out<- pQval(mod,eps=0.95)
#' item <- c(1,18)
#' plot(out,item=item,type="best",eps=0.95)#'
#'}
#'
#' @references
#' de la Torre, J., Qiu, X.-L., & Santos, C. K. (2021). An empirical Q-matrix validation method for the polytomous G-DINA model, \emph{Psychometrika, 87}, 693–724.
#'
#' de la Torre, J., & Ma, W. (2016, August). Cognitive diagnosis modeling: A general framework approach and its implementation in R. A Short Course at the Fourth Conference on Statistical Methods in Psychometrics, Columbia University, New York.
#' @export



plot.pQval <-
  function(x, item, type = "best", no.qvector = 10,
           data.label = TRUE,eps = "auto",
           original.q.label = FALSE,auto.ylim = TRUE,...)

  {
    if(eps=="auto") eps <- x$eps
    #Q <- extract.Qval(x,"Q")
    Q<- x$Q   # misQ
	Qhat<-x$sug.Q[,2:(K+1)]
    K <- ncol(Q)
    patt <- attributepattern(K,Q)  # contain zero patterns
    L <- (nrow(patt)-1)
    fullPVAF <- x$PVAF
    label.PVAF <- apply(patt,1,paste0,collapse = "")

    if (L<no.qvector) no.qvector <- L
      for (j in item){
	    Kj <- rowSums((patt != 0)*1)
	    Mj<- rowSums(patt*1)
		# split data according to Kj
		sp.PVAF.j=split(fullPVAF[,j],Kj[-1])
		# sorted with splited data
        bestPVAF.j <- lapply(sp.PVAF.j,function(x) {x[order(x,decreasing = TRUE)][1:3]})
		label.bestPVAF.j<- lapply(sp.PVAF.j,function(x) {names(x[order(x,decreasing = TRUE)])[1:3]})

		# convert results from a list to a matrix
        bestPVAFj <-matrix(unlist(bestPVAF.j,use.names=FALSE),ncol=1)
		label.bestPVAFj<- matrix(unlist(label.bestPVAF.j,use.names=FALSE),ncol=1)

		if(tolower(type)=="three"){
		ind<-3
		}else if (tolower(type)=="two"){
		ind<-2
		}else if (tolower(type)=="best"){
		ind<-1
		}

		bestPVAFj.tmp<- matrix(bestPVAFj,nrow=3)
		label.bestPVAFj.tmp<- matrix(label.bestPVAFj,nrow=3)
		bestPVAFj<- matrix(bestPVAFj.tmp[1:ind,],ncol=1)
		label.bestPVAFj<-matrix(label.bestPVAFj.tmp[1:ind,],ncol=1)

		bestPVAFj<- rbind(0,bestPVAFj)
        label.bestPVAFj<- rbind(paste0(patt[1,],collapse=""),label.bestPVAFj)

           #if (auto.ylim) ylim = c(max(0,round(min(bestPVAF[,j])-0.1,1)),1) else ylim=c(0,1)
		   if (auto.ylim) ylim = c(max(0,round(min(bestPVAFj)-0.1,1)),1) else ylim=c(0,1)

           graphics::plot(bestPVAFj,xaxt="n",type="o",ylab = "PVAF",xlab="q-vectors",
                       main = paste0("Mesa Plot for Item",j,"(",item.qly,",mis.p=",mis.p,",rep",r,")"),ylim = ylim)
           #graphics::axis(1,at=c(1:length(bestPVAFj)),labels = label.bestPVAFj,las=2)  # vergical axis labels
		   graphics::axis(1,at=c(1:length(bestPVAFj)),labels = label.bestPVAFj,cex.axis=1)

		   bestlocj<- data.frame(matrix(ncol=length(label.bestPVAFj),nrow=1))
	       for (jj in 1:length(label.bestPVAFj)) {
		      bestlocj[jj]<- which(label.bestPVAFj[jj]==label.PVAF)
		  }

        if (!is.null(eps)&&eps>0&&eps<1) abline(h=eps,lty=3);text(1.5,eps+0.03,paste("eps =",eps))
        yloc <- bestPVAFj-diff(ylim)/15
        yloc[yloc<=ylim[1]] <- yloc[yloc<=ylim[1]] + 2 * diff(ylim)/15
        #if (data.label) graphics::text(c(1:length(bestPVAFj)),yloc,round(bestPVAFj,3))
		if (data.label) graphics::text(c(1:length(bestPVAFj)),yloc, formatC( round(bestPVAFj, 3 ), format='f', digits=3 )            )


		locy0 <- which(apply(patt,1,function(x){
          all(x==Qhat[j,])}))
		locy1 <- which(apply(patt,1,function(x){
          all(x==trueQ[j,])}))
        locy2 <- which(apply(patt,1,function(x){
          all(x==Q[j,])}))
        if(locy0%in%bestlocj) graphics::points(which(bestlocj==locy0),bestPVAFj[which(bestlocj==locy0)],col="red",pch=19)
        if(locy1%in%bestlocj) graphics::points(which(bestlocj==locy1),bestPVAFj[which(bestlocj==locy1)],col="red",pch=3)
		#if (original.q.label) text(K-1,ylim[1]+diff(ylim)/6,paste("original q-vector:\n",names(fullPVAF[,j])[locy0]))
		text(length(bestPVAFj)*0.75,(ylim[1]+diff(ylim)/1.5)-0.5,paste("suggested q-vector:\n",label.bestPVAFj[which(bestlocj==locy0)]))
		text(length(bestPVAFj)*0.75,ylim[1]+diff(ylim)/1.5,paste("true q-vector:\n",label.PVAF[locy1]))
		text(length(bestPVAFj)*0.75,(ylim[1]+diff(ylim)/1.5)-0.25,paste("specified q-vector:\n",label.PVAF[locy2]))
      }
	}

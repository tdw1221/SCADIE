## In this file we put in the auxciliary functions used by SCADIE
##dependencies NMF, solveDampenedWLS



update_H <- function(W_input,Y){

  output<-NMF::fcnnls(x=W_input,y =Y )
  output_H <-output$x%*%diag(1/apply(output$x,MARGIN = 2,sum))
  return( output_H )
}


update_H_dwls <- function(W_input,Y){

  output_H <- matrix(NA,nrow=ncol(W_input),ncol=ncol(Y))
  for (i in 1:ncol(output_H)){
    output_H[,i] <- solveDampenedWLS(S =as.matrix(W_input),B=as.matrix(Y[,i],ncol=1) )
  }

  return( output_H )
}



update_W <-function(H_input,Y){

  if ( (length(which( duplicated(Y, MARGIN = 2) ))>0)|(length(which( duplicated(Y, MARGIN = 1) ))>0)|length(which(apply(Y,2,sd)==0))>0| length(which(apply(Y,1,sd)==0)) ){
    stop("Duplicated or constant cols/rows in Y matrix, please check!")
  }

  output<-NMF::fcnnls(x=t(H_input),y =t(Y))

  return( t(output$x) )
}

update_W_bs <-function(H_input,Y){

  try(output<-NMF::fcnnls(x=t(H_input),y =t(Y)), silent=TRUE)
  while (exists("output") == FALSE){
    n = nrow(t(Y)); p=ncol(t(Y)); yz = t(Y) + matrix(rnorm(n*p,0,0.1), nrow=n)
    try(output<-NMF::fcnnls(x=t(H_input),y =yz),silent=TRUE)
  }
  return( t(output$x) )
}






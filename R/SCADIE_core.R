## This file conains the core functions used in SCADIE



#' Implementation of separate W-update without SCAD Penalty
#'
#' This function should be called by parental function and is not supposed to be used individually.
#' @param ini_H_adjusted Initial H
#' @param ini_W_adjusted Initial W
#' @param bulk_expr_sub Bulk Y with only siganture gene rows
#' @param bulk_expr_full Full bulk Y
#' @param n_ct Number of cell types
#' @param itr Number of iteration
#' @param H_update_method ``NNLS" or ``DWLS"
#' @param H_update_gene ``all" or ``signature"
#' @param signature_gene_row_index row indexes of signature genes in full W or Y
#' @param duplicated_rows If there are duplicated rows present in Y
#'
#' @return Returns a list of W and H after iteration, as well as initial W and H
#' @export
#'
#' @examples
SIR_itr_general <- function(ini_H_adjusted,ini_W_adjusted,bulk_expr_sub,bulk_expr_full,n_ct,itr=200,H_update_method="NNLS",H_update_gene="all",signature_gene_row_index,duplicated_rows=F){
  H_tmp <-  ini_H_adjusted
  W_tmp <- ini_W_adjusted


  for (z in 2:itr){
    if (duplicated_rows==T)
    {
      W_tmp <- update_W_bs(H_input =H_tmp,Y = bulk_expr_full )
    }

    else {
      W_tmp <- update_W(H_input =H_tmp,Y = bulk_expr_full )
    }
    if (H_update_method=="NNLS"){

      if (H_update_gene=="all"){
        H_tmp <- update_H(W_input = W_tmp,Y = bulk_expr_full)
      }

      else if (H_update_gene=="signature"){
        H_tmp <- update_H(W_input = W_tmp[signature_gene_row_index,],Y = bulk_expr_full[signature_gene_row_index,])
      }

      stopifnot(H_update_gene=="all"|H_update_gene=="signature")
    }

    else if (H_update_method=="DWLS"){

      if (H_update_gene=="all"){
        H_tmp <- update_H_dwls(W_input = W_tmp,Y = bulk_expr_full)
      }

      else if (H_update_gene=="signature"){
        H_tmp <- update_H_dwls(W_input = W_tmp[signature_gene_row_index,],Y = bulk_expr_full[signature_gene_row_index,])
      }

      stopifnot(H_update_gene=="all"|H_update_gene=="signature")
    }


    stopifnot(H_update_method=="NNLS"|H_update_method=="DWLS")

  }
  return(list(W_end=W_tmp,H_end=H_tmp,W_ini=ini_W_adjusted,H_ini=ini_H_adjusted))
}



#' Implementation of SCAD-penalty itertion
#' This is a sub-function and is supposed to be called by its parental function
#' @param ini_H_adjusted1 Initial H1
#' @param ini_H_adjusted2 Initial H2
#' @param ini_W_adjusted1 Initial W1
#' @param ini_W_adjusted2 Initial W2
#' @param bulk_expr_sub1 Bulk Y1 with only signature gene rows
#' @param bulk_expr_sub2 Bulk Y2 with only signature gene rows
#' @param bulk_expr_full1 Bulk Y1
#' @param bulk_expr_full2 Bulk Y2
#' @param itr number of iteration
#' @param lambda the $\lambda$ controlling the penalty term in objective function
#' @param Weight_mat input weight matrix
#' @param H_update_method ``NNLS" or ``DWLS"
#' @param H_update_gene ``siganture" or ``all"
#' @param signature_gene_row_index signature genes' row indexes in full Y or W
#' @param weight_update_number update E every n times, not used
#' @param cutoff default 10^-6, see the document of ``Iterate_W_H_full_general" for more explanation
#' @param lambda_weight see the document of ``Iterate_W_H_full_general" for more explanation
#'
#' @return End Ws and Hs, initial Ws and Hs, as well as weight matrix E
#' @export
#'
#' @examples
SIR_itr_W_sim_updated_general <- function(ini_H_adjusted1, ini_H_adjusted2, ini_W_adjusted1, ini_W_adjusted2, bulk_expr_sub1,bulk_expr_sub2,bulk_expr_full1, bulk_expr_full2, itr=200, lambda, Weight_mat, H_update_method,H_update_gene,signature_gene_row_index,weight_update_number,cutoff,lambda_weight){
  H_tmp1 <-   ini_H_adjusted1
  H_tmp2 <-   ini_H_adjusted2
  W_tmp1 <- ini_W_adjusted1
  W_tmp2 <- ini_W_adjusted2


  W_dif_Weight <- matrix(NA,nrow=ncol(W_tmp1),ncol=nrow(W_tmp1))
  #browser()
  z = 1
  abs_cutoff = cutoff*((norm(W_tmp1,type="f")+norm(W_tmp2,type="f"))/2)

  W_out <- list()
  W_out$W1 <- 0
  W_out$W2 <- 0

  while (z <= itr & (norm(W_tmp1 - W_out$W1,type="f")>abs_cutoff | norm(W_tmp2 - W_out$W2,type="f")>abs_cutoff)  ){

    if (z > 1){
      W_tmp1 <- W_out$W1
      W_tmp2 <- W_out$W2
    }


    if (z%%weight_update_number == 0){
      scad_der = function(x, thres){
        (x <= thres) + (x>thres) * pmax(3.7*thres - x, 0) / (2.7*thres)
      }

      W_dif_Weight = scad_der(t(abs(W_tmp1 - W_tmp2)), lambda_weight)
    }



      ## prior to W_dif we use the input weight matrix
      if (  is.na(W_dif_Weight[1,1])  ){
        W_out <- update_W_sim_adaptive1_general(H_tmp1, H_tmp2, lambda, Y1=bulk_expr_full1, Y2=bulk_expr_full2, Weight_mat=Weight_mat)
      }
      else{
        W_out <- update_W_sim_adaptive1_general(H_tmp1, H_tmp2, lambda, Y1=bulk_expr_full1, Y2=bulk_expr_full2, Weight_mat=W_dif_Weight)}



    if (H_update_method=="NNLS"){

      if (H_update_gene=="all"){
        H_tmp1 <- update_H(W_input = W_tmp1,Y = bulk_expr_full1)
        H_tmp2 <- update_H(W_input = W_tmp2,Y = bulk_expr_full2)
      }

      else if (H_update_gene=="signature"){
        H_tmp1 <- update_H(W_input = W_tmp1[signature_gene_row_index,],Y = bulk_expr_full1[signature_gene_row_index,])
        H_tmp2 <- update_H(W_input = W_tmp2[signature_gene_row_index,],Y = bulk_expr_full2[signature_gene_row_index,])
      }

      stopifnot(H_update_gene=="all"|H_update_gene=="signature")
    }

    else if (H_update_method=="DWLS"){

      if (H_update_gene=="all"){
        H_tmp1 <- update_H_dwls(W_input = W_tmp1,Y = bulk_expr_full1)
        H_tmp2 <- update_H_dwls(W_input = W_tmp2,Y = bulk_expr_full2)
      }

      else if (H_update_gene=="signature"){
        H_tmp1 <- update_H_dwls(W_input = W_tmp1[signature_gene_row_index,],Y = bulk_expr_full1[signature_gene_row_index,])
        H_tmp2 <- update_H_dwls(W_input = W_tmp2[signature_gene_row_index,],Y = bulk_expr_full2[signature_gene_row_index,])
      }

      stopifnot(H_update_gene=="all"|H_update_gene=="signature")
    }


    stopifnot(H_update_method=="NNLS"|H_update_method=="DWLS")


    z = z + 1
  }

  return(list(W_end1=W_out$W1, W_end2=W_out$W2, H_end1=H_tmp1, H_end2=H_tmp2, W_ini1=ini_W_adjusted1,W_ini2=ini_W_adjusted2,H_ini1=ini_H_adjusted1,H_ini2=ini_H_adjusted2,weight_mat_updated=W_dif_Weight))
}


#' Implementation of SCAD-penalized W update
#'
#' @param H1 cell type proportion matrix for group1
#' @param H2 cell type proportion matrix for group2
#' @param lambda $\lambda$ controllling the SCAD penalty term in the objective function
#' @param Y1 full bulk matrix Y1
#' @param Y2 full bulk matrix Y2
#' @param Weight_mat weight matrix E
#'
#' @return the updated W1 and W2
#' @export
#'
#' @examples
update_W_sim_adaptive1_general = function(H1, H2, lambda, Y1, Y2, Weight_mat){
  r = nrow(H1)
  n1 = ncol(Y1)
  n2 = ncol(Y2)

  Y_w_t1 <- t(Y1)
  Y_w_t2 <- t(Y2)


  output_w <- matrix(nrow=nrow(Y1),ncol=2*r )
  for (i in 1:nrow(output_w)){
    X11 = t(H1)
    X22 = t(H2)
    X31 = sqrt(lambda) * diag(sqrt(Weight_mat[,i]))
    X32 = - sqrt(lambda) * diag(sqrt(Weight_mat[,i]))
    Xi = cbind(X11, array(0, c(n1,r)))
    Xi = rbind(Xi, cbind(array(0,c(n2,r)), X22))
    Xi = rbind(Xi, cbind(X31,X32))
    output_w[i,] <- NMF::fcnnls(y =rbind(matrix(Y_w_t1[,i],ncol=1), matrix(Y_w_t2[,i],ncol=1), array(0,c(r,1))), x = Xi)$x
  }

  return(list(W1=output_w[,1:r], W2 = output_w[, ((r+1):(2*r))]))
}




#' Main function for SCADIE point estimate
#'
#'This is the main interface function for users. Due to the large number of arguments it take in (W1, W2, H1, H2, etc.), it is recommended to wrap them all in a list. As the list will be extremely useful also in the standard error estimate function ``Estimate_sd_general"
#'
#' @param c an auxicilliary argument used only when this function is called by ``Estimate_sd_general", default value is NA and should be kept like that.
#' @param n_ct number of cell types
#' @param input_initial_H1 initial cell type proportion matrix $H_{1}$, each row represents a cell type and each column represents a sample
#' @param input_initial_H2 initial cell type proportion matrix $H_{2}$, each row represents a cell type and each column represents a sample
#' @param input_initial_W1 initial full $W_1$ matrix. Each row represents a gene and each column represents a cell type. Can be obtained by the ``update_W" function with $Y$ and initial $H$ as input. Please be noted that the $W$ is required to have rownames that are each gene's name
#' @param input_initial_W2 initial full $W_2$ matrix. Each row represents a gene and each column represents a cell type. Can be obtained by the ``update_W" function with $Y$ and initial $H$ as input. Please be noted that the $W$ is required to have rownames that are each gene's name
#' @param input_bulk_sub1 initial sub $Y_1$ matrix containing only the signature genes. Each row represents a gene and each column represents a sample.
#' @param input_bulk_sub2 initial sub $Y_2$ matrix containing only the signature genes. Each row represents a gene and each column represents a sample.
#' @param input_bulk_full1 initial full bulk $Y_1$ matrix. Each row represents a gene and each column represents a sample.
#' @param input_bulk_full2 initial full bulk $Y_2$ matrix. Each row represents a gene and each column represents a sample.
#' @param update_W_method the $W$-update method, can be ``NNLS" or ``SCAD", default is ``SCAD"
#' @param var_estimate an auxcilliary argument used when called by `Estimate_sd_general", default is NA and no need to change it when only running this function
#' @param H_update_method the $W$-update method, can be ``NNLS" or ``DWLS", default is ``NNLS"
#' @param H_update_gene in the $H$-update step, what set of genes to used for regression, this can be ``all" for using full $Y$ and $W$, or ``signature" for using only those signature genes, default is ``signature"
#' @param signature_gene_row_index a vector containing signature genes' row indexes in full $Y$ or $W$
#' @param max_itr max iteration number, default 50.
#' @param weight_update_number update the weight matrix every ``weight_update_number" rounds. Noted that this is a new feature that has not been included in our paper. In the current set-up, the weight matrix $E$ is obtained from initial ``warm-up" runs. But to make it more data-adaptive, we can update $E$ during the iteration.
#' @param cutoff the cutoff coefficient for algorithm convergence, the algorithm will stop when the norm difference between previous $W$ and current $W$ is smaller than cutoff*W_norm, default $10^-8$
#' @param eta_weight $\eta_{n}$ for SCAD penalty
#'
#' @return a list with length two, the first list field contains all the outputs related to group 1, including end $W_{1}, H_{1}$, the weight matrix $E$, the updated matrix $E$ (if any), as well as a record of initial $W_{1}$; the second list field contains all the outputs related to group 2, including end $W_{2}, H_{2}$, the weight matrix $E$, the updated matrix $E$ (if any), as well as a record of initial $W_{2}$
#' @export
#'
#' @examples Please refer to package vignette for a detailed walkthrough.
Iterate_W_H_full_general = function(c=NA, n_ct,input_initial_H1,input_initial_H2,input_initial_W1,input_initial_W2,input_bulk_sub1,input_bulk_sub2,input_bulk_full1,input_bulk_full2, update_W_method="SCAD" ,var_estimate=NA,H_update_method="NNLS",H_update_gene="all",signature_gene_row_index,max_itr=50,weight_update_number=999,cutoff=10^-6,eta_weight=4){
  if (is.na(c)){
    column_index1 <- seq(1:ncol(input_initial_H1))
    column_index2 <- seq(1:ncol(input_initial_H2))
  }
  else {
    if (var_estimate=="jackknife")
    {## Create index list
      column_index1 <- setdiff(seq(1:ncol(input_initial_H1) ),c)
      column_index2 <- setdiff(seq(1:ncol(input_initial_H2) ),c)
      input_initial_H1 <- input_initial_H1[,-c]
      input_initial_H2 <- input_initial_H2[,-c]
    }
    else if (var_estimate=="bootstrap")
    {
      ## Create index list
      set.seed(c)
      column_index1 <- sample(1:ncol(input_initial_H1), ncol(input_initial_H1), replace = T)
      column_index2 <- sample(1:ncol(input_initial_H2), ncol(input_initial_H2), replace = T)

      input_initial_H1 <- input_initial_H1[, column_index1]
      input_initial_H2 <- input_initial_H2[, column_index2]
    }

    stopifnot(var_estimate=="bootstrap"|var_estimate=="jackknife"|is.na(c) )

    }

  if (update_W_method=="NNLS")
  {
    ## Iteration by NNLS
    nnls_tmp1_0 = SIR_itr_general(ini_H_adjusted = input_initial_H1,ini_W_adjusted = input_initial_W1,bulk_expr_full = input_bulk_full1[,column_index1],bulk_expr_sub =input_bulk_sub1[,column_index1],n_ct = n_ct, itr = max_itr,H_update_method=H_update_method,H_update_gene=H_update_gene,signature_gene_row_index=signature_gene_row_index )

    nnls_tmp2_0 = SIR_itr_general(ini_H_adjusted = input_initial_H2,ini_W_adjusted = input_initial_W2,bulk_expr_full = input_bulk_full2[,column_index2],bulk_expr_sub =input_bulk_sub2[,column_index2],n_ct = n_ct ,itr = max_itr,H_update_method=H_update_method,H_update_gene=H_update_gene,signature_gene_row_index=signature_gene_row_index )

    rownames(nnls_tmp1_0$W_end) <- rownames(input_bulk_full1)
    rownames(nnls_tmp2_0$W_end) <- rownames(input_bulk_full2)



    nnls_output_vec_1 <- list(W_end=nnls_tmp1_0$W_end,H_end=nnls_tmp1_0$H_end )

    nnls_output_vec_2 <- list(W_end=nnls_tmp2_0$W_end,H_end=nnls_tmp2_0$H_end )


    return(list(output1=nnls_output_vec_1, output2=nnls_output_vec_2))


  }

  else if (update_W_method=="SCAD"){
    if ( is.na(c)  ){
      nnls_tmp1_0 = SIR_itr_general(ini_H_adjusted = input_initial_H1,ini_W_adjusted = input_initial_W1,bulk_expr_full = input_bulk_full1[,column_index1],bulk_expr_sub =input_bulk_sub1[,column_index1],n_ct = n_ct,itr = 3,H_update_method=H_update_method,H_update_gene=H_update_gene,signature_gene_row_index=signature_gene_row_index )

      nnls_tmp2_0 = SIR_itr_general(ini_H_adjusted = input_initial_H2,ini_W_adjusted = input_initial_W2,bulk_expr_full = input_bulk_full2[,column_index2],bulk_expr_sub =input_bulk_sub2[,column_index2],n_ct = n_ct,itr = 3,H_update_method=H_update_method,H_update_gene=H_update_gene,signature_gene_row_index=signature_gene_row_index )

    }
    else if ( var_estimate=="bootstrap" ){
      nnls_tmp1_0 = SIR_itr_general(ini_H_adjusted = input_initial_H1,ini_W_adjusted = input_initial_W1,bulk_expr_full = input_bulk_full1[,column_index1],bulk_expr_sub =input_bulk_sub1[,column_index1],n_ct = n_ct,itr = 3,H_update_method=H_update_method,H_update_gene=H_update_gene,signature_gene_row_index=signature_gene_row_index ,duplicated_rows=T)

      nnls_tmp2_0 = SIR_itr_general(ini_H_adjusted = input_initial_H2,ini_W_adjusted = input_initial_W2,bulk_expr_full = input_bulk_full2[,column_index2],bulk_expr_sub =input_bulk_sub2[,column_index2],n_ct = n_ct,itr = 3,H_update_method=H_update_method,H_update_gene=H_update_gene,signature_gene_row_index=signature_gene_row_index ,duplicated_rows=T)

    }

    else{
      nnls_tmp1_0 = SIR_itr_general(ini_H_adjusted = input_initial_H1,ini_W_adjusted = input_initial_W1,bulk_expr_full = input_bulk_full1[,column_index1],bulk_expr_sub =input_bulk_sub1[,column_index1],n_ct = n_ct,itr = 3,H_update_method=H_update_method,H_update_gene=H_update_gene,signature_gene_row_index=signature_gene_row_index )

      nnls_tmp2_0 = SIR_itr_general(ini_H_adjusted = input_initial_H2,ini_W_adjusted = input_initial_W2,bulk_expr_full = input_bulk_full2[,column_index2],bulk_expr_sub =input_bulk_sub2[,column_index2],n_ct = n_ct,itr = 3,H_update_method=H_update_method,H_update_gene=H_update_gene,signature_gene_row_index=signature_gene_row_index )

    }

    ##getting weight matrix ready

    scad_der = function(x, thres){
      return((x <= thres) + (x>thres) * pmax(3.7*thres - x, 0) / (2.7*thres) )
    }
    W_dif_Weight = scad_der(t(abs(nnls_tmp1_0$W_end - nnls_tmp2_0$W_end)), eta_weight)


    ## Run the scad program

    test_out_tmp <-SIR_itr_W_sim_updated_general(
      ini_H_adjusted1=nnls_tmp1_0$H_end, ini_H_adjusted2=nnls_tmp2_0$H_end,
      ini_W_adjusted1= nnls_tmp1_0$W_end,ini_W_adjusted2= nnls_tmp2_0$W_end,
      bulk_expr_sub1=input_bulk_sub1[,column_index1],bulk_expr_sub2=input_bulk_sub2[,column_index2],
      bulk_expr_full1=input_bulk_full1[,column_index1], bulk_expr_full2=input_bulk_full2[,column_index2],
      itr=max_itr, lambda=0.1, Weight_mat=W_dif_Weight,
     H_update_method=H_update_method,H_update_gene=H_update_gene,signature_gene_row_index=signature_gene_row_index,weight_update_number=weight_update_number,cutoff=cutoff,lambda_weight=eta_weight)

    rownames(test_out_tmp$W_end1) <- rownames(input_bulk_full1)
    rownames(test_out_tmp$W_end2) <- rownames(input_bulk_full2)


    output_vec_scad_1 <- list(W_end=test_out_tmp$W_end1,H_end=test_out_tmp$H_end1,weight=W_dif_Weight,updated_weight = test_out_tmp$weight_mat_updated,ini_M=nnls_tmp1_0$W_end )

    output_vec_scad_2 <- list(W_end=test_out_tmp$W_end2,H_end=test_out_tmp$H_end2,weight=W_dif_Weight,updated_weight = test_out_tmp$weight_mat_updated,ini_M=nnls_tmp2_0$W_end  )


    return(list(output1=output_vec_scad_1, output2=output_vec_scad_2))
  }
}









#' Standard Error Estimation for SCADIE through Jackknife (or Bootstrap)
#'
#'This function estiamtes the entry-wise standard error for W2-W1, which is used in combined with point estimates from `Iterate_W_H_full_general" for DEG identification. The default strategy is through leave-one-out Jackknife, but bootstrap is also implemented.
#'
#' @param input_list A list containing the following fields: bulk_full_1,bulk_full_2, initial_H_1, initial_H_2, initial_W_1, initial_W_2, bulk_sub_1, bulk_sub_2. For a detailed explanation of these variables, please refer to the document of ``Iterate_W_H_full_general"
#' @param update_W_method The method use for W-update, can be "SCAD" or "NNLS", default ``SCAD"
#' @param method ``jackknife" or ``bootstrap", default ``jackknife"
#' @param cores Number of cores used for computating
#' @param jk_subsample Number of jackknife resamplings, the default strategy for jackknife is to exhaustively apply levae-one-out through all samples, the only case we need to specify this argument is when sample size is extrmely large. Default value is NA, which equals the smaller number of samples in group1 and 2.
#' @param bs_num Number of bootstrap re-sampling, it is recommended to use at least 2x sample size
#' @param H_update_method the $W$-update method, can be ``NNLS" or ``DWLS", default is ``NNLS"
#' @param H_update_gene in the $H$-update step, what set of genes to used for regression, this can be ``all" for using full $Y$ and $W$, or ``signature" for using only those signature genes, default is ``signature"
#' @param signature_gene_row_index a vector containing signature genes' row indexes in full $Y$ or $W$
#' @param max_itr max iteration number, default 50.
#' @param weight_update_number update the weight matrix every ``weight_update_number" rounds. Noted that this is a new feature that has not been included in our paper (default 999, which does not update $E$ during iteration). In the current set-up, the weight matrix $E$ is obtained from initial ``warm-up" runs. But to make it more data-adaptive, we can update $E$ during the iteration.
#' @param cutoff the cutoff coefficient for algorithm convergence, the algorithm will stop when the norm difference between previous $W$ and current $W$ is smaller than cutoff*W_norm, default $10^-8$
#' @param eta_weight $\eta_{n}$ for SCAD penalty
#'
#' @return When method=``jackknife", the output is a list list of 4 fields: W1_vec, W2_vec,  w_diff_sd_jackknife and w_diff_sd_jackknife_raw. The W1_vec and W2_vec record all output W1s, W2s throughout the leave-one-out process. w_diff_sd_jackknife_raw.records the sample sd from all resamplings, and w_diff_sd_jackknife is the entry-wise sd. For DEG analysis, w_diff_sd_jackknife should be used. When method=``bootstrap", there would not be w_diff_sd_jackknife_raw.records in the output.
#' @export
#'
#' @examples Please refer to package vignette for a detailed walkthrough.
Estimate_sd_general <- function(input_list, update_W_method="SCAD", method="jackknife",cores,bs_num=500,H_update_method="NNLS",H_update_gene="all",signature_gene_row_index=NULL,jk_subsample=NA,max_itr=50,weight_update_number=999,cutoff=10^-6,eta_weight=4){

  n_sample = min(ncol(input_list$bulk_full_1),ncol(input_list$bulk_full_2))
  n_ct <- ncol(input_list$sig_matrix)
  if (method == "jackknife"){
    if (is.na(jk_subsample)){
      result_jackknife= mclapply(1:n_sample, Iterate_W_H_full_general, n_ct=n_ct,input_initial_H1=input_list$initial_H_1,input_initial_H2=input_list$initial_H_2,input_initial_W1=input_list$initial_W_1,input_initial_W2=input_list$initial_W_2,input_bulk_sub1=input_list$bulk_sub_1,input_bulk_sub2=input_list$bulk_sub_2,input_bulk_full1=input_list$bulk_full_1,input_bulk_full2=input_list$bulk_full_2, update_W_method=update_W_method,var_estimate=method ,H_update_method=H_update_method,H_update_gene=H_update_gene,signature_gene_row_index=signature_gene_row_index,max_itr=max_itr,weight_update_number=weight_update_number,cutoff=cutoff,eta_weight=eta_weight,mc.cores=cores)
      jackknife_vec1<-list(1)
      jackknife_vec2<-list(2)

      for (c in 1:n_sample){
        jackknife_vec1[[c]] = result_jackknife[[c]]$output1
        jackknife_vec2[[c]] = result_jackknife[[c]]$output2
      }

      w_diff_series_jf <- vector("list",length(jackknife_vec1))
      for (i in 1:length(w_diff_series_jf)){
        w_diff_series_jf[[i]] <- jackknife_vec1[[i]]$W_end - jackknife_vec2[[i]]$W_end
      }

      ## estimate element level standard deviation
      w_diff_sd_jackknife_raw <-apply(simplify2array(w_diff_series_jf), 1:2, sd)

      w_diff_sd_jackknife=w_diff_sd_jackknife_raw*sqrt( (n_sample-1)^2/n_sample )

      return(list(W1_vec=jackknife_vec1,W2_vec=jackknife_vec2,w_diff_sd_jackknife=w_diff_sd_jackknife,w_diff_sd_jackknife_raw=w_diff_sd_jackknife_raw))
    }

    else{
      sample_list = sample(seq(1,n_sample),jk_subsample)
      result_jackknife= mclapply(sample_list, Iterate_W_H_full_general, n_ct=n_ct,input_initial_H1=input_list$initial_H_1,input_initial_H2=input_list$initial_H_2,input_initial_W1=input_list$initial_W_1,input_initial_W2=input_list$initial_W_2,input_bulk_sub1=input_list$bulk_sub_1,input_bulk_sub2=input_list$bulk_sub_2,input_bulk_full1=input_list$bulk_full_1,input_bulk_full2=input_list$bulk_full_2, update_W_method=update_W_method,var_estimate=method ,H_update_method=H_update_method,H_update_gene=H_update_gene,signature_gene_row_index=signature_gene_row_index,max_itr=max_itr,weight_update_number=weight_update_number,cutoff=cutoff,eta_weight=eta_weight,mc.cores=cores)
      jackknife_vec1<-list(1)
      jackknife_vec2<-list(2)

      for (c in 1:length(sample_list)){
        jackknife_vec1[[c]] = result_jackknife[[c]]$output1
        jackknife_vec2[[c]] = result_jackknife[[c]]$output2
      }

      w_diff_series_jf <- vector("list",length(jackknife_vec1))
      for (i in 1:length(w_diff_series_jf)){
        w_diff_series_jf[[i]] <- jackknife_vec1[[i]]$W_end - jackknife_vec2[[i]]$W_end
      }

      ## estimate element level standard deviation
      w_diff_sd_jackknife_raw <-apply(simplify2array(w_diff_series_jf), 1:2, sd)

      w_diff_sd_jackknife=w_diff_sd_jackknife_raw*sqrt( (length(sample_list)-1)^2/length(sample_list) )

      return(list(W1_vec=jackknife_vec1,W2_vec=jackknife_vec2,w_diff_sd_jackknife=w_diff_sd_jackknife,w_diff_sd_jackknife_raw=w_diff_sd_jackknife_raw))


    }
  }

  else if (method == "bootstrap"){
    ## bootstrap sampling will create duplicated rows, be careful
    result_bootstrap= mclapply(1:bs_num, Iterate_W_H_full_general, n_ct=n_ct,input_initial_H1=input_list$initial_H_1,input_initial_H2=input_list$initial_H_2,input_initial_W1=input_list$initial_W_1,input_initial_W2=input_list$initial_W_2,input_bulk_sub1=input_list$bulk_sub_1,input_bulk_sub2=input_list$bulk_sub_2,input_bulk_full1=input_list$bulk_full_1,input_bulk_full2=input_list$bulk_full_2, update_W_method=update_W_method,var_estimate=method ,H_update_method=H_update_method,H_update_gene=H_update_gene,signature_gene_row_index=signature_gene_row_index,max_itr=30,weight_update_number=weight_update_number,cutoff=cutoff,eta_weight=eta_weight,mc.cores=cores)
    bootstrap_vec1<-list(1)
    bootstrap_vec2<-list(2)
    for (c in 1:bs_num){
      bootstrap_vec1[[c]] = result_bootstrap[[c]]$output1
      bootstrap_vec2[[c]] = result_bootstrap[[c]]$output2
    }

    w_diff_series_bs <- vector("list",length(bootstrap_vec1))
    for (i in 1:length(w_diff_series_bs)){
      w_diff_series_bs[[i]] <- bootstrap_vec1[[i]]$W_end - bootstrap_vec2[[i]]$W_end
    }

    ## estimate element level standard deviation
    w_diff_sd_bootstrap <-apply(simplify2array(w_diff_series_bs), 1:2, sd)

    w_diff_sd_bootstrap=w_diff_sd_bootstrap*1
    return(list(W1_vec=bootstrap_vec1,W2_vec=bootstrap_vec2,w_diff_sd_bootstrap=w_diff_sd_bootstrap) )
  }

  else{print("Warning: method can only be 'jackknife' or 'bootstrap'")}


}









#' Function to generate the MAP solution corresponding to patterns in Zp
#' Inputs:
#'@param C: 1xD string with data types, D = number of dimensions
#'@param Zp: PxK matrix of feature patterns for which it computes the MAP estimate
#'    (P is the number of feature patterns)
#' @param hidden: structure with latent variables learned by the model
#'@param B: latent feature matrix (D * K * maxR)  where
#'@param D: number of dimensions
#'@param K: number of latent variables
#'@param maxR: maximum number of categories across all dimensions
#'@param  mu: 1*D shift parameter
#' @param w:  1*D scale parameter
#' @param  theta: D*maxR matrix of auxiliary vars (for ordinal variables)
  # ----------------(optional) ------------------
  #' @param     - idxsD: dimensions to infer
#' @return X_map: P*Di matrix with MAP estimate where Di = length(idxsD)
    
GLFM_computeMAP<-function(C,Zp,hidden,params,varargin){
  # We need to call the transformation functions! 
  source("f_o.R")
  source("f_g.R")
  source("f_c.R")
  source("f_p.R")
  source("f_o.R")
  source("f_n.R")
  if (length(varargin) == 1){
    idxsD <- varargin[1]
  }
  else if(length(varargin) > 1){
    stop('Too many input arguments')
  }
  else{
    idxsD <- 1:dim(hidden$B)[1]
  }
  P <- dim(Zp)[1]
  K <- dim(hidden$B)[2]
  if(dim(Zp)[2]!= K){
    stop('Incongruent sizes between Zp and hidden.B: number of latent variables should not be different')
  }
  X_map<-matrix(0,nrow=P,ncol=length(idxsD))
  # For each dimension
  for(dd in 1:length(idxsD)){ # for each dimension
  d <- idxsD(dd)
  #('g','p','n','c','o')
  X_map[,dd]<-f_g(Zp*unlist(hidden$B[d]),hidden$mu[d],hidden$w[d])
  X_map[,dd]<-f_p(Zp*unlist(hidden$B[d]),hidden$mu[d],hidden$w[d])
  X_map[,dd]<-f_n(Zp*unlist(hidden$B[d]),hidden$mu[d],hidden$w[d])
  switch(C(d),'g'={X_map[,dd]<-f_g(Zp*unlist(hidden$B[d]),hidden$mu[d],hidden$w[d])},'p'={'goodbye'},'n'={'no'})
  }
}


  
    # d = idxsD(dd);
    # if isfield(params,'t') % if external transformations have been defined
    # if ~isempty(params.t{d}) % there is an external transform for data type d
    # C(d) = params.ext_dataType{d}; % set new type of data
    # end
    # end
    # switch C(d)
    # case 'g', X_map(:,dd) = f_g( Zp * squeeze(hidden.B(d,:,1))', hidden.mu(d), hidden.w(d) );
    #                              case 'p', X_map(:,dd) = f_p( Zp * squeeze(hidden.B(d,:,1))', hidden.mu(d), hidden.w(d) );
    # case 'n', X_map(:,dd) = f_n( Zp * squeeze(hidden.B(d,:,1))', hidden.mu(d), hidden.w(d) );
    #                              case 'c', X_map(:,dd) = f_c( Zp * squeeze(hidden.B(d,:,1:hidden.R(d))) );
    #                              case 'o', X_map(:,dd) = f_o( Zp * squeeze(hidden.B(d,:,1))', hidden.theta(d,1:(hidden.R(d)-1)) );
    # otherwise
    # error('Unknown data type');
    # end
    # if (sum(isnan(X_map(:,dd))) > 0)
    #   warning('Some values are nan!');
    # end
    # if isfield(params,'t')
    # if ~isempty(params.t{d})
    # X_map(:,dd) = params.t{d}( X_map(:,dd) );
    # 
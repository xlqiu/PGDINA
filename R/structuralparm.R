#' Update structural parameters of the attribute distribution (internal)
#'
#' Internal helper called in the M-step to update the joint attribute (latent class)
#' distribution under different assumptions. Supported attribute distribution models
#' include saturated, loglinear smoothing, independent, fixed, and higher-order
#' (delegated to \code{HO.est()}).
#'
#' @param AlphaPattern A binary matrix of all attribute patterns, of size \eqn{2^K \times K}.
#' @param no.mg Integer; number of groups \eqn{G}.
#' @param logprior A matrix of log prior probabilities of size \eqn{2^K \times G}.
#' @param att.dist Character vector of length 1 or \eqn{G} specifying the attribute distribution
#'   model per group. Values used include \code{"saturated"}, \code{"loglinear"},
#'   \code{"independent"}, \code{"fixed"}, and \code{"higher.order"}.
#' @param att.str Logical; whether an attribute structure is imposed. When \code{TRUE},
#'   only \code{att.dist = "saturated"} or \code{"fixed"} is allowed in this implementation.
#' @param saturated A list containing saturated prior information. Must include
#'   \code{saturated$prior}, expected to be a vector or matrix of size \eqn{2^K \times G}
#'   used as pseudo-count prior (NEED TO UPDATE).
#' @param initial.logprior A matrix of initial log prior probabilities of size \eqn{2^K \times G}.
#'   Used when \code{att.dist = "fixed"}.
#' @param lower.prior Numeric scalar or length-\eqn{G} vector giving the minimum allowed class prior
#'   probability (applied per group, then renormalized).
#' @param loglinear Integer scalar or length-\eqn{G} vector specifying maximum interaction order for
#'   \code{att.dist = "loglinear"}. Default is 2.
#' @param N Integer; total sample size (used in loglinear smoothing as \code{log(N * prior)}).
#' @param Ng Integer vector of length \eqn{G}; group sample sizes.
#' @param K Integer; number of attributes.
#' @param lambda List (length \eqn{G}) of current structural parameters; updated and returned.
#'   For example, it stores class priors for \code{"saturated"}, coefficients for \code{"loglinear"},
#'   marginal mastery probabilities for \code{"independent"}, or higher-order parameters
#'   when \code{"higher.order"} (NEED TO UPDATE).
#' @param higher.order List of higher-order model settings passed to \code{HO.est()} when
#'   \code{att.dist = "higher.order"} (NEED TO UPDATE structure).
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{logprior}: Updated log prior matrix of size \eqn{2^K \times G}.
#'   \item \code{lambda}: Updated structural parameter list (length \eqn{G}).
#'   \item \code{higher.order}: Updated higher-order settings (or \code{NULL}).
#'   \item \code{higher.order.npar}: Number of higher-order parameters (NEED TO UPDATE; \code{NA} if not applicable).
#' }
#'
#' @examples
#' \dontrun{
#' # Internal function called during EM updates.
#' }
#'
#' @seealso \code{\link{Est.inOptim}}
#' @keywords internal
structural.parm <- function(AlphaPattern,no.mg,logprior,att.dist,att.str,saturated,initial.logprior,
                            lower.prior,loglinear=2,N,Ng,K,lambda,higher.order){
  prior <- exp(logprior)
  parm <- NULL
  ho <- NULL
  if(att.str){
    if (all(att.dist=="saturated")){
      for(g in seq_len(no.mg)){
          prior[,g] <- (prior[,g]*Ng[g]+saturated$prior[,g])/sum(prior[,g]*Ng[g]+saturated$prior[,g])
          prior[which(prior[,g]<lower.prior[g]),g] <- lower.prior[g]
          prior[,g] <- prior[,g]/sum(prior[,g])
          logprior[,g] <- log(prior[,g])
          lambda[[g]] <- prior[,g]

        }
    }else if(all(att.dist=="fixed")){
      for(g in seq_len(no.mg)){
        logprior[,g] <- initial.logprior[,g]
        lambda[[g]] <- NA
      }

    }else{
      stop("Only att.dist = saturated or fixed is available when att.str = TRUE.",call. = FALSE)
    }
  }else{
    if (all(tolower(att.dist)=="higher.order"))
    {
      HO.out <- HO.est(lambda=lambda,AlphaPattern = AlphaPattern, HOgr = seq_len(no.mg), Rl = rowProd(prior,Ng),
                       higher.order = higher.order)
      logprior <- HO.out$logprior
      lambda <- HO.out$lambda
      ho <- HO.out$higher.order
      npar <- HO.out$npar
    }else{
      for(g in seq_len(no.mg)){
        if (att.dist[g]=="saturated"){
          prior[,g] <- (prior[,g]*Ng[g]+saturated$prior[,g])/sum(prior[,g]*Ng[g]+saturated$prior[,g])
          prior[which(prior[,g]<lower.prior[g]),g] <- lower.prior[g]
          prior[,g] <- prior[,g]/sum(prior[,g])
          logprior[,g] <- log(prior[,g])
          lambda[[g]] <- prior[,g]
        }else if (att.dist[g]=="loglinear"){
          Z <- designM(K,0)
          if(K<2) stop("loglinear smoothing is not available when K < 2.",call. = FALSE)
          if(loglinear[g]>K) stop("Argument 'loglinear' cannot be greater than K.",call. = FALSE)
          Z <- Z[,seq_len(1+sum(sapply(seq_len(loglinear[g]),choose,n=K)))]
          prior[prior[,g]<1e-9,g] <- 1e-9
          lambda[[g]] <- parm <- stats::lm.wfit(x=Z,log(N*prior[,g]),w=N*prior[,g])$coefficients
          logprior[,g] <- c(Z%*%parm)-log(sum(exp(Z%*%parm)))
          logprior[logprior[,g]<log(.Machine$double.eps),g] <- log(.Machine$double.eps)
          logprior[logprior[,g]>log(1-.Machine$double.eps),g] <- log(1-.Machine$double.eps)
        }else if(att.dist[g]=="independent"){
          lambda[[g]] <- pk <- colSums(rowProd(prior,Ng)[,g]*AlphaPattern)/Ng[g] #length of K
          logprior[,g] <- AlphaPattern%*%log(pk) + (1-AlphaPattern)%*%log(1-pk)
        }else if(att.dist[g]=="fixed"){
          logprior[,g] <- initial.logprior[,g]
          lambda[[g]] <- NA
        }
      }
    }
  }




  return(list(logprior=logprior,lambda=lambda,higher.order = ho,higher.order.npar=npar))
}

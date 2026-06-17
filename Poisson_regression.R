
# Dual Likelihood
poisson_dual_likelihood <- function(Y, gram_mat, lambda, u){
  
  
  
  m <- length(Y)
  XXT <- as.matrix(gram_mat) 
  
  return(  (sum((u - Y)*((log(Y - u)) - 1))/m - (t(u) %*% XXT %*% u)/(2*m*m*lambda))[1,1] )
}


# Dual Gradient
poisson_dual_grad <- function(Y, gram_mat, lambda, u){
  
  m <- length(Y)
  
  XXT <- as.matrix(gram_mat)
  
  log_term <- log(Y - u)
  regularization <- (XXT %*% u)/(m*m*lambda)
  
  g <- (log_term)/m - regularization
  
  
  return( g )
  
}

# Dual Hessian Matrix
poisson_dual_Hess <- function(Y, gram_mat, lambda, u){

  m <- length(Y)
  XXT <- as.matrix(gram_mat) 
  
  main_h <- (1/(u - Y))/m
  diag_main_h <- matrix(0, nrow = m, ncol = m) 
  diag(diag_main_h) <- main_h
  
  H <- diag_main_h - XXT/(m*m*lambda)
  
  return( H )
}


# Projected gradient descent algorithm
poisson_dual_gradient_ascent <- function(Y, gram_mat, lambda, u0 = NULL,
                                         eps = 1e-8,
                                         step0 = 1.0,
                                         max_iter = 50000 ,
                                         tol = 1e-6,
                                         backtrack_beta = 0.5,
                                         armijo_c = 1e-4,
                                         verbose = TRUE){
  
  Y <- as.numeric(Y)
  # X <- as.matrix(X)
  XXT <- as.matrix(gram_mat)
  m <- length(Y)
  
  # if (length(Y) != m) stop("length(Y) must equal nrow(X).")
  # if (!is.finite(lambda) || lambda <= 0) stop("lambda must be > 0.")
  
  # projection to the feasible set: u <= Y - eps (elementwise)
  proj_u <- function(u) pmin(u, Y - eps)
  
  # initialize u
  if (is.null(u0)) {
    u <- Y - 1.0  # a safe default if Y >= 1; we'll project anyway
  } else {
    u <- as.numeric(u0)
    if (length(u) != m) stop("length(u0) must equal length(Y).")
  }
  u <- proj_u(u)
  
  # helper: safe objective (returns -Inf if out of domain)
  f_safe <- function(u){
    if (any(Y - u <= 0) || any(!is.finite(u))) return(-Inf)
    poisson_dual_likelihood(Y, XXT, lambda, u)
  }
  
  f <- f_safe(u)
  if (!is.finite(f)) stop("Initial u is infeasible (Y - u must be > 0).")
  
  history <- data.frame(iter = integer(0), obj = numeric(0), grad_norm = numeric(0), step = numeric(0))
  
  for (iter in 1:max_iter){
    g <- poisson_dual_grad(Y, XXT, lambda, u)
    if (any(!is.finite(g))) stop("Gradient produced non-finite values; check feasibility and inputs.")
    
    gnorm <- sqrt(sum(g^2))
    if (gnorm < tol) {
      if (verbose) message(sprintf("Converged at iter %d: grad_norm=%.3e, obj=%.10f", iter, gnorm, f))
      history <- rbind(history, data.frame(iter=iter, obj=f, grad_norm=gnorm, step=0))
      break
    }
    
    step <- step0
    # ascent direction is +g
    # Armijo-style condition for ascent:
    # f(u_new) >= f(u) + c * step * ||g||^2
    repeat {
      u_new <- proj_u(u + step * g)  # projection enforces Y - u_new >= eps
      f_new <- f_safe(u_new)
      
      if (is.finite(f_new) && (f_new >= f + armijo_c * step * (gnorm^2))) {
        break
      }
      
      step <- step * backtrack_beta
      if (step < 1e-16) {
        if (verbose) message(sprintf("Step size collapsed at iter %d. Stopping.", iter))
        u_new <- u
        f_new <- f
        break
      }
    }
    
    u <- u_new
    f <- f_new
    
    
    
  }
  
  list(u = u, objective = f)
}


coef_recover <- function(XTu, u, lambda){
  
  b = XTu / (nrow(u)  * lambda)
  
  return(b)
  
}

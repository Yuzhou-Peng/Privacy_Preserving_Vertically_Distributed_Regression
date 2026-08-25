# Negative binomial dual objective D(u)
D_obj <- function(u, gram_mat, y, r, m, lambda) {
  # domain check to avoid invalid logs
  if (any(r + u <= 0) || any(y - u <= 0)) return(-Inf)
  
  XXT = as.matrix(gram_mat)
  
  
  h_term <- (r + u) * log((r + u) / r) +
    (y - u) * log(y - u) -
    (r + y) * log(r + y)
  D_val <- mean(-h_term) - 0.5 / (m^2 * lambda) * (t(u) %*% XXT %*% u)[1,1]
  as.numeric(D_val)
}

# Gradient of D(u)
grad_D <- function(u, gram_mat, y, r, m, lambda) {
  
  XXT = as.matrix(gram_mat)
  
  grad_h <- log((r + u) / (r * (y - u)))   # d/du h(-u)
  grad_val <- -grad_h / m - (1 / (m^2 * lambda)) * (XXT %*% u)
  as.vector(grad_val)
}

# Hessian of D(u)
hess_D <- function(u, gram_mat, y, r, m, lambda) {
  
  XXT = as.matrix(gram_mat)
  
  diag_h <- - (1 / (r + u) + 1 / (y - u)) / m   # -h''(-u_i)/n
  H <- diag(diag_h) - (1 / (m^2 * lambda)) * (XXT)
  H
}


NB_dual_gradient_ascent <- function(u_init, gram_mat, y, r, m, lambda, 
                                    lr = 0.01, tol = 1e-6, max_iter = 1e6) {
  
  XXT <- as.matrix(gram_mat)
  u <- u_init
  for (iter in 1:max_iter) {
    g <- grad_D(u, XXT, y, r, m, lambda)
    u_new <- u + lr * g
    if (any(r + u_new <= 0) || any(y - u_new <= 0)) {
      warning("Step out of domain, reducing step size")
      lr <- lr / 2
      next
    }
    if (sqrt(sum((u_new - u)^2)) < tol) break
    u <- u_new
  }
  
  
  list(u = u, iter = iter, value = D_obj(u, XXT, y, r, m, lambda))
}


NB_dual_optim <- function(Y, gram_mat, r, lambda,
                          u_init = NULL,
                          eps = 1e-8,
                          maxit = 10000) {
  
  Y <- as.numeric(Y)
  m <- length(Y)
  XXT <- as.matrix(gram_mat)
  
  if (nrow(XXT) != m || ncol(XXT) != m) {
    stop("gram_mat must be an m by m matrix where m = length(y).")
  }
  
  if (!is.finite(r) || r <= 0) {
    stop("r must be positive.")
  }
  
  if (!is.finite(lambda) || lambda <= 0) {
    stop("lambda must be positive.")
  }
  
  # Feasible domain:
  # r + u > 0  =>  u > -r
  # y - u > 0  =>  u < y
  lower_bound <- rep(-r + eps, m)
  upper_bound <- Y - eps
  
  if (is.null(u_init)) {
    u_init <- pmin(pmax(rep(0, m), lower_bound), upper_bound)
  } else {
    u_init <- as.numeric(u_init)
    if (length(u_init) != m) stop("u0 must have length equal to length(y).")
    u_init <- pmin(pmax(u_init, lower_bound), upper_bound)
  }
  
  neg_obj <- function(u) {
    -D_obj(
      u = u,
      gram_mat = XXT,
      y = Y,
      r = r,
      m = m,
      lambda = lambda
    )
  }
  
  neg_grad <- function(u) {
    -grad_D(
      u = u,
      gram_mat = XXT,
      y = Y,
      r = r,
      m = m,
      lambda = lambda
    )
  }
  
  fit <- optim(
    par = u_init,
    fn = neg_obj,
    gr = neg_grad,
    method = "L-BFGS-B",
    lower = lower_bound,
    upper = upper_bound,
    control = list(
      maxit = maxit,
      factr = 1e7
    )
  )
  
  list(
    u = fit$par,
    objective = -fit$value,
    convergence = fit$convergence,
    message = fit$message
  )
}


coef_recover <- function(XTu, u, lambda){
  
  b = XTu / (length(u)  * lambda)
  
  return(b)
  
}

NB_r_MoM <- function(Y,
                     gram_mat,
                     gamma,
                     lambda) {
  
  m <- length(Y)
  
  eta <- as.vector(gram_mat %*% gamma) /
    (m * lambda)
  
  mu <- exp(eta)
  
  numerator   <- sum(mu^2)
  denominator <- sum((Y - mu)^2 - mu)
  
  if (denominator <= 0) {
    warning("MoM estimate is not positive.")
    return(Inf)
  }
  
  numerator / denominator
}

NB_loglik_r <- function(Y, gram_mat, u, lambda, r) {
  m <- length(Y)
  
  eta <- as.vector(gram_mat %*% u) / (m * lambda)
  mu  <- exp(eta)
  
  sum(
    lgamma(Y + r) -
      lgamma(r) +
      r * log(r) +
      Y * eta -
      (Y + r) * log(r + mu)
  )
}


NB_r_optim <- function(
    Y,
    gram_mat,
    u,
    lambda,
    lower = 1e-6,
    upper = 1e5
) {
  # Dual feasibility requires u_i > -r
  lower <- max(lower, -min(u) + 1e-8)
  
  objective <- function(r) {
    -NB_loglik_r(
      Y = Y,
      gram_mat = gram_mat,
      u = u,
      lambda = lambda,
      r = r
    )
  }
  
  fit <- optimize(
    f = objective,
    interval = c(lower, upper)
  )
  
  fit$minimum
}

NB_alternating_optim <- function(
    Y,
    gram_mat,
    lambda,
    r_init,
    u_init = rep(0, length(Y)),
    max_iter = 100,
    tol = 1e-6,
    verbose = FALSE
) {
  Y <- as.numeric(Y)
  u_old <- as.numeric(u_init)
  r_old <- r_init
  
  history <- data.frame(
    iteration = integer(0),
    r = numeric(0),
    u_change = numeric(0),
    r_change = numeric(0)
  )
  
  converged <- FALSE
  
  for (iter in seq_len(max_iter)) {
    
    # Optimize u given the current r
    dual_fit <- NB_dual_optim(
      Y = Y,
      gram_mat = gram_mat,
      lambda = lambda,
      r = r_old,
      u_init = u_old
    )
    
    u_new <- dual_fit$u
    
    # Optimize r given the updated u
    r_new <- NB_r_optim(
      Y = Y,
      gram_mat = gram_mat,
      u = u_new,
      lambda = lambda
    )
    
    u_change <- sqrt(sum((u_new - u_old)^2)) /
      (1 + sqrt(sum(u_old^2)))
    
    r_change <- abs(r_new - r_old) /
      (1 + abs(r_old))
    
    history <- rbind(
      history,
      data.frame(
        iteration = iter,
        r = r_new,
        u_change = u_change,
        r_change = r_change
      )
    )
    
    if (verbose) {
      cat(
        sprintf(
          "Iteration %d: r = %.6f, u change = %.3e, r change = %.3e\n",
          iter,
          r_new,
          u_change,
          r_change
        )
      )
    }
    
    u_old <- u_new
    r_old <- r_new
    
    if (u_change < tol && r_change < tol) {
      converged <- TRUE
      break
    }
  }
  
  list(
    u = u_old,
    r = r_old,
    converged = converged,
    history = history
  )
}


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


# Projected Gradient Algorithm

# Helper function: project onto the feasible box: -r < u < y

project_u <- function(u, y, r, eps = 1e-8) {
  pmin(pmax(u, -r + eps), y - eps)
}

# Projected gradient ascent for maximizing D(u)
proj_grad_ascent_nb <- function(
    u_init,
    gram_mat,
    y,
    r,
    m,
    lambda,
    step_init = 1,
    beta = 0.5,
    sigma = 1e-4,
    max_iter = 5000,
    tol = 1e-6,
    eps = 1e-8,
    verbose = TRUE
) {
  # ensure vectors
  u <- as.vector(u_init)
  y <- as.vector(y)
  r <- as.vector(r)
  XXT <- as.matrix(gram_mat)
  
  
  

  # start from feasible point
  u <- project_u(u, y = y, r = r, eps = eps)
  
  obj_hist <- numeric(max_iter)
  
  for (iter in 1:max_iter) {
    obj_old <- D_obj(u, XXT, y, r, m, lambda)
    g <- grad_D(u, XXT, y, r, m, lambda)
    
    # if gradient has non-finite values, stop
    if (any(!is.finite(g))) {
      warning("Non-finite gradient encountered.")
      break
    }
    
    step <- step_init
    
    # projected gradient step
    repeat {
      u_trial <- project_u(u + step * g, y = y, r = r, eps = eps)
      obj_trial <- D_obj(u_trial, XXT, y = y, r = r, m = m, lambda)
      
      # Armijo-type ascent condition for projected gradient
      step_dir <- u_trial - u
      rhs <- obj_old + sigma * sum(g * step_dir)
      
      if (is.finite(obj_trial) && obj_trial >= rhs) {
        break
      }
      
      step <- beta * step
      
      if (step < 1e-16) {
        if (verbose) {
          message("Line search step became too small at iteration ", iter)
        }
        break
      }
    }
    
    obj_hist[iter] <- obj_trial
    
    # convergence check: projected step is tiny
    if (sqrt(sum((u_trial - u)^2)) < tol) {
      u <- u_trial
      if (verbose) {
        message("Converged at iteration ", iter)
      }
      obj_hist <- obj_hist[1:iter]
      return(list(
        u_opt = u,
        obj_value = D_obj(u, gram_mat, y, r, m, lambda),
        gradient = grad_D(u, gram_mat, y, r, m, lambda),
        iterations = iter,
        converged = TRUE,
        obj_history = obj_hist
      ))
    }
    
    u <- u_trial
    
    if (verbose && iter %% 100 == 0) {
      message("iter = ", iter,
              ", obj = ", round(obj_hist[iter], 8),
              ", step = ", signif(step, 3))
    }
  }
  
  obj_hist <- obj_hist[1:max_iter]
  
  list(
    u_opt = u,
    obj_value = D_obj(u, gram_mat, y, r, m, lambda),
    gradient = grad_D(u, gram_mat, y, r, m, lambda),
    iterations = max_iter,
    converged = FALSE,
    obj_history = obj_hist
  )
}


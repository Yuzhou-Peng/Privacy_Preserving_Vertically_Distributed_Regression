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



coef_recover <- function(XTu, u, lambda){
  
  b = XTu / (nrow(u)  * lambda)
  
  return(b)
  
}
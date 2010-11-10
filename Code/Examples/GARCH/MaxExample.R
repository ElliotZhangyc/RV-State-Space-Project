# An example script to find the maximum of a function.

# Some function.
f <- function(theta){
  theta %*% theta;
}

# The initial values to start our search.
iv = c( 0.5, 0.2);

# Let's see what happens.  The NULL parameter here is for calculating
# the gradient.  This is something I can actually do analytically for
# the GARCH problem.
output = optim(iv, f, NULL, method="L-BFGS-B",
               lower=c(-1,-1), upper=c(1,1), hessian=TRUE);

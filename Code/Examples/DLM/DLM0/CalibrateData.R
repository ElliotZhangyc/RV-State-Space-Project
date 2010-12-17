# We want to generate the posterior distribution through
# simulation for our x values and our parameter values.  We will use a
# Gibbs sampler.  Thus we break our draw from
#   p(phi, W, mu, x[0:N] | y[1:N])
# into steps
#   p(phi | everything else)
#   p(W | everything else)
#   p(mu | everything else)
#   p(x[0:n] | everything else).
# We chose conjugate priors for all of our parameters and a diffuse
# prior for x_0.  The derivation of the above priors can be found in
# my notes.

## ASSUMPTIONS ##

# We assume our data is stored in y.data

# The length of our time series.
N = length(y.data);

# The number of samples we take.
num.samples = 10000;

## DATA STRUCTURES ##

# Now set up our data structures.  We set the values to the ``true
# values'' so that we may check out our Gibbs steps.
phi = rep(phi.true, num.samples);
W = rep(W.true, num.samples);
mu = rep(mu.true, num.samples);
x = matrix(0, num.samples, N+1);

## PRIOR PARAMETERS ##

# Set up the prior parameters.  We need this information to draw from
# our conditional posterior distributions.

# We assume phi ~ N(m.phi, C.phi).
m.phi = 0.3;
C.phi = 0.5;

# We assume W ~ Inv-Gamma(a/2, b/2).
W.a = 2.0;
W.b = 2.0;

# We assume mu ~ N(m.mu, C.mu).
m.mu = 0.5;
C.mu = 1.0;

# We assume x_0 ~ N(m.0, C.0).
m.0 = 0.0;
C.0 = 5.0;

## SEED VALUES ##

# Let's seed our sampler with some values.
phi[1] = phi.true;
W[1] = W.true;
mu[1] = mu.true;
x[1,1] = 0.2;
# We need to generate the rest of the x values.
for(i in 2:(N+1)){
  x[1,i] = phi[1]*x[1,i-1] + rnorm(1, 0, sqrt(W[1]));
}

## TEMPORARY CALCULATIONS ##

for(i in 1:num.samples){
  x[i,] = x.data;
}

## CONDITIONAL DENSITIES ##

# Let's set up the function, which will create the draws from our
# conditional posteriors.

# To draw from x | everything else.
x.cond.post <- function(y.data, mu, V, phi, W, m.0, C.0){
  # We need to keep track of some information as we feed forward.  The
  # procedure we use here can be found in my notes.  We need an extra
  # spot for m.0 and C.0 in our m an C arrays.  Since R forces us to
  # start indexing at 1, we have that m[1] = m.0 and C[1] = C.0.  We
  # must adjust our indexing accordingly.
  n = length(y.data);
  R = rep(0, n);
  m = rep(0, n+1); m[1] = m.0;
  C = rep(0, n+1); C[1] = C.0;
  # Now let's feed forward.
  for(i in 1:n){
    R[i] = phi^2 * C[i] + W;
    Q = R[i] + V;
    A = R[i]/Q;
    m[i+1] = phi*m[i] + A*(y.data[i] - mu - phi*m[i]);
    C[i+1] = A*V;
  }
  # Now let's go backwards.
  # We need a data structure to hold our draw from x.
  x.draw = rep(0, n+1);
  # First draw x_n | y.data.
  x.draw[n+1] = rnorm(1, m[n+1], sqrt(C[n+1]) );
  # Now we draw backwards, conditionally on our previous values of x.
  for(i in n:1){
    A = phi*C[i]/R[i];
    m.x.draw = m[i] + A * (x.draw[i+1] - phi*m[i]);
    C.x.draw  = C[i]*W/R[i];
    # Now draw x_{i-1}, which is x[i] since we have an x_0 value.
    x.draw[i] = rnorm(1, m.x.draw, sqrt(C.x.draw) );
  }
  # Return x.draw.
  x.draw;
  # Return parameters.
  # matrix( c(m,C), N+1, 2 );
}

# To draw from phi | everything else.
# Remeber that x.data[1:(N+1)] = {x_0, \ldots, x_n}.
phi.cond.post <- function(x.data, W, m.phi, C.phi){
  n = length(x.data) - 1;
  g.0 = x.data[1:n] %*% x.data[1:n];
  g.1 = x.data[2:(n+1)] %*% x.data[1:n];
  # Some useful quantities.  Refer to notes.
  Q = g.0^2 * C.phi + g.0 * W;
  R = g.0 * C.phi;
  A = R/Q;
  m.phi.draw = m.phi + A * (g.1 - g.0 * m.phi);
  C.phi.draw = A*W;
  phi.draw = rnorm(1, m.phi.draw, sqrt(C.phi.draw) );
}

# To draw from mu | everything else.
mu.cond.post <- function(y.data, x.data, V, m.mu, C.mu){
  # Refer to notes for specifics.
  n = length(y.data);
  # We need the quantity a from the notes.
  a = mean(y.data-x.data[2:(n+1)]);
  Q = n*C.mu + V;
  m.mu.draw = V*m.mu/Q + n*C.mu*a/Q;
  C.mu.draw = C.mu*V/Q;
  mu.draw = rnorm(1, m.mu.draw, sqrt(C.mu.draw));
}

# To draw from W | everything else.
W.cond.post <- function(x.data, phi, W.a, W.b){
  # Again, review notes for details.
  n = length(x.data)-1;
  a.post = W.a + n;
  diff = x.data[2:(n+1)] - phi*x.data[1:n];
  b.post = W.b + diff %*% diff;
  W.recip.draw = rgamma(1, a.post/2, rate=b.post/2);
  W.draw = 1/W.recip.draw;
}

## GIBBS SAMPLING ##

# Now we can go ahead and do our Gibbs sampling.
# Again refer to the notes for a description.
for(i in 2:num.samples){
  phi[i] = phi.cond.post(x[i-1,], W[i-1], m.phi, C.phi);
  W[i] = W.cond.post(x[i-1,], phi[i], W.a, W.b);
  mu[i] = mu.cond.post(y.data, x[i-1,], V, m.mu, C.mu);
  x[i,] = x.cond.post(y.data, mu[i], V, phi[i], W[i], m.0, C.0);
}

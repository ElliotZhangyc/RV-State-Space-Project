# We are considering a dynamics linear model governed by
#  y_t = z_t + \nu_t, \nu_t \sim N(m.nu ,V)
#  z_t = \mu + \phi (z_{t-1} - \mu) + \omega_t, \omega_t \sim N(0,W).

# We want to generate the posterior distribution through
# simulation for our x values and our parameter values.  We will use a
# Gibbs sampler.  Thus we break our draw from
#   p(phi, W, mu, x[0:N] | y[1:N])
# into steps
#   p(phi | everything else)
#   p(W | everything else)
#   p(mu | everything else)
#   p(z[0:n] | everything else).
# We chose conjugate priors for all of our parameters and a diffuse
# prior for z_0.  The derivation of the above priors can be found in
# my notes.

## ASSUMPTIONS ##

# We assume our data is stored in y.data

# The length of our time series.
N = length(y.data);

# The number of samples we take.
num.samples = 3000;

# For the innovations nu we need to set the mean and variance.
m.nu = m.nu.true;
V = V.true;

## DATA STRUCTURES ##

# Now set up our data structures.  We set the values to the ``true
# values'' so that we may check out our Gibbs steps.
phi = rep(phi.true, num.samples);
W = rep(W.true, num.samples);
mu = rep(mu.true, num.samples);
z = matrix(0, num.samples, N+1);

## PRIOR PARAMETERS ##

# Set up the prior parameters.  We need this information to draw from
# our conditional posterior distributions.

# We assume phi ~ N(m.phi, C.phi).
m.phi = 0.5;
C.phi = 0.5;

# We assume W ~ Inv-Gamma(a/2, b/2).
W.a = 2.0;
W.b = 2.0;

# We assume mu ~ N(m.mu, C.mu).
m.mu = 0.0;
C.mu = 1.0;

# We assume z_0 ~ N(m.0, C.0).
m.0 = 0.0;
C.0 = 5.0;

## SEED VALUES ##

# Let's seed our sampler with some values.
# phi[1] = phi.true;
phi[1] = 0.5;
# W[1] = W.true;
W[1] = 0.6;
# mu[1] = mu.true;
mu[1] = 0;
# z[1,1] = 0.2;
z[1,1] = 0.2;
# We need to generate the rest of the z values.
for(i in 2:(N+1)){
  z[1,i] = mu[1] + phi[1]*(z[1,i-1]-mu[1]) + rnorm(1, 0, sqrt(W[1]));
}

## TEMPORARY CALCULATIONS ##

# Useful if we are fixing z.
#for(i in 1:num.samples){
#  z[i,] = z.data;
#}

## CONDITIONAL DENSITIES ##

# Let's set up the function, which will create the draws from our
# conditional posteriors.

# To draw from x | everything else.
z.cond.post <- function(y.data, mu, m.nu, V, phi, W, m.0, C.0){
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
    # NOTE: We include m.mu here to accomodate nu with nonzero mean.
    m[i+1] = mu + phi*(m[i]-mu) + A*( y.data[i] - (mu + phi*(m[i]-mu) + m.nu) );
    C[i+1] = A*V;
  }
  # Now let's go backwards.
  # We need a data structure to hold our draw from z.
  z.draw = rep(0, n+1);
  # First draw z_n | y.data.
  z.draw[n+1] = rnorm(1, m[n+1], sqrt(C[n+1]) );
  # Now we draw backwards, conditionally on our previous values of z.
  for(i in n:1){
    A = phi*C[i]/R[i];
    m.z.draw = m[i] + A * (z.draw[i+1] - (mu + phi*(m[i]-mu)) );
    C.z.draw  = C[i]*W/R[i];
    # Now draw z_{i-1}, which is z[i] since we have an z_0 value.
    z.draw[i] = rnorm(1, m.z.draw, sqrt(C.z.draw) );
  }
  # Return z.draw.
  z.draw;
  # Return parameters.
  # matrix( c(m,C), N+1, 2 );
}

# To draw from phi | everything else.
# Remeber that z.data[1:(N+1)] = {z_0, \ldots, z_n}.
# x.data = z.data - mu.
phi.cond.post <- function(x.data, W, m.phi, C.phi){
  n = length(x.data) - 1;
  g.0 = x.data[1:n] %*% x.data[1:n];
  g.1 = x.data[2:(n+1)] %*% x.data[1:n];
  # Some useful quantities.  Refer to notes.
  Q = C.phi + W / g.0;
  m.phi.draw = m.phi * (W/g.0) / Q + (g.1/g.0) * C.phi / Q;
  C.phi.draw = C.phi * (W/g.0) / Q;
  phi.draw = rnorm(1, m.phi.draw, sqrt(C.phi.draw) );
}

# To draw from mu | everything else.
mu.cond.post <- function(z.data, phi, W, m.mu, C.mu){
  # Refer to notes for specifics.
  n = length(z.data)-1;
  # We need the quantity a from the notes.
  # In particular, look at the notes for ``Model 2.''
  a = mean(z.data[2:(n+1)]- phi*z.data[1:n]);
  Q = (1-phi)^2 * n * C.mu + W;
  m.mu.draw = m.mu * W / Q + a * (1-phi) * n * C.mu / Q;
  # print(m.mu.draw);
  C.mu.draw = C.mu * W / Q;
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
  # A change of variables is useful.
  x.data = z[i-1,] - mu[i-1];
  phi[i] = phi.cond.post(x.data, W[i-1], m.phi, C.phi);
  W[i] = W.cond.post(x.data, phi[i], W.a, W.b);
  mu[i] = mu.cond.post(z[i-1,], phi[i], W[i], m.mu, C.mu);
  z[i,] = z.cond.post(y.data, mu[i], m.nu, V, phi[i], W[i], m.0, C.0);
  if (i %% 100 == 0) print(c(i,date()));
}

## CHECKING OUTPUT ##

# For normal posterior distributions.  QoI = quantity of interest.
#range = 3000:num.samples;
#QoI = mu[range];
#m.QoI = mean(QoI);
#sd.QoI = sd(QoI);
#hist(QoI, breaks=40, prob=TRUE);
#grid = seq(min(QoI), max(QoI), 0.01);
#dens = dnorm(grid, m.QoI, sd.QoI);
#lines(grid, dens);
#print(c(m.QoI, sd.QoI));

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

# For the innovations nu we need to set the mean and variance.
m.nu = nu$m;
V = nu$v;

## DATA STRUCTURES ##

# Now set up our data structures.  We set the values to the ``true
# values'' so that we may check out our Gibbs steps.
phi.gibbs = rep(true$phi, mcmc$samples);
W.gibbs = rep(true$W, mcmc$samples);
mu.gibbs = rep(true$mu, mcmc$samples);
z.gibbs = matrix(0, N+1, mcmc$samples);
gamma.gibbs = rep(0, mcmc$samples);

## SEED THE ARRAYS ##

# Let's seed our sampler with some values.

## SEED THE ARRAYS ##

# For mu
mu.gibbs[1] = seed$mu;
# For phi
phi.gibbs[1] = seed$phi;
# For W
W.gibbs[1] = seed$W;
# FOr gamma
gamma.gibbs[1] = seed$gamma;

# F or z.
z.gibbs[1,1] = seed$z.0;
for(i in 2:(T+1)){
  z.gibbs[i,1] = mu.gibbs[1] + phi.gibbs[1]*(z.gibbs[i-1,1]-mu.gibbs[1]) +
    gamma.gibbs[1]*ell[i-1] + rnorm(1, 0, sqrt(W.gibbs[1]));
}

## TEMPORARY CALCULATIONS ##

#mu.gibbs = rep(true$mu, mcmc$samples);
phi.gibbs = rep(true$phi, mcmc$samples);
W.gibbs = rep(true$phi, mcmc$samples);
#gamma.gibbs = rep(true$gamma, mcmc$samples);

# For q
#for(i in 1:mcmc$samples){
#  q.gibbs[,i] = q.data;
#}

for(i in 1:mcmc$samples){
  z.gibbs[,i] = z.data;
}

## CONDITIONAL DENSITIES ##

# Let's set up the function, which will create the draws from our
# conditional posteriors.

# To draw from x | everything else.
z.cond.post <- function(y.data, mu, m.nu, V, phi, W, gamma, m.0, C.0){
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
    # m.z.tprev = E[z_t | D_{t-1}].  In the notes we know that E[z_t |
    # D_{t-1}] = mu + phi*(m_{t-1} - mu) + gamma l_t.  Since m_{t-1} =
    # m[t] we have the following.
    m.z.tprev = mu + phi*(m[i] - mu) + gamma * ell[i];
    m[i+1] = m.z.tprev +
             A*( y.data[i] - (m.z.tprev + m.nu ));
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
    m.z.draw = m[i] + A * (z.draw[i+1] -
      (mu + phi*(m[i]-mu) + gamma*ell[i]) );
    C.z.draw  = C[i]*W/R[i];
    # Now draw z_{i-1}, which is z[i] since we have an z_0 value.
    z.draw[i] = rnorm(1, m.z.draw, sqrt(C.z.draw) );
  }
  # Return z.draw.
  z.draw;
}

# To draw from phi | everything else.
# Remeber that z.data[1:(N+1)] = {z_0, \ldots, z_n}.
# x.data = z.data - mu.
phi.cond.post <- function(x.data, W, m.phi, C.phi, gamma){
  n = length(x.data) - 1;
  g.0 = x.data[1:n] %*% x.data[1:n];
  g.1 = (x.data[2:(n+1)]-gamma*ell[1:n]) %*% x.data[1:n];
  # Some useful quantities.  Refer to notes.
  Q = C.phi + W / g.0;
  m.phi.draw = m.phi * (W/g.0) / Q + (g.1/g.0) * C.phi / Q;
  C.phi.draw = C.phi * (W/g.0) / Q;
  phi.draw = rnorm(1, m.phi.draw, sqrt(C.phi.draw) );
}

# To draw from mu | everything else.
mu.cond.post <- function(z.data, phi, W, m.mu, C.mu, gamma){
  # Refer to notes for specifics.
  n = length(z.data)-1;
  # We need the quantity a from the notes.
  # In particular, look at the notes for ``Model 2.''
  a = mean(z.data[2:(n+1)]- phi*z.data[1:n] - gamma*ell[1:n]);
  Q = (1-phi)^2 * n * C.mu + W;
  m.mu.draw = m.mu * W / Q + a * (1-phi) * n * C.mu / Q;
  # print(m.mu.draw);
  C.mu.draw = C.mu * W / Q;
  mu.draw = rnorm(1, m.mu.draw, sqrt(C.mu.draw));
}

# To draw from W | everything else.
W.cond.post <- function(x.data, phi, W.a, W.b, gamma){
  # Again, review notes for details.
  n = length(x.data)-1;
  a.post = W.a + n;
  diff = x.data[2:(n+1)] - phi*x.data[1:n] - gamma*ell[1:n];
  b.post = W.b + diff %*% diff;
  W.recip.draw = rgamma(1, a.post/2, rate=b.post/2);
  W.draw = 1/W.recip.draw;
}

# To draw from gamma | everything else
gamma.cond.post <- function(x.data, W, phi, m.gamma, C.gamma){
  n = length(x.data) - 1;
  # g.0 DOESN'T CHANGE as we iterate.
  g.0 = ell[1:n] %*% ell[1:n];
  g.1 = (x.data[2:(n+1)] - phi*x.data[1:n]) %*% ell[1:n];
  # Some useful quantities.  Refer to notes.
  # These values correspond to our ``intermediate'' density.
  Q = C.gamma + W / g.0;
  m.inter = m.gamma * (W/g.0) / Q + (g.1/g.0) * C.gamma / Q;
  C.inter = C.gamma * (W/g.0) / Q;
  # Draw from the posterior distribution.
  gamma.draw = rnorm(1, m.inter, sqrt(C.inter));
}
  
## GIBBS SAMPLING ##

# To keep track of the time.
the.time = rep(0, floor(mcmc$samples/100)+1);
delta.time = rep(0, length(the.time)-1);
the.time[1] = proc.time()[[1]];
  
# Now we can go ahead and do our Gibbs sampling.
# Again refer to the notes for a description.
for(i in 2:mcmc$samples){
  # A change of variables is useful.  Make sure to sample everything
  # else before you sample z and mu so that your x.data values are not
  # old.
  x.data = z.gibbs[,i-1] - mu.gibbs[i-1];
  # Sample phi | e.e.
#  phi.gibbs[i] = phi.cond.post(x.data, W.gibbs[i-1],
#             prior$m.phi, prior$v.phi, gamma.gibbs[i-1]);
  # Sample W | e.e.
#  W.gibbs[i] = W.cond.post(x.data, phi.gibbs[i],
#           prior$a.W, prior$b.W, gamma.gibbs[i-1]);
  # We want to make sure that we sample gamma before we sample mu.
  # This is necessary only because we are using x.data, which is based
  # on mu.gibbs[i-1].
  # Sample gamma | e.e.
  gamma.gibbs[i] = gamma.cond.post(x.data, W.gibbs[i], phi.gibbs[i],
               prior$m.gamma, prior$v.gamma);
  # Sample mu | e.e.
  mu.gibbs[i] = mu.cond.post(z.gibbs[,i-1], phi.gibbs[i], W.gibbs[i],
            prior$m.mu, prior$v.mu, gamma.gibbs[i]);
  # Sample z | e.e.
#  z.gibbs[,i] = z.cond.post(y.data, mu.gibbs[i], m.nu, V, phi.gibbs[i],
#           W.gibbs[i], gamma.gibbs[i], prior$m.0, prior$v.0);
  # Record the time if this is our 100th draw and output time estimates.
  if (i %% 100 == 0){
    idx = i/100 + 1; # We must increment by 1 because our array starts at 0.
    the.time[idx] = proc.time()[[1]];
    delta.time[idx-1] = the.time[idx] - the.time[idx-1];
    ave.time = mean(delta.time[1:(idx-1)]);
    min.to.complete = (mcmc$samples - i)/6000*ave.time;
    print(paste("Iter:", i, "|",
                round(delta.time[idx-1], 3), "s last 100", "|",
                round(ave.time, 3), "s per 100", "|",
                round(min.to.complete, 3), "min to done"));
  }
}

## RECORD AVERAGE TIME ##
mcmc$ave.time = ave.time;

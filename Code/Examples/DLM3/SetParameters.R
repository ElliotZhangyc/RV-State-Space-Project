# We set up the parameters of our system here for both the creation of
# synthetic data and for the estimation process.

## FOR SYNTHETIC DATA ##

# The length of our observed time series, (y_t).
T = 1000;

# Let's remember we are using synthetic data or not.
is.synthetic = TRUE;

# trial.id
trial.id = 1;

# The true values of our parameters.
# This is used when creating synthetic data.
true = data.frame(
  mu = -3.80,
  m.nu = -0.62,
  v.nu = 1.16,
  phi = 0.8,
  W = 0.5,
  z.0 = -3
);

## PRIOR PARAMETERS ##

# We need this information to draw from our conditional posterior
# distributions.

prior = data.frame(
  # We assume mu ~ N(m.mu, v.mu).
  m.mu = -3.50,
  v.mu = 4.0,
  # We assume phi ~ N(m.phi, v.phi).
  m.phi = 0.5,
  v.phi = 0.5,
  # We assume W ~ Inv-Gamma(a.W/2, b.W/2).
  a.W = 6.0,
  b.W = 2.0
);

# The values of nu, which have taken to be known.
nu = data.frame(
  m = true$m.nu,
  v = true$v.nu
);

## SEED VALUES ##

seed = data.frame(
  mu = -3.5,
  phi = 0.9,
  W = 0.5,
  z.0 = -3.5
);

## FOR THE MCMC ##

# Information for our MCMC.
mcmc = data.frame(
  model = 3,
  samples = 3000,
  burn.in = 1000,
  ave.time = 0 # ave.time will be set after we run the MCMC.
);

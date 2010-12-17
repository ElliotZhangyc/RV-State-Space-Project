# We set up the parameters of our system here for both the creation of
# synthetic data and for the estimation process for DLM3.

## FOR SYNTHETIC DATA ##

# The length of our observed time series, (y_t).
T = 1000;

# Let's remember we are using synthetic data or not.
run.info = data.frame(
  id = 18,
  is.synthetic = TRUE,
  time.stamp = NA
)

# The true values of our parameters.
# This is used when creating synthetic data.
true = data.frame(
  m.nu = -0.62,
  v.nu = 0.5,
  size = T,
  mu = -3.80,
  phi = 0.2,
  W = 0.2,
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
  a.W = 4.0,
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
  phi = 0.5,
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

## SYNTHETIC DATA STATISTICS ##

# We will set these values when we run Check.R.
synth.stats = data.frame(
  z.sm = NA,
  x.sv = NA,
  x.sacv = NA,
  x.sacr = NA,
  W.marg.m = NA
);

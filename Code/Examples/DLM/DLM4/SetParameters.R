# We set up the parameters of our system here for both the creation of
# synthetic data and for the estimation process for DLM4.

## FOR TRUE VALUES / SYNTHETIC DATA PARAMS ##

# The length of our observed time series, (y_t).
T = 1000;

# Let's remember we are using synthetic data or not.
run.info = data.frame(
  id = 24,
  is.synthetic = FALSE,
  time.stamp = NA
)

# The true values of our parameters.
# This is used when creating synthetic data.
#true = data.frame(
#  m.nu = NA,
#  v.nu = NA,
#  size = T,
#  mu = -3.80,
#  phi = 0.2,
#  W = 0.2,
#  z.0 = -3
#);

##  FOR REAL DATA ##
true = data.frame(
  m.nu = NA,
  v.nu = NA,
  size = NA,
  mu = NA,
  phi = NA,
  W = NA,
  z.0 = NA
);

## NORMAL MIXTURE ##

nmix = data.frame(
  q = c(0.0073, 0.1056, 0.2575, 0.34, 0.2456, 0.0440),
  b = c(-5.7002, -2.6216, -1.1793, -0.3255, 0.2624, 0.7537),
  v = c(1.4490, 0.6534, 0.3157, 0.16, 0.0851, 0.0418)
)

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
  model = 4,
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

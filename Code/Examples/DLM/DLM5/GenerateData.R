# We are considering a dynamics linear model governed by
#  y_t = z_t + \nu_t, \nu_t \sim N(m.nu ,V)
#  z_t = \mu + \phi (z_{t-1}-\mu) + \omega_t, \omega_t \sim N(0,W).

# Our data structures, y and x.  For plotting purposes we keep track
# of x.  y is indexed from 1 to n.  x is indexed from 0 to n.
y.true = rep(0,T);
z.true = rep(0,T+1);

# Set x_0 to something.  Since x_0 mean reverts to zero we chose
# something that is not too far away from 0.
z.true[1] = true$z.0;

# Now generate our data.
for (i in 1:T){
  z.true[i+1] = true$mu + true$phi*(z.true[i]-true$mu) +
    true$gamma * ell[i] +  rnorm(1, 0, sqrt(true$W));
  y.true[i] = z.true[i+1] + rnorm(1, true$m.nu, sqrt(true$v.nu));
}

# For purposes of calibrating.  Our CalibrateData.R script assumes
# there is a variable y.data with the observed y-values.
z.data = z.true;
y.data = y.true;

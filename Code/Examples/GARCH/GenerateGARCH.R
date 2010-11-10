# This script generates a time series from GARCH.
# Recall that a GARCH(1,1) model is
#   R_t = \sig_t \ep_t   (\ep_t ~ i.i.d. Normal or T or whatever)
# where \sig_t is the conditional variance and follows
#   \sig_t^2 = \omega + \alpha R_{t-1}^2 + \beta \sig_{t-1}^2.

# Let the user know what we are doing.
print("Generating GARCH(1,1) time series...");

# The length of the time series.
T = 200;

# Set aside space for the data.
# The returns:
R = rep(0, T);
# The conditional variance:
CV = rep(0,T);

# Seed the conditional variance with some value.
CV[1] = 0.3;

# We define a function for our innovations.  This way we can change it
# easily later.
innovation <- function(){
  rnorm(1, 0, 1);
}

# Set the parameters for our GARCH equation:
# CV[t] = omega + alpha * r[t-1]^2 + beta * CV[t-1].
omega.true = 0.1;
alpha.true = 0.2;
beta.true = 0.6;

# Generate the first return.
R[1] = sqrt( CV[1] ) * innovation();

# Now generate our GARCH(1,1) process.
for(i in 2:T){
  CV[i] = omega.true + alpha.true * R[i-1]^2 + beta.true * CV[i-1];
  R[i] = sqrt( CV[i] ) * innovation();
}

# To view the time series.
plot(1:T, R, type="l");

# Our time series should look like white noise.
# To calculate the autocorrelation:
# acf(R, lag.max=40);

# For use by other scripts.
r.data = R;
cv.data = CV;

# To write a plot.
# jpeg("filename.jpg");
# plot(1:T, R, type="l");
# title("Whatever");
# dev.off();

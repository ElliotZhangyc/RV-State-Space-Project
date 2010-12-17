## LOAD SP DATA ##
prices.df = read.table("SPXD.dat");
prices = prices.df$V2;
T = length(prices);
n = T-1;
R = log(prices[2:T])-log(prices[1:n]);
y = pmax(-20, 0.5*log(R^2));
y.data = y;

## LETS SIMULATE THE STATIONARY DISTRIBUTION ##
N1 = rnorm(100000, -4.6, sqrt(10));
N2 = rnorm(100000, 0, 1);
model.R = exp(N1)*N2;
trunc.model.R = pmin(0.1, pmax(-0.1,model.R));

## NOW PLOT ##

## OUTPUT FILE ##
#postscript(file="SP_Returns_Hist.eps", width=6,
#           height=3, horizontal = FALSE, onefile = FALSE, paper = "special");

## OPTIONS ##
par(mfrow=c(1,2));

## HISTOGRAMS ##
hist(R, breaks=40, main="Histogram of Returns");
hist(trunc.model.R, breaks=40, main="Histogram of Model Returns");

## TURN DEVICE OFF
# dev.off();

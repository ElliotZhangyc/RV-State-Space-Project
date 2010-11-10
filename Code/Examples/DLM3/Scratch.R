#phi.v3 = phi;
#W.v3 = W;
#mu.v3 = mu;
#z.v3 = z;

#phi.v3.hist = hist(phi.v3[3000:num.samples], breaks=40, plot=FALSE);
#W.v3.hist = hist(W.v3[3000:num.samples], breaks=40, plot=FALSE);
#mu.v3.hist = hist(mu.v3[3000:num.samples], breaks=40, plot=FALSE);
#z.v3.hist = hist(z.v3[3000:num.samples,50], breaks=40, plot=FALSE);

#hist(phi[3000:num.samples], breaks=40, prob=TRUE)
#lines(phi.v3.hist$mids, phi.v3.hist$density)
#hist(mu[3000:num.samples], breaks=40, prob=TRUE)
#lines(mu.v3.hist$mids, mu.v3.hist$density)
#hist(W[3000:num.samples], breaks=40, prob=TRUE)
#lines(W.v3.hist$mids, W.v3.hist$density)
#hist(z[3000:num.samples,50], breaks=40, prob=TRUE)
#lines(z.v3.hist$mids, z.v3.hist$density)

#jpeg(file="W_compare.jpg");
#plot(W.v2.hist$mids, W.v2.hist$density, col="red", type="l");
#lines(W.v3.hist$mids, W.v3.hist$density, col="blue");
#dev.off();

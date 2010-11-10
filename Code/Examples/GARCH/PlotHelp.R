# We want to produce a picture of distributions.  This will simply
# brute force calculate the posterior.  No Markov Chain techniques are
# used.  We may not even normalize the posterior.

# Number of points along an edge.
N = 51;

# The x and y grids.  They are equispaced on the intervals [0,1] and
# [0,1].
xgrid = seq(0,1,1/(N-1));
ygrid = xgrid;

# Set up the data structure to hold our posterior distribution.  We
# can easily visualize up to two dimensions.  We will use GnuPlot to
# visualize, which can hanle richer data by picking out specific
# columns, though we don't currently do that.
thedata = matrix(0, N*N, 3);

#for(i in 1:N){
#  for(j in 1:N){
#    thedata[(i-1)*N+j,] = c(xgrid[i],ygrid[j],
#         dmvnorm(c(xgrid[i],ygrid[j], output$par[3]), output$par, diag(3)*0.02));
#  }
#}

#for(i in 1:N){
#  for(j in 1:N){
    # When using mvt make sure that you tell it you do not want the
    # log dendsity.
#    thedata[(i-1)*N+j,] = c(xgrid[i],ygrid[j],
#         dmvt(c(xgrid[i],ygrid[j], output$par[3]),
#              output$par, C.cand, 1, log=F) );
#  }
#}

#for(i in 1:N){
#  for(j in 1:N){
#    thedata[(i-1)*N+j,] = c(xgrid[i],ygrid[j],
#         dmvt(c(xgrid[i]-output$par[1],ygrid[j]-output$par[2],0.06314999-output$par[3]),
#              rep(0,3), diag(3), 1, log=FALSE) );
#  }
#}

#for(i in 1:N){
#  for(j in 1:N){
#    thedata[(i-1)*N+j,] = c(xgrid[i],ygrid[j],
#         dmvt(c(xgrid[i],ygrid[j]),
#              output$par[1:2], C.cand[1:2,1:2], 1, log=FALSE) );
#  }
#}

for(i in 1:N){
  for(j in 1:N){
    thedata[(i-1)*N+j,] = c(xgrid[i],ygrid[j],
             likelihoodkernel(c(xgrid[i],ygrid[j], output$par[3])) );
  }
}

# Write the data to file.  The help file on write tells us that we
# must transpose the data for the file to match our row/column
# orientation from R.

# Normalize it for a comparison across different files.
normC = sum(thedata[,3]);
thedata[,3] = thedata[,3]/normC;

# MAKE SURE TO CHECK THE NAME OF THE FILE.
write(t(thedata), file="posterior.dat", ncolumns=3)

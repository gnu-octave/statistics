%!demo
%! ## Create a Log-uniform distribution with fixed parameters Lower=1 and
%! ## Upper=4 and plot its PDF.
%! pd = LoguniformDistribution (1, 4);
%! plot (pd)
%! title ("Fixed Log-uniform distribution with Lower=1 and Upper=4")

%!demo
%! ## Generate a data set of 100 random samples from a Log-uniform
%! ## distribution with parameters Lower=1 and Upper=10. Plot its CDF.
%! rand ("seed", 21);
%! data = exp (unifrnd (log (1), log (10), 100, 1));
%! pd = LoguniformDistribution (1, 10);
%! plot (pd, "PlotType", "cdf")
%! title ("Log-uniform distribution with Lower=1 and Upper=10")
%! 
%! ## Use this to visualize the CDF, useful for understanding cumulative
%! ## probabilities in the distribution.

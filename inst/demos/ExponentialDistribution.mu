%!demo
%! ## Create an Exponential distribution with default parameter
%! pd = makedist ("Exponential")
%!
%! ## Query parameter 'mu' (mean parameter)
%! pd.mu
%!
%! ## Set parameter 'mu'
%! pd.mu = 2
%!
%! ## Use this to initialize or modify the mean parameter of an Exponential
%! ## distribution. The mean parameter must be a positive real scalar.

%!demo
%! ## Create an Exponential distribution object by calling its constructor
%! pd = ExponentialDistribution (1.5)
%!
%! ## Query parameter 'mu'
%! pd.mu
%!
%! ## This demonstrates direct construction with a specific mean parameter,
%! ## useful for modeling waiting times with a known average.

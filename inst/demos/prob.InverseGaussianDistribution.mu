%!demo
%! ## Create an Inverse Gaussian distribution with default parameters
%! pd = makedist ("InverseGaussian")
%!
%! ## Query parameter 'mu' (mean parameter)
%! pd.mu
%!
%! ## Set parameter 'mu'
%! pd.mu = 2
%!
%! ## Use this to initialize or modify the mean parameter of an Inverse Gaussian
%! ## distribution. The mean parameter must be a positive real scalar.

%!demo
%! ## Create an Inverse Gaussian distribution object by calling its constructor
%! pd = InverseGaussianDistribution (1.5, 2)
%!
%! ## Query parameter 'mu'
%! pd.mu
%!
%! ## This demonstrates direct construction with a specific mean parameter,
%! ## useful for modeling non-negative skewed data with a known mean.

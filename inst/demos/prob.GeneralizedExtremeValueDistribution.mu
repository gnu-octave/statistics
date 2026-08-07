%!demo
%! ## Create a Generalized Extreme Value distribution with default parameters
%! pd = makedist ("GeneralizedExtremeValue")
%!
%! ## Query parameter 'mu' (location parameter)
%! pd.mu
%!
%! ## Set parameter 'mu'
%! pd.mu = 0.5
%!
%! ## Use this to initialize or modify the location parameter, which shifts the
%! ## Generalized Extreme Value distribution. It must be a real scalar.

%!demo
%! ## Create a Generalized Extreme Value distribution object by calling its constructor
%! pd = GeneralizedExtremeValueDistribution (0.1, 1, 0.5)
%!
%! ## Query parameter 'mu'
%! pd.mu
%!
%! ## This demonstrates setting the location parameter directly via the constructor,
%! ## useful for modeling the central tendency of extreme value data.

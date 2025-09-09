%!demo
%! ## Create an Extreme Value distribution with default parameters
%! pd = makedist ("ExtremeValue")
%!
%! ## Query parameter 'mu' (location parameter)
%! pd.mu
%!
%! ## Set parameter 'mu'
%! pd.mu = 2
%!
%! ## Use this to initialize or modify the location parameter of an Extreme Value
%! ## distribution. The location parameter must be a real scalar.

%!demo
%! ## Create an Extreme Value distribution object by calling its constructor
%! pd = ExtremeValueDistribution (1.5, 0.5)
%!
%! ## Query parameter 'mu'
%! pd.mu
%!
%! ## This demonstrates direct construction with a specific location parameter,
%! ## useful for modeling the position of maxima in data sets.

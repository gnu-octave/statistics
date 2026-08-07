%!demo
%! ## Create a Generalized Extreme Value distribution with default parameters
%! pd = makedist ("GeneralizedExtremeValue")
%!
%! ## Query parameter 'k' (shape parameter)
%! pd.k
%!
%! ## Set parameter 'k'
%! pd.k = 0.2
%!
%! ## Use this to initialize or modify the shape parameter of a Generalized Extreme
%! ## Value distribution. The shape parameter determines the tail behavior and must
%! ## be a real scalar.

%!demo
%! ## Create a Generalized Extreme Value distribution object by calling its constructor
%! pd = GeneralizedExtremeValueDistribution (0.1, 1, 0)
%!
%! ## Query parameter 'k'
%! pd.k
%!
%! ## This demonstrates direct construction with a specific shape parameter,
%! ## useful for modeling extreme value data with known tail behavior.

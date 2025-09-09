%!demo
%! ## Create a Half-normal distribution with default parameters
%! pd = makedist ("HalfNormal")
%!
%! ## Query parameter 'mu' (location parameter)
%! pd.mu
%!
%! ## Set parameter 'mu'
%! pd.mu = 1
%!
%! ## Use this to initialize or modify the location parameter of a Half-normal
%! ## distribution. The location parameter must be a real scalar, often set to 0
%! ## for modeling positive deviations from zero, such as measurement errors.

%!demo
%! ## Create a Half-normal distribution object by calling its constructor
%! pd = HalfNormalDistribution (1.5, 2)
%!
%! ## Query parameter 'mu'
%! pd.mu
%!
%! ## This demonstrates direct construction with a specific location parameter,
%! ## useful for modeling data shifted from zero, like distances or folded normals.

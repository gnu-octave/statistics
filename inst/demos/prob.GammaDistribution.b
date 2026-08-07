%!demo
%! ## Create a Gamma distribution with default parameters
%! pd = makedist ("Gamma")
%!
%! ## Query parameter 'b' (scale parameter)
%! pd.b
%!
%! ## Set parameter 'b'
%! pd.b = 2
%!
%! ## Use this to initialize or modify the scale parameter in a Gamma
%! ## distribution. The scale parameter must be a positive real scalar.

%!demo
%! ## Create a Gamma distribution object by calling its constructor
%! pd = GammaDistribution (2, 1)
%!
%! ## Query parameter 'b'
%! pd.b
%!
%! ## This shows how to set the scale parameter directly via the constructor,
%! ## ideal for scaling the distribution in lifetime models.

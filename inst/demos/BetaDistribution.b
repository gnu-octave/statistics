%!demo
%! ## Create a beta distribution with default parameters
%! pd = makedist ("Beta")
%!
%! ## Query parameter 'b' (second shape parameter)
%! pd.b
%!
%! ## Set parameter 'b'
%! pd.b = 2
%!
%! ## Use this to initialize or modify the second shape parameter of a beta
%! ## distribution. The parameter 'b' must be a positive real scalar.

%!demo
%! ## Create a beta distribution object by calling its constructor
%! pd = BetaDistribution (2, 3)
%!
%! ## Query parameter 'beta'
%! pd.b

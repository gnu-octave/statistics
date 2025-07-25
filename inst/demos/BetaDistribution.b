%!demo
%! ## Create a beta distribution with default parameters
%! pd = makedist ("beta")
%!
%! ## Query parameter 'beta'
%! pd.b
%!
%! ## Set parameter 'beta'
%! pd.b = 2

%!demo
%! ## Create a beta distribution object by calling its constructor
%! pd = BetaDistribution (2, 3)
%!
%! ## Query parameter 'beta'
%! pd.b

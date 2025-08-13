%!demo
%! ## Create a Burr distribution with default parameters
%! pd = makedist ("Burr")
%!
%! ## Query parameter 'k' (second shape parameter)
%! pd.k
%!
%! ## Set parameter 'k'
%! pd.k = 2
%!
%! ## Use this to initialize or modify the second shape parameter in a Burr
%! ## distribution. The second shape parameter must be a positive real scalar.

%!demo
%! ## Create a Burr distribution object by calling its constructor
%! pd = BurrDistribution (1, 2, 2)
%!
%! ## Query parameter 'k'
%! pd.k
%!
%! ## This shows how to set the second shape parameter directly via the constructor,
%! ## ideal for tailoring the tail behavior in distributions for income data.

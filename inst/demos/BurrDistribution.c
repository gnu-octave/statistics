%!demo
%! ## Create a Burr distribution with default parameters
%! pd = makedist ("Burr")
%!
%! ## Query parameter 'c' (first shape parameter)
%! pd.c
%!
%! ## Set parameter 'c'
%! pd.c = 3
%!
%! ## Use this to initialize or modify the first shape parameter in a Burr
%! ## distribution. The first shape parameter must be a positive real scalar.

%!demo
%! ## Create a Burr distribution object by calling its constructor
%! pd = BurrDistribution (1, 3, 1)
%!
%! ## Query parameter 'c'
%! pd.c
%!
%! ## This shows how to set the first shape parameter directly via the constructor,
%! ## useful for modeling variability in non-negative data like household income.

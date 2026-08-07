%!demo
%! ## Create a Triangular distribution with default parameters
%! pd = makedist ("Triangular", "A", 0, "B", 0.5, "C", 1);
%!
%! ## Query parameter 'C' (upper limit)
%! pd.C
%!
%! ## Set parameter 'C'
%! pd.C = 2
%!
%! ## Use this to initialize or modify the upper limit parameter of a Triangular
%! ## distribution. The upper limit must be a real scalar greater than A,
%! ## defining the maximum possible value in scenarios like budgeting or scheduling.

%!demo
%! ## Create a Triangular distribution object by calling its constructor
%! pd = TriangularDistribution (1, 2, 3);
%!
%! ## Query parameter 'C'
%! pd.C
%!
%! ## This demonstrates setting the upper limit directly via the constructor,
%! ## useful for modeling the maximum bound in applications like cost or time estimates.

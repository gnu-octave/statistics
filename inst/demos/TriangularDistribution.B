%!demo
%! ## Create a Triangular distribution with default parameters
%! pd = makedist ("Triangular", "A", 0, "B", 0.5, "C", 1);
%!
%! ## Query parameter 'B' (peak location)
%! pd.B
%!
%! ## Set parameter 'B'
%! pd.B = 0.7
%!
%! ## Use this to initialize or modify the peak location parameter of a Triangular
%! ## distribution. The peak location must be a real scalar between A and C,
%! ## representing the most likely value in applications like risk analysis.

%!demo
%! ## Create a Triangular distribution object by calling its constructor
%! pd = TriangularDistribution (1, 2, 3);
%!
%! ## Query parameter 'B'
%! pd.B
%!
%! ## This shows how to set the peak location directly via the constructor,
%! ## ideal for modeling the mode of triangularly distributed data, such as
%! ## expected task completion times.

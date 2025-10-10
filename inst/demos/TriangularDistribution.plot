%!demo
%! ## Create a Triangular distribution with fixed parameters A=0, B=1, C=2
%! ## and plot its PDF.
%! pd = TriangularDistribution (0, 1, 2);
%! plot (pd)
%! title ("Triangular distribution with A=0, B=1, C=2")
%!
%! ## Use this to visualize the PDF of a Triangular distribution with fixed
%! ## parameters, useful for understanding the shape of the distribution.

%!demo
%! ## Generate a data set of 100 random samples from a Triangular distribution
%! ## with parameters A=0, B=1, C=2. Plot its CDF.
%! rand ("seed", 21);
%! data = trirnd (0, 1, 2, 100, 1);
%! pd = makedist ("Triangular", "A", 0, "B", 1, "C", 2);
%! plot (pd, "PlotType", "cdf")
%! title ("Triangular distribution with A=0, B=1, C=2")
%! legend ({"Fitted CDF"}, "location", "southeast")
%!
%! ## Use this to visualize the CDF of a Triangular distribution, useful for
%! ## assessing cumulative probabilities in bounded data scenarios.

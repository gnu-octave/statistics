%!demo
%! ## Plot various CDFs from the Triangular distribution
%! x = -1:0.01:6;
%! pd1 = TriangularDistribution (0, 1, 2);
%! pd2 = TriangularDistribution (1, 2, 3);
%! pd3 = TriangularDistribution (2, 3, 4);
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"A=0, B=1, C=2", "A=1, B=2, C=3", "A=2, B=3, C=4"}, ...
%!         "location", "southeast")
%! title ("Triangular CDF")
%! xlabel ("Values in x (A <= x <= C)")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different Triangular distributions, showing how probability accumulates
%! ## within the bounds, useful in risk assessment or forecasting.

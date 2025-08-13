%!demo
%! ## Plot various CDFs from the Burr distribution
%! x = 0:0.01:5;
%! pd1 = makedist ("Burr", "alpha", 1, "c", 2, "k", 1);
%! pd2 = makedist ("Burr", "alpha", 1, "c", 3, "k", 1);
%! pd3 = makedist ("Burr", "alpha", 1, "c", 4, "k", 1);
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"alpha=1, c=2, k=1", "alpha=1, c=3, k=1", "alpha=1, c=4, k=1"}, ...
%!         "location", "southeast")
%! title ("Burr CDF")
%! xlabel ("values in x")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different Burr distributions, showing how probability accumulates
%! ## over non-negative values like income levels.

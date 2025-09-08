%!demo
%! ## Plot various CDFs from the Gamma distribution
%! x = 0:0.01:10;
%! pd1 = makedist ("Gamma", "a", 1, "b", 1);
%! pd2 = makedist ("Gamma", "a", 2, "b", 1);
%! pd3 = makedist ("Gamma", "a", 5, "b", 1);
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"a = 1, b = 1", "a = 2, b = 1", "a = 5, b = 1"}, ...
%!         "location", "southeast")
%! title ("Gamma CDF")
%! xlabel ("Values")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different Gamma distributions, showing how probability accumulates
%! ## over values.

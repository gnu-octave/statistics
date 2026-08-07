%!demo
%! ## Plot various CDFs from the Beta distribution
%! x = 0:0.01:1;
%! pd1 = makedist ("Beta", "a", 0.5, "b", 0.5);
%! pd2 = makedist ("Beta", "a", 2, "b", 2);
%! pd3 = makedist ("Beta", "a", 5, "b", 2);
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"a = 0.5, b = 0.5", "a = 2, b = 2", "a = 5, b = 2"}, ...
%!         "location", "southeast")
%! title ("Beta CDF")
%! xlabel ("Value")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different beta distributions, showing how probability accumulates
%! ## over the interval [0, 1].

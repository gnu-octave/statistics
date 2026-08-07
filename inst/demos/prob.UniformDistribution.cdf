%!demo
%! ## Plot various CDFs from the Uniform distribution
%! x = -1:0.01:12;
%! pd1 = makedist ("Uniform", "Lower", 0, "Upper", 5);
%! pd2 = makedist ("Uniform", "Lower", 2, "Upper", 8);
%! pd3 = makedist ("Uniform", "Lower", 4, "Upper", 10);
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"Lower = 0, Upper = 5", "Lower = 2, Upper = 8", "Lower = 4, Upper = 10"}, ...
%!         "location", "southeast")
%! title ("Uniform CDF")
%! xlabel ("Values in x")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different Uniform distributions, showing how probability accumulates
%! ## over the defined interval, useful in probability assessments or simulations.

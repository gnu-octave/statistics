%!demo
%! ## Plot various CDFs from the Binomial distribution
%! x = 0:10;
%! pd1 = makedist ("Binomial", "N", 10, "p", 0.2);
%! pd2 = makedist ("Binomial", "N", 10, "p", 0.5);
%! pd3 = makedist ("Binomial", "N", 10, "p", 0.8);
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "*b", x, p2, "*g", x, p3, "*r")
%! grid on
%! legend ({"N = 10, p = 0.2", "N = 10, p = 0.5", "N = 10, p = 0.8"}, ...
%!          "location", "southeast")
%! title ("Binomial CDF")
%! xlabel ("Number of successes")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different binomial distributions, showing how probability accumulates.

%!demo
%! ## Plot various CDFs from the Multinomial distribution
%! x = 1:0.1:6;
%! pd1 = MultinomialDistribution ([0.4, 0.3, 0.3]);
%! pd2 = MultinomialDistribution ([0.2, 0.2, 0.2, 0.2, 0.2]);
%! pd3 = MultinomialDistribution ([1/6, 1/6, 1/6, 1/6, 1/6, 1/6]);
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"3 outcomes", "5 outcomes", "6 outcomes (dice)"}, ...
%!         "location", "southeast")
%! title ("Multinomial CDF")
%! xlabel ("values in x (outcome index)")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different Multinomial distributions, showing how probability
%! ## accumulates across categorical outcomes, useful in decision analysis.

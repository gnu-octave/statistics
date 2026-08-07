%!demo
%! ## Plot various CDFs from the Negative Binomial distribution
%! x = 0:20;
%! data1 = nbinrnd(5, 0.3, 10000, 1);
%! data2 = nbinrnd(5, 0.5, 10000, 1);
%! data3 = nbinrnd(5, 0.7, 10000, 1);
%! pd1 = fitdist (data1, "NegativeBinomial");
%! pd2 = fitdist (data2, "NegativeBinomial");
%! pd3 = fitdist (data3, "NegativeBinomial");
%! p1 = cdf (pd1, x);
%! p2 = cdf (pd2, x);
%! p3 = cdf (pd3, x);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"R=5, P=0.3", "R=5, P=0.5", "R=5, P=0.7"}, "location", "southeast")
%! title ("Negative Binomial CDF")
%! xlabel ("values in x (number of failures)")
%! ylabel ("Cumulative probability")
%!
%! ## Use this to compute and visualize the cumulative distribution function
%! ## for different Negative Binomial distributions, showing how probability
%! ## accumulates for count data, useful in modeling overdispersed counts.

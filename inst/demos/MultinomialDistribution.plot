%!demo
%! ## Create a Multinomial distribution with fixed parameters and plot its PDF.
%! pd = MultinomialDistribution ([0.1, 0.2, 0.3, 0.2, 0.1, 0.1]);
%! plot (pd)
%! title ("Multinomial distribution PDF")

%!demo
%! ## Create a Multinomial distribution and plot its CDF.
%! pd = MultinomialDistribution ([0.1, 0.2, 0.3, 0.2, 0.1, 0.1]);
%! plot (pd, "PlotType", "cdf")
%! title ("Multinomial distribution CDF")
%!
%! ## Use this to visualize the cumulative distribution function,
%! ## useful for understanding probability accumulation across outcomes.

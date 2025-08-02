%!demo
%! ## Generate a data set of 100 random samples from a Birnbaum-Saunders
%! ## distribution with parameters beta = 1 and gamma = 0.5. Fit a Birnbaum-Saunders
%! ## distribution to this data and plot its CDF superimposed over an empirical CDF.
%!
%! pd = makedist ("BirnbaumSaunders", "beta", 1, "gamma", 0.5)
%! randg ("seed", 21);
%! data = random (pd, 100, 1);
%! pd = fitdist (data, "BirnbaumSaunders")
%! plot (pd, "PlotType", "cdf")
%! title (sprintf ("Fitted Birnbaum-Saunders distribution with beta = %0.2f and gamma = %0.2f", ...
%!                 pd.beta, pd.gamma))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## Use this to visualize the fitted CDF compared to the empirical CDF of the
%! ## data, useful for assessing model fit.

%!demo
%! ## Generate a data set of 200 random samples from a Birnbaum-Saunders
%! ## distribution with parameters beta = 1 and gamma = 0.5. Display a probability
%! ## plot for the Birnbaum-Saunders distribution fit to the data.
%!
%! pd = makedist ("BirnbaumSaunders", "beta", 1, "gamma", 0.5)
%! randg ("seed", 21);
%! data = random (pd, 200, 1);
%! pd = fitdist (data, "BirnbaumSaunders")
%! plot (pd, "PlotType", "probability")
%! title (sprintf ("Probability plot of fitted Birnbaum-Saunders distribution with beta = %0.2f and gamma = %0.2f", ...
%!                 pd.beta, pd.gamma))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## This creates a probability plot to compare the fitted distribution to the
%! ## data, useful for checking if the Birnbaum-Saunders model is appropriate.

%!demo
%! ## Generate a data set of 100 random samples from a Binomial distribution
%! ## with parameters N = 10 and p = 0.3. Fit a Binomial distribution to this
%! ## data and plot its CDF superimposed over an empirical CDF of the data
%!
%! pd = makedist ("Binomial", "N", 10, "p", 0.3)
%! rand ("seed", 22);
%! data = random (pd, 100, 1);
%! pd = fitdist (data, "Binomial", "ntrials", 10)
%! plot (pd, "PlotType", "cdf", "Discrete", true)
%! title (sprintf ("Fitted Binomial distribution with N = %d and p = %0.2f", ...
%!                  pd.N, pd.p))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## Use this to visualize the fitted CDF compared to the empirical CDF of
%! ## the data, useful for assessing model fit.

%!demo
%! ## Generate a data set of 200 random samples from a Binomial distribution
%! ## with parameters N = 10 and p = 0.3. Display a probability plot for the
%! ## Binomial distribution fit to the data.
%!
%! pd = makedist ("Binomial", "N", 10, "p", 0.3)
%! rand ("seed", 22);
%! data = random (pd, 200, 1);
%! pd = fitdist (data, "Binomial", "ntrials", 10)
%! plot (pd, "PlotType", "probability", "Discrete", true)
%! title (sprintf (["Probability plot of fitted Binomial distribution with " ...
%!                  "N = %d and p = %0.2f"], pd.N, pd.p));
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast");
%!
%! ## This creates a probability plot to compare the fitted distribution to the
%! ## data, useful for checking if the binomial model is appropriate.
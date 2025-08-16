%!demo
%! ## Create an Exponential distribution with fixed parameter mu = 2 and plot its PDF.
%!
%! pd = makedist ("Exponential", "mu", 2)
%! plot (pd)
%! title ("Fixed Exponential distribution with mu = 2")

%!demo
%! ## Generate a data set of 100 random samples from an Exponential
%! ## distribution with parameter mu = 2. Fit an Exponential
%! ## distribution to this data and plot its CDF superimposed over an empirical
%! ## CDF.
%!
%! pd_fixed = makedist ("Exponential", "mu", 2)
%! rand ("seed", 5);
%! data = random (pd_fixed, 100, 1);
%! pd_fitted = fitdist (data, "Exponential")
%! plot (pd_fitted, "PlotType", "cdf")
%! txt = "Fitted Exponential distribution with mu = %0.2f";
%! title (sprintf (txt, pd_fitted.mu))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## Use this to visualize the fitted CDF compared to the empirical CDF of the
%! ## data, useful for assessing model fit.

%!demo
%! ## Generate a data set of 200 random samples from an Exponential
%! ## distribution with parameter mu = 2. Display a probability
%! ## plot for the Exponential distribution fit to the data.
%!
%! pd_fixed = makedist ("Exponential", "mu", 2)
%! rand ("seed", 5);
%! data = random (pd_fixed, 200, 1);
%! pd_fitted = fitdist (data, "Exponential")
%! plot (pd_fitted, "PlotType", "probability")
%! txt = strcat ("Probability plot of fitted Exponential", ...
%!               " distribution with mu = %0.2f");
%! title (sprintf (txt, pd_fitted.mu))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## This creates a probability plot to compare the fitted distribution to the
%! ## data, useful for checking if the Exponential model is appropriate.

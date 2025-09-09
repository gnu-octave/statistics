%!demo
%! ## Create a Gamma distribution with fixed parameters a = 2 and
%! ## b = 1 and plot its PDF.
%!
%! pd = makedist ("Gamma", "a", 2, "b", 1)
%! plot (pd)
%! title ("Fixed Gamma distribution with a = 2 and b = 1")

%!demo
%! ## Generate a data set of 100 random samples from a Gamma
%! ## distribution with parameters a = 2 and b = 1. Fit a Gamma
%! ## distribution to this data and plot its CDF superimposed over an empirical
%! ## CDF.
%!
%! pd_fixed = makedist ("Gamma", "a", 2, "b", 1)
%! rand ("seed", 21);
%! data = gamrnd (2, 1, 100, 1);
%! pd_fitted = fitdist (data, "Gamma")
%! plot (pd_fitted, "PlotType", "cdf")
%! txt = "Fitted Gamma distribution with a = %0.2f and b = %0.2f";
%! title (sprintf (txt, pd_fitted.a, pd_fitted.b))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## Use this to visualize the fitted CDF compared to the empirical CDF of the
%! ## data, useful for assessing model fit.

%!demo
%! ## Generate a data set of 200 random samples from a Gamma
%! ## distribution with parameters a = 2 and b = 1. Display a probability
%! ## plot for the Gamma distribution fit to the data.
%!
%! pd_fixed = makedist ("Gamma", "a", 2, "b", 1)
%! rand ("seed", 21);
%! data = gamrnd (2, 1, 200, 1);
%! pd_fitted = fitdist (data, "Gamma")
%! plot (pd_fitted, "PlotType", "probability")
%! txt = strcat ("Probability plot of fitted Gamma", ...
%!               " distribution with a = %0.2f and b = %0.2f");
%! title (sprintf (txt, pd_fitted.a, pd_fitted.b))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## This creates a probability plot to compare the fitted distribution to the
%! ## data, useful for checking if the Gamma model is appropriate.

%!demo
%! ## Create a Burr distribution with fixed parameters alpha = 1, c = 2, k = 1
%! ## and plot its PDF.
%!
%! pd = makedist ("Burr", "alpha", 1, "c", 2, "k", 1)
%! plot (pd)
%! title ("Fixed Burr distribution with alpha = 1, c = 2, k = 1")

%!demo
%! ## Generate a data set of 100 random samples from a Burr distribution with
%! ## parameters alpha = 1, c = 2, k = 1. Fit a Burr distribution to this data
%! ## and plot its CDF superimposed over an empirical CDF.
%!
%! pd_fixed = makedist ("Burr", "alpha", 1, "c", 2, "k", 1)
%! rand ("seed", 21);
%! data = random (pd_fixed, 100, 1);
%! pd_fitted = fitdist (data, "Burr")
%! plot (pd_fitted, "PlotType", "cdf")
%! txt = "Fitted Burr distribution with alpha = %0.2f, c = %0.2f, k = %0.2f";
%! title (sprintf (txt, pd_fitted.alpha, pd_fitted.c, pd_fitted.k))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## Use this to visualize the fitted CDF compared to the empirical CDF of the
%! ## data, useful for assessing model fit in income data.

%!demo
%! ## Generate a data set of 200 random samples from a Burr distribution with
%! ## parameters alpha = 1, c = 2, k = 1. Display a probability plot for the
%! ## Burr distribution fit to the data.
%!
%! pd_fixed = makedist ("Burr", "alpha", 1, "c", 2, "k", 1)
%! rand ("seed", 21);
%! data = random (pd_fixed, 200, 1);
%! pd_fitted = fitdist (data, "Burr")
%! plot (pd_fitted, "PlotType", "probability")
%! txt = strcat ("Probability plot of fitted Burr distribution with ", ...
%!               "alpha = %0.2f, c = %0.2f, k = %0.2f");
%! title (sprintf (txt, pd_fitted.alpha, pd_fitted.c, pd_fitted.k))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## This creates a probability plot to compare the fitted distribution to the
%! ## data, useful for checking if the Burr model is appropriate for the dataset.

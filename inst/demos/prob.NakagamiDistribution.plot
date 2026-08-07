%!demo
%! ## Create a Nakagami distribution with fixed parameters mu = 1 and
%! ## omega = 1 and plot its PDF.
%!
%! data = nakarnd (1, 1, 10000, 1);
%! pd = fitdist (data, "Nakagami");
%! plot (pd)
%! title ("Fixed Nakagami distribution with mu = 1 and omega = 1")

%!demo
%! ## Generate a data set of 100 random samples from a Nakagami
%! ## distribution with parameters mu = 1 and omega = 1. Fit a Nakagami
%! ## distribution to this data and plot its CDF superimposed over an empirical
%! ## CDF.
%!
%! rand ("seed", 21);
%! data = nakarnd (1, 1, 100, 1);
%! pd_fitted = fitdist (data, "Nakagami");
%! plot (pd_fitted, "PlotType", "cdf")
%! txt = "Fitted Nakagami distribution with mu = %0.2f and omega = %0.2f";
%! title (sprintf (txt, pd_fitted.mu, pd_fitted.omega))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## Use this to visualize the fitted CDF compared to the empirical CDF of the
%! ## data, useful for assessing model fit.

%!demo
%! ## Generate a data set of 200 random samples from a Nakagami
%! ## distribution with parameters mu = 1 and omega = 1. Display a probability
%! ## plot for the Nakagami distribution fit to the data.
%!
%! rand ("seed", 21);
%! data = nakarnd (1, 1, 200, 1);
%! pd_fitted = fitdist (data, "Nakagami");
%! plot (pd_fitted, "PlotType", "probability")
%! txt = strcat ("Probability plot of fitted Nakagami", ...
%!               " distribution with mu = %0.2f and omega = %0.2f");
%! title (sprintf (txt, pd_fitted.mu, pd_fitted.omega))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## This creates a probability plot to compare the fitted distribution to the
%! ## data, useful for checking if the Nakagami model is appropriate.

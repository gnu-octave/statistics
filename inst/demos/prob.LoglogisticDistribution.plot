%!demo
%! ## Create a Log-logistic distribution with fixed parameters mu = 0 and
%! ## sigma = 1 and plot its PDF.
%!
%! data = loglrnd (0, 1, 10000, 1);
%! pd = fitdist (data, "Loglogistic");
%! x = linspace (0.01, 20, 1000);
%! y = pdf (pd, x);
%! plot (x, y, "b", "LineWidth", 2)
%! grid on
%! title ("Fixed Log-logistic distribution with mu = 0 and sigma = 1")
%! xlabel ("x")
%! ylabel ("PDF")

%!demo
%! ## Generate a data set of 100 random samples from a Log-logistic
%! ## distribution with parameters mu = 0 and sigma = 1. Fit a Log-logistic
%! ## distribution to this data and plot its CDF superimposed over an empirical
%! ## CDF.
%!
%! rand ("seed", 21);
%! data = loglrnd (0, 1, 100, 1);
%! pd_fitted = fitdist (data, "Loglogistic");
%! ecdf (data);
%! hold on;
%! x = linspace (icdf (pd_fitted, 0.01), icdf (pd_fitted, 0.99), 1000);
%! y = cdf (pd_fitted, x);
%! plot (x, y, "r", "LineWidth", 2);
%! txt = "Fitted Log-logistic distribution with mu = %0.2f and sigma = %0.2f";
%! title (sprintf (txt, pd_fitted.mu, pd_fitted.sigma))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%! xlabel ("x")
%! ylabel ("CDF")
%! grid on
%! hold off;
%!
%! ## Use this to visualize the fitted CDF compared to the empirical CDF of the
%! ## data, useful for assessing model fit in skewed positive data.

%!demo
%! ## Generate a data set of 200 random samples from a Log-logistic
%! ## distribution with parameters mu = 0 and sigma = 1. Display a probability
%! ## plot for the Log-logistic distribution fit to the data.
%!
%! rand ("seed", 21);
%! data = loglrnd (0, 1, 200, 1);
%! pd_fitted = fitdist (data, "Loglogistic");
%! plot (pd_fitted, "PlotType", "probability")
%! txt = strcat ("Probability plot of fitted Log-logistic", ...
%!               " distribution with mu = %0.2f and sigma = %0.2f");
%! title (sprintf (txt, pd_fitted.mu, pd_fitted.sigma))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## This creates a probability plot to compare the fitted distribution to the
%! ## data, useful for checking if the Log-logistic model captures the tail behavior.

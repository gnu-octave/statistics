%!demo
%! ## Create a Poisson distribution with fixed parameter lambda=5 and plot its PDF.
%! data = poissrnd (5, 10000, 1);
%! pd = fitdist (data, "Poisson");
%! plot (pd)
%! title ("Fixed Poisson distribution with lambda = 5")

%!demo
%! ## Generate a data set of 100 random samples from a Poisson distribution with
%! ## lambda=5. Fit a Poisson distribution to this data and plot its CDF superimposed
%! ## over an empirical CDF.
%! rand ("seed", 21);
%! data = poissrnd (5, 100, 1);
%! pd_fitted = fitdist (data, "Poisson");
%! plot (pd_fitted, "PlotType", "cdf")
%! txt = "Fitted Poisson distribution with lambda = %0.2f";
%! title (sprintf (txt, pd_fitted.lambda))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## Use this to visualize the fitted CDF compared to the empirical CDF of the data,
%! ## useful for assessing model fit in count data.

%!demo
%! ## Generate a data set of 200 random samples from a Poisson distribution with
%! ## lambda=5. Display a probability plot for the Poisson distribution fit to the data.
%! rand ("seed", 21);
%! data = poissrnd (5, 200, 1);
%! pd_fitted = fitdist (data, "Poisson");
%! plot (pd_fitted, "PlotType", "probability")
%! txt = strcat ("Probability plot of fitted Poisson distribution with lambda = %0.2f");
%! title (sprintf (txt, pd_fitted.lambda))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## This creates a probability plot to compare the fitted distribution to the data,
%! ## useful for checking if the Poisson model is appropriate for count data.

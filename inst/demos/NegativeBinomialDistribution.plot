%!demo
%! ## Create a Negative Binomial distribution with fixed parameters R=5 and
%! ## P=0.5 and plot its PDF.
%!
%! data = nbinrnd(5, 0.5, 10000, 1);
%! pd = fitdist (data, "NegativeBinomial");
%! plot (pd)
%! title ("Fixed Negative Binomial distribution with R=5 and P=0.5")

%!demo
%! ## Generate a data set of 100 random samples from a Negative Binomial
%! ## distribution with parameters R=5 and P=0.5. Fit a Negative Binomial
%! ## distribution to this data and plot its CDF superimposed over an empirical
%! ## CDF.
%!
%! rand ("seed", 21);
%! data = nbinrnd(5, 0.5, 100, 1);
%! pd_fitted = fitdist (data, "NegativeBinomial");
%! plot (pd_fitted, "PlotType", "cdf")
%! txt = "Fitted Negative Binomial distribution with R=%0.2f and P=%0.2f";
%! title (sprintf (txt, pd_fitted.R, pd_fitted.P))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## Use this to visualize the fitted CDF compared to the empirical CDF of the
%! ## data, useful for assessing model fit in count data.

%!demo
%! ## Generate a data set of 200 random samples from a Negative Binomial
%! ## distribution with parameters R=5 and P=0.5. Display a probability
%! ## plot for the Negative Binomial distribution fit to the data.
%!
%! rand ("seed", 21);
%! data = nbinrnd(5, 0.5, 200, 1);
%! pd_fitted = fitdist (data, "NegativeBinomial");
%! plot (pd_fitted, "PlotType", "probability")
%! txt = strcat ("Probability plot of fitted Negative Binomial", ...
%!               " distribution with R=%0.2f and P=%0.2f");
%! title (sprintf (txt, pd_fitted.R, pd_fitted.P))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")
%!
%! ## This creates a probability plot to compare the fitted distribution to the
%! ## data, useful for checking if the Negative Binomial model is appropriate.

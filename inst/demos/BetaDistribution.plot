%!demo
%! ## Create a Beta distribution with fixed parameters a = 2 and b = 5, and
%! plot its PDF.
%!
%! pd = makedist ("BirnbaumSaunders", "a", 2, "b", 5)
%! plot (pd)
%! title ("Fixed Beta distribution with a = 2 and b = 5")

%!demo
%! ## Generate a data set of 100 random samples from a Beta distribution with
%! ## parameters a = 2 and b = 4.  Fit a Beta distribution to this data and plot
%! ## its CDF superimposed over an empirical CDF of the data
%!
%! pd_fixed = makedist ("Beta", "a", 2, "b", 4)
%! randg ("seed", 21);
%! data = random (pd_fixed, 100, 1);
%! pd_fitted = fitdist (data, "Beta")
%! plot (pd_fitted, "plottype", "cdf")
%! txt = "Fitted Beta distribution with a = %0.2f and b = %0.2f";
%! title (sprintf (txt, pd_fitted.a, pd_fitted.b))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "east")

%!demo
%! ## Generate a data set of 200 random samples from a Beta distribution with
%! ## parameters a = 2 and b = 4.  Display a probability plot for the Beta
%! ## distribution fit to the data.
%!
%! pd_fixed = makedist ("Beta", "a", 2, "b", 4)
%! randg ("seed", 21);
%! data = random (pd_fixed, 200, 1);
%! pd_fitted = fitdist (data, "Beta")
%! plot (pd_fitted, "plottype", "probability")
%! txt = "Probability plot of a fitted Beta distribution with a = %0.2f and b = %0.2f";
%! title (sprintf (txt, pd_fitted.a, pd_fitted.b))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")

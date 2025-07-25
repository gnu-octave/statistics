%!demo
%! ## Generate a data set of 100 random samples from a Beta distribution with
%! ## parameters a = 2 and b = 4.  Fit a Beta distribution to this data and plot
%! ## its CDF superimposed over an empirical CDF of the data
%!
%! pd = makedist ("Beta", "a", 2, "b", 4)
%! randg ("seed", 21);
%! data = random (pd, 100, 1);
%! pd = fitdist (data, "Beta")
%! plot (pd, "plottype", "cdf")
%! title (sprintf ("Fitted Beta distribution with a = %0.2f and b = %0.2f", ...
%!                 pd.a, pd.b))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "east")

%!demo
%! ## Generate a data set of 200 random samples from a Beta distribution with
%! ## parameters a = 2 and b = 4.  Display a probability plot for the Beta
%! ## distribution fit to the data.
%!
%! pd = makedist ("Beta", "a", 2, "b", 4)
%! randg ("seed", 21);
%! data = random (pd, 200, 1);
%! pd = fitdist (data, "Beta")
%! plot (pd, "plottype", "probability")
%! title (sprintf ("Probability plot of a fitted Beta distribution with a = %0.2f and b = %0.2f", ...
%!                 pd.a, pd.b))
%! legend ({"empirical CDF", "fitted CDF"}, "location", "southeast")

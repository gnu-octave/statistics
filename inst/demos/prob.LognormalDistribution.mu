%!demo
%! ## Create a Lognormal distribution with default parameters
%! data = lognrnd (0, 1, 10000, 1);
%! pd = fitdist (data, "Lognormal");
%!
%! ## Query parameter 'mu' (mean of logarithmic values)
%! pd.mu
%!
%! ## Set parameter 'mu'
%! pd.mu = 1
%!
%! ## Use this to initialize or modify the mean parameter of the logarithmic values
%! ## in a Lognormal distribution. The mu parameter must be a real scalar, often
%! ## representing the log-mean of multiplicative processes, like growth rates.

%!demo
%! ## Create a Lognormal distribution object by calling its constructor
%! pd = LognormalDistribution (1.5, 0.5)
%!
%! ## Query parameter 'mu'
%! pd.mu
%!
%! ## This demonstrates direct construction with a specific mu parameter,
%! ## useful for modeling skewed positive data shifted in log-scale, such as incomes or sizes.

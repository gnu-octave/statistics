%!demo
%! ## Create a Log-logistic distribution with default parameters
%! data = loglrnd (0, 1, 10000, 1);
%! pd = fitdist (data, "Loglogistic");
%!
%! ## Query parameter 'mu' (mean of logarithmic values)
%! pd.mu
%!
%! ## Set parameter 'mu'
%! pd.mu = 1
%!
%! ## Use this to initialize or modify the mean of the logarithmic values in a
%! ## Log-logistic distribution. The mu parameter must be a nonnegative real
%! ## scalar, often representing the location in log-space for modeling
%! ## positive skewed data like survival times or income distributions.

%!demo
%! ## Create a Log-logistic distribution object by calling its constructor
%! pd = LoglogisticDistribution (1.5, 2);
%!
%! ## Query parameter 'mu'
%! pd.mu
%!
%! ## This demonstrates direct construction with a specific mu parameter,
%! ## useful for modeling data with a known log-mean, such as in reliability
%! ## engineering or financial modeling.

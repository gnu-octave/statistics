%!demo
%! ## Create a Logistic distribution with default parameters
%! pd = makedist ("Logistic")
%!
%! ## Query parameter 'mu' (location parameter)
%! pd.mu
%!
%! ## Set parameter 'mu'
%! pd.mu = 2
%!
%! ## Use this to initialize or modify the location parameter of a Logistic
%! ## distribution. The location parameter must be a finite real scalar and
%! ## represents the center of the distribution, often used in regression models.

%!demo
%! ## Create a Logistic distribution object by calling its constructor
%! pd = LogisticDistribution (1.5, 0.5)
%!
%! ## Query parameter 'mu'
%! pd.mu
%!
%! ## This demonstrates direct construction with a specific location parameter,
%! ## useful for modeling data centered around a known value, such as in logistic regression.

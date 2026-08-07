%!demo
%! ## Create a Weibull distribution with default parameters
%! data = wblrnd (1, 1, 10000, 1);
%! pd = fitdist (data, "Weibull");
%!
%! ## Query parameter 'k' (shape parameter)
%! pd.k
%!
%! ## Set parameter 'k'
%! pd.k = 2
%!
%! ## Use this to initialize or modify the shape parameter in a Weibull
%! ## distribution. The shape parameter must be a positive real scalar, controlling
%! ## the failure rate behavior (k<1 decreasing, k=1 constant, k>1 increasing).

%!demo
%! ## Create a Weibull distribution object by calling its constructor
%! pd = WeibullDistribution (1, 3)
%!
%! ## Query parameter 'k'
%! pd.k
%!
%! ## This shows how to set the shape parameter directly via the constructor,
%! ## ideal for modeling different hazard functions in reliability engineering.

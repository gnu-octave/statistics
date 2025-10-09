%!demo
%! ## Create a Piecewise Linear distribution with default parameters
%! load patients
%! [f, x] = ecdf (Weight);
%! f = f(1:5:end);
%! x = x(1:5:end);
%! pd = PiecewiseLinearDistribution (x, f);
%!
%! ## Query parameter 'x' (vector of x values)
%! pd.x
%!
%! ## Set parameter 'x'
%! pd.x = [100; 120; 140; 160; 180; 200]
%!
%! ## Use this to initialize or modify the vector of x values where the CDF changes
%! ## slope in a Piecewise Linear distribution. The x vector must be real, finite,
%! ## and increasing, often derived from empirical data points.

%!demo
%! ## Create a Piecewise Linear distribution object by calling its constructor
%! pd = PiecewiseLinearDistribution ([0; 1; 2; 3], [0; 0.3; 0.7; 1])
%!
%! ## Query parameter 'x'
%! pd.x
%!
%! ## This demonstrates direct construction with specific x values, useful for
%! ## modeling custom empirical distributions or interpolating between known points.

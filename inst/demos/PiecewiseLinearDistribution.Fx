%!demo
%! ## Create a Piecewise Linear distribution with default parameters
%! load patients
%! [f, x] = ecdf (Weight);
%! f = f(1:5:end);
%! x = x(1:5:end);
%! pd = PiecewiseLinearDistribution (x, f);
%!
%! ## Query parameter 'Fx' (vector of CDF values)
%! pd.Fx
%!
%! ## Set parameter 'Fx'
%! pd.Fx = [0; 0.2; 0.5; 0.6; 0.95; 1]
%!
%! ## Use this to initialize or modify the vector of CDF values corresponding to x.
%! ## Fx must be real, between 0 and 1, non-decreasing, starting at 0 and ending at 1.

%!demo
%! ## Create a Piecewise Linear distribution object by calling its constructor
%! pd = PiecewiseLinearDistribution ([0; 1; 2; 3], [0; 0.3; 0.7; 1])
%!
%! ## Query parameter 'Fx'
%! pd.Fx
%!
%! ## This shows how to set CDF values directly via the constructor, ideal for
%! ## defining piecewise linear approximations of arbitrary distributions.

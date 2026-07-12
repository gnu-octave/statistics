## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{mdl} =} fitnlm (@var{X}, @var{y}, @var{modelfun}, @var{beta0})
## @deftypefnx {statistics} {@var{mdl} =} fitnlm (@var{tbl}, @var{modelfun}, @var{beta0})
## @deftypefnx {statistics} {@var{mdl} =} fitnlm (@dots{}, @var{Name}, @var{Value})
##
## Fit a nonlinear regression model.
##
## @code{@var{mdl} = fitnlm (@var{X}, @var{y}, @var{modelfun}, @var{beta0})}
## fits the nonlinear regression model @code{@var{y} = @var{modelfun}
## (@var{beta}, @var{X})} to the response vector @var{y} and the
## @math{n}-by-@math{p} predictor matrix @var{X}, starting the iterative fit
## from the coefficient vector
## @var{beta0}, and returns a @code{NonLinearModel} object.  @var{modelfun} is a
## function handle @code{@@(@var{b}, @var{X})} returning the fitted responses.
##
## @code{@var{mdl} = fitnlm (@var{tbl}, @var{modelfun}, @var{beta0})} takes the
## predictors and response from the table @var{tbl}; the last column is the
## response unless overridden by @qcode{'ResponseVar'}.
##
## The following @var{Name}/@var{Value} pairs are accepted:
##
## @multitable @columnfractions 0.2 0.75
## @headitem Name @tab Value
## @item @qcode{'CoefficientNames'} @tab a cell array of names for the
## coefficients (default @qcode{'b1'}, @qcode{'b2'}, @dots{}).
## @item @qcode{'Weights'} @tab a vector of nonnegative observation weights.
## @item @qcode{'ErrorModel'} @tab the error-variance model: @qcode{'constant'}
## (default), @qcode{'proportional'}, or @qcode{'combined'}.
## @item @qcode{'RobustWgtFun'} @tab the name of a robust weight function,
## enabling robust fitting (see @code{nlinfit}).
## @item @qcode{'Options'} @tab a statset-style options structure controlling
## the iterative fit (@qcode{MaxIter}, @qcode{TolFun}, @qcode{TolX}).
## @item @qcode{'PredictorVars'}, @qcode{'ResponseVar'} @tab for table input,
## the predictor and response variable names.
## @item @qcode{'VarNames'} @tab a cell array of @math{p + 1} variable names
## (predictors followed by the response) for numeric @var{X}.
## @item @qcode{'Exclude'} @tab observations to exclude from the fit.
## @end multitable
##
## @seealso{NonLinearModel, nlinfit, nlparci, nlpredci, fitlm, fitglm}
## @end deftypefn

function mdl = fitnlm (varargin)

  if (nargin < 3)
    print_usage ();
  endif

  if (istable (varargin{1}))
    mdl = NonLinearModel (varargin{1}, [], varargin{2}, varargin{3}, ...
                          varargin{4:end});
  else
    if (nargin < 4)
      print_usage ();
    endif
    mdl = NonLinearModel (varargin{1}, varargin{2}, varargin{3}, ...
                          varargin{4}, varargin{5:end});
  endif

endfunction

%!demo
%! ## Fit an exponential growth model and inspect the summary.
%! x = [1:10]';
%! y = [2.1;2.9;4.2;5.3;7.1;9.4;12.8;16.5;22.1;29.8];
%! modelfun = @(b, x) b(1) .* exp (b(2) .* x);
%! mdl = fitnlm (x, y, modelfun, [1; 0.3])

%!demo
%! ## Predictions with 95% confidence intervals on the fitted curve.
%! x = [1:10]';
%! y = [2.1;2.9;4.2;5.3;7.1;9.4;12.8;16.5;22.1;29.8];
%! modelfun = @(b, x) b(1) .* exp (b(2) .* x);
%! mdl = fitnlm (x, y, modelfun, [1; 0.3]);
%! [ypred, yci] = predict (mdl, [2.5; 5.5; 8.5])

%!shared X, y, modelfun, beta0
%! X = [1;2;3;4;5;6;7;8;9;10];
%! y = [2.1;2.9;4.2;5.3;7.1;9.4;12.8;16.5;22.1;29.8];
%! modelfun = @(b, x) b(1) .* exp (b(2) .* x);
%! beta0 = [1; 0.3];

## Values verified against MATLAB's fitnlm.
%!test
%! mdl = fitnlm (X, y, modelfun, beta0);
%! assert_equal (mdl.Coefficients.Estimate, ...
%!   [1.683747025; 0.286911087], 1e-6);
%! assert_equal (mdl.Coefficients.SE, [0.035194899; 0.002350913], 1e-6);
%! assert_equal (mdl.Coefficients.tStat, [47.8406555; 122.042406], -1e-4);
%! assert_equal (mdl.RMSE, 0.170942956, 1e-7);
%! assert_equal (mdl.SSE, 0.233771954, 1e-7);
%! assert_equal (mdl.SST, 750.976, 1e-3);
%!test
%! mdl = fitnlm (X, y, modelfun, beta0);
%! assert_equal (mdl.Rsquared.Ordinary, 0.999688709, 1e-8);
%! assert_equal (mdl.Rsquared.Adjusted, 0.999649798, 1e-8);
%! assert_equal (mdl.LogLikelihood, 4.590586096, 1e-6);
%! assert_equal (mdl.ModelCriterion.AIC, -5.181172193, 1e-6);
%! assert_equal (mdl.ModelCriterion.BIC, -4.576002007, 1e-6);
%!test
%! ## coefCI matches beta +/- t * SE with the error degrees of freedom.
%! mdl = fitnlm (X, y, modelfun, beta0);
%! ci = coefCI (mdl);
%! b = mdl.Coefficients.Estimate;  se = mdl.Coefficients.SE;
%! t = tinv (0.975, mdl.DFE);
%! assert_equal (ci, [b - t .* se, b + t .* se], 1e-12);
%!test
%! ## Table input gives the same fit as matrix input.
%! tbl = table (X, y, "VariableNames", {"x", "y"});
%! mdl = fitnlm (tbl, modelfun, beta0);
%! assert_equal (mdl.Coefficients.Estimate, [1.683747025; 0.286911087], 1e-6);
%! assert_equal (mdl.CoefficientNames, {"b1", "b2"});
%!test
%! ## Custom coefficient names.
%! mdl = fitnlm (X, y, modelfun, beta0, "CoefficientNames", {"A", "k"});
%! assert_equal (mdl.CoefficientNames, {"A", "k"});
%!test
%! ## coefTest reports a Wald F statistic versus the zero model.
%! mdl = fitnlm (X, y, modelfun, beta0);
%! [p, F, df] = coefTest (mdl);
%! assert_equal (df, 2);
%! assert (F > 1e5);
%! assert (p < 1e-10);
%!test
%! ## predict returns the fitted values and confidence intervals.
%! mdl = fitnlm (X, y, modelfun, beta0);
%! [yp, yci] = predict (mdl, [2.5; 5.5; 8.5]);
%! assert_equal (yp, [3.449741842; 8.158274281; 19.293455074], 1e-6);
%! assert_equal (yci(:,1), [3.329126146; 7.997938483; 19.121921613], 1e-5);
%! assert_equal (yci(:,2), [3.570357538; 8.318610079; 19.464988535], 1e-5);
%!test  # plotting methods, feval, and random run without error
%! mdl = fitnlm (X, y, modelfun, beta0);
%! assert_equal (feval (mdl, [2.5; 5.5]), predict (mdl, [2.5; 5.5]), 1e-12);
%! assert_equal (numel (random (mdl)), 10);
%! hf = figure ("visible", "off");
%! unwind_protect
%!   plotResiduals (mdl);
%!   plotResiduals (mdl, "fitted");
%!   plotDiagnostics (mdl);
%!   plotSlice (mdl);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Test input validation
%!error<Invalid call> fitnlm ()
%!error<Invalid call> fitnlm ([1, 2], [1; 2])
%!error<NonLinearModel: MODELFUN must be a function handle.> ...
%! fitnlm ([1;2], [1;2], "bad", [1])
%!error<NonLinearModel: unknown parameter name 'foo'.> ...
%! fitnlm ([1;2;3], [1;2;3], @(b, x) b(1) * x, 1, "foo", 1)

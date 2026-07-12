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
## @deftypefn  {statistics} {@var{ci} =} nlparci (@var{beta}, @var{resid}, @qcode{'covar'}, @var{CovB})
## @deftypefnx {statistics} {@var{ci} =} nlparci (@var{beta}, @var{resid}, @qcode{'jacobian'}, @var{J})
## @deftypefnx {statistics} {@var{ci} =} nlparci (@dots{}, @qcode{'alpha'}, @var{alpha})
##
## Confidence intervals for the coefficients of a nonlinear regression.
##
## @code{@var{ci} = nlparci (@var{beta}, @var{resid}, @qcode{'covar'},
## @var{CovB})} returns the @math{100 (1 - @var{alpha})%} confidence intervals
## for the fitted coefficients @var{beta} of a nonlinear regression, given the
## residual vector @var{resid} and the estimated coefficient covariance matrix
## @var{CovB} (both produced by @code{nlinfit}).  @var{ci} is a
## @math{p}-by-@math{2} matrix whose rows are the lower and upper bounds for the
## corresponding coefficient.
##
## @code{@var{ci} = nlparci (@var{beta}, @var{resid}, @qcode{'jacobian'},
## @var{J})} instead derives the coefficient covariance from the Jacobian
## @var{J} and the residuals.  A legacy positional form @code{nlparci
## (@var{beta}, @var{resid}, @var{J})} is also accepted.
##
## The confidence level defaults to @math{95%}; pass @code{@qcode{'alpha'},
## @var{alpha}} for a @math{100 (1 - @var{alpha})%} interval.  The intervals use
## Student's @math{t} distribution with @code{numel (@var{resid}) - numel
## (@var{beta})} degrees of freedom.
##
## @seealso{nlinfit, nlpredci, fitnlm, NonLinearModel}
## @end deftypefn

function ci = nlparci (beta, resid, varargin)

  if (nargin < 3)
    print_usage ();
  endif

  if (! (isnumeric (beta) && isvector (beta) && isreal (beta)))
    error ("nlparci: BETA must be a real numeric vector.");
  endif
  if (! (isnumeric (resid) && isreal (resid)))
    error ("nlparci: RESID must be a real numeric vector.");
  endif

  beta  = beta(:);
  resid = resid(:);
  p     = numel (beta);
  dfe   = numel (resid) - p;
  if (dfe <= 0)
    error ("nlparci: not enough residuals to estimate the coefficients.");
  endif

  ## Parse the covariance source and the optional confidence level.  The third
  ## argument may be a bare Jacobian (legacy) or a 'covar'/'jacobian' keyword.
  alpha = 0.05;
  CovB  = [];
  J     = [];
  args  = varargin;
  if (! ischar (args{1}))
    J    = args{1};
    args = args(2:end);
  endif
  if (mod (numel (args), 2) != 0)
    error ("nlparci: Name/Value arguments must come in pairs.");
  endif
  for k = 1:2:numel (args)
    switch (lower (args{k}))
      case 'covar'
        CovB = args{k+1};
      case 'jacobian'
        J = args{k+1};
      case 'alpha'
        alpha = args{k+1};
      otherwise
        error ("nlparci: unknown parameter name '%s'.", args{k});
    endswitch
  endfor

  if (! (isscalar (alpha) && isreal (alpha) && alpha > 0 && alpha < 1))
    error ("nlparci: ALPHA must be a scalar in the range (0, 1).");
  endif

  ## Standard errors from either the supplied covariance or the Jacobian.
  if (! isempty (CovB))
    se = sqrt (diag (CovB));
  elseif (! isempty (J))
    rmse = sqrt (sum (resid .^ 2) / dfe);
    se   = rmse * sqrt (diag (pinv (J' * J)));
  else
    error ("nlparci: a covariance matrix or a Jacobian is required.");
  endif

  delta = se * tinv (1 - alpha / 2, dfe);
  ci    = [beta - delta, beta + delta];

endfunction

%!demo
%! ## 95% confidence intervals for the coefficients of an exponential fit.
%! x = [1:10]';
%! y = [2.1;2.9;4.2;5.3;7.1;9.4;12.8;16.5;22.1;29.8];
%! modelfun = @(b, x) b(1) .* exp (b(2) .* x);
%! [beta, R, J, CovB] = nlinfit (x, y, modelfun, [1; 0.3]);
%! ci = nlparci (beta, R, 'covar', CovB)

%!shared beta, R, J, CovB
%! x = [1;2;3;4;5;6;7;8;9;10];
%! y = [2.1;2.9;4.2;5.3;7.1;9.4;12.8;16.5;22.1;29.8];
%! modelfun = @(b, xx) b(1) .* exp (b(2) .* xx);
%! [beta, R, J, CovB] = nlinfit (x, y, modelfun, [1; 0.3]);

## Values verified against MATLAB's nlparci.
%!test
%! ci = nlparci (beta, R, "covar", CovB);
%! assert_equal (ci, [1.602588326, 1.764905576; ...
%!                    0.281490049, 0.292332135], 1e-6);
%!test
%! ## The Jacobian and covariance forms agree.
%! ci1 = nlparci (beta, R, "covar", CovB);
%! ci2 = nlparci (beta, R, "jacobian", J);
%! assert_equal (ci1, ci2, 1e-8);
%!test
%! ## The legacy positional Jacobian form matches the named form.
%! ci = nlparci (beta, R, J);
%! assert_equal (ci, nlparci (beta, R, "jacobian", J), 1e-12);
%!test
%! ## A 90% interval is narrower than the default 95% interval.
%! ci95 = nlparci (beta, R, "covar", CovB);
%! ci90 = nlparci (beta, R, "covar", CovB, "alpha", 0.10);
%! assert (all (diff (ci90, 1, 2) < diff (ci95, 1, 2)));

## Test input validation
%!error<Invalid call> nlparci (1, 2)
%!error<nlparci: unknown parameter name 'foo'.> ...
%! nlparci ([1;2], [1;2;3], "foo", 1)
%!error<nlparci: ALPHA must be a scalar in the range> ...
%! nlparci ([1;2], [1;2;3], "covar", eye (2), "alpha", 2)
%!error<nlparci: a covariance matrix or a Jacobian is required.> ...
%! nlparci ([1;2], [1;2;3], "alpha", 0.05)
%!error<nlparci: not enough residuals> nlparci ([1;2], 1, "covar", eye (2))

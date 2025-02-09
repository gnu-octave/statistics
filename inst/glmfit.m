## Copyright (C) 2024 Ruchika Sonagote <ruchikasonagote2003@gmail.com>
## Copyright (C) 2025 Swayam Shah <swayamshah66@gmail.com>
## Copyright (C) 2024-2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{b} =} glmfit (@var{X}, @var{y}, @var{distribution})
## @deftypefnx {statistics} {@var{b} =} glmfit (@var{X}, @var{y}, @var{distribution}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {[@var{b}, @var{dev}] =} glmfit (@dots{})
## @deftypefnx {statistics} {[@var{b}, @var{dev}, @var{stats}] =} glmfit (@dots{})
##
## Perform generalized linear model fitting.
##
## @code{@var{b} = glmfit (@var{X}, @var{y}, @var{distribution})} returns a
## vector @var{b} of coefficient estimates for a generalized linear regression
## model of the responses in @var{y} on the predictors in @var{X}, using the
## distribution defined in @var{distribution}.
##
## @itemize
## @item @var{X} is an @math{nxp} numeric matrix of predictor variables with
## @math{n} observations and @math{p} predictors.
## @item @var{y} is an @math{nx1} numeric vector of responses for all supported
## distributions, except for the 'binomial' distribution in which case @var{y}
## can be either a numeric or logical @math{nx1} vector or an @math{nx2}
## matrix, where the first column contains the number of successes and the
## second column contains the number of trials.
## @item @var{distribution} is a character vector specifying the distribution of
## the response variable. Supported distributions are @qcode{"normal"},
## @qcode{"binomial"}, @qcode{"poisson"}, @qcode{"gamma"}, and @qcode{"inverse
## gaussian"}.
## @end itemize
##
## @code{@var{b} = glmfit (@dots{}, @var{Name}, @var{Value})} specifies
## additional options using @qcode{Name-Value} pair arguments.
##
## @multitable @columnfractions 0.18 0.02 0.8
## @headitem @var{Name} @tab @tab @var{Value}
##
## @item @qcode{"B0"} @tab @tab A numeric vector specifying initial values for
## the coefficient estimates.  By default, the initial values are fitted values
## fitted from the data.
##
## @item @qcode{"Constant"} @tab @tab A character vector specifying whether to
## include a constant term in the model.  Valid options are @var{"on"} (default)
## and @var{"off"}.
##
## @item @qcode{"EstDisp"} @tab @tab A character vector specifying whether to
## compute dispersion parameter. Valid options are @var{"on"} and @var{"off"}.
## For @qcode{"binomial"} and @qcode{"poisson"} distributions the default is
## @var{"off"}, whereas for the @qcode{"normal"}, @qcode{"gamma"}, and
## @qcode{"inverse gaussian"} distributions the default is @var{"on"}.
##
## @item @qcode{"link"} @tab @tab A character vector specifying the name of a
## canonical link function or a numeric scalar for specifying a @qcode{"power"}
## link function.  Supported canonical link functions include @qcode{"identity"}
## (default for @qcode{"normal"} distribution), @qcode{"log"} (default for
## @qcode{"poisson"} distribution), @qcode{"logit"} (default for
## @qcode{"binomial"} distribution), @qcode{"probit"}, @qcode{"loglog"},
## @qcode{"comploglog"}, and @qcode{"reciprocal"} (default for the
## @qcode{"gamma"} distribution).  The @qcode{"power"} link function is the
## default for the @qcode{"inverse gaussian"} distribution with @math{p = -2}.
## For custom link functions, the user can provide cell array with three
## function handles: the link function, its derivative, and its inverse, or
## alternatively a structure @var{S} with three fields: @qcode{S.Link},
## @qcode{S.Derivative}, and @qcode{S.Inverse}.  Each field can either contain a
## function handle or a character vector with the name of an existing function.
## All custom link functions must accept a vector of inputs and return a vector
## of the same size.
##
## @item @qcode{"Offset"} @tab @tab A numeric vector of the same length as the
## response @var{y} specifying an offset variable in the fit. It is used as an
## additional predictor with a coefficient value fixed at 1.
##
## @item @qcode{"Options"} @tab @tab A scalar structure containing the fields
## @qcode{MaxIter} and @qcode{TolX}.  @qcode{MaxIter} must be a scalar positive
## integer specifying the maximum number of iteration allowed for fitting the
## model, and @qcode{TolX} must be a positive scalar value specifying the
## termination tolerance.
##
## @item @qcode{"Weights"} @tab @tab An @math{nx1} numeric vector of nonnegative
## values, where @math{n} is the number of observations in @var{X}.  By default,
## it is @code{ones (n, 1)}.
## @end multitable
##
## @code{[@var{b}, @var{dev}] = glmfit (@dots{})} also returns the deviance of
## the fit as a numeric value in @var{dev}.  Deviance is a generalization of the
## residual sum of squares.  It measures the goodness of fit compared to a
## saturated model.
##
## @code{[@var{b}, @var{dev}, @var{stats}] = glmfit (@dots{})} also returns the
## structure @var{stats}, which contains the model statistics in the following
## fields:
##
## @itemize
## @item @qcode{beta} - Coefficient estimates @var{b}
## @item @qcode{dfe} - Degrees of freedom for error
## @item @qcode{sfit} - Estimated dispersion parameter
## @item @qcode{s} - Theoretical or estimated dispersion parameter
## @item @qcode{estdisp} - @code{false} when @qcode{"EstDisp"} is @qcode{"off"}
## and @code{true} when @qcode{"EstDisp"} is @qcode{"on"}
## @item @qcode{covb} - Estimated covariance matrix for @var{b}
## @item @qcode{se} - Vector of standard errors of the coefficient estimates
## @var{b}
## @item @qcode{coeffcorr} - Correlation matrix for @var{b}
## @item @qcode{t} - @math{t} statistics for @var{b}
## @item @qcode{p} - @math{p}-values for @var{b}
## @item @qcode{resid} - Vector of residuals
## @item @qcode{residp} - Vector of Pearson residuals
## @item @qcode{residd} - Vector of deviance residuals
## @item @qcode{resida} - Vector of Anscombe residuals
## @end itemize
##
## @seealso{glmval}
## @end deftypefn

function [b, dev, stats] = glmfit (X, y, distribution, varargin)

  ## Check input arguments
  if (nargin < 3)
    error ("glmfit: too few input arguments.");
  elseif (mod (nargin - 3, 2) != 0)
    error ("glmfit: Name-Value arguments must be in pairs.");
  elseif (! isnumeric (X) || isempty (X))
    error ("glmfit: X must be a numeric matrix.");
  elseif (! (isnumeric (y) || islogical (y)) || isempty (y))
    error ("glmfit: Y must be either a numeric matrix or a logical vector.");
  elseif (size (X, 1) != size (y, 1))
    error ("glmfit: X and Y must have the same number of observations.");
  elseif (! ischar (distribution))
    error ("glmfit: DISTRIBUTION must be a character vector.");
  endif

  ## Remove missing values
  xymissing = any (isnan (y), 2) | any (isnan (X), 2);
  y(xymissing) = [];
  X(xymissing,:) = [];
  [ny, cy] = size (y);
  [nx, cx] = size (X);

  ## Check y dimensions based on distribution
  if (strcmpi (distribution, 'binomial'))
    if (cy > 2)
      error (["glmfit: for a 'binomial' distribution,", ...
              " Y must be an n-by-1 or n-by-2 matrix."]);
    ## Get y and N for binomial distribution
    elseif (cy == 2)
      if (! isnumeric (y))
        error (["glmfit: n-by-2 matrix Y for 'binomial' distribution", ...
                " must be numeric."]);
      endif
      N = y(:, 2);
      y = y(:, 1) ./ N;
    else
      if (islogical (y))
        y = double (y);
      endif
      N = ones (size (y));
    endif
  else
    if (cy != 1)
      error (["glmfit: for distributions other than 'binomial',", ...
              " Y must be an n-by-1 column vector."]);
    endif
  endif

  ## Set default link, variance, and deviance functions
  ## Set defaults for estimating dispersion parameter and limiting mu
  switch (tolower (distribution))
    case "normal"
      [flink, dlink, ilink] = getlinkfunctions ("identity");
      varFun = @(mu) ones (size (mu));
      devFun = @(mu, y) (y - mu) .^ 2;
      estDisp = true;
    case "binomial"
      [flink, dlink, ilink] = getlinkfunctions ("logit");
      varFun = @(mu, N) sqrt (mu) .* sqrt (1 - mu) ./ sqrt (N);
      devFun = @(mu, y, N) 2 * N .* (y .* log ((y + (y == 0)) ./ mu) + ...
                           (1 - y) .* log ((1 - y + (y == 1)) ./ (1 - mu)));
      estDisp = false;
      muLimits = [eps, 1-eps];
    case "poisson"
      [flink, dlink, ilink] = getlinkfunctions ("log");
      varFun = @(mu) sqrt (mu);
      devFun = @(mu, y) 2 * (y .* (log ((y + (y == 0)) ./ mu)) - (y - mu));
      estDisp = false;
      muLimits = realmin;
    case "gamma"
      [flink, dlink, ilink] = getlinkfunctions ("reciprocal");
      varFun = @(mu) mu;
      devFun = @(mu, y) 2 * (-log (y ./ mu) + (y - mu) ./ mu);
      estDisp = true;
      muLimits = realmin;
    case "inverse gaussian"
      [flink, dlink, ilink] = getlinkfunctions (-2);
      varFun = @(mu) mu .^ (3 / 2);
      devFun = @(mu, y) (((y - mu) ./ mu) .^ 2) ./  y;
      estDisp = true;
      muLimits = realmin;
    otherwise
      error ("glmfit: unsupported distribution.");
  endswitch

  ## Set defaults
  B0 = [];
  constant = true;
  offset = zeros (nx, 1);
  weight = ones (nx, 1);
  MaxIter = 100;
  TolX = 1e-6;

  ## Parse extra parameters
  while (numel (varargin) > 0)
    switch (tolower (varargin {1}))

      case "b0"
        B0 = varargin {2};
        if (! (isnumeric (B0) && isequal (size (B0), size (xymissing))))
          error ("glmfit: 'B0' must be a numeric vector of the same size as Y.");
        endif
        B0(xymissing) = [];

      case "constant"
        constant = tolower (varargin {2});
        if (strcmpi (constant, "on"))
          constant = true;
        elseif (strcmpi (constant, "off"))
          constant = false;
        else
          error ("glmfit: 'Constant' should be either 'on' or 'off'.");
        endif

      case "estdisp"
        estDisp = tolower (varargin {2});
        if (strcmpi (estDisp, "on"))
          estDisp = true;
        elseif (strcmpi (estDisp, "off"))
          estDisp = false;
        else
          error ("glmfit: 'EstDisp' should be either 'on' or 'off'.");
        endif

      case "link"
        linkArg = varargin {2};
        ## Input validation is performed in private function
        [flink, dlink, ilink, errmsg] = getlinkfunctions (linkArg);
        if (! isempty (errmsg))
          error ("glmfit: %s", errmsg);
        endif

      case "options"
        options = varargin {2};
        rf = {"MaxIter", "TolX"};
        if (! (isstruct (options) && all (ismember (rf, fieldnames (options)))))
          error (["glmfit: 'Options' must be a structure containing", ...
                  " the fields 'MaxIter', and 'TolX'."]);
        endif
        MaxIter = options.MaxIter;
        TolX = options.TolX;
        if (! isscalar (MaxIter) || MaxIter <= 0 || fix (MaxIter) != MaxIter)
          error (["glmfit: 'MaxIter' in 'Options' structure must", ...
                  " be a positive integer."]);
        endif
        if (! isscalar (TolX) || TolX <= 0)
          error (["glmfit: 'TolX' in 'Options' structure must", ...
                  " be a positive scalar."]);
        endif

      case "offset"
        offset = varargin {2};
        if (! (isnumeric (offset) && isequal (size (offset), size (xymissing))))
          error (["glmfit: 'Offset' must be a numeric vector", ...
                  " of the same size as Y."]);
        endif
        offset(xymissing) = [];

      case "weights"
        weight = varargin {2};
        if (! (isnumeric (weight) && isequal (size (weight), size (xymissing))))
          error (["glmfit: 'Weights' must be a numeric vector", ...
                  " of the same size as Y."]);
        endif
        weight(xymissing) = [];

      otherwise
        error ("glmfit: unknown parameter name.");
    endswitch
    varargin (1:2) = [];
  endwhile

  ## Adjust X based on constant
  if (constant)
    X = [ones(nx, 1), X];
    cx += 1;
  endif

  ## Check X for rank deficiency
  [~, R, P] = qr (X .* weight(:, ones (1, cx)), 0);
  if (isempty (R))
    rankX = 0;
  else
    rankX = sum (abs (diag (R)) > abs (R(1)) * max (nx, cx) * eps);
  endif
  if (rankX < cx)
    warning ("glmfit: X is ill-conditioned.");
    P = P(1:rankx);
    X = X(:,P);
  else
    P = [1:cx];
  endif

  ## Adjust number of observations for zero weights
  if (any (weight == 0))
    nx = nx - sum (weight == 0);
  endif

  ## Initialize mu and eta
  if (isempty (B0))  # from y
    switch (distribution)
      case "binomial"
        mu = (N .* y + 0.5) ./ (N + 1);
      case "poisson"
        mu = y + 0.25;
      case {"gamma", "inverse gaussian"}
        mu = max (y, eps);
      otherwise
        mu = y;
    endswitch
    eta = flink (mu);
  else               # from coefficient estimates
    eta = offset + X * B0(:);
    mu = ilink (eta);
  endif

  ## Initialize coefficient vector and iterations
  numc = size (X, 2);
  Bnew = zeros (numc, 1);
  seps = sqrt (eps);
  iter = 0;
  while (iter < MaxIter)
    iter += 1;
    Bold = Bnew;

    ## Compute iteratively reweighted least squares weights
    d_eta = dlink (mu);
    if (strcmpi (distribution, "binomial"))
      IRLS = abs (d_eta) .* varFun (mu, N);
    else
      IRLS = abs (d_eta) .* varFun (mu);
    endif
    squaredWeight = sqrt (weight) ./ IRLS;

    ## Estimate coefficients
    z_off = eta + (y - mu) .* d_eta - offset;
    yweighted = (z_off) .* squaredWeight;
    Xweighted = X .* squaredWeight(:, ones (1, numc));
    [Q, R] = qr (Xweighted, 0);
    Bnew = R \ (Q' * yweighted);

    ## Compute predicted mean using current linear predictor
    eta = offset + X * Bnew;
    mu = ilink (eta);

    ## Force predicted mean within distribution support limits
    if (strcmpi (distribution, "normal"))
      mu = mu;
    elseif (strcmpi (distribution, "binomial"))
      mu = max (min (mu, muLimits(2)), muLimits(1));
    else  # for "poisson", "gamma", and "inverse gaussian" distributions
      mu = max (mu, muLimits(1));
    endif

    ## Break if TolX is reached
    if (! any (abs (Bnew - Bold) > TolX * max (seps, abs (Bold))))
      break;
    endif
  endwhile

  ## Warn if iteration limit is reached
  if (iter == MaxIter)
    warning ("glmfit: maximum number of iterations has been reached.");
  endif

  ## Return estimated coefficients
  if (rankX < cx)
    b = zeros (cx, 1);
    b(P) = Bnew;
  else
    b = Bnew;
  endif

  ## Compute deviance
  if (nargout > 1)
    if (strcmpi (distribution, "binomial"))
      devn = devFun (mu, y, N);
      dev = sum (weight .* devn);
    else
      devn = devFun (mu, y);
      dev = sum (weight .* devn);
    endif
  endif

  ## Compute stats
  if (nargout > 2)
    ## Store coefficient estimates
    stats.beta = bb;

    ## Compute degrees of freedom for error
    stats.dfe = max (nx - numc, 0);

    ## Compute estimated dispersion parameter
    if (stats.dfe > 0)
      switch (tolower (distribution))
        case "normal"
          stats.sfit = sum (weight .* (y - mu) .^ 2) / stats.dfe;
        case "binomial"
          stats.sfit = sum (weight .* (y - mu) .^ 2 ./ ...
                            (mu .* (1 - mu) ./ N)) / stats.dfe;
        case "poisson"
          stats.sfit = sum (weight .* (y - mu) .^ 2 ./ mu) / stats.dfe;
        case "gamma"
          stats.sfit = sum (weight .* ((y - mu) ./ mu) .^ 2) / stats.dfe;
        case "inverse gaussian"
          stats.sfit = sum (weight .* ((y - mu) ./ mu .^ (3 / 2)) .^ 2) / ...
                       stats.dfe;
      endswitch
    else
      stats.sfit = NaN;
    endif

    ## Store theoretical or estimated dispersion parameter
    if (estDisp)
      stats.s = stats.sfit;
      stats.estdisp = estDisp;
    else
      stats.s = 1;
      stats.estdisp = estDisp;
    endif

    ## Compute covariance matrix, standard errors, correlation matrix,
    ## t-statistic, and p-value for coefficient estimates
    if (isnan (stats.s))
      stats.covb = NaN (numel (b));
      stats.se = NaN (size (b));
      stats.coeffcorr = NaN (numel (b));
      stats.t = NaN (size (b));
      stats.p = NaN (size (b));
    else
      stats.covb = zeros (cx, cx);
      stats.se = zeros (cx, 1);
      stats.coeffcorr = zeros (cx, cx);
      stats.t = NaN (cx, 1);
      stats.p = NaN (cx, 1);
      RI = R \ eye (numc);
      C = RI * RI';
      if (estDisp)
        C = C * stats.s ^ 2;
      endif
      stats.covb(P,P) = C;
      se = sqrt (diag (C));
      se = se(:);
      stats.se(P) = se;
      C = C ./ (se * se');
      stats.coeffcorr(P,P) = C;
      stats.t(P) = b ./ se;
      if (estDisp)
        stats.p = 2 * tcdf (-abs (stats.t), dfe);
      else
        stats.p = 2 * normcdf (-abs (stats.t));
      endif
    endif

    ## Compute residuals
    stats.resid = NaN (size (xymissing));   # Vector of residuals
    stats.residp = NaN (size (xymissing));  # Vector of Pearson residuals
    stats.residd = NaN (size (xymissing));  # Vector of deviance residuals
    stats.resida = NaN (size (xymissing));  # Vector of Anscombe residuals
    if (isequal (distribution, 'binomial'))
      stats.resid(! xymissing) = (y - mu) .* N;
      stats.residp(! xymissing) = (y - mu) ./ (varFun (mu, N) + (y == mu));
    else
      stats.resid(! xymissing)  = y - mu;
      stats.residp(! xymissing) = (y - mu) ./ (varFun (mu) + (y == mu));
    end
    stats.residd(! xymissing) = sign (y - mu) .* sqrt (max (0, divn));
    switch (tolower (distribution))
      case "normal"
        stats.resida(! xymissing) = y - mu;
      case "binomial"
        a = b = 2 / 3;
        stats.resida(! xymissing) = beta(a, b) ...
                                  * (betainc (y, a, b) - betainc (mu, a, b)) ...
                                  ./ ((mu .* (1 - mu)) .^ (1 / 6) ./ sqrt (N));
      case "poisson"
        stats.resida(! xymissing) = 1.5 * ((y .^ (2 / 3) - mu .^ (2 / 3)) ...
                                           ./ mu .^ (1 / 6));
      case "gamma"
        pwr = 1 / 3;
        stats.resida(! xymissing) = 3 * (y .^ pwr - mu .^ pwr) ./ mu .^ pwr;
      case "inverse gaussian"
        stats.resida(! xymissing) = (log (y) - log (mu)) ./ mu;
    endswitch
  endif

endfunction

%!demo
%! x = [210, 230, 250, 270, 290, 310, 330, 350, 370, 390, 410, 430]';
%! n = [48, 42, 31, 34, 31, 21, 23, 23, 21, 16, 17, 21]';
%! y = [1, 2, 0, 3, 8, 8, 14, 17, 19, 15, 17, 21]';
%! b = glmfit (x, [y n], "binomial", "Link", "probit");
%! yfit = glmval (b, x, "probit", "Size", n);
%! plot (x, y./n, 'o', x, yfit ./ n, '-')

%!demo
%! load fisheriris
%! X = meas (51:end, :);
%! y = strcmp ("versicolor", species(51:end));
%! b = glmfit (X, y, "binomial", "link", "logit")

## Test output
%!test
%! load fisheriris;
%! X = meas(51:end,:);
%! y = strcmp ("versicolor", species(51:end));
%! b = glmfit (X, y, "binomial", "link", "logit");
%! assert (b, [42.6379; 2.4652; 6.6809; -9.4294; -18.2861], 1e-4);

%!test
%! X = [1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9, 9.0, 10.1]';
%! y = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4]';
%! [Bnew, dev] = glmfit (X, y, "gamma", "link", "log");
%! b_matlab = [-0.7631; 0.1113];
%! dev_matlab = 0.0111;
%! assert (Bnew, b_matlab, 0.001);
%! assert (dev, dev_matlab, 0.001);

%!test
%! X = [1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9, 9.0, 10.1]';
%! y = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4]';
%! p_input = 1;
%! [Bnew, dev] = glmfit (X, y, "inverse gaussian", "link", p_input);
%! b_matlab = [0.3813; 0.0950];
%! dev_matlab = 0.0051;
%! assert (Bnew, b_matlab, 0.001);
%! assert (dev, dev_matlab, 0.001);

## Test input validation
%!error <glmfit: too few input arguments.> glmfit ()
%!error <glmfit: too few input arguments.> glmfit (1)
%!error <glmfit: too few input arguments.> glmfit (1, 2)
%!error <glmfit: Name-Value arguments must be in pairs.> ...
%! glmfit (rand (6, 1), rand (6, 1), 'poisson', 'link')
%!error <glmfit: X must be a numeric matrix.> ...
%! glmfit ('abc', rand (6, 1), 'poisson')
%!error <glmfit: X must be a numeric matrix.> ...
%! glmfit ([], rand (6, 1), 'poisson')
%!error <glmfit: Y must be either a numeric matrix or a logical vector.> ...
%! glmfit (rand (5, 2), 'abc', 'poisson')
%!error <glmfit: Y must be either a numeric matrix or a logical vector.> ...
%! glmfit (rand (5, 2), [], 'poisson')
%!error <glmfit: X and Y must have the same number of observations.> ...
%! glmfit (rand (5, 2), rand (6, 1), 'poisson')
%!error <glmfit: DISTRIBUTION must be a character vector.> ...
%! glmfit (rand (6, 2), rand (6, 1), 3)
%!error <glmfit: DISTRIBUTION must be a character vector.> ...
%! glmfit (rand (6, 2), rand (6, 1), {'poisson'})
%!error <glmfit: for a 'binomial' distribution, Y must be an n-by-1 or n-by-2 matrix.> ...
%! glmfit (rand (5, 2), rand (5, 3), 'binomial')
%!error <glmfit: n-by-2 matrix Y for 'binomial' distribution must be numeric.> ...
%! glmfit (rand (2, 2), [true, true; false, false], 'binomial')
%!error <glmfit: for distributions other than 'binomial', Y must be an n-by-1 column vector> ...
%! glmfit (rand (5, 2), rand (5, 2), 'normal')
%!error <glmfit: unsupported distribution.> ...
%! glmfit (rand (5, 2), rand (5, 1), 'chebychev')
%!error <glmfit: 'B0' must be a numeric vector of the same size as Y.> ...
%! glmfit (rand (5, 2), rand (5, 1), 'normal', 'B0', [1; 2; 3; 4])
%!error <glmfit: 'Constant' should be either 'on' or 'off'.> ...
%! glmfit (rand (5, 2), rand (5, 1), 'normal', 'constant', 1)
%!error <glmfit: 'Constant' should be either 'on' or 'off'.> ...
%! glmfit (rand (5, 2), rand (5, 1), 'normal', 'constant', 'o')
%!error <glmfit: 'Constant' should be either 'on' or 'off'.> ...
%! glmfit (rand (5, 2), rand (5, 1), 'normal', 'constant', true)
%!error <glmfit: 'EstDisp' should be either 'on' or 'off'.> ...
%! glmfit (rand (5, 2), rand (5, 1), 'normal', 'estdisp', 1)
%!error <glmfit: 'EstDisp' should be either 'on' or 'off'.> ...
%! glmfit (rand (5, 2), rand (5, 1), 'normal', 'estdisp', 'o')
%!error <glmfit: 'EstDisp' should be either 'on' or 'off'.> ...
%! glmfit (rand (5, 2), rand (5, 1), 'normal', 'estdisp', true)
%!error <glmfit: structure with custom link functions must be a scalar.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', struct ("Link", {1, 2}))
%!error <glmfit: structure with custom link functions requires the fields 'Link', 'Derivative', and 'Inverse'.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', struct ("Link", "norminv"))
%!error <glmfit: bad 'Link' function in custom link function structure.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', struct ("Link", "some", "Derivative", @(x)x, "Inverse", "normcdf"))
%!error <glmfit: bad 'Link' function in custom link function structure.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', struct ("Link", 1, "Derivative", @(x)x, "Inverse", "normcdf"))
%!error <glmfit: custom 'Link' function must return an output of the same size as input.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', struct ("Link", @(x) [x, x], "Derivative", @(x)x, "Inverse", "normcdf"))
%!error <glmfit: invalid custom 'Link' function.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', struct ("Link", "what", "Derivative", @(x)x, "Inverse", "normcdf"))
%!error <glmfit: bad 'Derivative' function in custom link function structure.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', struct ("Link", @(x)x, "Derivative", "some", "Inverse", "normcdf"))
%!error <glmfit: bad 'Derivative' function in custom link function structure.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', struct ("Link", @(x)x, "Derivative", 1, "Inverse", "normcdf"))
%!error <glmfit: custom 'Derivative' function must return an output of the same size as input.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', struct ("Link", @(x)x, "Derivative", @(x) [x, x], "Inverse", "normcdf"))
%!error <glmfit: invalid custom 'Derivative' function.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', struct ("Link", @(x)x, "Derivative", "what", "Inverse", "normcdf"))
%!error <glmfit: bad 'Inverse' function in custom link function structure.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', struct ("Link", @(x)x, "Derivative", "normcdf", "Inverse", "some"))
%!error <glmfit: bad 'Inverse' function in custom link function structure.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', struct ("Link", @(x)x, "Derivative", "normcdf", "Inverse", 1))
%!error <glmfit: custom 'Inverse' function must return an output of the same size as input.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', struct ("Link", @(x)x, "Derivative", "normcdf", "Inverse", @(x) [x, x]))
%!error <glmfit: invalid custom 'Inverse' function.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', struct ("Link", @(x)x, "Derivative", "normcdf", "Inverse", "what"))
%!error <glmfit: cell array with custom link functions must have three elements.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', {'log'})
%!error <glmfit: cell array with custom link functions must have three elements.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', {'log', 'hijy'})
%!error <glmfit: cell array with custom link functions must have three elements.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', {1, 2, 3, 4})
%!error <glmfit: bad 'Link' function in custom link function cell array.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', {"log", "dfv", "dfgvd"})
%!error <glmfit: custom 'Link' function must return an output of the same size as input.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', {@(x) [x, x], "dfv", "dfgvd"})
%!error <glmfit: invalid custom 'Link' function.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', {@(x) what (x), "dfv", "dfgvd"})
%!error <glmfit: bad 'Derivative' function in custom link function cell array.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', {@(x) x, "dfv", "dfgvd"})
%!error <glmfit: custom 'Derivative' function must return an output of the same size as input.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', {@(x) x, @(x) [x, x], "dfgvd"})
%!error <glmfit: invalid custom 'Derivative' function.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', {@(x) x, @(x) what (x), "dfgvd"})
%!error <glmfit: bad 'Inverse' function in custom link function cell array.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', {@(x) x, @(x) x, "dfgvd"})
%!error <glmfit: custom 'Inverse' function must return an output of the same size as input.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', {@(x) x, @(x) x, @(x) [x, x]})
%!error <glmfit: invalid custom 'Inverse' function.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', {@(x) x, @(x) x, @(x) what (x)})
%!error <glmfit: numeric input for custom link function must be a finite real scalar value.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', NaN)
%!error <glmfit: numeric input for custom link function must be a finite real scalar value.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', [1, 2])
%!error <glmfit: numeric input for custom link function must be a finite real scalar value.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', [1i])
%!error <glmfit: canonical link function name must be a character vector.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', ["log"; "log1"])
%!error <glmfit: canonical link function 'somelinkfunction' is not supported.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', 'somelinkfunction')
%!error <glmfit: invalid value for custom link function.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', true)
%!error <glmfit: 'Options' must be a structure containing the fields 'MaxIter', and 'TolX'.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'options', true)
%!error <glmfit: 'Options' must be a structure containing the fields 'MaxIter', and 'TolX'.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'options', struct ("MaxIter", 100))
%!error <glmfit: 'MaxIter' in 'Options' structure must be a positive integer.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'options', struct ("MaxIter", 4.5, "TolX", 1e-6))
%!error <glmfit: 'MaxIter' in 'Options' structure must be a positive integer.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'options', struct ("MaxIter", 0, "TolX", 1e-6))
%!error <glmfit: 'MaxIter' in 'Options' structure must be a positive integer.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'options', struct ("MaxIter", -100, "TolX", 1e-6))
%!error <glmfit: 'MaxIter' in 'Options' structure must be a positive integer.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'options', struct ("MaxIter", [50 ,50], "TolX", 1e-6))
%!error <glmfit: 'TolX' in 'Options' structure must be a positive scalar.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'options', struct ("MaxIter", 100, "TolX", 0))
%!error <glmfit: 'TolX' in 'Options' structure must be a positive scalar.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'options', struct ("MaxIter", 100, "TolX", -1e-6))
%!error <glmfit: 'TolX' in 'Options' structure must be a positive scalar.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'options', struct ("MaxIter", 100, "TolX", [1e-6, 1e-6]))
%!error <glmfit: 'Offset' must be a numeric vector of the same size as Y.> ...
%! glmfit (rand (5, 2), rand (5, 1), 'normal', 'offset', [1; 2; 3; 4])
%!error <glmfit: 'Offset' must be a numeric vector of the same size as Y.> ...
%! glmfit (rand (5, 2), rand (5, 1), 'normal', 'offset', 'asdfg')
%!error <glmfit: 'Weights' must be a numeric vector of the same size as Y.> ...
%! glmfit (rand (5, 2), rand (5, 1), 'normal', 'weights', [1; 2; 3; 4])
%!error <glmfit: 'Weights' must be a numeric vector of the same size as Y.> ...
%! glmfit (rand (5, 2), rand (5, 1), 'normal', 'weights', 'asdfg')

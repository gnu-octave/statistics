## Copyright (C) 2012-2019 Fernando Damian Nieuwveldt <fdnieuwveldt@gmail.com>
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {[@var{xload}, @var{yload}] =} plsregress (@var{X}, @var{Y})
## @deftypefnx {statistics} {[@var{xload}, @var{yload}] =} plsregress (@var{X}, @var{Y}, @var{NCOMP})
## @deftypefnx {statistics} {[@var{xload}, @var{yload}, @var{xscore}, @var{yscore}, @var{coef}, @var{pctVar}, @var{mse}, @var{stats}] =} plsregress (@var{X}, @var{Y}, @var{NCOMP})
## @deftypefnx {statistics} {[@var{xload}, @var{yload}, @var{xscore}, @var{yscore}, @var{coef}, @var{pctVar}, @var{mse}, @var{stats}] =} plsregress (@dots{}, @var{Name}, @var{Value})
##
## Calculate partial least squares regression using SIMPLS algorithm.
##
## @code{plsregress} uses the SIMPLS algorithm, and first centers @var{X} and
## @var{Y} by subtracting off column means to get centered variables, unless
## @qcode{'Intercept'} is @qcode{false}.  However, it does not rescale the
## columns.  To perform partial least squares regression with standardized
## variables, use @code{zscore} to normalize @var{X} and @var{Y}.
##
## @code{[@var{xload}, @var{yload}] = plsregress (@var{X}, @var{Y})} computes a
## partial least squares regression of  @var{Y} on @var{X}, using @var{NCOMP}
## PLS components, which by default are calculated as
## @qcode{min (size (@var{X}, 1) - 1, size(@var{X}, 2))}, and returns the
## the predictor and response loadings in @var{xload} and @var{yload},
## respectively.
## @itemize
## @item @var{X} is an @math{N*P} matrix of predictor variables, with rows
## corresponding to observations, and columns corresponding to variables.
## @item @var{Y} is an @math{N*M} response matrix.
## @item @var{xload} is a @math{P*NCOMP} matrix of predictor loadings, where
## each row of @var{xload} contains coefficients that define a linear
## combination of PLS components that approximate the original predictor
## variables.
## @item @var{yload} is an @math{M*NCOMP} matrix of response loadings, where
## each row of @var{yload} contains coefficients that define a linear
## combination of PLS components that approximate the original response
## variables.
## @end itemize
##
## @code{[@var{xload}, @var{yload}] = plsregress (@var{X}, @var{Y},
## @var{NCOMP})} defines the desired number of PLS components to use in the
## regression.  @var{NCOMP}, a scalar positive integer, must not exceed the
## default calculated value.
##
## @code{[@var{xload}, @var{yload}, @var{xscore}, @var{yscore}, @var{coef},
## @var{pctVar}, @var{mse}, @var{stats}] = plsregress (@var{X}, @var{Y},
## @var{NCOMP})} also returns the following arguments:
## @itemize
## @item @var{xscore} is an @math{N*NCOMP} orthonormal matrix with the predictor
## scores, i.e., the PLS components that are linear combinations of the
## variables in @var{X}, with rows corresponding to observations and columns
## corresponding to components.
## @item @var{yscore} is an @math{N*NCOMP} orthonormal matrix with the response
## scores, i.e., the linear combinations of the responses with which the PLS
## components @var{xscore} have maximum covariance, with rows corresponding to
## observations and columns corresponding to components.
## @item @var{coef} is a @math{(P+1)*M} matrix with the PLS regression
## coefficients, containing the intercepts in the first row.  If
## @qcode{'Intercept'} is @qcode{false}, @var{coef} is a @math{P*M} matrix
## without the intercepts.
## @item @var{pctVar} is a @math{2*NCOMP} matrix containing the percentage of
## the variance explained by the model with the first row containing the
## percentage of explained varianced in @var{X} by each PLS component and the
## second row containing the percentage of explained variance in @var{Y}.
## @item @var{mse} is a @math{2*(NCOMP+1)} matrix containing the estimated mean
## squared errors for PLS models with @qcode{0:@var{NCOMP}} components with the
## first row containing the squared errors for the predictor variables in
## @var{X} and the second row containing the mean squared errors for the
## response variable(s) in @var{Y}.
## @item @var{stats} is a structure with the following fields:
## @itemize
## @item @var{stats}@qcode{.W} is a @math{P*NCOMP} matrix of PLS weights.
## @item @var{stats}@qcode{.T2} is the @math{T^2} statistics for each point in
## @var{xscore}.
## @item @var{stats}@qcode{.Xresiduals} is an @math{N*P} matrix with the
## predictor residuals.
## @item @var{stats}@qcode{.Yresiduals} is an @math{N*M} matrix with the
## response residuals.
## @end itemize
## @end itemize
##
## @code{[@dots{}] = plsregress (@dots{}, @var{Name}, @var{Value}, @dots{})}
## specifies one or more of the following @var{Name}/@var{Value} pairs:
##
## @multitable @columnfractions 0.2 0.75
## @headitem @var{Name} @tab @var{Value}
## @item @qcode{'CV'} @tab The method used to compute @var{mse}.  When
## @var{Value} is a positive integer @math{K}, @code{plsregress} uses
## @math{K}-fold cross-validation.  Set @var{Value} to a cross-validation
## partition, created using @code{cvpartition}, to use other forms of
## cross-validation.  Set @var{Value} to @qcode{'resubstitution'} to use both
## @var{X} and @var{Y} to fit the model and to estimate the mean squared errors,
## without cross-validation. By default, @qcode{@var{Value} = "resubstitution"}.
## @item @qcode{'MCReps'} @tab A positive integer indicating the number of
## Monte-Carlo repetitions for cross-validation.  By default,
## @qcode{@var{Value} = 1}.  A different @qcode{'MCReps'} value is only
## meaningful when using the @qcode{'HoldOut'} method for cross-validation,
## previously set by a @code{cvpartition} object.  If no cross-validation method
## is used, then @qcode{'MCReps'} must be @qcode{1}.
## @item @qcode{'Intercept'} @tab A logical scalar, or @qcode{'on'} or
## @qcode{'off'}, specifying whether the model includes a constant term.  When
## the constant term is excluded, @var{X} and @var{Y} are not centered and
## @var{coef} contains no intercepts.  By default, @qcode{@var{Value} = true}.
## @end multitable
##
## Further information about the PLS regression can be found at
## @url{https://en.wikipedia.org/wiki/Partial_least_squares_regression}
##
## @subheading References
## @enumerate
## @item
## SIMPLS: An alternative approach to partial least squares regression.
## Chemometrics and Intelligent Laboratory Systems (1993)
##
## @end enumerate
## @end deftypefn

function [xload, yload, xscore, yscore, coef, pctVar, mse, stats] = ...
                                              plsregress (X, Y, NCOMP, varargin)

  ## Check input arguments and add defaults
  if (nargin < 2)
    error ("plsregress: function called with too few input arguments.");
  endif

  if (! isnumeric (X) || ! isnumeric (Y))
    error ("plsregress: X and Y must be real matrices.");
  endif

  ## Get size of predictor and response inputs
  [nobs, npred] = size (X);
  [Yobs, nresp] = size (Y);

  if (nobs != Yobs)
    error ("plsregress: X and Y observations mismatch.");
  endif

  ## Calculate max number of components
  NCOMPmax = min (nobs - 1, npred);

  if (nargin < 3)
    NCOMP = NCOMPmax;
  elseif (! isnumeric (NCOMP) || ! isscalar (NCOMP) ||
          NCOMP != fix (NCOMP) || NCOMP <= 0)
    error ("plsregress: invalid value for NCOMP.");
  elseif (NCOMP > NCOMPmax)
    error ("plsregress: NCOMP exceeds maximum components for X.");
  endif

  ## Add default optional arguments
  CV = false;
  mcreps = 1;
  intercept = true;

  ## Parse additional Name-Value pairs
  while (numel (varargin) > 0)
    if (strcmpi (varargin{1}, 'cv'))
      cvarg = varargin{2};
      if (isa (cvarg, 'cvpartition'))
        CV = true;
      elseif (isscalar (cvarg) && cvarg == fix (cvarg) && cvarg > 0)
        CV = true;
      elseif (! strcmpi (cvarg, 'resubstitution'))
        error ("plsregress: invalid VALUE for 'cv' optional argument.");
      endif
    elseif (strcmpi (varargin{1}, 'mcreps'))
      mcreps = varargin{2};
      if (! (isscalar (mcreps) && mcreps == fix (mcreps) && mcreps > 0))
        error ("plsregress: invalid VALUE for 'mcreps' optional argument.");
      endif
    elseif (strcmpi (varargin{1}, 'intercept'))
      intarg = varargin{2};
      if (ischar (intarg) && isrow (intarg) &&
          any (strcmpi (intarg, {'on', 'off'})))
        intercept = strcmpi (intarg, 'on');
      elseif ((isnumeric (intarg) || islogical (intarg)) && isscalar (intarg) &&
              (intarg == 0 || intarg == 1))
        intercept = logical (intarg);
      else
        error ("plsregress: invalid VALUE for 'intercept' optional argument.");
      endif
    else
      error ("plsregress: invalid NAME argument.");
    endif
    varargin(1:2) = [];
  endwhile

  ## Check MCREPS = 1 when "resubstitution" is set for cross validation
  if (! CV && mcreps != 1)
    error (strcat ("plsregress: 'mcreps' must be 1 when 'resubstitution'", ...
                   " is specified for cross validation."));
  endif


  ## Check number of output arguments
  if (nargout > 8)
    print_usage ();
  endif

  ## Mean centering Data matrix and responses, unless intercept is excluded
  if (intercept)
    Xmeans = mean (X);
    Ymeans = mean (Y);
  else
    Xmeans = zeros (1, npred);
    Ymeans = zeros (1, nresp);
  endif
  X0 = bsxfun (@minus, X, Xmeans);
  Y0 = bsxfun (@minus, Y, Ymeans);

  [P, Q, T, U, W] = simpls (X0, Y0, NCOMP, intercept);

  ## Store output arguments
  xload  = P;
  yload  = Q;
  xscore = T;
  yscore = U;

  ## Compute regression coefficients
  if (nargout > 4)
    coef = W * Q';
    if (intercept)
      coef = [Ymeans - Xmeans * coef; coef];
    endif
  endif

  ## Compute the percent of variance explained for X and Y
  if (nargout > 5)
    XVar = sum (abs (xload) .^ 2, 1) ./ sum (sum (abs (X0) .^ 2, 1));
    YVar = sum (abs (yload) .^ 2, 1) ./ sum (sum (abs (Y0) .^ 2, 1));
    pctVar = [XVar; YVar];
  endif

  ## Estimate the mean squared errors
  if (nargout > 6)
    ## Compute MSE by cross-validation
    if (CV)
      mse = NaN (2, NCOMP + 1);
      ## Check crossval method and recalculate max number of components
      if isa (cvarg, 'cvpartition')
        type = 'Partition';
        NCOMPmax = min (min (cvarg.TrainSize)-1,npred);
        ts = sum (cvarg.TestSize);
      else
        type = 'Kfold';
        NCOMPmax = min (floor ((nobs * (cvarg - 1) / cvarg) -1), npred);
        ts = nobs;
      endif
      if (NCOMP > NCOMPmax)
        warning (strcat ("plsregress: NCOMP exceeds maximum components", ...
                         " for cross validation."));
        NCOMP = NCOMPmax;
      endif
      ## Create function handle with NCOMP extra argument
      F = @(xtr, ytr, xte, yte) sseCV (xtr, ytr, xte, yte, NCOMP, intercept);
      ## Apply cross validation
      sse = crossval (F, X, Y, type, cvarg, 'mcreps', mcreps);
      ## Compute MSE from the SSEs collected from each cross validation set
      mse(:,1:NCOMP+1) = reshape (sum (sse, 1) / (ts * mcreps), [2, NCOMP+1]);

      ## Computed fitted if residuals are requested
      if (nargout > 7)
        xfitted = xscore * xload';
        yfitted = xscore * yload';
      endif

    ## Compute MSE by resubstitution
    else
      mse = zeros (2, NCOMP + 1);
      ## Model with 0 components
      mse(1,1) = sum (sum (abs (X0) .^ 2, 2));
      mse(2,1) = sum (sum (abs (Y0) .^ 2, 2));
      ## Models with 1:NCOMP components
      for i = 1:NCOMP
        xfitted = xscore(:,1:i) * xload(:,1:i)';
        yfitted = xscore(:,1:i) * yload(:,1:i)';
        mse(1,i+1) = sum (sum (abs (X0 - xfitted) .^ 2, 2));
        mse(2,i+1) = sum (sum (abs (Y0 - yfitted) .^ 2, 2));
      endfor
      ## Compute the mean of the sum of squares above
      mse = mse / nobs;
    endif
  endif

  ## Compute stats
  if (nargout > 7)
    ## Save weights
    stats.W = W;
    ## Compute T-squared
    stats.T2 = sum (bsxfun (@rdivide, abs (xscore) .^ 2, ...
                            var (xscore, [], 1)) , 2);
    ## Compute residuals for X and Y
    stats.Xresiduals = X0 - xfitted;
    stats.Yresiduals = Y0 - yfitted;
  endif

endfunction

## SIMPLS algorithm
function [P, Q, T, U, W] = simpls (X0, Y0, NCOMP, intercept)
  ## Get size of predictor and response inputs
  [nobs, npred] = size (X0);
  [Yobs, nresp] = size (Y0);
  ## Compute covariance
  S = X0' * Y0;
  ## Preallocate matrices
  W = P = V = zeros (npred, NCOMP);
  T = U = zeros (nobs, NCOMP);
  Q = zeros (nresp, NCOMP);
  ## Models with 1:NCOMP components
  for a = 1:NCOMP
    [eigvec, eigval] = eig (S' * S);      # Y factor weights
    ## Get max eigenvector, the first one when several tie
    [~, domindex] = max (diag (eigval));
    q  = eigvec(:,domindex);
    w  = S * q;                           # X block factor weights
    t  = X0 * w;                          # X block factor scores
    if (intercept)
      t  = t - mean (t);
    endif
    nt = sqrt (t' * t);                   # compute norm
    t  = t ./ nt;
    w  = w ./ nt;                         # normalize
    p  = X0' * t;                         # X block factor loadings
    q  = Y0' * t;                         # Y block factor loadings
    u  = Y0 * q;                          # Y block factor scores
    v  = p;
    ## Ensure orthogonality
    if (a > 1)
      v = v - V * (V' * p);
      u = u - T * (T' * u);
    endif
    v = v ./ sqrt (v' * v);   # normalize orthogonal loadings
    S = S - v * (v' * S);     # deflate S wrt loadings
    V(:,a) = v;
    ## Store data
    P(:,a) = p;     # Xloads
    Q(:,a) = q;     # Yloads
    T(:,a) = t;     # xscore
    U(:,a) = u;     # Yscores
    W(:,a) = w;     # Weights
  endfor
endfunction

## Helper function for SSE cross-validation
function sse = sseCV (XTR, YTR, XTE, YTE, NCOMP, intercept)
  ## Center train data, unless intercept is excluded
  if (intercept)
    XTRmeans = mean (XTR);
    YTRmeans = mean (YTR);
  else
    XTRmeans = zeros (1, columns (XTR));
    YTRmeans = zeros (1, columns (YTR));
  endif
  X0TR = bsxfun (@minus, XTR, XTRmeans);
  Y0TR = bsxfun (@minus, YTR, YTRmeans);
  ## Center test data
  X0TE = bsxfun (@minus, XTE, XTRmeans);
  Y0TE = bsxfun (@minus, YTE, YTRmeans);
  ## Fit the full model
  [xload, yload, ~, ~, W] = simpls (X0TR, Y0TR, NCOMP, intercept);
  XTEscore = X0TE * W;
  ## Preallocate SSE matrix
  sse = zeros (2, NCOMP + 1);
  ## Model with 0 components
  sse(1,1) = sum (sum (abs (X0TE) .^ 2, 2));
  sse(2,1) = sum (sum (abs (Y0TE) .^ 2, 2));
  ## Models with 1:NCOMP components
  for i = 1:NCOMP
    X0fitted = XTEscore(:,1:i) * xload(:,1:i)';
    sse(1,i+1) = sum (sum (abs (X0TE - X0fitted) .^ 2, 2));
    Y0fitted = XTEscore(:,1:i) * yload(:,1:i)';
    sse(2,i+1) = sum (sum (abs (Y0TE - Y0fitted) .^ 2, 2));
  endfor
  ## crossval collects a row per test set, and the caller reshapes the summed
  ## rows back to [2, NCOMP+1]
  sse = sse(:)';
endfunction

%!demo
%! ## Perform Partial Least-Squares Regression
%!
%! ## Load the spectra data set and use the near infrared (NIR) spectral
%! ## intensities (NIR) as the predictor and the corresponding octave
%! ## ratings (octave) as the response.
%! load spectra
%!
%! ## Perform PLS regression with 10 components
%! [xload, yload, xscore, yscore, coef, ptcVar] = plsregress (NIR, octane, 10);
%!
%! ## Plot the percentage of explained variance in the response variable
%! ## (PCTVAR) as a function of the number of components.
%! plot (1:10, cumsum (100 * ptcVar(2,:)), '-ro');
%! xlim ([1, 10]);
%! xlabel ('Number of PLS components');
%! ylabel ('Percentage of Explained Variance in octane');
%! title ('Explained Variance per PLS components');
%!
%! ## Compute the fitted response and display the residuals.
%! octane_fitted = [ones(size(NIR,1),1), NIR] * coef;
%! residuals = octane - octane_fitted;
%! figure
%! stem (residuals, 'color', 'r', 'markersize', 4, 'markeredgecolor', 'r')
%! xlabel ('Observations');
%! ylabel ('Residuals');
%! title ('Residuals in octane''s fitted response');

%!demo
%! ## Calculate Variable Importance in Projection (VIP) for PLS Regression
%!
%! ## Load the spectra data set and use the near infrared (NIR) spectral
%! ## intensities (NIR) as the predictor and the corresponding octave
%! ## ratings (octave) as the response.  Variables with a VIP score greater than
%! ## 1 are considered important for the projection of the PLS regression model.
%! load spectra
%!
%! ## Perform PLS regression with 10 components
%! [xload, yload, xscore, yscore, coef, pctVar, mse, stats] = ...
%!                                                 plsregress (NIR, octane, 10);
%!
%! ## Calculate the normalized PLS weights
%! W0 = stats.W ./ sqrt (sum (stats.W.^2,1));
%!
%! ## Calculate the VIP scores for 10 components
%! nobs = size (xload, 1);
%! SS = sum (xscore .^ 2, 1) .* sum (yload .^ 2, 1);
%! VIPscore = sqrt (nobs * sum (SS .* (W0 .^ 2), 2) ./ sum (SS, 2));
%!
%! ## Find variables with a VIP score greater than or equal to 1
%! VIPidx = find (VIPscore >= 1);
%!
%! ## Plot the VIP scores
%! scatter (1:length (VIPscore), VIPscore, 'xb');
%! hold on
%! scatter (VIPidx, VIPscore(VIPidx), 'xr');
%! plot ([1, length(VIPscore)], [1, 1], '--k');
%! hold off
%! axis ('tight');
%! xlabel ('Predictor Variables');
%! ylabel ('VIP scores');
%! title ('VIP scores for each predictor variable with 10 components');

## Test output
%!test
%! load spectra
%! [xload, yload, xscore, yscore, coef, pctVar] = plsregress (NIR, octane, 10);
%! xload1_out = [-0.0170, 0.0039, 0.0095,  0.0258, 0.0025, ...
%!               -0.0075, 0.0000, 0.0018, -0.0027, 0.0020];
%! yload_out = [6.6384, 9.3106, 2.0505, 0.6471, 0.9625, ...
%!              0.5905, 0.4244, 0.2437, 0.3516, 0.2548];
%! xscore1_out = [-0.0401, -0.1764, -0.0340, 0.1669,  0.1041, ...
%!                -0.2067,  0.0457,  0.1565, 0.0706, -0.1471];
%! yscore1_out = [-12.4635, -15.0003,  0.0638,  0.0652, -0.0070, ...
%!                 -0.0634,   0.0062, -0.0012, -0.0151, -0.0173];
%! assert_equal (xload(1,:), xload1_out, 1e-4);
%! assert_equal (yload, yload_out, 1e-4);
%! assert_equal (xscore(1,:), xscore1_out, 1e-4);
%! assert_equal (yscore(1,:), yscore1_out, 1e-4);
%!test
%! load spectra
%! [xload, yload, xscore, yscore, coef, pctVar] = plsregress (NIR, octane, 5);
%! xload1_out = [-0.0170, 0.0039, 0.0095, 0.0258, 0.0025];
%! yload_out = [6.6384, 9.3106, 2.0505, 0.6471, 0.9625];
%! xscore1_out = [-0.0401, -0.1764, -0.0340, 0.1669, 0.1041];
%! yscore1_out = [-12.4635, -15.0003, 0.0638, 0.0652, -0.0070];
%! assert_equal (xload(1,:), xload1_out, 1e-4);
%! assert_equal (yload, yload_out, 1e-4);
%! assert_equal (xscore(1,:), xscore1_out, 1e-4);
%! assert_equal (yscore(1,:), yscore1_out, 1e-4);

## Test Intercept option
%!shared Xi, Yi
%! Xi = [1 2 4; 2 1 3; 3 5 2; 4 3 6; 5 6 5; 6 4 9; 7 8 7; 8 7 10];
%! Yi = [3 1; 2 4; 6 2; 5 7; 9 5; 8 9; 12 8; 11 11];
%!test
%! [xload, yload, xscore, yscore, coef] = plsregress (Xi, Yi, 2, ...
%!                                                    'Intercept', false);
%! xload_out = [-14.2220866289987, -0.382908357308142; ...
%!              -13.7943470180713, -3.6736441811456; ...
%!              -17.5828395894702, 3.22646836023154];
%! yload_out = [-21.7155246617638, -3.37516200500198; ...
%!              -18.5549450196154, 2.92602803399922];
%! xscore_out = [-0.159463240281499, 0.176166759987928; ...
%!               -0.135201662297805, 0.219123900595871; ...
%!               -0.209658303362604, -0.554109815239913; ...
%!               -0.290146968405309, 0.275009629902996; ...
%!               -0.344859965660409, -0.334985914644042; ...
%!               -0.425348630703114, 0.494133530498867; ...
%!               -0.474906494400693, -0.403685011986161; ...
%!               -0.550240025885877, 0.137611435218759];
%! yscore_out = [-83.7015190049069, -4.16947710556566; ...
%!               -117.650829401989, 7.52277173856108; ...
%!               -167.403038009814, -10.4151724456727; ...
%!               -238.462238446127, 9.11950487145006; ...
%!               -288.214447053951, -9.19359074611281; ...
%!               -340.718702470649, 7.41505853701072; ...
%!               -409.025856098089, -8.06996238088916; ...
%!               -442.975166495171, 5.51470553389892];
%! coef_out = [0.42805469436003, 0.390758452559943; ...
%!             0.979698857562371, -0.111297240033807; ...
%!             0.120196040659335, 0.826534140491024];
%! assert_equal (xload, xload_out, 1e-10);
%! assert_equal (yload, yload_out, 1e-10);
%! assert_equal (xscore, xscore_out, 1e-10);
%! assert_equal (yscore, yscore_out, 1e-10);
%! assert_equal (coef, coef_out, 1e-10);
%!test
%! [~, ~, ~, ~, ~, pctVar, mse] = plsregress (Xi, Yi, 2, 'Intercept', false);
%! pctVar_out = [0.963884623305986, 0.0330389814001876; ...
%!               0.965502953865825, 0.0236134421488264];
%! mse_out = [91, 3.28649927915526, 0.27995197173819; ...
%!            105.625, 3.64375049792233, 1.14958067095253];
%! assert_equal (pctVar, pctVar_out, 1e-10);
%! assert_equal (mse, mse_out, 1e-10);
%!test
%! [~, ~, ~, ~, ~, ~, ~, stats] = plsregress (Xi, Yi, 2, 'Intercept', false);
%! W_out = [-0.0203808431932239, 0.00430379596844149; ...
%!          -0.0197436438096987, -0.163238171288745; ...
%!          -0.0248987773672195, 0.124584826649244];
%! T2_out = [1.31872433173103; 1.12791712302996; 4.05334038679012; ...
%!           4.17604263860742; 5.93709890387945; 9.54610820786475; ...
%!           10.9102185426489; 13.247290390847];
%! Xres_out = [-1.20044429274505, 0.447482719671833, 0.627786948487707; ...
%!             0.161054619240101, -0.0600354047820661, -0.0842254726551177; ...
%!             -0.193951832033687, 0.0722983097243972, 0.101429470340183; ...
%!             -0.021191834171586, 0.00789955822793997, 0.0110825275173898; ...
%!             -0.0328972127935856, 0.0122629049423306, 0.0172039976850655; ...
%!             0.139862785068515, -0.0521358465541262, -0.0731429451377288; ...
%!             0.0912643311496959, -0.0340200801975641, -0.0477277923719566; ...
%!             0.227131253714974, -0.0846664120518552, -0.118781052596658];
%! Yres_out = [0.131763427877889, -2.47430053445642; ...
%!             -0.196396368260435, 0.85017791344085; ...
%!             -0.423050452210615, -0.268857438492841; ...
%!             -0.372491594083243, 0.811653066915507; ...
%!             0.380553179528951, -0.41867952505199; ...
%!             0.431112037656321, -0.338169019643642; ...
%!             0.32465399436548, 0.369329768913692; ...
%!             -0.484289964410188, 0.387671654846869];
%! assert_equal (stats.W, W_out, 1e-10);
%! assert_equal (stats.T2, T2_out, 1e-10);
%! assert_equal (stats.Xresiduals, Xres_out, 1e-10);
%! assert_equal (stats.Yresiduals, Yres_out, 1e-10);
%!test
%! [~, ~, ~, ~, coef] = plsregress (Xi, Yi, 2);
%! coef_out = [-0.0955463955587028, -0.810241454727286; ...
%!             0.502555885912014, 0.483878040247745; ...
%!             0.950692448939054, -0.0936986428112594; ...
%!             0.0966832849963297, 0.857292898480539];
%! assert_equal (coef, coef_out, 1e-10);
%!test
%! Xc = Xi - mean (Xi);
%! Yc = Yi - mean (Yi);
%! [xl1, yl1, xs1, ys1, c1] = plsregress (Xc, Yc, 2, 'Intercept', false);
%! [xl2, yl2, xs2, ys2, c2] = plsregress (Xc, Yc, 2);
%! assert_equal (xl1, xl2, 1e-10);
%! assert_equal (yl1, yl2, 1e-10);
%! assert_equal (xs1, xs2, 1e-10);
%! assert_equal (ys1, ys2, 1e-10);
%! assert_equal (c1, c2(2:end,:), 1e-10);
%! assert_equal (c2(1,:), [0, 0], 1e-10);
%!test
%! [~, ~, ~, ~, coef] = plsregress (Xi, Yi, 2, 'Intercept', 'OFF');
%! assert_equal (size (coef), [3, 2]);
%! [~, ~, ~, ~, coef] = plsregress (Xi, Yi, 2, 'Intercept', 'On');
%! assert_equal (size (coef), [4, 2]);
%!test
%! [~, ~, ~, ~, coef] = plsregress (Xi, Yi, 2, 'Intercept', int8 (0));
%! assert_equal (size (coef), [3, 2]);
%! [~, ~, ~, ~, coef] = plsregress (Xi, Yi, 2, 'Intercept', uint8 (1));
%! assert_equal (size (coef), [4, 2]);
%!test
%! [~, ~, ~, ~, ~, ~, mse] = plsregress (Xi, Yi, 2, 'Intercept', false, ...
%!                                       'CV', 8);
%! mse_out = [91, 4.32192776866963, 0.340549626995997; ...
%!            105.625, 4.71960207316224, 1.60619414427168];
%! assert_equal (mse, mse_out, 1e-10);
%!test
%! [~, ~, ~, ~, ~, ~, mse] = plsregress (Xi, Yi, 2, 'CV', 8);
%! mse_out = [22.7755102040816, 5.53413798133651, 0.437067443370619; ...
%!            28.8775510204082, 6.35600933981599, 2.45515345309489];
%! assert_equal (mse, mse_out, 1e-10);
%!warning<plsregress: NCOMP exceeds maximum components for cross validation.>
%! [~, ~, ~, ~, ~, ~, mse] = plsregress (Xi(1:4,:), Yi(1:4,:), 3, ...
%!                                       'Intercept', false, 'CV', 4);
%! mse_out = [33.5, 6.11816888105823, 4.40871978441352; ...
%!            36, 8.82425775321653, 7.84649635299314];
%! assert_equal (mse(:,1:3), mse_out, 1e-10);
%! assert_equal (isnan (mse(:,4)), [true; true]);
%!test
%! [~, yload, ~, ~, coef, ~, ~, stats] = plsregress (Xi, Yi(:,1), 2, ...
%!                                                   'Intercept', false);
%! yload_out = [21.7303444547984, 3.27301621090877];
%! coef_out = [0.490175387381509; 0.948569023203085; 0.0943735846723177];
%! Yres_out = [0.235192227523052; -0.212040551983055; -0.402118447504584; ...
%!             -0.372650127169194; 0.385841000512361; 0.415309320847751; ...
%!             0.319605009998538; -0.505122108196838];
%! assert_equal (yload, yload_out, 1e-10);
%! assert_equal (coef, coef_out, 1e-10);
%! assert_equal (stats.Yresiduals, Yres_out, 1e-10);
%!test
%! Y = [1:6; (1:6) .^ 2]';
%! [~, ~, ~, ~, b, pv, mse, st] = plsregress (ones (6, 2), Y, 1);
%! assert_equal (b, NaN (3, 2));
%! assert_equal (pv, [NaN; NaN]);
%! assert_equal (mse, [0, NaN; 152.055555555556, NaN], 1e-12);
%! assert_equal (st.T2, NaN (6, 1));
%!test
%! X = [ones(6, 1), (1:6)'];
%! Y = [1:6; (1:6) .^ 2]';
%! [xl, yl, ~, ~, b, pv, mse, st] = plsregress (X, Y, 2);
%! assert_equal (abs (xl), [0, NaN; 4.18330013267038, NaN], 1e-12);
%! assert_equal (abs (yl), [4.18330013267038, NaN; 29.2831009286927, NaN], ...
%!               1e-12);
%! assert_equal (b, NaN (3, 2));
%! assert_equal (pv, [1, NaN; 0.959079283887468, NaN], 1e-12);
%! assert_equal (mse, [2.91666666666667, 0, NaN; ...
%!                     152.055555555556, 6.22222222222222, NaN], 1e-12);
%! assert_equal (st.T2, NaN (6, 1));

## Test input validation
%!error<plsregress: function called with too few input arguments.>
%! plsregress (1)
%!error<plsregress: X and Y must be real matrices.> plsregress (1, 'asd')
%!error<plsregress: X and Y must be real matrices.> plsregress (1, {1,2,3})
%!error<plsregress: X and Y must be real matrices.> plsregress ('asd', 1)
%!error<plsregress: X and Y must be real matrices.> plsregress ({1,2,3}, 1)
%!error<plsregress: X and Y observations mismatch.> ...
%! plsregress (ones (20,3), ones (15,1))
%!error<plsregress: invalid value for NCOMP.> ...
%! plsregress (ones (20,3), ones (20,1), 0)
%!error<plsregress: invalid value for NCOMP.> ...
%! plsregress (ones (20,3), ones (20,1), -5)
%!error<plsregress: invalid value for NCOMP.> ...
%! plsregress (ones (20,3), ones (20,1), 3.2)
%!error<plsregress: invalid value for NCOMP.> ...
%! plsregress (ones (20,3), ones (20,1), [2, 3])
%!error<plsregress: NCOMP exceeds maximum components for X.> ...
%! plsregress (ones (20,3), ones (20,1), 4)
%!error<plsregress: invalid VALUE for 'cv' optional argument.> ...
%! plsregress (ones (20,3), ones (20,1), 3, 'cv', 4.5)
%!error<plsregress: invalid VALUE for 'cv' optional argument.> ...
%! plsregress (ones (20,3), ones (20,1), 3, 'cv', -1)
%!error<plsregress: invalid VALUE for 'cv' optional argument.> ...
%! plsregress (ones (20,3), ones (20,1), 3, 'cv', 'somestring')
%!error<plsregress: invalid VALUE for 'mcreps' optional argument.> ...
%! plsregress (ones (20,3), ones (20,1), 3, 'cv', 3, 'mcreps', 2.2)
%!error<plsregress: invalid VALUE for 'mcreps' optional argument.> ...
%! plsregress (ones (20,3), ones (20,1), 3, 'cv', 3, 'mcreps', -2)
%!error<plsregress: invalid VALUE for 'mcreps' optional argument.> ...
%! plsregress (ones (20,3), ones (20,1), 3, 'cv', 3, 'mcreps', [1, 2])
%!error<plsregress: invalid NAME argument.> ...
%! plsregress (ones (20,3), ones (20,1), 3, 'Name', 3, 'mcreps', 1)
%!error<plsregress: invalid NAME argument.> ...
%! plsregress (ones (20,3), ones (20,1), 3, 'cv', 3, 'Name', 1)
%!error<plsregress: 'mcreps' must be 1 when 'resubstitution' is specified> ...
%! plsregress (ones (20,3), ones (20,1), 3, 'mcreps', 2)
%!error<plsregress: 'mcreps' must be 1 when 'resubstitution' is specified> ...
%! plsregress (ones (20,3), ones (20,1), 3, 'cv', 'resubstitution', 'mcreps', 2)
%!error<plsregress: invalid VALUE for 'intercept' optional argument.> ...
%! plsregress (ones (20,3), ones (20,1), 3, 'Intercept', 2)
%!error<plsregress: invalid VALUE for 'intercept' optional argument.> ...
%! plsregress (ones (20,3), ones (20,1), 3, 'Intercept', 'true')
%!error<plsregress: invalid VALUE for 'intercept' optional argument.> ...
%! plsregress (ones (20,3), ones (20,1), 3, 'Intercept', [true, false])
%!error<plsregress: invalid VALUE for 'intercept' optional argument.> ...
%! plsregress (ones (20,3), ones (20,1), 3, 'Intercept', {false})
%!error<plsregress: invalid VALUE for 'intercept' optional argument.> ...
%! plsregress (ones (20,3), ones (20,1), 3, 'Intercept', ['off'; 'off'])
%!test
%! ## A single output returns the predictor loadings, as MATLAB allows;
%! ## fewer than two outputs used to be refused outright.
%! X = [1, 2, 3; 2, 3, 5; 3, 5, 8; 4, 7, 11; 5, 11, 16; 6, 13, 19];
%! Y = [1; 2; 3; 4; 5; 6];
%! xl = plsregress (X, Y, 2);
%! [xl2, yl2] = plsregress (X, Y, 2);
%! assert_equal (xl, xl2);
%! assert_equal (size (xl), [3, 2]);

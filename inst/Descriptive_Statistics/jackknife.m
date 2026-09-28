## Copyright (C) 2011 Alexander Klein <alexander.klein@math.uni-giessen.de>
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
## @deftypefn  {statistics} {@var{jackstat} =} jackknife (@var{E}, @var{x})
## @deftypefnx {statistics} {@var{jackstat} =} jackknife (@var{E}, @var{x}, @dots{})
## @deftypefnx {statistics} {@var{jackstat} =} jackknife (@dots{}, @qcode{'Options'}, @var{options})
##
## Compute jackknife estimates of a parameter taking one or more given samples
## as parameters.
##
## In particular, @var{E} is the estimator to be jackknifed as a function name,
## handle, or inline function, and @var{x} is the sample for which the estimate
## is to be taken.  Each row of @var{x} is an observation, and a row vector is
## taken as a column of them.  The @var{i}-th row of @var{jackstat} holds the
## estimator's value on @var{x} with its @var{i}-th row omitted, reshaped to a
## row, so an estimator returning @math{m} values gives an @math{n*m} matrix
## for @math{n} observations, and a scalar estimator a column vector.  Every
## sample must give the estimator's result the same number of elements as the
## whole of @var{x} does.  With fewer than two observations there is nothing
## to omit, and @var{jackstat} is the estimator's value on @var{x} itself.
##
## @example
## @group
## jackstat (@var{i},:) = @var{E} (@var{x}([1:@var{i}-1, @var{i}+1:end],:))
## @end group
## @end example
##
## Depending on the number of samples to be used, the estimator must have the
## appropriate form:
## @itemize
## @item
## If only one sample is used, then the estimator need not be concerned with
## cell arrays, for example jackknifing the standard deviation of a sample can
## be performed with @code{@var{jackstat} = jackknife (@@std, rand (100, 1))}.
## @item
## If, however, more than one sample is to be used, the samples must all have
## the same number of rows, and the estimator must address them as elements of
## a cell-array, in which they are aggregated in their order of appearance:
## @end itemize
##
## @example
## @group
## @var{jackstat} = jackknife (@@(x) std(x@{1@})/var(x@{2@}),
## rand (100, 1), randn (100, 1))
## @end group
## @end example
##
## @qcode{'Options'} takes an options structure, as made by @code{statset},
## for MATLAB compatibility; it changes nothing, the samples being evaluated
## one after another.
##
## If all goes well, a theoretical value @var{P} for the parameter is already
## known, @var{n} is the sample size,
##
## @code{@var{t} = @var{n} * @var{E}(@var{x}) - (@var{n} - 1) *
## mean(@var{jackstat})}
##
## and
##
## @code{@var{v} = sumsq(@var{n} * @var{E}(@var{x}) - (@var{n} - 1) *
## @var{jackstat} - @var{t}) / (@var{n} * (@var{n} - 1))}
##
## then
##
## @code{(@var{t}-@var{P})/sqrt(@var{v})} should follow a t-distribution with
## @var{n}-1 degrees of freedom.
##
## Jackknifing is a well known method to reduce bias.
## Further details can be found in:
## @subheading References
##
## @enumerate
## @item
## Rupert G. Miller. The jackknife - a review. Biometrika (1974), 61(1):1-15.
## doi:10.1093/biomet/61.1.1
## @item
## Rupert G. Miller. Jackknifing Variances. Ann. Math. Statist. (1968),
## Volume 39, Number 2, 567-582. doi:10.1214/aoms/1177698418
## @end enumerate
## @end deftypefn

function jackstat = jackknife (anEstimator, varargin)

  if (nargin < 2)
    print_usage ();
  endif

  ## Convert function name to handle if necessary, or throw an error.
  if (ischar (anEstimator))
    anEstimator = str2func (anEstimator);
  elseif (! is_function_handle (anEstimator))
    error (strcat ("jackknife: estimators must be passed as function", ...
                   " names or handles."));
  endif

  ## An options structure is taken and changes nothing: the samples are
  ## evaluated one after another and nothing is drawn at random.
  if (numel (varargin) >= 3 && ischar (varargin{end-1})
      && strcmpi (varargin{end-1}, 'Options'))
    if (! (isstruct (varargin{end}) || isempty (varargin{end})))
      error ("jackknife: 'Options' must be a structure.");
    endif
    varargin(end-1:end) = [];
  endif

  ## Each row is an observation, and a row vector is a column of them.
  for i = 1:numel (varargin)
    if (isrow (varargin{i}))
      varargin{i} = varargin{i}(:);
    endif
  endfor
  n = cellfun (@rows, varargin);
  if (any (n != n(1)))
    error ("jackknife: all passed data must have the same number of rows.");
  endif
  n = n(1);

  ## A single sample reaches the estimator as it is, several as a cell array.
  if (numel (varargin) == 1 && (isnumeric (varargin{1})
                                || islogical (varargin{1})))
    estimate = @(idx) anEstimator (varargin{1}(idx,:));
  else
    estimate = @(idx) anEstimator (cellfun (@(x) x(idx,:), varargin, ...
                                            'UniformOutput', false));
  endif

  ## One row per jackknife sample, holding the estimate as a row.  With
  ## fewer than two observations there is no sample to leave one out of,
  ## and the estimate on the whole data is returned.
  jackstat = estimate (1:n);
  jackstat = jackstat(:).';
  if (n < 2)
    return;
  endif
  m = numel (jackstat);
  jackstat = repmat (jackstat, n, 1);
  for k = 1:n
    s = estimate ([1:k-1, k+1:n]);
    if (numel (s) != m)
      error (strcat ("jackknife: the estimator returned %d values for a", ...
                     " jackknife sample and %d for the whole data."), ...
             numel (s), m);
    endif
    jackstat(k,:) = s(:).';
  endfor

endfunction


%!demo
%! rng (42);
%! for k = 1:1000
%!   x = rand (10, 1);
%!   s(k) = std (x);
%!   jackstat = jackknife (@std, x);
%!   j(k) = 10 * std (x) - 9 * mean (jackstat);
%! endfor
%! figure ();
%! hist ([s', j'], 0:sqrt (1/12)/10:2*sqrt (1/12))

%!demo
%! rng (42);
%! for k = 1:1000
%!   x = randn (1, 50);
%!   y = rand (1, 50);
%!   jackstat = jackknife (@(x) std (x{1})/std (x{2}), y, x);
%!   j(k) = 50 * std (y) / std (x) - 49 * mean (jackstat);
%!   v(k) = sumsq ((50 * std (y) / std (x) - 49 * jackstat) - j(k)) / (50 * 49);
%! endfor
%! t = (j - sqrt (1 / 12)) ./ sqrt (v);
%! figure ();
%! plot (sort (tcdf (t, 49)), ...
%!       '-;Almost linear mapping indicates good fit with t-distribution.;')

## Test output
%!test
%! ##Example from Quenouille, Table 1
%! d=[0.18 4.00 1.04 0.85 2.14 1.01 3.01 2.33 1.57 2.19];
%! jackstat = jackknife ( @(x) 1/mean (x), d );
%! assert_equal ( 10 / mean (d) - 9 * mean (jackstat), 0.5240, 1e-5 );

%!test
%! ## Empty input
%! assert_equal (jackknife (@mean, []), NaN);

%!test
%! ## Single-element input
%! assert_equal (jackknife (@mean, 5), 5);

%!test
%! ## Estimator returning multiple values
%! expected = [2.5, sqrt(0.5);
%!             2.0, sqrt(2);
%!             1.5, sqrt(0.5)];
%! jackstat = jackknife (@(x) [mean(x); std(x)], [1 2 3]);
%! assert_equal (jackstat, expected, 1e-5);

## Expected values below are MATLAB R2024a's.
%!test
%! ## Each row of a matrix is an observation
%! assert_equal (jackknife (@mean, magic (3)), ...
%!               [3.5, 7, 4.5; 6, 5, 4; 5.5, 3, 6.5]);

%!test
%! assert_equal (jackknife (@mean, [1, 2; 3, 4; 5, 7]), ...
%!               [4, 5.5; 3, 4.5; 2, 3]);

%!test
%! ## A row vector is a column of observations
%! assert_equal (jackknife (@mean, [1, 2, 3]), [2.5; 2; 1.5]);

%!test
%! ## An estimate is reshaped to a row
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 10];
%! assert_equal (jackknife (@(v) [v(1), v(end)], x), [4, 10; 1, 10; 1, 6]);

%!test
%! assert_equal (jackknife (@(v) sum (v(:)) * ones (2, 2), [1, 2, 3, 4]), ...
%!               [9, 9, 9, 9; 8, 8, 8, 8; 7, 7, 7, 7; 6, 6, 6, 6]);

%!test
%! assert_equal (jackknife (@mean, zeros (1, 0)), NaN);

%!test
%! assert_equal (jackknife (@(x) [mean(x); std(x)], []), [NaN, NaN]);

%!test
%! assert_equal (jackknife (@(x) [mean(x); std(x)], 5), [5, 0]);

%!test
%! assert_equal (jackknife (@(x) mean (x) > 2, [1, 2, 3]), ...
%!               [true; false; false]);

%!test
%! assert_equal (jackknife (@mean, single ([1, 2, 3])), single ([2.5; 2; 1.5]));

%!test
%! assert_equal (jackknife ('mean', [1, 2, 3]), [2.5; 2; 1.5]);

%!test
%! assert_equal (jackknife (@mean, [1, 2, 3], 'Options', statset ()), ...
%!               [2.5; 2; 1.5]);

%!test
%! ## Several samples, each omitting the same row
%! assert_equal (jackknife (@(x) mean (x{1}) - mean (x{2}), [1, 2, 3], ...
%!                          [4, 5, 7]), [-3.5; -3.5; -3]);

%!test
%! assert_equal (jackknife (@(x) mean (x{1}) + mean (x{2}), 5, 6), 11);

## Test input validation
%!error<Invalid call to jackknife> jackknife (@mean)
%!error<jackknife: estimators must be passed as function names or handles.> ...
%! jackknife (1, [1, 2, 3])
%!error<jackknife: 'Options' must be a structure.> ...
%! jackknife (@mean, [1, 2, 3], 'Options', 1)
%!error<jackknife: all passed data must have the same number of rows.> ...
%! jackknife (@(x) mean (x{1}), [1, 2, 3], [4, 5])
%!error<jackknife: the estimator returned 2 values for a jackknife sample and 3 for the whole data.> ...
%! jackknife (@(v) v, [1, 2, 3])

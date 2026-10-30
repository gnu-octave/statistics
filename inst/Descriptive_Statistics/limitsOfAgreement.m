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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{LoA} =} limitsOfAgreement (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{LoA} =} limitsOfAgreement (@var{X}, @var{Y}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{LoA}, @var{stats}] =} limitsOfAgreement (@dots{})
##
## Bland-Altman limits of agreement between two measurement methods.
##
## @code{@var{LoA} = limitsOfAgreement (@var{X}, @var{Y})} compares two
## methods that measured the same subjects, the measurements of the first in
## the vector @var{X} and of the second, in the same order, in @var{Y}.  With
## @code{d = @var{X} - @var{Y}} the differences, @code{dbar} their mean and
## @code{s} their standard deviation, the bias of the first method against
## the second is @code{dbar} and the limits of agreement are
## @code{dbar - 1.96 * s} and @code{dbar + 1.96 * s}, between which about 95%
## of the differences between the methods fall when the differences are
## normal.  Two methods
## agree well enough to be used in place of one another when the limits are
## too close to matter in practice, a judgement the data cannot make.
##
## @var{LoA} is a table with the rows @qcode{Bias}, @qcode{LowerLimit} and
## @qcode{UpperLimit} and the variables @qcode{Estimate} and
## @qcode{ConfidenceIntervals}, each row's confidence interval at the level
## @math{100 (1 - alpha)} percent: @code{dbar + [-1, 1] * t * s / sqrt (n)}
## for the bias and @code{[-1, 1] * t * sqrt (3 * s^2 / n)} about each limit,
## @code{t} being the @math{1 - alpha / 2} quantile of Student's t
## distribution on @math{n - 1} degrees of freedom and @math{n} the number of
## subjects (Bland and Altman, 1986).  @var{X} and @var{Y} are vectors of
## type double or single, of the same length, and a subject with a missing
## value (@qcode{NaN}) in either is left out.  Where either is single, so are
## the results.
##
## The following @qcode{Name-Value} pairs are supported:
##
## @multitable @columnfractions 0.2 0.8
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'Multiplier'} @tab The number of standard deviations from
## the bias to each limit, a positive scalar.  The default is 1.96, which
## Bland and Altman use.
##
## @item @qcode{'Alpha'} @tab The significance level of the confidence
## intervals, a scalar in the range @math{(0,1)}.  The default is
## @math{0.05}.
## @end multitable
##
## @code{[@var{LoA}, @var{stats}] = limitsOfAgreement (@dots{})} also returns
## a structure with the fields @qcode{sd}, the standard deviation of the
## differences, @qcode{n}, the number of subjects, @qcode{differences} and
## @qcode{means}, the difference and the mean of the two measurements of each
## subject, which @code{blandAltmanPlot} draws.
##
## The bias with its interval is the paired mean difference of
## @code{meanEffectSize}.  MATLAB has no limits of agreement, so
## @code{limitsOfAgreement} is specific to Octave.
##
## Reference: J. M. Bland and D. G. Altman.  Statistical methods for
## assessing agreement between two methods of clinical measurement.
## @emph{The Lancet}, 327(8476):307--310, 1986.
##
## @seealso{blandAltmanPlot, meanEffectSize, cohenKappa}
## @end deftypefn

function [LoA, stats] = limitsOfAgreement (X, Y, varargin)

  if (nargin < 2)
    print_usage ();
  endif
  if (! (isfloat (X) && isreal (X) && isvector (X)))
    error (strcat ("limitsOfAgreement: X must be a real vector of type", ...
                   " double or single."));
  endif
  if (! (isfloat (Y) && isreal (Y) && isvector (Y)))
    error (strcat ("limitsOfAgreement: Y must be a real vector of type", ...
                   " double or single."));
  endif
  if (numel (X) != numel (Y))
    error ("limitsOfAgreement: X and Y must have the same length.");
  endif

  ## Parse optional Name-Value paired arguments
  optNames = {'Multiplier', 'Alpha'};
  dfValues = {1.96, 0.05};
  [mult, alpha, rest] = parsePairedArguments (optNames, dfValues, varargin(:));
  if (! isempty (rest))
    if (mod (numel (rest), 2) == 0 && ischar (rest{1}))
      error ("limitsOfAgreement: unknown parameter name '%s'.", rest{1});
    endif
    error ("limitsOfAgreement: optional arguments must be Name-Value pairs.");
  endif
  if (! (isfloat (mult) && isreal (mult) && isscalar (mult) && mult > 0
         && isfinite (mult)))
    error ("limitsOfAgreement: 'Multiplier' must be a positive finite scalar.");
  endif
  if (! (isfloat (alpha) && isreal (alpha) && isscalar (alpha)
         && alpha > 0 && alpha < 1))
    error ("limitsOfAgreement: 'Alpha' must be a scalar between 0 and 1.");
  endif

  outclass = 'double';
  if (isa (X, 'single') || isa (Y, 'single'))
    outclass = 'single';
  endif
  x = double (X(:));
  y = double (Y(:));
  keep = ! (isnan (x) | isnan (y));
  x = x(keep);
  y = y(keep);
  if (any (isinf (x)) || any (isinf (y)))
    error ("limitsOfAgreement: X and Y must not hold infinite values.");
  endif
  n = numel (x);
  if (n < 2)
    error ("limitsOfAgreement: X and Y must hold at least two subjects.");
  endif

  ## Bias, limits and their intervals (Bland and Altman, 1986)
  d = x - y;
  bias = mean (d);
  s = std (d);
  t = tinv (1 - alpha / 2, n - 1);
  E = bias + [0; -1; 1] * mult * s;
  h = t * [s / sqrt(n); sqrt(3 * s ^ 2 / n); sqrt(3 * s ^ 2 / n)];
  CI = E + [-1, 1] .* h;

  LoA = table (cast (E, outclass), cast (CI, outclass), 'VariableNames', ...
               {'Estimate', 'ConfidenceIntervals'}, 'RowNames', ...
               {'Bias'; 'LowerLimit'; 'UpperLimit'});
  stats = struct ('sd', cast (s, outclass), 'n', n, 'differences', ...
                  cast (d, outclass), 'means', cast ((x + y) / 2, outclass));

endfunction

%!demo
%! ## Two meters measured the peak expiratory flow of 17 subjects
%! x = [494, 395, 516, 434, 476, 557, 413, 442, 650, 433, 417, 656, 267, ...
%!      478, 178, 423, 427];
%! y = [512, 430, 520, 428, 500, 600, 364, 380, 658, 445, 432, 626, 260, ...
%!      477, 259, 350, 451];
%! LoA = limitsOfAgreement (x, y)

## Test output
## Values checked against R 4.5.0, BlandAltmanLeh 0.3.1 (bland.altman.stats)
%!shared x, y
%! x = [494, 395, 516, 434, 476, 557, 413, 442, 650, 433, 417, 656, 267, ...
%!      478, 178, 423, 427];
%! y = [512, 430, 520, 428, 500, 600, 364, 380, 658, 445, 432, 626, 260, ...
%!      477, 259, 350, 451];
%!test
%! [L, s] = limitsOfAgreement (x, y);
%! assert_equal (L.Properties.RowNames, {'Bias'; 'LowerLimit'; 'UpperLimit'});
%! assert_equal (L.Estimate, [-2.11764705882353; -78.097301611094; ...
%!                            73.8620074934469], 1e-12);
%! assert_equal (L.ConfidenceIntervals, ...
%!               [-22.0488376966452, 17.8135435789981; ...
%!                -112.619136451142, -43.5754667710458; ...
%!                39.3401726533987, 108.383842333495], 1e-12);
%! assert_equal (s.sd, 38.7651298736074, 1e-12);
%! assert_equal (s.n, 17);
%! assert_equal (s.differences, (x - y)');
%! assert_equal (s.means, ((x + y) / 2)');
%!test
%! L = limitsOfAgreement ([x, NaN], [y, 3], 'Multiplier', 2, 'Alpha', 0.1);
%! assert_equal (L.Estimate, [-2.11764705882353; -79.6479068060383; ...
%!                            75.4126126883912], 1e-12);
%! assert_equal (L.ConfidenceIntervals, ...
%!               [-18.5323144480251, 14.2970203303781; ...
%!                -108.078944713479, -51.2168688985971; ...
%!                46.9815747809501, 103.843650595832], 1e-12);
%!test
%! L = limitsOfAgreement (x, y);
%! T = meanEffectSize (x, y, 'Paired', true);
%! assert_equal (L.Estimate(1), T.Effect, 1e-12);
%! assert_equal (L.ConfidenceIntervals(1,:), T.ConfidenceIntervals, 1e-12);
%!test
%! [L, s] = limitsOfAgreement (single (x), y);
%! assert_equal (class (L.Estimate), 'single');
%! assert_equal (class (L.ConfidenceIntervals), 'single');
%! assert_equal (class (s.sd), 'single');
%! assert_equal (double (L.Estimate), [-2.11764705882353; -78.097301611094; ...
%!                                     73.8620074934469], 1e-4);

## Test input validation
%!error <Invalid call to limitsOfAgreement> limitsOfAgreement (1)
%!error <limitsOfAgreement: X must be a real vector of type double or single.> ...
%! limitsOfAgreement (int8 ([1, 2]), [1, 2])
%!error <limitsOfAgreement: Y must be a real vector of type double or single.> ...
%! limitsOfAgreement ([1, 2], [1, 2] + 1i)
%!error <limitsOfAgreement: X and Y must have the same length.> ...
%! limitsOfAgreement ([1, 2, 3], [1, 2])
%!error <limitsOfAgreement: unknown parameter name 'Foo'.> ...
%! limitsOfAgreement (x, y, 'Foo', 1)
%!error <limitsOfAgreement: optional arguments must be Name-Value pairs.> ...
%! limitsOfAgreement (x, y, 'Alpha')
%!error <limitsOfAgreement: 'Multiplier' must be a positive finite scalar.> ...
%! limitsOfAgreement (x, y, 'Multiplier', 0)
%!error <limitsOfAgreement: 'Alpha' must be a scalar between 0 and 1.> ...
%! limitsOfAgreement (x, y, 'Alpha', 0)
%!error <limitsOfAgreement: X and Y must not hold infinite values.> ...
%! limitsOfAgreement ([1, Inf, 3], [1, 2, 3])
%!error <limitsOfAgreement: X and Y must hold at least two subjects.> ...
%! limitsOfAgreement ([1, NaN], [2, 3])

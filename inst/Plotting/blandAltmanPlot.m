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
## @deftypefn  {statistics} {} blandAltmanPlot (@var{X}, @var{Y})
## @deftypefnx {statistics} {} blandAltmanPlot (@var{ax}, @var{X}, @var{Y})
## @deftypefnx {statistics} {} blandAltmanPlot (@dots{}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {@var{H} =} blandAltmanPlot (@dots{})
##
## Bland-Altman plot of the agreement between two measurement methods.
##
## @code{blandAltmanPlot (@var{X}, @var{Y})} draws, for each subject measured
## by both methods, the difference @code{@var{X} - @var{Y}} against the mean
## of the two measurements, with a solid line at the bias, the mean
## difference, dashed lines at the limits of agreement, and dotted lines at
## the bounds of the confidence interval of each.  @code{limitsOfAgreement}
## computes them, and the options it takes, @qcode{'Multiplier'} and
## @qcode{'Alpha'}, are passed on to it, as are @var{X} and @var{Y}, whose
## requirements it states.  A pattern in the points, the differences growing
## with the mean for instance, is a sign that the limits, which assume one
## spread of differences throughout, do not describe the methods well.
##
## @code{blandAltmanPlot (@var{ax}, @dots{})} draws into the axes @var{ax}
## rather than the current axes.
##
## @code{@var{H} = blandAltmanPlot (@dots{})} returns the graphics objects
## drawn, in a row vector: the points, the lines at the bias, the lower limit
## and the upper limit, and the dotted pairs of lines of the intervals of the
## bias, the lower limit and the upper limit, each pair one line object.  Each
## carries a @qcode{DisplayName}, so @code{legend} labels them.
##
## @seealso{limitsOfAgreement, gardnerAltmanPlot}
## @end deftypefn

function varargout = blandAltmanPlot (varargin)

  ## Input validation
  ax = [];
  if (numel (varargin) > 0 && isscalar (varargin{1})
      && ishghandle (varargin{1})
      && strcmp (get (varargin{1}, 'type'), 'axes'))
    ax = varargin{1};
    varargin(1) = [];
  endif
  if (numel (varargin) < 2 || ischar (varargin{2})
      || isstring (varargin{2}))
    error ("blandAltmanPlot: two measurements X and Y are required.");
  endif

  ## The bias, the limits and their intervals, the data and options checked
  ## by limitsOfAgreement
  [LoA, stats] = limitsOfAgreement (varargin{:});
  E = double (LoA.Estimate);
  CI = double (LoA.ConfidenceIntervals);
  m = double (stats.means);
  d = double (stats.differences);

  if (isempty (ax))
    ax = gca ();
  endif
  ax = newplot (ax);
  held = ishold (ax);
  hold (ax, 'on');
  co = get (ax, 'colororder');
  xr = [min(m), max(m)];
  if (xr(1) == xr(2))
    xr += [-0.5, 0.5];
  endif
  xs = [xr, NaN, xr];
  names = {'Bias', 'Lower limit', 'Upper limit'};
  h = line (ax, m, d, 'linestyle', 'none', 'marker', 'o', ...
            'color', co(1,:), 'displayname', 'Subjects');
  styles = {'-', '--', '--'};
  for i = 1:3
    h(end+1) = line (ax, xr, [E(i), E(i)], 'color', [0, 0, 0], ...
                     'linestyle', styles{i}, 'displayname', names{i});
  endfor
  for i = 1:3
    h(end+1) = line (ax, xs, [CI(i,[1, 1]), NaN, CI(i,[2, 2])], ...
                     'color', [0, 0, 0], 'linestyle', ':', ...
                     'displayname', [names{i}, ' interval']);
  endfor
  if (! held)
    hold (ax, 'off');
  endif
  xlabel (ax, 'Mean of X and Y');
  ylabel (ax, 'X - Y');

  ## The handles are handed back only where they were asked for, so a call
  ## made for the drawing alone prints nothing
  if (nargout > 0)
    varargout{1} = h;
  endif

endfunction

%!demo
%! ## Two meters measured the peak expiratory flow of 17 subjects
%! x = [494, 395, 516, 434, 476, 557, 413, 442, 650, 433, 417, 656, 267, ...
%!      478, 178, 423, 427];
%! y = [512, 430, 520, 428, 500, 600, 364, 380, 658, 445, 432, 626, 260, ...
%!      477, 259, 350, 451];
%! H = blandAltmanPlot (x, y);
%! legend (H(1:4), 'location', 'northeastoutside');

## Test output
%!shared x, y
%! x = [494, 395, 516, 434, 476, 557, 413, 442, 650, 433, 417, 656, 267, ...
%!      478, 178, 423, 427];
%! y = [512, 430, 520, 428, 500, 600, 364, 380, 658, 445, 432, 626, 260, ...
%!      477, 259, 350, 451];
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   H = blandAltmanPlot (x, y);
%!   assert_equal (size (H), [1, 7]);
%!   assert_equal (get (H(1), 'xdata'), (x + y) / 2);
%!   assert_equal (get (H(1), 'ydata'), x - y);
%!   L = limitsOfAgreement (x, y);
%!   assert_equal (get (H(2), 'ydata'), L.Estimate([1, 1])');
%!   assert_equal (get (H(3), 'ydata'), L.Estimate([2, 2])');
%!   assert_equal (get (H(4), 'ydata'), L.Estimate([3, 3])');
%!   assert_equal (get (H(5), 'ydata'), [L.ConfidenceIntervals(1,[1, 1]), ...
%!                 NaN, L.ConfidenceIntervals(1,[2, 2])]);
%!   assert_equal (get (H(2:4), 'linestyle'), {'-'; '--'; '--'});
%!   assert_equal (get (H(1), 'displayname'), 'Subjects');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   ax = axes ();
%!   plot (ax, 1:3, 1:3);
%!   hold (ax, 'on');
%!   H = blandAltmanPlot (ax, x, y, 'Multiplier', 2);
%!   assert_equal (ishold (ax), true);
%!   assert_equal (numel (get (ax, 'children')), 8);
%!   L = limitsOfAgreement (x, y, 'Multiplier', 2);
%!   assert_equal (get (H(3), 'ydata'), L.Estimate([2, 2])');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   blandAltmanPlot (x, y);
%!   assert_equal (ishold (gca ()), false);
%!   assert_equal (get (get (gca (), 'xlabel'), 'string'), 'Mean of X and Y');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Test input validation
%!error <blandAltmanPlot: two measurements X and Y are required.> ...
%! blandAltmanPlot ([1, 2, 3])
%!error <blandAltmanPlot: two measurements X and Y are required.> ...
%! blandAltmanPlot ([1, 2, 3], 'Alpha', 0.1)
%!error <limitsOfAgreement: X and Y must have the same length.> ...
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   blandAltmanPlot ([1, 2, 3], [1, 2]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

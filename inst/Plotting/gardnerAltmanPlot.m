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
## @deftypefn  {statistics} {} gardnerAltmanPlot (@var{X}, @var{Y})
## @deftypefnx {statistics} {} gardnerAltmanPlot (@var{ax}, @var{X}, @var{Y})
## @deftypefnx {statistics} {} gardnerAltmanPlot (@dots{}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {@var{H} =} gardnerAltmanPlot (@dots{})
##
## Gardner-Altman plot of two samples and the effect size between them.
##
## @code{gardnerAltmanPlot (@var{X}, @var{Y})} draws the observations of the
## samples @var{X} and @var{Y} at the horizontal positions 1 and 2, and the
## difference between their means, with its confidence interval, at position
## 3, measured against a second vertical axis on the right.  That axis is
## placed so that its zero is level with the centre of @var{Y} and the effect
## size level with the centre of @var{X}, so the effect reads off the data it
## comes from.  @var{X} and @var{Y} are vectors of type double or single, and
## a missing value (@qcode{NaN}) is left out.
##
## Two unpaired samples are drawn as swarms, the points of each spread
## sideways where they crowd, with a black line from the centre of each
## sample to the effect size: solid for @var{X}, dashed for @var{Y}.  The
## centre is the median for @qcode{'mediandiff'} and the mean otherwise.
## Paired samples are drawn as a line from each observation of @var{X} to its
## partner in @var{Y}, the pairs that increase, decrease and stay equal each
## in a colour of their own.  For @qcode{'cliff'} and @qcode{'kstest'}, whose
## values do not share the scale of the data, the right axis spans the values
## the effect can take, @math{[-1.1, 1.1]} and @math{[-0.1, 1.1]}, and two
## unpaired samples get a reference line at zero instead of the lines from
## their centres.
##
## @code{gardnerAltmanPlot (@var{ax}, @dots{})} draws into the axes @var{ax}
## rather than the current axes.
##
## @code{gardnerAltmanPlot (@dots{}, @var{Name}, @var{Value})} takes the
## options @qcode{'Effect'}, @qcode{'Paired'}, @qcode{'VarianceType'},
## @qcode{'Alpha'}, @qcode{'ConfidenceIntervalType'}, @qcode{'NumBootstraps'},
## @qcode{'BootstrapOptions'} and @qcode{'Resampling'} of
## @code{meanEffectSize}, which computes the effect size and its interval;
## @qcode{'Effect'} names a single effect.  @qcode{'Mean'} does not apply,
## there always being two samples.
##
## @code{@var{H} = gardnerAltmanPlot (@dots{})} returns the graphics objects
## drawn, in a row vector.  For two unpaired samples they are the two swarms,
## the error bar, and the lines from the centres of @var{X} and @var{Y}, or
## the reference line.  For paired samples they are the lines of the pairs
## that increase, decrease and stay equal, each where there is any, and the
## error bar.  Where @qcode{'ConfidenceIntervalType'} is @qcode{'none'} the
## error bar is a line holding a single marker.
##
## MATLAB draws the effect size on a second vertical axis of the same axes;
## Octave has no such axis, so here it is a second axes laid over the first,
## which follows the first's position and vertical limits.  The error bar is
## Octave's @code{errorbar} object, whose bounds are its @qcode{'ldata'} and
## @qcode{'udata'} properties.
##
## Reference: M. J. Gardner and D. G. Altman (1986).  Confidence intervals
## rather than P values: estimation rather than hypothesis testing.  British
## Medical Journal, 292(6522), 746-750.
##
## @seealso{meanEffectSize, swarmchart, errorbar}
## @end deftypefn

function varargout = gardnerAltmanPlot (varargin)

  ## Input validation
  ax = [];
  if (numel (varargin) > 0 && isscalar (varargin{1})
      && ishghandle (varargin{1})
      && strcmp (get (varargin{1}, 'type'), 'axes'))
    ax = varargin{1};
    varargin(1) = [];
  endif
  if (numel (varargin) < 2 || ischar (varargin{2})
      || isa (varargin{2}, 'string'))
    error ("gardnerAltmanPlot: two samples X and Y are required.");
  endif
  X = varargin{1};
  Y = varargin{2};
  args = varargin(3:end);
  optNames = {'Effect', 'Paired', 'VarianceType', 'Alpha', ...
              'ConfidenceIntervalType', 'NumBootstraps', ...
              'BootstrapOptions', 'Resampling'};
  dfValues = {'meandiff', false, [], [], [], [], [], []};
  [effect, paired, ~, ~, ~, ~, ~, ~, rest] = ...
                       parsePairedArguments (optNames, dfValues, args(:));
  if (! isempty (rest))
    error ("gardnerAltmanPlot: invalid optional paired argument.");
  endif
  if ((iscell (effect) || isa (effect, 'string')) && numel (effect) != 1)
    error ("gardnerAltmanPlot: 'Effect' must name a single effect.");
  endif

  ## The effect size and its interval, the options checked by meanEffectSize
  ## and its messages raised under this name
  try
    T = meanEffectSize (X, Y, args{:});
  catch err
    error ("gardnerAltmanPlot: %s", ...
           regexprep (err.message, '^meanEffectSize: ', ''));
  end_try_catch
  name = T.Properties.RowNames{1};
  e = double (T.Effect);
  noCI = ! any (strcmp (T.Properties.VariableNames, 'ConfidenceIntervals'));
  if (! noCI)
    ci = double (T.ConfidenceIntervals);
  endif
  ## A switch meanEffectSize accepted: true, false, 0, 1, 'on' or 'off'
  if (ischar (paired) || isa (paired, 'string'))
    paired = strcmpi (char (paired), 'on');
  else
    paired = logical (paired);
  endif

  ## Missing values left out, of both samples where they are paired
  x = double (X(:));
  y = double (Y(:));
  if (paired)
    keep = ! (isnan (x) | isnan (y));
    x = x(keep);
    y = y(keep);
  else
    x = x(! isnan (x));
    y = y(! isnan (y));
  endif
  if (strcmp (name, 'MedianDifference'))
    cx = median (x);
    cy = median (y);
    centre = 'Median';
  else
    cx = mean (x);
    cy = mean (y);
    centre = 'Mean';
  endif
  bounded = any (strcmp (name, {'CliffsDelta', ...
                                'KolmogorovSmirnovStatistic'}));

  ## The data on the left axes
  if (isempty (ax))
    ax = gca ();
  endif
  ## The right axes and listeners of an earlier plot into these axes go
  removeOverlay (ax);
  ax = newplot (ax);
  held = ishold (ax);
  hold (ax, 'on');
  co = get (ax, 'colororder');
  nc = 0;
  if (paired)
    kinds = {y > x, y < x, y == x};
    h = [];
    for k = 1:3
      take = kinds{k};
      if (any (take))
        nc += 1;
        xd = repmat ([1; 2; NaN], 1, sum (take));
        yd = [x(take)'; y(take)'; NaN(1, sum (take))];
        h(end+1) = line (ax, xd(:)', yd(:)', 'color', ...
                         co(mod (nc - 1, rows (co)) + 1,:));
      endif
    endfor
  else
    h = zeros (1, 2);
    h(1) = swarmchart (ax, ones (numel (x), 1), x, 36, co(1,:), ...
                       'XJitterWidth', 0.5, 'DisplayName', 'X');
    h(2) = swarmchart (ax, 2 * ones (numel (y), 1), y, 36, co(2,:), ...
                       'XJitterWidth', 0.5, 'DisplayName', 'Y');
    nc = 2;
  endif
  ebColor = co(mod (nc, rows (co)) + 1,:);
  set (ax, 'xlim', [0.4, 3.3], 'xtick', 1:3, ...
       'xticklabel', {'X', 'Y', effectLabel(name)}, 'box', 'on');
  title (ax, 'Gardner-Altman Plot');

  ## The effect on a second axes laid over the first, with the axis on the
  ## right; MATLAB puts it on a second axis of the same axes
  ax2 = axes ('parent', get (ax, 'parent'), 'position', ...
              get (ax, 'position'), 'color', 'none', 'yaxislocation', ...
              'right', 'xlim', [0.4, 3.3], 'xtick', [], 'box', 'off', ...
              'ycolor', ebColor);
  hold (ax2, 'on');
  if (noCI)
    eb = line (ax2, 3, e, 'marker', 'square', 'linestyle', 'none', ...
               'color', ebColor, 'displayname', effectLabel (name));
  else
    eb = errorbar (ax2, 3, e, e - ci(1), ci(2) - e);
    set (eb, 'marker', 'square', 'color', ebColor, 'displayname', ...
         effectLabel (name));
  endif
  hold (ax2, 'off');
  h(end+1) = eb;

  ## The lines from each centre to the effect, or the reference line
  if (! paired)
    if (bounded)
      h(end+1) = line (ax2, [2.5, 3.5], [0, 0], 'color', [0, 0, 0]);
    else
      h(end+1) = line (ax, [1, 3], [cx, cx], 'color', [0, 0, 0], ...
                       'marker', '*', 'displayname', ['X ', centre]);
      h(end+1) = line (ax, [2, 3], [cy, cy], 'color', [0, 0, 0], ...
                       'marker', '*', 'linestyle', '--', ...
                       'displayname', ['Y ', centre]);
    endif
  endif
  if (! held)
    hold (ax, 'off');
  endif

  ## The right axis follows the left one: its limits, unless they are
  ## fixed, and its position
  if (bounded)
    k = [];
    if (strcmp (name, 'CliffsDelta'))
      set (ax2, 'ylim', [-1.1, 1.1]);
    else
      set (ax2, 'ylim', [-0.1, 1.1]);
    endif
  else
    k = e / (cx - cy);
    if (! isfinite (k) || k == 0)
      k = 1;
    endif
  endif
  onYlim = @(~, ~) followLimits (ax, ax2, cy, k);
  onPosition = @(~, ~) followPosition (ax, ax2);
  onYlim ();
  addlistener (ax, 'ylim', onYlim);
  addlistener (ax, 'position', onPosition);
  setappdata (ax, 'gardnerAltmanPlot', struct ('axes', ax2, 'ylim', ...
                                               onYlim, 'position', onPosition));
  ## The right axes lives as long as the data on the left: clearing or
  ## deleting the left axes deletes that data, and with it the right axes
  set (h(1), 'deletefcn', @(~, ~) removeOverlay (ax));
  set (ancestor (ax, 'figure'), 'currentaxes', ax);

  ## The handles are handed back only where they were asked for, so a call
  ## made for the drawing alone prints nothing
  if (nargout > 0)
    varargout{1} = h;
  endif

endfunction

## Delete the right axes laid over the axes AX, and the listeners that keep
## it in step
function removeOverlay (ax)
  if (! ishghandle (ax) || ! isappdata (ax, 'gardnerAltmanPlot'))
    return;
  endif
  old = getappdata (ax, 'gardnerAltmanPlot');
  rmappdata (ax, 'gardnerAltmanPlot');
  dellistener (ax, 'ylim', old.ylim);
  dellistener (ax, 'position', old.position);
  if (ishghandle (old.axes))
    delete (old.axes);
  endif
endfunction

## The limits of the right axes AX2 from those of the left axes AX, so that
## zero is level with CY and the scale is K times the data's; none if K is
## empty
function followLimits (ax, ax2, cy, k)
  if (! isempty (k) && ishghandle (ax) && ishghandle (ax2))
    set (ax2, 'ylim', (get (ax, 'ylim') - cy) * k);
  endif
endfunction

## The position of the right axes AX2 kept on that of the left axes AX
function followPosition (ax, ax2)
  if (ishghandle (ax) && ishghandle (ax2))
    set (ax2, 'position', get (ax, 'position'));
  endif
endfunction

## The label of each effect on the horizontal axis
function out = effectLabel (name)
  switch (name)
    case 'MeanDifference'
      out = 'Mean Difference';
    case 'CohensD'
      out = 'Cohen''s d';
    case 'GlasssDelta'
      out = 'Glass''s delta';
    case 'CliffsDelta'
      out = 'Cliff''s Delta';
    case 'MedianDifference'
      out = 'Median Difference';
    case 'RobustCohensD'
      out = 'Robust Cohen''s d';
    case 'AKPCohensD'
      out = 'AKP Cohen''s d';
    case 'KolmogorovSmirnovStatistic'
      out = 'Kolmogorov-Smirnov Statistic';
  endswitch
endfunction

%!demo
%! ## Petal length of two iris species.  The swarms show every flower, and
%! ## the error bar on the right the difference between the mean lengths,
%! ## read on the right axis, whose zero is level with the mean of the
%! ## second species.
%! load fisheriris
%! x = meas(51:100,3);
%! y = meas(101:150,3);
%! rng (42);
%! gardnerAltmanPlot (x, y);

%!demo
%! ## The same comparison as Cohen's d, the difference in units of the
%! ## pooled standard deviation.  The right axis is scaled so that the
%! ## effect still sits level with the mean of the first species.
%! load fisheriris
%! x = meas(51:100,3);
%! y = meas(101:150,3);
%! rng (42);
%! gardnerAltmanPlot (x, y, 'Effect', 'cohen');

%!demo
%! ## Paired samples: two exam grades of the same students.  Each line joins
%! ## one student's two grades, coloured by whether the second is higher,
%! ## lower or equal.
%! load examgrades
%! gardnerAltmanPlot (grades(:,1), grades(:,2), 'Paired', true);

%!shared x, y, yp, isScatter
%! x = [2.1; 3.4; 1.9; 5.6; 4.4; 3.8; 2.7; 6.1; 3.3; 4.9];
%! y = [1.2; 2.8; 0.9; 2.2; 3.1; 1.7; 2.5; 0.4; 1.9; 3.6; 2.0; 1.1];
%! yp = [1.8; 2.9; 2.2; 4.1; 4.9; 2.6; 3.0; 4.8; 2.1; 4.0];
%! ## Octave's scatter builds an hggroup under the gnuplot toolkit
%! isScatter = ! strcmp (graphics_toolkit (), 'gnuplot');
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   H = gardnerAltmanPlot (x, y);
%!   assert_equal (size (H), [1, 5]);
%!   if (isScatter)
%!     assert_equal (get (H(1:2), 'type'), {'scatter'; 'scatter'});
%!   endif
%!   assert_equal (get (H(3:5), 'type'), {'hggroup'; 'line'; 'line'});
%!   assert_equal (get (H(1), 'ydata'), x);
%!   assert_equal (get (H(2), 'ydata'), y);
%!   assert_equal (get (H(3), 'xdata'), 3);
%!   assert_equal (get (H(3), 'ydata'), 1.87, -1e-14);
%!   assert_equal (get (H(3), 'ldata'), 1.0606770309586725, -1e-13);
%!   assert_equal (get (H(3), 'udata'), 1.0606770309586709, -1e-13);
%!   assert_equal (get (H(4), 'xdata'), [1, 3]);
%!   assert_equal (get (H(4), 'ydata'), [3.82, 3.82], -1e-14);
%!   assert_equal (get (H(5), 'xdata'), [2, 3]);
%!   assert_equal (get (H(5), 'ydata'), [1.95, 1.95], -1e-14);
%!   assert_equal (get (H(5), 'linestyle'), '--');
%!   ax = gca ();
%!   assert_equal (get (ax, 'xticklabel'), {'X'; 'Y'; 'Mean Difference'});
%!   assert_equal (get (get (ax, 'title'), 'string'), 'Gardner-Altman Plot');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   H = gardnerAltmanPlot (x, y);
%!   ax = gca ();
%!   ax2 = get (H(3), 'parent');
%!   assert_equal (get (ax2, 'ylim'), get (ax, 'ylim') - 1.95, -1e-14);
%!   set (ax, 'ylim', [-2, 9]);
%!   assert_equal (get (ax2, 'ylim'), [-3.95, 7.05], -1e-14);
%!   assert_equal (get (ax2, 'position'), get (ax, 'position'));
%!   assert_equal (get (ax2, 'yaxislocation'), 'right');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   H = gardnerAltmanPlot (x, y, 'Effect', 'cohen');
%!   ax = gca ();
%!   ax2 = get (H(3), 'parent');
%!   assert_equal (get (H(3), 'ydata'), 1.51473229411906, -1e-13);
%!   assert_equal (get (ax2, 'ylim'), (get (ax, 'ylim') - 1.95) ...
%!                 * 1.51473229411906 / 1.87, -1e-13);
%!   assert_equal (get (ax, 'xticklabel'), {'X'; 'Y'; 'Cohen''s d'});
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   H = gardnerAltmanPlot (x, y, 'Effect', 'cliff');
%!   assert_equal (size (H), [1, 4]);
%!   assert_equal (get (H(3), 'ydata'), 0.725, -1e-14);
%!   assert_equal (get (H(4), 'xdata'), [2.5, 3.5]);
%!   assert_equal (get (H(4), 'ydata'), [0, 0]);
%!   assert_equal (get (get (H(3), 'parent'), 'ylim'), [-1.1, 1.1]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   rand ('state', 1);
%!   H = gardnerAltmanPlot (x, y, 'Effect', 'kstest', 'NumBootstraps', 100);
%!   assert_equal (size (H), [1, 4]);
%!   assert_equal (get (get (H(3), 'parent'), 'ylim'), [-0.1, 1.1]);
%!   assert_equal (get (gca (), 'xticklabel'), ...
%!                 {'X'; 'Y'; 'Kolmogorov-Smirnov Statistic'});
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   rand ('state', 1);
%!   H = gardnerAltmanPlot (x, y, 'Effect', 'mediandiff', 'NumBootstraps', 100);
%!   assert_equal (get (H(4), 'ydata'), [3.6, 3.6], -1e-14);
%!   assert_equal (get (H(5), 'ydata'), [1.95, 1.95], -1e-14);
%!   assert_equal (get (H(3), 'ydata'), 1.65, -1e-14);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   H = gardnerAltmanPlot (x, yp, 'Paired', true);
%!   assert_equal (get (H, 'type'), {'line'; 'line'; 'hggroup'});
%!   assert_equal (get (H(1), 'xdata'), [1, 2, NaN, 1, 2, NaN, 1, 2, NaN]);
%!   assert_equal (get (H(1), 'ydata'), ...
%!                 [1.9, 2.2, NaN, 4.4, 4.9, NaN, 2.7, 3, NaN]);
%!   assert_equal (numel (get (H(2), 'ydata')), 21);
%!   assert_equal (get (H(3), 'ydata'), 0.58, -1e-13);
%!   ax2 = get (H(3), 'parent');
%!   assert_equal (get (ax2, 'ylim'), get (gca (), 'ylim') - 3.24, -1e-14);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   H = gardnerAltmanPlot ([1; 2; 3; 4], [1; 3; 2; 4], 'Paired', true);
%!   assert_equal (get (H, 'type'), {'line'; 'line'; 'line'; 'hggroup'});
%!   assert_equal (get (H(3), 'ydata'), [1, 1, NaN, 4, 4, NaN]);
%!   co = get (gca (), 'colororder');
%!   assert_equal (get (H(4), 'color'), co(4,:));
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   H = gardnerAltmanPlot ([1; 2; 3; 4], [2; 3; 4; 5], 'Paired', true);
%!   assert_equal (get (H, 'type'), {'line'; 'hggroup'});
%!   co = get (gca (), 'colororder');
%!   assert_equal (get (H(2), 'color'), co(2,:));
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   H = gardnerAltmanPlot (x, y, 'ConfidenceIntervalType', 'none');
%!   assert_equal (get (H(3), 'type'), 'line');
%!   assert_equal (get (H(3), 'marker'), 'square');
%!   assert_equal (get (H(3), 'ydata'), 1.87, -1e-14);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   H = gardnerAltmanPlot (x, [y; NaN]);
%!   assert_equal (get (H(2), 'ydata'), y);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   ax = subplot (1, 2, 2);
%!   H = gardnerAltmanPlot (ax, x, y);
%!   assert_equal (get (H(1), 'parent'), ax);
%!   assert_equal (get (get (H(3), 'parent'), 'position'), ...
%!                 get (ax, 'position'));
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   gardnerAltmanPlot (x, y);
%!   gardnerAltmanPlot (x, yp, 'Paired', true);
%!   assert_equal (numel (findobj (hf, 'type', 'axes')), 2);
%!   plot (1:3);
%!   assert_equal (numel (findobj (hf, 'type', 'axes')), 1);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!error<gardnerAltmanPlot: two samples X and Y are required.> ...
%! gardnerAltmanPlot ([1; 2; 3])
%!error<gardnerAltmanPlot: two samples X and Y are required.> ...
%! gardnerAltmanPlot ([1; 2; 3], 'Effect', 'cohen')
%!error<gardnerAltmanPlot: invalid optional paired argument.> ...
%! gardnerAltmanPlot ([1; 2; 3], [2; 3; 4], 'Mean', 3)
%!error<gardnerAltmanPlot: 'Effect' must name a single effect.> ...
%! gardnerAltmanPlot ([1; 2; 3], [2; 3; 4], 'Effect', {'cohen', 'glass'})
%!error<gardnerAltmanPlot: X must be a nonempty vector of type double or single.> ...
%! gardnerAltmanPlot ([1, 2; 3, 4], [2; 3; 4])
%!error<gardnerAltmanPlot: effect 'glass' does not apply to paired samples.> ...
%! gardnerAltmanPlot ([1; 2; 3], [2; 3; 4], 'Paired', true, 'Effect', 'glass')

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
## @deftypefn  {statistics} {} plotPartialDependence (@var{Mdl}, @var{Vars})
## @deftypefnx {statistics} {} plotPartialDependence (@var{Mdl}, @var{Vars}, @var{Labels})
## @deftypefnx {statistics} {} plotPartialDependence (@dots{}, @var{Data})
## @deftypefnx {statistics} {} plotPartialDependence (@var{fun}, @var{Vars}, @var{Data})
## @deftypefnx {statistics} {} plotPartialDependence (@dots{}, @var{name}, @var{value})
## @deftypefnx {statistics} {@var{ax} =} plotPartialDependence (@dots{})
##
## Plot partial dependence and individual conditional expectation.
##
## @code{plotPartialDependence (@var{Mdl}, @var{Vars})} draws the partial
## dependence of the model @var{Mdl} on the predictors named by @var{Vars}:
## a line where one is named and a surface where two are.  The arguments are
## those of @code{partialDependence}, which computes what is drawn, and
## @var{Labels} is required for a classification model in the same way.
##
## @var{ax} is the axes drawn into.
##
## @multitable @columnfractions 0.28 0.02 0.7
## @headitem @var{Name} @tab @tab @var{Value}
##
## @item @qcode{'Conditional'} @tab @tab What to draw:
## @qcode{'none'} (default) draws the partial dependence alone;
## @qcode{'absolute'} draws a curve per observation, a marker where each one
## sits on its own curve, and the mean of those curves over them;
## @qcode{'centered'} draws the same with every curve shifted to start at
## zero.  Only one predictor may be varied, and a classification model may
## name only one class.
##
## @item @qcode{'Parent'} @tab @tab The axes to draw into.  The default is
## @code{gca}.
## @end multitable
##
## Every other name-value pair is passed to @code{partialDependence}; see
## there for @qcode{'QueryPoints'}, @qcode{'NumObservationsToSample'} and the
## rest.
##
## The line a conditional plot draws is the mean of the curves drawn with it,
## which is not always what @code{partialDependence} returns.  A decision
## tree, an ensemble of bagged trees and a generalized additive model are
## answered over the distribution they were fitted on, while a curve belongs
## to one observation and must come from @code{predict}, so the two part
## company for those models whenever the observations are not the ones the
## model was fitted on.  MATLAB R2024a draws it this way and so does this.
##
## @seealso{partialDependence, PredictiveModel}
## @end deftypefn

function ax = plotPartialDependence (Mdl, Vars, varargin)

  if (nargin < 2)
    error ("plotPartialDependence: too few input arguments.");
  endif
  [pd, x, y, opts] = pdCompute ('plotPartialDependence', Mdl, Vars, ...
                                varargin{:});
  F = opts.Frame;

  cond = opts.Conditional;
  if (isempty (cond))
    cond = 'none';
  endif
  if (! (ischar (cond) && isrow (cond)
         && any (strcmpi (cond, {'none', 'absolute', 'centered'}))))
    error (strcat ("plotPartialDependence: 'Conditional' must be 'none',", ...
                   " 'absolute' or 'centered'."));
  endif
  cond = lower (cond);
  two = ! isempty (y);
  if (! strcmp (cond, 'none'))
    if (two)
      error (strcat ("plotPartialDependence: a conditional plot varies", ...
                     " one predictor, not two."));
    endif
    if (F.IsClass && numel (F.LabelIdx) != 1)
      error (strcat ("plotPartialDependence: a conditional plot of a", ...
                     " classification model names one class."));
    endif
  endif

  ax = opts.Parent;
  if (isempty (ax))
    ax = gca ();
  elseif (! (isscalar (ax) && ishghandle (ax) && isaxes (ax)))
    error ("plotPartialDependence: 'Parent' must be an axes handle.");
  endif

  if (strcmp (cond, 'none'))
    pdPlotPlain (ax, F, pd, x, y, two);
  else
    pdPlotIce (ax, Mdl, F, pd, x, cond, opts.PredArgs);
  endif

  xlabel (ax, F.PredictorNames{F.Vars(1)});
  if (two)
    ylabel (ax, F.PredictorNames{F.Vars(2)});
    zlabel (ax, pdResponseLabel (Mdl, F));
  else
    ylabel (ax, pdResponseLabel (Mdl, F));
  endif
  if (strcmp (cond, 'none'))
    title (ax, 'Partial Dependence Plot');
  else
    title (ax, 'Individual Conditional Expectation Plot');
  endif

endfunction

## The partial dependence on its own: a line per class, or a surface.
function pdPlotPlain (ax, F, pd, x, y, two)

  if (two)
    surf (ax, x, y, pd);
    return;
  endif
  if (! F.IsClass)
    plot (ax, x, pd(1,:));
    return;
  endif
  k = numel (F.LabelIdx);
  co = get (ax, 'ColorOrder');
  hold (ax, 'on');
  names = cell (1, k);
  for j = 1:k
    plot (ax, x, pd(j,:), 'Color', co(mod (j - 1, rows (co)) + 1, :));
    names{j} = pdClassName (F, j);
  endfor
  hold (ax, 'off');
  if (k > 1)
    legend (ax, names{:});
  endif

endfunction

## The curves, a marker where each observation sits, and their mean over them.
function pdPlotIce (ax, Mdl, F, pd, x, cond, predArgs)

  ice = pdIce (Mdl, F, predArgs);
  ice = ice(:,:,1);
  at = pdScoreAt (Mdl, F, predArgs);
  if (strcmp (cond, 'centered'))
    off = ice(:,1);
    ice -= off;
    at -= off;
  endif

  hold (ax, 'on');
  for i = 1:rows (ice)
    plot (ax, x, ice(i,:), 'Color', [0.5, 0.5, 0.5]);
  endfor
  scatter (ax, F.X(:,F.Vars(1)), at, 'o');
  plot (ax, x, mean (ice, 1), 'Color', [1, 0, 0]);
  hold (ax, 'off');

endfunction

## What the model answers for each observation as it stands.
function s = pdScoreAt (Mdl, F, predArgs)

  if (is_function_handle (Mdl))
    s = Mdl (F.X);
    if (! isempty (F.OutCols))
      s = s(:, F.OutCols);
    endif
    s = s(:,1);
  elseif (F.IsClass)
    [~, sc] = predict (Mdl, F.X, predArgs{:});
    s = sc(:, F.LabelIdx(1));
  else
    s = predict (Mdl, F.X, predArgs{:});
    s = s(:);
  endif

endfunction

## What the vertical axis measures.
function txt = pdResponseLabel (Mdl, F)

  if (F.IsClass)
    if (numel (F.LabelIdx) == 1)
      txt = sprintf ('Score of class %s', pdClassName (F, 1));
    else
      txt = 'Scores';
    endif
    return;
  endif
  txt = 'Y';
  if (! is_function_handle (Mdl)
      && any (strcmp (properties (Mdl), 'ResponseName'))
      && ! isempty (Mdl.ResponseName))
    r = Mdl.ResponseName;
    if (ischar (r) && isrow (r))
      txt = r;
    endif
  endif

endfunction

## One class name as text, whatever type ClassNames holds it in.
function txt = pdClassName (F, j)

  lab = labelsFromIndex (F.ClassNames, F.LabelIdx(j));
  if (iscellstr (lab))
    txt = lab{1};
  elseif (ischar (lab))
    txt = lab;
  elseif (isnumeric (lab) || islogical (lab))
    txt = num2str (lab);
  else
    txt = char (lab);
    if (! isrow (txt))
      txt = txt(1,:);
    endif
  endif

endfunction

%!demo
%! ## The partial dependence of a fitted tree on one predictor.  The line is
%! ## what the tree answers on average as that predictor is varied over its
%! ## range, the others left as the data holds them.
%! load fisheriris
%! Mdl = fitrtree (meas(:,2:4), meas(:,1));
%! plotPartialDependence (Mdl, 1);

%!demo
%! ## The same tree, with one curve per observation.  Each grey curve is what
%! ## the model answers for a single flower as the predictor is varied; the
%! ## red line is their mean and the circles mark where each flower sits.
%! load fisheriris
%! Mdl = fitrtree (meas(:,2:4), meas(:,1));
%! plotPartialDependence (Mdl, 1, 'Conditional', 'absolute');

%!demo
%! ## Two predictors give a surface, and a classification model is told which
%! ## class to answer for.
%! load fisheriris
%! Mdl = fitcsvm (meas(51:end,3:4), species(51:end));
%! plotPartialDependence (Mdl, [1, 2], 'versicolor');

%!shared ppX, ppYr, ppY3, ppQ, ppLo
%! x1 = repmat ([1;2;3], 4, 1);
%! x2 = reshape (repmat (1:4, 3, 1), 12, 1);
%! ppX = [x1, x2];
%! ppYr = 10 * x2 + x1;
%! ppY3 = repmat ({'a'; 'b'; 'c'}, 4, 1);
%! ppQ = [1; 2; 3];
%! lo1 = repmat ([1;2;3], 2, 1);
%! ppLo = [lo1, ones(6, 1)];

## MATLAB parity: the plain plot draws what partialDependence returns
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   Mdl = fitrsvm (ppX, ppYr);
%!   ax = plotPartialDependence (Mdl, 1, 'QueryPoints', ppQ);
%!   ch = get (ax, 'children');
%!   assert_equal (numel (ch), 1);
%!   assert_equal (get (ch(1), 'ydata'), ...
%!                 partialDependence (Mdl, 1, 'QueryPoints', ppQ), 1e-12);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # two predictors give one surface, sized by the second then the first
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   Mdl = fitrsvm (ppX, ppYr);
%!   ax = plotPartialDependence (Mdl, [1, 2], 'QueryPoints', {ppQ, [2;4]});
%!   ch = get (ax, 'children');
%!   assert_equal (get (ch(1), 'type'), 'surface');
%!   assert_equal (size (get (ch(1), 'zdata')), [2, 3]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # a classifier draws one line per class named
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   Mdl = fitctree (ppX, ppY3);
%!   ax = plotPartialDependence (Mdl, 1, {'a', 'b'}, 'QueryPoints', ppQ);
%!   assert_equal (numel (get (ax, 'children')), 2);
%!   assert_equal (get (get (ax, 'ylabel'), 'string'), 'Scores');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # MATLAB parity: one class names the axis after it
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   Mdl = fitctree (ppX, ppY3);
%!   ax = plotPartialDependence (Mdl, 1, 'a', 'QueryPoints', ppQ, ...
%!                               'Conditional', 'absolute');
%!   assert_equal (get (get (ax, 'ylabel'), 'string'), 'Score of class a');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # a conditional plot draws a curve per observation, a marker, a line
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   Mdl = fitrsvm (ppX, ppYr);
%!   ax = plotPartialDependence (Mdl, 1, 'QueryPoints', ppQ, ...
%!                               'Conditional', 'absolute');
%!   ch = get (ax, 'children');
%!   assert_equal (numel (ch), rows (ppX) + 2);
%!   assert_equal (sum (strcmp (get (ch, 'type'), 'scatter')), 1);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # MATLAB parity: centering shifts every curve to start at zero
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   Mdl = fitrsvm (ppX, ppYr);
%!   ax = plotPartialDependence (Mdl, 1, 'QueryPoints', ppQ, ...
%!                               'Conditional', 'centered');
%!   ch = get (ax, 'children');
%!   for k = 1:numel (ch)
%!     if (strcmp (get (ch(k), 'type'), 'line'))
%!       yd = get (ch(k), 'ydata');
%!       assert_equal (yd(1), 0, 1e-12);
%!     endif
%!   endfor
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## MATLAB parity: the line of a conditional plot is the mean of the curves it
## draws, which for a tree over other observations is not what
## partialDependence returns
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   Mdl = fitrtree (ppX, ppYr);
%!   ax = plotPartialDependence (Mdl, 1, ppLo, 'QueryPoints', ppQ, ...
%!                               'Conditional', 'absolute');
%!   ch = get (ax, 'children');
%!   red = [];
%!   grey = [];
%!   for k = 1:numel (ch)
%!     if (! strcmp (get (ch(k), 'type'), 'line'))
%!       continue;
%!     endif
%!     if (isequal (get (ch(k), 'color'), [1, 0, 0]))
%!       red = get (ch(k), 'ydata');
%!     else
%!       grey(end+1,:) = get (ch(k), 'ydata');
%!     endif
%!   endfor
%!   assert_equal (red, mean (grey, 1), 1e-12);
%!   assert_equal (red, [17, 17, 17], 1e-12);
%!   assert_equal (partialDependence (Mdl, 1, ppLo, 'QueryPoints', ppQ), ...
%!                 [27, 27, 27], 1e-12);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # the titles say which plot it is
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   Mdl = fitrsvm (ppX, ppYr);
%!   ax = plotPartialDependence (Mdl, 1, 'QueryPoints', ppQ);
%!   assert_equal (get (get (ax, 'title'), 'string'), ...
%!                 'Partial Dependence Plot');
%!   clf (hf);
%!   ax = plotPartialDependence (Mdl, 1, 'QueryPoints', ppQ, ...
%!                               'Conditional', 'centered');
%!   assert_equal (get (get (ax, 'title'), 'string'), ...
%!                 'Individual Conditional Expectation Plot');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # the horizontal axis is named after the predictor varied
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   Mdl = fitrsvm (ppX, ppYr, 'PredictorNames', {'alpha', 'beta'});
%!   ax = plotPartialDependence (Mdl, 'beta', 'QueryPoints', [2;4]);
%!   assert_equal (get (get (ax, 'xlabel'), 'string'), 'beta');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # MATLAB parity: Parent draws into the axes given and returns it
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   a1 = subplot (1, 2, 1);
%!   a2 = subplot (1, 2, 2);
%!   Mdl = fitrsvm (ppX, ppYr);
%!   ax = plotPartialDependence (Mdl, 1, 'QueryPoints', ppQ, 'Parent', a2);
%!   assert_equal (ax, a2);
%!   assert_equal (isempty (get (a2, 'children')), false);
%!   assert_equal (isempty (get (a1, 'children')), true);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Input validation
%!error<plotPartialDependence: too few input arguments.> ...
%! plotPartialDependence (1)

%!error<plotPartialDependence: a conditional plot varies one predictor, not two.> ...
%! plotPartialDependence (fitrsvm (ppX, ppYr), [1, 2], ...
%!                        'Conditional', 'absolute')

%!error<plotPartialDependence: a conditional plot of a classification model names one class.> ...
%! plotPartialDependence (fitctree (ppX, ppY3), 1, {'a', 'b'}, ...
%!                        'Conditional', 'absolute')

%!error<plotPartialDependence: 'Conditional' must be 'none', 'absolute' or 'centered'.> ...
%! plotPartialDependence (fitrsvm (ppX, ppYr), 1, 'Conditional', 'sideways')

%!error<plotPartialDependence: 'Parent' must be an axes handle.> ...
%! plotPartialDependence (fitrsvm (ppX, ppYr), 1, 'Parent', 42)

%!error<partialDependence: 'Conditional' and 'Parent' belong to plotPartialDependence.> ...
%! partialDependence (fitrsvm (ppX, ppYr), 1, 'Conditional', 'absolute')

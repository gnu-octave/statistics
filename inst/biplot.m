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
## @deftypefn  {statistics} {} biplot (@var{coefs})
## @deftypefnx {statistics} {} biplot (@var{coefs}, @var{name}, @var{value}, @dots{})
## @deftypefnx {statistics} {} biplot (@var{ax}, @dots{})
## @deftypefnx {statistics} {@var{h} =} biplot (@dots{})
##
## Create a biplot of the coefficients in @var{coefs}.
##
## @code{biplot (@var{coefs})} plots the rows of @var{coefs}, typically the
## principal component coefficients returned by @code{pca} or @code{pcacov} or
## the factor loadings returned by @code{factoran}, as vectors from the origin.
## @var{coefs} has one row per observed variable and either two columns (for a
## 2-D biplot) or three columns (for a 3-D biplot).
##
## The following name/value pairs are accepted:
##
## @table @asis
## @item @qcode{"Scores"}
## A matrix of scores with the same number of columns as @var{coefs} (one row
## per observation).  The scores are plotted as points, scaled to fit within the
## unit circle relative to the maximum coefficient length: each score is divided
## by the maximum absolute score value and multiplied by the length of the
## longest coefficient vector.
##
## @item @qcode{"VarLabels"}
## A character array or cell array of strings labeling each variable vector.
##
## @item @qcode{"ObsLabels"}
## A character array or cell array of strings labeling each observation.
##
## @item @qcode{"Positive"}
## If @qcode{true}, the reference axes are drawn only over the positive quadrant
## (2-D) or octant (3-D).  The default is @qcode{false}.
## @end table
##
## Any additional name/value pairs are treated as line properties and applied to
## the variable vectors.
##
## For readability the sign of each column of @var{coefs} is chosen so that its
## largest-magnitude element is positive; the same sign change is applied to the
## corresponding column of the scores.
##
## @code{biplot (@var{ax}, @dots{})} plots into the axes @var{ax}.
##
## The optional output @var{h} is a column vector of handles to the plotted
## graphics objects, ordered as the variable vector lines, the variable markers,
## the variable text labels (if any), the observation markers, the observation
## text labels (if any), and finally the reference axis lines.
##
## @seealso{pca, pcacov, factoran, rotatefactors}
## @end deftypefn

function h = biplot (varargin)

  ## Optional leading axes handle
  hax = [];
  if (numel (varargin) > 0 && isaxes (varargin{1}))
    hax = varargin{1};
    varargin(1) = [];
  endif

  if (numel (varargin) < 1)
    print_usage ();
  endif

  coefs = varargin{1};
  varargin(1) = [];
  if (! isnumeric (coefs) || ! isreal (coefs) || ndims (coefs) > 2 ...
      || (columns (coefs) != 2 && columns (coefs) != 3))
    error ("biplot: COEFS must be a real matrix with 2 or 3 columns.");
  endif
  p = rows (coefs);
  d = columns (coefs);

  ## Parse name/value options
  scores = [];
  varlabels = {};
  obslabels = {};
  positive = false;
  lineprops = {};
  if (mod (numel (varargin), 2) != 0)
    error ("biplot: name/value arguments must come in pairs.");
  endif
  for i = 1:2:numel (varargin)
    name = varargin{i};
    value = varargin{i+1};
    if (! ischar (name))
      error ("biplot: property names must be strings.");
    endif
    switch (lower (name))
      case "scores"
        scores = value;
      case "varlabels"
        varlabels = value;
      case "obslabels"
        obslabels = value;
      case "positive"
        positive = value;
      otherwise
        lineprops(end+1:end+2) = {name, value};
    endswitch
  endfor

  have_scores = ! isempty (scores);
  if (have_scores)
    if (! isnumeric (scores) || ! isreal (scores) || columns (scores) != d)
      error ("biplot: SCORES must be a real matrix with the same number of columns as COEFS.");
    endif
  endif
  if (ischar (varlabels))
    varlabels = cellstr (varlabels);
  endif
  if (ischar (obslabels))
    obslabels = cellstr (obslabels);
  endif

  ## Sign convention: force the largest-magnitude element of each column positive
  for j = 1:d
    [~, idx] = max (abs (coefs(:,j)));
    if (coefs(idx,j) < 0)
      coefs(:,j) = -coefs(:,j);
      if (have_scores)
        scores(:,j) = -scores(:,j);
      endif
    endif
  endfor

  ## Scale the scores to the coefficient space
  if (have_scores)
    maxlen = max (sqrt (sum (coefs .^ 2, 2)));
    maxscore = max (abs (scores(:)));
    if (maxscore > 0)
      scores = scores * (maxlen / maxscore);
    endif
  endif

  ## Reference axis extent
  axlim = 1.1 * max (abs (coefs(:)));
  if (axlim == 0 || isnan (axlim))
    axlim = 1;
  endif

  if (isempty (hax))
    hax = newplot ();
  else
    newplot (hax);
  endif
  old_hold = ishold (hax);
  hold (hax, "on");

  vcol = [0 0 1];      # variable vectors in blue
  ocol = [1 0 0];      # observations in red
  acol = [0.5 0.5 0.5]; # reference axes in gray

  varlines = zeros (p, 1);
  varmarks = zeros (p, 1);
  vartext = zeros (numel (varlabels), 1);
  obsmarks = zeros (0, 1);
  obstext = zeros (numel (obslabels), 1);

  ## Variable vector lines (origin to each coefficient)
  for i = 1:p
    if (d == 2)
      varlines(i) = line (hax, [0, coefs(i,1)], [0, coefs(i,2)], ...
                          "color", vcol, "marker", "none", lineprops{:});
    else
      varlines(i) = line (hax, [0, coefs(i,1)], [0, coefs(i,2)], ...
                          [0, coefs(i,3)], "color", vcol, "marker", "none", ...
                          lineprops{:});
    endif
  endfor

  ## Variable markers at the vector tips
  for i = 1:p
    if (d == 2)
      varmarks(i) = line (hax, coefs(i,1), coefs(i,2), "linestyle", "none", ...
                          "marker", "o", "color", vcol);
    else
      varmarks(i) = line (hax, coefs(i,1), coefs(i,2), coefs(i,3), ...
                          "linestyle", "none", "marker", "o", "color", vcol);
    endif
  endfor

  ## Variable text labels
  for i = 1:numel (varlabels)
    if (d == 2)
      vartext(i) = text (hax, coefs(i,1), coefs(i,2), varlabels{i}, ...
                         "color", vcol);
    else
      vartext(i) = text (hax, coefs(i,1), coefs(i,2), coefs(i,3), ...
                         varlabels{i}, "color", vcol);
    endif
  endfor

  ## Observation markers
  if (have_scores)
    m = rows (scores);
    obsmarks = zeros (m, 1);
    for j = 1:m
      if (d == 2)
        obsmarks(j) = line (hax, scores(j,1), scores(j,2), ...
                            "linestyle", "none", "marker", ".", "color", ocol);
      else
        obsmarks(j) = line (hax, scores(j,1), scores(j,2), scores(j,3), ...
                            "linestyle", "none", "marker", ".", "color", ocol);
      endif
    endfor
    ## Observation text labels
    for j = 1:numel (obslabels)
      if (d == 2)
        obstext(j) = text (hax, scores(j,1), scores(j,2), obslabels{j}, ...
                           "color", ocol);
      else
        obstext(j) = text (hax, scores(j,1), scores(j,2), scores(j,3), ...
                           obslabels{j}, "color", ocol);
      endif
    endfor
  endif

  ## Reference axis lines through the origin
  lo = ifelse (positive, 0, -axlim);
  if (d == 2)
    axline = line (hax, [lo, axlim, NaN, 0, 0], [0, 0, NaN, lo, axlim], ...
                   "color", acol);
  else
    axline = line (hax, [lo, axlim, NaN, 0, 0, NaN, 0, 0], ...
                   [0, 0, NaN, lo, axlim, NaN, 0, 0], ...
                   [0, 0, NaN, 0, 0, NaN, lo, axlim], "color", acol);
  endif

  if (d == 3)
    view (hax, 3);
    zlabel (hax, "Component 3");
  endif
  xlabel (hax, "Component 1");
  ylabel (hax, "Component 2");

  if (! old_hold)
    hold (hax, "off");
  endif

  if (nargout > 0)
    h = [varlines; varmarks; vartext; obsmarks; obstext; axline];
  endif

endfunction

## Return a if cond is true, else b.
function r = ifelse (cond, a, b)
  if (cond)
    r = a;
  else
    r = b;
  endif
endfunction

%!demo
%! ## Biplot of the first two principal components of Fisher's iris data.
%!
%! load fisheriris;
%! [coefs, score] = pca (zscore (meas));
%! biplot (coefs(:,1:2), "Scores", score(:,1:2), ...
%!         "VarLabels", {"SL", "SW", "PL", "PW"});

%!demo
%! ## Three-component biplot of the same data.
%!
%! load fisheriris;
%! coefs = pca (zscore (meas));
%! biplot (coefs(:,1:3), "VarLabels", {"SL", "SW", "PL", "PW"});

## Test output
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   coefs = [0.6 -0.3; -0.2 0.7; 0.5 0.5];
%!   score = [1 2; -3 1; 0.5 -2; 4 0];
%!   h = biplot (coefs, "Scores", score);
%!   assert (numel (h), 11);
%!   ## variable vector lines (origin to coefficient)
%!   assert (get (h(1), "xdata"), [0 0.6], 1e-12);
%!   assert (get (h(1), "ydata"), [0 -0.3], 1e-12);
%!   assert (get (h(2), "xdata"), [0 -0.2], 1e-12);
%!   ## variable tip markers
%!   assert (get (h(4), "xdata"), 0.6, 1e-12);
%!   assert (get (h(4), "ydata"), -0.3, 1e-12);
%!   ## observation markers (scaled scores)
%!   assert (get (h(7), "xdata"), 0.182, 1e-4);
%!   assert (get (h(7), "ydata"), 0.364, 1e-4);
%!   assert (get (h(10), "xdata"), 0.728, 1e-4);
%!   ## reference axis extent
%!   assert (get (h(11), "xdata"), [-0.77 0.77 NaN 0 0], 1e-12);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test  # 3-D biplot without scores
%! hf = figure ("visible", "off");
%! unwind_protect
%!   coefs = [0.6 -0.3 0.2; -0.2 0.7 0.1; 0.5 0.5 -0.6];
%!   h = biplot (coefs);
%!   assert (numel (h), 7);
%!   ## column 3 has its largest-magnitude element (-0.6) forced positive,
%!   ## so the whole column is negated: [0.2 0.1 -0.6] -> [-0.2 -0.1 0.6]
%!   assert (get (h(1), "zdata"), [0 -0.2], 1e-12);
%!   assert (get (h(6), "zdata"), 0.6, 1e-12);
%!   assert (get (h(7), "zdata"), [0 0 NaN 0 0 NaN -0.77 0.77], 1e-12);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test  # sign convention flips a column with a negative largest element
%! hf = figure ("visible", "off");
%! unwind_protect
%!   coefs = [-0.8 0.1; 0.2 0.9; -0.5 0.5];
%!   h = biplot (coefs);
%!   assert (get (h(1), "xdata"), [0 0.8], 1e-12);
%!   assert (get (h(2), "xdata"), [0 -0.2], 1e-12);
%!   assert (get (h(3), "xdata"), [0 0.5], 1e-12);
%!   assert (get (h(1), "ydata"), [0 0.1], 1e-12);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test  # labels add text handles in the documented order
%! hf = figure ("visible", "off");
%! unwind_protect
%!   coefs = [0.6 -0.3; -0.2 0.7];
%!   score = [1 2; -3 1];
%!   h = biplot (coefs, "Scores", score, "VarLabels", {"a", "b"}, ...
%!               "ObsLabels", {"x", "y"});
%!   ## 2 varlines + 2 varmarks + 2 vartext + 2 obsmarks + 2 obstext + 1 axis
%!   assert (numel (h), 11);
%!   assert (strcmp (get (h(5), "string"), "a"));
%!   assert (strcmp (get (h(9), "string"), "x"));
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Test input validation
%!error <Invalid call to biplot> biplot ()
%!error <biplot: COEFS must be a real matrix with 2 or 3 columns.> ...
%! biplot (ones (3, 4))
%!error <biplot: COEFS must be a real matrix with 2 or 3 columns.> biplot ({1})
%!error <biplot: SCORES must be a real matrix with the same number of columns> ...
%! biplot ([0.6 0.3; 0.2 0.7], "Scores", [1 2 3])
%!error <biplot: name/value arguments must come in pairs.> ...
%! biplot ([0.6 0.3; 0.2 0.7], "Scores")

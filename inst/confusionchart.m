## Copyright (C) 2020-2021 Stefano Guidoni <ilguido@users.sf.net>
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
## @deftypefn {} {} confusionchart (@var{trueLabels}, @var{predictedLabels})
## @deftypefnx {} {} confusionchart (@var{m})
## @deftypefnx {} {} confusionchart (@var{m}, @var{classLabels})
## @deftypefnx {} {} confusionchart (@var{parent}, @dots{})
## @deftypefnx {} {} confusionchart (@dots{}, @var{prop}, @var{val}, @dots{})
## @deftypefnx {} {@var{cm} =} confusionchart (@dots{})
##
## Display a chart of a confusion matrix.
##
## The two vectors of values @var{trueLabels} and @var{predictedLabels}, which
## are used to compute the confusion matrix, must be defined with the same
## format as the inputs of @code{confusionmat}.
## Otherwise a confusion matrix @var{m} as computed by @code{confusionmat} can
## be given.
##
## @var{classLabels} is an array of labels, i.e. the list of the class names.
##
## If the first argument is a handle to a @code{figure} or to a @code{uipanel},
## then the confusion matrix chart is displayed inside that object.
##
## Optional property/value pairs are passed directly to the underlying objects,
## e.g. @qcode{"xlabel"}, @qcode{"ylabel"}, @qcode{"title"}, @qcode{"fontname"},
## @qcode{"fontsize"} etc.
##
## The optional return value @var{cm} is a @code{ConfusionMatrixChart} object.
## Specific properties of a @code{ConfusionMatrixChart} object are:
## @itemize @bullet
## @item @qcode{"DiagonalColor"}
## The color of the patches on the diagonal, default is [0.0, 0.4471, 0.7412].
##
## @item @qcode{"OffDiagonalColor"}
## The color of the patches off the diagonal, default is [0.851, 0.3255, 0.098].
##
## @item @qcode{"GridVisible"}
## Available values: @qcode{on} (default), @qcode{off}.
##
## @item @qcode{"Normalization"}
## Available values: @qcode{absolute} (default), @qcode{column-normalized},
## @qcode{row-normalized}, @qcode{total-normalized}.
##
## @item @qcode{"ColumnSummary"}
## Available values: @qcode{off} (default), @qcode{absolute},
## @qcode{column-normalized},@qcode{total-normalized}.
##
## @item @qcode{"RowSummary"}
## Available values: @qcode{off} (default), @qcode{absolute},
## @qcode{row-normalized}, @qcode{total-normalized}.
## @end itemize
##
## Run @code{demo confusionchart} to see some examples.
##
## @end deftypefn
##
## @seealso{confusionmat, sortClasses}

## Author: Stefano Guidoni <ilguido@users.sf.net>

function cm = confusionchart (varargin)

  ## check the input parameters
  if (nargin < 1)
    print_usage ();
  endif

  p_i = 1;

  if (ishghandle (varargin{p_i}))
    ## parameter is a parent figure
    handle_type = get (varargin{p_i}, "type");
    if (strcmp (handle_type, "figure"))
      h = figure (varargin{p_i});
      hax = axes ("parent", h);
    elseif (strcmp (handle_type, "uipanel"))
      h = varargin{p_i};
      hax = axes ("parent", varargin{p_i});
    else
      ## MATLAB compatibility: on MATLAB are also available Tab objects,
      ## TiledChartLayout objects, GridLayout objects
      error ("confusionchart: invalid handle to parent object");
    endif
    p_i++;
  else
    h = figure ();
    hax = axes ("parent", h);
  endif

  if (ismatrix (varargin{p_i}) && rows (varargin{p_i}) == ...
      columns (varargin{p_i}))
    ## parameter is a confusion matrix
    conmat = varargin{p_i};
    p_i++;

    if (p_i <= nargin && ((isvector (varargin{p_i}) && ...
      length (varargin{p_i}) == rows (conmat)) || ...
      (ischar ( varargin{p_i}) && rows (varargin{p_i}) == rows (conmat)) ...
      || iscellstr (varargin{p_i})))
      ## parameter is an array of labels
      labarr = varargin{p_i};

      if (isrow (labarr))
        labarr = vec (labarr);
      endif

      p_i++;
    else
      labarr = [1 : (rows (conmat))]';
    endif
  elseif (isvector (varargin{p_i}))
    ## parameter must be a group for confusionmat
    [conmat, labarr] = confusionmat (varargin{p_i}, varargin{p_i + 1});
    p_i = p_i + 2;
  else
    close (h);
    error ("confusionchart: invalid argument");
  endif

  ## remaining parameters are stored
  i = p_i;
  args = {};
  while (i <= nargin)
    args{end + 1} = varargin{i++};
  endwhile

  ## prepare the labels
  if (! iscellstr (labarr))
    if (! ischar (labarr))
      labarr = num2str (labarr);
    endif

    labarr = cellstr (labarr);
  endif

  ## MATLAB compatibility: labels are sorted
  [labarr, I] = sort (labarr);
  conmat = conmat(I, :);
  conmat = conmat(:, I);

  cm = ConfusionMatrixChart (hax, conmat, labarr, args);

endfunction


## Test input validation

## Get current figure visibility so it can be restored after tests
%!shared visibility_setting
%! visibility_setting = get (0, "DefaultFigureVisible");

%!test
%! set (0, "DefaultFigureVisible", "off");
%! fail ("confusionchart ()", "Invalid call");
%! set (0, "DefaultFigureVisible", visibility_setting);

%!test
%! set (0, "DefaultFigureVisible", "off");
%! fail ("confusionchart ([1 1; 2 2; 3 3])", "invalid argument");
%! set (0, "DefaultFigureVisible", visibility_setting);

%!test
%! set (0, "DefaultFigureVisible", "off");
%! fail ("confusionchart ([1 2], [0 1], 'xxx', 1)", "invalid property");
%! set (0, "DefaultFigureVisible", visibility_setting);

%!test
%! set (0, "DefaultFigureVisible", "off");
%! fail ("confusionchart ([1 2], [0 1], 'XLabel', 1)", "XLabel .* string");
%! set (0, "DefaultFigureVisible", visibility_setting);

%!test
%! set (0, "DefaultFigureVisible", "off");
%! fail ("confusionchart ([1 2], [0 1], 'YLabel', [1 0])", ".* YLabel .* string");
%! set (0, "DefaultFigureVisible", visibility_setting);

%!test
%! set (0, "DefaultFigureVisible", "off");
%! fail ("confusionchart ([1 2], [0 1], 'Title', .5)", ".* Title .* string");
%! set (0, "DefaultFigureVisible", visibility_setting);

%!test
%! set (0, "DefaultFigureVisible", "off");
%! fail ("confusionchart ([1 2], [0 1], 'FontName', [])", ".* FontName .* string");
%! set (0, "DefaultFigureVisible", visibility_setting);

%!test
%! set (0, "DefaultFigureVisible", "off");
%! fail ("confusionchart ([1 2], [0 1], 'FontSize', 'b')", ".* FontSize .* numeric");
%! set (0, "DefaultFigureVisible", visibility_setting);

%!test
%! set (0, "DefaultFigureVisible", "off");
%! fail ("confusionchart ([1 2], [0 1], 'DiagonalColor', 'h')", ".* DiagonalColor .* color");
%! set (0, "DefaultFigureVisible", visibility_setting);

%!test
%! set (0, "DefaultFigureVisible", "off");
%! fail ("confusionchart ([1 2], [0 1], 'OffDiagonalColor', [])", ".* OffDiagonalColor .* color");
%! set (0, "DefaultFigureVisible", visibility_setting);

%!test
%! set (0, "DefaultFigureVisible", "off");
%! fail ("confusionchart ([1 2], [0 1], 'Normalization', '')", ".* invalid .* Normalization");
%! set (0, "DefaultFigureVisible", visibility_setting);

%!test
%! set (0, "DefaultFigureVisible", "off");
%! fail ("confusionchart ([1 2], [0 1], 'ColumnSummary', [])", ".* invalid .* ColumnSummary");
%! set (0, "DefaultFigureVisible", visibility_setting);

%!test
%! set (0, "DefaultFigureVisible", "off");
%! fail ("confusionchart ([1 2], [0 1], 'RowSummary', 1)", ".* invalid .* RowSummary");
%! set (0, "DefaultFigureVisible", visibility_setting);

%!test
%! set (0, "DefaultFigureVisible", "off");
%! fail ("confusionchart ([1 2], [0 1], 'GridVisible', .1)", ".* invalid .* GridVisible");
%! set (0, "DefaultFigureVisible", visibility_setting);

%!test
%! set (0, "DefaultFigureVisible", "off");
%! fail ("confusionchart ([1 2], [0 1], 'HandleVisibility', .1)", ".* invalid .* HandleVisibility");
%! set (0, "DefaultFigureVisible", visibility_setting);

%!test
%! set (0, "DefaultFigureVisible", "off");
%! fail ("confusionchart ([1 2], [0 1], 'OuterPosition', .1)", ".* invalid .* OuterPosition");
%! set (0, "DefaultFigureVisible", visibility_setting);

%!test
%! set (0, "DefaultFigureVisible", "off");
%! fail ("confusionchart ([1 2], [0 1], 'Position', .1)", ".* invalid .* Position");
%! set (0, "DefaultFigureVisible", visibility_setting);

%!test
%! set (0, "DefaultFigureVisible", "off");
%! fail ("confusionchart ([1 2], [0 1], 'Units', .1)", ".* invalid .* Units");
%! set (0, "DefaultFigureVisible", visibility_setting);

## Demonstration using the confusion matrix example from
## R.Bonnin, "Machine Learning for Developers", pp. 55-56
%!demo "Setting the chart properties"
%! Yt = [8 5 6 8 5 3 1 6 4 2 5 3 1 4]';
%! Yp = [8 5 6 8 5 2 3 4 4 5 5 7 2 6]';
%! confusionchart (Yt, Yp, "Title", ...
%!   "Demonstration with summaries","Normalization",...
%!   "absolute","ColumnSummary", "column-normalized","RowSummary",...
%!   "row-normalized")

## example: confusion matrix and class labels
%!demo "Cellstr as inputs"
%! Yt = {'Positive', 'Positive', 'Positive', 'Negative', 'Negative' };
%! Yp = {'Positive', 'Positive', 'Negative', 'Negative', 'Negative' };
%! m = confusionmat ( Yt, Yp );
%! confusionchart ( m, { 'Positive', 'Negative' } );

## example: editing the properties of an existing ConfusionMatrixChart object
%!demo "Editing the object properties"
%! Yt = {'Positive', 'Positive', 'Positive', 'Negative', 'Negative' };
%! Yp = {'Positive', 'Positive', 'Negative', 'Negative', 'Negative' };
%! cm = confusionchart ( Yt, Yp );
%! cm.Title = "This is an example with a green diagonal";
%! cm.DiagonalColor = [0.4660, 0.6740, 0.1880];

## example: drawing the chart inside a uipanel
%!demo "Confusion chart in a uipanel"
%! h = uipanel ();
%! Yt = {'Positive', 'Positive', 'Positive', 'Negative', 'Negative' };
%! Yp = {'Positive', 'Positive', 'Negative', 'Negative', 'Negative' };
%! cm = confusionchart ( h, Yt, Yp );

## example: sortClasses
%!demo "Sorting classes"
%! Yt = [8 5 6 8 5 3 1 6 4 2 5 3 1 4]';
%! Yp = [8 5 6 8 5 2 3 4 4 5 5 7 2 6]';
%! cm = confusionchart (Yt, Yp, "Title", ...
%!   "Classes are sorted according to clusters");
%! sortClasses (cm, "cluster");

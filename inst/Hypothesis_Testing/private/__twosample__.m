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
## @deftypefn  {Private Function} {[@var{Z}, @var{iscat}, @var{mx}, @var{my}, @var{errmsg}, @var{names}, @var{labels}] =} __twosample__ (@var{X}, @var{Y}, @var{vnames}, @var{catvars})
## @deftypefnx {Private Function} {[@dots{}] =} __twosample__ (@dots{}, @var{argnames})
##
## Read the two samples of @code{knntest} and @code{mmdtest}.
##
## @var{X} and @var{Y} are numeric matrices or tables, @var{vnames} the
## @qcode{'VariableNames'} given and @var{catvars} the
## @qcode{'CategoricalVariables'} given, each empty where not.  @var{Z} holds
## the pooled sample, the @var{mx} rows of @var{X} over the @var{my} rows of
## @var{Y}, one column per variable used: a continuous variable as numbers,
## one holding levels as codes @code{1}, @code{2}, @dots{} shared across the
## samples.  A row holding a missing value is left out.  @var{iscat} says
## which columns hold levels.  @var{errmsg} is empty, or the body of the
## message the caller raises under its own name.  @var{names} holds the
## variable names of two tables and is empty for matrices, and @var{labels}
## holds, for each variable holding levels, the name of each code in turn.
## @var{argnames}, @code{@{'X', 'Y'@}} by default, names the two samples in
## @var{errmsg}.
##
## @end deftypefn

function [Z, iscat, mx, my, errmsg, names, labels] = __twosample__ (X, Y, ...
                                                  vnames, catvars, argnames)

  [Z, iscat, mx, my, errmsg, names, labels] = tsRead (X, Y, vnames, catvars);
  if (! isempty (errmsg) && nargin > 4)
    errmsg = regexprep (errmsg, '(?<![A-Za-z])X(?![A-Za-z])', argnames{1});
    errmsg = regexprep (errmsg, '(?<![A-Za-z])Y(?![A-Za-z])', argnames{2});
  endif

endfunction

## The samples read with the messages naming them X and Y.
function [Z, iscat, mx, my, errmsg, names, labels] = tsRead (X, Y, ...
                                                             vnames, catvars)

  Z = [];
  mx = 0;
  my = 0;
  labels = {};
  [xc, yc, names, iscat, errmsg] = tsVariables (X, Y, vnames);
  if (! isempty (errmsg))
    return;
  endif
  [iscat, errmsg] = tsCategorical (catvars, iscat, names, istable (X));
  if (! isempty (errmsg))
    return;
  endif
  [xc, mx] = tsComplete (xc);
  [yc, my] = tsComplete (yc);
  if (mx < 1 || my < 1)
    errmsg = strcat ("X and Y must each hold an observation with no", ...
                     " missing value.");
    return;
  endif
  K = numel (xc);
  Z = zeros (mx + my, K);
  labels = cell (1, K);
  for j = 1:K
    if (iscat(j))
      ## Codes 1, 2, ... over the levels present, each with its name
      [g, gn] = grp2idx (vertcat (xc{j}, yc{j}));
      [u, ~, Z(:,j)] = unique (g);
      labels{j} = gn(u);
    else
      Z(:,j) = double (vertcat (xc{j}, yc{j}));
    endif
  endfor

endfunction

## The columns of X and Y, one per variable used, the variable names, and
## which variables hold levels by their type.
function [xc, yc, names, iscat, errmsg] = tsVariables (X, Y, vnames)

  xc = {};
  yc = {};
  names = {};
  iscat = [];
  errmsg = '';
  if (istable (X) != istable (Y))
    errmsg = "X and Y must both be matrices or both be tables.";
    return;
  endif

  if (! istable (X))
    if (! isempty (vnames))
      errmsg = "'VariableNames' applies only where X and Y are tables.";
      return;
    endif
    if (! (isnumeric (X) && isreal (X) && ismatrix (X)
           && isnumeric (Y) && isreal (Y) && ismatrix (Y)))
      errmsg = "X and Y must be real numeric matrices or tables.";
      return;
    endif
    if (columns (X) != columns (Y) || columns (X) < 1)
      errmsg = "X and Y must have the same number of columns.";
      return;
    endif
    K = columns (X);
    xc = cell (1, K);
    yc = cell (1, K);
    for j = 1:K
      xc{j} = X(:,j);
      yc{j} = Y(:,j);
    endfor
    names = cell (1, K);
    iscat = false (1, K);
    return;
  endif

  xn = X.Properties.VariableNames;
  yn = Y.Properties.VariableNames;
  if (all (ismember (yn, xn)))
    shared = yn;
  elseif (all (ismember (xn, yn)))
    shared = xn;
  else
    errmsg = strcat ("the variable names of one of X and Y must all be", ...
                     " names of the other.");
    return;
  endif
  names = shared;
  if (! isempty (vnames))
    if (isa (vnames, 'string'))
      vnames = cellstr (vnames);
    elseif (ischar (vnames) && isrow (vnames))
      vnames = {vnames};
    endif
    if (! (iscellstr (vnames) && all (ismember (vnames, shared))))
      errmsg = "'VariableNames' must name variables X and Y share.";
      return;
    endif
    names = vnames(:)';
  endif
  K = numel (names);
  xc = cell (1, K);
  yc = cell (1, K);
  iscat = false (1, K);
  for j = 1:K
    a = X.(names{j});
    c = Y.(names{j});
    if (columns (a) != 1 || columns (c) != 1)
      errmsg = sprintf ("variable '%s' must be one column.", names{j});
      return;
    endif
    [ka, oka] = tsKind (a);
    [kc, okc] = tsKind (c);
    if (! (oka && okc))
      errmsg = sprintf (strcat ("variable '%s' must hold numbers, logical", ...
                                " values, categories or text."), names{j});
      return;
    endif
    if (ka != kc)
      errmsg = sprintf (strcat ("variable '%s' must hold the same kind of", ...
                                " values in X and Y."), names{j});
      return;
    endif
    iscat(j) = ka;
    if (isa (a, 'categorical') && ! ka)
      ## An ordinal category is read as its position among the categories
      a = double (grp2idx (a));
      c = double (grp2idx (c));
    endif
    xc{j} = a;
    yc{j} = c;
  endfor

endfunction

## Whether a table variable holds levels by its type, and whether its type is
## one this test reads at all.
function [iscat, ok] = tsKind (v)
  ok = true;
  if (islogical (v) || isa (v, 'string') || iscellstr (v))
    iscat = true;
  elseif (isa (v, 'categorical'))
    iscat = ! isordinal (v);
  elseif (isnumeric (v) && isreal (v))
    iscat = false;
  else
    iscat = false;
    ok = false;
  endif
endfunction

## The variables named in 'CategoricalVariables', added to those holding levels
## by their type.
function [iscat, errmsg] = tsCategorical (cv, iscat, names, istab)

  errmsg = '';
  if (isempty (cv))
    return;
  endif
  K = numel (iscat);
  if (isa (cv, 'string'))
    cv = cellstr (cv);
  endif
  if (iscellstr (cv) && numel (cv) == 1)
    cv = cv{1};
  endif
  if (ischar (cv) && strcmpi (cv, 'all'))
    iscat(:) = true;
  elseif (ischar (cv) || iscellstr (cv))
    if (! istab)
      errmsg = strcat ("'CategoricalVariables' can name variables only", ...
                       " where X and Y are tables.");
      return;
    endif
    cv = cellstr (cv);
    if (! all (ismember (cv, names)))
      errmsg = "'CategoricalVariables' must name variables in use.";
      return;
    endif
    iscat(ismember (names, cv)) = true;
  elseif (islogical (cv))
    if (numel (cv) != K)
      errmsg = sprintf (strcat ("'CategoricalVariables' must hold one", ...
                                " logical value for each of the %d", ...
                                " variables."), K);
      return;
    endif
    iscat(cv(:)') = true;
  elseif (isnumeric (cv) && isreal (cv) && all (cv(:) == fix (cv(:)))
          && all (cv(:) >= 1) && all (cv(:) <= K))
    iscat(cv(:)') = true;
  else
    errmsg = sprintf (strcat ("'CategoricalVariables' must be 'all',", ...
                              " indices from 1 to %d, a logical vector", ...
                              " or variable names."), K);
  endif

endfunction

## The columns with every observation holding a missing value removed.
function [cols, n] = tsComplete (cols)
  keep = true (rows (cols{1}), 1);
  for j = 1:numel (cols)
    keep &= ! ismissing (cols{j});
  endfor
  for j = 1:numel (cols)
    cols{j} = cols{j}(keep);
  endfor
  n = sum (keep);
endfunction

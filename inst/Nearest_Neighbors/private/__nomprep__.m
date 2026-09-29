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
## @deftypefn {Private Function} {[@var{C}, @var{measure}, @var{w}, @var{errmsg}] =} __nomprep__ (@var{X}, @var{Y}, @var{hasY}, @var{args})
##
## Validate and code the input of @code{nomdist} and @code{nomdist2}.
##
## @var{C} holds the level codes of @var{X}, @var{Y} and the reference sample,
## positive integers shared across the three and one column per variable,
## the last empty where no @qcode{'Reference'} was given.  @var{measure} is
## the measure named in lower case, @qcode{'goodall3'} where none was, and
## @var{w} a row of one weight per variable.  @var{errmsg} is empty, or the
## body of the message the caller raises under its own name.
##
## @var{args} holds what followed @var{X} and @var{Y}: a measure where their
## number is odd, then @qcode{'Weights'} and @qcode{'Reference'} pairs.
##
## @end deftypefn

function [C, measure, w, errmsg] = __nomprep__ (X, Y, hasY, args)

  C = {};
  w = [];
  errmsg = '';
  names = {'anderberg', 'burnaby', 'eskin', 'gambaryan', 'goodall1', ...
           'goodall2', 'goodall3', 'goodall4', 'iof', 'lin', 'lin1', 'of', ...
           'sm', 'smirnov', 've', 'vm'};

  ## The measure, which an odd number of arguments leads with
  measure = 'goodall3';
  if (mod (numel (args), 2) == 1)
    m = args{1};
    if (isa (m, 'string') && isscalar (m))
      m = char (m);
    endif
    if (! (ischar (m) && isrow (m) && any (strcmpi (m, names))))
      errmsg = sprintf ("MEASURE must be one of %s.", ...
                        strjoin (strcat ("'", names, "'"), ', '));
      return;
    endif
    measure = lower (m);
    args(1) = [];
  endif

  ## An explicit empty value counts as not given
  [w, R, rem] = parsePairedArguments ({'Weights', 'Reference'}, {[], []}, ...
                                      args(:));
  if (! isempty (rem))
    errmsg = "invalid optional paired argument.";
    return;
  endif
  hasR = ! isempty (R);

  ## One column of levels per variable, a table read by X's variable names
  [xc, vars, errmsg] = nomColumns (X, 'X', X, {});
  if (! isempty (errmsg))
    return;
  endif
  if (size (X, 1) < 1 || isempty (xc))
    errmsg = "X must not be empty.";
    return;
  endif
  K = numel (xc);
  yc = {};
  rc = {};
  if (hasY)
    [yc, ~, errmsg] = nomColumns (Y, 'Y', X, vars);
    if (! isempty (errmsg))
      return;
    endif
    if (numel (yc) != K)
      errmsg = "Y must have as many columns as X.";
      return;
    endif
  endif
  if (hasR)
    [rc, ~, errmsg] = nomColumns (R, "'Reference'", X, vars);
    if (! isempty (errmsg))
      return;
    endif
    if (numel (rc) != K)
      errmsg = "'Reference' must have as many columns as X.";
      return;
    endif
  endif

  ## Codes shared across the samples, one variable at a time
  nx = rows (xc{1});
  ny = 0;
  nr = 0;
  if (hasY)
    ny = rows (yc{1});
  endif
  if (hasR)
    nr = rows (rc{1});
  endif
  Xc = zeros (nx, K);
  Yc = zeros (ny, K);
  Rc = zeros (nr, K);
  for k = 1:K
    f = levelType (xc{k});
    parts = {xc{k}};
    if (hasY)
      if (! strcmp (levelType (yc{k}), f))
        errmsg = "Y must hold the same type as X in every variable.";
        return;
      endif
      parts{end+1} = yc{k};
    endif
    if (hasR)
      if (! strcmp (levelType (rc{k}), f))
        errmsg = "'Reference' must hold the same type as X in every variable.";
        return;
      endif
      parts{end+1} = rc{k};
    endif
    if (strcmp (f, 'numeric'))
      parts = cellfun (@double, parts, 'UniformOutput', false);
    endif
    g = grp2idx (vertcat (parts{:}));
    Xc(:,k) = g(1:nx);
    Yc(:,k) = g(nx+1:nx+ny);
    Rc(:,k) = g(nx+ny+1:end);
  endfor
  if (any (isnan (Xc(:))))
    errmsg = "X must not hold missing values.";
    return;
  endif
  if (any (isnan (Yc(:))))
    errmsg = "Y must not hold missing values.";
    return;
  endif
  if (any (isnan (Rc(:))))
    errmsg = "'Reference' must not hold missing values.";
    return;
  endif
  if (hasR)
    for k = 1:K
      if (! all (ismember (Xc(:,k), Rc(:,k))))
        errmsg = "'Reference' must hold every level that X holds.";
        return;
      endif
    endfor
  else
    Rc = [];
  endif
  C = {Xc, Yc, Rc};

  ## The weights, read by every measure but three
  if (isempty (w))
    w = ones (1, K);
    return;
  endif
  if (any (strcmp (measure, {'anderberg', 'gambaryan', 'smirnov'})))
    errmsg = sprintf ("'Weights' does not apply to the '%s' measure.", ...
                      measure);
    return;
  endif
  if (! (isnumeric (w) && isreal (w) && isvector (w) && numel (w) == K
         && all (isfinite (w)) && all (w >= 0 & w <= 1)))
    errmsg = sprintf (strcat ("'Weights' must be a vector of %d values", ...
                              " between 0 and 1."), K);
    return;
  endif
  if (! any (w > 0))
    errmsg = "'Weights' must hold at least one positive value.";
    return;
  endif
  w = double (w(:)');

endfunction

## The columns of one sample as a cell of column vectors.  X decides whether
## the samples are tables, and a table is read by X's variable names.
function [cols, vars, errmsg] = nomColumns (A, name, X, vars)

  cols = {};
  errmsg = '';
  if (istable (X))
    if (! istable (A))
      errmsg = sprintf ("%s must be a table, as X is.", name);
      return;
    endif
    if (isempty (vars))
      vars = A.Properties.VariableNames;
    elseif (! all (ismember (vars, A.Properties.VariableNames)))
      errmsg = sprintf ("%s must hold every variable of X.", name);
      return;
    endif
    cols = cell (1, numel (vars));
    for k = 1:numel (vars)
      v = A.(vars{k});
      if (! (isLevels (v) && columns (v) == 1 && ndims (v) == 2))
        errmsg = sprintf (strcat ("every variable of %s must be one", ...
                                  " column of levels."), name);
        return;
      endif
      cols{k} = v;
    endfor
    return;
  endif

  if (istable (A))
    errmsg = sprintf ("%s must be a matrix, as X is.", name);
    return;
  endif
  if (! isLevels (A) || ndims (A) != 2)
    errmsg = sprintf (strcat ("%s must be a numeric, logical, categorical,", ...
                              " string or cellstr matrix, or a table."), name);
    return;
  endif
  cols = cell (1, columns (A));
  for k = 1:columns (A)
    cols{k} = A(:,k);
  endfor

endfunction

## Whether an array can hold levels.  A char matrix cannot: its rows would
## be the observations and its columns meaningless.
function tf = isLevels (v)
  tf = (isnumeric (v) && isreal (v)) || islogical (v) || iscellstr (v) ...
       || isa (v, 'categorical') || isa (v, 'string');
endfunction

## The kind of level a column holds, which must agree across the samples.
function f = levelType (v)
  if (isnumeric (v))
    f = 'numeric';
  elseif (iscellstr (v))
    f = 'cellstr';
  else
    f = class (v);
  endif
endfunction

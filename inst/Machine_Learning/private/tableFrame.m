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
## @deftypefn {private} {[@var{X}, @var{Y}, @var{args}, @var{lev}, @var{errmsg}] =} tableFrame (@var{X}, @var{Y}, @var{args})
##
## Resolve table input into the predictors and the response a learner takes.
##
## @var{X} is the predictor data, a table or a numeric matrix, and @var{Y}
## either the response, the name of the variable holding it, a model
## formula naming the response and the predictors together, or empty where
## there is no response and every variable is a predictor,
## @qcode{'Y ~ x1 + x2'}.  A formula holds main effects only: no wildcard,
## no products and no powers, as R2024a takes none of them.  @var{args} are
## the name-value arguments the call carried.
##
## Where @var{X} is not a table nothing is resolved and the three come back as
## they went in, with @var{lev} empty.  Where it is, @var{X} comes back as a
## real numeric matrix, @var{Y} as the response, and @var{args} carrying
## @qcode{'PredictorNames'}, @qcode{'ResponseName'} and
## @qcode{'CategoricalPredictors'} as the table gives them.  @var{lev} is one
## cell per predictor, holding the levels of a predictor read from text and
## empty for one read from numbers, which is what lets a level map to the code
## at prediction that it carried at fitting.
##
## This runs before any model exists, which is why it is here rather than on
## @code{PredictiveModel}: mapping a table onto a model that has already been
## fitted is that class's to do.
##
## @end deftypefn

function [X, Y, args, lev, errmsg] = tableFrame (X, Y, args)

  lev = {};
  errmsg = '';
  if (! istable (X))
    return;
  endif

  names = X.Properties.VariableNames;

  ## The response is named where it is one of the table's own variables,
  ## and given outright otherwise
  respname = 'Y';
  keep = true (1, numel (names));
  formula = {};
  if (ischar (Y) && isrow (Y) || (isa (Y, 'string') && isscalar (Y)))
    spec = char (Y);
    if (any (spec == '~'))
      [respname, formula, errmsg] = tfFormula (spec, names);
      if (! isempty (errmsg))
        return;
      endif
    else
      respname = spec;
      j = find (strcmp (names, respname), 1);
      if (isempty (j))
        ## A formula missing its tilde reads as a response name that no
        ## column carries, which is a poor way to be told of a typo
        if (any (ismember (spec, '+*:^')))
          errmsg = sprintf (strcat ("'%s' is not a model formula: a", ...
                                    " formula names its response before", ...
                                    " a '~'."), spec);
        else
          errmsg = sprintf ("the table holds no variable '%s'.", respname);
        endif
        return;
      endif
    endif
    Y = X.(respname);
    keep(strcmp (names, respname)) = false;
  elseif (isempty (Y))
    ## No response at all: every variable is a predictor.  An explainer
    ## works over predictor data alone, there being nothing to predict.
    Y = [];
  elseif (rows (Y) != height (X))
    errmsg = "the table must have one row per response.";
    return;
  endif

  ## Whatever is left is a predictor, unless the formula or the call named
  ## a subset; the two together are a conflict, as R2024a calls it one
  pnames = names(keep);
  [pn, args] = tfTakeArg (args, 'PredictorNames');
  if (! isempty (formula))
    if (! isempty (pn))
      errmsg = strcat ("'PredictorNames' cannot be given beside a model", ...
                       " formula, which names the predictors itself.");
      return;
    endif
    pnames = formula;
  elseif (! isempty (pn))
    [pnames, errmsg] = tfNameSubset (pn, pnames);
    if (! isempty (errmsg))
      return;
    endif
  endif
  if (isempty (pnames))
    errmsg = "the table holds no predictors.";
    return;
  endif

  ## A column of levels is a categorical predictor whether or not the call
  ## said so, and what the call said is added to that rather than replacing
  ## it, as R2024a adds it
  [cp, args] = tfTakeArg (args, 'CategoricalPredictors');
  M = numel (pnames);
  iscat = false (1, M);
  lev = cell (1, M);
  Xn = zeros (height (X), M);
  for k = 1:M
    v = X.(pnames{k});
    if (columns (v) != 1)
      errmsg = sprintf (strcat ("the table variable '%s' must hold one", ...
                                " column."), pnames{k});
      return;
    endif
    [Xn(:,k), lev{k}, iscat(k), errmsg] = tfColumn (v, pnames{k});
    if (! isempty (errmsg))
      return;
    endif
  endfor
  idx = find (iscat);
  if (! isempty (cp))
    [named, errmsg] = tfCatSpec (cp, pnames);
    if (! isempty (errmsg))
      return;
    endif
    idx = unique ([idx, named]);
  endif

  X = Xn;
  args = [{'PredictorNames', pnames, 'ResponseName', respname, ...
           'CategoricalPredictors', idx}, args];

endfunction

## One table variable as the numbers a learner takes, and the levels it was
## read through where it holds text.
function [v, lev, iscat, errmsg] = tfColumn (col, name)

  v = [];
  lev = [];
  iscat = false;
  errmsg = '';

  if (isnumeric (col) && isreal (col))
    v = double (col);
    return;
  endif
  if (islogical (col))
    v = double (col);
    iscat = true;
    return;
  endif
  if (isa (col, 'categorical'))
    lev = categories (col)(:)';
    v = double (col);
    iscat = true;
    return;
  endif
  if (ischar (col))
    col = cellstr (col);
  elseif (isa (col, 'string'))
    col = cellstr (col);
  endif
  if (iscellstr (col))
    lev = unique (col(:)')';
    lev = lev(:)';
    v = zeros (numel (col), 1);
    for k = 1:numel (lev)
      v(strcmp (col(:), lev{k})) = k;
    endfor
    iscat = true;
    return;
  endif
  errmsg = sprintf (strcat ("the table variable '%s' must hold numbers,", ...
                            " levels or text."), name);

endfunction

## The response and the predictors a model formula names, read by
## parseWilkinsonFormula rather than by anything of our own.  A learner
## takes main effects only, so a term over more than one variable is a
## product and a variable carrying a caret is a power, and R2024a refuses
## both.  The parser is asked for the order the formula wrote the names in,
## which is the order a learner reports them in, rather than the sorted
## order it gives by default.
function [respname, terms, errmsg] = tfFormula (spec, names)

  respname = '';
  terms = {};
  errmsg = '';
  try
    F = parseWilkinsonFormula (spec, 'matrix', 'stable');
  catch
    errmsg = sprintf ("'%s' is not a model formula.", spec);
    return;
  end_try_catch
  if (! (isfield (F, 'ResponseIdx') && isscalar (F.ResponseIdx)))
    errmsg = sprintf ("'%s' names no response.", spec);
    return;
  endif
  respname = F.VariableNames{F.ResponseIdx};

  ## A term over more than one variable is a product
  if (any (sum (F.Terms, 2) > 1))
    errmsg = strcat ("a model formula holds main effects only, so no", ...
                     " products, powers or wildcards.");
    return;
  endif
  used = F.VariableNames(any (F.Terms == 1, 1));
  if (any (cellfun (@(v) any (v == '^'), used)))
    errmsg = strcat ("a model formula holds main effects only, so no", ...
                     " products, powers or wildcards.");
    return;
  endif

  for nm = [{respname}, used]
    if (! any (strcmp (names, nm{1})))
      errmsg = sprintf (strcat ("the model formula names '%s', which the", ...
                                " table does not hold."), nm{1});
      return;
    endif
  endfor
  if (any (strcmp (used, respname)))
    errmsg = strcat ("a model formula cannot name the response among its", ...
                     " predictors.");
    return;
  endif

  terms = used;

endfunction

## Take one name-value pair out of the argument list, where it is there.
function [val, args] = tfTakeArg (args, name)

  val = [];
  for k = 1:2:numel (args) - 1
    if (ischar (args{k}) && strcmpi (args{k}, name))
      val = args{k+1};
      args(k:k+1) = [];
      return;
    endif
  endfor

endfunction

## The predictors named, in the order the call named them.
function [sub, errmsg] = tfNameSubset (pn, pnames)

  sub = {};
  errmsg = '';
  if (ischar (pn) && isrow (pn))
    pn = {pn};
  elseif (isa (pn, 'string'))
    pn = cellstr (pn);
  endif
  if (! iscellstr (pn))
    errmsg = "'PredictorNames' must name variables of the table.";
    return;
  endif
  for k = 1:numel (pn)
    if (! any (strcmp (pnames, pn{k})))
      errmsg = sprintf (strcat ("'PredictorNames' names no variable", ...
                                " '%s' of the table."), pn{k});
      return;
    endif
  endfor
  sub = pn(:)';

endfunction

## The categorical predictors a call named, as indices into the predictors.
function [idx, errmsg] = tfCatSpec (cp, pnames)

  idx = [];
  errmsg = '';
  M = numel (pnames);
  if (ischar (cp) && isrow (cp) && strcmpi (cp, 'all'))
    idx = 1:M;
    return;
  endif
  if (islogical (cp))
    if (numel (cp) != M)
      errmsg = strcat ("'CategoricalPredictors' as a logical vector must", ...
                       " have one value per predictor.");
      return;
    endif
    idx = find (cp(:)');
    return;
  endif
  if (isnumeric (cp))
    if (! (isreal (cp) && all (cp == fix (cp)) && all (cp >= 1)
           && all (cp <= M)))
      errmsg = strcat ("'CategoricalPredictors' must index the", ...
                       " predictors of the table.");
      return;
    endif
    idx = unique (double (cp(:)'));
    return;
  endif
  if (ischar (cp))
    cp = cellstr (cp);
  elseif (isa (cp, 'string'))
    cp = cellstr (cp);
  endif
  if (! iscellstr (cp))
    errmsg = strcat ("'CategoricalPredictors' must name or index the", ...
                     " predictors of the table.");
    return;
  endif
  for k = 1:numel (cp)
    j = find (strcmp (pnames, cp{k}), 1);
    if (isempty (j))
      errmsg = sprintf (strcat ("'CategoricalPredictors' names no", ...
                                " predictor '%s' of the table."), cp{k});
      return;
    endif
    idx(end+1) = j;
  endfor
  idx = unique (idx);

endfunction

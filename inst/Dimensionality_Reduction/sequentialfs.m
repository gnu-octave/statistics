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
## @deftypefn  {statistics} {@var{tf} =} sequentialfs (@var{fun}, @var{X}, @var{y})
## @deftypefnx {statistics} {@var{tf} =} sequentialfs (@var{fun}, @var{X1}, @dots{}, @var{XN})
## @deftypefnx {statistics} {@var{tf} =} sequentialfs (@dots{}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{tf}, @var{history}] =} sequentialfs (@dots{})
##
## Sequential feature selection.
##
## @code{@var{tf} = sequentialfs (@var{fun}, @var{X}, @var{y})} selects the
## columns of @var{X} that best predict @var{y}, adding one column at a time,
## each time the one that lowers the cross-validated criterion computed by
## @var{fun} the most, until no column lowers it enough.  @var{tf} is a
## logical row vector, true for each selected column of @var{X}.  By default
## the criterion is estimated by 10-fold cross-validation, on a random
## partition of the rows that is not stratified.
##
## @var{fun} is a function handle called as
## @code{@var{criterion} = @var{fun} (@var{XTRAIN}, @var{ytrain}, @var{XTEST},
## @var{ytest})}, where @var{XTRAIN} and @var{ytrain} are the training rows of
## @var{X}, restricted to the candidate columns, and of @var{y}, and
## @var{XTEST} and @var{ytest} the test rows.  It returns a real scalar, which
## is typically a loss summed over the test rows: the values @var{fun}
## returns for every test set are summed and divided by the total number of
## test rows, so a loss averaged over the test rows must be multiplied by
## their number.  Lower values are better.
##
## @code{@var{tf} = sequentialfs (@var{fun}, @var{X1}, @dots{}, @var{XN})}
## partitions the rows of every data argument alike, all of which must have
## the same number of rows, at least two, and selects columns of @var{X1}
## only.  @var{fun} is then called with the training parts of all of them
## followed by the test parts, @code{@var{fun} (@var{X1TRAIN}, @dots{},
## @var{XNTRAIN}, @var{X1TEST}, @dots{}, @var{XNTEST})}.  The data arguments
## are the leading arguments up to the first option name.
##
## The following @qcode{Name-Value} pairs are supported:
##
## @multitable @columnfractions 0.2 0.8
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'CV'} @tab How the criterion is validated: a positive
## integer @math{k}, greater than 1 and not greater than the number of rows,
## for @math{k}-fold cross-validation on a random partition that is not
## stratified; a @code{cvpartition} object; @qcode{'resubstitution'}, calling
## @code{@var{fun} (@var{X}, @var{y}, @var{X}, @var{y})} on all rows and
## dividing by their number; or @qcode{'none'}, calling @code{@var{fun}
## (@var{X}, @var{y})} on all rows and taking its value as it is.  The
## default is 10.
##
## @item @qcode{'MCReps'} @tab The number of Monte Carlo repetitions of the
## cross-validation, a positive integer, each on a new random partition, the
## criterion being summed over all of them before it is divided by the total
## number of test rows.  The default is 1.  It is set to 1, with a warning,
## when @qcode{'CV'} is @qcode{'none'} or @qcode{'resubstitution'}, or a
## @code{cvpartition} object for leave-one-out, resubstitution or a custom
## partition, none of which can be drawn anew.
##
## @item @qcode{'Direction'} @tab @qcode{'forward'}, the default, starts from
## no column and adds one at each step; @qcode{'backward'} starts from all
## columns and removes one at each step.
##
## @item @qcode{'KeepIn'} @tab Columns always selected, given as a logical
## vector with one element per column of @var{X} or as column indices.  The
## default is none.
##
## @item @qcode{'KeepOut'} @tab Columns never selected, given as
## @qcode{'KeepIn'}.  A column cannot be in both.
##
## @item @qcode{'NFeatures'} @tab The number of columns to select, a positive
## integer.  The search then goes on until that many columns are selected,
## whatever the criterion and @qcode{MaxIter}.  The default selects by the
## criterion alone.
##
## @item @qcode{'NullModel'} @tab @qcode{true} to let the selection hold no
## column, in which case @var{fun} must accept data arguments with no column;
## the default is @qcode{false}.  A forward search then starts by evaluating
## the empty selection, and a backward search may remove the last column.
##
## @item @qcode{'Options'} @tab A structure, such as @code{statset
## ('sequentialfs')} returns, with the fields @qcode{Display}, @qcode{'off'}
## (the default), @qcode{'final'} or @qcode{'iter'}; @qcode{MaxIter}, the
## largest number of steps, @code{Inf} by default; @qcode{TolFun}, the
## tolerance on the criterion, 1e-6 for a forward search and 0 for a
## backward search by default; and @qcode{TolTypeFun}, @qcode{'rel'} (the
## default) or @qcode{'abs'}.  Missing or empty fields take their defaults.
## @end multitable
##
## At each step every remaining column is tried and the one giving the
## lowest criterion is taken, the first of them in the order of the columns
## on a tie; a criterion of @qcode{NaN} is never taken over a number.  The
## search stops when no column is left, after @qcode{MaxIter} steps, or when
## the new criterion @math{c} is greater than @math{c_0 - TolFun (|c_0| +
## sqrt (eps))} for a relative tolerance, or than @math{c_0 - TolFun} for an
## absolute one, @math{c_0} being the criterion before the step, so a step
## that leaves the criterion unchanged is taken when @qcode{TolFun} is 0.  A
## criterion of @qcode{NaN} on either side never stops the search.  The
## first column of a forward search is always added, unless
## @qcode{'NullModel'} or @qcode{'KeepIn'} gives a starting selection to
## compare it with.
##
## @code{[@var{tf}, @var{history}] = sequentialfs (@dots{})} also returns a
## structure with the fields @qcode{In}, a logical matrix whose row @math{i}
## marks the columns selected after step @math{i}, and @qcode{Crit}, a row
## vector of the criterion after each step.  The starting selection is the
## first step when it is evaluated: always for a backward search, and for a
## forward search when @qcode{'KeepIn'} or @qcode{'NullModel'} is given.  Both
## fields are empty when nothing is evaluated.
##
## The option names are not case sensitive but must be given in full, where
## MATLAB accepts any shortened name.  @code{sequentialfs} refuses an
## @qcode{'NFeatures'} of 0, a @qcode{'CV'} larger than the number of rows,
## a @qcode{'NullModel'} other than true or false, and a criterion that is
## not a real number, all of which MATLAB accepts.  Parallel computation is
## not supported, and the @qcode{UseParallel}, @qcode{UseSubstreams} and
## @qcode{Streams} options are ignored.
##
## @seealso{relieff, cvpartition, crossval, statset, stepwiselm}
## @end deftypefn

function [tf, history] = sequentialfs (fun, varargin)

  if (nargin < 2)
    print_usage ();
  endif
  if (! is_function_handle (fun))
    error ("sequentialfs: FUN must be a function handle.");
  endif

  ## Data arguments, up to the first option name
  nd = numel (varargin);
  for i = 1:numel (varargin)
    if ((ischar (varargin{i}) && rows (varargin{i}) <= 1)
        || (isstring (varargin{i}) && isscalar (varargin{i})))
      nd = i - 1;
      break;
    endif
  endfor
  if (nd == 0)
    print_usage ();
  endif
  data = varargin(1:nd);
  n = size (data{1}, 1);
  if (any (cellfun (@(d) size (d, 1), data) != n))
    error (strcat ("sequentialfs: the data arguments must have the same", ...
                   " number of rows."));
  endif
  if (n < 2)
    error ("sequentialfs: the data arguments must have at least two rows.");
  endif
  p = size (data{1}, 2);

  ## Parse optional Name-Value paired arguments
  optNames = {'CV', 'MCReps', 'Direction', 'KeepIn', 'KeepOut', ...
              'NFeatures', 'NullModel', 'Options'};
  dfValues = {10, 1, 'forward', [], [], [], false, []};
  [cv, mcreps, direction, keepin, keepout, nfeat, nullmodel, options, ...
   args] = parsePairedArguments (optNames, dfValues, varargin(nd+1:end));
  if (! isempty (args))
    if (mod (numel (args), 2) == 0 && ischar (args{1}))
      error ("sequentialfs: unknown parameter name '%s'.", args{1});
    endif
    error ("sequentialfs: optional arguments must be Name-Value pairs.");
  endif

  if (! ((ischar (direction) && isrow (direction))
         || (isstring (direction) && isscalar (direction))))
    error ("sequentialfs: 'Direction' must be a character vector.");
  endif
  direction = lower (char (direction));
  if (! any (strcmp (direction, {'forward', 'backward'})))
    error ("sequentialfs: unknown direction '%s'.", direction);
  endif
  forward = strcmp (direction, 'forward');

  if (! (isnumeric (mcreps) && isreal (mcreps) && isscalar (mcreps)
         && isfinite (mcreps) && mcreps >= 1 && mcreps == fix (mcreps)))
    error ("sequentialfs: 'MCReps' must be a positive integer.");
  endif

  if (isstring (cv) && isscalar (cv))
    cv = char (cv);
  endif
  if (isa (cv, 'cvpartition'))
    if (cv.NumObservations != n)
      error (strcat ("sequentialfs: the cvpartition object must have as", ...
                     " many observations as the data arguments have rows."));
    endif
    if (mcreps > 1 && (cv.IsCustom
                       || any (strcmpi (cv.Type, {'leaveout', ...
                                                  'resubstitution'}))))
      warning (strcat ("sequentialfs: 'MCReps' is set to 1 for a", ...
                       " leave-one-out, resubstitution or custom", ...
                       " partition."));
      mcreps = 1;
    endif
    validation = 'partition';
  elseif (ischar (cv) && any (strcmpi (cv, {'resubstitution', 'none'})))
    validation = lower (cv);
    if (mcreps > 1)
      warning (strcat ("sequentialfs: 'MCReps' is set to 1 when 'CV' is", ...
                       " 'none' or 'resubstitution'."));
      mcreps = 1;
    endif
  elseif (isnumeric (cv) && isreal (cv) && isscalar (cv) && isfinite (cv)
          && cv > 1 && cv == fix (cv))
    if (cv > n)
      error (strcat ("sequentialfs: 'CV' must not be greater than the", ...
                     " number of rows."));
    endif
    validation = 'partition';
    cv = cvpartition (n, 'KFold', cv);
  else
    error (strcat ("sequentialfs: 'CV' must be a cvpartition object, an", ...
                   " integer greater than 1, 'resubstitution' or 'none'."));
  endif

  [keepin, errmsg] = column_mask (keepin, p, 'KeepIn');
  if (! isempty (errmsg))
    error (errmsg);
  endif
  [keepout, errmsg] = column_mask (keepout, p, 'KeepOut');
  if (! isempty (errmsg))
    error (errmsg);
  endif
  if (any (keepin & keepout))
    error ("sequentialfs: a column cannot be in both 'KeepIn' and 'KeepOut'.");
  endif

  if (! isempty (nfeat))
    if (! (isnumeric (nfeat) && isreal (nfeat) && isscalar (nfeat)
           && nfeat >= 1 && nfeat <= p && nfeat == fix (nfeat)))
      error (strcat ("sequentialfs: 'NFeatures' must be a positive integer", ...
                     " not greater than the number of columns of X."));
    endif
    if (nfeat < sum (keepin))
      error (strcat ("sequentialfs: 'NFeatures' must not be smaller than", ...
                     " the number of columns in 'KeepIn'."));
    endif
    if (nfeat > p - sum (keepout))
      error (strcat ("sequentialfs: 'NFeatures' must not be greater than", ...
                     " the number of columns not in 'KeepOut'."));
    endif
  endif

  if (! ((islogical (nullmodel) || isnumeric (nullmodel))
         && isscalar (nullmodel) && (nullmodel == 0 || nullmodel == 1)))
    error ("sequentialfs: 'NullModel' must be true or false.");
  endif
  nullmodel = logical (nullmodel);

  ## Search options, missing or empty fields taking their defaults
  display = 'off';
  maxiter = Inf;
  tolfun = 0;
  if (forward)
    tolfun = 1e-6;
  endif
  reltol = true;
  if (! isempty (options))
    if (! (isstruct (options) && isscalar (options)))
      error ("sequentialfs: 'Options' must be a structure.");
    endif
    if (isfield (options, 'Display') && ! isempty (options.Display))
      display = options.Display;
      if (isstring (display) && isscalar (display))
        display = char (display);
      endif
      if (! (ischar (display)
             && any (strcmpi (display, {'off', 'final', 'iter'}))))
        error ("sequentialfs: unknown 'Display' value.");
      endif
      display = lower (display);
    endif
    if (isfield (options, 'MaxIter') && ! isempty (options.MaxIter))
      maxiter = options.MaxIter;
      if (! (isnumeric (maxiter) && isreal (maxiter) && isscalar (maxiter)
             && maxiter >= 0 && maxiter == fix (maxiter)))
        error (strcat ("sequentialfs: 'MaxIter' must be a non-negative", ...
                       " integer or Inf."));
      endif
    endif
    if (isfield (options, 'TolFun') && ! isempty (options.TolFun))
      tolfun = options.TolFun;
      if (! (isnumeric (tolfun) && isreal (tolfun) && isscalar (tolfun)
             && tolfun >= 0 && isfinite (tolfun)))
        error ("sequentialfs: 'TolFun' must be a non-negative scalar.");
      endif
    endif
    if (isfield (options, 'TolTypeFun') && ! isempty (options.TolTypeFun))
      toltype = options.TolTypeFun;
      if (isstring (toltype) && isscalar (toltype))
        toltype = char (toltype);
      endif
      if (! (ischar (toltype) && any (strcmpi (toltype, {'abs', 'rel'}))))
        error ("sequentialfs: 'TolTypeFun' must be 'abs' or 'rel'.");
      endif
      reltol = strcmpi (toltype, 'rel');
    endif
  endif

  ## Training and test rows of every repetition of the partition
  trainrows = testrows = {};
  if (strcmp (validation, 'partition'))
    for r = 1:mcreps
      if (r > 1)
        cv = repartition (cv);
      endif
      for i = 1:cv.NumTestSets
        trainrows{end+1} = training (cv, i);
        testrows{end+1} = test (cv, i);
      endfor
    endfor
  endif

  ## Starting selection
  if (forward)
    in = keepin;
  else
    in = ! keepout;
  endif
  H = false (0, p);
  C = zeros (1, 0);
  verbose = ! strcmp (display, 'off');
  if (verbose)
    if (forward)
      printf ("Start forward sequential feature selection:\n");
      printf ("Initial columns included:  %s\n", column_list (in));
      printf ("Columns that can not be included:  %s\n", ...
              column_list (keepout));
    else
      printf ("Start backward sequential feature selection:\n");
      printf ("Initial columns included:  %s\n", column_list (in));
      printf ("Columns that must be included:  %s\n", column_list (keepin));
    endif
  endif
  crit = NaN;
  if ((forward && (any (in) || nullmodel)) || (! forward && any (in)))
    crit = criterion (fun, data, in, validation, trainrows, testrows, n);
    H(end+1,:) = in;
    C(end+1) = crit;
    if (verbose)
      printf ("Step %d, used initial columns, criterion value %g\n", ...
              numel (C), crit);
    endif
  endif

  ## Add or remove one column at each step
  steps = 0;
  while (! isempty (nfeat) || steps < maxiter)
    if (forward)
      cand = find (! in & ! keepout);
      if (! isempty (nfeat) && sum (in) >= nfeat)
        break;
      endif
    else
      cand = find (in & ! keepin);
      if ((! nullmodel && sum (in) <= 1)
          || (! isempty (nfeat) && sum (in) <= nfeat))
        break;
      endif
    endif
    if (isempty (cand))
      break;
    endif
    c = zeros (1, numel (cand));
    for j = 1:numel (cand)
      trial = in;
      trial(cand(j)) = forward;
      c(j) = criterion (fun, data, trial, validation, trainrows, testrows,
                        n);
    endfor
    [newcrit, best] = min (c);
    if (isempty (nfeat))
      if (reltol)
        limit = crit - tolfun * (abs (crit) + sqrt (eps));
      else
        limit = crit - tolfun;
      endif
      if (newcrit > limit)
        break;
      endif
    endif
    in(cand(best)) = forward;
    crit = newcrit;
    steps += 1;
    H(end+1,:) = in;
    C(end+1) = crit;
    if (strcmp (display, 'iter'))
      if (forward)
        printf ("Step %d, added column %d, criterion value %g\n", ...
                numel (C), cand(best), crit);
      else
        printf ("Step %d, removed column %d, criterion value %g\n", ...
                numel (C), cand(best), crit);
      endif
    endif
  endwhile

  if (verbose)
    printf ("Final columns included:  %s\n", column_list (in));
  endif
  tf = in;
  if (isempty (C))
    history = struct ('In', false (0, 0), 'Crit', []);
  else
    history = struct ('In', H, 'Crit', C);
  endif

endfunction

## A column selection as a logical row, from a logical vector or indices
function [mask, errmsg] = column_mask (v, p, name)
  mask = false (1, p);
  errmsg = '';
  if (isempty (v))
    return;
  endif
  if (islogical (v) && isvector (v) && numel (v) == p)
    mask = v(:)';
  elseif (isnumeric (v) && isreal (v) && isvector (v)
          && all (v == fix (v)) && all (v >= 1) && all (v <= p))
    mask(v) = true;
  else
    errmsg = sprintf (strcat ("sequentialfs: '%s' must be a logical", ...
                              " vector with one element per column of X,", ...
                              " or a vector of column indices."), name);
  endif
endfunction

## The selected columns as text, as the display shows them
function txt = column_list (mask)
  if (! any (mask))
    txt = 'none';
  elseif (all (mask))
    txt = 'all';
  else
    txt = sprintf ('%d ', find (mask));
  endif
endfunction

## The criterion of the columns in MASK, validated as asked
function crit = criterion (fun, data, mask, validation, trainrows, testrows, n)
  sel = data;
  sel{1} = data{1}(:, mask);
  switch (validation)
    case 'none'
      crit = checked (fun (sel{:}));
    case 'resubstitution'
      crit = checked (fun (sel{:}, sel{:})) / n;
    otherwise
      total = 0;
      ntest = 0;
      for i = 1:numel (trainrows)
        tr = cellfun (@(d) d(trainrows{i},:), sel, 'UniformOutput', false);
        te = cellfun (@(d) d(testrows{i},:), sel, 'UniformOutput', false);
        total += checked (fun (tr{:}, te{:}));
        ntest += sum (testrows{i});
      endfor
      crit = total / ntest;
  endswitch
endfunction

## The value FUN returned, which must be a real scalar
function v = checked (v)
  if (! ((isnumeric (v) || islogical (v)) && isreal (v) && isscalar (v)))
    error ("sequentialfs: FUN must return a real scalar.");
  endif
  v = double (v);
endfunction

%!demo
%! ## Select the ingredients of cement that predict the heat it gives off,
%! ## judging each selection by its sum of squared errors under
%! ## leave-one-out cross-validation
%! load hald
%! sse = @(XT, yT, Xt, yt) sum ((yt - [ones(rows (Xt), 1), Xt] * ...
%!                               ([ones(rows (XT), 1), XT] \ yT)) .^ 2);
%! cv = cvpartition (13, 'LeaveOut');
%! opts = statset ('Display', 'iter');
%! tf = sequentialfs (sse, ingredients, heat, 'CV', cv, 'Options', opts)

%!demo
%! ## Remove ingredients one at a time from the full model instead
%! load hald
%! sse = @(XT, yT, Xt, yt) sum ((yt - [ones(rows (Xt), 1), Xt] * ...
%!                               ([ones(rows (XT), 1), XT] \ yT)) .^ 2);
%! cv = cvpartition (13, 'LeaveOut');
%! [tf, history] = sequentialfs (sse, ingredients, heat, 'CV', cv, ...
%!                               'Direction', 'backward')

## Test output
## Values checked against MATLAB R2024a.  Column j of X holds the value j,
## so the criterion F looks up the selected columns in the table T.
%!shared X, y, m
%! X = repmat (1:3, 10, 1);
%! y = (1:10)';
%! m = @(Z) 1 + sum (2 .^ (Z(1,:) - 1));
%!test
%! T = [100, 10, 20, 9, 30, 9.5, 25, 8];
%! [tf, h] = sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV', 'none');
%! assert_equal (tf, true (1, 3));
%! assert_equal (h.In, logical ([1, 0, 0; 1, 1, 0; 1, 1, 1]));
%! assert_equal (h.Crit, [10, 9, 8]);
%!test
%! T = [100, 10, 20, 9, 30, 9.5, 25, 8];
%! [tf, h] = sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV', 'none', ...
%!                         'Options', statset ('TolFun', 0.1));
%! assert_equal (tf, logical ([1, 0, 0]));
%! assert_equal (h.Crit, 10);
%!test
%! T = [100, 10, 50, 8.99999999];
%! [~, h] = sequentialfs (@(Z, w) T(m (Z)), X(:,1:2), y, 'CV', 'none', ...
%!                        'Options', statset ('TolFun', 0.1));
%! assert_equal (h.Crit, [10, 8.99999999]);
%!test
%! T = [100, 10, 50, 8.999999999];
%! [~, h] = sequentialfs (@(Z, w) T(m (Z)), X(:,1:2), y, 'CV', 'none', ...
%!                        'Options', statset ('TolFun', 0.1));
%! assert_equal (h.Crit, 10);
%!test
%! T = [100, 0, 50, -0.05];
%! [~, h] = sequentialfs (@(Z, w) T(m (Z)), X(:,1:2), y, 'CV', 'none', ...
%!                        'Options', statset ('TolFun', 0.1));
%! assert_equal (h.Crit, [0, -0.05]);
%!test
%! T = [100, 10, 20, 9.05, 30, 9.5, 25, 9.05];
%! [tf, h] = sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV', 'none', ...
%!                         'Options', struct ('TolFun', 0));
%! assert_equal (tf, true (1, 3));
%! assert_equal (h.Crit, [10, 9.05, 9.05]);
%!test
%! T = [100, 10, 20, 9.05, 30, 9.5, 25, 8];
%! [~, h] = sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV', 'none', ...
%!                        'Options', statset ('TolFun', 0.95, ...
%!                                            'TolTypeFun', 'abs'));
%! assert_equal (h.Crit, [10, 9.05, 8]);
%!test
%! T = [100, 10, 20, 9.05, 30, 9.5, 25, 8];
%! [~, h] = sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV', 'none', ...
%!                        'Options', statset ('TolFun', 0.96, ...
%!                                            'TolTypeFun', 'abs'));
%! assert_equal (h.Crit, 10);
%!test
%! T = [100, 10, 10, 9, 20, 9, 25, 8];
%! [~, h] = sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV', 'none');
%! assert_equal (h.In, logical ([1, 0, 0; 1, 1, 0; 1, 1, 1]));
%!test
%! T = [100, NaN, 20, 9, 30, 9.5, 25, 8];
%! [~, h] = sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV', 'none');
%! assert_equal (h.In, logical ([0, 1, 0; 1, 1, 0; 1, 1, 1]));
%! assert_equal (h.Crit, [20, 9, 8]);
%!test
%! T = [100, NaN, NaN, 9, NaN, 9.5, 25, 8];
%! [~, h] = sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV', 'none');
%! assert_equal (h.In, logical ([1, 0, 0; 1, 1, 0; 1, 1, 1]));
%! assert_equal (h.Crit, [NaN, 9, 8]);
%!test
%! T = [100, 10, 50, NaN];
%! [~, h] = sequentialfs (@(Z, w) T(m (Z)), X(:,1:2), y, 'CV', 'none');
%! assert_equal (h.Crit, [10, NaN]);
%!test
%! T = [5, 10, 20, 9, 30, 9.5, 25, 8];
%! [tf, h] = sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV', 'none', ...
%!                         'NullModel', true);
%! assert_equal (tf, false (1, 3));
%! assert_equal (h.In, false (1, 3));
%! assert_equal (h.Crit, 5);
%!test
%! T = [10.5, 10, 20, 9, 30, 9.5, 25, 8];
%! [~, h] = sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV', 'none', ...
%!                        'NullModel', true);
%! assert_equal (h.In, logical ([0, 0, 0; 1, 0, 0; 1, 1, 0; 1, 1, 1]));
%! assert_equal (h.Crit, [10.5, 10, 9, 8]);
%!test
%! T = [100, 10, 20, 9, 30, 9.5, 25, 8];
%! [~, h] = sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV', 'none', 'KeepIn', 3);
%! assert_equal (h.In, logical ([0, 0, 1; 1, 0, 1; 1, 1, 1]));
%! assert_equal (h.Crit, [30, 9.5, 8]);
%!test
%! T = [100, 10, 20, 9, 30, 9.5, 25, 8];
%! [tf, h] = sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV', 'none', ...
%!                         'KeepOut', [false, true, false]);
%! assert_equal (tf, logical ([1, 0, 1]));
%! assert_equal (h.Crit, [10, 9.5]);
%!test
%! T = [100, 10, 20, 9, 30, 9.5, 25, 8];
%! [tf, h] = sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV', 'none', ...
%!                         'KeepOut', [1, 2, 3]);
%! assert_equal (tf, false (1, 3));
%! assert_equal (size (h.In), [0, 0]);
%! assert_equal (size (h.Crit), [0, 0]);
%!test
%! T = [100, 10, 20, 9, 30, 9.5, 25, 8];
%! [~, h] = sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV', 'none', ...
%!                        'NFeatures', 3, 'Options', statset ('TolFun', 0.5));
%! assert_equal (h.Crit, [10, 9, 8]);
%!test
%! T = [100, 10, 20, 9, 30, 9.5, 25, 8];
%! [~, h] = sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV', 'none', ...
%!                        'NFeatures', 2, 'Options', statset ('MaxIter', 1));
%! assert_equal (h.Crit, [10, 9]);
%!test
%! T = [100, 10, 20, 9, 30, 9.5, 25, 8];
%! [~, h] = sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV', 'none', ...
%!                        'Options', statset ('MaxIter', 2));
%! assert_equal (h.Crit, [10, 9]);
%!test
%! T = [100, 10, 20, 9, 30, 9.5, 25, 8];
%! [~, h] = sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV', 'none', ...
%!                        'NullModel', true, 'Options', statset ('MaxIter', 1));
%! assert_equal (h.In, logical ([0, 0, 0; 1, 0, 0]));
%! assert_equal (h.Crit, [100, 10]);
%!test
%! T = [100, 9, 8, 10, 7, 12, 10, 10];
%! [tf, h] = sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV', 'none', ...
%!                         'Direction', 'backward');
%! assert_equal (tf, logical ([0, 0, 1]));
%! assert_equal (h.In, logical ([1, 1, 1; 0, 1, 1; 0, 0, 1]));
%! assert_equal (h.Crit, [10, 10, 7]);
%!test
%! T = [100, 50, 10.000000000000002, 10];
%! [~, h] = sequentialfs (@(Z, w) T(m (Z)), X(:,1:2), y, 'CV', 'none', ...
%!                        'Direction', 'backward');
%! assert_equal (h.Crit, 10);
%!test
%! T = [100, 50, 9, 10];
%! [~, h] = sequentialfs (@(Z, w) T(m (Z)), X(:,1:2), y, 'CV', 'none', ...
%!                        'Direction', 'backward', ...
%!                        'Options', statset ('TolFun', 0.1));
%! assert_equal (h.Crit, 10);
%!test
%! T = [100, 50, 8.99, 10];
%! [~, h] = sequentialfs (@(Z, w) T(m (Z)), X(:,1:2), y, 'CV', 'none', ...
%!                        'Direction', 'backward', ...
%!                        'Options', statset ('TolFun', 0.1));
%! assert_equal (h.Crit, [10, 8.99]);
%!test
%! T = [100, 50, 9, 10];
%! [~, h] = sequentialfs (@(Z, w) T(m (Z)), X(:,1:2), y, 'CV', 'none', ...
%!                        'Direction', 'backward', ...
%!                        'Options', statset ('TolFun', 1, ...
%!                                            'TolTypeFun', 'abs'));
%! assert_equal (h.Crit, [10, 9]);
%!test
%! T = [100, 50, NaN, 10];
%! [~, h] = sequentialfs (@(Z, w) T(m (Z)), X(:,1:2), y, 'CV', 'none', ...
%!                        'Direction', 'backward');
%! assert_equal (h.Crit, 10);
%!test
%! T = [6, 50, 7, 10];
%! [tf, h] = sequentialfs (@(Z, w) T(m (Z)), X(:,1:2), y, 'CV', 'none', ...
%!                         'Direction', 'backward', 'NullModel', true);
%! assert_equal (tf, false (1, 2));
%! assert_equal (h.In, logical ([1, 1; 0, 1; 0, 0]));
%! assert_equal (h.Crit, [10, 7, 6]);
%!test
%! T = [100, 9, 8, 10.5, 7, 12, 10, 10];
%! [~, h] = sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV', 'none', ...
%!                        'Direction', 'backward', 'KeepIn', 2);
%! assert_equal (h.In, logical ([1, 1, 1; 0, 1, 1; 0, 1, 0]));
%! assert_equal (h.Crit, [10, 10, 8]);
%!test
%! T = [100, 9, 8, 10.5, 7, 12, 10, 10];
%! [~, h] = sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV', 'none', ...
%!                        'Direction', 'backward', 'KeepOut', 3);
%! assert_equal (h.In, logical ([1, 1, 0; 0, 1, 0]));
%! assert_equal (h.Crit, [10.5, 8]);
%!test
%! T = [100, 10, 20, 9, 30, 9.5, 25, 8];
%! [~, h] = sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV', 'none', ...
%!                        'Direction', 'backward', 'KeepIn', 3, ...
%!                        'NFeatures', 1);
%! assert_equal (h.Crit, [8, 9.5, 30]);
%!test
%! T = [100, 10, 20, 9, 30, 9.5, 25, 8];
%! [tf, h] = sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV', 'none', ...
%!                         'Direction', 'backward', ...
%!                         'Options', struct ('MaxIter', 0));
%! assert_equal (tf, true (1, 3));
%! assert_equal (h.Crit, 8);
%!test
%! T = [100, 10, 20, 9, 30, 9.5, 25, 8];
%! [~, h] = sequentialfs (@(A, b, C, d) T(m (A)), X, y, 'CV', 5);
%! assert_equal (h.Crit, [5, 4.5, 4]);
%!test
%! T = [100, 10, 20, 9, 30, 9.5, 25, 8];
%! [~, h] = sequentialfs (@(A, b, C, d) T(m (A)), X, y, 'CV', 5, 'MCReps', 3);
%! assert_equal (h.Crit, [5, 4.5, 4]);
%!test
%! T = [100, 10, 20, 9, 30, 9.5, 25, 8];
%! [~, h] = sequentialfs (@(A, b, C, d) T(m (A)), X, y);
%! assert_equal (h.Crit, [10, 9, 8]);
%!test
%! T = [100, 10, 20, 9, 30, 9.5, 25, 8];
%! [~, h] = sequentialfs (@(A, b, C, d) T(m (A)), X, y, ...
%!                        'CV', 'resubstitution');
%! assert_equal (h.Crit, [1, 0.9, 0.8], eps);
%!test
%! T = [100, 10, 20, 9, 30, 9.5, 25, 8];
%! [~, h] = sequentialfs (@(A, b, C, d) T(m (A)), X, y, ...
%!                        'CV', cvpartition (10, 'Holdout', 3), 'MCReps', 4);
%! assert_equal (h.Crit, [10, 9, 8] / 3, eps);
%!test
%! T = [100, 10, 20, 9, 30, 9.5, 25, 8];
%! q = @(A, B, c, D, E, z) T(m (A)) + columns (B) / 1000 + rows (E) / 1e6;
%! [~, h] = sequentialfs (q, X, [y, y], y, 'CV', 'resubstitution');
%! assert_equal (h.Crit, [1.000201, 0.900201, 0.800201], 1e-15);
%!test
%! T = [100, 10, 20, 9, 30, 9.5, 25, 8];
%! [tf, h] = sequentialfs (@(A) T(m (A)), X, 'CV', 'none');
%! assert_equal (h.Crit, [10, 9, 8]);
%!test
%! T = [100, 9, 8, 10.5, 7, 12, 10, 10];
%! out = evalc (strcat ("sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV',", ...
%!                      " 'none', 'Options', statset ('Display', 'iter'));"));
%! assert_equal (out, ["Start forward sequential feature selection:\n", ...
%!                     "Initial columns included:  none\n", ...
%!                     "Columns that can not be included:  none\n", ...
%!                     "Step 1, added column 3, criterion value 7\n", ...
%!                     "Final columns included:  3 \n"]);
%!test
%! T = [100, 10, 20, 9, 30, 9.5, 25, 8];
%! out = evalc (strcat ("sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV',", ...
%!                      " 'none', 'KeepIn', [1, 3], 'Options',", ...
%!                      " statset ('Display', 'iter'));"));
%! assert_equal (out, ["Start forward sequential feature selection:\n", ...
%!                     "Initial columns included:  1 3 \n", ...
%!                     "Columns that can not be included:  none\n", ...
%!                     "Step 1, used initial columns, ", ...
%!                     "criterion value 9.5\n", ...
%!                     "Step 2, added column 2, criterion value 8\n", ...
%!                     "Final columns included:  all\n"]);
%!test
%! T = [100, 10, 20, 9, 30, 9.5, 25, 8];
%! out = evalc (strcat ("sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV',", ...
%!                      " 'none', 'KeepOut', [1, 3], 'Options',", ...
%!                      " statset ('Display', 'iter'));"));
%! assert_equal (out, ["Start forward sequential feature selection:\n", ...
%!                     "Initial columns included:  none\n", ...
%!                     "Columns that can not be included:  1 3 \n", ...
%!                     "Step 1, added column 2, criterion value 20\n", ...
%!                     "Final columns included:  2 \n"]);
%!test
%! T = [100, 10, 20, 9, 30, 9.5, 25, 8];
%! out = evalc (strcat ("sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV',", ...
%!                      " 'none', 'NullModel', true, 'Options',", ...
%!                      " statset ('Display', 'final'));"));
%! assert_equal (out, ["Start forward sequential feature selection:\n", ...
%!                     "Initial columns included:  none\n", ...
%!                     "Columns that can not be included:  none\n", ...
%!                     "Step 1, used initial columns, ", ...
%!                     "criterion value 100\n", ...
%!                     "Final columns included:  all\n"]);
%!test
%! T = [100, 9, 8, 10.5, 7, 12, 10, 10];
%! out = evalc (strcat ("sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV',", ...
%!                      " 'none', 'Direction', 'backward', 'Options',", ...
%!                      " statset ('Display', 'iter'));"));
%! assert_equal (out, ["Start backward sequential feature selection:\n", ...
%!                     "Initial columns included:  all\n", ...
%!                     "Columns that must be included:  none\n", ...
%!                     "Step 1, used initial columns, criterion value 10\n", ...
%!                     "Step 2, removed column 1, criterion value 10\n", ...
%!                     "Step 3, removed column 2, criterion value 7\n", ...
%!                     "Final columns included:  3 \n"]);
%!test
%! T = [100, 10, 20, 9, 30, 9.5, 25, 8];
%! out = evalc (strcat ("sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV',", ...
%!                      " 'none', 'Direction', 'backward', 'KeepOut', 2,", ...
%!                      " 'Options', statset ('Display', 'final'));"));
%! assert_equal (out, ["Start backward sequential feature selection:\n", ...
%!                     "Initial columns included:  1 3 \n", ...
%!                     "Columns that must be included:  none\n", ...
%!                     "Step 1, used initial columns, ", ...
%!                     "criterion value 9.5\n", ...
%!                     "Final columns included:  1 3 \n"]);
%!test
%! T = [100, 10, 20, 9, 30, 9.5, 25, 8];
%! out = evalc (strcat ("sequentialfs (@(Z, w) T(m (Z)), X, y, 'CV',", ...
%!                      " 'none', 'Direction', 'backward', 'KeepOut',", ...
%!                      " [1, 2, 3], 'Options',", ...
%!                      " statset ('Display', 'iter'));"));
%! assert_equal (out, ["Start backward sequential feature selection:\n", ...
%!                     "Initial columns included:  none\n", ...
%!                     "Columns that must be included:  none\n", ...
%!                     "Final columns included:  none\n"]);
%!shared sse, cv
%! sse = @(XT, yT, Xt, yt) sum ((yt - [ones(rows (Xt), 1), Xt] * ...
%!                               ([ones(rows (XT), 1), XT] \ yT)) .^ 2);
%! cv = cvpartition (13, 'LeaveOut');
%!test
%! load hald
%! [tf, h] = sequentialfs (sse, ingredients, heat, 'CV', cv);
%! assert_equal (tf, logical ([1, 1, 0, 1]));
%! assert_equal (h.In, logical ([0, 0, 0, 1; 1, 0, 0, 1; 1, 1, 0, 1]));
%! assert_equal (h.Crit, [91.86293870821909, 9.3249533115070999, ...
%!                        6.5654708620637523], 1e-10);
%!test
%! load hald
%! [tf, h] = sequentialfs (sse, ingredients, heat, 'CV', cv, ...
%!                         'Direction', 'backward');
%! assert_equal (tf, logical ([1, 1, 0, 1]));
%! assert_equal (h.Crit, [8.4881966849412365, 6.5654708620637523], 1e-10);
%!test
%! load hald
%! [~, h] = sequentialfs (sse, ingredients, heat, 'CV', cv, 'NullModel', true);
%! assert_equal (h.Crit, [245.17305555555561, 91.86293870821909, ...
%!                        9.3249533115070999, 6.5654708620637523], 1e-10);
%!test
%! load hald
%! [tf, h] = sequentialfs (sse, ingredients, heat, 'CV', 'resubstitution', ...
%!                         'Direction', 'backward', 'NullModel', true);
%! assert_equal (tf, true (1, 4));
%! assert_equal (h.Crit, 3.6818184115768386, 1e-12);
%!warning <sequentialfs: 'MCReps' is set to 1 when 'CV' is 'none' or 'resubstitution'.> ...
%! sequentialfs (@(A, b) 1, [1, 2; 3, 4], [1; 2], 'CV', 'none', 'MCReps', 3);
%!warning <sequentialfs: 'MCReps' is set to 1 for a leave-one-out, resubstitution or custom partition.> ...
%! sequentialfs (@(A, b, C, d) 1, [1, 2; 3, 4], [1; 2], ...
%!               'CV', cvpartition (2, 'LeaveOut'), 'MCReps', 3);

## Test input validation
%!shared X, y, f
%! X = repmat (1:3, 10, 1);
%! y = (1:10)';
%! f = @(Z, w) sum (Z(1,:));
%!error <Invalid call to sequentialfs> sequentialfs (f)
%!error <Invalid call to sequentialfs> sequentialfs (f, 'CV', 'none')
%!error <sequentialfs: FUN must be a function handle.> sequentialfs (3, X, y)
%!error <sequentialfs: the data arguments must have the same number of rows.> ...
%! sequentialfs (f, X, y(1:9), 'CV', 'none')
%!error <sequentialfs: the data arguments must have at least two rows.> ...
%! sequentialfs (f, X(1,:), y(1), 'CV', 'none')
%!error <sequentialfs: unknown parameter name 'Foo'.> ...
%! sequentialfs (f, X, y, 'CV', 'none', 'Foo', 1)
%!error <sequentialfs: unknown parameter name 'dir'.> ...
%! sequentialfs (f, X, y, 'CV', 'none', 'dir', 'backward')
%!error <sequentialfs: optional arguments must be Name-Value pairs.> ...
%! sequentialfs (f, X, y, 'CV')
%!error <sequentialfs: 'Direction' must be a character vector.> ...
%! sequentialfs (f, X, y, 'CV', 'none', 'Direction', 1)
%!error <sequentialfs: unknown direction 'foo'.> ...
%! sequentialfs (f, X, y, 'CV', 'none', 'Direction', 'foo')
%!error <sequentialfs: 'MCReps' must be a positive integer.> ...
%! sequentialfs (f, X, y, 'MCReps', 0)
%!error <sequentialfs: 'MCReps' must be a positive integer.> ...
%! sequentialfs (f, X, y, 'MCReps', 1.5)
%!error <sequentialfs: the cvpartition object must have as many observations as the data arguments have rows.> ...
%! sequentialfs (f, X, y, 'CV', cvpartition (9, 'LeaveOut'))
%!error <sequentialfs: 'CV' must not be greater than the number of rows.> ...
%! sequentialfs (f, X, y, 'CV', 11)
%!error <sequentialfs: 'CV' must be a cvpartition object, an integer greater than 1, 'resubstitution' or 'none'.> ...
%! sequentialfs (f, X, y, 'CV', 'foo')
%!error <sequentialfs: 'CV' must be a cvpartition object, an integer greater than 1, 'resubstitution' or 'none'.> ...
%! sequentialfs (f, X, y, 'CV', 1)
%!error <sequentialfs: 'CV' must be a cvpartition object, an integer greater than 1, 'resubstitution' or 'none'.> ...
%! sequentialfs (f, X, y, 'CV', 2.5)
%!error <sequentialfs: 'KeepIn' must be a logical vector with one element per column of X, or a vector of column indices.> ...
%! sequentialfs (f, X, y, 'CV', 'none', 'KeepIn', 4)
%!error <sequentialfs: 'KeepIn' must be a logical vector with one element per column of X, or a vector of column indices.> ...
%! sequentialfs (f, X, y, 'CV', 'none', 'KeepIn', [true, false])
%!error <sequentialfs: 'KeepIn' must be a logical vector with one element per column of X, or a vector of column indices.> ...
%! sequentialfs (f, X, y, 'CV', 'none', 'KeepIn', 1.5)
%!error <sequentialfs: 'KeepOut' must be a logical vector with one element per column of X, or a vector of column indices.> ...
%! sequentialfs (f, X, y, 'CV', 'none', 'KeepOut', 0)
%!error <sequentialfs: a column cannot be in both 'KeepIn' and 'KeepOut'.> ...
%! sequentialfs (f, X, y, 'CV', 'none', 'KeepIn', 1, 'KeepOut', 1)
%!error <sequentialfs: 'NFeatures' must be a positive integer not greater than the number of columns of X.> ...
%! sequentialfs (f, X, y, 'CV', 'none', 'NFeatures', 4)
%!error <sequentialfs: 'NFeatures' must be a positive integer not greater than the number of columns of X.> ...
%! sequentialfs (f, X, y, 'CV', 'none', 'NFeatures', 0)
%!error <sequentialfs: 'NFeatures' must be a positive integer not greater than the number of columns of X.> ...
%! sequentialfs (f, X, y, 'CV', 'none', 'NFeatures', 1.5)
%!error <sequentialfs: 'NFeatures' must not be smaller than the number of columns in 'KeepIn'.> ...
%! sequentialfs (f, X, y, 'CV', 'none', 'KeepIn', [1, 2], 'NFeatures', 1)
%!error <sequentialfs: 'NFeatures' must not be greater than the number of columns not in 'KeepOut'.> ...
%! sequentialfs (f, X, y, 'CV', 'none', 'KeepOut', [1, 2], 'NFeatures', 2)
%!error <sequentialfs: 'NullModel' must be true or false.> ...
%! sequentialfs (f, X, y, 'CV', 'none', 'NullModel', 2)
%!error <sequentialfs: 'NullModel' must be true or false.> ...
%! sequentialfs (f, X, y, 'CV', 'none', 'NullModel', 'on')
%!error <sequentialfs: 'Options' must be a structure.> ...
%! sequentialfs (f, X, y, 'CV', 'none', 'Options', 1)
%!error <sequentialfs: unknown 'Display' value.> ...
%! sequentialfs (f, X, y, 'CV', 'none', 'Options', struct ('Display', 'on'))
%!error <sequentialfs: 'MaxIter' must be a non-negative integer or Inf.> ...
%! sequentialfs (f, X, y, 'CV', 'none', 'Options', struct ('MaxIter', -1))
%!error <sequentialfs: 'TolFun' must be a non-negative scalar.> ...
%! sequentialfs (f, X, y, 'CV', 'none', 'Options', struct ('TolFun', -1))
%!error <sequentialfs: 'TolTypeFun' must be 'abs' or 'rel'.> ...
%! sequentialfs (f, X, y, 'CV', 'none', 'Options', struct ('TolTypeFun', 'foo'))
%!error <sequentialfs: FUN must return a real scalar.> ...
%! sequentialfs (@(Z, w) [1, 2], X, y, 'CV', 'none')
%!error <sequentialfs: FUN must return a real scalar.> ...
%! sequentialfs (@(Z, w) 'a', X, y, 'CV', 'none')
%!error <sequentialfs: FUN must return a real scalar.> ...
%! sequentialfs (@(Z, w) 1i, X, y, 'CV', 'none')
%!error <sequentialfs: FUN must return a real scalar.> ...
%! sequentialfs (@(Z, w) [], X, y, 'CV', 'none')

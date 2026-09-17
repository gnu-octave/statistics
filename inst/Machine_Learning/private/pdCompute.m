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
## @deftypefn {Private Function} {[@var{pd}, @var{x}, @var{y}] =} pdCompute (@var{fname}, @var{Mdl}, @var{Vars}, @dots{})
##
## Compute partial dependence, for whichever entry point was called.
##
## @code{partialDependence} and the method of the same name on
## @code{PredictiveModel} both come here, so that the two cannot drift apart.
## @var{fname} names the caller in every error message, and the rest is what
## that caller was given.  @code{plotPartialDependence} comes here too, for
## the values it draws.
##
## @end deftypefn

function [pd, x, y] = pdCompute (fname, Mdl, Vars, varargin)

  if (! (isa (Mdl, 'PredictiveModel') || is_function_handle (Mdl)))
    error (strcat ("%s: MDL must be a fitted model that", ...
                   " predicts, or a function handle."), fname);
  endif

  ## A classifier is told which classes to answer for, and the class names
  ## come before the data, as MATLAB orders them
  isclass = (! is_function_handle (Mdl)
             && any (strcmp (properties (Mdl), 'ClassNames'))
             && ! isempty (Mdl.ClassNames));
  args = varargin;
  Labels = [];
  if (isclass)
    if (isempty (args))
      error (strcat ("%s: LABELS is required for a", ...
                     " classification model."), fname);
    endif
    Labels = args{1};
    args(1) = [];
  endif
  Data = [];
  if (! isempty (args) && isnumeric (args{1}))
    Data = args{1};
    args(1) = [];
  elseif (! isempty (args)
          && ! (ischar (args{1}) || isa (args{1}, 'string')))
    ## Neither the data nor the name of an option: for a regression model the
    ## likeliest mistake is classes it has none of, and for a classifier,
    ## whose classes were taken above, it is the data
    if (isclass)
      error ("%s: DATA must be a real numeric matrix.", fname);
    else
      error (strcat ("%s: LABELS applies only to a classification", ...
                     " model."), fname);
    endif
  endif

  optNames = {'QueryPoints', 'NumObservationsToSample', ...
              'CategoricalPredictors', 'IncludeInteractions', ...
              'IncludeIntercept', 'OutputColumns', 'UseParallel', ...
              'PredictionForMissingValue'};
  dfValues = {[], [], [], [], true, 'all', [], []};
  [QP, NumObs, CatPred, Inter, Icept, OutCols, Par, Miss, rem] = ...
                              parsePairedArguments (optNames, dfValues, args);
  if (! isempty (rem))
    error (strcat ("%s: unknown optional argument or", ...
                   " misplaced value."), fname);
  endif
  if (! isempty (Par))
    error ("%s: 'UseParallel' is not implemented.", fname);
  endif
  if (! isempty (Miss))
    error (strcat ("%s: 'PredictionForMissingValue' is not", ...
                   " implemented."), fname);
  endif

  [predArgs, subIcept, OutCols, errmsg] = pdOptions (Mdl, Inter, Icept, ...
                                                     OutCols, CatPred);
  if (! isempty (errmsg))
    error ("%s: %s", fname, errmsg);
  endif

  [F, errmsg] = pdFrame (Mdl, Vars, Labels, Data, QP, NumObs, CatPred, ...
                         OutCols);
  if (! isempty (errmsg))
    error ("%s: %s", fname, errmsg);
  endif

  pd = pdValues (Mdl, F, predArgs, subIcept);
  x = F.QP{1};
  if (numel (F.Vars) == 2)
    y = F.QP{2};
  else
    y = [];
  endif

endfunction

## The options that belong to a kind of model rather than to every one.
function [predArgs, subIcept, OutCols, errmsg] = pdOptions (Mdl, Inter, ...
                                                            Icept, OutCols, ...
                                                            CatPred)

  predArgs = {};
  subIcept = false;
  errmsg = '';
  isfh = is_function_handle (Mdl);
  if (isfh)
    props = {};
  else
    props = properties (Mdl);
  endif
  isgam = any (strcmp (props, 'Interactions'));

  if (! isempty (Inter))
    if (! isgam)
      errmsg = strcat ("'IncludeInteractions' applies only to a", ...
                       " generalized additive model.");
      return;
    endif
    if (! (islogical (Inter) && isscalar (Inter)))
      errmsg = "'IncludeInteractions' must be true or false.";
      return;
    endif
    predArgs = {'IncludeInteractions', Inter};
  endif

  if (! (islogical (Icept) && isscalar (Icept)))
    errmsg = "'IncludeIntercept' must be true or false.";
    return;
  endif
  if (! Icept)
    if (! isgam)
      errmsg = strcat ("'IncludeIntercept' applies only to a generalized", ...
                       " additive model.");
      return;
    endif
    subIcept = true;
  endif

  if (! isfh && ! isempty (CatPred))
    errmsg = strcat ("'CategoricalPredictors' applies only to a function", ...
                     " handle, a model carrying its own.");
    return;
  endif

  if (ischar (OutCols) && strcmpi (OutCols, 'all'))
    OutCols = [];
  elseif (isnumeric (OutCols) && isreal (OutCols) && isvector (OutCols)
          && all (OutCols == fix (OutCols)) && all (OutCols >= 1))
    OutCols = double (OutCols(:)');
  else
    errmsg = strcat ("'OutputColumns' must be a vector of positive", ...
                     " integers or 'all'.");
    return;
  endif
  if (! isfh && ! isempty (OutCols))
    errmsg = "'OutputColumns' applies only to a function handle.";
  endif

endfunction

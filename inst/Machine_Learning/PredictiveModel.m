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

classdef (Abstract) PredictiveModel
  ## -*- texinfo -*-
  ## @deftp {statistics} PredictiveModel
  ##
  ## A fitted model that answers for new observations.
  ##
  ## @code{PredictiveModel} is the superclass of every classification and
  ## regression model in this package that carries a @code{predict} method,
  ## and of nothing else.  The models fitted by @code{fitlm} and its
  ## relatives derive from it as the learners do.  The cross-validated
  ## classes answer through @code{kfoldPredict} over their folds rather than
  ## through @code{predict}, and so do not derive from it.
  ##
  ## The class is abstract and cannot be instantiated.  It is where
  ## behaviour shared by all those models is written once rather than once
  ## per class, and @code{isa (@var{obj}, @qcode{'PredictiveModel'})} is how
  ## to ask whether an object is a model that can be told to predict.  The
  ## only data it holds is the levels a predictor read from a table was
  ## coded through, which every model needs to read a table again at
  ## prediction and which MATLAB carries by keeping the table itself.
  ##
  ## @end deftp

  ## This file must live in a directory that sorts before every directory
  ## holding a subclass of it: Machine_Learning before Regression, today.
  ## Octave below 11.2.0 rebuilds every doc-cache when the package is
  ## installed, one directory at a time and with only that directory on the
  ## path, so a subclass in another directory can resolve this class only
  ## where the walk has already met it.  genpath sorts, and a class stays
  ## resolved once it has been, so the earlier directory settles it.  The
  ## last test in this file fails if that ever stops being true.
  ##
  ## The shared behaviour of the models goes here.  The block stays while it
  ## is empty: a class body holding no block at all gets no help text, help
  ## synthesising a default constructor instead of reading the class block
  ## above.
  properties (GetAccess = public, SetAccess = protected, Hidden)
    ## One cell per predictor, holding the levels a predictor read from a
    ## table was coded through and empty for one read from numbers.  It is
    ## hidden because MATLAB's models carry no such property, keeping the
    ## table itself instead.
    PredictorLevels = {};
  endproperties

  methods (Access = protected)

    ## Resolve table input into the predictors and the response, keeping
    ## the levels the coding used.  A constructor calls this before it
    ## validates anything, and a matrix passes through untouched.
    function [this, X, Y, args] = resolveTable (this, caller, X, Y, args)

      if (! istable (X))
        return;
      endif
      [X, Y, args, lev, errmsg] = tableFrame (X, Y, args);
      if (! isempty (errmsg))
        error ("%s: %s", caller, errmsg);
      endif
      this.PredictorLevels = lev;

    endfunction

    ## Map a table onto the predictors this model was fitted on.  The work
    ## is a private function rather than this method alone, because an
    ## explainer holds a model without being one and needs the same mapping.
    function Z = tableColumns (this, caller, T)

      Z = tableToMatrix (this.PredictorNames, this.PredictorLevels, T, caller);

    endfunction

    ## Resolve table input into the predictors and the response that loss,
    ## margin and edge take.  METH is the method's own name, GIVEN says
    ## whether the call gave a third argument at all, and ARGS are the ones
    ## after it.  A matrix passes through untouched.
    function [X, Y, args] = tableResponse (this, meth, X, Y, args, given)

      if (! istable (X))
        return;
      endif
      caller = sprintf ('%s.%s', class (this), meth);

      ## The argument after the table is the response where an odd number
      ## of arguments follows the table, and a name-value name where an
      ## even number does, pairs being even.  Then the response is the
      ## variable the model was fitted on.
      if (mod (given + numel (args), 2) == 0)
        if (given)
          args = [{Y}, args];
        endif
        Y = this.ResponseName;
      endif

      ## A name is looked up among the table's variables; anything else is
      ## the response itself, given beside the table
      if ((ischar (Y) && isrow (Y)) || (isa (Y, 'string') && isscalar (Y)))
        name = char (Y);
        if (! any (strcmp (X.Properties.VariableNames, name)))
          error ("%s: the table holds no variable '%s'.", caller, name);
        endif
        Y = X.(name);
      endif

      X = tableColumns (this, caller, X);

    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {PredictiveModel} {@var{pd} =} partialDependence (@var{obj}, @var{Vars})
    ## @deftypefnx {PredictiveModel} {@var{pd} =} partialDependence (@var{obj}, @var{Vars}, @var{Labels})
    ## @deftypefnx {PredictiveModel} {@var{pd} =} partialDependence (@dots{}, @var{Data})
    ## @deftypefnx {PredictiveModel} {@var{pd} =} partialDependence (@dots{}, @var{name}, @var{value})
    ## @deftypefnx {PredictiveModel} {[@var{pd}, @var{x}, @var{y}] =} partialDependence (@dots{})
    ##
    ## Compute partial dependence.
    ##
    ## @code{partialDependence (@var{obj}, @var{Vars})} is the same call as
    ## @code{partialDependence} the function, written on the model instead of
    ## before it, and takes and returns exactly what it does.  See
    ## @code{partialDependence} for the arguments and for what @var{pd},
    ## @var{x} and @var{y} hold.
    ##
    ## @seealso{partialDependence, plotPartialDependence}
    ## @end deftypefn
    function varargout = partialDependence (this, varargin)

      [varargout{1:max(nargout, 1)}] = pdCompute ('partialDependence', ...
                                                  this, varargin{:});

    endfunction

  endmethods

endclassdef


%!test  # the class is abstract and cannot be instantiated
%! fail ('PredictiveModel ()', 'abstract');

%!test  # a model that predicts derives from it, a cross-validated one does not
%! load fisheriris
%! Mdl = fitctree (meas, species);
%! assert_equal (isa (Mdl, 'PredictiveModel'), true);
%! assert_equal (isa (compact (Mdl), 'PredictiveModel'), true);
%! assert_equal (isa (crossval (Mdl, 'KFold', 2), 'PredictiveModel'), false);

%!test  # a model fitted by fitlm derives from it too, across directories
%! X = [1, 2; 2, 3; 3, 4; 1, 5; 2, 6; 3, 7];
%! y = [2.5; 3.1; 4.8; 2.2; 3.9; 5.1];
%! assert_equal (isa (fitlm (X, y), 'PredictiveModel'), true);

%!test  # no subclass may sit in a directory sorting before this class's own
%! base = fileparts (which ('PredictiveModel'));
%! root = fileparts (base);
%! here = base(numel (root) + 2:end);
%! d = dir (root);
%! subs = {d([d.isdir]).name};
%! subs = subs(! ismember (subs, {'.', '..', here}));
%! for k = 1:numel (subs)
%!   f = dir (fullfile (root, subs{k}, '*.m'));
%!   for j = 1:numel (f)
%!     src = fileread (fullfile (root, subs{k}, f(j).name));
%!     pat = '^classdef[^\n]*<[^\n]*PredictiveModel';
%!     if (! isempty (regexp (src, pat, 'once', 'lineanchors')))
%!       order = sort ({here, subs{k}});
%!       assert_equal (order{1}, here);
%!     endif
%!   endfor
%! endfor

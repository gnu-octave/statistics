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
## @deftp {statistics} LinearFormula
##
## Model formula of a linear or generalized linear regression.
##
## A @qcode{LinearFormula} object describes the terms of a fitted model: which
## variables the model draws on, how they combine into terms, and how the whole
## thing reads back as a formula.  It is the class of the @code{Formula}
## property of a @code{LinearModel} and of a @code{GeneralizedLinearModel}, and
## is normally obtained from a fitted model rather than built directly.
##
## The object is defined by its terms matrix and the names of the variables that
## matrix is written over; every other property is derived from those two.  Each
## row of @code{Terms} is one term of the model and each column is one variable,
## the entry giving the power that variable carries in that term.  An all-zero
## row is the intercept.  The response variable occupies a column of its own,
## which is always zero.
##
## Converting the object with @code{char} renders the whole formula, response
## included, as @qcode{"y ~ 1 + x1 + x2"}; the @code{LinearPredictor} property
## holds the right-hand side on its own.  For a generalized linear model the
## response carries its link function, as in @qcode{"logit(y) ~ 1 + x1"}.
##
## @end deftp

classdef LinearFormula

  properties (SetAccess = private)

    ## -*- texinfo -*-
    ## @deftp {LinearFormula} {property} ResponseName
    ##
    ## Name of the response variable
    ##
    ## A character vector naming the variable on the left-hand side of the
    ## formula.  This property is read-only.
    ##
    ## @end deftp
    ResponseName = '';

    ## -*- texinfo -*-
    ## @deftp {LinearFormula} {property} VariableNames
    ##
    ## Names of all variables available to the model
    ##
    ## A cell array of character vectors naming every variable the model was
    ## given, whether or not it is used, in the order the data lists them.  The
    ## response is included.  This property is read-only.
    ##
    ## @end deftp
    VariableNames = {};

    ## -*- texinfo -*-
    ## @deftp {LinearFormula} {property} PredictorNames
    ##
    ## Names of the variables the model actually uses
    ##
    ## A cell array of character vectors holding those elements of
    ## @code{VariableNames} that appear in at least one term.  This property is
    ## read-only.
    ##
    ## @end deftp
    PredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {LinearFormula} {property} TermNames
    ##
    ## Name of each term of the model
    ##
    ## A column cell array of character vectors, one per row of @code{Terms},
    ## naming the term over the model's @emph{variables}: the intercept is
    ## @qcode{"(Intercept)"}, an interaction joins its factors with a colon, and
    ## a power is written with a caret.  A categorical variable contributes one
    ## term under its own name however many indicator columns it expands to, so
    ## these are not the coefficient names.  This property is read-only.
    ##
    ## @end deftp
    TermNames = {};

    ## -*- texinfo -*-
    ## @deftp {LinearFormula} {property} Terms
    ##
    ## Terms matrix of the model
    ##
    ## A numeric matrix with one row per term and one column per variable of
    ## @code{VariableNames}, each entry giving the power that variable carries
    ## in that term.  An all-zero row is the intercept.  This property is
    ## read-only.
    ##
    ## @end deftp
    Terms = [];

    ## -*- texinfo -*-
    ## @deftp {LinearFormula} {property} InModel
    ##
    ## Which variables take part in the model
    ##
    ## A logical row vector with one element per variable of
    ## @code{VariableNames}, true where that variable appears in at least one
    ## term.  The response is false.  This property is read-only.
    ##
    ## @end deftp
    InModel = logical ([]);

    ## -*- texinfo -*-
    ## @deftp {LinearFormula} {property} HasIntercept
    ##
    ## Whether the model carries an intercept
    ##
    ## A logical scalar, true when @code{Terms} holds an all-zero row.  This
    ## property is read-only.
    ##
    ## @end deftp
    HasIntercept = false;

    ## -*- texinfo -*-
    ## @deftp {LinearFormula} {property} LinearPredictor
    ##
    ## Right-hand side of the formula
    ##
    ## A character vector rendering the model's terms, such as
    ## @qcode{"1 + x1 + x2"}.  A pair of variables appearing both on their own
    ## and as an interaction is written as a product, so that
    ## @qcode{"x1 + x2 + x1:x2"} reads @qcode{"x1*x2"}.  This property is
    ## read-only.
    ##
    ## @end deftp
    LinearPredictor = '';

    ## -*- texinfo -*-
    ## @deftp {LinearFormula} {property} Link
    ##
    ## Link function applied to the response
    ##
    ## The link of a generalized linear model, given as its name, its numeric
    ## exponent, or a structure of function handles.  It is @qcode{"identity"}
    ## for a linear model.  This property is read-only.
    ##
    ## @end deftp
    Link = 'identity';

    ## -*- texinfo -*-
    ## @deftp {LinearFormula} {property} ModelFun
    ##
    ## Function computing the linear predictor
    ##
    ## A function handle taking the coefficient vector and the design matrix.
    ## This property is read-only.
    ##
    ## @end deftp
    ModelFun = @(b, X) X * b;

    ## -*- texinfo -*-
    ## @deftp {LinearFormula} {property} FunctionCalls
    ##
    ## Functions called from within the formula
    ##
    ## A cell array of character vectors, empty unless the formula applies a
    ## function to a variable.  This property is read-only.
    ##
    ## @end deftp
    FunctionCalls = cell (1, 0);

    ## -*- texinfo -*-
    ## @deftp {LinearFormula} {property} NTerms
    ##
    ## Number of terms in the model
    ##
    ## A non-negative integer, the number of rows of @code{Terms}.  This
    ## property is read-only.
    ##
    ## @end deftp
    NTerms = 0;

    ## -*- texinfo -*-
    ## @deftp {LinearFormula} {property} NVars
    ##
    ## Number of variables available to the model
    ##
    ## A non-negative integer, the number of elements of
    ## @code{VariableNames}.  This property is read-only.
    ##
    ## @end deftp
    NVars = 0;

    ## -*- texinfo -*-
    ## @deftp {LinearFormula} {property} NPredictors
    ##
    ## Number of variables the model uses
    ##
    ## A non-negative integer, the number of true elements of @code{InModel}.
    ## This property is read-only.
    ##
    ## @end deftp
    NPredictors = 0;

  endproperties

  methods (Hidden)

    function disp (this)
      fprintf ("%s\n", char (this));
    endfunction

    function display (this)
      name = inputname (1);
      if (isempty (name))
        name = 'ans';
      endif
      fprintf ("%s = %s\n", name, char (this));
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {LinearFormula} {@var{obj} =} LinearFormula ()
    ## @deftypefnx {LinearFormula} {@var{obj} =} LinearFormula (@var{terms}, @var{varnames})
    ## @deftypefnx {LinearFormula} {@var{obj} =} LinearFormula (@var{terms}, @var{varnames}, @var{name}, @var{value}, @dots{})
    ##
    ## Create a model formula from a terms matrix.
    ##
    ## @code{@var{obj} = LinearFormula ()} returns an empty formula.
    ##
    ## @code{@var{obj} = LinearFormula (@var{terms}, @var{varnames})} builds a
    ## formula whose terms matrix is @var{terms} and whose variables are named
    ## by the cell array of character vectors @var{varnames}.  @var{terms} must
    ## have one column per element of @var{varnames}; each row is one term and
    ## each entry the power its variable carries in that term.  An all-zero row
    ## is the intercept.
    ##
    ## The remaining properties are derived from these two arguments, except
    ## those given as @var{name}-@var{value} pairs:
    ##
    ## @multitable @columnfractions 0.25 0.75
    ## @headitem @var{name} @tab @var{value}
    ##
    ## @item @qcode{"ResponseName"} @tab A character vector naming the response
    ## variable.  It must be one of @var{varnames}.
    ##
    ## @item @qcode{"Link"} @tab The link function applied to the response, as a
    ## name, a numeric exponent, or a structure of function handles.  It
    ## defaults to @qcode{"identity"}.
    ##
    ## @item @qcode{"ModelFun"} @tab A function handle computing the linear
    ## predictor from the coefficients and the design matrix.
    ##
    ## @item @qcode{"FunctionCalls"} @tab A cell array of character vectors
    ## naming functions the formula applies to its variables.
    ## @end multitable
    ##
    ## @end deftypefn

    function this = LinearFormula (varargin)

      if (nargin == 0)
        return;
      endif

      if (nargin < 2)
        error ("LinearFormula: TERMS and VARNAMES are both required.");
      endif

      terms    = varargin{1};
      varnames = varargin{2};

      if (! isnumeric (terms) || ! isreal (terms) || ndims (terms) > 2)
        error ("LinearFormula: TERMS must be a real matrix.");
      endif
      if (! iscellstr (varnames))
        error (strcat ("LinearFormula: VARNAMES must be a cell array of", ...
                       " character vectors."));
      endif
      if (! isempty (terms) && columns (terms) != numel (varnames))
        error (strcat ("LinearFormula: TERMS must have one column per", ...
                       " element of VARNAMES."));
      endif

      if (mod (nargin - 2, 2) != 0)
        error ("LinearFormula: optional arguments must be name-value pairs.");
      endif

      respname = '';
      for i = 3:2:nargin
        name = varargin{i};
        if (! (ischar (name) && isrow (name)))
          error ("LinearFormula: option names must be character vectors.");
        endif
        switch (lower (name))
          case 'responsename'
            respname = varargin{i+1};
            if (! (ischar (respname) && (isrow (respname) ...
                                         || isempty (respname))))
              error (strcat ("LinearFormula: 'ResponseName' must be a", ...
                             " character vector."));
            endif
          case 'link'
            this.Link = varargin{i+1};
          case 'modelfun'
            this.ModelFun = varargin{i+1};
          case 'functioncalls'
            this.FunctionCalls = varargin{i+1};
          otherwise
            error ("LinearFormula: unknown option '%s'.", name);
        endswitch
      endfor

      varnames = varnames(:)';
      if (! isempty (respname) && ! any (strcmp (varnames, respname)))
        error ("LinearFormula: '%s' is not one of VARNAMES.", respname);
      endif

      this.ResponseName  = respname;
      this.VariableNames = varnames;
      this.Terms         = terms;
      this.NVars         = numel (varnames);
      this.NTerms        = rows (terms);

      if (isempty (terms))
        this.InModel = false (1, this.NVars);
      else
        this.InModel = any (terms != 0, 1);
      endif
      this.PredictorNames = varnames(this.InModel);
      this.NPredictors    = sum (this.InModel);
      this.HasIntercept   = ! isempty (terms) && any (all (terms == 0, 2));
      this.TermNames      = term_names (terms, varnames);
      this.LinearPredictor = predictor_string (terms, varnames, ...
                                               this.TermNames);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {LinearFormula} {@var{str} =} char (@var{obj})
    ##
    ## Render a model formula as a character vector.
    ##
    ## @code{@var{str} = char (@var{obj})} returns the whole formula, response
    ## included, as in @qcode{"y ~ 1 + x1 + x2"}.  The response carries the
    ## link function of a generalized linear model, as in
    ## @qcode{"logit(y) ~ 1 + x1"}.
    ##
    ## @end deftypefn

    function str = char (this)

      if (isempty (this.LinearPredictor) && isempty (this.ResponseName))
        str = '';
        return;
      endif
      str = sprintf ("%s ~ %s", linked_response (this.ResponseName, ...
                                                 this.Link), ...
                     this.LinearPredictor);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {LinearFormula} {@var{str} =} string (@var{obj})
    ##
    ## Render a model formula as a string scalar.
    ##
    ## @code{@var{str} = string (@var{obj})} is the @code{string} counterpart of
    ## @code{char}.
    ##
    ## @end deftypefn

    function str = string (this)
      str = string (char (this));
    endfunction

  endmethods

endclassdef

## Name each row of a terms matrix over the model's variables.
function names = term_names (terms, varnames)

  n_terms = rows (terms);
  names   = cell (n_terms, 1);
  for t = 1:n_terms
    idx = find (terms(t, :) != 0);
    if (isempty (idx))
      names{t} = '(Intercept)';
    else
      parts = cell (1, numel (idx));
      for k = 1:numel (idx)
        parts{k} = factor_string (varnames{idx(k)}, terms(t, idx(k)));
      endfor
      names{t} = strjoin (parts, ':');
    endif
  endfor

endfunction

## One factor of a term, carrying its power when that is not 1.
function str = factor_string (name, power)
  if (power == 1)
    str = name;
  else
    str = sprintf ("%s^%d", name, power);
  endif
endfunction

## Render the right-hand side of a formula.  Two variables that appear both on
## their own and as their interaction are written as a product -- 'x1 + x2 +
## x1:x2' becomes 'x1*x2' -- and the main effects it absorbs are then left out.
## Only pairs collapse: a three-way interaction is written out as it stands even
## when every one of its sub-terms is present, and a power never collapses.
function str = predictor_string (terms, varnames, names)

  n_terms  = rows (terms);
  absorbed = false (n_terms, 1);
  rendered = names;

  for t = 1:n_terms
    idx = find (terms(t, :) != 0);
    if (numel (idx) != 2 || any (terms(t, idx) != 1))
      continue;
    endif
    main = zeros (1, 2);
    for k = 1:2
      row = zeros (1, columns (terms));
      row(idx(k)) = 1;
      hit = find (all (terms == row, 2), 1);
      if (isempty (hit))
        main = [];
        break;
      endif
      main(k) = hit;
    endfor
    if (isempty (main))
      continue;
    endif
    rendered{t}    = strjoin (varnames(idx), '*');
    absorbed(main) = true;
  endfor

  parts = {};
  for t = 1:n_terms
    if (absorbed(t))
      continue;
    elseif (all (terms(t, :) == 0))
      parts{end+1} = '1';
    else
      parts{end+1} = rendered{t};
    endif
  endfor

  if (isempty (parts))
    str = '';
  else
    str = strjoin (parts, ' + ');
  endif

endfunction

## The response as the formula shows it, wrapped in its link function.
function str = linked_response (respname, link)

  if (isempty (link) || (ischar (link) && strcmpi (link, 'identity')))
    str = respname;
  elseif (ischar (link))
    str = sprintf ("%s(%s)", link, respname);
  elseif (isnumeric (link) && isscalar (link))
    ## The identity link written as an exponent leaves the response alone.
    if (link == 1)
      str = respname;
    else
      str = sprintf ("power(%s,%g)", respname, link);
    endif
  else
    str = sprintf ("link(%s)", respname);
  endif

endfunction

%!demo
%! ## The Formula property of a fitted model is a LinearFormula object.
%! x1 = [1 2 3 4 5 6 7 8]';
%! x2 = [2 1 4 3 6 5 8 7]';
%! y  = [3.1 4.2 5.3 6.4 7.5 8.6 9.7 10.8]';
%! mdl = fitlm ([x1, x2], y, 'interactions');
%! f = mdl.Formula
%! ## Rendering it as text gives back the whole formula, response included,
%! ## while the LinearPredictor property holds the right-hand side alone.
%! char (f)
%! f.LinearPredictor
%! ## A pair of variables present both on their own and as their interaction
%! ## is written as a product.
%! f.TermNames
%! f.Terms

%!demo
%! ## A formula can also be built directly from a terms matrix.  Each column is
%! ## a variable and each entry the power it carries in that term; the all-zero
%! ## row is the intercept and the response keeps a column of its own.
%! terms = [0 0 0; 1 0 0; 0 1 0; 1 1 0];
%! f = LinearFormula (terms, {'dose', 'age', 'score'}, ...
%!                    'ResponseName', 'score')
%! f.PredictorNames
%! f.InModel

## Properties derived from the terms matrix
%!test
%! f = LinearFormula ([0 0 0; 1 0 0; 0 1 0], {'u', 'g', 'resp'}, ...
%!                    'ResponseName', 'resp');
%! assert_equal (f.ResponseName, 'resp');
%! assert_equal (f.VariableNames, {'u', 'g', 'resp'});
%! assert_equal (f.PredictorNames, {'u', 'g'});
%! assert_equal (f.TermNames, {'(Intercept)'; 'u'; 'g'});
%! assert_equal (f.Terms, [0 0 0; 1 0 0; 0 1 0]);
%! assert_equal (f.InModel, [true, true, false]);
%! assert_equal (f.HasIntercept, true);
%! assert_equal (f.LinearPredictor, '1 + u + g');
%! assert_equal (f.NTerms, 3);
%! assert_equal (f.NVars, 3);
%! assert_equal (f.NPredictors, 2);
%! assert_equal (f.Link, 'identity');
%! assert_equal (f.FunctionCalls, cell (1, 0));

%!test  # the response takes no part in the model
%! f = LinearFormula ([0 0 0 0; 1 0 0 0; 0 0 1 0], {'u', 'v', 'w', 'resp'}, ...
%!                    'ResponseName', 'resp');
%! assert_equal (f.InModel, [true, false, true, false]);
%! assert_equal (f.PredictorNames, {'u', 'w'});
%! assert_equal (f.NPredictors, 2);
%! assert_equal (f.NVars, 4);

%!test  # an empty formula is renderable
%! f = LinearFormula ();
%! assert_equal (char (f), '');
%! assert_equal (f.NTerms, 0);
%! assert_equal (f.HasIntercept, false);

## Rendering
%!test  # char gives the whole formula, LinearPredictor the right-hand side
%! f = LinearFormula ([0 0 0; 1 0 0; 0 1 0], {'u', 'g', 'resp'}, ...
%!                    'ResponseName', 'resp');
%! assert_equal (char (f), 'resp ~ 1 + u + g');
%! assert_equal (f.LinearPredictor, '1 + u + g');

%!test  # a pair present on its own and as an interaction reads as a product
%! f = LinearFormula ([0 0 0; 1 0 0; 0 1 0; 1 1 0], {'u', 'v', 'resp'}, ...
%!                    'ResponseName', 'resp');
%! assert_equal (char (f), 'resp ~ 1 + u*v');

%!test  # a product needs both main effects; without one the term stands alone
%! f = LinearFormula ([0 0 0; 1 0 0; 1 1 0], {'u', 'v', 'resp'}, ...
%!                    'ResponseName', 'resp');
%! assert_equal (char (f), 'resp ~ 1 + u + u:v');

%!test  # only pairs collapse, never a three-way interaction
%! terms = [0 0 0 0; 1 0 0 0; 0 1 0 0; 0 0 1 0; ...
%!          1 1 0 0; 1 0 1 0; 0 1 1 0; 1 1 1 0];
%! f = LinearFormula (terms, {'u', 'v', 'w', 'resp'}, ...
%!                    'ResponseName', 'resp');
%! assert_equal (char (f), 'resp ~ 1 + u*v + u*w + v*w + u:v:w');

%!test  # a power is never written as a product
%! f = LinearFormula ([0 0 0; 1 0 0; 0 1 0; 1 1 0; 2 0 0; 0 2 0], ...
%!                    {'u', 'v', 'resp'}, 'ResponseName', 'resp');
%! assert_equal (char (f), 'resp ~ 1 + u*v + u^2 + v^2');
%! assert_equal (f.TermNames, ...
%!               {'(Intercept)'; 'u'; 'v'; 'u:v'; 'u^2'; 'v^2'});

%!test  # dropping the intercept drops the leading 1
%! f = LinearFormula ([1 0 0; 0 1 0], {'u', 'g', 'resp'}, ...
%!                    'ResponseName', 'resp');
%! assert_equal (char (f), 'resp ~ u + g');
%! assert_equal (f.HasIntercept, false);

%!test  # an intercept-only model
%! f = LinearFormula ([0 0 0], {'u', 'g', 'resp'}, 'ResponseName', 'resp');
%! assert_equal (char (f), 'resp ~ 1');
%! assert_equal (f.TermNames, {'(Intercept)'});
%! assert_equal (f.NPredictors, 0);
%! assert_equal (f.PredictorNames, cell (1, 0));

%!test  # string renders the same text as char
%! f = LinearFormula ([0 0; 1 0], {'u', 'resp'}, 'ResponseName', 'resp');
%! s = string (f);
%! assert_equal (class (s), 'string');
%! assert_equal (char (s), char (f));

%!test  # disp writes the rendered formula
%! f = LinearFormula ([0 0; 1 0], {'u', 'resp'}, 'ResponseName', 'resp');
%! s = evalc ('disp (f)');
%! assert_equal (strtrim (s), 'resp ~ 1 + u');

## The link wraps the response
%!test  # a named link
%! f = LinearFormula ([0 0 0; 1 0 0; 0 1 0], {'u', 'g', 'yb'}, ...
%!                    'ResponseName', 'yb', 'Link', 'logit');
%! assert_equal (char (f), 'logit(yb) ~ 1 + u + g');
%! assert_equal (f.Link, 'logit');

%!test  # the identity link leaves the response alone, named or as an exponent
%! f = LinearFormula ([0 0; 1 0], {'u', 'y'}, 'ResponseName', 'y', ...
%!                    'Link', 'identity');
%! assert_equal (char (f), 'y ~ 1 + u');
%! f = LinearFormula ([0 0; 1 0], {'u', 'y'}, 'ResponseName', 'y', 'Link', 1);
%! assert_equal (char (f), 'y ~ 1 + u');

%!test  # a numeric link is shown as a power
%! f = LinearFormula ([0 0; 1 0], {'u', 'y'}, 'ResponseName', 'y', 'Link', -2);
%! assert_equal (char (f), 'power(y,-2) ~ 1 + u');

%!test  # a link given as a structure of handles is named generically
%! lk = struct ('Link', @(x) x, 'Derivative', @(x) 1, 'Inverse', @(x) x);
%! f = LinearFormula ([0 0; 1 0], {'u', 'y'}, 'ResponseName', 'y', ...
%!                    'Link', lk);
%! assert_equal (char (f), 'link(y) ~ 1 + u');

## Remaining options
%!test  # ModelFun and FunctionCalls are carried through
%! fun = @(b, X) X * b;
%! f = LinearFormula ([0 0; 1 0], {'u', 'y'}, 'ResponseName', 'y', ...
%!                    'ModelFun', fun, 'FunctionCalls', {'log'});
%! assert_equal (func2str (f.ModelFun), func2str (fun));
%! assert_equal (f.FunctionCalls, {'log'});

%!test  # option names are matched without regard to case
%! f = LinearFormula ([0 0; 1 0], {'u', 'y'}, 'responsename', 'y', ...
%!                    'LINK', 'log');
%! assert_equal (char (f), 'log(y) ~ 1 + u');

## Input validation
%!error <LinearFormula: TERMS and VARNAMES are both required.> ...
%! LinearFormula ([0 0])
%!error <LinearFormula: TERMS must be a real matrix.> ...
%! LinearFormula ('abc', {'u', 'y'})
%!error <LinearFormula: TERMS must be a real matrix.> ...
%! LinearFormula ([1+2i, 0], {'u', 'y'})
%!error <LinearFormula: VARNAMES must be a cell array of character vectors.> ...
%! LinearFormula ([0 0], 'uy')
%!error <LinearFormula: VARNAMES must be a cell array of character vectors.> ...
%! LinearFormula ([0 0], {1, 2})
%!error <LinearFormula: TERMS must have one column per element of VARNAMES.> ...
%! LinearFormula ([0 0 0], {'u', 'y'})
%!error <LinearFormula: optional arguments must be name-value pairs.> ...
%! LinearFormula ([0 0], {'u', 'y'}, 'Link')
%!error <LinearFormula: option names must be character vectors.> ...
%! LinearFormula ([0 0], {'u', 'y'}, 5, 'log')
%!error <LinearFormula: unknown option 'bogus'.> ...
%! LinearFormula ([0 0], {'u', 'y'}, 'bogus', 1)
%!error <LinearFormula: 'ResponseName' must be a character vector.> ...
%! LinearFormula ([0 0], {'u', 'y'}, 'ResponseName', 5)
%!error <LinearFormula: 'z' is not one of VARNAMES.> ...
%! LinearFormula ([0 0], {'u', 'y'}, 'ResponseName', 'z')

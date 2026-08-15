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
## @deftypefn  {statistics} {@var{options} =} statset ()
## @deftypefnx {statistics} {@var{options} =} statset (@var{funcname})
## @deftypefnx {statistics} {@var{options} =} statset (@var{name}, @var{value}, @dots{})
## @deftypefnx {statistics} {@var{options} =} statset (@var{oldopts}, @var{name}, @var{value}, @dots{})
## @deftypefnx {statistics} {@var{options} =} statset (@var{oldopts}, @var{newopts})
## @deftypefnx {statistics} {} statset ()
##
## Create or modify an options structure for iterative statistics algorithms.
##
## @code{@var{options} = statset ()} returns a structure carrying every
## recognized option name, each set to an empty value.  An empty option means
## @qcode{"use the calling function's own default"}, so an all-empty structure
## changes nothing wherever it is passed.
##
## @code{@var{options} = statset (@var{funcname})} returns the options that
## @var{funcname} uses by default, with the remaining fields left empty.
## @var{funcname} must name a function of this package that documents an
## @qcode{"Options"} argument; see the list below.  Unlike the name/value
## forms, this form takes no further arguments.
##
## @code{@var{options} = statset (@var{name}, @var{value}, @dots{})} returns an
## otherwise empty structure with the named options set.  Option names are
## matched case-insensitively and must be given in full.
##
## @code{@var{options} = statset (@var{oldopts}, @var{name}, @var{value},
## @dots{})} copies @var{oldopts} and applies the given name/value pairs to the
## copy.  @var{oldopts} is left unchanged.
##
## @code{@var{options} = statset (@var{oldopts}, @var{newopts})} merges two
## structures: every @emph{non-empty} field of @var{newopts} overrides its
## counterpart in @var{oldopts}, while an empty field of @var{newopts} leaves
## the @var{oldopts} value in place.  Fields that are not recognized option
## names are ignored in both structures.
##
## @code{statset ()} called with no output argument displays the recognized
## option names together with their permitted values, marking each default
## in braces.
##
## The recognized options are:
##
## @multitable @columnfractions 0.25 0.75
## @headitem Option @tab Description
## @item @qcode{"Display"} @tab Level of reporting: @qcode{"off"},
## @qcode{"final"}, or @qcode{"iter"}.
## @item @qcode{"MaxFunEvals"} @tab Maximum number of objective function
## evaluations, a positive scalar.
## @item @qcode{"MaxIter"} @tab Maximum number of iterations, a positive scalar.
## @item @qcode{"TolBnd"} @tab Positive scalar tolerance on parameter bounds.
## @item @qcode{"TolFun"} @tab Positive scalar tolerance on the objective
## function value.
## @item @qcode{"TolTypeFun"} @tab Whether @qcode{"TolFun"} is absolute,
## @qcode{"abs"}, or relative, @qcode{"rel"}.
## @item @qcode{"TolX"} @tab Positive scalar tolerance on the parameters.
## @item @qcode{"TolTypeX"} @tab Whether @qcode{"TolX"} is absolute,
## @qcode{"abs"}, or relative, @qcode{"rel"}.
## @item @qcode{"GradObj"} @tab Whether the objective function returns a
## gradient, @qcode{"off"} or @qcode{"on"}.
## @item @qcode{"Jacobian"} @tab Whether the model function returns a Jacobian,
## @qcode{"off"} or @qcode{"on"}.
## @item @qcode{"DerivStep"} @tab Relative step size for finite-difference
## derivatives, a positive scalar or vector.
## @item @qcode{"FunValCheck"} @tab Whether to check the objective function for
## invalid values, @qcode{"off"} or @qcode{"on"}.
## @item @qcode{"Robust"} @tab Whether to invoke a robust fit, @qcode{"off"} or
## @qcode{"on"}.  Superseded by @qcode{"RobustWgtFun"}.
## @item @qcode{"RobustWgtFun"} @tab Weight function for robust fitting: one of
## @qcode{"andrews"}, @qcode{"bisquare"}, @qcode{"cauchy"}, @qcode{"fair"},
## @qcode{"huber"}, @qcode{"logistic"}, @qcode{"talwar"}, @qcode{"welsch"}, a
## function handle, or empty for a non-robust fit.
## @item @qcode{"WgtFun"} @tab Weight function used with @qcode{"Robust"}.
## Superseded by @qcode{"RobustWgtFun"}.
## @item @qcode{"Tune"} @tab Positive tuning constant for the robust weight
## function.  Set automatically for a named weight function; required for a
## function handle.
## @item @qcode{"UseParallel"} @tab Logical flag requesting parallel
## computation.
## @item @qcode{"UseSubstreams"} @tab Logical flag requesting reproducible
## random substreams.
## @item @qcode{"Streams"} @tab A random stream or a cell array of them.
## @item @qcode{"OutputFcn"} @tab A function handle, or a cell array of them,
## called after each iteration.
## @end multitable
##
## @var{funcname} may name any of the following functions, each of which
## documents an @qcode{"Options"} argument: @code{copulafit}, @code{coxphfit},
## @code{crossval}, @code{evfit}, @code{factoran}, @code{fitglm},
## @code{fitglme}, @code{fitlme},
## @code{fitlmematrix}, @code{fitnlm}, @code{gamfit}, @code{gevfit},
## @code{glmfit}, @code{gmdistribution}, @code{gpfit}, @code{kmeans},
## @code{kmedoids}, @code{lasso}, @code{lassoglm}, @code{lognfit},
## @code{mdscale}, @code{mlecov}, @code{mlecustom}, @code{mvncdf},
## @code{mvtcdf}, @code{nbinfit}, @code{nlinfit}, @code{nnmf}, @code{normfit},
## @code{pca}, @code{plsregress}, @code{ppca}, @code{rocmetrics}, @code{tsne},
## @code{wblfit}, @code{GeneralizedLinearMixedModel}, and
## @code{LinearMixedModel}.
##
## Any function accepting an @qcode{"Options"} argument also accepts a plain
## structure carrying only the fields it needs, so @code{statset} is a
## convenience rather than a requirement.
##
## MATLAB's @code{statset} additionally accepts the names of functions this
## package does not provide.  Those names are rejected here rather than
## answered, since returning options for an absent function would assert a
## capability that does not exist.
##
## @seealso{statget, nlinfit, fitnlm, nnmf, mdscale, ppca, tsne, kmedoids}
## @end deftypefn

function options = statset (varargin)

  ## The recognized option names, in MATLAB's own field order.  'Robust' and
  ## 'WgtFun' are present in MATLAB's structure but absent from its published
  ## table; they are kept because nlinfit and fitnlm consume them.
  names = {'Display', 'MaxFunEvals', 'MaxIter', 'TolBnd', 'TolFun', ...
           'TolTypeFun', 'TolX', 'TolTypeX', 'GradObj', 'Jacobian', ...
           'DerivStep', 'FunValCheck', 'Robust', 'RobustWgtFun', 'WgtFun', ...
           'Tune', 'UseParallel', 'UseSubstreams', 'Streams', 'OutputFcn'};

  ## No output and no input: report the option names and their values.
  if (nargin == 0 && nargout == 0)
    display_options ();
    return;
  endif

  ## Start from an all-empty structure; 'Streams' is an empty cell.
  options = cell2struct (repmat ({[]}, numel (names), 1), names, 1);
  options.Streams = {};

  if (nargin == 0)
    return;
  endif

  ## A single character or string argument names a function.
  if (nargin == 1 && (ischar (varargin{1}) || isstring_scalar (varargin{1})))
    options = func_defaults (options, char (varargin{1}));
    return;
  endif

  ## An optional leading structure supplies the starting values, and an
  ## optional second structure overrides them where it is not empty.
  args = varargin;
  if (isstruct (args{1}))
    options = merge_struct (options, args{1}, names);
    args(1) = [];
    if (numel (args) > 0 && isstruct (args{1}))
      options = merge_struct (options, args{1}, names);
      args(1) = [];
    endif
  elseif (! ischar (args{1}) && ! isstring_scalar (args{1}))
    error (strcat ("statset: first argument must be a function name,", ...
                   " an option name, or an options structure."));
  endif

  if (mod (numel (args), 2) != 0)
    error ("statset: arguments must occur in NAME/VALUE pairs.");
  endif

  ## Apply the name/value pairs, remembering whether 'Tune' was given here.
  tune_given = false;
  for i = 1:2:numel (args)
    name = args{i};
    if (! ischar (name) && ! isstring_scalar (name))
      error ("statset: option name must be a character vector or a string.");
    endif
    idx = find (strcmpi (char (name), names));
    if (isempty (idx))
      error ("statset: unrecognized option name '%s'.", char (name));
    endif
    field = names{idx};
    options.(field) = check_value (field, args{i+1});
    if (strcmp (field, 'Tune'))
      tune_given = true;
    endif
  endfor

  ## A named weight function carries a default tuning constant, unless the
  ## caller set 'Tune' explicitly in this same call.
  if (! tune_given && ! isempty (options.RobustWgtFun) ...
      && ischar (options.RobustWgtFun))
    options.Tune = default_tune (options.RobustWgtFun);
  endif

endfunction

## Copy the recognized, non-empty fields of S over those of OPTIONS.
function options = merge_struct (options, s, names)
  if (! isscalar (s))
    error ("statset: an options structure must be a scalar structure.");
  endif
  fn = fieldnames (s);
  for i = 1:numel (fn)
    idx = find (strcmpi (fn{i}, names));
    if (isempty (idx))
      continue;                       # unrecognized fields are ignored
    endif
    value = s.(fn{i});
    if (! isempty (value))
      options.(names{idx}) = value;
    endif
  endfor
endfunction

## Validate VALUE for the option FIELD and return it in canonical form.
function value = check_value (field, value)

  ## An empty value always means "unset", whatever the option.
  if (isempty (value) && ! iscell (value))
    value = [];
    return;
  endif

  switch (field)

    case 'Display'
      value = check_string (field, value, {'off', 'final', 'iter'});

    case {'GradObj', 'Jacobian', 'FunValCheck', 'Robust'}
      value = check_string (field, value, {'off', 'on'});

    case {'TolTypeFun', 'TolTypeX'}
      value = check_string (field, value, {'abs', 'rel'});

    case {'MaxFunEvals', 'MaxIter'}
      if (! (isnumeric (value) && isscalar (value) && isreal (value)
             && isfloat (value) && value > 0))
        error (strcat ("statset: option '%s' must be a real positive", ...
                       " scalar."), field);
      endif

    case {'TolBnd', 'TolFun', 'TolX', 'Tune'}
      if (! (isnumeric (value) && isscalar (value) && isreal (value)
             && isfloat (value) && value > 0))
        error (strcat ("statset: option '%s' must be a real positive", ...
                       " scalar."), field);
      endif

    case 'DerivStep'
      if (! (isnumeric (value) && isreal (value) && isfloat (value)
             && all (value(:) > 0)))
        error (strcat ("statset: option 'DerivStep' must be a real", ...
                       " positive scalar or vector."));
      endif

    case {'RobustWgtFun', 'WgtFun'}
      if (is_function_handle (value))
        return;
      endif
      value = check_string (field, value, {'andrews', 'bisquare', 'cauchy', ...
                                           'fair', 'huber', 'logistic', ...
                                           'talwar', 'welsch'});

    case {'UseParallel', 'UseSubstreams'}
      if (! (isscalar (value) && (islogical (value)
             || (isnumeric (value) && (value == 0 || value == 1)))))
        error (strcat ("statset: option '%s' must be a logical", ...
                       " scalar."), field);
      endif
      value = logical (value);

    case 'OutputFcn'
      if (! (is_function_handle (value)
             || (iscell (value) && all (cellfun (@is_function_handle, value)))))
        error (strcat ("statset: option 'OutputFcn' must be a function", ...
                       " handle or a cell array of function handles."));
      endif

    case 'Streams'
      if (! (isa (value, 'RandStream') || iscell (value)))
        error (strcat ("statset: option 'Streams' must be a random stream", ...
                       " or a cell array of random streams."));
      endif

  endswitch

endfunction

## Validate a character option against the list of values it accepts.
function value = check_string (field, value, valid)
  if (! (ischar (value) || isstring_scalar (value)))
    error (strcat ("statset: option '%s' must be a character vector or", ...
                   " a string scalar."), field);
  endif
  value = char (value);
  idx = find (strcmpi (value, valid));
  if (isempty (idx))
    quoted = strcat ("'", strjoin (valid, "', '"), "'");
    error ("statset: option '%s' must be one of %s.", field, quoted);
  endif
  value = valid{idx};
endfunction

## The default tuning constant of each named robust weight function.
function t = default_tune (wgtfun)
  switch (lower (wgtfun))
    case 'andrews';     t = 1.339;
    case 'bisquare';    t = 4.685;
    case 'cauchy';      t = 2.385;
    case 'fair';        t = 1.400;
    case 'huber';       t = 1.345;
    case 'logistic';    t = 1.205;
    case 'talwar';      t = 2.795;
    case 'welsch';      t = 2.985;
    otherwise;          t = [];
  endswitch
endfunction

## The options each supported function uses by default.
function options = func_defaults (options, fname)

  ## The finite-difference step is eps^(1/3) throughout, and mlecov's is 2^-13.
  dstep = eps ^ (1/3);

  switch (lower (fname))

    case {'crossval', 'lasso', 'lassoglm', 'plsregress', 'rocmetrics'}
      options.Display = 'off';
      options.UseParallel = false;
      options.UseSubstreams = false;

    case {'kmeans', 'kmedoids'}
      options.Display = 'off';
      options.MaxIter = 100;
      options.UseParallel = false;
      options.UseSubstreams = false;

    case 'nnmf'
      options.Display = 'off';
      options.MaxIter = 100;
      options.TolFun = 1e-4;
      options.TolX = 1e-4;
      options.UseParallel = false;
      options.UseSubstreams = false;

    case {'fitglme', 'fitlme', 'fitlmematrix', ...
          'generalizedlinearmixedmodel', 'linearmixedmodel'}
      options.Display = 'off';
      options.MaxIter = 10000;
      options.TolFun = 1e-6;
      options.TolX = 1e-12;

    case {'fitglm', 'glmfit'}
      options.Display = 'off';
      options.MaxIter = 100;
      options.TolX = 1e-6;

    case {'fitnlm', 'nlinfit'}
      options.Display = 'off';
      options.MaxIter = 200;
      options.TolFun = 1e-8;
      options.TolX = 1e-8;
      options.DerivStep = dstep;
      options.FunValCheck = 'on';
      options.Robust = 'off';
      options.WgtFun = 'bisquare';

    case {'gamfit', 'lognfit', 'normfit'}
      options.Display = 'off';
      options.MaxFunEvals = 200;
      options.MaxIter = 100;
      options.TolBnd = 1e-6;
      options.TolFun = 1e-8;
      options.TolX = 1e-8;

    case {'gevfit', 'gpfit', 'nbinfit'}
      options.Display = 'off';
      options.MaxFunEvals = 400;
      options.MaxIter = 200;
      options.TolBnd = 1e-6;
      options.TolFun = 1e-6;
      options.TolX = 1e-6;

    case {'evfit', 'wblfit'}
      options.Display = 'off';
      options.TolX = 1e-6;

    case {'mvncdf', 'mvtcdf'}
      options.Display = 'off';
      options.MaxFunEvals = 1e7;

    case {'pca', 'ppca'}
      options.Display = 'off';
      options.MaxIter = 1000;
      options.TolFun = 1e-6;
      options.TolX = 1e-6;

    case 'coxphfit'
      options.Display = 'off';
      options.MaxFunEvals = 200;
      options.MaxIter = 100;
      options.TolFun = 1e-8;
      options.TolX = 1e-8;

    case 'copulafit'
      options.Display = 'off';
      options.MaxFunEvals = 200;
      options.MaxIter = 100;
      options.TolBnd = 1e-6;
      options.TolX = 1e-6;

    case 'factoran'
      options.Display = 'off';
      options.MaxFunEvals = 400;
      options.MaxIter = 100;
      options.TolFun = 1e-8;
      options.TolX = 1e-8;

    case 'gmdistribution'
      options.Display = 'off';
      options.MaxIter = 100;
      options.TolFun = 1e-6;

    case 'mdscale'
      options.Display = 'off';
      options.MaxIter = 200;
      options.TolFun = 1e-6;
      options.TolX = 1e-6;

    case 'mlecov'
      options.Display = 'off';
      options.GradObj = 'off';
      options.DerivStep = 2 ^ -13;

    case 'mlecustom'
      options.Display = 'off';
      options.MaxFunEvals = 400;
      options.MaxIter = 200;
      options.TolBnd = 1e-6;
      options.TolFun = 1e-6;
      options.TolX = 1e-6;
      options.GradObj = 'off';
      options.DerivStep = dstep;
      options.FunValCheck = 'on';

    case 'tsne'
      options.Display = 'off';
      options.MaxIter = 1000;
      options.TolFun = 1e-10;
      options.OutputFcn = '';

    otherwise
      error ("statset: no default options available for the function '%s'.", ...
             fname);

  endswitch

endfunction

## Print the option names and the values each accepts.
function display_options ()
  tbl = {'Display', '[ {off} | final | iter ]'; ...
         'MaxFunEvals', '[ positive scalar ]'; ...
         'MaxIter', '[ positive scalar ]'; ...
         'TolBnd', '[ positive scalar ]'; ...
         'TolFun', '[ positive scalar ]'; ...
         'TolTypeFun', '[ abs | rel ]'; ...
         'TolX', '[ positive scalar ]'; ...
         'TolTypeX', '[ abs | rel ]'; ...
         'GradObj', '[ {off} | on ]'; ...
         'Jacobian', '[ {off} | on ]'; ...
         'DerivStep', '[ positive scalar or vector ]'; ...
         'FunValCheck', '[ off | {on} ]'; ...
         'Robust', '[ {off} | on ]'; ...
         'RobustWgtFun', ...
         '[ {[]} | andrews | bisquare | cauchy | fair | huber | logistic | talwar | welsch | function handle ]'; ...
         'WgtFun', '[ {[]} | andrews | bisquare | cauchy | fair | huber | logistic | talwar | welsch | function handle ]'; ...
         'Tune', '[ positive scalar ]'; ...
         'UseParallel', '[ {false} | true ]'; ...
         'UseSubstreams', '[ {false} | true ]'; ...
         'Streams', '[ {} | RandStream or cell array ]'; ...
         'OutputFcn', '[ {[]} | function handle or cell array ]'};
  for i = 1:size (tbl, 1)
    printf ("%20s: %s\n", tbl{i,1}, tbl{i,2});
  endfor
endfunction

## Return true for a string scalar.  isa returns false when the class is
## absent, so this does not require the datatypes package to be loaded.
function tf = isstring_scalar (x)
  tf = isa (x, 'string') && isscalar (x);
endfunction

%!demo
%! ## The default options of a given function
%! options = statset ('nlinfit')

%!demo
%! ## Raise the iteration limit of an existing options structure
%! options = statset ('nlinfit');
%! options = statset (options, 'MaxIter', 500);
%! [options.MaxIter, options.TolFun]

## Structure shape
%!test
%! options = statset ();
%! assert_equal (isstruct (options), true);
%! assert_equal (numel (fieldnames (options)), 20);

%!test
%! assert_equal (fieldnames (statset ())', {'Display', 'MaxFunEvals', ...
%! 'MaxIter', 'TolBnd', 'TolFun', 'TolTypeFun', 'TolX', 'TolTypeX', ...
%! 'GradObj', 'Jacobian', 'DerivStep', 'FunValCheck', 'Robust', ...
%! 'RobustWgtFun', 'WgtFun', 'Tune', 'UseParallel', 'UseSubstreams', ...
%! 'Streams', 'OutputFcn'});

%!test
%! options = statset ();
%! assert_equal (options.MaxIter, []);
%! assert_equal (options.Streams, {});

## Per-function defaults, measured against MATLAB R2024a
%!test
%! options = statset ('nlinfit');
%! assert_equal (options.Display, 'off');
%! assert_equal (options.MaxIter, 200);
%! assert_equal (options.TolFun, 1e-8);
%! assert_equal (options.TolX, 1e-8);
%! assert_equal (options.DerivStep, 6.0554544523933429e-06, 1e-20);
%! assert_equal (options.FunValCheck, 'on');
%! assert_equal (options.Robust, 'off');
%! assert_equal (options.WgtFun, 'bisquare');

%!test
%! assert_equal (statset ('fitnlm'), statset ('nlinfit'));

%!test
%! options = statset ('factoran');
%! assert_equal (options.MaxFunEvals, 400);
%! assert_equal (options.MaxIter, 100);
%! assert_equal (options.TolFun, 1e-8);
%! assert_equal (options.TolX, 1e-8);
%! assert_equal (options.TolBnd, []);

%!test
%! options = statset ('fitlme');
%! assert_equal (options.MaxIter, 10000);
%! assert_equal (options.TolFun, 1e-6);
%! assert_equal (options.TolX, 1e-12);

%!test
%! options = statset ('tsne');
%! assert_equal (options.MaxIter, 1000);
%! assert_equal (options.TolFun, 1e-10);
%! assert_equal (options.OutputFcn, '');

%!test
%! options = statset ('nnmf');
%! assert_equal (options.TolFun, 1e-4);
%! assert_equal (options.UseParallel, false);

%!test
%! options = statset ('mlecov');
%! assert_equal (options.GradObj, 'off');
%! assert_equal (options.DerivStep, 0.0001220703125);

%!test
%! options = statset ('coxphfit');
%! assert_equal (options.MaxFunEvals, 200);
%! assert_equal (options.MaxIter, 100);
%! assert_equal (options.TolFun, 1e-8);
%! assert_equal (options.TolX, 1e-8);

%!test
%! options = statset ('mvncdf');
%! assert_equal (options.MaxFunEvals, 1e7);

%!test
%! options = statset ('kmedoids');
%! assert_equal (options.MaxIter, 100);
%! assert_equal (options.UseSubstreams, false);

## A function name is matched case-insensitively
%!test
%! assert_equal (statset ('NLINFIT'), statset ('nlinfit'));

%!test
%! assert_equal (statset ('LinearMixedModel'), statset ('fitlme'));

## Name/value pairs
%!test
%! options = statset ('MaxIter', 100);
%! assert_equal (options.MaxIter, 100);
%! assert_equal (options.TolX, []);

%!test
%! options = statset ('maxiter', 50);
%! assert_equal (options.MaxIter, 50);

%!test
%! options = statset ('MaxIter', 10, 'TolX', 1e-3);
%! assert_equal (options.MaxIter, 10);
%! assert_equal (options.TolX, 1e-3);

## An empty value leaves the option unset
%!test
%! assert_equal (statset ('MaxIter', []), statset ());

## A non-integer iteration count is accepted, as it is by MATLAB
%!test
%! assert_equal (statset ('MaxIter', 2.5).MaxIter, 2.5);

%!test
%! assert_equal (statset ('MaxIter', Inf).MaxIter, Inf);

## Modifying an existing structure leaves the original alone
%!test
%! old = statset ('nlinfit');
%! new = statset (old, 'MaxIter', 999);
%! assert_equal (new.MaxIter, 999);
%! assert_equal (old.MaxIter, 200);
%! assert_equal (new.TolFun, 1e-8);

## Merging two structures: non-empty fields of the second win
%!test
%! old = statset ('nlinfit');
%! new = statset (old, statset ('MaxIter', 777));
%! assert_equal (new.MaxIter, 777);
%! assert_equal (new.TolFun, 1e-8);

## An all-empty second structure overrides nothing
%!test
%! old = statset ('nlinfit');
%! assert_equal (statset (old, statset ()), old);

## Unrecognized fields of a supplied structure are ignored
%!test
%! options = statset (struct ('NotAnOption', 1), 'MaxIter', 5);
%! assert_equal (options.MaxIter, 5);
%! assert_equal (numel (fieldnames (options)), 20);

%!test
%! assert_equal (statset (struct ('NotAnOption', 1)), statset ());

## A named weight function brings its own tuning constant
%!test
%! options = statset ('RobustWgtFun', 'bisquare');
%! assert_equal (options.RobustWgtFun, 'bisquare');
%! assert_equal (options.Tune, 4.685);

%!test
%! assert_equal (statset ('RobustWgtFun', 'huber').Tune, 1.345);

## An explicit tuning constant is not overwritten
%!test
%! options = statset ('RobustWgtFun', 'bisquare', 'Tune', 3);
%! assert_equal (options.Tune, 3);

## A function handle carries no default tuning constant
%!test
%! options = statset ('RobustWgtFun', @(r) 1 ./ (1 + r .^ 2));
%! assert_equal (is_function_handle (options.RobustWgtFun), true);
%! assert_equal (options.Tune, []);

## Enumerated options are matched case-insensitively and stored canonically
%!test
%! assert_equal (statset ('Display', 'ITER').Display, 'iter');

%!test
%! assert_equal (statset ('TolTypeFun', 'REL').TolTypeFun, 'rel');

## A logical option is stored as a logical
%!test
%! options = statset ('UseParallel', 1);
%! assert_equal (islogical (options.UseParallel), true);
%! assert_equal (options.UseParallel, true);

## Error conditions
%!error<statset: no default options available for the function 'nosuchfun'.> ...
%! statset ('nosuchfun')

%!error<statset: no default options available for the function 'MaxIter'.> ...
%! statset ('MaxIter')

%!error<statset: no default options available for the function 'TreeBagger'.> ...
%! statset ('TreeBagger')

%!error<statset: first argument must be a function name, an option name, or an options structure.> ...
%! statset (1, 2)

%!error<statset: arguments must occur in NAME/VALUE pairs.> ...
%! statset ('MaxIter', 5, 'TolX')

%!error<statset: arguments must occur in NAME/VALUE pairs.> ...
%! statset (statset (), 'MaxIter')

%!error<statset: unrecognized option name 'NoSuchOption'.> ...
%! statset ('NoSuchOption', 1)

%!error<statset: option 'MaxIter' must be a real positive scalar.> ...
%! statset ('MaxIter', 'abc')

%!error<statset: option 'MaxIter' must be a real positive scalar.> ...
%! statset ('MaxIter', -1)

%!error<statset: option 'MaxIter' must be a real positive scalar.> ...
%! statset ('MaxIter', int32 (5))

%!error<statset: option 'TolX' must be a real positive scalar.> ...
%! statset ('TolX', -1)

%!error<statset: option 'Tune' must be a real positive scalar.> ...
%! statset ('Tune', 0)

%!error<statset: option 'Display' must be one of 'off', 'final', 'iter'.> ...
%! statset ('Display', 'bogus')

%!error<statset: option 'GradObj' must be one of 'off', 'on'.> ...
%! statset ('GradObj', 'maybe')

%!error<statset: option 'UseParallel' must be a logical scalar.> ...
%! statset ('UseParallel', 'maybe')

%!error<statset: option 'DerivStep' must be a real positive scalar or vector.> ...
%! statset ('DerivStep', -1)

%!error<statset: option 'OutputFcn' must be a function handle or a cell array of function handles.> ...
%! statset ('OutputFcn', 'notafunction')

%!error<statset: option 'Streams' must be a random stream or a cell array of random streams.> ...
%! statset ('Streams', 5)

%!error<statset: option 'RobustWgtFun' must be one of 'andrews', 'bisquare', 'cauchy', 'fair', 'huber', 'logistic', 'talwar', 'welsch'.> ...
%! statset ('RobustWgtFun', 'nosuchweight')

%!error<statset: an options structure must be a scalar structure.> ...
%! statset (struct ('MaxIter', {1, 2}), 'TolX', 1)

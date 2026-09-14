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
## @deftypefn  {statistics} {@var{T} =} templateEnsemble (@var{Method}, @var{NLearn}, @var{Learners})
## @deftypefnx {statistics} {@var{T} =} templateEnsemble (@dots{}, @var{name}, @var{value})
##
## Create a template for an ensemble learner.
##
## @code{@var{T} = templateEnsemble (@var{Method}, @var{NLearn},
## @var{Learners})} returns a template for an ensemble grown by @var{Method}
## from @var{NLearn} learners, each fitted as @var{Learners} says.  A template
## names a learner and the options it is to be fitted with, without fitting
## anything: it is given to @code{fitcecoc}, which grows one such ensemble for
## every binary learner it trains.
##
## @var{Method} is one of the methods of @code{fitcensemble} or
## @code{fitrensemble}, in any letter case.  @var{NLearn} is the number of
## learning cycles and @var{Learners} a learner name, such as
## @qcode{'tree'}, or a template of one, such as @code{templateTree} returns;
## they are the @qcode{'NumLearningCycles'} and @qcode{'Learners'} options of
## @code{fitcensemble}.
##
## @code{@var{T} = templateEnsemble (@dots{}, @var{name}, @var{value})}
## also stores the given options.  They are the name-value arguments of
## @code{fitcensemble}, such as @qcode{'LearnRate'}.
##
## @example
## @group
## T = templateEnsemble ('GentleBoost', 50, templateTree ('MaxNumSplits', 1));
## Mdl = fitcecoc (X, Y, 'Learners', T);
## @end group
## @end example
##
## @var{T} is a structure carrying @qcode{Method}, @qcode{Type},
## @qcode{LearnerTemplates}, @qcode{NLearn} and one field per option given.
## @qcode{Type} is @qcode{'regression'} for LSBoost and
## @qcode{'classification'} for every other method; a @qcode{'Type'} option
## may choose it for Bag, which serves both.
##
## @subheading Deviation from MATLAB
##
## MATLAB returns an object of a class whose name we cannot use, which has no
## public properties and one method this package declines package wide, so a
## structure carries everything a user can observe.  This is what
## @qcode{ModelParameters} already does throughout the package.
##
## Only the method is checked here.  The ensemble owns the list of options it
## takes and checks them, @var{NLearn} and @var{Learners} included, when the
## template is used, so a value it refuses is refused then rather than now.
## MATLAB checks @qcode{'LearnRate'} already here.
##
## @seealso{fitcecoc, fitcensemble, templateTree}
## @end deftypefn

function T = templateEnsemble (Method, NLearn, Learners, varargin)

  if (nargin < 3)
    print_usage ();
  endif

  names = {'AdaBoostM1', 'AdaBoostM2', 'GentleBoost', 'LogitBoost', ...
           'LPBoost', 'TotalBoost', 'RobustBoost', 'RUSBoost', ...
           'Subspace', 'Bag', 'LSBoost'};
  if (! (ischar (Method) && isrow (Method)))
    error ("templateEnsemble: METHOD must be a character vector.");
  endif
  k = find (strcmpi (Method, names));
  if (isempty (k))
    error ("templateEnsemble: '%s' is not an ensemble method.", Method);
  endif
  Method = names{k};

  [opts, errmsg] = templateStruct (Method, varargin);
  if (! isempty (errmsg))
    error ("templateEnsemble: %s", errmsg);
  endif

  Type = methodType (Method);
  given = find (strcmpi ('Type', varargin(1:2:end)));
  if (! isempty (given))
    Type = varargin{2 * given(end)};
    if (! (ischar (Type) && any (strcmp (Type, {'classification', ...
                                                 'regression'}))))
      error (strcat ("templateEnsemble: 'Type' must be 'classification'", ...
                     " or 'regression'."));
    endif
    if (! strcmp (Method, 'Bag') && ! strcmp (Type, methodType (Method)))
      error ("templateEnsemble: the '%s' method is for %s.", Method, ...
             methodType (Method));
    endif
  endif

  T = struct ('Method', Method, 'Type', Type, ...
              'LearnerTemplates', {Learners}, 'NLearn', NLearn);
  for [val, name] = opts
    if (! any (strcmpi (name, {'Method', 'Type'})))
      T.(name) = val;
    endif
  endfor

endfunction

## The kind of model a method grows, Bag growing either.
function type = methodType (Method)

  if (strcmp (Method, 'LSBoost'))
    type = 'regression';
  else
    type = 'classification';
  endif

endfunction

## Tests
%!test  # a template names its method, type, learners and cycles
%! T = templateEnsemble ('AdaBoostM1', 20, 'tree');
%! assert_equal (class (T), 'struct');
%! assert_equal (fieldnames (T), ...
%!               {'Method'; 'Type'; 'LearnerTemplates'; 'NLearn'});
%! assert_equal (T.Method, 'AdaBoostM1');
%! assert_equal (T.Type, 'classification');
%! assert_equal (T.LearnerTemplates, 'tree');
%! assert_equal (T.NLearn, 20);

%!test  # the method is matched in any letter case
%! T = templateEnsemble ('gentleboost', 10, 'tree');
%! assert_equal (T.Method, 'GentleBoost');

%!test  # LSBoost grows a regression ensemble
%! T = templateEnsemble ('LSBoost', 10, 'tree');
%! assert_equal (T.Type, 'regression');

%!test  # Bag takes its type from the 'Type' option
%! T = templateEnsemble ('Bag', 10, 'tree', 'Type', 'regression');
%! assert_equal (T.Type, 'regression');
%! assert_equal (isfield (T, 'Type'), true);
%! assert_equal (numfields (T), 4);

%!test  # a learner template and the options are stored as they stand
%! S = templateTree ('MaxNumSplits', 1);
%! T = templateEnsemble ('GentleBoost', 5, S, 'LearnRate', 0.5);
%! assert_equal (T.LearnerTemplates, S);
%! assert_equal (T.LearnRate, 0.5);

%!test  # a value the ensemble would refuse is not refused here
%! T = templateEnsemble ('AdaBoostM1', 0, 'svm');
%! assert_equal (T.NLearn, 0);

## Test input validation
%!error<Invalid call to templateEnsemble> templateEnsemble ('AdaBoostM1', 10)
%!error<templateEnsemble: METHOD must be a character vector.> ...
%! templateEnsemble (1, 10, 'tree')
%!error<templateEnsemble: 'Foo' is not an ensemble method.> ...
%! templateEnsemble ('Foo', 10, 'tree')
%!error<templateEnsemble: name-value arguments must be in pairs.> ...
%! templateEnsemble ('AdaBoostM1', 10, 'tree', 'LearnRate')
%!error<templateEnsemble: 'Type' must be 'classification' or 'regression'.> ...
%! templateEnsemble ('Bag', 10, 'tree', 'Type', 'both')
%!error<templateEnsemble: the 'LSBoost' method is for regression.> ...
%! templateEnsemble ('LSBoost', 10, 'tree', 'Type', 'classification')

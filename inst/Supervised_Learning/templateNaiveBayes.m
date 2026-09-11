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
## @deftypefn  {statistics} {@var{T} =} templateNaiveBayes ()
## @deftypefnx {statistics} {@var{T} =} templateNaiveBayes (@var{name}, @var{value})
##
## Create a template for a naive Bayes classifier.
##
## @code{@var{T} = templateNaiveBayes ()} returns a template carrying the
## default options of @code{ClassificationNaiveBayes}.  A template names a
## learner and the options it is to be fitted with, without fitting anything: it
## is given to a function that fits many models, such as @code{fitcecoc}, which
## uses it for every binary learner it trains.
##
## @code{@var{T} = templateNaiveBayes (@var{name}, @var{value})} also stores
## the given options.  They are the name-value arguments of
## @code{fitcnb}, and any of them may be given here instead.
##
## @example
## @group
## T = templateNaiveBayes ('DistributionNames', 'kernel');
## Mdl = fitcecoc (X, Y, 'Learners', T);
## @end group
## @end example
##
## @var{T} is a structure carrying @qcode{Method}, @qcode{Type} and one
## field per option given.
##
## @subheading Deviation from MATLAB
##
## MATLAB returns an object of a class whose name we cannot use, which has no
## public properties and one method this package declines package wide, so a
## structure carries everything a user can observe.  This is what
## @qcode{ModelParameters} already does throughout the package.
##
## An option name is not checked here.  @code{ClassificationNaiveBayes} owns
## the list of options it takes and checks it when the template is used, so a
## name it does not know is refused then rather than now.
##
## @seealso{fitcecoc, ClassificationNaiveBayes, fitcnb}
## @end deftypefn

function T = templateNaiveBayes (varargin)

  [T, errmsg] = templateStruct ('NaiveBayes', varargin);
  if (! isempty (errmsg))
    error ("templateNaiveBayes: %s", errmsg);
  endif

endfunction

## Tests
%!test  # the default template names its learner and nothing else
%! T = templateNaiveBayes ();
%! assert_equal (class (T), 'struct');
%! assert_equal (T.Method, 'NaiveBayes');
%! assert_equal (T.Type, 'classification');
%! assert_equal (numfields (T), 2);

%!test  # an option given is stored under its own name, as it stands
%! T = templateNaiveBayes ('DistributionNames', 'kernel');
%! assert_equal (numfields (T), 3);
%! assert_equal (T.DistributionNames, 'kernel');

%!test  # a name the learner does not know is not refused here
%! T = templateNaiveBayes ('NoSuchOption', 42);
%! assert_equal (T.NoSuchOption, 42);

## Test input validation
%!error<templateNaiveBayes: name-value arguments must be in pairs.> ...
%! templateNaiveBayes ('KernelScale')
%!error<templateNaiveBayes: parameter names must be character vectors.> ...
%! templateNaiveBayes (42, 1)
%!error<templateNaiveBayes: 'not a name' is not a valid parameter name.> ...
%! templateNaiveBayes ('not a name', 1)

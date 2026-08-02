## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {ClassificationSVM} {@var{obj} = } loadmodel (@var{filename})
##
## Load a Classification or Regression model from a file.
##
## @code{@var{obj} = loadmodel (@var{filename})} loads a Classification or
## Regression object, @var{obj}, from a file defined in @var{filename}.
##
## @seealso{savemodel, ClassificationDiscriminant, ClassificationGAM,
## ClassificationKNN, ClassificationNeuralNetwork,
## ClassificationPartitionedModel, ClassificationSVM, RegressionGAM}
## @end deftypefn

function obj = loadmodel (filename)

  ## Check input parameters
  if (nargin < 1)
    error ("loadmodel: too few arguments.");
  endif

  ## Supported Classification and Regression objects
  supported = {'ClassificationKNN'};

  ## Read file into a structure
  data = load (filename);

  ## Check that 'classdef_name' variable exists and that it
  ## contains a valid Classification or Regression object
  if (! isfield (data, 'classdef_name'))
    msg = ' ''%s'' does not contain a Classification or Regression object.';
    error (strcat ("loadmodel:", msg), filename);
  endif

  ## Remove 'classdef_name' field from data structure
  classdef_name = data.classdef_name;
  data = rmfield (data, 'classdef_name');

  ## Parse data structure to the static load method of specified classdef
  switch (classdef_name)

    case 'ClassificationDiscriminant'
      obj = ClassificationDiscriminant.load_model (filename, data);

    case 'CompactClassificationDiscriminant'
      obj = CompactClassificationDiscriminant.load_model (filename, data);

    case 'ClassificationGAM'
      obj = ClassificationGAM.load_model (filename, data);

    case 'CompactClassificationGAM'
      obj = CompactClassificationGAM.load_model (filename, data);

    case 'ClassificationKNN'
      obj = ClassificationKNN.load_model (filename, data);

    case 'ClassificationNeuralNetwork'
      obj = ClassificationNeuralNetwork.load_model (filename, data);

    case 'CompactClassificationNeuralNetwork'
      obj = CompactClassificationNeuralNetwork.load_model (filename, data);

    case 'ClassificationSVM'
      obj = ClassificationSVM.load_model (filename, data);

    case 'CompactClassificationSVM'
      obj = CompactClassificationSVM.load_model (filename, data);

    case 'RegressionGAM'
      obj = RegressionGAM.load_model (filename, data);

    otherwise
      error ("loadmodel: '%s' is not supported.", classdef_name);

  endswitch

endfunction

## A saved model must load back as the same class with the same state.  Every
## loader but the neural-network one compared the saved fieldnames against
## fieldnames (mdl) for exact equality; a private property such as STname is
## saved but never reported by fieldnames, so the comparison could not match
## and the load always failed.
%!test
%! load fisheriris
%! Yb = strcmp (species, 'setosa');
%! mdls = {fitcknn(meas, species, 'ScoreTransform', 'logit'), ...
%!         fitcdiscr(meas, species, 'ScoreTransform', 'logit'), ...
%!         fitcsvm(meas(1:100,:), Yb(1:100), 'ScoreTransform', 'logit'), ...
%!         fitcgam(meas(1:100,:), Yb(1:100))};
%! for i = 1:numel (mdls)
%!   m = mdls{i};
%!   fn = tempname ();
%!   unwind_protect
%!     savemodel (m, fn);
%!     m2 = loadmodel (fn);
%!     assert_equal (class (m2), class (m));
%!     for f = fieldnames (m)'
%!       a = m.(f{1});
%!       b = m2.(f{1});
%!       if (is_function_handle (a))
%!         assert_equal (func2str (b), func2str (a));
%!       else
%!         assert_equal (b, a);
%!       endif
%!     endfor
%!   unwind_protect_cleanup
%!     if (exist (fn, 'file'))
%!       delete (fn);
%!     endif
%!   end_unwind_protect
%! endfor

## A compact model must come back compact.  CompactClassificationSVM was
## dispatched to ClassificationSVM.load_model, which builds the wrong class.
%!test
%! load fisheriris
%! Yb = strcmp (species, 'setosa');
%! c = compact (fitcsvm (meas(1:100,:), Yb(1:100)));
%! fn = tempname ();
%! unwind_protect
%!   savemodel (c, fn);
%!   c2 = loadmodel (fn);
%!   assert_equal (class (c2), 'CompactClassificationSVM');
%!   assert_equal (predict (c2, meas(1:10,:)), predict (c, meas(1:10,:)));
%! unwind_protect_cleanup
%!   if (exist (fn, 'file'))
%!     delete (fn);
%!   endif
%! end_unwind_protect

## The score transform's name is private state and must survive.  It was not
## saved at all by ClassificationSVM, so it silently reverted to 'none'.
%!test
%! load fisheriris
%! Yb = strcmp (species, 'setosa');
%! m = fitcsvm (meas(1:100,:), Yb(1:100), 'ScoreTransform', 'logit');
%! fn = tempname ();
%! unwind_protect
%!   savemodel (m, fn);
%!   m2 = loadmodel (fn);
%!   assert_equal (strfind (evalc ('disp (m2)'), "'logit'") > 0, true);
%! unwind_protect_cleanup
%!   if (exist (fn, 'file'))
%!     delete (fn);
%!   endif
%! end_unwind_protect

## ClassificationGAM did not save LearningRate or NumIterations at all.
%!test
%! load fisheriris
%! Yb = strcmp (species, 'setosa');
%! m = fitcgam (meas(1:100,:), Yb(1:100));
%! fn = tempname ();
%! unwind_protect
%!   savemodel (m, fn);
%!   m2 = loadmodel (fn);
%!   assert_equal (m2.LearningRate, m.LearningRate);
%!   assert_equal (m2.NumIterations, m.NumIterations);
%! unwind_protect_cleanup
%!   if (exist (fn, 'file'))
%!     delete (fn);
%!   endif
%! end_unwind_protect

## Test input validation
%!error<loadmodel: too few arguments.> loadmodel ()
%!error<loadmodel: 'fisheriris.mat' does not contain a Classification or Regression object.> ...
%! loadmodel ('fisheriris.mat')
%!error<loadmodel: 'ClassificationModel' is not supported.> ...
%! loadmodel ('fail_loadmodel.mdl')
%!error<ClassificationKNN.load_model: invalid model in 'fail_load_model.mdl'.> ...
%! loadmodel ('fail_load_model.mdl')


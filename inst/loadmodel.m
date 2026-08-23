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

    case 'ClassificationKernel'
      obj = ClassificationKernel.load_model (filename, data);

    case 'ClassificationKNN'
      obj = ClassificationKNN.load_model (filename, data);

    case 'ClassificationLinear'
      obj = ClassificationLinear.load_model (filename, data);

    case 'ClassificationNeuralNetwork'
      obj = ClassificationNeuralNetwork.load_model (filename, data);

    case 'CompactClassificationNeuralNetwork'
      obj = CompactClassificationNeuralNetwork.load_model (filename, data);

    case 'ClassificationSVM'
      obj = ClassificationSVM.load_model (filename, data);

    case 'CompactClassificationSVM'
      obj = CompactClassificationSVM.load_model (filename, data);

    case 'RegressionGP'
      obj = RegressionGP.load_model (filename, data);

    case 'CompactRegressionGP'
      obj = CompactRegressionGP.load_model (filename, data);

    case 'RegressionGAM'
      obj = RegressionGAM.load_model (filename, data);

    case 'RegressionKernel'
      obj = RegressionKernel.load_model (filename, data);

    case 'RegressionLinear'
      obj = RegressionLinear.load_model (filename, data);

    case 'CompactRegressionGAM'
      obj = CompactRegressionGAM.load_model (filename, data);

    case 'RegressionNeuralNetwork'
      obj = RegressionNeuralNetwork.load_model (filename, data);

    case 'RegressionSVM'
      obj = RegressionSVM.load_model (filename, data);

    case 'CompactRegressionSVM'
      obj = CompactRegressionSVM.load_model (filename, data);

    case 'CompactRegressionNeuralNetwork'
      obj = CompactRegressionNeuralNetwork.load_model (filename, data);

    otherwise
      error ("loadmodel: '%s' is not supported.", classdef_name);

  endswitch

endfunction

## A saved model must load back as the same class with the same state.  Every
## loader but the neural-network one compared the saved fieldnames against
## fieldnames (mdl) for exact equality; a private property such as STname is
## saved but never reported by fieldnames, so the comparison could not match
## and the load always failed.
## Every property of a saved ClassificationKNN comes back as it was.
%!test
%! load fisheriris
%! Yb = strcmp (species, 'setosa');
%! m = fitcknn (meas, species, 'ScoreTransform', 'logit');
%! fn = tempname ();
%! unwind_protect
%!   savemodel (m, fn);
%!   m2 = loadmodel (fn);
%!   assert_equal (class (m2), class (m));
%! assert_equal (m2.W, m.W);
%! assert_equal (m2.X, m.X);
%! assert_equal (m2.Y, m.Y);
%! assert_equal (m2.NumObservations, m.NumObservations);
%! assert_equal (m2.RowsUsed, m.RowsUsed);
%! assert_equal (m2.NumPredictors, m.NumPredictors);
%! assert_equal (m2.PredictorNames, m.PredictorNames);
%! assert_equal (m2.ResponseName, m.ResponseName);
%! assert_equal (m2.ClassNames, m.ClassNames);
%! assert_equal (m2.Sigma, m.Sigma);
%! assert_equal (m2.Mu, m.Mu);
%! assert_equal (m2.BreakTies, m.BreakTies);
%! assert_equal (m2.NumNeighbors, m.NumNeighbors);
%! assert_equal (m2.Distance, m.Distance);
%! assert_equal (m2.DistanceWeight, m.DistanceWeight);
%! assert_equal (m2.DistParameter, m.DistParameter);
%! assert_equal (m2.NSMethod, m.NSMethod);
%! assert_equal (m2.IncludeTies, m.IncludeTies);
%! assert_equal (m2.BucketSize, m.BucketSize);
%! assert_equal (m2.Cost, m.Cost);
%! assert_equal (m2.Prior, m.Prior);
%! assert_equal (m2.ScoreTransform, m.ScoreTransform);
%! unwind_protect_cleanup
%!   if (exist (fn, 'file'))
%!     delete (fn);
%!   endif
%! end_unwind_protect

## Every property of a saved ClassificationDiscriminant comes back as it was.
%!test
%! load fisheriris
%! Yb = strcmp (species, 'setosa');
%! m = fitcdiscr (meas, species, 'ScoreTransform', 'logit');
%! fn = tempname ();
%! unwind_protect
%!   savemodel (m, fn);
%!   m2 = loadmodel (fn);
%!   assert_equal (class (m2), class (m));
%! assert_equal (m2.W, m.W);
%! assert_equal (m2.X, m.X);
%! assert_equal (m2.Y, m.Y);
%! assert_equal (m2.NumObservations, m.NumObservations);
%! assert_equal (m2.RowsUsed, m.RowsUsed);
%! assert_equal (m2.NumPredictors, m.NumPredictors);
%! assert_equal (m2.PredictorNames, m.PredictorNames);
%! assert_equal (m2.ResponseName, m.ResponseName);
%! assert_equal (m2.ClassNames, m.ClassNames);
%! assert_equal (m2.Sigma, m.Sigma);
%! assert_equal (m2.Mu, m.Mu);
%! assert_equal (m2.Coeffs, m.Coeffs);
%! assert_equal (m2.Delta, m.Delta);
%! assert_equal (m2.DiscrimType, m.DiscrimType);
%! assert_equal (m2.Gamma, m.Gamma);
%! assert_equal (m2.MinGamma, m.MinGamma);
%! assert_equal (m2.LogDetSigma, m.LogDetSigma);
%! assert_equal (m2.XCentered, m.XCentered);
%! assert_equal (m2.Cost, m.Cost);
%! assert_equal (m2.Prior, m.Prior);
%! assert_equal (m2.ScoreTransform, m.ScoreTransform);
%! unwind_protect_cleanup
%!   if (exist (fn, 'file'))
%!     delete (fn);
%!   endif
%! end_unwind_protect

## Every property of a saved ClassificationSVM comes back as it was.
%!test
%! load fisheriris
%! Yb = strcmp (species, 'setosa');
%! m = fitcsvm (meas(1:100,:), Yb(1:100), 'ScoreTransform', 'logit');
%! fn = tempname ();
%! unwind_protect
%!   savemodel (m, fn);
%!   m2 = loadmodel (fn);
%!   assert_equal (class (m2), class (m));
%! assert_equal (m2.X, m.X);
%! assert_equal (m2.Y, m.Y);
%! assert_equal (m2.NumObservations, m.NumObservations);
%! assert_equal (m2.RowsUsed, m.RowsUsed);
%! assert_equal (m2.NumPredictors, m.NumPredictors);
%! assert_equal (m2.PredictorNames, m.PredictorNames);
%! assert_equal (m2.ResponseName, m.ResponseName);
%! assert_equal (m2.ClassNames, m.ClassNames);
%! assert_equal (m2.Sigma, m.Sigma);
%! assert_equal (m2.Mu, m.Mu);
%! assert_equal (m2.ModelParameters, m.ModelParameters);
%! assert_equal (m2.Model, m.Model);
%! assert_equal (m2.Alpha, m.Alpha);
%! assert_equal (m2.Beta, m.Beta);
%! assert_equal (m2.Bias, m.Bias);
%! assert_equal (m2.IsSupportVector, m.IsSupportVector);
%! assert_equal (m2.SupportVectorLabels, m.SupportVectorLabels);
%! assert_equal (m2.SupportVectors, m.SupportVectors);
%! assert_equal (m2.Prior, m.Prior);
%! assert_equal (m2.Cost, m.Cost);
%! assert_equal (m2.W, m.W);
%! assert_equal (m2.CategoricalPredictors, m.CategoricalPredictors);
%! assert_equal (m2.ExpandedPredictorNames, m.ExpandedPredictorNames);
%! assert_equal (m2.ScoreTransform, m.ScoreTransform);
%! unwind_protect_cleanup
%!   if (exist (fn, 'file'))
%!     delete (fn);
%!   endif
%! end_unwind_protect

## Every property of a saved ClassificationGAM comes back as it was.
%!test
%! load fisheriris
%! Yb = strcmp (species, 'setosa');
%! m = fitcgam (meas(1:100,:), Yb(1:100));
%! fn = tempname ();
%! unwind_protect
%!   savemodel (m, fn);
%!   m2 = loadmodel (fn);
%!   assert_equal (class (m2), class (m));
%! assert_equal (m2.X, m.X);
%! assert_equal (m2.Y, m.Y);
%! assert_equal (m2.NumObservations, m.NumObservations);
%! assert_equal (m2.RowsUsed, m.RowsUsed);
%! assert_equal (m2.NumPredictors, m.NumPredictors);
%! assert_equal (m2.PredictorNames, m.PredictorNames);
%! assert_equal (m2.ResponseName, m.ResponseName);
%! assert_equal (m2.ClassNames, m.ClassNames);
%! assert_equal (m2.Prior, m.Prior);
%! assert_equal (m2.Formula, m.Formula);
%! assert_equal (m2.Interactions, m.Interactions);
%! assert_equal (m2.Knots, m.Knots);
%! assert_equal (m2.Order, m.Order);
%! assert_equal (m2.DoF, m.DoF);
%! assert_equal (m2.LearningRate, m.LearningRate);
%! assert_equal (m2.NumIterations, m.NumIterations);
%! assert_equal (m2.Intercept, m.Intercept);
%! assert_equal (m2.W, m.W);
%! assert_equal (m2.CategoricalPredictors, m.CategoricalPredictors);
%! assert_equal (m2.ExpandedPredictorNames, m.ExpandedPredictorNames);
%! assert_equal (m2.BaseModel, m.BaseModel);
%! assert_equal (m2.ModelwInt, m.ModelwInt);
%! assert_equal (m2.IntMatrix, m.IntMatrix);
%! assert_equal (m2.Cost, m.Cost);
%! assert_equal (m2.ScoreTransform, m.ScoreTransform);
%! unwind_protect_cleanup
%!   if (exist (fn, 'file'))
%!     delete (fn);
%!   endif
%! end_unwind_protect

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


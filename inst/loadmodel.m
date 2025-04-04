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
  supported = {"ClassificationKNN"};

  ## Read file into a structure
  data = load (filename);

  ## Check that 'classdef_name' variable exists and that it
  ## contains a valid Classification or Regression object
  if (! isfield (data, "classdef_name"))
    msg = " '%s' does not contain a Classification or Regression object.";
    error (strcat ("loadmodel:", msg), filename);
  endif

  ## Remove 'classdef_name' field from data structure
  classdef_name = data.classdef_name;
  data = rmfield (data, "classdef_name");

  ## Parse data structure to the static load method of specified classdef
  switch (classdef_name)

    case "ClassificationDiscriminant"
      obj = ClassificationDiscriminant.load_model (filename, data);

    case "CompactClassificationDiscriminant"
      obj = CompactClassificationDiscriminant.load_model (filename, data);

    case "ClassificationGAM"
      obj = ClassificationGAM.load_model (filename, data);

    case "CompactClassificationGAM"
      obj = CompactClassificationGAM.load_model (filename, data);

    case "ClassificationKNN"
      obj = ClassificationKNN.load_model (filename, data);

    case "ClassificationNeuralNetwork"
      obj = ClassificationNeuralNetwork.load_model (filename, data);

    case "CompactClassificationNeuralNetwork"
      obj = CompactClassificationNeuralNetwork.load_model (filename, data);

    case "ClassificationSVM"
      obj = ClassificationSVM.load_model (filename, data);

    case "CompactClassificationSVM"
      obj = ClassificationSVM.load_model (filename, data);

    case "RegressionGAM"
      obj = RegressionGAM.load_model (filename, data);

    otherwise
      error ("loadmodel: '%s' is not supported.", classdef_name);

  endswitch

endfunction

## Test input validation
%!error<loadmodel: too few arguments.> loadmodel ()
%!error<loadmodel: 'fisheriris.mat' does not contain a Classification or Regression object.> ...
%! loadmodel ("fisheriris.mat")
%!error<loadmodel: 'ClassificationModel' is not supported.> ...
%! loadmodel ("fail_loadmodel.mdl")
%!error<ClassificationKNN.load_model: invalid model in 'fail_load_model.mdl'.> ...
%! loadmodel ("fail_load_model.mdl")


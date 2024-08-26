/*
Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>

This file is part of the statistics package for GNU Octave.

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program; if not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <memory>
#include <cmath>
#include <octave/oct.h>
#include <octave/parse.h>
#include <octave/ov-struct.h>
#include "fcnn.cpp"

using namespace std;

DEFUN_DLD(fcnnpredict, args, nargout,
          "-*- texinfo -*-\n\
 @deftypefn  {statistics} {@var{pred_Y} =} fcnnpredict (@var{Mdl}, @var{XC})\n\
 @deftypefnx {statistics} {[@var{pred_Y}, @var{scores}] =} fcnnpredict @\
 (@var{Mdl}, @var{XC})\n\
\n\
\n\
Make predictions from a fully connected Neural Network. \n\
\n\n\
\n\n\
@code{@var{pred_Y} = fcnnpredict (@var{Mdl}, @var{XC})} requires the following \
input arguments.\
\n\n\
@itemize \n\
@item @var{Mdl} : A structure containing the trained model parameters as \
generated by the @code{fcnntrain} function. \
\n\
\n\
@item @var{X} : An @math{NxM} matrix containing the data set to be predicted \
upon.  Rows @math{N} correspond to individual samples and columns @math{M} \
correspond to features (dimensions).  Type of @var{X} must be double and the \
number of features must correspond to those of the trained model. \n\
@end itemize \n\
\n\
@code{fcnnpredict} returns the predicted labels, @var{pred_Y}, and if a second \
output argument is requested, it also returns the corresponding values of the \
neural networks output in @var{scores}. \
\n\
\n\
@seealso{fcnntrain, fitcnet, ClassificationNeuralNetwork} \n\
@end deftypefn")
{
  // Check for correct number of input/output arguments
  if (args.length () < 2)
  {
    error ("fcnnpredict: too few input arguments.");
  }
  if (nargout > 2)
  {
    error ("fcnnpredict: too many output arguments.");
  }

  // Do some input validation while loading the trained model
  if (! (args(0).isstruct () && args(0).numel () == 1))
  {
    error ("fcnnpredict: first argument must be a scalar structure.");
  }
  octave_scalar_map fcnn_model = args(0).scalar_map_value ();
  if (! fcnn_model.isfield ("LayerWeights"))
  {
    error ("fcnnpredict: model does not have a 'LayerWeights' field.");
  }
  if (! fcnn_model.getfield("LayerWeights").iscell () ||
      ! (fcnn_model.getfield("LayerWeights").rows () == 1 &&
         fcnn_model.getfield("LayerWeights").columns () > 1))
  {
    error ("fcnnpredict: 'LayerWeights' must be a cell row vector.");
  }
  Cell LayerWeights = fcnn_model.getfield("LayerWeights").cell_value();
  if (! fcnn_model.isfield ("Activations"))
  {
    error ("fcnnpredict: model does not have an 'Activations' field.");
  }
  if (! fcnn_model.getfield("Activations").isnumeric () ||
      ! (fcnn_model.getfield("Activations").rows () == 1 &&
         fcnn_model.getfield("Activations").columns () > 1))
  {
    error ("fcnnpredict: 'Activations' must be a numeric row vector.");
  }
  RowVector ActiveCode = fcnn_model.getfield("Activations").row_vector_value();

  // Do some input validation while loading the testing data
  if (! args(1).isnumeric () || args(1).iscomplex ())
  {
    error ("fcnnpredict: XC must be a real numeric matrix.");
  }
  if (args(1).isempty ())
  {
    error ("fcnnpredict: XC cannot be empty.");
  }
  if (args(1).columns () != LayerWeights.elem(0).columns () - 1)
  {
    error ("fcnnpredict: the features in XC do not match the trained model.");
  }
  Matrix X = args(1).matrix_value ();
  int n = args(1).rows ();
  int d = args(1).columns ();

  // Construct 2D vector from data in XC
  vector<vector<double>> data (n, vector<double>(d, 0));
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < d; j++)
    {
      data[i][j] = X(i,j);
    }
  }

  // Create a vector of layers sized appropriately
  vector<DenseLayer> WeightBias;
  vector<ActivationLayer> Activation;
  int numlayers = LayerWeights.numel ();
  int input_size = d;
  int output_size;
  for (int i = 0; i < numlayers; i++)
  {
    Matrix WB = LayerWeights.elem(i).matrix_value ();
    output_size = (int) WB.rows ();
    // Construct 2D vector with Weights-Biases
    vector<vector<double>> Wb_matrix;
    for (int r = 0; r < WB.rows (); r++)
    {
      vector<double> WB_vector;
      for (int c = 0; c < WB.columns (); c++)
      {
        WB_vector.push_back (WB(r,c));
      }
      Wb_matrix.push_back (WB_vector);
    }
    // Create dense layer and set its values
    DenseLayer DL = DenseLayer (input_size, output_size);
    DL.set_layer (Wb_matrix);
    // Create activation layer
    ActivationLayer AL = ActivationLayer (ActiveCode(i), 0.01);
    WeightBias.push_back (DL);
    Activation.push_back (AL);
    input_size = output_size;
  }

  // Initialize Prediction and Score
  vector<int> predictions = vector<int> ();
  vector<vector<double>> scores;

  // Go through all testing samples
  for (int sample_idx = 0; sample_idx < n; sample_idx++)
  {
    vector<double> sample = data[sample_idx];

    // Forward pass
    for (int layer_idx = 0; layer_idx < numlayers; layer_idx++)
    {
      sample = WeightBias[layer_idx].forward (sample);
      sample = Activation[layer_idx].forward (sample);
    }
    // Save scores
    scores.push_back (sample);

    // Get the prediction for this sample and store it to vector
    int prediction = 0;           // Search for highest value
    for (int j = 0; j < output_size; j++)
    {
      if (sample[j] > sample[prediction])
      {
        prediction = j;
      }
    }
    predictions.push_back (prediction + 1);
  }

  // Store predicted labels in ColumnVector
  ColumnVector Y_pred(n);
  for (int sample_idx = 0; sample_idx < n; sample_idx++)
  {
    Y_pred(sample_idx) = predictions[sample_idx];
  }

  // Store predicted scores in Matrix
  Matrix Y_scores(n,output_size);
  for (int r = 0; r < n; r++)
  {
    for (int c = 0; c < output_size; c++)
    {
      Y_scores(r,c) = scores[r][c];
    }
  }

  // Prepare returning arguments
  octave_value_list retval (nargout);
  retval(0) = Y_pred;
  if (nargout > 0)
  {
    retval(1) = Y_scores;
  }
  return retval;
}

/*
%!shared X, Y, MODEL
%! load fisheriris
%! X = meas;
%! Y = grp2idx (species);
%! MODEL = fcnntrain (X, Y, 10, [1, 1], 0.025, 100, false);
%!test
%! [Y_pred, Y_scores] = fcnnpredict (MODEL, X);
%! assert (numel (Y_pred), numel (Y));
%! assert (isequal (size (Y_pred), size (Y)), true);
%! assert (columns (Y_scores), numel (unique (Y)));
%! assert (rows (Y_scores), numel (Y));
%!error <fcnnpredict: too few input arguments.> ...
%! fcnnpredict (MODEL);
%!error <fcnnpredict: too many output arguments.> ...
%! [Q, W, E] = fcnnpredict (MODEL, X);
%!error <fcnnpredict: first argument must be a scalar structure.> ...
%! fcnnpredict (1, X);
%!error <fcnnpredict: first argument must be a scalar structure.> ...
%! fcnnpredict (struct ("L", {1, 2, 3}), X);
%!error <fcnnpredict: model does not have a 'LayerWeights' field.> ...
%! fcnnpredict (struct ("L", 1), X);
%!error <fcnnpredict: 'LayerWeights' must be a cell row vector.> ...
%! fcnnpredict (struct ("LayerWeights", 1), X);
%!error <fcnnpredict: 'LayerWeights' must be a cell row vector.> ...
%! fcnnpredict (struct ("LayerWeights", {1}), X);
%!error <fcnnpredict: 'LayerWeights' must be a cell row vector.> ...
%! fcnnpredict (struct ("LayerWeights", {{1; 2; 3}}), X);
%!error <fcnnpredict: model does not have an 'Activations' field.> ...
%! fcnnpredict (struct ("LayerWeights", {[{ones(3)},{ones(3)}]}, "R", 2), X);
%!error <fcnnpredict: 'Activations' must be a numeric row vector.> ...
%! fcnnpredict (struct ("LayerWeights", {[{ones(3)},{ones(3)}]}, ...
%!                      "Activations", [2]), X);
%!error <fcnnpredict: 'Activations' must be a numeric row vector.> ...
%! fcnnpredict (struct ("LayerWeights", {[{ones(3)},{ones(3)}]}, ...
%!                      "Activations", [2; 2]), X);
%!error <fcnnpredict: 'Activations' must be a numeric row vector.> ...
%! fcnnpredict (struct ("LayerWeights", {[{ones(3)},{ones(3)}]}, ...
%!                      "Activations", {{2, 2}}), X);
%!error <fcnnpredict: 'Activations' must be a numeric row vector.> ...
%! fcnnpredict (struct ("LayerWeights", {[{ones(3)},{ones(3)}]}, ...
%!                      "Activations", {{"sigmoid", "softmax"}}), X);
%!error <fcnnpredict: 'Activations' must be a numeric row vector.> ...
%! fcnnpredict (struct ("LayerWeights", {[{ones(3)},{ones(3)}]}, ...
%!                      "Activations", "sigmoid"), X);
%!error <fcnnpredict: XC must be a real numeric matrix.> ...
%! fcnnpredict (MODEL, complex (X));
%!error <fcnnpredict: XC must be a real numeric matrix.> ...
%! fcnnpredict (MODEL, {1, 2, 3, 4});
%!error <fcnnpredict: XC must be a real numeric matrix.> ...
%! fcnnpredict (MODEL, "asd");
%!error <fcnnpredict: XC cannot be empty.> ...
%! fcnnpredict (MODEL, []);
%!error <fcnnpredict: the features in XC do not match the trained model.> ...
%! fcnnpredict (MODEL, X(:,[1:3]));
*/
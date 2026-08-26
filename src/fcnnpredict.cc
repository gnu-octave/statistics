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
 @deftypefn  {statistics} {@var{pred_Y} =} fcnnpredict (@var{LayerWeights}, @\
 @var{LayerBiases}, @var{Activations}, @var{OutputLayerActivation}, @var{XC})\n\
 @deftypefnx {statistics} {@var{pred_Y} =} fcnnpredict @\
 (@dots{}, @var{NumThreads})\n\
 @deftypefnx {statistics} {[@var{pred_Y}, @var{scores}] =} fcnnpredict (@dots{})\n\
\n\
\n\
Make predictions from a fully connected Neural Network. \n\
\n\n\
\n\n\
@code{@var{pred_Y} = fcnnpredict (@var{LayerWeights}, @var{LayerBiases}, \
@var{Activations}, @var{OutputLayerActivation}, @var{XC})} requires the \
following input arguments.\
\n\n\
@itemize \n\
@item @var{LayerWeights} : A cell row vector holding one matrix per layer, \
each with one row per neuron of that layer and one column per input to it. \
\n\
\n\
@item @var{LayerBiases} : A cell row vector holding one bias column per \
layer, matching @var{LayerWeights} layer for layer and row for row. \
\n\
\n\
@item @var{Activations} : The activation function of the hidden layers, \
named as a character vector applying to all of them or as a cellstring naming \
them one by one.  The supported names are listed under @code{fcnntrain}. \
\n\
\n\
@item @var{OutputLayerActivation} : The activation function of the output \
layer, named as a character vector. \
\n\
\n\
@item @var{XC} : An @math{NxM} matrix containing the data set to be predicted \
upon.  Rows @math{N} correspond to individual samples and columns @math{M} \
correspond to features (dimensions).  Type of @var{XC} must be double and the \
number of features must correspond to those of the trained model. \n\
@end itemize \n\
@code{fcnnpredict} can also be called with a sixth input argument, in which \
case, @var{NumThreads}, a positive scalar integer value, defines the number \
of threads to be used when computing the activation layers.  For layers with \
less than 1000 neurons, @var{NumThreads} always defaults to 1. \
\n\
@code{fcnnpredict} returns the predicted labels, @var{pred_Y}, and if a second \
output argument is requested, it also returns the corresponding values of the \
neural networks output in @var{scores}. \
\n\
\n\
Installation Note: in order to support parallel processing on MacOS, users \
have to manually add support for OpenMP by adding the following flags to \
@qcode{CFLAGS} and @qcode{CXXFLAGS} prior to installing the statistics \
package:\n\n\
@code{setenv (\"CPPFLAGS\", \"-I/opt/homebrew/opt/libomp/include -Xclang -fopenmp\")} \
\n\
\n\
@seealso{fcnntrain, fitcnet, ClassificationNeuralNetwork} \n\
@end deftypefn")
{
  // Check for correct number of input/output arguments
  if (args.length () < 5)
  {
    error ("fcnnpredict: too few input arguments.");
  }
  if (nargout > 2)
  {
    error ("fcnnpredict: too many output arguments.");
  }

  // The weights and the biases arrive as the classes hold them, one cell per
  // layer each, rather than as the augmented [W b] matrices fcnntrain emits.
  if (! args(0).iscell () ||
      ! (args(0).rows () == 1 && args(0).columns () > 1))
  {
    error ("fcnnpredict: 'LayerWeights' must be a cell row vector.");
  }
  Cell LayerWeights = args(0).cell_value ();
  if (! args(1).iscell () ||
      ! (args(1).rows () == 1 && args(1).columns () > 1))
  {
    error ("fcnnpredict: 'LayerBiases' must be a cell row vector.");
  }
  Cell LayerBiases = args(1).cell_value ();
  if (LayerBiases.numel () != LayerWeights.numel ())
  {
    error ("fcnnpredict: 'LayerBiases' must match 'LayerWeights'.");
  }
  // The layer count comes from the weights, so the names alone say what each
  // layer computes: the hidden layers from 'Activations' and the last from
  // 'OutputLayerActivation'.
  RowVector ActiveCode = activation_codes (args(2), args(3),
                                           LayerWeights.numel () - 1,
                                           "fcnnpredict");

  // Do some input validation while loading the testing data
  if (! args(4).isnumeric () || args(4).iscomplex ())
  {
    error ("fcnnpredict: XC must be a real numeric matrix.");
  }
  if (args(4).isempty ())
  {
    error ("fcnnpredict: XC cannot be empty.");
  }
  if (args(4).columns () != LayerWeights.elem(0).columns ())
  {
    error ("fcnnpredict: the features in XC do not match the trained model.");
  }
  Matrix X = args(4).matrix_value ();
  int n = args(4).rows ();
  int d = args(4).columns ();

  // Check for optional sixth input argument to set number of threads
  int NumThreads = 1;
  if (args.length () == 6)
  {
    if (! args(5).is_scalar_type () || ! args(5).isnumeric () ||
        args(5).scalar_value () < 1 || args(5).iscomplex ())
    {
      error ("fcnnpredict: NumThreads must be a positive integer scalar value.");
    }
    NumThreads = args(5).scalar_value ();
  }

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
    Matrix W = LayerWeights.elem(i).matrix_value ();
    Matrix B = LayerBiases.elem(i).matrix_value ();
    output_size = (int) W.rows ();
    if (B.numel () != output_size)
    {
      error ("fcnnpredict: 'LayerBiases' must match 'LayerWeights'.");
    }
    // set_layer wants each neuron as [weights, bias], which is how the layer
    // packs itself; the two cells are joined a row at a time here rather than
    // as a matrix in the m-code.
    vector<vector<double>> Wb_matrix;
    for (int r = 0; r < W.rows (); r++)
    {
      vector<double> WB_vector;
      for (int c = 0; c < W.columns (); c++)
      {
        WB_vector.push_back (W(r,c));
      }
      WB_vector.push_back (B(r));
      Wb_matrix.push_back (WB_vector);
    }
    // Create dense layer and set its values
    DenseLayer DL = DenseLayer (input_size, output_size);
    DL.set_layer (Wb_matrix);
    // Create activation layer
    ActivationLayer AL = ActivationLayer (ActiveCode(i), NumThreads);
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
%!shared X, Y, MODEL, W, B
%! load fisheriris
%! X = meas;
%! Y = grp2idx (species);
%! MODEL = fcnntrain (X, Y, 10, "sigmoid", "sigmoid", 1, 0.025, 100, false);
%! W = cellfun (@(m) m(:,1:end-1), MODEL.LayerWeights, "UniformOutput", false);
%! B = cellfun (@(m) m(:,end), MODEL.LayerWeights, "UniformOutput", false);
%!test
%! [Y_pred, Y_scores] = fcnnpredict (W, B, "sigmoid", "sigmoid", X);
%! assert_equal (numel (Y_pred), numel (Y));
%! assert_equal (isequal (size (Y_pred), size (Y)), true);
%! assert_equal (columns (Y_scores), numel (unique (Y)));
%! assert_equal (rows (Y_scores), numel (Y));

## A trained network drives its outputs to the targets, not to a fraction of
## them.  A gradient of 2*y-t rather than 2*(y-t) settles at y = t/2, which
## leaves every label right and every score halved, so the scores are what
## has to be checked.
%!test
%! rand ("seed", 42);
%! randn ("seed", 42);
%! Xs = [randn(40,2)*0.3 + 3; randn(40,2)*0.3 - 3];
%! Ys = [ones(40,1); 2*ones(40,1)];
%! M = fcnntrain (Xs, Ys, [8, 8], "sigmoid", "sigmoid", 1, 0.05, 400, false);
%! Wm = cellfun (@(m) m(:,1:end-1), M.LayerWeights, "UniformOutput", false);
%! Bm = cellfun (@(m) m(:,end), M.LayerWeights, "UniformOutput", false);
%! [pred, scores] = fcnnpredict (Wm, Bm, "sigmoid", "sigmoid", [3, 3; -3, -3]);
%! assert_equal (pred, [1; 2]);
%! assert_equal (max (scores(1,:)) > 0.8, true);
%! assert_equal (max (scores(2,:)) > 0.8, true);
%! assert_equal (all (abs (sum (scores, 2) - 1) < 0.1), true);

## Cross entropy with a softmax output trains at least as confidently as the
## mean squared error does, and its rows are a probability by construction.
%!test
%! rand ("seed", 42);
%! randn ("seed", 42);
%! Xs = [randn(40,2)*0.3 + 3; randn(40,2)*0.3 - 3];
%! Ys = [ones(40,1); 2*ones(40,1)];
%! Mm = fcnntrain (Xs, Ys, [8, 8], "sigmoid", "softmax", 1, 0.05, 400, false, 0);
%! Mc = fcnntrain (Xs, Ys, [8, 8], "sigmoid", "softmax", 1, 0.05, 400, false, 1);
%! Wm = cellfun (@(m) m(:,1:end-1), Mm.LayerWeights, "UniformOutput", false);
%! Bm = cellfun (@(m) m(:,end), Mm.LayerWeights, "UniformOutput", false);
%! Wc = cellfun (@(m) m(:,1:end-1), Mc.LayerWeights, "UniformOutput", false);
%! Bc = cellfun (@(m) m(:,end), Mc.LayerWeights, "UniformOutput", false);
%! [pm, sm] = fcnnpredict (Wm, Bm, "sigmoid", "softmax", [3, 3; -3, -3]);
%! [pc, sc] = fcnnpredict (Wc, Bc, "sigmoid", "softmax", [3, 3; -3, -3]);
%! assert_equal (pm, [1; 2]);
%! assert_equal (pc, [1; 2]);
%! assert_equal (all (abs (sum (sc, 2) - 1) < 1e-8), true);
%! assert_equal (min (max (sc, [], 2)) >= min (max (sm, [], 2)), true);

## Omitting the loss selector keeps the mean squared error.
%!test
%! rand ("seed", 42);
%! randn ("seed", 42);
%! Xs = [randn(20,2)*0.3 + 3; randn(20,2)*0.3 - 3];
%! Ys = [ones(20,1); 2*ones(20,1)];
%! rand ("seed", 1);
%! M9 = fcnntrain (Xs, Ys, 6, "sigmoid", "sigmoid", 1, 0.05, 100, false);
%! rand ("seed", 1);
%! M0 = fcnntrain (Xs, Ys, 6, "sigmoid", "sigmoid", 1, 0.05, 100, false, 0);
%! assert_equal (M9.Loss, M0.Loss, 0);

%!error <fcnnpredict: too few input arguments.> ...
%! fcnnpredict (W, B, "sigmoid", "sigmoid");
%!error <fcnnpredict: too many output arguments.> ...
%! [Q, E, R] = fcnnpredict (W, B, "sigmoid", "sigmoid", X);
%!error <fcnnpredict: 'LayerWeights' must be a cell row vector.> ...
%! fcnnpredict (1, B, "sigmoid", "sigmoid", X);
%!error <fcnnpredict: 'LayerWeights' must be a cell row vector.> ...
%! fcnnpredict ({1}, B, "sigmoid", "sigmoid", X);
%!error <fcnnpredict: 'LayerWeights' must be a cell row vector.> ...
%! fcnnpredict ({1; 2; 3}, B, "sigmoid", "sigmoid", X);
%!error <fcnnpredict: 'LayerBiases' must be a cell row vector.> ...
%! fcnnpredict (W, 1, "sigmoid", "sigmoid", X);
%!error <fcnnpredict: 'LayerBiases' must be a cell row vector.> ...
%! fcnnpredict (W, {1}, "sigmoid", "sigmoid", X);
%!error <fcnnpredict: 'LayerBiases' must match 'LayerWeights'.> ...
%! fcnnpredict (W, [B, B], "sigmoid", "sigmoid", X);
%!error <fcnnpredict: 'LayerBiases' must match 'LayerWeights'.> ...
%! fcnnpredict (W, {B{1}(1:end-1), B{2}}, "sigmoid", "sigmoid", X);
%!error <fcnnpredict: 'Activations' must be a character vector or a cellstring.> ...
%! fcnnpredict (W, B, 2, "sigmoid", X);
%!error <fcnnpredict: 'Activations' must be a character vector or a cellstring.> ...
%! fcnnpredict (W, B, {2, 2}, "sigmoid", X);
%!error <fcnnpredict: 'Activations' does not match the number of layers.> ...
%! fcnnpredict (W, B, {"sigmoid", "relu"}, "sigmoid", X);
%!error <fcnnpredict: unsupported 'Activations' function: 'sgmoid'.> ...
%! fcnnpredict (W, B, "sgmoid", "sigmoid", X);
%!error <fcnnpredict: 'OutputLayerActivation' must be a character vector.> ...
%! fcnnpredict (W, B, "sigmoid", 4, X);
%!error <fcnnpredict: unsupported 'OutputLayerActivation' function: 'softmx'.> ...
%! fcnnpredict (W, B, "sigmoid", "softmx", X);
%!error <fcnnpredict: XC must be a real numeric matrix.> ...
%! fcnnpredict (W, B, "sigmoid", "sigmoid", complex (X));
%!error <fcnnpredict: XC must be a real numeric matrix.> ...
%! fcnnpredict (W, B, "sigmoid", "sigmoid", {1, 2, 3, 4});
%!error <fcnnpredict: XC must be a real numeric matrix.> ...
%! fcnnpredict (W, B, "sigmoid", "sigmoid", "asd");
%!error <fcnnpredict: XC cannot be empty.> ...
%! fcnnpredict (W, B, "sigmoid", "sigmoid", []);
%!error <fcnnpredict: the features in XC do not match the trained model.> ...
%! fcnnpredict (W, B, "sigmoid", "sigmoid", X(:,[1:3]));
%!error <fcnnpredict: NumThreads must be a positive integer scalar value.> ...
%! fcnnpredict (W, B, "sigmoid", "sigmoid", X, 0);
*/

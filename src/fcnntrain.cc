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

DEFUN_DLD(fcnntrain, args, nargout,
          "-*- texinfo -*-\n\
 @deftypefn  {statistics} {@var{Mdl} =} fcnntrain (@var{X}, @var{Y}, @\
 @var{LayerSizes}, @var{Activations}, @var{NumThreads}, @var{Alpha}, @\
 @var{LearningRate}, @var{Epochs}, @var{DisplayInfo})\n\
\n\
\n\
Train a fully connected Neural Network. \n\
\n\n\
@code{@var{Mdl} = fcnntrain (@dots{})} requires the following input arguments.\
\n\n\
@itemize \n\
@item @var{X} : An @math{NxM} matrix containing the data set to be trained \
upon.  Rows @math{N} correspond to individual samples and columns @math{M} \
correspond to features (dimensions).  Type of @var{X} must be double. \
\n\
\n\
@item @var{Y} : An @math{Nx1} column vector containing the labels of the \
training dataset.  The labels must be natural numbers (positive integers) \
starting from 1 up to the number of classes, similarily as returned by the \
`grp2idx` function. Type of @var{Y} must be double. \
\n\
\n\
@item @var{LayerSizes} : A numeric row vector of integer values defining the \
size of the hidden layers of the network.  Input and output layers are \
automatically determined by the training data and their labels. \
\n\
\n\
@item @var{Activations} : A numeric row vector of integer values defining the \
activation functions to be used at each layer including the output layer.  The \
corresponding codes to activation functions is: \n\
@itemize \n\
@item @code{0} : @qcode{'Linear'} \n\
@item @code{1} : @qcode{'Sigmoid'} \n\
@item @code{2} : @qcode{'Rectified Linear Unit (ReLU)'} \n\
@item @code{3} : @qcode{'Hyperbolic tangent (tanh)'} \n\
@item @code{4} : @qcode{'Softmax'} \n\
@item @code{5} : @qcode{'Parametric or Leaky ReLU'} \n\
@item @code{6} : @qcode{'Exponential Linear Unit (ELU)'} \n\
@item @code{7} : @qcode{'Gaussian Error Linear Unit (GELU)'} \n\
@end itemize \n\
\n\
\n\
@item @var{NumThreads} : A positive scalar integer value defining the number \
of threads used for computing the activation layers.  For layers with less \
than 1000 neurons, @var{NumThreads} always defaults to 1. \
\n\
\n\
@item @var{Alpha} : A positive scalar value defining the parameter \
@qcode{alpha} used in @qcode{ReLU} and @qcode{ELU} activation layers. \
\n\
\n\
@item @var{LearningRate} : A positive scalar value defining the learning rate \
used by the gradient descend algorithm during training. \
\n\
\n\
@item @var{Epochs} : A positive scalar value defining the number of epochs for \
training the model. \
\n\
\n\
@item @var{DisplayInfo} : A boolean scalar indicating whether to print \
information during training. \n\
@end itemize \n\
\n\
\n\
@code{fcnntrain} returns the trained model, @var{Mdl}, as a structure \
containing the following fields: \
\n\
\n\
@itemize \n\
@item @code{LayerWeights} : A cell array with each element containing a matrix \
with the Weights and Biases of each layer including the output layer.\n\
\n\
\n\
@item @code{Activations} : A numeric row vector of integer values defining the \
activation functions to be used at each layer including the output layer. \
\n\
\n\
@item @code{Accuracy} : The prediction accuracy at each iteration during the \
neural network model's training process. \
\n\
\n\
@item @code{Loss} : The loss value recorded at each iteration during the \
neural network model's training process. \
\n\
\n\
@item @code{Alpha} : The value of the Alpha parameter used in @qcode{ReLU} and \
@qcode{ELU} activation layers. \
\n\
\n\
@end itemize \
\n\
\n\
Installation Note: in order to support parallel processing on MacOS, users \
have to manually add support for OpenMP by adding the following flags to \
@qcode{CFLAGS} and @qcode{CXXFLAGS} prior to installing the statistics \
package:\n\n\
@code{setenv (\"CPPFLAGS\", \"-I/opt/homebrew/opt/libomp/include -Xclang -fopenmp\")} \
\n\
\n\
@seealso{fcnnpredict, fitcnet, ClassificationNeuralNetwork} \n\
@end deftypefn")
{
  // Check for correct number of input/output arguments
  if (args.length () < 9)
  {
    error ("fcnntrain: too few input arguments.");
  }
  if (nargout > 1)
  {
    error ("fcnntrain: too many output arguments.");
  }
  // int seed = time(NULL);
  int seed = 0;
  srand (seed);

  // Do some input validation while loading trainind data and labels
  if (! args(0).isnumeric () || args(0).iscomplex ())
  {
    error ("fcnntrain: X must be a real numeric matrix.");
  }
  if (args(0).isempty ())
  {
    error ("fcnntrain: X cannot be empty.");
  }
  if (! args(1).isnumeric () || args(1).iscomplex ())
  {
    error ("fcnntrain: Y must be a real numeric matrix.");
  }
  if (args(1).isempty ())
  {
    error ("fcnntrain: Y cannot be empty.");
  }
  if (args(0).rows () != args(1).rows ())
  {
    error ("fcnntrain: X and Y must have the same number of rows.");
  }

  // Construct 2D vector from data in X
  int n = args(0).rows ();
  int d = args(0).columns ();
  vector<vector<double>> data (n, vector<double>(d, 0));
  Matrix X = args(0).matrix_value ();
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < d; j++)
    {
      data[i][j] = X(i,j);
    }
  }

  // Construct 1D vector from labels in Y
  vector<int> labels(n, 0);
  ColumnVector Y = args(1).column_vector_value ();
  for (int i = 0; i < n; i++)
  {
    labels[i] = Y(i);
    if (labels[i] < 1)
    {
      error ("fcnntrain: labels in Y must be positive integers.");
    }
  }

  // Check LayerSizes and Activations input arguments
  if (! args(2).isnumeric () || args(2).iscomplex () || args(2).isempty () ||
        args(2).rows () != 1)
  {
    error ("fcnntrain: 'LayerSizes' must be a row vector of integer values.");
  }
  if (! args(3).isnumeric () || args(3).iscomplex () || args(3).isempty () ||
        args(3).rows () != 1 || args(3).columns () < 2)
  {
    error ("fcnntrain: 'Activations' must be a row vector of integer values.");
  }
  if (args(2).numel () != args(3).numel () - 1)
  {
    error ("fcnntrain: 'Activations' do not match LayerSizes.");
  }

  // Check NumThreads and Alpha input arguments
  if (! args(4).is_scalar_type () || ! args(4).isnumeric () ||
        args(4).scalar_value () < 1 || args(4).iscomplex ())
  {
    error ("fcnntrain: 'NumThreads' must be a positive integer scalar value.");
  }
  int NumThreads = args(4).scalar_value ();
  if (! args(5).is_scalar_type () || ! args(5).isnumeric () ||
        args(5).scalar_value () <= 0 || args(5).iscomplex ())
  {
    error ("fcnntrain: 'Alpha' must be a positive scalar value.");
  }
  double Alpha = args(5).scalar_value ();

  // Create a vector of layers sized appropriately
  vector<DenseLayer> WeightBias;
  vector<ActivationLayer> Activation;
  RowVector LayerSizes = args(2).row_vector_value ();
  RowVector ActiveCode = args(3).row_vector_value ();
  int numlayers = args(2).numel () + 1;
  int input_size = d;
  for (int i = 0; i < args(2).numel (); i++)
  {
    int output_size = (int) LayerSizes(i);
    if (output_size < 1)
    {
      error ("fcnntrain: cannot have a layer of zero size.");
    }
    DenseLayer DL = DenseLayer (input_size, output_size);
    int code = ActiveCode(i);
    if (code < 0 || code > 7)
    {
      error ("fcnntrain: invalid 'Activations' code.");
    }
    ActivationLayer AL = ActivationLayer (code, NumThreads, Alpha);
    WeightBias.push_back (DL);
    Activation.push_back (AL);
    input_size = output_size;
  }
  // Push back last dense layer
  int output_size = set<int> (labels.begin (), labels.end ()).size ();
  DenseLayer DL = DenseLayer (input_size, output_size);
  int last_AC = args(2).numel ();
  ActivationLayer AL = ActivationLayer (ActiveCode(last_AC), NumThreads, Alpha);
  WeightBias.push_back (DL);
  Activation.push_back (AL);

  // Input validation on LearningRate, Epochs, and DisplayInfo
  if (! args(6).is_scalar_type () || ! args(6).isnumeric ())
  {
    error ("fcnntrain: 'LearningRate' must be a positive scalar value.");
  }
  double learning_rate = args(6).scalar_value ();
  if (learning_rate <= 0)
  {
    error ("fcnntrain: 'LearningRate' must be a positive scalar value.");
  }
  if (! args(7).is_scalar_type () || ! args(7).isnumeric ())
  {
    error ("fcnntrain: 'Epochs' must be a positive scalar value.");
  }
  if (args(7).scalar_value () < 1)
  {
    error ("fcnntrain: 'Epochs' must be a positive scalar value.");
  }
  if (! args(8).is_bool_scalar ())
  {
    error ("fcnntrain: 'DisplayInfo' must be a boolean scalar.");
  }

  // Initialize return variables
  octave_idx_type max_epochs = args(7).idx_type_value ();
  vector<double> Accuracy (max_epochs);
  vector<double> Loss (max_epochs);

  // Start training
  octave_idx_type epoch = 0;
  for (; epoch < max_epochs; epoch++)
  {
    // Initialize Loss and Prediction
    double sum_loss = 0.0;
    vector<int> predictions = vector<int> ();

    // Go through all training samples
    for (int sample_idx = 0; sample_idx < n; sample_idx++)
    {
      vector<double> sample = data[sample_idx];
      int label = labels[sample_idx];

      // Forward pass
      for (int layer_idx = 0; layer_idx < numlayers; layer_idx++)
      {
        sample = WeightBias[layer_idx].forward (sample);
        sample = Activation[layer_idx].forward (sample);
      }

      // Get the prediction for this sample
      vector<double> label_vector = vector<double> (output_size);
      label_vector[label-1] = 1.0;  // Labels in Y start from 1
      int prediction = 0;           // Search for highest value
      for (int j = 0; j < output_size; j++)
      {
        if (sample[j] > sample[prediction])
        {
          prediction = j;
        }
      }
      predictions.push_back (prediction);

      // Compute loss
      MeanSquaredErrorLoss loss = MeanSquaredErrorLoss ();
      double loss_output = loss.forward (sample, label_vector);
      sum_loss += loss_output;

      // Print output
      if (args(8).scalar_value () != 0)
      {
        if (sample_idx % 500 == 0)
        {
          cout << setprecision(4) << "i:" << sample_idx << " | Mean Loss: ";
          cout << (sum_loss / (sample_idx + 1)) << "\r" << flush;
        }
      }

      // Backward pass
      for (int layer_idx = 0; layer_idx < numlayers; layer_idx++)
      {
        WeightBias[layer_idx].zero_gradient (); // Reset gradients to zero
      }

      loss.backward(1.0);

      // Compute gradients
      Activation[numlayers-1].backward(loss.grad);

      for (int layer_idx = numlayers; layer_idx > 0; layer_idx--)
      {
        WeightBias[layer_idx-1].backward (Activation[layer_idx-1].grad);
        if (layer_idx > 1)
        {
          Activation[layer_idx-2].backward (WeightBias[layer_idx-1]);
        }
      }

      // Update weights
      for (int layer_idx = 0; layer_idx < numlayers; layer_idx++)
      {
        WeightBias[layer_idx].descend (learning_rate);
      }
    }

    // Compute Loss and Accuracy on the training set
    double A = accuracy (predictions, labels);
    double L = sum_loss / n;
    Accuracy.push_back (A);
    Loss.push_back (L);

    // Print output
    if (args(8).scalar_value () != 0)
    {
      cout << "                                              \r"
           << "Epoch: " << epoch + 1 << " | Loss: "
           << L << " | Train Accuracy: " << A << endl;
    }
  }

  // Get weights and biases from each layer and store them in a cell array
  Cell LayerWeights(1, numlayers);
  for (int layer_idx = 0; layer_idx < numlayers; layer_idx++)
  {
    DenseLayer DL = WeightBias[layer_idx];
    vector<vector<double>> Wb_matrix = DL.get_layer ();
    vector<double> WB_vector = Wb_matrix[0];
    octave_idx_type row = Wb_matrix.size ();
    octave_idx_type col = WB_vector.size ();
    Matrix WB (row, col);
    for (int r = 0; r < row; r++)
    {
      for (int c = 0; c < col; c++)
      {
        WB(r,c) = Wb_matrix[r][c];
      }
    }
    LayerWeights.elem(layer_idx) = WB;
  }

  // Store accuracy and loss vectors in RowVector
  RowVector A(max_epochs);
  RowVector L(max_epochs);
  for (epoch = 0; epoch < max_epochs; epoch++)
  {
    A(epoch) = Accuracy[epoch];
    L(epoch) = Loss[epoch];
  }

  // Prepare returning arguments
  octave_scalar_map fcnn_model;
  fcnn_model.assign ("LayerWeights", LayerWeights);
  fcnn_model.assign ("Activations", ActiveCode);
  fcnn_model.assign ("Accuracy", A);
  fcnn_model.assign ("Loss", L);
  fcnn_model.assign ("Alpha", Alpha);
  octave_value_list retval (1);
  retval(0) = fcnn_model;
  return retval;
}

/*
%!shared X, Y, MODEL
%! load fisheriris
%! X = meas;
%! Y = grp2idx (species);
%!error <fcnntrain: too few input arguments.> ...
%! model = fcnntrain (X, Y);
%!error <fcnntrain: too many output arguments.> ...
%! [Q, W] = fcnntrain (X, Y, 10, [1, 1], 1, 0.01, 0.025, 50, false);
%!error <fcnntrain: X must be a real numeric matrix.> ...
%! fcnntrain (complex (X), Y, 10, [1, 1], 1, 0.01, 0.025, 50, false);
%!error <fcnntrain: X must be a real numeric matrix.> ...
%! fcnntrain ({X}, Y, 10, [1, 1], 1, 0.01, 0.025, 50, false);
%!error <fcnntrain: X cannot be empty.> ...
%! fcnntrain ([], Y, 10, [1, 1], 1, 0.01, 0.025, 50, false);
%!error <fcnntrain: Y must be a real numeric matrix.> ...
%! fcnntrain (X, complex (Y), 10, 1, 0.01, [1, 1], 0.025, 50, false);
%!error <fcnntrain: Y must be a real numeric matrix.> ...
%! fcnntrain (X, {Y}, 10, [1, 1], 1, 0.01, 0.025, 50, false);
%!error <fcnntrain: Y cannot be empty.> ...
%! fcnntrain (X, [], 10, [1, 1], 1, 0.01, 0.025, 50, false);
%!error <fcnntrain: X and Y must have the same number of rows.> ...
%! fcnntrain (X, Y([1:50]), 10, [1, 1], 1, 0.01, 0.025, 50, false);
%!error <fcnntrain: labels in Y must be positive integers.> ...
%! fcnntrain (X, Y - 1, 10, [1, 1], 1, 0.01, 0.025, 50, false);
%!error <fcnntrain: 'LayerSizes' must be a row vector of integer values.> ...
%! fcnntrain (X, Y, [10; 5], [1, 1, 1], 1, 0.01, 0.025, 50, false);
%!error <fcnntrain: 'LayerSizes' must be a row vector of integer values.> ...
%! fcnntrain (X, Y, "10", [1, 1], 1, 0.01, 0.025, 50, false);
%!error <fcnntrain: 'LayerSizes' must be a row vector of integer values.> ...
%! fcnntrain (X, Y, {10}, [1, 1], 1, 0.01, 0.025, 50, false);
%!error <fcnntrain: 'LayerSizes' must be a row vector of integer values.> ...
%! fcnntrain (X, Y, complex (10), [1, 1], 1, 0.01, 0.025, 50, false);
%!error <fcnntrain: 'Activations' must be a row vector of integer values.> ...
%! fcnntrain (X, Y, 10, [1; 1], 1, 0.01, 0.025, 50, false);
%!error <fcnntrain: 'Activations' must be a row vector of integer values.> ...
%! fcnntrain (X, Y, 10, {1, 1}, 1, 0.01, 0.025, 50, false);
%!error <fcnntrain: 'Activations' must be a row vector of integer values.> ...
%! fcnntrain (X, Y, 10, "1", 1, 0.01, 0.025, 50, false);
%!error <fcnntrain: 'Activations' must be a row vector of integer values.> ...
%! fcnntrain (X, Y, 10, complex ([1, 1]), 1, 0.01, 0.025, 50, false);
%!error <fcnntrain: 'Activations' do not match LayerSizes.> ...
%! fcnntrain (X, Y, 10, [1, 1, 1], 1, 0.01, 0.025, 50, false);
%!error <fcnntrain: cannot have a layer of zero size.> ...
%! fcnntrain (X, Y, [10, 0, 5], [1, 1, 1, 1], 1, 0.01, 0.025, 50, false);
%!error <fcnntrain: invalid 'Activations' code.> ...
%! fcnntrain (X, Y, 10, [-1, 1], 1, 0.01, 0.025, 50, false);
%!error <fcnntrain: invalid 'Activations' code.> ...
%! fcnntrain (X, Y, 10, [8, 1], 1, 0.01, 0.025, 50, false);
%!error <fcnntrain: 'NumThreads' must be a positive integer scalar value.> ...
%! fcnntrain (X, Y, 10, [1, 1], 0, 0.01, 0.025, 50, false);
%!error <fcnntrain: 'Alpha' must be a positive scalar value.> ...
%! fcnntrain (X, Y, 10, [1, 1], 1, -0.01, 0.025, 50, false);
%!error <fcnntrain: 'LearningRate' must be a positive scalar value.> ...
%! fcnntrain (X, Y, 10, [1, 1], 1, 0.01, -0.025, 50, false);
%!error <fcnntrain: 'LearningRate' must be a positive scalar value.> ...
%! fcnntrain (X, Y, 10, [1, 1], 1, 0.01, 0, 50, false);
%!error <fcnntrain: 'LearningRate' must be a positive scalar value.> ...
%! fcnntrain (X, Y, 10, [1, 1], 1, 0.01, [0.025, 0.001], 50, false);
%!error <fcnntrain: 'LearningRate' must be a positive scalar value.> ...
%! fcnntrain (X, Y, 10, [1, 1], 1, 0.01, {0.025}, 50, false);
%!error <fcnntrain: 'Epochs' must be a positive scalar value.> ...
%! fcnntrain (X, Y, 10, [1, 1], 1, 0.01, 0.025, 0, false);
%!error <fcnntrain: 'Epochs' must be a positive scalar value.> ...
%! fcnntrain (X, Y, 10, [1, 1], 1, 0.01, 0.025, [50, 25], false);
%!error <fcnntrain: 'DisplayInfo' must be a boolean scalar.> ...
%! fcnntrain (X, Y, 10, [1, 1], 1, 0.01, 0.025, 50, 0);
%!error <fcnntrain: 'DisplayInfo' must be a boolean scalar.> ...
%! fcnntrain (X, Y, 10, [1, 1], 1, 0.01, 0.025, 50, 1);
%!error <fcnntrain: 'DisplayInfo' must be a boolean scalar.> ...
%! fcnntrain (X, Y, 10, [1, 1], 1, 0.01, 0.025, 50, [false, false]);
*/

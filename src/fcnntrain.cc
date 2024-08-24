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
 @var{LayerSizes}, @var{Activations}, @var{LearningRate}, @var{Epochs}, @\
 @var{DisplayInfo})\n\
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
@end itemize \n\n\
@seealso{fcnnpredict, fitcnet, ClassificationNeuralNetwork} \n\
@end deftypefn")
{
  // int seed = time(NULL);
  int seed = 0;
  srand (seed);
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
  }
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
    DenseLayer DL = DenseLayer (input_size, output_size);
    ActivationLayer AL = ActivationLayer (ActiveCode(i), 0.01);
    WeightBias.push_back (DL);
    Activation.push_back (AL);
    input_size = output_size;
  }
  // Push back last dense layer
  int output_size = set<int> (labels.begin (), labels.end ()).size ();
  DenseLayer DL = DenseLayer (input_size, output_size);
  int last_AC = args(2).numel ();
  ActivationLayer AL = ActivationLayer (ActiveCode(last_AC), 0.01);
  WeightBias.push_back (DL);
  Activation.push_back (AL);

  vector<double> Accuracy;
  vector<double> Loss;
  double learning_rate = args(4).scalar_value ();

  // Start training
  int epoch = 0;
  for (; epoch < args(5).uint_value (); epoch++)
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
      if (args(6).scalar_value () != 0)
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
    if (args(6).scalar_value () != 0)
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

  // Store loss vector in RowVector
  RowVector A(n);
  RowVector L(n);
  for (int sample_idx = 0; sample_idx < n; sample_idx++)
  {
    A(sample_idx) = Accuracy[sample_idx];
    L(sample_idx) = Loss[sample_idx];
  }

  // Prepare returning arguments
  octave_scalar_map fcnn_model;
  fcnn_model.assign ("LayerWeights", LayerWeights);
  fcnn_model.assign ("Activations", ActiveCode);
  fcnn_model.assign ("Accuracy", A);
  fcnn_model.assign ("Loss", L);
  octave_value_list retval (1);
  retval(0) = fcnn_model;
  return retval;
}


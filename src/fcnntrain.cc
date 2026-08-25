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
#include "lbfgs.h"

using namespace std;

// fitcnet's wording for each stopping criterion, measured on R2024a.  It
// lives here and not in lbfgs.h because MATLAB words the same criterion
// differently per function: fitclinear reports "Tolerance on gradient
// satisfied." where this reports "Relative gradient tolerance reached.",
// so the engine hands out a token and every caller supplies its own prose.
static const char *
criterion_message (int crit)
{
  switch (crit)
  {
    case lbfgs::CRIT_GRADIENT:
      return "Relative gradient tolerance reached.";
    case lbfgs::CRIT_STEP:
      return "Step size tolerance reached.";
    case lbfgs::CRIT_LOSS:
      return "Loss tolerance reached.";
    case lbfgs::CRIT_ITERATION_LIMIT:
      return "Iteration limit reached.";
    default:
      return "Line search could not improve the objective.";
  }
}

// Full-batch objective: the mean loss over every training sample together
// with its gradient against the flat parameter vector.  One call is one sweep
// forward and back over the whole training set, which is what each trial step
// of a line search costs.
//
// This is the piece the epoch loop cannot supply.  That loop clears the
// gradient before every sample and consumes it immediately, so no gradient of
// the summed loss ever exists; here it is cleared once and accumulated across
// the batch.
class fcnn_objective
{
public:

  fcnn_objective (vector<DenseLayer>& wb, vector<ActivationLayer>& act,
                  const vector<vector<double>>& data,
                  const vector<vector<double>>& targets,
                  const vector<int>& labels, int output_size, int lossfcn,
                  bool regression)
    : m_wb (wb), m_act (act), m_data (data), m_targets (targets),
      m_labels (labels), m_output_size (output_size), m_lossfcn (lossfcn),
      m_regression (regression)
  { }

  double operator () (const vector<double>& w, vector<double>& g)
  {
    octave_quit ();

    const int numlayers = m_wb.size ();
    const int n = m_data.size ();

    const double *p = w.data ();
    for (int l = 0; l < numlayers; l++)
    {
      m_wb[l].unpack (p);
      p += m_wb[l].nparams ();
    }

    for (int l = 0; l < numlayers; l++)
    {
      m_wb[l].zero_gradient ();
    }

    double total = 0.0;
    for (int s = 0; s < n; s++)
    {
      vector<double> sample = m_data[s];
      for (int l = 0; l < numlayers; l++)
      {
        sample = m_wb[l].forward (sample);
        sample = m_act[l].forward (sample);
      }

      vector<double> label_vector;
      if (m_regression)
      {
        label_vector = m_targets[s];
      }
      else
      {
        label_vector = vector<double> (m_output_size, 0.0);
        label_vector[m_labels[s]-1] = 1.0;   // labels in Y start from 1
      }

      vector<double> loss_grad;
      if (m_lossfcn == 1)
      {
        CrossEntropyLoss loss = CrossEntropyLoss ();
        total += loss.forward (sample, label_vector);
        // The objective is the mean, so each sample's gradient carries 1/n.
        loss.backward (1.0 / n);
        loss_grad = loss.grad;
      }
      else
      {
        MeanSquaredErrorLoss loss = MeanSquaredErrorLoss ();
        total += loss.forward (sample, label_vector);
        loss.backward (1.0 / n);
        loss_grad = loss.grad;
      }

      m_act[numlayers-1].backward (loss_grad);
      for (int l = numlayers; l > 0; l--)
      {
        m_wb[l-1].backward (m_act[l-1].grad);
        if (l > 1)
        {
          m_act[l-2].backward (m_wb[l-1]);
        }
      }
    }

    double *q = &g[0];
    for (int l = 0; l < numlayers; l++)
    {
      m_wb[l].pack_grad (q);
      q += m_wb[l].nparams ();
    }

    return total / n;
  }

private:

  vector<DenseLayer>& m_wb;
  vector<ActivationLayer>& m_act;
  const vector<vector<double>>& m_data;
  const vector<vector<double>>& m_targets;
  const vector<int>& m_labels;
  int m_output_size;
  int m_lossfcn;
  bool m_regression;
};

DEFUN_DLD(fcnntrain, args, nargout,
          "-*- texinfo -*-\n\
 @deftypefn  {statistics} {@var{Mdl} =} fcnntrain (@var{X}, @var{Y}, @\
 @var{LayerSizes}, @var{Activations}, @var{OutputLayerActivation}, @\
 @var{NumThreads}, @var{LearningRate}, @var{Epochs}, @var{DisplayInfo})\n\
 @deftypefnx  {statistics} {@var{Mdl} =} fcnntrain (@dots{}, @\
 @var{LossFunction})\n\
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
starting from 1 up to the number of classes, similarly as returned by the \
`grp2idx` function. Type of @var{Y} must be double.  Under regression, \
selected by @var{LossFunction} 2, @var{Y} is instead an @math{NxR} matrix of \
response values, which may take any finite value, and the output layer is \
sized to its @math{R} columns rather than to a number of classes. \
\n\
\n\
@item @var{LayerSizes} : A numeric row vector of integer values defining the \
size of the hidden layers of the network.  Input and output layers are \
automatically determined by the training data and their labels. \
\n\
\n\
@item @var{Activations} : The activation function of the hidden layers, named \
as a character vector applying to all of them or as a cellstring naming them \
one by one, in which case it must have one name per hidden layer.  The \
supported names are: \n\
@itemize \n\
@item @qcode{'linear'} or @qcode{'none'} : the identity \n\
@item @qcode{'sigmoid'} \n\
@item @qcode{'relu'} : rectified linear unit \n\
@item @qcode{'tanh'} : hyperbolic tangent \n\
@item @qcode{'softmax'} \n\
@item @qcode{'lrelu'} or @qcode{'prelu'} : leaky rectified linear unit, whose \
negative slope is a constant 0.01 \n\
@item @qcode{'elu'} : exponential linear unit, saturating at @math{-1} \n\
@item @qcode{'gelu'} : Gaussian error linear unit \n\
@end itemize \n\
\n\
\n\
@item @var{OutputLayerActivation} : The activation function of the output \
layer, named as a character vector from the same list. \
\n\
\n\
@item @var{NumThreads} : A positive scalar integer value defining the number \
of threads used for computing the activation layers.  For layers with less \
than 1000 neurons, @var{NumThreads} always defaults to 1. \
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
@code{@var{Mdl} = fcnntrain (@dots{}, @var{LossFunction})} also selects the \
loss the network is trained against.  @var{LossFunction} is a scalar: 0 for \
mean squared error over a one-hot target, which is the default, 1 for cross \
entropy, and 2 for mean squared error over a continuous response.  Cross \
entropy expects the output layer to report a probability over the classes, \
so it belongs with a softmax output; paired that way the two gradients \
compose to @math{y - t}.  Its loss is undefined where the predicted \
probability of the true class is zero, so both the logarithm and its \
derivative are floored.  Code 2 is regression: @var{Y} holds response values \
rather than labels, the output layer belongs with the identity activation, \
and the returned model carries no @code{Accuracy} field, there being no \
labels to count. \n\
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
neural network model's training process.  Absent under regression. \
\n\
\n\
@item @code{Loss} : The loss value recorded at each iteration during the \
neural network model's training process. \
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
  if (args.length () > 11)
  {
    error ("fcnntrain: too many input arguments.");
  }
  if (nargout > 1)
  {
    error ("fcnntrain: too many output arguments.");
  }

  // Optional tenth argument selecting the loss: 0 for mean squared error over
  // a one-hot target, which is the default and what earlier releases always
  // used, 1 for cross entropy, which expects the output layer to report a
  // probability, and 2 for mean squared error over a continuous response,
  // which is regression.  It is read first because it decides whether Y holds
  // class labels or response values.
  int lossfcn = 0;
  if (args.length () >= 10)
  {
    if (! args(9).isnumeric () || ! args(9).is_scalar_type ()
        || args(9).iscomplex ())
    {
      error ("fcnntrain: 'LossFunction' must be a numeric scalar value.");
    }
    lossfcn = args(9).int_value ();
    if (lossfcn < 0 || lossfcn > 2)
    {
      error ("fcnntrain: invalid 'LossFunction' code.");
    }
  }
  bool regression = (lossfcn == 2);
  // Do some input validation while loading training data and labels
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

  // Construct the training targets from Y.  Under regression Y holds the
  // response itself, one column per response, and is taken as it stands;
  // otherwise it holds class labels, which index a one-hot target built below.
  vector<int> labels;
  vector<vector<double>> targets;
  int n_response = 0;
  if (regression)
  {
    n_response = args(1).columns ();
    Matrix R = args(1).matrix_value ();
    targets = vector<vector<double>> (n, vector<double> (n_response, 0.0));
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n_response; j++)
      {
        if (! octave::math::isfinite (R(i,j)))
        {
          error ("fcnntrain: Y must be finite.");
        }
        targets[i][j] = R(i,j);
      }
    }
  }
  else
  {
    labels = vector<int> (n, 0);
    ColumnVector Y = args(1).column_vector_value ();
    for (int i = 0; i < n; i++)
    {
      labels[i] = Y(i);
      if (labels[i] < 1)
      {
        error ("fcnntrain: labels in Y must be positive integers.");
      }
    }
  }

  // Check LayerSizes and Activations input arguments
  if (! args(2).isnumeric () || args(2).iscomplex () || args(2).isempty () ||
        args(2).rows () != 1)
  {
    error ("fcnntrain: 'LayerSizes' must be a row vector of integer values.");
  }

  // Check the NumThreads input argument
  if (! args(5).is_scalar_type () || ! args(5).isnumeric () ||
        args(5).scalar_value () < 1 || args(5).iscomplex ())
  {
    error ("fcnntrain: 'NumThreads' must be a positive integer scalar value.");
  }
  int NumThreads = args(5).scalar_value ();

  // Create a vector of layers sized appropriately.  The initial weights are
  // drawn from Octave's generator, so rand ('seed', s) reproduces a fit.
  uniform_scope draw_uniform;
  vector<DenseLayer> WeightBias;
  vector<ActivationLayer> Activation;
  RowVector LayerSizes = args(2).row_vector_value ();
  RowVector ActiveCode = activation_codes (args(3), args(4), args(2).numel (),
                                          "fcnntrain");
  int numlayers = args(2).numel () + 1;
  int input_size = d;
  for (int i = 0; i < args(2).numel (); i++)
  {
    int output_size = (int) LayerSizes(i);
    if (output_size < 1)
    {
      error ("fcnntrain: cannot have a layer of zero size.");
    }
    int code = ActiveCode(i);
    DenseLayer DL = DenseLayer (input_size, output_size, code);
    ActivationLayer AL = ActivationLayer (code, NumThreads);
    WeightBias.push_back (DL);
    Activation.push_back (AL);
    input_size = output_size;
  }
  // Push back last dense layer
  int output_size = regression ? n_response
                    : (int) set<int> (labels.begin (), labels.end ()).size ();
  int last_AC = args(2).numel ();
  DenseLayer DL = DenseLayer (input_size, output_size,
                              (int) ActiveCode(last_AC));
  ActivationLayer AL = ActivationLayer (ActiveCode(last_AC), NumThreads);
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
  // Reserved, not sized: these are filled with push_back below, and sizing
  // them here would leave max_epochs zeros in front of the values and the
  // reported history reading back as all zero.
  // Optional eleventh argument: the solver and its tolerances, as a scalar
  // struct.  Absent, the epoch loop below runs exactly as it always has, and
  // every default in this file is unchanged.
  bool use_lbfgs = false;
  lbfgs::options lbopt;
  lbopt.iteration_limit = max_epochs;
  lbfgs::result lbres;
  if (args.length () > 10)
  {
    if (! args(10).isstruct () || args(10).numel () != 1)
    {
      error ("fcnntrain: 'SolverOptions' must be a scalar struct.");
    }
    octave_scalar_map so = args(10).scalar_map_value ();
    if (so.isfield ("Solver"))
    {
      std::string sv = so.contents ("Solver").string_value ();
      if (sv == "lbfgs")
      {
        use_lbfgs = true;
      }
      else if (sv != "sgd")
      {
        error ("fcnntrain: 'Solver' must be 'sgd' or 'lbfgs'.");
      }
    }
    if (so.isfield ("GradientTolerance"))
    {
      lbopt.gradient_tolerance
        = so.contents ("GradientTolerance").double_value ();
    }
    if (so.isfield ("LossTolerance"))
    {
      lbopt.loss_tolerance = so.contents ("LossTolerance").double_value ();
    }
    if (so.isfield ("StepTolerance"))
    {
      lbopt.step_tolerance = so.contents ("StepTolerance").double_value ();
    }
    if (so.isfield ("HistorySize"))
    {
      lbopt.history_size = so.contents ("HistorySize").int_value ();
    }
  }
  vector<double> Accuracy;
  vector<double> Loss;
  Accuracy.reserve (max_epochs);
  Loss.reserve (max_epochs);

  // Order the samples are visited in, reshuffled every epoch below.  Visiting
  // them in a fixed order makes every update for one class precede the first
  // update for the next whenever the labels arrive sorted, which is what
  // grp2idx returns, and the weights then swing between the classes instead of
  // settling.
  vector<int> order (n);
  for (int i = 0; i < n; i++)
  {
    order[i] = i;
  }

  // The epoch loop is the stochastic solver.  LBFGS replaces it whole rather
  // than adjusting it: it needs the gradient of the summed loss, which no
  // per-sample update can produce, so the two paths share the network and
  // nothing else.  The sample order is not shuffled here, the fit being
  // deterministic once the weights are drawn.
  if (use_lbfgs)
  {
    int nparams = 0;
    for (int layer_idx = 0; layer_idx < numlayers; layer_idx++)
    {
      nparams += WeightBias[layer_idx].nparams ();
    }

    vector<double> w (nparams, 0.0);
    double *p = &w[0];
    for (int layer_idx = 0; layer_idx < numlayers; layer_idx++)
    {
      WeightBias[layer_idx].pack (p);
      p += WeightBias[layer_idx].nparams ();
    }

    fcnn_objective fobj (WeightBias, Activation, data, targets, labels,
                         output_size, lossfcn, regression);
    lbres = lbfgs::minimize (fobj, w, lbopt);

    // minimize leaves the network holding whatever the last trial step set,
    // which is the accepted point only by construction of the line search;
    // writing the returned vector back makes that explicit.
    const double *q = w.data ();
    for (int layer_idx = 0; layer_idx < numlayers; layer_idx++)
    {
      WeightBias[layer_idx].unpack (q);
      q += WeightBias[layer_idx].nparams ();
    }

    if (args(8).scalar_value () != 0)
    {
      cout << "Iterations: " << lbres.iterations << " | Loss: "
           << lbres.fval << " | " << criterion_message (lbres.crit)
           << endl;
    }
  }
  else
  {
    // Start training
    octave_idx_type epoch = 0;
    for (; epoch < max_epochs; epoch++)
    {
      // Fisher-Yates over Octave's generator, so that a seeded fit stays
      // reproducible
      for (int i = n - 1; i > 0; i--)
      {
        int j = (int) (octave::rand::scalar () * (i + 1));
        if (j > i)          // scalar () is documented on [0, 1); guard the end
        {
          j = i;
        }
        int keep = order[i];
        order[i] = order[j];
        order[j] = keep;
      }

      // Running loss, for the progress line only: the weights move under it, so
      // it is not the loss of any one network and is never recorded.
      double running_loss = 0.0;

      // Go through all training samples
      for (int visit = 0; visit < n; visit++)
      {
        int sample_idx = order[visit];
        vector<double> sample = data[sample_idx];

        // Forward pass
        for (int layer_idx = 0; layer_idx < numlayers; layer_idx++)
        {
          sample = WeightBias[layer_idx].forward (sample);
          sample = Activation[layer_idx].forward (sample);
        }

        vector<double> label_vector;
        if (regression)
        {
          label_vector = targets[sample_idx];
        }
        else
        {
          label_vector = vector<double> (output_size);
          label_vector[labels[sample_idx]-1] = 1.0;  // Labels in Y start from 1
        }

        // Compute loss and the gradient it hands back
        double loss_output;
        vector<double> loss_grad;
        if (lossfcn == 1)
        {
          CrossEntropyLoss loss = CrossEntropyLoss ();
          loss_output = loss.forward (sample, label_vector);
          loss.backward (1.0);
          loss_grad = loss.grad;
        }
        else
        {
          MeanSquaredErrorLoss loss = MeanSquaredErrorLoss ();
          loss_output = loss.forward (sample, label_vector);
          loss.backward (1.0);
          loss_grad = loss.grad;
        }
        running_loss += loss_output;

        // Print output
        if (args(8).scalar_value () != 0)
        {
          if (visit % 500 == 0)
          {
            cout << setprecision(4) << "i:" << visit << " | Mean Loss: ";
            cout << (running_loss / (visit + 1)) << "\r" << flush;
          }
        }

        // Backward pass
        for (int layer_idx = 0; layer_idx < numlayers; layer_idx++)
        {
          WeightBias[layer_idx].zero_gradient (); // Reset gradients to zero
        }

        // Compute gradients
        Activation[numlayers-1].backward (loss_grad);

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

      // Loss and accuracy of the network as it stands at the end of the
      // epoch, measured in one forward-only pass with the weights held
      // still.  Summing them inside the loop above instead scores every
      // sample against different weights, so the figure belongs to no
      // network that ever existed: on class-interleaved data a network stuck
      // on a constant output reports an accuracy of exactly zero, each sample
      // being scored against weights just pulled toward the one before it.
      double sum_loss = 0.0;
      vector<int> predictions = vector<int> ();
      predictions.reserve (n);
      for (int sample_idx = 0; sample_idx < n; sample_idx++)
      {
        vector<double> sample = data[sample_idx];
        for (int layer_idx = 0; layer_idx < numlayers; layer_idx++)
        {
          sample = WeightBias[layer_idx].forward (sample);
          sample = Activation[layer_idx].forward (sample);
        }

        vector<double> label_vector;
        if (regression)
        {
          label_vector = targets[sample_idx];
        }
        else
        {
          int prediction = 0;         // Search for highest value
          for (int j = 0; j < output_size; j++)
          {
            if (sample[j] > sample[prediction])
            {
              prediction = j;
            }
          }
          predictions.push_back (prediction);
          label_vector = vector<double> (output_size);
          label_vector[labels[sample_idx]-1] = 1.0;
        }

        if (lossfcn == 1)
        {
          CrossEntropyLoss loss = CrossEntropyLoss ();
          sum_loss += loss.forward (sample, label_vector);
        }
        else
        {
          MeanSquaredErrorLoss loss = MeanSquaredErrorLoss ();
          sum_loss += loss.forward (sample, label_vector);
        }
      }

      // Accuracy counts correct labels, which regression has none of.
      double A = regression ? 0.0 : accuracy (predictions, labels);
      double L = sum_loss / n;
      Accuracy.push_back (A);
      Loss.push_back (L);

      // Print output
      if (args(8).scalar_value () != 0)
      {
        cout << "                                              \r"
             << "Epoch: " << epoch + 1 << " | Loss: " << L;
        if (! regression)
        {
          cout << " | Train Accuracy: " << A;
        }
        cout << endl;
      }
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

  // The recorded history: one row per epoch under the stochastic solver, one
  // per iteration under LBFGS, which also reports the two quantities it
  // measured to decide it had converged and which of them stopped it.
  // Accuracy is not among them, MATLAB not reporting it either, and computing
  // it would cost a forward pass over the whole set at every iteration.
  octave_idx_type nrec = use_lbfgs
                         ? (octave_idx_type) lbres.history.size ()
                         : (octave_idx_type) Loss.size ();
  RowVector A(nrec), L(nrec), G(nrec), S(nrec);
  for (octave_idx_type i = 0; i < nrec; i++)
  {
    if (use_lbfgs)
    {
      L(i) = lbres.history[i].fval;
      G(i) = lbres.history[i].gradient;
      S(i) = lbres.history[i].step;
    }
    else
    {
      A(i) = Accuracy[i];
      L(i) = Loss[i];
    }
  }

  // Prepare returning arguments
  octave_scalar_map fcnn_model;
  fcnn_model.assign ("LayerWeights", LayerWeights);
  fcnn_model.assign ("Activations", ActiveCode);
  if (use_lbfgs)
  {
    fcnn_model.assign ("Loss", L);
    fcnn_model.assign ("Gradient", G);
    fcnn_model.assign ("Step", S);
    fcnn_model.assign ("Criterion", criterion_message (lbres.crit));
  }
  else
  {
    if (! regression)
    {
      fcnn_model.assign ("Accuracy", A);
    }
    fcnn_model.assign ("Loss", L);
  }
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
%! [Q, W] = fcnntrain (X, Y, 10, "sigmoid", "sigmoid", 1, 0.025, 50, false);
%!error <fcnntrain: X must be a real numeric matrix.> ...
%! fcnntrain (complex (X), Y, 10, "sigmoid", "sigmoid", 1, 0.025, 50, false);
%!error <fcnntrain: X must be a real numeric matrix.> ...
%! fcnntrain ({X}, Y, 10, "sigmoid", "sigmoid", 1, 0.025, 50, false);
%!error <fcnntrain: X cannot be empty.> ...
%! fcnntrain ([], Y, 10, "sigmoid", "sigmoid", 1, 0.025, 50, false);
%!error <fcnntrain: Y must be a real numeric matrix.> ...
%! fcnntrain (X, complex (Y), 10, "sigmoid", "sigmoid", 0.01, 0.025, 50, false);
%!error <fcnntrain: Y must be a real numeric matrix.> ...
%! fcnntrain (X, {Y}, 10, "sigmoid", "sigmoid", 1, 0.025, 50, false);
%!error <fcnntrain: Y cannot be empty.> ...
%! fcnntrain (X, [], 10, "sigmoid", "sigmoid", 1, 0.025, 50, false);
%!error <fcnntrain: X and Y must have the same number of rows.> ...
%! fcnntrain (X, Y([1:50]), 10, "sigmoid", "sigmoid", 1, 0.025, 50, false);
%!error <fcnntrain: labels in Y must be positive integers.> ...
%! fcnntrain (X, Y - 1, 10, "sigmoid", "sigmoid", 1, 0.025, 50, false);
%!error <fcnntrain: 'LayerSizes' must be a row vector of integer values.> ...
%! fcnntrain (X, Y, [10; 5], "sigmoid", "sigmoid", 1, 0.025, 50, false);
%!error <fcnntrain: 'LayerSizes' must be a row vector of integer values.> ...
%! fcnntrain (X, Y, "10", "sigmoid", "sigmoid", 1, 0.025, 50, false);
%!error <fcnntrain: 'LayerSizes' must be a row vector of integer values.> ...
%! fcnntrain (X, Y, {10}, "sigmoid", "sigmoid", 1, 0.025, 50, false);
%!error <fcnntrain: 'LayerSizes' must be a row vector of integer values.> ...
%! fcnntrain (X, Y, complex (10), "sigmoid", "sigmoid", 1, 0.025, 50, false);
%!error <fcnntrain: 'Activations' must be a character vector or a cellstring.> ...
%! fcnntrain (X, Y, 10, [1; 1], "sigmoid", 1, 0.025, 50, false);
%!error <fcnntrain: 'Activations' must be a character vector or a cellstring.> ...
%! fcnntrain (X, Y, 10, {1, 1}, "sigmoid", 1, 0.025, 50, false);
%!error <fcnntrain: 'Activations' must be a character vector or a cellstring.> ...
%! fcnntrain (X, Y, 10, complex ([1, 1]), "sigmoid", 1, 0.025, 50, false);
%!error <fcnntrain: 'Activations' does not match the number of layers.> ...
%! fcnntrain (X, Y, 10, {"sigmoid", "relu"}, "sigmoid", 1, 0.025, 50, false);
%!error <fcnntrain: cannot have a layer of zero size.> ...
%! fcnntrain (X, Y, [10, 0, 5], "sigmoid", "sigmoid", 1, 0.025, 50, false);
%!error <fcnntrain: unsupported 'Activations' function: 'sgmoid'.> ...
%! fcnntrain (X, Y, 10, "sgmoid", "sigmoid", 1, 0.025, 50, false);
%!error <fcnntrain: unsupported 'Activations' function: 'bogus'.> ...
%! fcnntrain (X, Y, 10, {"bogus"}, "sigmoid", 1, 0.025, 50, false);
%!error <fcnntrain: 'OutputLayerActivation' must be a character vector.> ...
%! fcnntrain (X, Y, 10, "sigmoid", 4, 1, 0.025, 50, false);
%!error <fcnntrain: unsupported 'OutputLayerActivation' function: 'softmx'.> ...
%! fcnntrain (X, Y, 10, "sigmoid", "softmx", 1, 0.025, 50, false);
%!error <fcnntrain: 'NumThreads' must be a positive integer scalar value.> ...
%! fcnntrain (X, Y, 10, "sigmoid", "sigmoid", 0, 0.025, 50, false);
%!error <fcnntrain: 'LearningRate' must be a positive scalar value.> ...
%! fcnntrain (X, Y, 10, "sigmoid", "sigmoid", 1, -0.025, 50, false);
%!error <fcnntrain: 'LearningRate' must be a positive scalar value.> ...
%! fcnntrain (X, Y, 10, "sigmoid", "sigmoid", 1, 0, 50, false);
%!error <fcnntrain: 'LearningRate' must be a positive scalar value.> ...
%! fcnntrain (X, Y, 10, "sigmoid", "sigmoid", 1, [0.025, 0.001], 50, false);
%!error <fcnntrain: 'LearningRate' must be a positive scalar value.> ...
%! fcnntrain (X, Y, 10, "sigmoid", "sigmoid", 1, {0.025}, 50, false);
%!error <fcnntrain: 'Epochs' must be a positive scalar value.> ...
%! fcnntrain (X, Y, 10, "sigmoid", "sigmoid", 1, 0.025, 0, false);
%!error <fcnntrain: 'Epochs' must be a positive scalar value.> ...
%! fcnntrain (X, Y, 10, "sigmoid", "sigmoid", 1, 0.025, [50, 25], false);
%!error <fcnntrain: 'DisplayInfo' must be a boolean scalar.> ...
%! fcnntrain (X, Y, 10, "sigmoid", "sigmoid", 1, 0.025, 50, 0);
%!error <fcnntrain: 'DisplayInfo' must be a boolean scalar.> ...
%! fcnntrain (X, Y, 10, "sigmoid", "sigmoid", 1, 0.025, 50, 1);
%!error <fcnntrain: 'DisplayInfo' must be a boolean scalar.> ...
%! fcnntrain (X, Y, 10, "sigmoid", "sigmoid", 1, 0.025, 50, [false, false]);
## The reported training history has one entry per epoch and holds the values
## actually recorded, not the zeros a pre-sized vector would leave in front of
## them.
%!test
%! rand ('seed', 42);
%! randn ('seed', 42);
%! Xs = [randn(30,2)*0.4 + 2; randn(30,2)*0.4 - 2];
%! Ys = [ones(30,1); 2*ones(30,1)];
%! M = fcnntrain (Xs, Ys, 8, "relu", "softmax", 1, 0.05, 60, false, 1);
%! assert_equal (numel (M.Loss), 60);
%! assert_equal (numel (M.Accuracy), 60);
%! assert_equal (any (M.Loss != 0), true);
%! assert_equal (M.Loss(end) < M.Loss(1), true);
%! assert_equal (M.Accuracy(end) >= M.Accuracy(1), true);

%!error <fcnntrain: too many input arguments.> ...
%! fcnntrain (X, Y, 10, "sigmoid", "sigmoid", 1, 0.025, 50, false, 0, struct (), 0);
%!error <fcnntrain: 'SolverOptions' must be a scalar struct.> ...
%! fcnntrain (X, Y, 10, "sigmoid", "sigmoid", 1, 0.025, 50, false, 0, 0);
%!error <fcnntrain: 'Solver' must be 'sgd' or 'lbfgs'.> ...
%! fcnntrain (X, Y, 10, "sigmoid", "sigmoid", 1, 0.025, 50, false, 0, ...
%!            struct ("Solver", "bogus"));
%!error <fcnntrain: 'LossFunction' must be a numeric scalar value.> ...
%! fcnntrain (X, Y, 10, "sigmoid", "sigmoid", 1, 0.025, 50, false, 'ce');
%!error <fcnntrain: 'LossFunction' must be a numeric scalar value.> ...
%! fcnntrain (X, Y, 10, "sigmoid", "sigmoid", 1, 0.025, 50, false, [0, 1]);
%!error <fcnntrain: invalid 'LossFunction' code.> ...
%! fcnntrain (X, Y, 10, "sigmoid", "sigmoid", 1, 0.025, 50, false, 3);
%!error <fcnntrain: invalid 'LossFunction' code.> ...
%! fcnntrain (X, Y, 10, "sigmoid", "sigmoid", 1, 0.025, 50, false, -1);

## Loss function 2 is regression: Y holds the response, the output layer is
## sized to its columns, and no Accuracy is reported because there are no
## labels to count.
%!test
%! rand ('seed', 42);
%! randn ('seed', 42);
%! Xr = linspace (-2, 2, 60)';
%! Yr = 3 * Xr - 1;
%! M = fcnntrain (Xr, Yr, [8, 8], "relu", "linear", 1, 0.005, 300, false, 2);
%! assert_equal (fieldnames (M), {'LayerWeights'; 'Activations'; 'Loss'});
%! assert_equal (rows (M.LayerWeights{end}), 1);
%! assert_equal (M.Loss(end) < M.Loss(1), true);

## The recorded loss is the mean squared error of the network it belongs to.
%!test
%! rand ('seed', 42);
%! Xr = linspace (0, 1, 40)';
%! Yr = 5 * Xr + 2;
%! M = fcnntrain (Xr, Yr, 10, "relu", "linear", 1, 0.005, 200, false, 2);
%! Wm = cellfun (@(m) m(:,1:end-1), M.LayerWeights, "UniformOutput", false);
%! Bm = cellfun (@(m) m(:,end), M.LayerWeights, "UniformOutput", false);
%! [~, yFit] = fcnnpredict (Wm, Bm, "relu", "linear", Xr);
%! assert_equal (M.Loss(end), mean ((Yr - yFit) .^ 2), 1e-12);

## A response of several columns gets one output unit per column.
%!test
%! rand ('seed', 42);
%! Xr = linspace (0, 1, 30)';
%! M = fcnntrain (Xr, [Xr, 2 * Xr, 3 * Xr], 6, "relu", "linear", 1, 0.005, ...
%!                50, false, 2);
%! assert_equal (rows (M.LayerWeights{end}), 3);
%! Wm = cellfun (@(m) m(:,1:end-1), M.LayerWeights, "UniformOutput", false);
%! Bm = cellfun (@(m) m(:,end), M.LayerWeights, "UniformOutput", false);
%! [~, yFit] = fcnnpredict (Wm, Bm, "relu", "linear", Xr);
%! assert_equal (columns (yFit), 3);

## Regression takes a response the classification path would refuse.
%!test
%! rand ('seed', 42);
%! Xr = linspace (0, 1, 20)';
%! Yr = linspace (-3.5, 2.25, 20)';
%! M = fcnntrain (Xr, Yr, 6, "relu", "linear", 1, 0.005, 50, false, 2);
%! assert_equal (all (isfinite (M.Loss)), true);

%!error <fcnntrain: Y must be finite.> ...
%! fcnntrain ([1; 2; 3], [1; Inf; 3], 4, "relu", "linear", 1, 0.01, 10, false, 2);
%!error <fcnntrain: Y must be finite.> ...
%! fcnntrain ([1; 2; 3], [1; NaN; 3], 4, "relu", "linear", 1, 0.01, 10, false, 2);

## The full-batch solver drives the loss down and reports what it measured to
## decide it had stopped.  Accuracy is not among them: MATLAB does not report
## it either, and it would cost a forward pass over the whole set per
## iteration.
%!test
%! so = struct ("Solver", "lbfgs");
%! M = fcnntrain (X, Y, 10, "relu", "softmax", 1, 0.005, 100, false, 1, so);
%! assert_equal (fieldnames (M), ...
%!               {'LayerWeights'; 'Activations'; 'Loss'; 'Gradient'; ...
%!                'Step'; 'Criterion'});
%! assert_equal (M.Loss(end) < M.Loss(1), true);
%! assert_equal (numel (M.Gradient), numel (M.Loss));
%! assert_equal (numel (M.Step), numel (M.Loss));

## It reaches a lower training loss than the epoch loop does, in fewer passes
## over the data, which is the whole reason for offering it.
%!test
%! rand ("state", 3); randn ("state", 3);
%! Ms = fcnntrain (X, Y, 10, "relu", "softmax", 1, 0.005, 200, false, 1);
%! rand ("state", 3); randn ("state", 3);
%! so = struct ("Solver", "lbfgs");
%! Ml = fcnntrain (X, Y, 10, "relu", "softmax", 1, 0.005, 200, false, 1, so);
%! assert_equal (Ml.Loss(end) < Ms.Loss(end), true);
%! assert_equal (numel (Ml.Loss) < numel (Ms.Loss), true);

## An explicit 'sgd' is the epoch loop, unchanged.
%!test
%! rand ("state", 5); randn ("state", 5);
%! Ma = fcnntrain (X, Y, 10, "relu", "softmax", 1, 0.005, 30, false, 1);
%! rand ("state", 5); randn ("state", 5);
%! so = struct ("Solver", "sgd");
%! Mb = fcnntrain (X, Y, 10, "relu", "softmax", 1, 0.005, 30, false, 1, so);
%! assert_equal (Mb.Loss, Ma.Loss);
%! assert_equal (Mb.LayerWeights, Ma.LayerWeights);

## The tolerances reach the solver: from the same starting weights, a loose
## gradient tolerance stops sooner than a tight one.
%!test
%! rand ("state", 9); randn ("state", 9);
%! so = struct ("Solver", "lbfgs", "GradientTolerance", 1e3);
%! Ma = fcnntrain (X, Y, 10, "relu", "softmax", 1, 0.005, 100, false, 1, so);
%! rand ("state", 9); randn ("state", 9);
%! so = struct ("Solver", "lbfgs", "GradientTolerance", 1e-8);
%! Mb = fcnntrain (X, Y, 10, "relu", "softmax", 1, 0.005, 100, false, 1, so);
%! assert_equal (Ma.Criterion, "Relative gradient tolerance reached.");
%! assert_equal (numel (Ma.Loss) < numel (Mb.Loss), true);
*/

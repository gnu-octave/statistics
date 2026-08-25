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
#include <algorithm>
#include <string>
#include <octave/oct-rand.h>
#if defined (_OPENMP)
  #include <omp.h>
  #define MY_OMP_SET_THREADS (omp_set_num_threads (this->n_threads))
#else
  #define MY_OMP_SET_THREADS
#endif

using namespace std;

// Helper functions

// Return random number between -1 and 1, drawn from Octave's generator rather
// than the C library's.  The two are unconnected, so a network seeded from
// rand () ignored rand ('seed', s) and rng () entirely and every fit on the
// same data returned the same model.
double get_random ()
{
  return octave::rand::scalar () * 2 - 1;
}

// Select the uniform distribution for as long as the object is in scope, and
// put back whatever the caller had selected.  The distribution is global state
// owned by the caller, and an error raised while a network is being built must
// not leave it changed.
class uniform_scope
{
public:
  uniform_scope () : saved (octave::rand::distribution ())
  {
    octave::rand::uniform_distribution ();
  }
  ~uniform_scope ()
  {
    octave::rand::distribution (this->saved);
  }

private:
  std::string saved;
};

// Compute accuracy of predicted samples during training
double accuracy (vector<int> predictions, vector<int> labels)
{
  double correct = 0.0;
  for (int i = 0; i < predictions.size (); i++)
  {
    if (predictions[i] == labels[i] - 1)
    {
      correct += 1.0;
    }
  }
  double accuracy = correct / (double) predictions.size ();
  return accuracy;
}

// Class definitions

// Negative-side parameters of the leaky family.  MathWorks uses the same two
// values, leakyReluLayer defaulting to a scale of 0.01 and eluLayer to an
// alpha of 1, the latter being what makes an ELU saturate at -alpha.  Neither
// is a fitted quantity and neither is reachable from the m-code.
static const double LRELU_ALPHA = 0.01;
static const double ELU_ALPHA = 1.0;

// Half-width of the uniform range a layer's weights are drawn from.  A
// constant range ignores how many inputs each neuron sums, so the variance of
// a pre-activation grows with the width of the layer feeding it: deep networks
// then saturate a sigmoid or drive every ReLU unit negative, where it stays
// with a gradient of exactly zero.  Scaling by fan-in holds that variance
// steady from layer to layer.  A rectifier passes only half its input, so it
// needs the wider He range; the symmetric activations take Glorot, which
// accounts for the backward pass through the layer as well.
double init_scale (int activation, int fan_in, int fan_out)
{
  bool rectifier = (activation == 2 || activation == 5 || activation == 6
                    || activation == 7);
  if (rectifier)
  {
    return sqrt (6.0 / fan_in);
  }
  return sqrt (6.0 / (fan_in + fan_out));
}

class Neuron
{
public:
  // constructors
  Neuron (int input_size, double scale);
  Neuron (int input_size);
  // destructor
  ~Neuron ();

  // methods
  double forward (vector<double> inputs);
  void backward (vector<double> last_input, double grad);
  void descend (double learning_rate);
  vector<double> get_neuron ();
  void set_neuron (vector<double> Wb_vector);
  void zero_gradient ();

  // Flat view, laid out as [weights, bias] to match get_neuron.  It exists
  // so that a whole network can be handed to an optimizer as one vector and
  // set back from one, which the per-neuron accessors above cannot do
  // without copying the model on every trial step of a line search.
  int nparams () const;
  void pack (double *p) const;
  void unpack (const double *p);
  void pack_grad (double *p) const;

  // data
  vector<double> weights;
  vector<double> wgrad;
  double bias;
  double bgrad;
  // The gradient this neuron received from the sample being processed, as
  // distinct from bgrad, which accumulates it.  Backpropagation into the
  // layer below needs the one and full-batch training needs the other, and
  // reading bgrad for both is only correct while the gradient is cleared
  // between samples.
  double delta;
};

Neuron::Neuron (int input_size, double scale)
{
  this->weights = vector<double> (input_size);
  this->wgrad = vector<double> (input_size, 0.0);
  this->bias = 0.01 * get_random ();
  this->bgrad = 0.0;
  this->delta = 0.0;
  for (int i = 0; i < input_size; i++)
  {
    this->weights[i] = scale * get_random ();
  }
}

// Parameters left at zero, for a model that is about to have them loaded.
// Drawing them would take numbers from Octave's generator and move the
// caller's random stream on every call to fcnnpredict.
Neuron::Neuron (int input_size)
{
  this->weights = vector<double> (input_size);
  this->wgrad = vector<double> (input_size, 0.0);
  this->bias = 0.0;
  this->bgrad = 0.0;
  this->delta = 0.0;
}
Neuron::~Neuron () {}

double Neuron::forward (vector<double> inputs)
{
  double total = this->bias;
  for (int i = 0; i < inputs.size (); i++)
  {
    total += inputs[i] * this->weights[i];
  }
  return total;
}

void Neuron::backward (vector<double> last_input, double grad)
{
  this->delta = grad;
  this->bgrad += grad;
  for (int i = 0; i < this->wgrad.size (); i++)
  {
    this->wgrad.at (i) = this->wgrad.at (i) + grad * last_input.at (i);
  }
}

void Neuron::descend (double learning_rate)
{
  this->bias -= this->bgrad * learning_rate;
  for (int i = 0; i < this->weights.size (); i++)
  {
    this->weights.at (i) -= this->wgrad.at (i) * learning_rate;
  }
}

vector<double> Neuron::get_neuron ()
{
  vector<double> Wb_vector = this->weights;
  Wb_vector.push_back (this->bias);
  return Wb_vector;
}

void Neuron::set_neuron (vector<double> Wb_vector)
{
  int w_len = Wb_vector.size () - 1;
  for (int i = 0; i < w_len; i++)
  {
    this->weights.at (i) = Wb_vector[i];
  }
  this->bias = Wb_vector[w_len];
}

void Neuron::zero_gradient ()
{
  this->wgrad = vector<double> (this->weights.size ());
  this->bgrad = 0.0;
}

int Neuron::nparams () const
{
  return this->weights.size () + 1;
}

void Neuron::pack (double *p) const
{
  int n = this->weights.size ();
  for (int i = 0; i < n; i++)
  {
    p[i] = this->weights[i];
  }
  p[n] = this->bias;
}

void Neuron::unpack (const double *p)
{
  int n = this->weights.size ();
  for (int i = 0; i < n; i++)
  {
    this->weights[i] = p[i];
  }
  this->bias = p[n];
}

void Neuron::pack_grad (double *p) const
{
  int n = this->wgrad.size ();
  for (int i = 0; i < n; i++)
  {
    p[i] = this->wgrad[i];
  }
  p[n] = this->bgrad;
}

class DenseLayer
{
public:
  // constructors
  DenseLayer (int input_size, int output_size, int activation);
  DenseLayer (int input_size, int output_size);
  // destructor
  ~DenseLayer ();

  // methods
  vector<double> forward (vector<double> inputs);
  void backward (vector<double> grad);
  void descend (double learning_rate);
  vector<vector<double>> get_layer ();
  void set_layer (vector<vector<double>> Wb_matrix);
  void zero_gradient ();

  // Flat view over every neuron of the layer, in neuron order.
  int nparams () const;
  void pack (double *p) const;
  void unpack (const double *p);
  void pack_grad (double *p) const;

  // data
  vector<Neuron> neurons;
  vector<double> last_input;
};

DenseLayer::DenseLayer (int input_size, int output_size, int activation)
{
  // initialize neurons on a range set by the fan-in and the activation
  double scale = init_scale (activation, input_size, output_size);
  this->neurons = vector<Neuron> ();
  for (int i = 0; i < output_size; i++)
  {
    Neuron to_add = Neuron (input_size, scale);
    this->neurons.push_back (to_add);
  }
}

// Layer of zeroed neurons, for a model whose weights are about to be loaded.
DenseLayer::DenseLayer (int input_size, int output_size)
{
  this->neurons = vector<Neuron> ();
  for (int i = 0; i < output_size; i++)
  {
    this->neurons.push_back (Neuron (input_size));
  }
}
DenseLayer::~DenseLayer () {}

vector<double> DenseLayer::forward (vector<double> inputs)
{
  this->last_input = inputs;
  vector<double> outputs = vector<double> (this->neurons.size());
  for (int i = 0; i < this->neurons.size (); i++)
  {
    outputs[i] = this->neurons[i].forward (inputs);
  }
  return outputs;
}

void DenseLayer::backward (vector<double> grad)
{
  for (int i = 0; i < this->neurons.size (); i++)
  {
    this->neurons[i].backward (last_input, grad[i]);
  }
}

void DenseLayer::descend (double learning_rate)
{
  for (int i = 0; i < this->neurons.size (); i++)
  {
    this->neurons[i].descend (learning_rate);
  }
}

vector<vector<double>> DenseLayer::get_layer ()
{
  vector<vector<double>> Wb_matrix;
  for (int i = 0; i < this->neurons.size (); i++)
  {
    vector<double> WB_vector = this->neurons[i].get_neuron ();
    Wb_matrix.push_back (WB_vector);
  }
  return Wb_matrix;
}

void DenseLayer::set_layer (vector<vector<double>> Wb_matrix)
{
  for (int i = 0; i < Wb_matrix.size (); i++)
  {
    this->neurons[i].set_neuron (Wb_matrix[i]);
  }
}

void DenseLayer::zero_gradient ()
{
  for (int i = 0; i < this->neurons.size (); i++)
  {
    this->neurons[i].zero_gradient ();
  }
}

int DenseLayer::nparams () const
{
  int total = 0;
  for (size_t i = 0; i < this->neurons.size (); i++)
  {
    total += this->neurons[i].nparams ();
  }
  return total;
}

void DenseLayer::pack (double *p) const
{
  for (size_t i = 0; i < this->neurons.size (); i++)
  {
    this->neurons[i].pack (p);
    p += this->neurons[i].nparams ();
  }
}

void DenseLayer::unpack (const double *p)
{
  for (size_t i = 0; i < this->neurons.size (); i++)
  {
    this->neurons[i].unpack (p);
    p += this->neurons[i].nparams ();
  }
}

void DenseLayer::pack_grad (double *p) const
{
  for (size_t i = 0; i < this->neurons.size (); i++)
  {
    this->neurons[i].pack_grad (p);
    p += this->neurons[i].nparams ();
  }
}

// The activation functions a layer can carry, and the names the m-code spells
// them with.  The mapping lives here because this is the only place the codes
// mean anything: nothing outside the engine needs to know the numbers.
enum activation_kind
{
  ACT_LINEAR  = 0,
  ACT_SIGMOID = 1,
  ACT_RELU    = 2,
  ACT_TANH    = 3,
  ACT_SOFTMAX = 4,
  ACT_LRELU   = 5,
  ACT_ELU     = 6,
  ACT_GELU    = 7
};

struct activation_entry
{
  const char *name;
  int code;
};

// 'none' is MATLAB's name for the identity, which fitcnet has always
// documented as available; 'linear' is this package's older spelling of the
// same map.  'prelu' and 'lrelu' are one function here, the negative slope
// being a constant rather than a fitted parameter.
static const activation_entry ACTIVATION_TABLE[] =
{
  { "linear",  ACT_LINEAR  },
  { "none",    ACT_LINEAR  },
  { "sigmoid", ACT_SIGMOID },
  { "relu",    ACT_RELU    },
  { "tanh",    ACT_TANH    },
  { "softmax", ACT_SOFTMAX },
  { "lrelu",   ACT_LRELU   },
  { "prelu",   ACT_LRELU   },
  { "elu",     ACT_ELU     },
  { "gelu",    ACT_GELU    }
};

// The code for NAME, matched without regard to case, or -1 if there is none.
static int
activation_code (const string& name)
{
  string key = name;
  for (size_t i = 0; i < key.size (); i++)
  {
    key[i] = tolower (key[i]);
  }
  const int n = sizeof (ACTIVATION_TABLE) / sizeof (ACTIVATION_TABLE[0]);
  for (int i = 0; i < n; i++)
  {
    if (key == ACTIVATION_TABLE[i].name)
    {
      return ACTIVATION_TABLE[i].code;
    }
  }
  return -1;
}

// One code per layer, the output layer last.  ACTS names the hidden layers,
// either as one character vector applying to all NLAYERS of them or as a
// cellstring naming them one by one; OUT_ACT names the output layer.  CALLER
// prefixes every message, so each front end reports under its own name.
static RowVector
activation_codes (const octave_value& acts, const octave_value& out_act,
                  int nlayers, const string& caller)
{
  RowVector codes (nlayers + 1);
  if (acts.is_string () && acts.rows () == 1)
  {
    int code = activation_code (acts.string_value ());
    if (code < 0)
    {
      error ("%s: unsupported 'Activations' function: '%s'.",
             caller.c_str (), acts.string_value ().c_str ());
    }
    for (int i = 0; i < nlayers; i++)
    {
      codes(i) = code;
    }
  }
  else if (acts.iscellstr ())
  {
    Cell c = acts.cell_value ();
    if (c.numel () != nlayers)
    {
      error ("%s: 'Activations' does not match the number of layers.",
             caller.c_str ());
    }
    for (int i = 0; i < nlayers; i++)
    {
      int code = activation_code (c.elem(i).string_value ());
      if (code < 0)
      {
        error ("%s: unsupported 'Activations' function: '%s'.",
               caller.c_str (), c.elem(i).string_value ().c_str ());
      }
      codes(i) = code;
    }
  }
  else
  {
    error ("%s: 'Activations' must be a character vector or a cellstring.",
           caller.c_str ());
  }

  if (! (out_act.is_string () && out_act.rows () == 1))
  {
    error ("%s: 'OutputLayerActivation' must be a character vector.",
           caller.c_str ());
  }
  int code = activation_code (out_act.string_value ());
  if (code < 0)
  {
    error ("%s: unsupported 'OutputLayerActivation' function: '%s'.",
           caller.c_str (), out_act.string_value ().c_str ());
  }
  codes(nlayers) = code;
  return codes;
}

class ActivationLayer
{
public:
  // constructor
  ActivationLayer (int activation, int n_threads);
  // destructor
  ~ActivationLayer ();

  // methods
  vector<double> forward (vector<double> inputs);
  void backward (vector<double> grad);
  void backward (DenseLayer &prev_layer);

  // data
  vector<double> last_input;
  vector<double> grad;
  vector<double> last_output;
private:
  int activation;
  int n_threads;
};

ActivationLayer::ActivationLayer (int activation, int n_threads)
{
  this->activation = activation;
  this->n_threads = n_threads;
}
ActivationLayer::~ActivationLayer () {}

vector<double> ActivationLayer::forward (vector<double> inputs)
{
  this->last_input = inputs;
  int layer_size = inputs.size ();
  if (layer_size < 1000)
  {
    this->n_threads = 1;
  }
  vector<double> outputs = vector<double> (layer_size);
  if (this->activation == 0) // 'Linear'
  {
    outputs = inputs;
  }
  else if (this->activation == 1) // Sigmoid function
  {
    MY_OMP_SET_THREADS;
    #pragma omp parallel
    {
      #pragma omp parallel for
      for (int i = 0; i < layer_size; i++)
      {
        outputs[i] = 1 / (1 + exp (-inputs[i]));
      }
    }
  }
  else if (this->activation == 2) // Rectified Linear Unit (ReLU)
  {
    MY_OMP_SET_THREADS;
    #pragma omp parallel
    {
      #pragma omp parallel for
      for (int i = 0; i < layer_size; i++)
      {
        outputs[i] = inputs[i] > 0 ? inputs[i] : 0;
      }
    }
  }
  else if (this->activation == 3) // Hyperbolic tangent (tanh)
  {
    MY_OMP_SET_THREADS;
    #pragma omp parallel
    {
      #pragma omp parallel for
      for (int i = 0; i < layer_size; i++)
      {
        double ex = exp (inputs[i]);
        double e_x = exp (-inputs[i]);
        outputs[i] = (ex - e_x) / (ex + e_x);
      }
    }
  }
  else if (this->activation == 4) // Softmax activation
  {
    double total = 0.0;
    double maxel = *max_element (inputs.begin (), inputs.end ());
    for (int i = 0; i < layer_size; i++)
    {
        outputs[i] = exp (inputs[i] - maxel);
        total += outputs[i];
    }
    for (int i = 0; i < layer_size; i++)
    {
        outputs[i] /= total;
    }
  }
  else if (this->activation == 5) // Parametric or Leaky ReLU
  {
    MY_OMP_SET_THREADS;
    #pragma omp parallel
    {
      #pragma omp parallel for
      for (int i = 0; i < layer_size; i++)
      {
        outputs[i] = inputs[i] >= 0 ? inputs[i] : inputs[i] * LRELU_ALPHA;
      }
    }
  }
  else if (this->activation == 6) // Exponential Linear Unit (ELU)
  {
    MY_OMP_SET_THREADS;
    #pragma omp parallel
    {
      #pragma omp parallel for
      for (int i = 0; i < layer_size; i++)
      {
        outputs[i] = inputs[i] >= 0 ? inputs[i] : (exp (inputs[i]) - 1)
                                                  * ELU_ALPHA;
      }
    }
  }
  else if (this->activation == 7) // Gaussian Error Linear Unit (GELU)
  {
    MY_OMP_SET_THREADS;
    #pragma omp parallel
    {
      #pragma omp parallel for
      for (int i = 0; i < layer_size; i++)
      {
        // x * Phi(x), the standard normal CDF taken from erfc so that the
        // tails are accurate; the tanh form is an approximation of this.
        outputs[i] = inputs[i] * 0.5 * erfc (-inputs[i] * M_SQRT1_2);
      }
    }
  }
  this->last_output = outputs;
  return outputs;
}

void ActivationLayer::backward (vector<double> chain_grad)
{
  int layer_size = this->last_input.size ();
  if (layer_size < 1000)
  {
    this->n_threads = 1;
  }
  this->grad = vector<double> (layer_size);
  if (this->activation == 0) // 'Linear'
  {
    this->grad = chain_grad;
  }
  else if (this->activation == 1) // Sigmoid function
  {
    MY_OMP_SET_THREADS;
    #pragma omp parallel
    {
      #pragma omp parallel for
      for (int i = 0; i < layer_size; i++)
      {
        this->grad[i] = this->last_output[i] * (1 - this->last_output[i])
                                             * chain_grad[i];
      }
    }
  }
  else if (this->activation == 2) // Rectified Linear Unit (ReLU)
  {
    MY_OMP_SET_THREADS;
    #pragma omp parallel
    {
      #pragma omp parallel for
      for (int i = 0; i < layer_size; i++)
      {
        this->grad[i] = this->last_input[i] > 0 ? chain_grad[i] : 0;
      }
    }
  }
  else if (this->activation == 3) // Hyperbolic tangent (tanh)
  {
    MY_OMP_SET_THREADS;
    #pragma omp parallel
    {
      #pragma omp parallel for
      for (int i = 0; i < layer_size; i++)
      {
        this->grad[i] = (1 - pow (this->last_output[i], 2)) * chain_grad[i];
      }
    }
  }
  else if (this->activation == 4) // Softmax activation
  {
    // The Jacobian is diag (y) - y * y', so the product with the incoming
    // gradient is y .* (g - y' * g).  Passing g through unchanged is only
    // right when softmax is paired with a cross-entropy loss, which this
    // implementation does not provide.
    double dot = 0.0;
    for (int i = 0; i < layer_size; i++)
    {
      dot += this->last_output[i] * chain_grad[i];
    }
    for (int i = 0; i < layer_size; i++)
    {
      this->grad[i] = this->last_output[i] * (chain_grad[i] - dot);
    }
  }
  else if (this->activation == 5) // Parametric or Leaky ReLU
  {
    MY_OMP_SET_THREADS;
    #pragma omp parallel
    {
      #pragma omp parallel for
      for (int i = 0; i < layer_size; i++)
      {
        this->grad[i] = this->last_input[i] >= 0 ? chain_grad[i] :
                                                   chain_grad[i] * LRELU_ALPHA;
      }
    }
  }
  else if (this->activation == 6) // Exponential Linear Unit (ELU)
  {
    MY_OMP_SET_THREADS;
    #pragma omp parallel
    {
      #pragma omp parallel for
      for (int i = 0; i < layer_size; i++)
      {
        this->grad[i] = this->last_input[i] >= 0 ? chain_grad[i] : chain_grad[i]
                                      * exp (this->last_input[i]) * ELU_ALPHA;
      }
    }
  }
  else if (this->activation == 7) // Gaussian Error Linear Unit (GELU)
  {
    // d/dx [x * Phi(x)] = Phi(x) + x * phi(x), both taken at the input
    static const double inv_sqrt_2pi = 0.3989422804014327;
    MY_OMP_SET_THREADS;
    #pragma omp parallel
    {
      #pragma omp parallel for
      for (int i = 0; i < layer_size; i++)
      {
        double x = this->last_input[i];
        double cdf = 0.5 * erfc (-x * M_SQRT1_2);
        double pdf = inv_sqrt_2pi * exp (-0.5 * x * x);
        this->grad[i] = (cdf + x * pdf) * chain_grad[i];
      }
    }
  }
}

void ActivationLayer::backward (DenseLayer &prev_layer)
{
  // The gradient arriving here is W' * delta, where delta is the gradient at
  // each neuron of the layer that follows.  Every neuron stores the delta of
  // the sample being processed, so it is read from there rather than
  // recovered by dividing the weight gradient by the input, which is
  // undefined wherever an activation output is zero.  Reading the bias
  // gradient instead is only correct while it is cleared between samples,
  // which full-batch training does not do.  The local derivative is
  // then applied by the overload above, so both paths share one definition of
  // every activation.
  int layer_size = this->last_input.size ();
  vector<double> chain_grad = vector<double> (layer_size, 0.0);
  for (int n = 0; n < prev_layer.neurons.size (); n++)
  {
    double delta = prev_layer.neurons[n].delta;
    for (int i = 0; i < layer_size; i++)
    {
      chain_grad[i] += prev_layer.neurons[n].weights[i] * delta;
    }
  }
  this->backward (chain_grad);
}

class MeanSquaredErrorLoss
{
public:
  MeanSquaredErrorLoss ();
  ~MeanSquaredErrorLoss ();

  double forward (vector<double> inputs, vector<double> targets);
  void backward (double grad);

  // data
  vector<double> last_input;
  vector<double> last_target;
  vector<double> grad;
};

MeanSquaredErrorLoss::MeanSquaredErrorLoss () {}
MeanSquaredErrorLoss::~MeanSquaredErrorLoss () {}

double MeanSquaredErrorLoss::forward (vector<double> inputs,
                                      vector<double> targets)
{
  // we only need to calculate the loss for the target class
  this->last_input = inputs;
  this->last_target = targets;

  double total = 0;

  for (int i = 0; i < inputs.size (); i++)
  {
    total += pow (inputs[i] - targets[i], 2);
  }

  double loss = total;
  return loss;
}

void MeanSquaredErrorLoss::backward (double grad)
{
  this->grad = vector<double> (this->last_input.size ());
  for (int i = 0; i < this->last_input.size (); i++)
  {
    // d/dy of sum (y - t)^2 is 2 * (y - t); the bracket matters, since
    // 2 * y - t is a different function wherever the target is not zero.
    this->grad.at(i) = 2 * (this->last_input[i] - this->last_target[i]);
    this->grad.at(i) *= grad;
  }
}

// Cross-entropy loss, for an output layer that reports a probability over the
// classes.  Paired with softmax the composition of the two gradients reduces
// to y - t, which is the pairing this loss exists for.
class CrossEntropyLoss
{
public:
  CrossEntropyLoss ();
  ~CrossEntropyLoss ();

  double forward (vector<double> inputs, vector<double> targets);
  void backward (double grad);

  // data
  vector<double> last_input;
  vector<double> last_target;
  vector<double> grad;
};

CrossEntropyLoss::CrossEntropyLoss () {}
CrossEntropyLoss::~CrossEntropyLoss () {}

// A predicted probability of zero for the true class carries infinite loss, so
// the logarithm and the division below are floored rather than allowed to
// diverge.
static const double CE_FLOOR = 1e-15;

double CrossEntropyLoss::forward (vector<double> inputs, vector<double> targets)
{
  this->last_input = inputs;
  this->last_target = targets;

  double total = 0.0;

  for (int i = 0; i < inputs.size (); i++)
  {
    if (targets[i] != 0.0)
    {
      total -= targets[i] * log (inputs[i] > CE_FLOOR ? inputs[i] : CE_FLOOR);
    }
  }
  return total;
}

void CrossEntropyLoss::backward (double grad)
{
  this->grad = vector<double> (this->last_input.size ());
  for (int i = 0; i < this->last_input.size (); i++)
  {
    double y = this->last_input[i] > CE_FLOOR ? this->last_input[i] : CE_FLOOR;
    this->grad.at (i) = -this->last_target[i] / y;
    this->grad.at (i) *= grad;
  }
}


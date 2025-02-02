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
#include <omp.h>

using namespace std;

// Helper functions

// Return random number between -1 and 1
double get_random ()
{
    // return random number between -1 and 1
    return (rand () / (double) RAND_MAX) * 2 - 1;
}

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

class Neuron
{
public:
  // constructor
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

  // data
  vector<double> weights;
  vector<double> wgrad;
  double bias;
  double bgrad;
};

Neuron::Neuron (int input_size)
{
  this->weights = vector<double> (input_size);
  this->bias = 0.01 * get_random ();
  for (int i = 0; i < input_size; i++)
  {
    this->weights[i] = get_random ();
  }
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
  for (int i; i < w_len; i++)
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

class DenseLayer
{
public:
  // constructor
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

  // data
  vector<Neuron> neurons;
  vector<double> last_input;
};

DenseLayer::DenseLayer (int input_size, int output_size)
{
  // initialize neurons
  this->neurons = vector<Neuron> ();
  for (int i = 0; i < output_size; i++)
  {
    Neuron to_add = Neuron (input_size);
    this->neurons.push_back (to_add);
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
  for (int i; i < Wb_matrix.size (); i++)
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

class ActivationLayer
{
public:
  // constructor
  ActivationLayer (int activation, int n_threads, double alpha);
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
  double alpha;
};

ActivationLayer::ActivationLayer (int activation, int n_threads, double alpha)
{
  this->activation = activation;
  this->n_threads = n_threads;
  this->alpha = alpha;
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
    omp_set_num_threads (this->n_threads);
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
    omp_set_num_threads (this->n_threads);
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
    omp_set_num_threads (this->n_threads);
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
    omp_set_num_threads (this->n_threads);
    #pragma omp parallel
    {
      #pragma omp parallel for
      for (int i = 0; i < layer_size; i++)
      {
        outputs[i] = inputs[i] >= 0 ? inputs[i] : inputs[i] * this->alpha;
      }
    }
  }
  else if (this->activation == 6) // Exponential Linear Unit (ELU)
  {
    omp_set_num_threads (this->n_threads);
    #pragma omp parallel
    {
      #pragma omp parallel for
      for (int i = 0; i < layer_size; i++)
      {
        outputs[i] = inputs[i] >= 0 ? inputs[i] : (exp (inputs[i]) - 1)
                                                  * this->alpha;
      }
    }
  }
  else if (this->activation == 7) // Gaussian Error Linear Unit (GELU)
  {
    omp_set_num_threads (this->n_threads);
    #pragma omp parallel
    {
      #pragma omp parallel for
      for (int i = 0; i < layer_size; i++)
      {
        double x_3 = pow (inputs[i], 3) * 0.044715;
        outputs[i] = 0.5 * inputs[i] * (tanh (sqrt (2 / M_PI)
                                     * (inputs[i] + x_3)));
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
    omp_set_num_threads (this->n_threads);
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
    omp_set_num_threads (this->n_threads);
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
    omp_set_num_threads (this->n_threads);
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
    // WARNING: this code may be incorrect
    this->grad = chain_grad;
  }
  else if (this->activation == 5) // Parametric or Leaky ReLU
  {
    omp_set_num_threads (this->n_threads);
    #pragma omp parallel
    {
      #pragma omp parallel for
      for (int i = 0; i < layer_size; i++)
      {
        this->grad[i] = this->last_input[i] >= 0 ? chain_grad[i] :
                                                   chain_grad[i] * this->alpha;
      }
    }
  }
  else if (this->activation == 6) // Exponential Linear Unit (ELU)
  {
    omp_set_num_threads (this->n_threads);
    #pragma omp parallel
    {
      #pragma omp parallel for
      for (int i = 0; i < layer_size; i++)
      {
        this->grad[i] = this->last_input[i] >= 0 ? chain_grad[i] : chain_grad[i]
                                     * exp (this->last_output[i]) * this->alpha;
      }
    }
  }
  else if (this->activation == 7) // Gaussian Error Linear Unit (GELU)
  {
    // WARNING: this code may be incorrect
    static const double inv_sqrt_2pi = 0.3989422804014327;
    omp_set_num_threads (this->n_threads);
    #pragma omp parallel
    {
      #pragma omp parallel for
      for (int i = 0; i < layer_size; i++)
      {
        double pdf_val = inv_sqrt_2pi * exp (-0.5 * this->last_output[i]
                                                  * this->last_output[i]);
        this->grad[i] = this->last_output[i] * pdf_val * chain_grad[i];
      }
    }
  }
}

void ActivationLayer::backward (DenseLayer &prev_layer)
{
  this->grad = vector<double> (this->last_input.size ());
  if (this->activation == 0) // 'Linear'
  {
    for (int i = 0; i < prev_layer.last_input.size (); i++)
    {
      for (int n = 0; n < prev_layer.neurons.size (); n++)
      {
        double curr_grad = this->last_output[n];
        double chain_grad = prev_layer.neurons[n].weights[i] *
                            prev_layer.neurons[n].wgrad[i] /
                            prev_layer.last_input[i];
        this->grad[i] += curr_grad * chain_grad;
      }
    }
  }
  else if (this->activation == 1) // Sigmoid function
  {
    for (int i = 0; i < prev_layer.last_input.size (); i++)
    {
      for (int n = 0; n < prev_layer.neurons.size (); n++)
      {
        double curr_grad = this->last_output[n] * (1 - this->last_output[n]);
        double chain_grad = prev_layer.neurons[n].weights[i] *
                            prev_layer.neurons[n].wgrad[i] /
                            prev_layer.last_input[i];
        this->grad[i] += curr_grad * chain_grad;
      }
    }
  }
  else if (this->activation == 2) // Rectified Linear Unit (ReLU)
  {
    for (int i = 0; i < prev_layer.last_input.size(); i++)
    {
      for (int n = 0; n < prev_layer.neurons.size(); n++)
      {
        double grad = prev_layer.neurons[n].wgrad[i] / prev_layer.last_input[i];
        if (this->last_input[i] < 0)
        {
          this->grad[i] = 0;
        }
        else
        {
          this->grad[i] += prev_layer.neurons[n].weights[i] * grad;
        }
      }
    }
  }
  else if (this->activation == 3) // Hyperbolic tangent (tanh)
  {
    for (int i = 0; i < prev_layer.last_input.size (); i++)
    {
      for (int n = 0; n < prev_layer.neurons.size (); n++)
      {
        double curr_grad = this->last_output[n] *
                           (1 - pow (this->last_output[i], 2));
        double chain_grad = prev_layer.neurons[n].weights[i] *
                            prev_layer.neurons[n].wgrad[i] /
                            prev_layer.last_input[i];
        this->grad[i] += curr_grad * chain_grad;
      }
    }
  }
  else if (this->activation == 4) // Softmax activation
  {
    // WARNING: this code may be incorrect
    for (int i = 0; i < prev_layer.last_input.size (); i++)
    {
      for (int n = 0; n < prev_layer.neurons.size (); n++)
      {
        double curr_grad = this->last_output[n];
        double chain_grad = prev_layer.neurons[n].weights[i] *
                            prev_layer.neurons[n].wgrad[i] /
                            prev_layer.last_input[i];
        this->grad[i] += curr_grad * chain_grad;
      }
    }
  }
  else if (this->activation == 5) // Parametric or Leaky ReLU
  {
    for (int i = 0; i < prev_layer.last_input.size(); i++)
    {
      for (int n = 0; n < prev_layer.neurons.size(); n++)
      {
        double grad = prev_layer.neurons[n].wgrad[i] / prev_layer.last_input[i];
        this->grad[i] += prev_layer.neurons[n].weights[i] * grad;
      }
    }
    for (int i = 0; i < this->last_input.size(); i++)
    {
      this->grad[i] = this->last_input[i] >= 0 ? this->grad[i] :
                                                 this->grad[i] * this->alpha;
    }
  }
  else if (this->activation == 6) // Exponential Linear Unit (ELU)
  {
    for (int i = 0; i < prev_layer.last_input.size(); i++)
    {
      for (int n = 0; n < prev_layer.neurons.size(); n++)
      {
        double grad = prev_layer.neurons[n].wgrad[i] / prev_layer.last_input[i];
        this->grad[i] += prev_layer.neurons[n].weights[i] * grad;
      }
    }
    for (int i = 0; i < this->last_input.size(); i++)
    {
      this->grad[i] = this->last_input[i] >= 0 ? this->grad[i] : this->grad[i]
                                   * exp (this->last_output[i]) * this->alpha;
    }
  }
  else if (this->activation == 7) // Gaussian Error Linear Unit (GELU)
  {
    // WARNING: this code may be incorrect
    static const double inv_sqrt_2pi = 0.3989422804014327;
    for (int i = 0; i < prev_layer.last_input.size (); i++)
    {
      for (int n = 0; n < prev_layer.neurons.size (); n++)
      {
        double pdf_val = inv_sqrt_2pi * exp (-0.5 * this->last_output[i]
                                                  * this->last_output[i]);
        double curr_grad = this->last_output[n] * pdf_val * this->last_output[i];
        double chain_grad = prev_layer.neurons[n].weights[i] *
                            prev_layer.neurons[n].wgrad[i] /
                            prev_layer.last_input[i];
        this->grad[i] += curr_grad * chain_grad;
      }
    }
  }
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
    this->grad.at(i) = 2 * this->last_input[i] - this->last_target[i];
    this->grad.at(i) *= grad;
  }
}

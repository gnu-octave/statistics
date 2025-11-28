/*
Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
Copyright (C) 2025 Avanish Salunke <avanishsalunke16@gmail.com>

Based on the Octave LIBSVM wrapper adapted by Alan Meeson (2014) based on an
earlier version of the LIBSVM (3.18) library for MATLAB. Current implementation
is based on LIBSVM 3.36 (2025) by Chih-Chung Chang and Chih-Jen Lin.

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

#include <octave/oct.h>
#include <octave/dMatrix.h>
#include <octave/dColVector.h>
#include <octave/dRowVector.h>
#include <octave/ov.h>
#include <octave/ov-base.h>
#include <octave/ov-struct.h>
#include <string.h>
#include "svm.h"

#include "svm_model_octave.h"

#define CMD_LEN 2048

int print_null(const char *s,...) {return 0;}
int (*info)(const char *fmt,...) = &printf;

void read_sparse_instance(const SparseMatrix &args, int index, struct svm_node *x)
{
	int i, j, low, high;
	octave_idx_type *ir, *jc;
	double *samples;

	ir = (octave_idx_type*)args.ridx();
	jc = (octave_idx_type*)args.cidx();
	samples = (double*)args.data();

	// each column is one instance
	j = 0;
	low = (int)jc[index], high = (int)jc[index+1];
	for(i=low;i<high;i++)
	{
		x[j].index = (int)ir[i] + 1;
		x[j].value = samples[i];
		j++;
	}
	x[j].index = -1;
}

static void fake_answer(int nlhs, octave_value_list &plhs)
{
	int i;
	for(i=0;i<nlhs;i++) plhs(i) = Matrix(0,0);
}

void predict(int nlhs, octave_value_list &plhs, const octave_value_list &args,
             struct svm_model *model, const int predict_probability)
{
	int label_vector_row_num, label_vector_col_num;
	int feature_number, testing_instance_number;
	int instance_index;
	double *ptr_instance, *ptr_label, *ptr_predict_label;
	double *ptr_prob_estimates, *ptr_dec_values, *ptr;
	struct svm_node *x;
	SparseMatrix pplhs(0,0); // transposed instance sparse matrix
	octave_value_list tplhs(3); // temporary storage for plhs[]

	int correct = 0;
	int total = 0;
	double error = 0;
	double sump = 0, sumt = 0, sumpp = 0, sumtt = 0, sumpt = 0;

	int svm_type = svm_get_svm_type(model);
	int nr_class = svm_get_nr_class(model);
	double *prob_estimates = NULL;

	// args[1] = testing instance matrix
	feature_number = (int)args(1).columns();
	testing_instance_number = (int)args(1).rows();
	label_vector_row_num = (int)args(0).rows();
	label_vector_col_num = (int)args(0).columns();

	if(label_vector_row_num != testing_instance_number)
	{
		printf("svmpredict: length of label vector does not match # of instances.\n");
		fake_answer(nlhs, plhs);
		return;
	}
	if(label_vector_col_num != 1)
	{
		printf("svmpredict: label (1st argument) should be a vector (# of column is 1).\n");
		fake_answer(nlhs, plhs);
		return;
	}

  ColumnVector label_vec = args(0).matrix_value();
	ptr_label    = (double*)label_vec.data();

	// transpose instance matrix
	Matrix t_data(0,0);
	if(args(1).issparse())
	{
		if(model->param.kernel_type == PRECOMPUTED)
		{
			// precomputed kernel requires dense matrix, so we make one
			t_data = args(1).matrix_value();
		}
		else
		{
			//If it's a sparse matrix with a non PRECOMPUTED kernel, transpose it
			pplhs = args(1).sparse_matrix_value().transpose();
		}
  }
  else
  {
		t_data = args(1).matrix_value();
	}
	ptr_instance = (double*)t_data.data();
	if(predict_probability)
	{
		if(svm_type==NU_SVR || svm_type==EPSILON_SVR)
    {
			info("Prob. model for test data: target value = predicted value + z,\nz: Laplace distribution e^(-|z|/sigma)/(2sigma),sigma=%g\n", svm_get_svr_probability(model));
    }
		else
    {
			prob_estimates = (double*) malloc(nr_class*sizeof(double));
    }
	}

	ColumnVector cv_predictions(testing_instance_number);
	tplhs(0) = cv_predictions;
	if(predict_probability)
	{
		// prob estimates are in plhs[2]
		if(svm_type == C_SVC || svm_type == NU_SVC || svm_type == ONE_CLASS)
    {
			Matrix m_pe(testing_instance_number, nr_class);
			tplhs(2) = m_pe;
		}
    else
    {
			Matrix m_pe(0,0);
			tplhs(2) = m_pe;
		}
	}
	else
	{
		// decision values are in plhs[2]
		if(svm_type == ONE_CLASS ||
		   svm_type == EPSILON_SVR ||
		   svm_type == NU_SVR ||
		   nr_class == 1)
    {
			 Matrix m_pe(testing_instance_number, 1);
			 tplhs(2) = m_pe;
		}
    else
    {
			Matrix m_pe(testing_instance_number, nr_class*(nr_class-1)/2);
			tplhs(2) = m_pe;
		}
	}

	ptr_predict_label = (double*)tplhs(0).column_vector_value().data();
	ptr_prob_estimates = (double*)tplhs(2).matrix_value().data();
	ptr_dec_values = (double*)tplhs(2).matrix_value().data();
	x = (struct svm_node*)malloc((feature_number+1)*sizeof(struct svm_node));
	for(instance_index=0;instance_index<testing_instance_number;instance_index++)
	{
		int i;
		double target_label, predict_label;

		target_label = ptr_label[instance_index];

		if(args(1).issparse() && model->param.kernel_type != PRECOMPUTED)
    {
			read_sparse_instance(pplhs, instance_index, x);
    }
		else
		{
			for(i = 0; i < feature_number; i++)
			{
				x[i].index = i+1;
				x[i].value = ptr_instance[testing_instance_number*i+instance_index];
			}
			x[feature_number].index = -1;
		}

		if(predict_probability)
		{
			if(svm_type == C_SVC || svm_type == NU_SVC || svm_type == ONE_CLASS)
			{
				predict_label = svm_predict_probability(model, x, prob_estimates);
				ptr_predict_label[instance_index] = predict_label;
				for(i = 0; i < nr_class; i++)
        {
					ptr_prob_estimates[instance_index + i * testing_instance_number] = prob_estimates[i];
        }
			}
      else
      {
				predict_label = svm_predict(model,x);
				ptr_predict_label[instance_index] = predict_label;
			}
		}
		else
		{
			if(svm_type == ONE_CLASS ||
			   svm_type == EPSILON_SVR ||
			   svm_type == NU_SVR)
			{
				double res;
				predict_label = svm_predict_values(model, x, &res);
				ptr_dec_values[instance_index] = res;
			}
			else
			{
				double *dec_values = (double *) malloc(sizeof(double) * nr_class*(nr_class-1)/2);
				predict_label = svm_predict_values(model, x, dec_values);
				if(nr_class == 1)
        {
					ptr_dec_values[instance_index] = 1;
        }
				else
        {
					for(i = 0; i < (nr_class * (nr_class - 1)) / 2; i++)
          {
						ptr_dec_values[instance_index + i * testing_instance_number] = dec_values[i];
          }
        }
				free(dec_values);
			}
			ptr_predict_label[instance_index] = predict_label;
		}

		if(predict_label == target_label)
    {
			++correct;
    }
		error += (predict_label-target_label)*(predict_label-target_label);
		sump += predict_label;
		sumt += target_label;
		sumpp += predict_label*predict_label;
		sumtt += target_label*target_label;
		sumpt += predict_label*target_label;
		++total;
	}
	if(svm_type==NU_SVR || svm_type==EPSILON_SVR)
	{
		info("Mean squared error = %g (regression)\n",error/total);
		info("Squared correlation coefficient = %g (regression)\n",
			  ((total*sumpt-sump*sumt)*(total*sumpt-sump*sumt))/
			  ((total*sumpp-sump*sump)*(total*sumtt-sumt*sumt)));
	}
	else
  {
		info("Accuracy = %g%% (%d/%d) (classification)\n",
			  (double)correct/total*100,correct,total);
  }
	// return accuracy, mean squared error, squared correlation coefficient
	ColumnVector cv_acc(3);
	ptr = (double*)cv_acc.data();
	ptr[0] = (double)correct/total*100;
	ptr[1] = error/total;
	ptr[2] = ((total*sumpt-sump*sumt)*(total*sumpt-sump*sumt))/
				((total*sumpp-sump*sump)*(total*sumtt-sumt*sumt));
	tplhs(1) = cv_acc;
	free(x);
	if(prob_estimates != NULL)
  {
		free(prob_estimates);
  }
	switch(nlhs)
	{
		case 3:
			plhs(2) = tplhs(2);
			plhs(1) = tplhs(1);
		case 1:
		case 0:
			plhs(0) = tplhs(0);
	}
}


DEFUN_DLD (svmpredict, args, nargout,
           "-*- texinfo -*- \n\n\
 @deftypefn  {statistics} {@var{predicted_label} =} svmpredict (@var{labels}, @var{data}, @var{model})\n\
 @deftypefnx {statistics} {@var{predicted_label} =} svmpredict (@var{labels}, @var{data}, @var{model}, ""libsvm_options"")\n\
 @deftypefnx {statistics} {[@var{predicted_label}, @var{accuracy}, @var{decision_values}] =} svmpredict (@var{labels}, @var{data}, @var{model}, ""libsvm_options"")\n\
 @deftypefnx {statistics} {[@var{predicted_label}, @var{accuracy}, @var{prob_estimates}] =} svmpredict (@var{labels}, @var{data}, @var{model}, ""libsvm_options"")\n\
\n\
\n\
This function predicts new labels from a testing instance matrix based on an \
SVM @var{model} created with @code{svmtrain}. \
\n\
\n\
@itemize \n\
@item @var{labels} : An m by 1 vector of prediction labels. If labels \
of test data are unknown, simply use any random values. (type must be double) \
\n\
\n\
@item @var{data} : An m by n matrix of m testing instances with n features. \
It can be dense or sparse. (type must be double) \
\n\
\n\
@item @var{model} : The output of @code{svmtrain} function. \
\n\
\n\
@item @code{libsvm_options} : A string of testing options in the same format \
as that of LIBSVM. \
\n\
\n\
@end itemize \
\n\
\n\
@code{libsvm_options} :\n\
\n\
@itemize \n\
@item @code{-b} : probability_estimates; whether to predict probability \
estimates.\n\
\n\
@end itemize \
@multitable @columnfractions 0.1 0.1 0.8 \n\
@item @tab 0 @tab return decision values. (default) \n\
\n\
@item @tab 1 @tab return probability estimates. \
\n\
@end multitable \
\n\
\n\
@itemize \n\
@item @code{-q} : quiet mode. (no outputs) \
\n\
@end itemize \
\n\
\n\
The @code{svmpredict} function has three outputs.  The first one, \
@var{predicted_label}, is a vector of predicted labels.  The second output, \
@var{accuracy}, is a vector including accuracy (for classification), mean \
squared error, and squared correlation coefficient (for regression).  The \
third is a matrix containing decision values or probability estimates \
(if @code{-b 1}' is specified).  If @math{k} is the number of classes in \
training data, for decision values, each row includes results of predicting \
@math{k(k-1)/2} binary-class SVMs.  For classification, @math{k = 1} is a \
special case.  Decision value +1 is returned for each testing instance, \
instead of an empty vector.  For probabilities, each row contains @math{k} \
values indicating the probability that the testing instance is in each class.  \
Note that the order of classes here is the same as @code{Label} field in the \
@var{model} structure. \
\n\
\n\
\\\n\\\
@emph{Note on LIBSVM 3.36 Update}: This implementation is based on LIBSVM 3.36 (2025) and now supports probability estimates for One-Class SVM (@code{-s 2}) when combined with the probability flag (@code{-b 1}). For One-Class SVM, the @var{prob_estimates} output is a single column vector containing the probability of the instance being an inlier. \
\\\n\\\
@end deftypefn")
{
{
	int nlhs = nargout;
	int nrhs = args.length();
	octave_value_list plhs(nlhs);
	int prob_estimate_flag = 0;
	struct svm_model *model;
	info = &print_null;

	if(nlhs == 2 || nlhs > 3)
	{
		error ("svmpredict: wrong number of output arguments.");
	}
  if(nrhs > 4 || nrhs < 3)

	{
		error ("svmpredict: wrong number of input arguments.");
	}

	if(!args(0).is_double_type() || !args(1).is_double_type())
  {
		error ("svmpredict: label vector and instance matrix must be double.");
	}

	if(args(2).isstruct())
	{
		const char *error_msg;

		// parse options
		if(nrhs==4)
		{
			int i, argc = 1;
			char cmd[CMD_LEN], *argv[CMD_LEN/2];

			// put options in argv[]
			strncpy(cmd, args(3).string_value().c_str(), CMD_LEN);
			if((argv[argc] = strtok(cmd, " ")) != NULL)
      {
				while((argv[++argc] = strtok(NULL, " ")) != NULL);
      }
			for(i=1;i<argc;i++)
			{
				if(argv[i][0] != '-') break;
				if((++i>=argc) && argv[i-1][1] != 'q')
				{
					fake_answer(nlhs, plhs);
					return plhs;
				}
				switch(argv[i-1][1])
				{
					case 'b':
						prob_estimate_flag = atoi(argv[i]);
						break;
					case 'q':
						i--;
						info = &print_null;
						break;
					default:
						printf("svmpredict: unknown option: -%c\n", argv[i-1][1]);
						fake_answer(nlhs, plhs);
						return plhs;
				}
			}
		}
		octave_scalar_map osm_model = args(2).scalar_map_value();
		model = octave_matrix_to_model(osm_model, &error_msg);
		if (model == NULL)
		{
			printf("svmpredict: can't read model: %s\n", error_msg);
			fake_answer(nlhs, plhs);
			return plhs;
		}

		if(prob_estimate_flag)
		{
			// Check if the SVM type supports probability, new support for ONE_CLASS
			if (model->param.svm_type != C_SVC && model->param.svm_type != NU_SVC && model->param.svm_type != ONE_CLASS)
			{
				svm_free_and_destroy_model(&model);
				error ("svmpredict: probability estimates are not supported for this SVM type (only C-SVC, NU-SVC, and ONE-CLASS).\n");
			}

			// Check if the model itself was trained with probability info (-b 1)
			if(svm_check_probability_model(model)==0)
			{
				svm_free_and_destroy_model(&model);
				error ("svmpredict: model does not support probability estimates. Train with '-b 1'.\n");
			}
		}
		else
		{
			if(svm_check_probability_model(model)!=0)
				info("Model supports probability estimates, but disabled in prediction.\n");
		}

		predict(nlhs, plhs, args, model, prob_estimate_flag);
		// destroy model
		svm_free_and_destroy_model(&model);
	}
	else
	{
		error ("svmpredict: model should be a struct array.");
	}
	return plhs;
}

/*
%!test
%! # Test 1: Standard C-SVC Prediction (Original Regression Test)
%! [L, D] = libsvmread (file_in_loadpath ("heart_scale.dat"));
%! model = svmtrain (L, D, '-c 1 -g 0.07');
%! [predict_label, accuracy, dec_values] = svmpredict (L, D, model);
%! assert (size (predict_label), size (dec_values));
%! assert (accuracy, [86.666, 0.533, 0.533]', [1e-3, 1e-3, 1e-3]');
%! assert (dec_values(1), 1.225836001973273, 1e-14);
%! assert (dec_values(2), -0.3212992933043805, 1e-14);
%! assert (predict_label(1), 1);
%!
%!test
%! # Test 2: One-Class Probability (NEW LIBSVM 3.36 FEATURE)
%! [L, D] = libsvmread (file_in_loadpath ("heart_scale.dat"));
%! # Train One-Class (-s 2) with Probability (-b 1)
%! model_oc = svmtrain (L, D, '-s 2 -n 0.1 -g 0.07 -b 1');
%! # FIX: Changed // to # below to fix syntax error
%! assert (isstruct(model_oc), true, "svmtrain failed to return a valid struct model."); # <-- FIXED COMMENT HERE
%! # Predict with Probability (-b 1)
%! [pred, acc, probs] = svmpredict (L, D, model_oc, '-b 1');
%! 
%! # Detail Check A: Output must be N x 2 (Column 1: Normal, Column 2: Outlier)
%! assert (size (probs), [length(L), 2]);
%! 
%! # Detail Check B: Probabilities must sum to 1.0 for every instance
%! assert (sum (probs, 2), ones (length(L), 1), 1e-5);
%! 
%! # Detail Check C: Values must be valid probabilities [0, 1]
%! assert (all (all (probs >= 0 & probs <= 1)));
%! clear model_oc
%!
%!test
%! # Test 3: One-Class Decision Values (Standard Check)
%! # Verifies that the upgrade didn't break standard One-Class prediction (-b 0)
%! [L, D] = libsvmread (file_in_loadpath ("heart_scale.dat"));
%! model_oc = svmtrain (L, D, '-s 2 -n 0.1 -g 0.07');
%! [pred, acc, dec] = svmpredict (L, D, model_oc);
%! # Standard One-Class output is N x 1 (Scalar decision values)
%! assert (size (dec), [length(L), 1]);
%! clear model_oc
%!
%!shared L, D, model
%! # Test 4: Error Handling (Original Checks)
%! [L, D] = libsvmread (file_in_loadpath ("heart_scale.dat"));
%! model = svmtrain (L, D, '-c 1 -g 0.07');
%!
%!error <svmpredict: wrong number of output arguments.> ...
%! [p, a] = svmpredict (L, D, model);
%!error <svmpredict: wrong number of input arguments.> p = svmpredict (L, D);
%!error <svmpredict: label vector and instance matrix must be double.> ...
%! p = svmpredict (single (L), D, model);
%!error <svmpredict: model should be a struct array.> p = svmpredict (L, D, 123);
*/

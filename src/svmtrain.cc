/*
Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
Based on the Octave LIBSVM wrapper adapted by Alan Meeson (2014) based on an
earlier version of the LIBSVM (3.18) library for MATLAB. Current implementation
is based on LIBSVM 3.25 (2021) by Chih-Chung Chang and Chih-Jen Lin.

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
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

void print_null(const char *s) {}
//void print_string_octave(const char *s) {printf(s);}

// svm arguments
struct svm_parameter param;		// set by parse_command_line
struct svm_problem prob;		// set by read_problem
struct svm_model *model;
struct svm_node *x_space;
int cross_validation;
int nr_fold;


double do_cross_validation()
{
	int i;
	int total_correct = 0;
	double total_error = 0;
	double sumv = 0, sumy = 0, sumvv = 0, sumyy = 0, sumvy = 0;
	double *target = Malloc(double,prob.l);
	double retval = 0.0;

	svm_cross_validation(&prob,&param,nr_fold,target);
	if(param.svm_type == EPSILON_SVR ||
	   param.svm_type == NU_SVR)
	{
		for(i=0;i<prob.l;i++)
		{
			double y = prob.y[i];
			double v = target[i];
			total_error += (v-y)*(v-y);
			sumv += v;
			sumy += y;
			sumvv += v*v;
			sumyy += y*y;
			sumvy += v*y;
		}
		printf("Cross Validation Mean squared error = %g\n",total_error/prob.l);
		printf("Cross Validation Squared correlation coefficient = %g\n",
			((prob.l*sumvy-sumv*sumy)*(prob.l*sumvy-sumv*sumy))/
			((prob.l*sumvv-sumv*sumv)*(prob.l*sumyy-sumy*sumy))
			);
		retval = total_error/prob.l;
	}
	else
	{
		for(i=0;i<prob.l;i++)
			if(target[i] == prob.y[i])
				++total_correct;
		printf("Cross Validation Accuracy = %g%%\n",100.0*total_correct/prob.l);
		retval = 100.0*total_correct/prob.l;
	}
	free(target);
	return retval;
}

// nrhs should be 3
int parse_command_line(int nrhs, const octave_value_list args, char *model_file_name)
{
	int i, argc = 1;
	char cmd[CMD_LEN];
	char *argv[CMD_LEN/2];
	void (*print_func)(const char *) = print_null;

	// default values
	param.svm_type = C_SVC;
	param.kernel_type = RBF;
	param.degree = 3;
	param.gamma = 0;	// 1/num_features
	param.coef0 = 0;
	param.nu = 0.5;
	param.cache_size = 100;
	param.C = 1;
	param.eps = 1e-3;
	param.p = 0.1;
	param.shrinking = 1;
	param.probability = 0;
	param.nr_weight = 0;
	param.weight_label = NULL;
	param.weight = NULL;
	cross_validation = 0;

	if(nrhs <= 1)
		return 1;

	if(nrhs > 2)
	{
		// put options in argv[]
		strncpy(cmd, args(2).string_value().c_str(), CMD_LEN);
		//mxGetString(args[2], cmd, mxGetN(args[2]) + 1);
		if((argv[argc] = strtok(cmd, " ")) != NULL)
			while((argv[++argc] = strtok(NULL, " ")) != NULL)
				;
	}

	// parse options
	for(i=1;i<argc;i++)
	{
		if(argv[i][0] != '-') break;
		++i;
		if(i>=argc && argv[i-1][1] != 'q')	// since option -q has no parameter
			return 1;
		switch(argv[i-1][1])
		{
			case 's':
				param.svm_type = atoi(argv[i]);
				break;
			case 't':
				param.kernel_type = atoi(argv[i]);
				break;
			case 'd':
				param.degree = atoi(argv[i]);
				break;
			case 'g':
				param.gamma = atof(argv[i]);
				break;
			case 'r':
				param.coef0 = atof(argv[i]);
				break;
			case 'n':
				param.nu = atof(argv[i]);
				break;
			case 'm':
				param.cache_size = atof(argv[i]);
				break;
			case 'c':
				param.C = atof(argv[i]);
				break;
			case 'e':
				param.eps = atof(argv[i]);
				break;
			case 'p':
				param.p = atof(argv[i]);
				break;
			case 'h':
				param.shrinking = atoi(argv[i]);
				break;
			case 'b':
				param.probability = atoi(argv[i]);
				break;
			case 'q':
				print_func = &print_null;
				i--;
				break;
			case 'v':
				cross_validation = 1;
				nr_fold = atoi(argv[i]);
				if(nr_fold < 2)
				{
					printf("n-fold cross validation: n must >= 2\n");
					return 1;
				}
				break;
			case 'w':
				++param.nr_weight;
				param.weight_label = (int *)realloc(param.weight_label,sizeof(int)*param.nr_weight);
				param.weight = (double *)realloc(param.weight,sizeof(double)*param.nr_weight);
				param.weight_label[param.nr_weight-1] = atoi(&argv[i-1][2]);
				param.weight[param.nr_weight-1] = atof(argv[i]);
				break;
			default:
				printf("svmtrain: unknown option -%c\n", argv[i-1][1]);
				return 1;
		}
	}
	svm_set_print_string_function(print_func);
	return 0;
}

// read in a problem (in svmlight format)
int read_problem_dense(ColumnVector &label_vec, Matrix &instance_mat)
{
	int i, j, k;
	int elements, max_index, sc, label_vector_row_num;
	double *samples, *labels;

	prob.x = NULL;
	prob.y = NULL;
	x_space = NULL;

	labels = (double*)label_vec.mex_get_data();//mxGetPr(label_vec);
	samples = (double*)instance_mat.mex_get_data();
	sc = (int)instance_mat.cols();

	elements = 0;
	// the number of instance
	prob.l = (int)instance_mat.rows();
	label_vector_row_num = (int)label_vec.rows();

	if(label_vector_row_num!=prob.l)
	{
		printf("svmtrain: length of label vector does not match # of instances.\n");
		return -1;
	}

	if(param.kernel_type == PRECOMPUTED)
		elements = prob.l * (sc + 1);
	else
	{
		for(i = 0; i < prob.l; i++)
		{
			for(k = 0; k < sc; k++)
				if(samples[k * prob.l + i] != 0)
					elements++;
			// count the '-1' element
			elements++;
		}
	}

	prob.y = Malloc(double,prob.l);
	prob.x = Malloc(struct svm_node *,prob.l);
	x_space = Malloc(struct svm_node, elements);

	max_index = sc;
	j = 0;
	for(i = 0; i < prob.l; i++)
	{
		prob.x[i] = &x_space[j];
		prob.y[i] = labels[i];

		for(k = 0; k < sc; k++)
		{
			if(param.kernel_type == PRECOMPUTED || samples[k * prob.l + i] != 0)
			{
				x_space[j].index = k + 1;
				x_space[j].value = samples[k * prob.l + i];
				j++;
			}
		}
		x_space[j++].index = -1;
	}

	if(param.gamma == 0 && max_index > 0)
		param.gamma = 1.0/max_index;

	if(param.kernel_type == PRECOMPUTED)
		for(i=0;i<prob.l;i++)
		{
			if((int)prob.x[i][0].value <= 0 || (int)prob.x[i][0].value > max_index)
			{
				printf("svmtrain: wrong input format: sample_serial_number out of range\n");
				return -1;
			}
		}
	return 0;
}

int read_problem_sparse(ColumnVector &label_vec, SparseMatrix &instance_mat)
{
	int i, j, k, low, high;
	octave_idx_type *ir, *jc;
	int elements, max_index, num_samples, label_vector_row_num;
	double *samples, *labels;
  // transposed instance sparse matrix
	SparseMatrix instance_mat_col = instance_mat.transpose();

	prob.x = NULL;
	prob.y = NULL;
	x_space = NULL;

	// each column is one instance
	labels = (double*)label_vec.mex_get_data();
	samples = (double*)instance_mat_col.mex_get_data();
	ir = instance_mat_col.mex_get_ir();
	jc = instance_mat_col.mex_get_jc();

	num_samples = (int)instance_mat_col.nzmax();

	// the number of instance
	prob.l = (int)instance_mat_col.cols();
	label_vector_row_num = (int)label_vec.rows();

	if(label_vector_row_num!=prob.l)
	{
		printf("svmtrain: length of label vector does not match # of instances.\n");
		return -1;
	}

	elements = num_samples + prob.l;
	max_index = (int)instance_mat_col.rows();

	prob.y = Malloc(double,prob.l);
	prob.x = Malloc(struct svm_node *,prob.l);
	x_space = Malloc(struct svm_node, elements);

	j = 0;
	for(i=0;i<prob.l;i++)
	{
		prob.x[i] = &x_space[j];
		prob.y[i] = labels[i];
		low = (int)jc[i], high = (int)jc[i+1];
		for(k=low;k<high;k++)
		{
			x_space[j].index = (int)ir[k] + 1;
			x_space[j].value = samples[k];
			j++;
	 	}
		x_space[j++].index = -1;
	}

	if(param.gamma == 0 && max_index > 0)
  {
		param.gamma = 1.0/max_index;
  }
	return 0;
}

static void fake_answer(int nlhs, octave_value_list &plhs)
{
	int i;
	for(i=0;i<nlhs;i++) plhs(i) = Matrix(0,0);
}


DEFUN_DLD (svmtrain, args, nargout,
           "-*- texinfo -*-\n\
@deftypefn{Function} @var{model} = svmtrain (@var{labels}, @var{data}, 'libsvm_options')\n\
\n\
\n\
This function trains an SVM @var{model} based on known @var{labels} and their \
corresponding @var{data} which comprise an instance matrtix. \
\n\
\n\
---  @var{labels} : An m by 1 vector of prediction labels. (type must be double) \
\n\
\n\
---  @var{data} : An m by n matrix of m testing instances with n features. \
It can be dense or sparse. (type must be double) \
\n\
\n\
---  'libsvm_options' : A string of testing options in the same format as that \
of LIBSVM. \
\n\
\n\
* libsvm_options:\n\
\n\
--s svm_type : set type of SVM (default 0)\n\
\n\
    	0 -- C-SVC		(multi-class classification)\n\
\n\
    	1 -- nu-SVC		(multi-class classification)\n\
\n\
    	2 -- one-class SVM\n\
\n\
    	3 -- epsilon-SVR	(regression)\n\
\n\
    	4 -- nu-SVR		(regression)\n\
\n\
--t kernel_type : set type of kernel function (default 2)\n\
\n\
    	0 -- linear: u'*v\n\
\n\
    	1 -- polynomial: (gamma*u'*v + coef0)^degree\n\
\n\
    	2 -- radial basis function: exp(-gamma*|u-v|^2)\n\
\n\
    	3 -- sigmoid: tanh(gamma*u'*v + coef0)\n\
\n\
    	4 -- precomputed kernel (kernel values in training_instance_matrix)\n\
\n\
--d degree : set degree in kernel function (default 3)\n\
\n\
--g gamma : set gamma in kernel function (default 1/num_features)\n\
\n\
--r coef0 : set coef0 in kernel function (default 0)\n\
\n\
--c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)\n\
\n\
--n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)\n\
\n\
--p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)\n\
\n\
--m cachesize : set cache memory size in MB (default 100)\n\
\n\
--e epsilon : set tolerance of termination criterion (default 0.001)\n\
\n\
--h shrinking : whether to use the shrinking heuristics, 0 or 1 (default 1)\n\
\n\
--b probability_estimates : whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)\n\
\n\
--wi weight : set the parameter C of class i to weight*C, for C-SVC (default 1)\n\
\n\
--v n : n-fold cross validation mode\n\
\n\
--q : quiet mode (no outputs)\n\
\n\
\n\
The function 'svmtrain' function returns a @var{model} which can be used for \
future prediction.  It is a structure and is organized as [Parameters, \
nr_class, totalSV, rho, Label, ProbA, ProbB, nSV, sv_coef, SVs]: \
\n\
\n\
---  Parameters: parameters \
\n\
\n\
---  nr_class: number of classes; = 2 for regression/one-class svm \
\n\
\n\
---  totalSV: total #SV \
\n\
\n\
---  rho: -b of the decision function(s) wx+b \
\n\
\n\
---  Label: label of each class; empty for regression/one-class SVM \
\n\
\n\
---  sv_indices: values in [1,...,num_traning_data] to indicate SVs in the \
training set \
\n\
\n\
---  ProbA: pairwise probability information; empty if -b 0 or in one-class SVM \
\n\
\n\
---  ProbB: pairwise probability information; empty if -b 0 or in one-class SVM \
\n\
\n\
---  nSV: number of SVs for each class; empty for regression/one-class SVM \
\n\
\n\
---  sv_coef: coefficients for SVs in decision functions \
\n\
\n\
---  SVs: support vectors \
\n\
\n\
If you do not use the option '-b 1', ProbA and ProbB are empty \
matrices. If the '-v' option is specified, cross validation is \
conducted and the returned model is just a scalar: cross-validation \
accuracy for classification and mean-squared error for regression. \
\n\
\n\
@end deftypefn")
{
	const char *error_msg;
	octave_value_list plhs(nargout);
	// fix random seed to have same results for each run
	// (for cross validation and probability estimation)
	srand(1);
	int nlhs = nargout;
	int nrhs = args.length();
	if(nlhs > 1)
	{
    error ("svmtrain: wrong number of output arguments.");
	}

	// Transform the input Matrix to libsvm format
	if(nrhs > 1 && nrhs < 4)
	{
		int err;

		if(!args(0).is_double_type() || !args(1).is_double_type())
    {
			error ("svmtrain: label vector and instance matrix must be double.");
		}

		if(parse_command_line(nrhs, args, NULL))
		{
			svm_destroy_param(&param);
      error ("svmtrain: wrong values in parameter string.");
		}

		if(args(1).issparse())
		{
			if(param.kernel_type == PRECOMPUTED)
			{
				// precomputed kernel requires dense matrix, so we make one
        ColumnVector cv_lab = args(0).column_vector_value();
        Matrix m_dat = args(1).matrix_value();
				err = read_problem_dense(cv_lab, m_dat);
			}
			else {
        ColumnVector cv_lab = args(0).column_vector_value();
        SparseMatrix m_dat = args(1).sparse_matrix_value();
				err = read_problem_sparse(cv_lab, m_dat);
      }
		}
		else {
      ColumnVector cv_lab = args(0).column_vector_value();
      Matrix m_dat = args(1).matrix_value();
			err = read_problem_dense(cv_lab, m_dat);
    }

		// svmtrain's original code
		error_msg = svm_check_parameter(&prob, &param);

		if(err || error_msg)
		{
			if (error_msg != NULL)
      {
				printf("svmtrain: %s\n", error_msg);
      }
			svm_destroy_param(&param);
			free(prob.y);
			free(prob.x);
			free(x_space);
			fake_answer(nlhs, plhs);
			return plhs;
		}

		if(cross_validation)
		{
			double ptr = do_cross_validation();
			plhs(0) = octave_value(ptr);
		}
		else
		{
			int nr_feat = (int)args(1).matrix_value().cols();
			const char *error_msg;
			model = svm_train(&prob, &param);
			error_msg = model_to_octave_structure(plhs, nr_feat, model);
			if(error_msg)
      {
				printf("svmtrain: can't convert libsvm model to matrix structure: %s\n",
                error_msg);
      }
			svm_free_and_destroy_model(&model);
		}
		svm_destroy_param(&param);
		free(prob.y);
		free(prob.x);
		free(x_space);
    return plhs;
	}
	else
	{
    error ("svmtrain: wrong number of input arguments.");
	}
}

/*
%!test
%! [L, D] = libsvmread ("heart_scale.dat");
%! model = svmtrain(L, D, '-c 1 -g 0.07');
%! [predict_label, accuracy, dec_values] = svmpredict(L, D, model);
%! assert (isstruct (model), true);
%! assert (isfield (model, "Parameters"), true);
%! assert (model.totalSV, 130);
%! assert (model.nr_class, 2);
%! assert (size (model.Label), [2, 1]);
%!shared L, D
%! [L, D] = libsvmread ("heart_scale.dat");
%!error <svmtrain: wrong number of output arguments.> [L, D] = svmtrain (L, D);
%!error <svmtrain: label vector and instance matrix must be double.> model = svmtrain (single (L), D);
%!error <svmtrain: wrong number of input arguments.> model = svmtrain (L, D, "", "");
*/
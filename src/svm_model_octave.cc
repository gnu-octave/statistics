/*
Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
Based on the Octave LIBSVM wrapper created by Alan Meeson (2014) based on an
earlier version of the LIBSVM (3.18) library for MATLAB. Current implementation
is based on LIBSVM 3.25 (2021) by Chih-Chung Chang and Chih-Jen Lin.

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

#include <octave/oct.h>
#include <octave/ov-struct.h>
#include <octave/dMatrix.h>
#include <octave/dColVector.h>
#include <string.h>
#include "svm.h"

#define NUM_OF_RETURN_FIELD 12

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

static const char *field_names[] = {
	"Parameters",
	"nr_class",
	"totalSV",
	"rho",
	"Label",
	"sv_indices",
	"ProbA",
	"ProbB",
	"nSV",
	"sv_coef",
	"SVs",
	"ProbDensityMarks"
};

const char *model_to_octave_structure(octave_value_list &plhs, int num_of_feature, struct svm_model *model)
{
	int i, j, n;
	double *ptr;

	octave_scalar_map osm_model;

	// Parameters
	ColumnVector cm_parameters (5);
	cm_parameters(0) = model->param.svm_type;
	cm_parameters(1) = model->param.kernel_type;
	cm_parameters(2) = model->param.degree;
	cm_parameters(3) = model->param.gamma;
	cm_parameters(4) = model->param.coef0;
	osm_model.assign("Parameters", cm_parameters);

	// nr_class
	osm_model.assign("nr_class", octave_value(model->nr_class));

	// total SV
	osm_model.assign("totalSV", octave_value(model->l));

	// rho
	n = model->nr_class*(model->nr_class-1)/2;
	ColumnVector cm_rho(n);
	for(i = 0; i < n; i++)
  {
    cm_rho(i) = model->rho[i];
  }
	osm_model.assign("rho", cm_rho);

	// Label
	if(model->label)
  {
		ColumnVector cm_label(model->nr_class);
		for(i = 0; i < model->nr_class; i++)
    {
      cm_label(i) = model->label[i];
    }
		osm_model.assign("Label", cm_label);
	}
  else
  {
		osm_model.assign("Label", ColumnVector(0));
	}

	// sv_indices
	if(model->sv_indices)
	{
		ColumnVector cm_sv_indices(model->l);
		for(i = 0; i < model->l; i++)
    {
      cm_sv_indices(i) = model->sv_indices[i];
    }
		osm_model.assign("sv_indices", cm_sv_indices);
	}
  else
  {
		osm_model.assign("sv_indices", ColumnVector(0));
	}

	// probA
	if(model->probA != NULL)
	{
		ColumnVector cm_proba(n);
		for(i = 0; i < n; i++)
    {
      cm_proba(i) = model->probA[i];
    }
		osm_model.assign("ProbA", cm_proba);
	}
  else
  {
		osm_model.assign("ProbA", ColumnVector(0));
	}

	// probB
	if(model->probB != NULL)
	{
		ColumnVector cm_probb(n);
		for(i = 0; i < n; i++)
    {
      cm_probb(i) = model->probB[i];
    }
		osm_model.assign("ProbB", cm_probb);
	}
  else
  {
		osm_model.assign("ProbB", ColumnVector(0));
	}

	// nSV
	if(model->nSV)
	{
		ColumnVector cm_nsv(model->nr_class);
		for(i = 0; i < model->nr_class; i++)
    {
      cm_nsv(i) = model->nSV[i];
    }
		osm_model.assign("nSV", cm_nsv);
	}
  else
  {
		osm_model.assign("nSV", ColumnVector(0));
	}

	// sv_coef
	Matrix m_sv_coef(model->l, model->nr_class-1);
	for (i = 0; i < model->nr_class-1; i++)
  {
		for(j = 0; j < model->l; j++)
    {
			m_sv_coef(j,i) = model->sv_coef[i][j];
    }
  }
	osm_model.assign("sv_coef", m_sv_coef);

	// SVs
	{
		int ir_index, nonzero_element;
		octave_idx_type *ir, *jc;
		//mxArray *pprhs[1], *pplhs[1];

		if(model->param.kernel_type == PRECOMPUTED)
		{
			nonzero_element = model->l;
			num_of_feature = 1;
		}
		else
		{
			nonzero_element = 0;
			for(i = 0; i < model->l; i++)
      {
				j = 0;
				while(model->SV[i][j].index != -1)
				{
					nonzero_element++;
					j++;
				}
			}
		}

		// SV in column, easier accessing
		SparseMatrix sm_rhs = SparseMatrix((octave_idx_type)num_of_feature, (octave_idx_type)model->l, (octave_idx_type)nonzero_element);
		ir = sm_rhs.ridx();
		jc = sm_rhs.cidx();
		ptr = (double*) sm_rhs.data();
		jc[0] = ir_index = 0;
		for(i = 0; i < model->l; i++)
		{
			if(model->param.kernel_type == PRECOMPUTED)
			{
				// make a (1 x model->l) matrix
				ir[ir_index] = 0;
				ptr[ir_index] = model->SV[i][0].value;
				ir_index++;
				jc[i+1] = jc[i] + 1;
			}
			else
			{
				int x_index = 0;
				while (model->SV[i][x_index].index != -1)
				{
					ir[ir_index] = model->SV[i][x_index].index - 1;
					ptr[ir_index] = model->SV[i][x_index].value;
					ir_index++, x_index++;
				}
				jc[i+1] = jc[i] + x_index;
			}
		}
		// transpose back to SV in row
		sm_rhs = sm_rhs.transpose();
		osm_model.assign("SVs", sm_rhs);
	}

	// changes from libsvm 3.36
	if(model->prob_density_marks)
    {
        int nr_marks = 10; 
        Matrix m_marks(nr_marks, 1);
        for(int i = 0; i < nr_marks; i++)
        {
            m_marks(i) = model->prob_density_marks[i];
        }
        osm_model.setfield("ProbDensityMarks", m_marks);
    }
    else
    {
        osm_model.setfield("ProbDensityMarks", Matrix(0, 0));
    }

	/* return */
	plhs(0) = osm_model;
	return NULL;
}

struct svm_model *octave_matrix_to_model(octave_scalar_map &octave_model, const char **msg)
{
	int i, j, n, num_of_fields;
	double *ptr;
	int id = 0;
	struct svm_node *x_space;
	struct svm_model *model;

	model = Malloc(struct svm_model, 1);
	model->rho = NULL;
	model->probA = NULL;
	model->probB = NULL;
	model->label = NULL;
	model->sv_indices = NULL;
	model->nSV = NULL;
	model->free_sv = 1; // XXX

	//Parameters
	ColumnVector cm_parameters = octave_model.getfield("Parameters").column_vector_value();
	model->param.svm_type = (int)cm_parameters(0);
	model->param.kernel_type = (int)cm_parameters(1);
	model->param.degree = (int)cm_parameters(2);
	model->param.gamma = cm_parameters(3);
	model->param.coef0 = cm_parameters(4);

	//nr_class
	model->nr_class = (int)octave_model.getfield("nr_class").int_value();

	//total SV
	model->l = (int)octave_model.getfield("totalSV").int_value();

	//rho
	n = model->nr_class * (model->nr_class-1)/2;
	model->rho = (double*) malloc(n*sizeof(double));
	ColumnVector cm_rho = octave_model.getfield("rho").column_vector_value();
	for(i = 0; i < n; i++)
  {
    model->rho[i] = cm_rho(i);
  }

	//label
	if (!octave_model.getfield("Label").isempty())
  {
		model->label = (int*) malloc(model->nr_class*sizeof(int));
		ColumnVector cm_label = octave_model.getfield("Label").column_vector_value();
		for(i = 0; i < model->nr_class; i++)
    {
      model->label[i] = (int)cm_label(i);
    }
	}

	//sv_indices
	if (!octave_model.getfield("sv_indices").isempty())
  {
		model->sv_indices = (int*) malloc(model->l*sizeof(int));
		ColumnVector cv_svi = octave_model.getfield("sv_indices").column_vector_value();
		for(i = 0; i < model->l; i++)
    {
      model->sv_indices[i] = (int)cv_svi(i);
    }
	}

	// probA
	if(!octave_model.getfield("ProbA").isempty())
	{
		model->probA = (double*) malloc(n*sizeof(double));
		ColumnVector cv_proba = octave_model.getfield("ProbA").column_vector_value();
		for(i = 0; i < n; i++)
    {
      model->probA[i] = cv_proba(i);
    }
	}

	// probB
	if(!octave_model.getfield("ProbB").isempty())
	{
		model->probB = (double*) malloc(n*sizeof(double));
		ColumnVector cv_probb = octave_model.getfield("ProbB").column_vector_value();
		for(i = 0; i < n; i++)
    {
      model->probB[i] = cv_probb(i);
    }
	}

	// nSV
	if(!octave_model.getfield("nSV").isempty())
	{
		model->nSV = (int*) malloc(model->nr_class*sizeof(int));
		ColumnVector cv_nsv = octave_model.getfield("nSV").column_vector_value();
		for(i = 0; i < model->nr_class; i++)
    {
      model->nSV[i] = (int)cv_nsv(i);
    }
	}

	// sv_coef
	Matrix m_sv_coef = octave_model.getfield("sv_coef").matrix_value();
	ptr = (double*) m_sv_coef.data();
	model->sv_coef = (double**) malloc((model->nr_class-1)*sizeof(double));
	for(i = 0; i < model->nr_class - 1; i++ )
  {
		model->sv_coef[i] = (double*) malloc((model->l)*sizeof(double));
  }
	for(i = 0; i < model->nr_class - 1; i++)
  {
		for(j = 0; j < model->l; j++)
    {
			model->sv_coef[i][j] = ptr[i*(model->l)+j];//m_sv_coef(i,j);
    }
  }

	// SV
	{
		int sr, sc, elements;
		int num_samples;
		octave_idx_type *ir, *jc;

		// transpose SV
		SparseMatrix sm_sv = octave_model.getfield("SVs").sparse_matrix_value();
		sm_sv = sm_sv.transpose();

		sr = (int)sm_sv.cols();
		sc = (int)sm_sv.rows();

		ptr = (double*)sm_sv.data();
		ir = sm_sv.ridx();
		jc = sm_sv.cidx();

		num_samples = (int)sm_sv.nzmax();

		elements = num_samples + sr;

		model->SV = (struct svm_node **) malloc(sr * sizeof(struct svm_node *));
		x_space = (struct svm_node *)malloc(elements * sizeof(struct svm_node));

		// SV is in column
		for(i = 0; i < sr; i++)
		{
			int low = (int)jc[i], high = (int)jc[i+1];
			int x_index = 0;
			model->SV[i] = &x_space[low+i];
			for(j = low; j < high; j++)
			{
				model->SV[i][x_index].index = (int)ir[j] + 1;
				model->SV[i][x_index].value = ptr[j];
				x_index++;
			}
			model->SV[i][x_index].index = -1;
		}

		id++;
	}

	// changes from libsvm 3.36
	if (octave_model.isfield("ProbDensityMarks"))
	{
		Matrix m_marks = octave_model.getfield("ProbDensityMarks").matrix_value();
		if (m_marks.numel() > 0)
		{
			int nr_marks = 10;
			model->prob_density_marks = (double*) malloc(nr_marks * sizeof(double));
			for(int i = 0; i < nr_marks; i++)
			{
				model->prob_density_marks[i] = m_marks(i);
			}
		}
		else
		{
			model->prob_density_marks = NULL;
		}
	}
	else
	{
		model->prob_density_marks = NULL;
	}
	return model;
}

/*
Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
Adapted from MATLAB libsvmwrite.c file from the LIBSVM 3.25 (2021) library
by Chih-Chung Chang and Chih-Jen Lin.

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
#include <octave/dMatrix.h>
#include <octave/dColVector.h>
#include <octave/dRowVector.h>
#include <octave/ov.h>
#include <octave/ov-base.h>
#include <octave/ov-struct.h>
#include <fstream>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <string.h>

using namespace std;

void write(string filename, ColumnVector label_vec, SparseMatrix instance_mat)
{
	// open file
  FILE *fp = fopen(filename.c_str(),"w");
  if (fp == NULL)
  {
    error ("libsvmwrite: error opening file for write.");
  }
  else
  {
    // check for equal number of instances and labels
    int im_rows = (int)instance_mat.rows();
    int lv_rows = (int)label_vec.rows();
    if(im_rows != lv_rows)
    {
      // close file
      fclose (fp);
      remove (filename.c_str());
      error ("libsvmwrite: length of label vector does not match instances.");
    }
    // transpose instance sparse matrix in column format
    SparseMatrix instance_mat_col = instance_mat.transpose();

    octave_idx_type *ir, *jc, k, low, high;
    size_t i, l, label_vector_row_num;
    double *samples, *labels;
    // each column is one instance
    labels = (double*)label_vec.data();
    samples = (double*)instance_mat_col.data();
    ir = (octave_idx_type*)instance_mat_col.ridx();
    jc = (octave_idx_type*)instance_mat_col.cidx();
    for(int i = 0; i < lv_rows; i++)
    {
      fprintf(fp, "%.17g", labels[i]);
      low = jc[i], high = jc[i+1];

      for(k=low;k<high;k++)
      {
        fprintf(fp ," %lu:%g", (size_t)ir[k]+1, samples[k]);
      }
		  fprintf(fp, "\n");
    }
    // close file
    fclose(fp);
  }
}


DEFUN_DLD (libsvmwrite, args, nargout,
           "-*- texinfo -*-\n\
@deftypefn  {statistics} {} libsvmwrite (@var{filename}, @var{labels}, @var{data})\n\
\n\
\n\
This function saves the labels and the corresponding instance_matrix in a file \
specified by @var{filename}.  @var{data} must be a sparse matrix.  Both \
@var{labels}, @var{data} must be of double type. \
\n\
\n\
@end deftypefn")
{
	if(nargout > 0)
	{
		error ("libsvmwrite: wrong number of output arguments.");
		return octave_value_list();
	}

	// Transform the input Matrix to libsvm format
	if(args.length() == 3)
	{
		if(!args(1).is_double_type() || !args(2).is_double_type())
		{
			error ("libsvmwrite: label vector and instance matrix must be double.");
		}
    if(!args(0).is_string())
    {
      error ("libsvmwrite: filename must be a string.");
    }
    string filename = args(0).string_value();

		if(args(2).issparse())
    {
      ColumnVector label_vec = args(1).column_vector_value();
      SparseMatrix instance_mat = args(2).sparse_matrix_value();
			write(filename, label_vec, instance_mat);
    }
		else
		{
			error ("libsvmwrite: instance_matrix must be sparse.");
		}
	}
	else
	{
		error ("libsvmwrite: wrong number of input arguments.");
	}
  return octave_value_list();
}

/*
%!shared L, D
%! [L, D] = libsvmread (file_in_loadpath ("heart_scale.dat"));
%!error <libsvmwrite: error opening file for write.> libsvmwrite ("", L, D);
%!error <libsvmwrite: length of label vector does not match instances.> ...
%! libsvmwrite (tempname (), [L;L], D);
%!error <libsvmwrite: wrong number of output arguments.> ...
%! OUT = libsvmwrite (tempname (), L, D);
%!error <libsvmwrite: label vector and instance matrix must be double.> ...
%! libsvmwrite (tempname (), single (L), D);
%!error <libsvmwrite: filename must be a string.> libsvmwrite (13412, L, D);
%!error <libsvmwrite: instance_matrix must be sparse.> ...
%! libsvmwrite (tempname (), L, full (D));
%!error <libsvmwrite: wrong number of input arguments.> ...
%! libsvmwrite (tempname (), L, D, D);
*/

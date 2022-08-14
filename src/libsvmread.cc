/*
Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
Adapted from MATLAB libsvmwrite.c file from the LIBSVM 3.25 (2021) library
by Chih-Chung Chang and Chih-Jen Lin.

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
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

#ifndef max
#define max(x,y) (((x)>(y))?(x):(y))
#endif
#ifndef min
#define min(x,y) (((x)<(y))?(x):(y))
#endif

using namespace std;

typedef int64_t mwIndex;

static char *line;
static int max_line_len;

static char* readline(FILE *input)
{
	int len;

	if(fgets(line,max_line_len,input) == NULL)
  {
		return NULL;
  }
	while(strrchr(line,'\n') == NULL)
	{
		max_line_len *= 2;
		line = (char *) realloc(line, max_line_len);
		len = (int) strlen(line);
		if(fgets(line+len,max_line_len-len,input) == NULL)
    {
			break;
    }
	}
	return line;
}

// read the file in libsvm format
void read(string filename, octave_value_list &retval)
{
	int max_index, min_index, inst_max_index;
	size_t elements, k, i, l=0;
	FILE *fp = fopen(filename.c_str(),"r");
	char *endptr;
	octave_idx_type *ir, *jc;
	double *labels, *samples;

	if(fp == NULL)
	{
		printf("can't open input file %s\n",filename.c_str());
		return;
	}

	max_line_len = 1024;
	line = (char *) malloc(max_line_len*sizeof(char));

	max_index = 0;
	min_index = 1; // our index starts from 1
	elements = 0;
	while(readline(fp) != NULL)
	{
		char *idx, *val;
		// features
		int index = 0;

		inst_max_index = -1;
		strtok(line," \t");
		while (1)
		{
			idx = strtok(NULL,":");
			val = strtok(NULL," \t");
			if(val == NULL)
				break;

			errno = 0;
			index = (int) strtol(idx,&endptr,10);
			if(endptr == idx || errno != 0 || *endptr != '\0' || index <= inst_max_index)
			{
				printf("libsvmread: wrong input format at line %d.\n", (int)l+1);
				return;
			}
			else
				inst_max_index = index;

			min_index = min(min_index, index);
			elements++;
		}
		max_index = max(max_index, inst_max_index);
		l++;
	}
	rewind(fp);
  
	// y
	retval(0) = ColumnVector(l, 1);
	// x^T
	if (min_index <= 0)
  {
    octave_idx_type r = max_index-min_index+1;
    octave_idx_type c = l;
    octave_idx_type val = elements;
		retval(1) = SparseMatrix(r, c, val);
  }
	else
  {
    octave_idx_type r = max_index-min_index+1;
    octave_idx_type c = l;
    octave_idx_type val = elements;
		retval(1) = SparseMatrix(r, c, val);
  }
	labels = (double*)retval(0).mex_get_data();
	samples = (double*)retval(1).mex_get_data();
	ir = (octave_idx_type*)retval(1).mex_get_ir();
	jc = (octave_idx_type*)retval(1).mex_get_jc();

	k=0;
	for(i=0;i<l;i++)
	{
		char *idx, *val, *label;
		jc[i] = k;

		readline(fp);

		label = strtok(line," \t\n");
		if(label == NULL)
		{
			printf("libsvmread: empty line at line %d.\n", (int)i+1);
			return;
		}
		labels[i] = strtod(label,&endptr);
		if(endptr == label || *endptr != '\0')
		{
			printf("Wrong input format at line %d\n", (int)i+1);
			return;
		}

		// features
		while(1)
		{
			idx = strtok(NULL,":");
			val = strtok(NULL," \t");
			if(val == NULL)
      {
				break;
      }
      // precomputed kernel has <index> start from 0
			ir[k] = (mwIndex) (strtol(idx,&endptr,10) - min_index);

			errno = 0;
			samples[k] = strtod(val,&endptr);
			if (endptr == val || errno != 0 || (*endptr != '\0' && !isspace(*endptr)))
			{
				printf("libsvmread: wrong input format at line %d.\n", (int)i+1);
				return;
			}
			++k;
		}
	}
	jc[l] = k;

	fclose(fp);
	free(line);
  
  // transpose instance sparse matrix in row format
  SparseMatrix instance_mat_row = retval(1).sparse_matrix_value();
	retval(1) = instance_mat_row.transpose();
}


DEFUN_DLD (libsvmread, args, nargout,
           "-*- texinfo -*-\n\
@deftypefn{Function} [@var{labels}, @var{data}] = libsvmread (@var{filename})\n\
\n\
\n\
This function reads the labels and the corresponding instance_matrix from a \
LIBSVM data file and stores them in @var{labels} and @var{data} respectively. \
These can then be used as inputs to @code{svmtrain} or @code{svmpredict} \
function. \
\n\
\n\
@end deftypefn")
{
	if(args.length() != 1 || nargout != 2)
	{
		error ("libsvmread: wrong number of input or output arguments.");
	}
  if(!args(0).is_string())
  {
    error ("libsvmread: filename must be a string.");
  }
  string filename = args(0).string_value();
  octave_value_list retval(nargout);
	read(filename, retval);
	return retval;
}

/*
%!error <libsvmread: filename must be a string.> [L, D] = libsvmread (24);
%!error <libsvmread: wrong number of input or output arguments.> D = libsvmread ("filename");
%!test
%! [L, D] = libsvmread (file_in_loadpath ("heart_scale.dat"));
%! assert (size (L), [270, 1]);
%! assert (size (D), [270, 13]);
%!test
%! [L, D] = libsvmread (file_in_loadpath ("heart_scale.dat"));
%! assert (issparse (L), false);
%! assert (issparse (D), true);
*/
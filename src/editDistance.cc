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

#include <algorithm>
#include <iostream>
#include <cmath>
#include <cstring>
#include <fstream>
#include <vector>
#include <octave/oct.h>
#include <octave/Cell.h>
#include <octave/parse.h>

using namespace std;

struct UniqueVecOut
{
  ColumnVector IA;
  ColumnVector IC;
  Cell IA_c;
};

// Function for computing the minimum of three integer values
int minimum (int a, int b, int c)
{
	int min = a;
	if (b < min) min = b;
	if (c < min) min = c;
	return min;
}

// Function for computing the Levenshtein distance between two strings
int LevensDistStr (const string& s1, const string& s2)
{
  const size_t rows = s1.length();
	const size_t cols = s2.length();
  vector<int> curr(cols+1, 0);
  int prev;
  // Prepopulate 1st column
	for (size_t j = 0; j <= cols; j++)
	{
		curr[j] = j;
	}
  // Compute all other elements in distance matrix
	for (size_t i = 1; i <= rows; i++)
	{
    prev = curr[0];
    curr[0] = i;
		for (size_t j = 1; j <= cols; j++)
		{
      int temp = curr[j];
			if (s1[i - 1] == s2[j - 1])
			{
				curr[j] = prev;
			}
			else
			{
				curr[j] = 1 + minimum (prev, curr[j - 1], curr[j]);
			}
      prev = temp;
		}
	}
	return curr[cols];
}

// Function for computing the Levenshtein distance between two documents
int LevensDistDoc (const Cell& d1, const Cell& d2)
{
  const size_t rows = d1.numel();
	const size_t cols = d2.numel();
  vector<int> curr(cols+1, 0);
  int prev;
  // Prepopulate 1st column
	for (size_t j = 0; j <= cols; j++)
	{
		curr[j] = j;
	}
  // Compute all other elements in distance matrix
	for (size_t i = 1; i <= rows; i++)
	{
    prev = curr[0];
    curr[0] = i;
		for (size_t j = 1; j <= cols; j++)
		{
      int temp = curr[j];
			if (d1(i - 1).string_value() == d2(i - 1).string_value())
			{
				curr[j] = prev;
			}
			else
			{
				curr[j] = 1 + minimum (prev, curr[j - 1], curr[j]);
			}
      prev = temp;
		}
	}
	return curr[cols];
}

// Transform a distance triu matrix to a boolean triu matrix
boolMatrix double2bool (const Matrix& D, const int& minDist)
{
  const size_t sz = D.rows();
  boolMatrix Bmat(sz, sz);
  for (size_t i = 0; i < sz - 1; i++)
  {
    Bmat(i,i) = true;
    for (size_t j = i + 1; j < sz; j++)
    {
      if (D(i,j) <= minDist)
      {
        Bmat(i,j) = true;
        Bmat(j,i) = false;
      }
      else
      {
        Bmat(i,j) = false;
        Bmat(j,i) = false;
      }
    }
  }
  Bmat(sz - 1,sz - 1) = true;
  return Bmat;
}

// Transform a distance triu matrix to a distance vector
Matrix triu2Dvec (const Matrix& D)
{
  const size_t szA = D.rows();
  const size_t sz = szA * (szA - 1) / 2;
  octave_idx_type idx = 0;
  Matrix Dvec(sz, 1);
  for (size_t i = 0; i < szA - 1; i++)
  {
    for (size_t j = i + 1; j < szA; j++)
    {
      Dvec(idx++,0) = D(i,j);
    }
  }
  return Dvec;
}

// Compute unique indexing IA cell of vectors
vector<vector<int>> IAcellvec (const boolMatrix& B)
{
  const size_t rows = B.rows();
  const size_t cols = B.columns();
  vector<vector<int>> IAcell;
  for (size_t i = 0; i < rows; i++)
  {
    vector<int> IA_cIdx;
    IA_cIdx.push_back(i);
    for (size_t j = i + 1; j < cols; j++)
    {
      if (B(i,j))
      {
        IA_cIdx.push_back(j);
      }
    }
    IAcell.push_back(IA_cIdx);
  }
  return IAcell;
}

// Transform to IAcellvec to Cell
Cell IA2cell (const vector<vector<int>>& IAc, const vector<int>& IAv)
{
  const size_t sz_v = IAv.size();
  Cell IA(sz_v, 1);
  for (size_t i = 0; i < sz_v; i++)
  {
    int idx = IAv[i];
    const size_t sz_c = IAc[idx].size();
    Matrix IAidx(sz_c, 1);
    for (size_t j = 0; j < sz_c; j++)
    {
      IAidx(j,0) = IAc[idx][j] + 1;
    }
    IA(i,0) = IAidx;
  }
  return IA;
}

// Compute unique indexing IA vector
vector<int> IAvector (const vector<vector<int>>& IAcell)
{
  vector<int> IA;
  vector<int> IA_done;
  for (size_t i = 0; i < IAcell.size(); i++)
  {
    for (size_t j = 0; j < IAcell[i].size(); j++)
    {
      if (binary_search(IA_done.begin(), IA_done.end(), IAcell[i][j]))
      {
        break;
      }
      else
      {
        if (j == 0)
        {
          IA.push_back(IAcell[i][j]);
        }
        else
        {
          IA_done.push_back(IAcell[i][j]);
        }
      }
    }
    sort (IA_done.begin(), IA_done.end());
  }
  return IA;
}

// Transform to IAvector to Matrix
Matrix IA2mat (const vector<int>& IAv)
{
  const size_t sz_v = IAv.size();
  Matrix IA(sz_v, 1);
  for (size_t i = 0; i < sz_v; i++)
  {
    IA(i,0) = IAv[i] + 1;
  }
  return IA;
}

// Compute unique indexing IA vector
Matrix ICvector (const vector<vector<int>>& IAc, const size_t& szA)
{
  Matrix IC(szA, 1);
  vector<int> IC_done;
  for (size_t i = 0; i < IAc.size(); i++)
  {
    for (size_t j = 0; j < IAc[i].size(); j++)
    {
      if (binary_search(IC_done.begin(), IC_done.end(), IAc[i][j]))
      {
        break;
      }
      else
      {
        octave_idx_type idx = IAc[i][j];
        IC(idx,0) = i + 1;
        IC_done.push_back(IAc[i][j]);
      }
    }
  }
  return IC;
}

// Functionality for uniquetol
octave_value_list uniquetol (const int& nargout, const Cell& A, const Matrix& D,
                             const int& minDist, const bool& OutputAllIndices)
{
  octave_value_list retval (nargout);
  boolMatrix B = double2bool (D, minDist);
  vector<vector<int>> IAc = IAcellvec (B);
  vector<int> IAv = IAvector (IAc);
  // Build cellstr with unique elements
  Cell C(IAv.size(), 1);
  if (A.iscellstr())
  {
    for (size_t i = 0; i < IAv.size(); i++)
    {
      C(i,0) = A(IAv[i]).string_value();
    }
  }
  else
  {
    for (size_t i = 0; i < IAv.size(); i++)
    {
      C(i,0) = A.elem(IAv[i]);
    }
  }
  retval(0) = C;
  // Build IA vector output
  if (nargout > 1 && OutputAllIndices)
  {
    retval(1) = IA2cell (IAc, IAv);
  }
  else if (nargout > 1)
  {
    retval(1) = IA2mat (IAv);
  }
  // Build IC vector output
  if (nargout > 2)
  {
    retval(2) = ICvector (IAc, A.numel());
  }
  return retval;
}

// Expand a cell scalar to a cell vector
Cell expand (const Cell& IN, const size_t& sz)
{
  //octave_idx_type sz = static_cast<int>(sz_out);
  Cell OUT(sz, 1);
  for (size_t i = 0; i < sz; i++)
  {
    OUT(i,0) = IN.elem(0);
  }
  return OUT;
}

DEFUN_DLD(editDistance, args, nargout,
          "-*- texinfo -*-\n\
 @deftypefn  {statistics} {@var{d} =} editDistance (@var{str})\n\
 @deftypefnx {statistics} {@var{d} =} editDistance (@var{doc})\n\
 @deftypefnx {statistics} {@var{C} =} editDistance (@dots{}, @var{minDist})\n\
 @deftypefnx {statistics} {[@var{C}, @var{IA}, @var{IC}] =} editDistance @\
 (@dots{}, @var{minDist})\n\
 @deftypefnx {statistics} {[@var{C}, @var{IA}, @var{IC}] =} editDistance @\
 (@dots{}, @var{minDist}, @qcode{\"OutputAllIndices\"}, @var{value})\n\
 @deftypefnx {statistics} {@var{d} =} editDistance (@var{str1}, @var{str2})\n\
 @deftypefnx {statistics} {@var{d} =} editDistance (@var{doc1}, @var{doc2})\n\
\n\
\n\
Compute the edit (Levenshtein) distance between strings or documents. \
\n\n\
@code{@var{d} = editDistance (@var{str})} takes a cell array of character \
vectors and computes the Levenshtein distance between each pair of strings in \
@var{str} as the lowest number of grapheme insertions, deletions, and \
substitutions required to convert string @qcode{@var{str}@{1@}} to string \
@qcode{@var{str}@{2@}}.  If @var{str} is a @qcode{cellstr} vector with \
@math{N} elements, the returned distance @var{d} is an @math{(N * (N-1)) / 2)} \
column vector of doubles.  If @var{str} is an array (that is @code{all (size \
(str) > 1) = true}), then it is transformed to a column vector as in \
@code{str = str(:)}.  @code{editDistance} expects @var{str} to be a column \
vector, if it is row vector, it is transformed to a column vector.\n\n\
\
@code{@var{d} = editDistance (@var{doc})} can also take a cell array \
containing cell arrays of character vectors, in which case each element of \
@var{doc} is regarded as a document, and the character vector in each element \
of the cell string array is regarded a token.  @code{editDistance} computes \
the Levenshtein distance between each pair of cell elements in @var{doc} as \
the lowest number of token insertions, deletions, and substitutions required \
to convert document @qcode{@var{doc}@{1@}} to document @qcode{@var{doc}@{2@}}. \
If @var{doc} is a @qcode{cell} vector with @math{N} elements, the distance \
@var{d} is an @math{(N * (N-1)) / 2)} column vector of doubles.  If @var{doc} \
is an array (that is @code{all (size (doc) > 1) = true}), then it is converted \
to a column vector as in @code{doc = doc(:)}.\n\n\
\
@code{@var{C} = editDistance (@dots{}, @var{minDist})} specifies a minimum \
distance, @var{minDist}, which is regarded as a similarity threshold between \
each pair of strings or documents, defined in the previous syntaces.  In this \
case, @code{editDistance} resembles the functionality of the @code{uniquetol} \
function and returns the unique strings or documents that are similar up to \
@var{minDist} distance.  @var{C} is either a cellstring array or a cell array \
of cellstrings, depending on the first input argument.\n\n\
\
@code{[@var{C}, @var{IA}, @var{IC}] = editDistance (@dots{}, @var{minDist})} \
also returns index vectors @var{IA} and @var{IC}.  Assuming @var{A} contains \
either strings @var{str} or documents @var{doc} as defined above, @var{IA} \
is a column vector of indices to the first occurrence of similar elements such \
that @qcode{@var{C} = @var{A}(@var{IA})}, and @var{IC} is a column vector of \
indices such that @qcode{@var{A} ~ @var{C}(@var{IC})} where @qcode{~} means \
that the strings or documents are within the specified distance @var{minDist} \
of each other.\n\n\
\
@code{[@var{C}, @var{IA}, @var{IC}] = editDistance (@dots{}, @var{minDist}, \
@qcode{\"OutputAllIndices\"}, @var{value})} specifies the type of the second \
output index @var{IA}.  @var{value} must be a logical scalar.  When set to \
@code{true}, @var{IA} is a cell array containing the vectors of indices for \
ALL elements in @var{A} that are within the specified distance @var{minDist} \
of each other.  Each cell in @var{IA} corresponds to a value in @var{C} and \
the values in each cell correspond to locations in @var{A}.  If @var{value} is \
set to @code{false}, then @var{IA} is returned as an index vector described in \
the previous syntax.\n\n\
\
@code{@var{d} = editDistance (@var{str1}, @var{str2})} can also take two \
character vectors, @var{str1} and @var{str2} and compute the Levenshtein \
distance @var{d} as the lowest number of grapheme insertions, deletions, and \
substitutions required to convert @var{str1} to @var{str2}.  @var{str1} and \
@var{str2} may also be cellstring arrays, in which case the pairwise distance \
is computed between @qcode{@var{str1}@{n@}} and @qcode{@var{str1}@{n@}}.  The \
cellstring arrays must be of the same size or scalars, in which case the \
scalar is expanded to the size of the other cellstring input.  The returned \
distance @var{d} is a column vector with the same number of elements as the \
cellstring arrays.  If @var{str1}  or @var{str2} is an array, then it is \
transformed to a column vector.  @code{editDistance} expects both @var{str1} \
and @var{str2} to be a column vectors, if not, they are transformed into \
column vectors.\n\n\
\
@code{@var{d} = editDistance (@var{doc1}, @var{doc2})} can also take two cell \
array containing cell arrays of character vectors, in which case each element \
of @var{doc1} and @var{dos2} is regarded as a document, and the character \
vector in each element of the cell string array is regarded a token.  \
@code{editDistance} computes the pairwise Levenshtein distance between the \
of cell elements in @var{doc1} and @var{doc2} as the lowest number of token \
insertions, deletions, and substitutions required to convert document \
@qcode{@var{doc1}@{n@}} to document @qcode{@var{doc1}@{n@}}.\n\n\
@end deftypefn")
{
  int nargin = args.length();
  // Add default options
  bool OutputAllIndices = false;
  // Parse Name-Value paired arguments
  if (nargin > 2 && args(nargin-2).is_string())
  {
    string ParamName = "OutputAllIndices";
    if (args(nargin - 2).string_value() == ParamName)
    {
      const int idx = nargin - 1;
      if (args(idx).islogical() && args(idx).numel() == 1)
      {
        const boolMatrix tmp = args(idx).bool_matrix_value();
        OutputAllIndices = tmp(0);
      }
      else
      {
        error ("editDistance: value for OutputAllIndices "
               "must be a logical scalar.");
      }
      nargin--;
      nargin--;
    }
  }
  // Check for invalid number of input arguments
  if (nargin > 3)
  {
    error ("editDistance: too many input arguments.");
  }
  // Check for last argument being numeric (minDist)
  int minDist;
  bool doMinDist;
  if (nargin > 1 && args(nargin-1).isnumeric())
  {
    // Check minDist input argument
    if (args(nargin -1 ).numel() != 1)
    {
      error ("editDistance: minDist must be a scalar value.");
    }
    Matrix tmp = args(nargin - 1).matrix_value();
    if (tmp(0,0) < 0 || floor (tmp(0,0)) != tmp(0,0))
    {
      error ("editDistance: minDist must be a nonnegative integer.");
    }
    minDist = static_cast<int>(tmp(0,0));
    doMinDist = true;
    nargin--;
  }
  else
  {
    doMinDist = false;
  }
  // Check for invalid number of output arguments
  if ((nargout > 3 && doMinDist) || (nargout > 1 && ! doMinDist))
  {
    error ("editDistance: too many output arguments.");
  }
  // Check cases of string arguments
  octave_value_list retval (nargout);
  if (nargin == 1)
  {
    if (args(0).iscellstr())
    {
      // Get cellstr input argument
      const Cell strA = args(0).cellstr_value();
      size_t szA = strA.numel();
      // For scalar input return distance to itself, i.e. 0
      if (szA == 1)
      {
        retval(0) = double (0);
        return retval;
      }
      // Compute the edit distance
      Matrix D(szA, szA);
      #pragma omp parallel
      {
        #pragma omp parallel for
        for (size_t i = 0; i < szA - 1; i++)
        {
          D(i,i) = 0;
          string s1 = strA(i).string_value();
          for (size_t j = i + 1; j < szA; j++)
          {
            D(i,j) = LevensDistStr (s1, strA(j).string_value());
          }
        }
        D(szA - 1,szA - 1) = 0;
      }
      // If minDist is given, change functionality from 'pdist' to 'uniquetol'
      if (doMinDist)
      {
        retval = uniquetol (nargout, strA, D, minDist, OutputAllIndices);
      }
      else
      {
        // Transform to distance vector
        retval(0) = triu2Dvec (D);
      }
      return retval;
    }
    else if (args(0).iscell())
    {
      // Get cell input argument
      const Cell docA = args(0).cell_value();
      size_t szA = docA.numel();
      // Check that all cell elements contain cellstring arrays
      for (size_t i = 0; i < szA; i++)
      {
        Cell tmp = docA.elem(i);
        if (! tmp.iscellstr())
        {
          error ("editDistance: tokenizedDocument "
                 "must contain cellstr arrays.");
        }
      }
      // For scalar input return distance to itself, i.e. 0
      if (szA == 1)
      {
        retval(0) = double (0);
        return retval;
      }
      // Compute the edit distance
      Matrix D(szA, szA);
      #pragma omp parallel
      {
        #pragma omp parallel for
        for (size_t i = 0; i < szA - 1; i++)
        {
          D(i,i) = 0;
          Cell d1 = docA.elem(i);
          for (size_t j = i + 1; j < szA; j++)
          {
            D(i,j) = LevensDistDoc (d1, docA.elem(j));
          }
        }
        D(szA - 1,szA - 1) = 0;
      }
      // If minDist is given, change functionality from 'pdist' to 'uniquetol'
      if (doMinDist)
      {
        retval = uniquetol (nargout, docA, D, minDist, OutputAllIndices);
      }
      else
      {
        // Transform to distance vector
        retval(0) = triu2Dvec (D);
      }
      return retval;
    }
    else
    {
      error ("editDistance: STR1 must be a cellstr.");
    }
  }
  else if (nargin == 2)
  {
    if (args(0).iscellstr() && args(1).iscellstr())
    {
      // Get cellstr input arguments
      Cell strA = args(0).cellstr_value();
      Cell strB = args(1).cellstr_value();
      // Check cellstr sizes match
      size_t szA = strA.numel();
      size_t szB = strB.numel();
      if (szA != 1 && szB != 1 && szA != szB)
      {
        error ("editDistance: cellstr input arguments size mismatch.");
      }
      // Preallocate the distance vector and expand as necessary
      size_t sz = szA;
      if (szA == 1 && szB != 1)
      {
        sz = szB;
        strA = expand (strA, sz);
      }
      else if (szA != 1 && szB == 1)
      {
        strB = expand (strB, sz);
      }
      Matrix D(sz, 1);
      // Compute the distance vector
      for (size_t i = 0; i < sz; i++)
      {
        D(i,0) = LevensDistStr (strA(i).string_value(), strB(i).string_value());
      }
      retval(0) = D;
    }
    else if (args(0).iscell() && args(1).iscell())
    {
      // Get cell input arguments
      Cell docA = args(0).cell_value();
      Cell docB = args(1).cell_value();
      // Check cell sizes match
      size_t szA = docA.numel();
      size_t szB = docB.numel();
      if (szA != 1 && szB != 1 && szA != szB)
      {
        error ("editDistance: cellstr input arguments size mismatch.");
      }
      // Check both cell arrays contain cellstring arrays
      for (size_t i = 0; i < szA; i++)
      {
        Cell tmp = docA.elem(i);
        if (! tmp.iscellstr())
        {
          error ("editDistance: first tokenizedDocument "
                 "does not contain cellstr arrays.");
        }
      }
      for (size_t i = 0; i < szB; i++)
      {
        Cell tmp = docB.elem(i);
        if (! tmp.iscellstr())
        {
          error ("editDistance: second tokenizedDocument "
                 "does not contain cellstr arrays.");
        }
      }
      // Preallocate the distance vector and expand as necessary
      int sz = szA;
      if (szA == 1 && szB != 1)
      {
        sz = szB;
        docA = expand (docA, sz);
      }
      else if (szA != 1 && szB == 1)
      {
        docB = expand (docB, sz);
      }
      Matrix D(sz, 1);
      // Compute the distance vector
      for (size_t i = 0; i < sz; i++)
      {
        D(i,0) = LevensDistDoc (docA.elem(i), docB.elem(i));
      }
      retval(0) = D;
    }
    else if (args(0).is_string() && args(1).is_string())
    {
      retval(0) = LevensDistStr (args(0).string_value(),args(1).string_value());
    }
    else
    {
      error ("editDistance: STR1 and STR2 must be either strings or cellstr.");
    }
  }
  return retval;
}

/*
%!error <editDistance: too many input arguments.> d = editDistance (1, 2, 3, 4);
%!error <editDistance: too many output arguments.> ...
%! [C, IA, IC, I] = editDistance ({"AS","SD","AD"}, 1);
%!error <editDistance: too many output arguments.> ...
%! [C, IA] = editDistance ({"AS","SD","AD"});
%!error <editDistance: minDist must be a scalar value.> ...
%! d = editDistance ({"AS","SD","AD"}, [1, 2]);
%!error <editDistance: minDist must be a nonnegative integer.> ...
%! d = editDistance ({"AS","SD","AD"}, -2);
%!error <editDistance: minDist must be a nonnegative integer.> ...
%! d = editDistance ({"AS","SD","AD"}, 1.25);
%!error <editDistance: minDist must be a scalar value.> ...
%! d = editDistance ({"AS","SD","AD"}, {"AS","SD","AD"}, [1, 2]);
%!error <editDistance: minDist must be a nonnegative integer.> ...
%! d = editDistance ({"AS","SD","AD"}, {"AS","SD","AD"}, -2);
%!error <editDistance: minDist must be a nonnegative integer.> ...
%! d = editDistance ({"AS","SD","AD"}, {"AS","SD","AD"}, 1.25);
%!error <editDistance: minDist must be a scalar value.> ...
%! d = editDistance ("string1", "string2", [1, 2]);
%!error <editDistance: minDist must be a nonnegative integer.> ...
%! d = editDistance ("string1", "string2", -2);
%!error <editDistance: minDist must be a nonnegative integer.> ...
%! d = editDistance ("string1", "string2", 1.25);
%!error <editDistance: tokenizedDocument must contain cellstr arrays.> ...
%! d = editDistance ({{"string1", "string2"}, 2});
%!error <editDistance: tokenizedDocument must contain cellstr arrays.> ...
%! d = editDistance ({{"string1", "string2"}, 2}, 2);
%!error <editDistance: STR1 must be a cellstr.> ...
%! d = editDistance ([1, 2, 3]);
%!error <editDistance: STR1 must be a cellstr.> ...
%! d = editDistance (["AS","SD","AD","AS"]);
%!error <editDistance: STR1 must be a cellstr.> ...
%! d = editDistance (["AS","SD","AD"], 2);
%!error <editDistance: STR1 and STR2 must be either strings or cellstr.> ...
%! d = editDistance (logical ([1,2,3]), {"AS","AS","AD"});
%!error <editDistance: STR1 and STR2 must be either strings or cellstr.> ...
%! d = editDistance ({"AS","SD","AD"}, logical ([1,2,3]));
%!error <editDistance: STR1 and STR2 must be either strings or cellstr.> ...
%! d = editDistance ([1,2,3], {"AS","AS","AD"});
%!error <editDistance: first tokenizedDocument does not contain cellstr arrays.> ...
%! d = editDistance ({1,2,3}, {"AS","SD","AD"});
%!error <editDistance: second tokenizedDocument does not contain cellstr arrays.> ...
%! d = editDistance ({"AS","SD","AD"}, {1,2,3});
%!error <editDistance: cellstr input arguments size mismatch.> ...
%! d = editDistance ({"AS","SD","AD"}, {"AS", "AS"});
%!test
%! d = editDistance ({"AS","SD","AD"});
%! assert (d, [2; 1; 1]);
%! assert (class (d), "double");
%!test
%! C = editDistance ({"AS","SD","AD"}, 1);
%! assert (iscellstr (C), true);
%! assert (C, {"AS";"SD"});
%!test
%! [C, IA] = editDistance ({"AS","SD","AD"}, 1);
%! assert (class (IA), "double");
%! assert (IA, [1;2]);
%!test
%! A = {"ASS"; "SDS"; "FDE"; "EDS"; "OPA"};
%! [C, IA] = editDistance (A, 2, "OutputAllIndices", false);
%! assert (class (IA), "double");
%! assert (A(IA), C);
%!test
%! A = {"ASS"; "SDS"; "FDE"; "EDS"; "OPA"};
%! [C, IA] = editDistance (A, 2, "OutputAllIndices", true);
%! assert (class (IA), "cell");
%! assert (C, {"ASS"; "FDE"; "OPA"});
%! assert (A(IA{1}), {"ASS"; "SDS"; "EDS"});
%! assert (A(IA{2}), {"FDE"; "EDS"});
%! assert (A(IA{3}), {"OPA"});
%!test
%! A = {"ASS"; "SDS"; "FDE"; "EDS"; "OPA"};
%! [C, IA, IC] = editDistance (A, 2);
%! assert (class (IA), "double");
%! assert (A(IA), C);
%! assert (IC, [1; 1; 3; 1; 5]);
%!test
%! d = editDistance ({"AS","SD","AD"}, {"AS", "AD", "SE"});
%! assert (d, [0; 1; 2]);
%! assert (class (d), "double");
%!test
%! d = editDistance ({"AS","SD","AD"}, {"AS"});
%! assert (d, [0; 2; 1]);
%! assert (class (d), "double");
%!test
%! d = editDistance ({"AS"}, {"AS","SD","AD"});
%! assert (d, [0; 2; 1]);
%! assert (class (d), "double");
%!test
%! b = editDistance ("Octave", "octave");
%! assert (b, 1);
%! assert (class (b), "double");
*/

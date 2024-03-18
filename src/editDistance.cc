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
  const int rows = s1.length();
	const int cols = s2.length();
  vector<int> curr(cols+1, 0);
  int prev;
  // Prepopulate 1st column
	for (int j = 0; j <= cols; j++)
	{
		curr[j] = j;
	}
  // Compute all other elements in distance matrix
	for (int i = 1; i <= rows; i++)
	{
    prev = curr[0];
    curr[0] = i;
		for (int j = 1; j <= cols; j++)
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
  const int rows = d1.numel();
	const int cols = d2.numel();
  vector<int> curr(cols+1, 0);
  int prev;
  // Prepopulate 1st column
	for (int j = 0; j <= cols; j++)
	{
		curr[j] = j;
	}
  // Compute all other elements in distance matrix
	for (int i = 1; i <= rows; i++)
	{
    prev = curr[0];
    curr[0] = i;
		for (int j = 1; j <= cols; j++)
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

// Transform a distance triu matric to a boolean triu matrix
boolMatrix double2bool (const Matrix& D, const int& minDist)
{
  const int sz = D.rows();
  boolMatrix Bmat(sz, sz);
  for (octave_idx_type i = 0; i < sz - 1; i++)
  {
    Bmat(i,i) = true;
    for (octave_idx_type j = i + 1; j < sz; j++)
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

// Transform a distance triu matric to a distance vector
Matrix triu2Dvec (const Matrix& D)
{
  const int szA = D.rows();
  const int sz = szA * (szA - 1) / 2;
  octave_idx_type idx = 0;
  Matrix Dvec(sz, 1);
  for (octave_idx_type i = 0; i < szA - 1; i++)
  {
    for (octave_idx_type j = i + 1; j < szA; j++)
    {
      Dvec(idx++,0) = D(i,j);
    }
  }
  return Dvec;
}

// Compute unique indexing IA cell of vectors
vector<vector<int>> IAcellvec (const boolMatrix& B)
{
  int rows = B.rows();
  int cols = B.columns();
  vector<vector<int>> IAcell;
  for (int i = 0; i < rows; i++)
  {
    vector<int> IA_cIdx;
    IA_cIdx.push_back(i);
    for (int j = i + 1; j < cols; j++)
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

// Compute unique indexing IA vector
vector<int> IAvector (const vector<vector<int>>& IAcell)
{
  vector<int> IA;
  vector<int> IA_done;
  for (int i = 0; i < IAcell.size(); i++)
  {
    for (int j = 0; j < IAcell[i].size(); j++)
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

DEFUN_DLD(editDistance, args, nargout,
          "-*- texinfo -*-\n\
 @deftypefn  {statistics} {@var{d} =} editDistance (@var{str1})\n\
 @deftypefnx {statistics} {@var{d} =} editDistance (@var{doc1})\n\
 @deftypefnx {statistics} {@var{C} =} editDistance (@dots{}, @var{minDist})\n\
 @deftypefnx {statistics} {[@var{C}, @var{IA}, @var{IC}] =} editDistance @\
 (@dots{}, @var{minDist})\n\
 @deftypefnx {statistics} {[@var{C}, @var{IA}, @var{IC}] =} editDistance @\
 (@dots{}, @var{minDist}, @qcode{\"OutputAllIndices\"}, @qcode{true})\n\
 @deftypefnx {statistics} {[@var{C}, @var{IA}, @var{IC}] =} editDistance @\
 (@dots{}, @var{minDist}, @qcode{\"OutputAllIndices\"}, @qcode{false})\n\
 @deftypefnx {statistics} {@var{d} =} editDistance (@var{str1}, @var{str2})\n\
 @deftypefnx {statistics} {@var{d} =} editDistance (@dots{}, @var{minDist})\n\
\n\
\n\
Compute the edit (Levenshtein) distance between character vectors. \
\n\n\
@code{@var{d} = editDistance (@var{str1})} takes a cell array of character \
vectors and computes the Levenshtein distance between each pair of elements in \
@var{str1}.  If @var{str1} is a @qcode{cellstr} vector with @math{N} elements, \
the distance @var{d} is an @math{(N * (N-1)) / 2)} column vector of doubles. \
If @var{str1} is an array (that is @code{all (size (str1) > 1) = true}), then \
it is converted to a column vector as in @code{str1 = str1(:)}. \n\n\
@code{@var{d} = editDistance (@var{str1}, @var{str2})} can also take two \
character vectors, @var{str1} and @var{str2}, and compute the Levenshtein \
distance between them or two cell arrays of character vectors and compute the \
pairwise distance between them.  In the latter case, @var{str1} and @var{str2} \
must have the same number of cell elements.  If @var{str1} and/or @var{str2} \
are arrays (that is @code{all (size (str1) > 1) = true}), then they are \
converted to a column vector as in @code{str1 = str1(:)}. \n\n\
@code{@var{d} = editDistance (@dots{}, @var{minDist})} can also take an \
optional argument, @var{minDist}, which defines an upper threshold of \
Levenshtein distance to compare with.  This optional argument supports both \
previous syntaxes, but the returned argument @var{d} is a logical vector \
identifying as @qcode{true} the pairs that DO NOT exceed @var{minDist}. \n\n\
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
        error ("editDistance: value for OutputAllIndices must be a logical scalar.");
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
  // Check cases of string arguments
  octave_value_list retval (nargout);
  if (nargin == 1)
  {
    if (args(0).iscellstr())
    {
      // Get cellstr input argument
      const Cell strA = args(0).cellstr_value();
      int szA = strA.numel();
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
        for (octave_idx_type i = 0; i < szA - 1; i++)
        {
          D(i,i) = 0;
          string s1 = strA(i).string_value();
          for (octave_idx_type j = i + 1; j < szA; j++)
          {
            D(i,j) = LevensDistStr (s1, strA(j).string_value());
          }
        }
        D(szA - 1,szA - 1) = 0;
      }
      // If minDist is given, revert functionality from 'pdist' to 'uniquetol'
      // and return a vector indexing similar strings
      if (doMinDist)
      {
        boolMatrix B = double2bool (D, minDist);
        vector<vector<int>> IAc;
        vector<int> IAv;
        IAc = IAcellvec (B);
        IAv = IAvector (IAc);
        // Build cellstr with unique elements
        Cell C(IAv.size(), 1);
        for (octave_idx_type i = 0; i < IAv.size(); i++)
        {
          C(i,0) = strA(IAv[i]).string_value();
        }
        retval(0) = C;
        // Build IA vector output
        if (nargout > 1)
        {
          if (OutputAllIndices)
          {
            Cell IA(IAv.size(), 1);
            for (octave_idx_type i = 0; i < IAv.size(); i++)
            {
              int idx = IAv[i];
              Matrix IAidx(IAc[idx].size(), 1);
              for (octave_idx_type j = 0; j < IAc[idx].size(); j++)
              {
                IAidx(j,0) = IAc[idx][j] + 1;
              }
              IA(i,0) = IAidx;
            }
            retval(1) = IA;
          }
          else
          {
            Matrix IA(IAv.size(), 1);
            for (octave_idx_type i = 0; i < IAv.size(); i++)
            {
              IA(i,0) = IAv[i] + 1;
            }
            retval(1) = IA;
          }
        }
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
      // Get cellstr input argument
      const Cell docA = args(0).cell_value();
      int szA = docA.numel();
      // Check that all cell elements contain cellstr arrays
      for (octave_idx_type i = 0; i < szA - 1; i++)
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
        for (octave_idx_type i = 0; i < szA - 1; i++)
        {
          D(i,i) = 0;
          Cell d1 = docA.elem(i);
          for (octave_idx_type j = i + 1; j < szA; j++)
          {
            D(i,j) = LevensDistDoc (d1, docA.elem(j));
          }
        }
        D(szA - 1,szA - 1) = 0;
      }
      // If minDist is given, revert functionality from 'pdist' to 'uniquetol'
      // and return a vector indexing similar documents
      if (doMinDist)
      {
        boolMatrix B = double2bool (D, minDist);
        vector<vector<int>> IAc;
        vector<int> IAv;
        IAc = IAcellvec (B);
        IAv = IAvector (IAc);
        // Build cellstr with unique elements
        Cell C(IAv.size(), 1);
        for (octave_idx_type i = 0; i < IAv.size(); i++)
        {
          C(i,0) = docA(IAv[i]).cell_value();
        }
        retval(0) = C;
        // Build IA vector output
        if (nargout > 1)
        {
          if (OutputAllIndices)
          {
            Cell IA(IAv.size(), 1);
            for (octave_idx_type i = 0; i < IAv.size(); i++)
            {
              int idx = IAv[i];
              Matrix IAidx(IAc[idx].size(), 1);
              for (octave_idx_type j = 0; j < IAc[idx].size(); j++)
              {
                IAidx(j,0) = IAc[idx][j] + 1;
              }
              IA(i,0) = IAidx;
            }
            retval(1) = IA;
          }
          else
          {
            Matrix IA(IAv.size(), 1);
            for (octave_idx_type i = 0; i < IAv.size(); i++)
            {
              IA(i,0) = IAv[i] + 1;
            }
            retval(1) = IA;
          }
        }
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
      const Cell strA = args(0).cellstr_value();
      const Cell strB = args(1).cellstr_value();
      // Check cellstr sizes match
      int szA = strA.numel();
      int szB = strB.numel();
      if (szA != 1 && szB != 1 && szA != szB)
      {
        error ("editDistance: cellstr input arguments size mismatch.");
      }
      // Preallocate the distance vector
      int szV = szA;
      if (szA == 1)
      {
        szV = szB;
      }
      Matrix D(szV, 1);
      // Compute the distance vector
      if (szA == 1 && szB != 1)
      {
        for (octave_idx_type i = 0; i < szB; i++)
        {
          D(i,0) = LevensDistStr (strA(0).string_value(),
                                  strB(i).string_value());
        }
        retval(0) = D;
      }
      else if (szA != 1 && szB == 1)
      {
        for (octave_idx_type i = 0; i < szA; i++)
        {
          D(i,0) = LevensDistStr (strA(i).string_value(),
                                  strB(0).string_value());
        }
        retval(0) = D;
      }
      else
      {
        for (octave_idx_type i = 0; i < szA; i++)
        {
          D(i,0) = LevensDistStr (strA(i).string_value(),
                                  strB(i).string_value());
        }
        retval(0) = D;
      }
      if (doMinDist)
      {
        boolMatrix Bvec(szV, 1);
        for (octave_idx_type i = 0; i < szV; i++)
        {
          if (D(i,0) <= minDist)
          {
            Bvec(i,0) = true;
          }
          else
          {
            Bvec(i,0) = false;
          }
        }
        retval(0) = Bvec;
      }
      else
      {
        retval(0) = D;
      }
    }
    else if (args(0).is_string() && args(1).is_string())
    {
      int d = LevensDistStr (args(0).string_value(), args(1).string_value());
      if (doMinDist)
      {
        if (d <= minDist)
        {
          retval(0) = true;
        }
        else
        {
          retval(0) = false;
        }
      }
      else
      {
        retval(0) = d;
      }
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
%!error <editDistance: STR1 and STR2 must be either strings or cellstr.> ...
%! d = editDistance ({"AS","SD","AD"}, {1,2,3});
%!error <editDistance: cellstr input arguments size mismatch.> ...
%! d = editDistance ({"AS","SD","AD"}, {"AS", "AS"});
%!error <editDistance: cellstr input arguments size mismatch.> ...
%! d = editDistance ({"AS","SD","AD"}, {"AS","SD","AD","AF"}, 2);
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
%! b = editDistance ({"AS","SD","AD"}, {"AS", "AD", "SE"}, 1);
%! assert (b, logical ([1; 1; 0]));
%! assert (class (b), "logical");
%!test
%! b = editDistance ({"AS","SD","AD"}, {"AS"}, 1);
%! assert (b, logical ([1; 0; 1]));
%! assert (class (b), "logical");
%!test
%! b = editDistance ({"AS"}, {"AS", "AD", "SE"}, 1);
%! assert (b, logical ([1; 1; 0]));
%! assert (class (b), "logical");
%!test
%! b = editDistance ("Octave", "octave");
%! assert (b, 1);
%! assert (class (b), "double");
%!test
%! b = editDistance ("Octave", "octave", 1);
%! assert (b, true);
%! assert (class (b), "logical");
*/

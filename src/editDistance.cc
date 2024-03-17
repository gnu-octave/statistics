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
#include <cmath>
#include <cstring>
#include <fstream>
#include <vector>
#include <octave/oct.h>
#include <octave/Cell.h>
#include <octave/parse.h>

using namespace std;

// Function for computing the minimum of three integer values
int minimum (int a, int b, int c)
{
	int min = a;
	if (b < min) min = b;
	if (c < min) min = c;
	return min;
}


// Function for computing the Levenshtein distance
int LevensDist (const string& s1, const string s2)
{
  const int rows = s1.length() + 1;
	const int cols = s2.length() + 1;
  // Dynamically allocate distance matrix
	int **dist;
  dist = new int*[rows];
  for (int i = 0; i <rows; i++)
  {
    dist[i] = new int[cols];
  }
  // Prepopulate 1st row and column
	for (int i = 0; i < rows; i++)
	{
		dist[i][0] = i;
	}
	for (int j = 0; j < cols; j++)
	{
		dist[0][j] = j;
	}
  // Compute all other elements in distance matrix
	for (int i = 1; i < rows; i++)
	{
		for (int j = 1; j < cols; j++)
		{
			if (s1[i - 1] == s2[j - 1])
			{
				dist[i][j] = dist[i - 1][j - 1];
			}
			else
			{
				dist[i][j] = minimum (dist[i - 1][j - 1] + 1, dist[i - 1][j] + 1,
						                  dist[i][j - 1] + 1);
			}
		}
	}
  int d = dist[rows - 1][cols - 1];
  // Free memory
  for(int i = 0; i < rows; i++)
  {
    delete dist[i];
  }
  delete dist;
	return d;
}

// Function for comparing the Levenshtein distance against a minimum distance
bool minLevensDist (const string& s1, const string& s2, int const minDist)
{
  int rows = s1.length() + 1;
	int cols = s2.length() + 1;
  // Dynamically allocate distance matrix
	int **dist;
  dist = new int*[rows];
  for (int i = 0; i <rows; i++)
  {
    dist[i] = new int[cols];
  }
  // Prepopulate 1st row and column
	for (int i = 0; i < rows; i++)
	{
		dist[i][0] = i;
	}
	for (int j = 0; j < cols; j++)
	{
		dist[0][j] = j;
	}
  // Compute all other elements in distance matrix until minDist is exceeded
	for (int i = 1; i < rows; i++)
	{
		for (int j = 1; j < cols; j++)
		{
			if (s1[i - 1] == s2[j - 1])
			{
				dist[i][j] = dist[i - 1][j - 1];
			}
			else
			{
				dist[i][j] = minimum (dist[i - 1][j - 1] + 1, dist[i - 1][j] + 1,
						                  dist[i][j - 1] + 1);
			  int d = dist[rows - 1][cols - 1];
        if (d > minDist)
        {
          // Free memory
          for(int i = 0; i < rows; i++)
          {
            delete dist[i];
          }
          delete dist;
          return false;
        }
      }
		}
	}
  // Free memory
  for(int i = 0; i < rows; i++)
  {
    delete dist[i];
  }
  delete dist;
	return true;
}

DEFUN_DLD(editDistance, args, nargout,
          "-*- texinfo -*-\n\
 @deftypefn  {statistics} {@var{d} =} editDistance (@var{str1})\n\
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
    if (args(nargin-1).numel() != 1)
    {
      error ("editDistance: minDist must be a scalar value.");
    }
    Matrix tmp = args(nargin-1).matrix_value();
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
  octave_value_list retval;
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
      // Compute the distance vector
      int sz = szA * (szA - 1) / 2;
      octave_idx_type idx = 0;
      if (doMinDist)
      {
        boolMatrix D(sz, 1);
        #pragma omp parallel for
        for (octave_idx_type i = 0; i < szA - 1; i++)
        {
          string s1 = strA(i).string_value();
          for (octave_idx_type j = i + 1; j < szA; j++)
          {
            D(idx++,0) = minLevensDist (s1, strA(j).string_value(), minDist);
          }
        }
        retval(0) = D;
      }
      else
      {
        Matrix D(sz, 1);
        #pragma omp parallel for
        for (octave_idx_type i = 0; i < szA - 1; i++)
        {
          string s1 = strA(i).string_value();
          for (octave_idx_type j = i + 1; j < szA; j++)
          {
            D(idx++, 0) = LevensDist (s1, strA(j).string_value());
          }
        }
        retval(0) = D;
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
      // Compute the distance vector
      if (szA == 1 && szB != 1)
      {
        if (doMinDist)
        {
          boolMatrix D (szB, 1);
          for (octave_idx_type i = 0; i < szB; i++)
          {
            D(i,0) = minLevensDist (strA(0).string_value(),
                                    strB(i).string_value(), minDist);
          }
          retval(0) = D;
        }
        else
        {
          Matrix D (szB, 1);
          for (octave_idx_type i = 0; i < szB; i++)
          {
            D(i,0) = LevensDist (strA(0).string_value(),
                                 strB(i).string_value());
          }
          retval(0) = D;
        }
      }
      else if (szA != 1 && szB == 1)
      {
        if (doMinDist)
        {
          boolMatrix D (szA, 1);
          for (octave_idx_type i = 0; i < szA; i++)
          {
            D(i,0) = minLevensDist (strA(i).string_value(),
                                    strB(0).string_value(), minDist);
          }
          retval(0) = D;
        }
        else
        {
          Matrix D (szA, 1);
          for (octave_idx_type i = 0; i < szA; i++)
          {
            D(i,0) = LevensDist (strA(i).string_value(),
                                 strB(0).string_value());
          }
          retval(0) = D;
        }
      }
      else
      {
        if (doMinDist)
        {
          boolMatrix D (szA, 1);
          for (octave_idx_type i = 0; i < szA; i++)
          {
            D(i,0) = minLevensDist (strA(i).string_value(),
                                    strB(i).string_value(), minDist);
          }
          retval(0) = D;
        }
        else
        {
          Matrix D (szA, 1);
          for (octave_idx_type i = 0; i < szA; i++)
          {
            D(i,0) = LevensDist (strA(i).string_value(),
                                 strB(i).string_value());
          }
          retval(0) = D;
        }
      }

    }
    else if (args(0).is_string() && args(1).is_string())
    {
      if (doMinDist)
      {
        retval(0) = minLevensDist (args(0).string_value(),
                                   args(1).string_value(), minDist);
      }
      else
      {
        retval(0) = LevensDist (args(0).string_value(), args(1).string_value());
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
%! d = editDistance ({"AS","SD","AD"}, 1);
%! assert (d, logical ([0; 1; 1]));
%! assert (class (d), "logical");
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
%! d = editDistance ({"AS","SD","AD"}, {"AS", "AD", "SE"}, 1);
%! assert (d, logical ([1; 1; 0]));
%! assert (class (d), "logical");
%!test
%! d = editDistance ({"AS","SD","AD"}, {"AS"}, 1);
%! assert (d, logical ([1; 0; 1]));
%! assert (class (d), "logical");
%!test
%! d = editDistance ({"AS"}, {"AS", "AD", "SE"}, 1);
%! assert (d, logical ([1; 1; 0]));
%! assert (class (d), "logical");
%!test
%! d = editDistance ("Octave", "octave");
%! assert (d, 1);
%! assert (class (d), "double");
%!test
%! d = editDistance ("Octave", "octave", 1);
%! assert (d, true);
%! assert (class (d), "logical");
*/

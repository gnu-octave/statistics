/*
Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>

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

#include <cmath>
#include <vector>

// Barnes-Hut summation of the t-SNE repulsive term.
//
// The gradient of the t-SNE cost splits into an attraction over the pairs
// carrying a nonzero P, which is cheap and stays in the interpreter, and a
// repulsion over EVERY pair, which is what makes the exact algorithm quadratic.
// This file computes the repulsion in O(N log N) by the method of van der
// Maaten (2014): build a 2^D-ary tree over the embedding, then walk it once per
// point, collapsing a whole cell into its centre of mass whenever the cell is
// far enough away that its internal structure cannot matter.
//
// Returned unnormalized, because the normalizer is itself one of the sums the
// walk accumulates and the caller needs both:
//
//   FREP(i,:) = sum over cells  n * q^2 * (y_i - com),   q = 1 / (1 + d2)
//   Z         = sum over cells  n * q                    (over every i)
//
// so the repulsive force on point i is FREP(i,:) / Z and the low-dimensional
// affinities are q / Z.

// The subdivision stops here whatever the criterion says.  Coincident points
// never separate, so a tree that splits until every cell holds one point does
// not terminate on duplicated rows, which a t-SNE embedding produces readily
// in its first iterations while every point still sits on the origin.
static const int BHT_MAX_DEPTH = 40;

struct BHNode
{
  // Geometry of the cell, as centre and half-width per dimension.
  std::vector<double> centre;
  std::vector<double> halfw;

  // What the cell carries: how many points, and where their centre of mass is.
  octave_idx_type count;
  std::vector<double> com;

  // The longest edge, which is the length scale the opening criterion uses.
  double width;

  // A leaf names its points; an internal node names its children.  Exactly one
  // of the two is ever non-empty.
  std::vector<octave_idx_type> points;
  std::vector<octave_idx_type> children;
};

// Build the subtree rooted at NODE over the point indices IDX, subdividing
// until a cell holds a single point or the depth cap stops it.
static void
bht_build (std::vector<BHNode>& tree, octave_idx_type node,
           const Matrix& Y, std::vector<octave_idx_type>& idx, int depth)
{
  const octave_idx_type d = Y.columns ();
  const octave_idx_type n = idx.size ();

  tree[node].count = n;
  tree[node].com.assign (d, 0.0);
  for (octave_idx_type k = 0; k < n; k++)
  {
    for (octave_idx_type j = 0; j < d; j++)
    {
      tree[node].com[j] += Y(idx[k], j);
    }
  }
  for (octave_idx_type j = 0; j < d; j++)
  {
    tree[node].com[j] /= (double) n;
  }

  tree[node].width = 0.0;
  for (octave_idx_type j = 0; j < d; j++)
  {
    double w = 2.0 * tree[node].halfw[j];
    if (w > tree[node].width)
    {
      tree[node].width = w;
    }
  }

  if (n <= 1 || depth >= BHT_MAX_DEPTH)
  {
    tree[node].points = idx;
    return;
  }

  // Sort the points into the 2^d octants of this cell, one bit per dimension.
  const octave_idx_type nsub = (octave_idx_type) 1 << d;
  std::vector<std::vector<octave_idx_type>> bucket (nsub);
  for (octave_idx_type k = 0; k < n; k++)
  {
    octave_idx_type b = 0;
    for (octave_idx_type j = 0; j < d; j++)
    {
      if (Y(idx[k], j) > tree[node].centre[j])
      {
        b |= ((octave_idx_type) 1 << j);
      }
    }
    bucket[b].push_back (idx[k]);
  }

  // Every point landing in one octant means the cell cannot be split usefully
  // at this level, but the halved cell still tightens, so recursion continues
  // and the depth cap is what ends it.
  for (octave_idx_type b = 0; b < nsub; b++)
  {
    if (bucket[b].empty ())
    {
      continue;
    }

    BHNode child;
    child.centre.resize (d);
    child.halfw.resize (d);
    for (octave_idx_type j = 0; j < d; j++)
    {
      double h = 0.5 * tree[node].halfw[j];
      child.halfw[j] = h;
      child.centre[j] = (b & ((octave_idx_type) 1 << j))
                        ? tree[node].centre[j] + h
                        : tree[node].centre[j] - h;
    }
    child.count = 0;

    tree.push_back (child);
    octave_idx_type c = tree.size () - 1;
    tree[node].children.push_back (c);
    bht_build (tree, c, Y, bucket[b], depth + 1);
  }
}

// Accumulate the repulsion felt by point I from the subtree at NODE.  A leaf is
// summed pair by pair, skipping I itself; an internal node is collapsed into
// its centre of mass when the opening criterion allows and descended otherwise.
static void
bht_forces (const std::vector<BHNode>& tree, octave_idx_type node,
            const Matrix& Y, octave_idx_type i, double theta,
            double *frep, double& Z)
{
  const octave_idx_type d = Y.columns ();
  const BHNode& nd = tree[node];

  if (nd.count == 0)
  {
    return;
  }

  if (! nd.points.empty ())
  {
    for (std::size_t k = 0; k < nd.points.size (); k++)
    {
      octave_idx_type j = nd.points[k];
      if (j == i)
      {
        continue;
      }
      double d2 = 0.0;
      for (octave_idx_type m = 0; m < d; m++)
      {
        double t = Y(i, m) - Y(j, m);
        d2 += t * t;
      }
      double q = 1.0 / (1.0 + d2);
      Z += q;
      double qq = q * q;
      for (octave_idx_type m = 0; m < d; m++)
      {
        frep[m] += qq * (Y(i, m) - Y(j, m));
      }
    }
    return;
  }

  // Distance from the point to the cell's centre of mass.
  double d2 = 0.0;
  for (octave_idx_type m = 0; m < d; m++)
  {
    double t = Y(i, m) - nd.com[m];
    d2 += t * t;
  }

  // The opening criterion.  A cell counts as distant when its width is small
  // against that distance; theta = 0 therefore opens every cell and the walk
  // degenerates to the exact pairwise sum, which is what pins this file's
  // correctness without an oracle.
  bool distant = (nd.width < theta * std::sqrt (d2));

  if (distant)
  {
    double q = 1.0 / (1.0 + d2);
    double n = (double) nd.count;
    Z += n * q;
    double nqq = n * q * q;
    for (octave_idx_type m = 0; m < d; m++)
    {
      frep[m] += nqq * (Y(i, m) - nd.com[m]);
    }
    return;
  }

  for (std::size_t k = 0; k < nd.children.size (); k++)
  {
    bht_forces (tree, nd.children[k], Y, i, theta, frep, Z);
  }
}

DEFUN_DLD(__bhtsne__, args, ,
"-*- texinfo -*-\n\
@deftypefn {statistics} {[@var{Frep}, @var{Z}] =} __bhtsne__ (@var{Y}, @var{theta})\n\
\n\
Barnes-Hut summation of the t-SNE repulsive term.  Internal; called by\n\
@code{tsne} and not meant to be used directly.\n\
\n\
@var{Y} is the @math{N*D} embedding, @math{D} being 1, 2 or 3, and @var{theta}\n\
the opening criterion, a non-negative scalar.  @var{Frep} is the @math{N*D}\n\
matrix of unnormalized repulsive forces and @var{Z} the normalizer, so the\n\
force on a point is its row of @var{Frep} divided by @var{Z}.\n\
\n\
@var{theta} of zero opens every cell and reproduces the exact pairwise sum.\n\
\n\
@end deftypefn")
{
  if (args.length () != 2)
  {
    print_usage ();
  }

  if (! args(0).isnumeric () || args(0).iscomplex () || args(0).isempty ())
  {
    error ("__bhtsne__: Y must be a real numeric matrix.");
  }
  if (! args(1).is_scalar_type () || ! args(1).isnumeric ()
      || args(1).iscomplex ())
  {
    error ("__bhtsne__: THETA must be a real scalar.");
  }

  Matrix Y = args(0).matrix_value ();
  double theta = args(1).scalar_value ();

  if (theta < 0.0 || ! octave::math::isfinite (theta))
  {
    error ("__bhtsne__: THETA must be non-negative and finite.");
  }

  const octave_idx_type n = Y.rows ();
  const octave_idx_type d = Y.columns ();

  if (d < 1 || d > 3)
  {
    error ("__bhtsne__: Y must have 1, 2, or 3 columns.");
  }

  // The root cell, sized to contain every point with room to spare so that a
  // point never sits exactly on a boundary.
  BHNode root;
  root.centre.resize (d);
  root.halfw.resize (d);
  for (octave_idx_type j = 0; j < d; j++)
  {
    double lo = Y(0, j);
    double hi = Y(0, j);
    for (octave_idx_type i = 1; i < n; i++)
    {
      if (Y(i, j) < lo)
      {
        lo = Y(i, j);
      }
      if (Y(i, j) > hi)
      {
        hi = Y(i, j);
      }
    }
    root.centre[j] = 0.5 * (lo + hi);
    root.halfw[j] = 0.5 * (hi - lo) + 1e-5;
  }
  root.count = 0;

  std::vector<BHNode> tree;
  tree.push_back (root);

  std::vector<octave_idx_type> idx (n);
  for (octave_idx_type i = 0; i < n; i++)
  {
    idx[i] = i;
  }
  bht_build (tree, 0, Y, idx, 0);

  Matrix Frep (n, d, 0.0);
  double Z = 0.0;
  std::vector<double> frep (d);
  for (octave_idx_type i = 0; i < n; i++)
  {
    for (octave_idx_type m = 0; m < d; m++)
    {
      frep[m] = 0.0;
    }
    bht_forces (tree, 0, Y, i, theta, &frep[0], Z);
    for (octave_idx_type m = 0; m < d; m++)
    {
      Frep(i, m) = frep[m];
    }
  }

  return ovl (Frep, Z);
}

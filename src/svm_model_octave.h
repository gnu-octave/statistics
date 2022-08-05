/*
Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
Based on the Octave LIBSVM wrapper created by Alan Meeson (2014) based on an
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
#include <octave/ov-struct.h>
#include "svm.h"

const char *model_to_octave_structure(octave_value_list &plhs, int num_of_feature, struct svm_model *model);
struct svm_model *octave_matrix_to_model(octave_scalar_map &octave_struct, const char **error_message);
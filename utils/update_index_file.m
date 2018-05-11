## Copyright (C) 2018 John Donoghue
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {} update_index_file (@var{indexfilename}, @var{altindexfilename}, @var{altfunctionlist})
## insert function from altfunctionlist into INDEX file indexfilename, using index info 
## from INDEX file altindexfilename.
##
## ie: a conditional merge of two INDEX files
##
## @end deftypefn

function update_index_file (indexfilename, infoindexname, installed_functions)

  index_info = read_index_file (indexfilename);
  tmp_info = read_index_file (infoindexname);

  for i = 1:numel(installed_functions)
    func = get_index_file_function (tmp_info, installed_functions{i});
    if !isempty (func)
      index_info = add_index_file_function (index_info, func.category, func.name);
    else
      printf ("warning: didnt find func %s in index\n", installed_functions{i});
    endif
  endfor

  write_index_file (indexfilename, index_info);
endfunction

function index_info = read_index_file (filename)
  index_info = [];
  index_info.name = "<not set>";
  index_info.categories = {};

  categoryname = '';
  functions = {};

  fd = fopen (filename, "rt");

  if fd != -1
    index_info.name = fgetl (fd);
  
    while !feof (fd),
      tmp = fgetl (fd);
      if ischar (tmp)
        if length (tmp) > 0 && tmp (1) != " "
          % start of category
          if numel (functions) > 0
            c = {};
            c.name = categoryname;
            c.functions = functions;
            index_info.categories{end+1} = c;
            functions = {};
          endif
        
          categoryname = tmp;
        
        elseif length (strtrim (tmp)) > 0
          % starts with ' ', add to current category
          tmp_functions = strsplit (strtrim (tmp), {' ', '\t'}, "collapsedelimiters", true);
          for i=1:numel (tmp_functions)
            functions{end+1} = tmp_functions{i};
          endfor
        endif
      endif
    endwhile
  
    % last block of functions not yet added
    if numel (functions) > 0
      c = {};
      c.name = categoryname;
      c.functions = functions;
      index_info.categories{end+1} = c;
    endif
  endif
  
  fclose(fd);
endfunction

function func_info = get_index_file_function(index_info, functioname)
  % search everything until find the matching namelengthmax
  func_info = [];
  found = false;
  
  for c = 1:numel (index_info.categories)
      category = index_info.categories{c};
      
      for f = 1:numel (category.functions)
        func = category.functions{f};
        if strcmp (func, functioname)
          found = true;      
          func_info.name = func;
          func_info.category = category.name;
        endif
      endfor
      
      if found
        break;
      endif
    endfor
endfunction

function [index_info, found] = add_index_file_function (index_info, functioncategory, functionname)
  found = false;
  
  for c = 1:numel (index_info.categories)
    category = index_info.categories{c};
    if strcmpi (category.name, functioncategory)
      found = true;
      category.functions{end+1} = functionname;
      index_info.categories{c} = category;
    endif
  endfor  
  
  if !found
    % category not found, so add one
    c = {};
    c.name = functioncategory;
    c.functions = {};
    c.functions{end+1} = functionname;
    index_info.categories{end+1} = c;
  endif
endfunction

function ok = write_index_file (filename, index_info)
  ok = true;
  
  fd = fopen (filename, "w+t");
  
  if fd != -1
    
    fprintf (fd, "%s\n", index_info.name);
    
    for c = 1:numel (index_info.categories)
      category = index_info.categories{c};
      fprintf (fd, "%s\n", category.name);
      
      for f = 1:numel (category.functions)
        func = category.functions{f};
        fprintf (fd, " %s\n", func);
      endfor
    endfor
    
    fclose (fd);
    
  endif
  
endfunction

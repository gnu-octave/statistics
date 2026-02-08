classdef myanova

  properties
    % ---------------- Input ----------------
    Y               % response vector
    Groups          % cell array of grouping variables
    ModelSpec       % struct describing the model

    % ---------------- Results ----------------
    SS              % sums of squares
    DF              % degrees of freedom
    MS              % mean squares
    F               % F statistics
    P               % p-values

    Table           % myanova
    
    Stats           % stats struct
  endproperties


  %========================================================
  % PUBLIC METHODS
  %========================================================
  methods

    %---------------- Constructor ----------------
    function obj = myanova(Y, Groups, ModelSpec)

      if (nargin < 3)
        error("myanova : too few input arguments.");
      endif

      if (! isvector (Y))
        error ("myanova : Y must be a vector.");
      endif

      if (! iscell (Groups))
        error ("myanova: Groups must be a cell array.");
      endif

      obj.Y = Y(:);
      obj.Groups = Groups;
      obj.ModelSpec = ModelSpec;

      obj = obj.compute ();
    endfunction


    %---------------- Main pipeline ----------------
    function obj = compute (obj)
      obj = obj.computeSS ();
      obj = obj.computeDF ();
      obj = obj.computeMS ();
      obj = obj.computeF ();
      obj = obj.computeP ();
      obj = obj.buildTable ();
      obj = obj.buildStats ();
    endfunction

  endmethods


  %========================================================
  % PRIVATE METHODS
  %========================================================
  methods (Access = private)

    %---------------- Sum of Squares ----------------
    function obj = computeSS (obj)

      Y = obj.Y;
      G = obj.Groups;
      N = numel (Y);

      grand_mean = mean (Y);
      SST = sum ((Y - grand_mean).^2);

      switch (obj.ModelSpec.Design)

        %================ ONE-WAY =================
        case 'oneway'

          g = G{1}(:);
          levels = unique (g);
          SSA = 0;

          for i = 1:numel (levels)
            yi = Y(g == levels(i));
            SSA += numel (yi) * (mean (yi) - grand_mean)^2;
          endfor

          SSE = SST - SSA;
          obj.SS = [SSA; SSE; SST];


        %================ TWO-WAY =================
        case 'twoway'

          g1 = G{1}(:);
          g2 = G{2}(:);

          lv1 = unique (g1);
          lv2 = unique (g2);

          SSA = 0;
          SSB = 0;
          SSAB = 0;

          % Main effect A
          for i = 1:numel (lv1)
            yi = Y(g1 == lv1(i));
            SSA += numel (yi) * (mean (yi) - grand_mean)^2;
          endfor

          % Main effect B
          for j = 1:numel (lv2)
            yj = Y(g2 == lv2(j));
            SSB += numel (yj) * (mean (yj) - grand_mean)^2;
          endfor

          % Interaction A×B
          if (isfield (obj.ModelSpec, "Interactions") && ...
              obj.ModelSpec.Interactions)

            for i = 1:numel (lv1)
              for j = 1:numel (lv2)
                idx = (g1 == lv1(i)) & (g2 == lv2(j));
                yij = Y(idx);
                if (! isempty (yij))
                  SSAB += numel (yij) * ...
                    ( mean (yij) ...
                    - mean (Y(g1 == lv1(i))) ...
                    - mean (Y(g2 == lv2(j))) ...
                    + grand_mean )^2;
                endif
              endfor
            endfor

          endif

          SSE = SST - SSA - SSB - SSAB;

          if (isfield (obj.ModelSpec, "Interactions") && ...
              obj.ModelSpec.Interactions)
            obj.SS = [SSA; SSB; SSAB; SSE; SST];
          else
            obj.SS = [SSA; SSB; SSE; SST];
          endif


        otherwise
          error ("myanova: unknown design type.");
      endswitch

    endfunction


    %---------------- Degrees of Freedom ----------------
    function obj = computeDF (obj)

      N = numel (obj.Y);
      G = obj.Groups;

      switch (obj.ModelSpec.Design)

        case 'oneway'
          a = numel (unique (G{1}));
          obj.DF = [a-1; N-a; N-1];

        case 'twoway'
          a = numel (unique (G{1}));
          b = numel (unique (G{2}));

          dfA = a - 1;
          dfB = b - 1;
          dfAB = (a - 1) * (b - 1);

          if (isfield (obj.ModelSpec, "Interactions") && ...
              obj.ModelSpec.Interactions)
            dfE = N - a*b;
            obj.DF = [dfA; dfB; dfAB; dfE; N-1];
          else
            dfE = N - a - b + 1;
            obj.DF = [dfA; dfB; dfE; N-1];
          endif

      endswitch

    endfunction


    %---------------- Mean Squares ----------------
    function obj = computeMS (obj)
      obj.MS = obj.SS ./ obj.DF;
    endfunction


    %---------------- F statistics ----------------
    function obj = computeF (obj)
      MSerr = obj.MS(end-1);
      obj.F = obj.MS(1:end-2) ./ MSerr;
    endfunction


    %---------------- P values ----------------
    function obj = computeP (obj)
      df1 = obj.DF(1:end-2);
      df2 = obj.DF(end-1);
      obj.P = 1 - fcdf (obj.F, df1, df2);
    endfunction


    %---------------- Build myanovatable ----------------
  function obj = buildTable (obj)

  switch (obj.ModelSpec.Design)

    case 'oneway'
      src = {'Groups'; 'Error'; 'Total'};

    case 'twoway'
      if (isfield (obj.ModelSpec, "Interactions") && ...
          obj.ModelSpec.Interactions)
        src = {'Factor1'; 'Factor2'; 'Factor1:Factor2'; 'Error'; 'Total'};
      else
        src = {'Factor1'; 'Factor2'; 'Error'; 'Total'};
      endif

  endswitch

  % Preallocate table (rows = sources + header)
  obj.Table = cell (numel (src) + 1, 6);

  % Header
  obj.Table(1,:) = {'Source', 'SS', 'df', 'MS', 'F', 'Prob>F'};

  % Data rows
  for i = 1:numel (src)

    obj.Table{i+1,1} = src{i};        % Source name
    obj.Table{i+1,2} = obj.SS(i);     % SS
    obj.Table{i+1,3} = obj.DF(i);     % df
    obj.Table{i+1,4} = obj.MS(i);     % MS

    if (i <= numel (obj.F))
      obj.Table{i+1,5} = obj.F(i);    % F
      obj.Table{i+1,6} = obj.P(i);    % p-value
    else
      obj.Table{i+1,5} = NaN;
      obj.Table{i+1,6} = NaN;
    endif

  endfor

endfunction




    %---------------- Stats struct ----------------
    function obj = buildStats (obj)
      obj.Stats = struct ();
      obj.Stats.df = obj.DF(end-1);
      obj.Stats.s = sqrt (obj.MS(end-1));
      obj.Stats.gnames = obj.Groups;
      obj.Stats.modelspec = obj.ModelSpec;
    endfunction

  endmethods
endclassdef

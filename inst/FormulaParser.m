function varargout = FormulaParser(varargin)

    if (nargin < 1)
        error("FormulaParser : Input formula string is required.");
    endif

    formula_str = varargin{1};
    
    if (nargin > 1)
        mode = varargin{2};
    else
        mode = "tokenize";
    endif

    switch (lower(mode))
        case "tokenize"
            varargout{1} = run_lexer(formula_str);
            
        case "parse"
            
            
        otherwise
            error("FormulaParser : Unknown mode: %s", mode);
    endswitch
endfunction

function tokens = run_lexer(formula_str)
    if (isempty(formula_str))
        tokens = struct("type", {}, "value", {}, "pos", {});
        return;
    endif
    
    str = char(formula_str);
    n = length(str);
    
    tokens(n) = struct("type", "", "value", "", "pos", 0);
    tok_idx = 0;
    
    i = 1;
    while (i <= n)
        c = str(i);
        start_pos = i;
        
        if (isspace(c))
            i = i + 1;
            continue;
        endif
        
        ## Factors and Variables
        if (isletter(c))
            val = c;
            i = i + 1;
            while (i <= n)
                next_c = str(i);
                if (isletter(next_c) || (next_c >= '0' && next_c <= '9') || next_c == '_')
                    val = [val, next_c];
                    i = i + 1;
                else
                    break;
                endif
            endwhile
            tok_idx = tok_idx + 1;
            tokens(tok_idx) = create_token("IDENTIFIER", val, start_pos);
            continue;
        endif
        
        ## Numbers
        if (c >= '0' && c <= '9')
            val = c;
            i = i + 1;
            while (i <= n)
                next_c = str(i);
                if (next_c >= '0' && next_c <= '9')
                    val = [val, next_c];
                    i = i + 1;
                elseif (next_c == '.')
                    ## Distinguish decimal from interaction operator.
                    if (i + 1 <= n && str(i+1) >= '0' && str(i+1) <= '9')
                        val = [val, next_c];
                        i = i + 1;
                    else
                        break; 
                    endif
                else
                    break;
                endif
            endwhile
            tok_idx = tok_idx + 1;
            tokens(tok_idx) = create_token("NUMBER", val, start_pos);
            continue;
        endif
        
        ## Operators
        type = "";
        val = c;
        skip = 0;
        
        switch (c)
            case '-'
                ## Check for Deletion operators.
                if (i + 1 <= n)
                    if (str(i+1) == '/')
                        type = "OP_MINUS_MARGIN"; 
                        val = "-/";
                        skip = 1;
                    elseif (str(i+1) == '*')
                        type = "OP_MINUS_CLEAN"; 
                        val = "-*";
                        skip = 1;
                    else
                        type = "OP_MINUS"; 
                    endif
                else
                    type = "OP_MINUS";
                endif
                
            case '*'
                ## Check for Exponentiation.
                if (i + 1 <= n && str(i+1) == '*')
                    type = "OP_POWER";
                    val = "**";
                    skip = 1;
                else
                    type = "OP_CROSS";
                endif
                
            case '/'
                type = "OP_NEST";
                
            case '+'
                type = "OP_PLUS";
                
            case {'.', ':'}
                type = "OP_DOT";
                
            case '^'
                type = "OP_POWER";
                
            case '~'
                type = "SEPARATOR";
                
            case '('
                type = "LPAREN";
                
            case ')'
                type = "RPAREN";
                
            case ','
                type = "COMMA";
                
            otherwise
                error("FormulaParser : Unexpected character '%s' at position %d", c, i);
        endswitch
        
        tok_idx = tok_idx + 1;
        tokens(tok_idx) = create_token(type, val, start_pos);
        i = i + 1 + skip;
    endwhile
    
    tokens = tokens(1:tok_idx);
    tokens(end+1) = create_token("EOF", "EOF", i);
endfunction

function t = create_token(type, val, pos)
    t.type = type;
    t.value = val;
    t.pos = pos;
endfunction
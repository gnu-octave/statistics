## fullfact([n1 n2 n3])
##    full factorial design with choices 1 through n_i for each factor i
##
## fullfact(n)
##    full factorial design with n binary choices, 0 and 1
function A = fullfact(n)
  if length(n) == 1
    % combinatorial design with n either/or choices
    A = fullfact(2*ones(1,n))-1;
  else
    % combinatorial design with n(i) choices per level
    A = [1:n(end)]';
    for i=length(n)-1:-1:1
      A = [kron([1:n(i)]',ones(rows(A),1)), repmat(A,n(i),1)];
    end
  end

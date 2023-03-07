function med = find_median(x)
    n = length(x);
    x = sort(x);
    if mod(n, 2) == 1
        med = x((n+1)/2);
    else
        med = (x(n/2) + x(n/2+1)) / 2;
    end
end

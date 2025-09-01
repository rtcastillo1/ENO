%% ENO coefficients with uniform stencil
% This program generates uniform coefficients for ENO approximation, r is
% the offset, k is the number of points and k+1 is the order. Note that r
% can accept both positive and negative values (to the right of the
% reference index i).

function output = eno_unicoeff(r0,k0)
r = r0;
k = k0;

mvals = zeros(1, k+1);
for j = 0:k
    total_sum = 0;
    for m = j+1:k
        % Numerator
        num_sum = 0;
        for l = 0:k
            if l == m, continue; end
            prod_q = 1;
            for q = 0:k
                if q == l || q == m, continue; end
                prod_q = prod_q * (r - q + 1);
            end
            num_sum = num_sum + prod_q;
        end
        
        % Denominator
        prod_l = 1;
        for l = 0:k
            if l == m, continue; end
            prod_l = prod_l * (m - l);
        end
        
        total_sum = total_sum + num_sum/prod_l;
    end
    mvals(j+1) = total_sum;
end

output = mvals;
end




%% UNDIVIDED DIFFERENCE. 
%{ 
This program computes the 0th to nth order undivided difference, where n is
the order, x is the domain, y is function values over x. Actually, x0 is
not neeed. But we include it anyway.

The first column corresponds to the 0th order undivided difference.
%}

function output = eno_undivdiff(y0)

n = length(y0); % Number of points
% x = x0;        % Not really needed
y = y0;         % Function values

dd = zeros(n,n);
dd(:,1) = y; % Initialize

for j = 2:n
    for i = 1:n-(j-1)
        dd(i,j) = (dd(i+1,j-1)-dd(i,j-1));
    end
end

output = dd;

end

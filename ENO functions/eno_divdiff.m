%% DIVIDED DIFFERENCE. 
%{ 
This program computes the 0th to nth order divided difference, where n is
the order, x is the domain, y is function values over x. 

The first column corresponds to the 0th order undivided difference.
%}

function output = eno_divdiff(x0,y0)

n = length(x0); % Number of points
x = x0;         % x values     
y = y0;         % Function values

dd = zeros(n,n);
dd(:,1) = y; % Initialize

for j = 2:n
    for i = 1:n-(j-1)
        dd(i,j) = (dd(i+1,j-1)-dd(i,j-1))/(x(i+j-1)-x(i));
    end
end

output = dd;

end

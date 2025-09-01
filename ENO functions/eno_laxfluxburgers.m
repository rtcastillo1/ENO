%% LAX-FRIEDRICH FLUX - BURGERS' EQUATION
% This function implements the Lax-Friedrich Flux for ENO algorithm. 
% NOTE: THIS FLUX FUNCTION IS ONLY FOR THE BURGERS' EQUATION. 
% In this function, a0 is the lower cell boundary, b0 is the higher cell
% boundary, u0 is the 

function output = eno_laxfluxburgers(a0,b0,u0)

a = a0;
b = b0;
u = u0;

% We introduce flux. This only applies for Burgers' equation.
f =@(u) (1/2)*u.^2;
alpha = max(u);

h =@(a,b) (1/2)*(f(a) + f(b) - alpha.*(b-a));

output = h(a,b);

end
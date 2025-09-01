%% RK4 SOLVER
%{
This is a generic RK4 solver for ODEs or systems of ODEs. 

rk4_solve(x,init)
where x is the differential equation
init is the initial value of x
%}

function ret = rk4_solve(x,init)

%tm=1:.1:10; %time horizon

y(1) = init;

N = 10;
h = 1/N;
t = linspace(0,1,N+1);


for i = 1:N;
    k(1) = x(i,y(i));
    k(2) = x(i + h/2, x(i,y(i))+(h/2)*k(1));
    k(3) = x(i + h/2, x(t,y(i))+(h/2)*k(2));
    k(4) = x(i + h, x(i,y(i))+h*k(3));    
    y(i+1) = y(i) + h/6*(k(1) + 2*k(2) + 2*k(3) + k(4));
end

ret = y;

end
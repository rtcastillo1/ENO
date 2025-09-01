%% Solution of Inviscid Burger's equation using Upwind vs Lax-Wendroff method

clear
close all hidden
clc


%% SETUP 

xpoints = 100;
tpoints = 200;
tee = linspace(0, 1, tpoints);
ecks = linspace(0, 1, xpoints);
dx = (ecks(2)-ecks(1));
dt = (tee(2)-tee(1));
mu = dt/dx;

%% UPWIND SCHEME


% Initialize solution matrix
uu = zeros(tpoints,xpoints);

% Initial values

uu(1, :) = zeros(1, xpoints);  % Start with all zeros

% Define ramp region
x_start = 0.2;
x_end = 0.8;
ramp_idx = (ecks >= x_start) & (ecks <= x_end);
uu(1, ecks < x_start) = 1;  % flat at 1 before the ramp
uu(1, ramp_idx) = 1 - (ecks(ramp_idx) - x_start) / (x_end - x_start);  % decreasing ramp
% Already 0 for ecks > x_end due to initialization

% Boundary condition 
uu(:,1) = 1;
uu(:, end) = uu(:, end-1);

% SOLUTION 
for n = 1:tpoints-1
    for j = 2:xpoints-1
        aa = uu(n, j); 
        if aa<0
            uu(n+1,j) = uu(n,j) - aa*mu*(uu(n,j+1) - uu(n,j));
        else
            uu(n+1,j) = uu(n,j) - aa*mu*(uu(n,j) - uu(n,j-1));
        end
    end
    uu(n+1, end) = uu(n+1, end-1);
end

%% LAX WENDROFF 

% Initialize solution matrix
uul = zeros(tpoints,xpoints);

% Initial values
uul(1, :) = zeros(1, xpoints);  % Start with all zeros

% Define ramp region
x_start = 0.2;
x_end = 0.8;
ramp_idx = (ecks >= x_start) & (ecks <= x_end);
uul(1, ecks < x_start) = 1;  % flat at 1 before the ramp
uul(1, ramp_idx) = 1 - (ecks(ramp_idx) - x_start) / (x_end - x_start);  % decreasing ramp
% Already 0 for ecks > x_end due to initialization

% Boundary condition 
uul(:,1) = 1;
uul(:, end) = uul(:, end-1);

F =@(zz) ((1/2)*zz^2); % See the conservative form

for n = 1:tpoints-1
    for j = 2:xpoints-1
        a1 = (uul(n,j) + uul(n,j-1))/2;
        a2 = (uul(n,j) + uul(n,j+1))/2;

        uul(n+1,j) = uul(n,j) - (1/2)*mu*( (1-a2*mu)*( F(uul(n,j+1)) - F(uul(n,j)) ) + ...
            (1 + a1*mu)* ( F(uul(n,j)) - F(uul(n,j-1)) ) );
    end
    uul(n+1, end) = uul(n+1, end-1);
end




%% PLOT
gif_filename = 'burgers_solution.gif';
figure
for n = 1:tpoints
    plot(ecks, uu(n,:), 'b-', ecks, uul(n,:), 'r--', 'LineWidth', 2)
    axis([0 1 0 1.2])
    xlabel('x')
    ylabel('u(x,t)')
    legend('Upwind','Lax-Wendroff')
    title(sprintf('Numerical solution of Inviscid Burgers'' Equation at t = %.5f', tee(n)))
    grid on

    drawnow
    frame = getframe(gcf);
    img = frame2im(frame);
    [A,map] = rgb2ind(img,256);
    
    if n == 1
        imwrite(A,map,gif_filename,'gif','LoopCount',Inf,'DelayTime',0.05);
    else
        imwrite(A,map,gif_filename,'gif','WriteMode','append','DelayTime',0.05);
    end
end
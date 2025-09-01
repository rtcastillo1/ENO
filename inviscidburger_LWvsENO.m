%% Solution of Inviscid Burger's equation: Comparison of Lax-Wendroff method and ENO

clear
close all hidden
clc

%% LAX-WENDROFF ===========================================================

xpoints = 100;
tpoints = 200;
tee = linspace(0, 1, tpoints);
ecks = linspace(0, 1, xpoints);
dx = (ecks(2)-ecks(1));
dt = (tee(2)-tee(1));
mu = dt/dx;

% Initialize solution matrix
uu = zeros(tpoints,xpoints);

% Initial values with discontinuity
uu(1,:) = 1;
uu(1,26:end) = .5;

% Boundary condition 
uu(:,1) = 1;
uu(:, end) = uu(:, end-1);

% SOLUTION 
F =@(zz) ((1/2)*zz^2); % See the conservative form
for n = 1:tpoints-1
    for j = 2:xpoints-1
        a1 = (uu(n,j) + uu(n,j-1))/2;
        a2 = (uu(n,j) + uu(n,j+1))/2;

        uu(n+1,j) = uu(n,j) - (1/2)*mu*( (1-a2*mu)*( F(uu(n,j+1)) - F(uu(n,j)) ) + ...
            (1 + a1*mu)* ( F(uu(n,j)) - F(uu(n,j-1)) ) );
    end
    uu(n+1, end) = uu(n+1, end-1);
end



% PLOT
gif_filename = 'burgers_solution.gif';
figure
for n = 1:tpoints
    plot(ecks, uu(n,:), 'b.-', 'LineWidth', 2)
    axis([0 1 0 1.2])
    xlabel('x')
    ylabel('u(x,t)')
    title(sprintf('Inviscid Burgers'' Equation at t = %.5f', tee(n)))
    grid on
    
    drawnow
    frame = getframe(gcf);
    img = frame2im(frame);
    [A,map] = rgb2ind(img,256);
    
    if n == 1
        imwrite(A,map,gif_filename,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,gif_filename,'gif','WriteMode','append','DelayTime',0.1);
    end
end

%% ENO ====================================================================
% Add folder of ENO functions to the path
addpath('C:\Users\czjca\Documents\000 - Math\Math --- - Numerical PDE\Numerics\ENO\ENO\ENO functions')

xpoints = 20;
tpoints = 500;
ecks = linspace(0,1,xpoints);
tee = linspace(0,1,tpoints);
mu = (tee(2)-tee(1))/(ecks(2)-ecks(1));

% Initialize solution matrix
uu1 = zeros(tpoints,xpoints);
% Initial values with discontinuity
discont = floor(xpoints/2);
uu1(1,:) = 1;
uu1(1,discont:end) = .5;

%uu1(1,:) = sin(10*ecks);
%uu1(1,discont:end) = .5;

% Boundary condition 
uu1(:,1) = 1;
uu1(:, end) = uu1(:, end-1);

% Step 1: Compute undivided difference because the stencil is uniform.

% Step 2: Create a stencil. The initial stencil is just the single average.
% In step 3, we will expand this stencil by adding k points to the left or
% to the right. 
k = 3; % number of points, also order of approximation
i = k+1; % trial cell. To be looped later. The loop will go from k+1 to n-k.

sol = zeros(tpoints,xpoints-2*(k+1));

for t = 1:tpoints;
y0 = uu1(t,:);

% Step 3: Compare undivided differences since we deal with uniform stencil

for i = k+1:xpoints-k
    initstencil = y0(i);
    label0 = ["center"];
    trackr = 0; % We track the offset created by adding points.

    for kk = 1:k
        temp0 = initstencil;
        label1 = label0;
        tempbck = [y0(i-kk) temp0]; % Append to the left
        tempfwd = [temp0 y0(i+kk)]; % Append to the right
        n = length(temp0)+1; % This is also the order of differencing
        
        % We take the sum because this should just be a scalar because the
        % the order of differencing equals the length of the stencil.
        diffbck = eno_undivdiff(tempbck);
        diffbck = diffbck(1,n); % The first entry of the last column is the desired difference.
        difffwd = eno_undivdiff(tempfwd);
        difffwd = difffwd(1,n);
    
        % We compare the differences and decide where to add points.
        
        if abs(diffbck) < abs(difffwd)
            temp1 = tempbck;
            label1 = ['left' label1];
            trackr = trackr+1;
        else
            temp1 = tempfwd;
            label1 = [label1 'right'];
        end
        initstencil = temp1;
        label0 = label1;
    end
    r = find(label0=="center")-1; % This is the offset from the initial i (center).
    if trackr == r
        %display("r is correct")
    else
        %display("r is wrong")
    end
    % Step 4: We interpolate using these points to get the left and right
    % points to y_i.
    
    % Obtain coefficients
    coef = eno_unicoeff(r,k+1);
    coef = coef(1:k+1);
    
    v_i(i-k) = coef*initstencil'; % This is v_i+1/2. We loop over i to get all these.
end

% Setup fhatplus and fhatminus
uminus = v_i(1:end-1);
uplus = v_i(2:end);

fhat = eno_laxfluxburgers(uminus,uplus,y0);
fhatplus = fhat(2:end);
fhatminus = fhat(1:end-1);

% Setup the conservative scheme
uu1(t+1,k+2:xpoints-(k+1)) =  uu1(t,k+2:xpoints-(k+1)) - mu*(fhatplus - fhatminus);
end

final = uu1(:,k+2:xpoints-(k+1));

% PLOT 
figure;
for t = 1:tpoints
    plot(0:11,final(t,:),'LineWidth',1.5);
    xlabel('x'); ylabel('u');
    title(['Burgers Equation, t = ', num2str(tee(t))]);
    ylim([-0.1 1.1]); % adjust to match your solution range
    grid on;
    drawnow;
end

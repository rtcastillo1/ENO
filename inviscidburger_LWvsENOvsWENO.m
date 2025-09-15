%% Solution of Inviscid Burger's equation: Comparison of Lax-Wendroff method and ENO and WENO

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

xpoints = 100;
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

% We just fill in these parts with the boundary condition because these are 
% left out in the algorithm. 
uu1(t+1,1:k+1) = 1; 
uu1(t+1,end-k:end) = 0.5;
end

final = uu1(:,k+2:xpoints-(k+1));

% PLOT 
figure;
for t = 1:tpoints
    plot(0:xpoints-(2*(k+1)+1),final(t,:),'LineWidth',1.5);
    xlabel('x'); ylabel('u');
    title(['Burgers Equation, t = ', num2str(tee(t))]);
    ylim([-0.1 1.1]); % adjust to match your solution range
    grid on;
    drawnow;
end

%% WENO ===================================================================
% Add folder of ENO functions to the path. We need the reconstruction
% algorithm that we used in ENO. 
addpath('C:\Users\czjca\Documents\000 - Math\Math --- - Numerical PDE\Numerics\ENO\ENO\ENO functions')

xpoints = 50;
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

% Step 1: We obtain the reconstruction.
k = 3; % number of points, also order of approximation
i = k+1; % trial cell. To be looped.
v0 = uu1(1,:); % To be updated later.

for t=1:tpoints
v0 = uu1(t,:);

% We loop everything together.
stencilset = zeros(k,k,xpoints-2*(k-1)); % This is for checking of stencil.
vstencil = zeros(k,k,xpoints-2*(k-1)); % We extract the function value corresponding to the stencil.
coefstencil = zeros(k,k,xpoints-2*(k-1)); % We obtain the coefficient for reconstruction.
for i = k:xpoints-(k-1)
    for r = 0:k-1
       stencilset(r+1,:,i-(k-1)) = i-r:i-r+k-1;
       vstencil(r+1,:,i-(k-1)) = v0(i-r:i-r+k-1);
       temp = eno_unicoeff(r,k);
       coefstencil(r+1,:,i-(k-1)) = temp(1:k);
    end
end
% Note: The rows is each vstencil correspond to the value of v as we change
% r. There are 3 points in each batch.

% Obtain the reconstructed values for all i, r. 
recon = zeros(k,xpoints-2*(k-1));
for j = 1:k
    for i = 1:xpoints-2*(k-1)
        recon(j,i) = vstencil(j,:,i)*coefstencil(j,:,i)';
    end
end

% Step 2: We obtain the constants dr and dr_bar. For now, we will just use
% pre-computed values. 
dr = flip([1/10 3/5 3/10])';

% Step 3: We compute the weights omega. For now, we will just use
% the formula for k=3 for betas, omegas and alphas.
% Initialize variables. 
alpha = zeros(k,xpoints-2*(k-1));
omega = zeros(k,xpoints-2*(k-1));
beta = zeros(k,xpoints-2*(k-1));

% Functions for beta
betafunc0 =@(vv) (13/12)*(vv(1)-2*vv(2)+vv(3))^2 + ... 
    (1/4)*(3*vv(1) - 4*vv(2) + vv(3))^2;
betafunc1 =@(vv) (13/12)*(vv(1)-2*vv(2)+vv(3))^2 + ...
    (1/4)*(vv(1) - vv(3))^2;
betafunc2 =@(vv) (13/12)*(vv(1)-2*vv(2)+vv(3))^2 + ...
    (1/4)*(vv(1)-4*vv(2)+3*vv(3))^2;
% We compute beta.
for i = 1:xpoints-2*(k-1)
        beta(:,i) = [betafunc0(vstencil(1,:,i));...
                     betafunc1(vstencil(2,:,i));...
                     betafunc2(vstencil(3,:,i))];
end
% Compute alpha and omega. 
epsilon = 10^(-6);
alpha = dr./(epsilon + beta).^2;

alphasum = sum(alpha,1);
omega = alpha./alphasum;

% Step 5: Compute (2k-1)th order reconstruction.
v_i = sum(omega.*recon, 1);


% Setup fhatplus and fhatminus
uminus = v_i(1:end-1);
uplus = v_i(2:end);

fhat = eno_laxfluxburgers(uminus,uplus,v0);
fhatplus = fhat(2:end);
fhatminus = fhat(1:end-1);

% Setup the conservative scheme
uu1(t+1,k+1:xpoints-(k)) =  uu1(t,k+1:xpoints-(k)) - mu*(fhatplus - fhatminus);

% We just fill in these parts with the boundary condition because these are 
% left out in the algorithm. 
uu1(t+1,1:k+1) = 1; 
uu1(t+1,end-k:end) = 0.5;
end

final = uu1(:,k+2:xpoints-(k+1));

% PLOT 
figure;
for t = 1:tpoints
    plot(0:xpoints-(2*(k+1)+1),final(t,:),'LineWidth',1.5);
    xlabel('x'); ylabel('u');
    title('Burgers Equation using WENO');
    ylim([-0.1 1.1]); % adjust to match your solution range
    grid on;
    drawnow;
end

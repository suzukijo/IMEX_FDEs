%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ex3_Case2.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Yongtao Zhou, Jorge Suzuki, Mohsen Zayernouri
% 
% Description: Solution of a nonlinear FDE (Example 3, case 2) from
% "Implicit-explicit time-integration of nonlinar fractional differential
% equations", Applied Numerical Mathematics 156 (2020). The problem is
% defined as:
%
%   C_Du(t) = lambda*u(t) + f(u(t),t),
%
%   lambda = -1, f(u(t),t) = -0.1 u^2 + g(t)
%   g(t) given for a fabricated u(t) = 1 + t + t^2/2 + t^3/3 + t^4/4
%
% Revision history:
%
% Rev1 (September 2019) - Initial file
% Rev2 (November 2021)  - Added selection for first/second-order methods
%

clear all;
clc;

addpath('../src/');
format short e

%% Domain Parameters %%
t0 = 0;              % Initial time
T  = 1;              % Final time
y0 = 1;              % Initial condition  
alpha = 0.2;         % Fractional order value
epsil = 0.5*10^-8;   % Epsilon value for fast Topelitz inversion
tol = 1e-7;          % Numerical tolerance for the Picard iteration (nonlinear cases)

Ntests = 5;          % Number of runs to be performed at different dt
lambdamrho = -1;     % linear term coefficient (lambda <= 0)
imex_type = 2;       % Type of IMEX solver (1 = first order, 2 = second order)

%% Correction terms %%

% linear term
Mlu =2;
sigma1 = [1-alpha, 1, 3 * alpha, 4 * alpha, 5 * alpha];

% force term
Mf = 2;
sigma2 = [1-alpha, 1, 3 * alpha, 4 * alpha, 5 * alpha];

% history term
Mu = 1;
sigma3 = [1, 2 - alpha, 3 * alpha, 4 * alpha, 5 * alpha];

% force term linearization
Me = 2;
sigma4 = [1-alpha, 1, 3 * alpha, 4 * alpha, 5 * alpha];

% Analytical solution powers
Nterms = 4;
sigma = [1, 2, 3, 4];

% Initializing arrays
N     = zeros(Ntests, 1);    % Number of time steps
dt    = zeros(Ntests, 1);    % Time-step size
Err   = zeros(Ntests, 1);    % Error array
pErr  = zeros(Ntests, 1);    % Convergence rate array
ctime = zeros(Ntests, 1);    % CPU time array

% Loop over the number of tests
for k=1:Ntests
        
    % Computing time-step size and number of time steps
    dt(k) = 2^-(5+k);
    N(k)  = T/dt(k);
    
    % IMEX Solver selection
    if imex_type == 1
        [t, y11, ctime(k)] = IMEX_I(@fun_f1,t0,T, y0,N(k),alpha,lambdamrho, ...
            Mlu,sigma1,Mf,sigma2,Mu,sigma3,Me, sigma4,tol,epsil);

    elseif imex_type == 2
        [t, y11, ctime(k)] = IMEX_II(@fun_f1,t0,T, y0,N(k),alpha,lambdamrho, ...
            Mlu,sigma1,Mf,sigma2,Mu,sigma3,Me, sigma4,tol,epsil);

    end

    % Evaluating known analytical solution
    y1 = y0*ones(1, size(y11,2));
    % Power-law form of analytical solution
    for i=1:Nterms
        y1 = y1 + t.^sigma(i);
    end
    
    % Computing the global L^infinity error
    Err(k) = norm(y1-y11, Inf)/norm(y1,Inf);
    
    % Computing the convergence rate
    if k > 1
        pErr(k) = log2(Err(k-1)/Err(k));    
    end
    
end

% Printing the results
Data_Fast = [N, dt, Err, pErr, ctime]

%% Right-hand-side definition %%
function f = fun_f1(t,y,alpha, lambda)

    % Fabricated RHS from defined analytical solution
    % Analytical solution powers
    Nterms = 4;
    sigma = [1, 2, 3, 4];  
    % Coefficient from nonlinear term
    rho = 0.1;
                
    % Fractional derivative of fabricated solution
    dut = 0;    
    for i=1:Nterms
        dut = dut + (gamma(sigma(i)+1)/gamma(sigma(i)+1-alpha))*t^(sigma(i)-alpha);
    end
    
    % Linear term of fabricated solution
    ut = 1;
    for i=1:Nterms
        ut = ut + t^sigma(i);
    end
    
    % Analytical RHS
    ft = dut - lambda * ut + rho*ut^2;
    
    % Constructed RHS
    f = -rho*y^2 + ft;

end






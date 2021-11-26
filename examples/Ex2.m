%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ex2.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Yongtao Zhou, Jorge Suzuki, Mohsen Zayernouri
% 
% Description: Solution of a stiff FDE system (Example 2) from
% "Implicit-explicit time-integration of nonlinar fractional differential
% equations", Applied Numerical Mathematics 156 (2020). The problem is
% defined as the following system of FDEs:
%
%   C_Du(t) = A*u(t) + B*u(t) + g(t),
%
%   g(t) given for a fabricated u(t) = [1 + a1*t^s1+ a2*t^s2;
%                                       1 + a3*t^s3+ a4*t^s4;
%                                       1 + a5*t^s5+ a6*t^s6];
%
%   where A and B are stiff/nonstiff coefficient matrices.
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
T  = 10;              % Final time
y0 = [1, 1, 1];      % Initial condition   
alpha = 0.3;         % Fractional order
epsil = 0.5*10^(-8); % Epsilon value for fast Topelitz inversion
tol = 5e-7;          % Numerical tolerance for the Picard iteration (nonlinear cases)

Ntests = 5;          % Number of runs to be performed at different dt
lambdamrho = 0;      % linear term coefficient (lambda <= 0)
imex_type = 1;       % Type of IMEX solver (1 = first order, 2 = second order)

%% Correction terms %%
% Number of correction terms
NCorr = 3;

% linear term
Mlu = NCorr;
sigma1 = [2 * alpha, 1 + alpha, 5 * alpha, 2, 2 + alpha] - alpha;

% force term
Mf = NCorr;
sigma2 = [2 * alpha, 1 + alpha, 5 * alpha, 2, 2 + alpha] - alpha;

% history term
Mu = NCorr;
sigma3 = [alpha, 2 * alpha, 1 + alpha, 5 * alpha, 2, 2 + alpha];

% force term linearization
Me = NCorr;
sigma4 = [2 * alpha, 1 + alpha, 5 * alpha, 2, 2 + alpha] - alpha;

%% Fabricated solution information %%
% Coefficients
a = [0.5, 0.8, 1, 1, 1, 1];
% Powers
sigma = [alpha, 2 * alpha, 1 + alpha, 5 * alpha, 2, 2 + alpha];

% Coefficient matrix for linear system
A = 0.001*[-1000,    0,    1 ;
         -0.5 ,-0.8, -0.2 ;
         1     ,    0,   -1];

% Initializing arrays
N     = zeros(Ntests, 1);    % Number of time steps
dt    = zeros(Ntests, 1);    % Time-step size
Err   = zeros(Ntests, 1);    % Error array
pErr  = zeros(Ntests, 1);    % Convergence rate array
ctime = zeros(Ntests, 1);    % CPU time array

% Loop over the number of tests
for k=1:Ntests
    
    % Computing time-step size and number of time steps
    dt(k) = 2^-(2+k);
    N(k)  = T/dt(k);
        
    % IMEX Solver selection
    if imex_type == 1
        [t, y11, ctime(k)] = Fast_FAMMI_E_Correction_A_rev2_new(@fun_f1,...
            t0,T,y0,N(k),alpha,lambdamrho,A,Mlu,sigma1,Mf,sigma2,Mu,sigma3, ...
            Me,sigma4,tol,epsil);

    elseif imex_type == 2
        [t, y11, ctime(k)] = Fast_FAMMII_E_Correction_A_rev2(@fun_f1, ...
            t0,T,y0,N(k),alpha,lambdamrho,A,Mlu,sigma1,Mf,sigma2,Mu,sigma3, ...
            Me,sigma4,tol,epsil);
    end
    
    % Analytical solution
    y1 = [a(1)*t.^sigma(1) + a(2)*t.^sigma(2) + 1, ...
          a(3)*t.^sigma(3) + a(4)*t.^sigma(4) + 1, ...
          a(5)*t.^sigma(5) + a(6)*t.^sigma(6) + 1];

    % Computing the global L^infinity error
    Err(k) = norm(y1-y11, Inf)/norm(y1,Inf);
        
    % Computing the convergence rate
    if k > 1
        pErr(k) = log2(Err(k-1)/Err(k));    
    end
    
end

% Printing the results
Data = [N, dt, Err, pErr, ctime]

%% Right-hand-side
function f=fun_f1(t,y,alpha,lambda)
     
    A = 0.001*[-1000,    0,    1 ;
         -0.5 ,-0.8, -0.2 ;
         1     ,    0,   -1];
     
    B = 0.01*[-0.6  ,    0,  0.2 ;
         -0.1  , -0.2,    0 ;
            0  , -0.5, -0.8];
        
    sigma = [alpha, 2 * alpha, 1 + alpha, 5 * alpha, 2, 2 + alpha];
    a = [0.5, 0.8, 1, 1, 1, 1];
    
    g1 = [a(1)*Gammak(sigma(1),alpha)*t^(sigma(1)-alpha) + a(2)*Gammak(sigma(2),alpha)*t^(sigma(2)-alpha);
          a(3)*Gammak(sigma(3),alpha)*t^(sigma(3)-alpha) + a(4)*Gammak(sigma(4),alpha)*t^(sigma(4)-alpha);
          a(5)*Gammak(sigma(5),alpha)*t^(sigma(5)-alpha) + a(6)*Gammak(sigma(6),alpha)*t^(sigma(6)-alpha)];
    
    g2 = [a(1)*t^sigma(1) + a(2)*t^sigma(2) + 1;
          a(3)*t^sigma(3) + a(4)*t^sigma(4) + 1;
          a(5)*t^sigma(5) + a(6)*t^sigma(6) + 1];
      
    g = g1 - (A + B) * g2;
    
    f = (B*y' + g)';

end

% auxiliary function
function [Gk] = Gammak(sigma, alpha)

    Gk = gamma(sigma+1)/gamma(sigma+1-alpha);

end





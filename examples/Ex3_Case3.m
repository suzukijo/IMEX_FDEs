%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ex3_Case3.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Yongtao Zhou, Jorge Suzuki, Mohsen Zayernouri
% 
% Description: Solution of a nonlinear FDE (Example 3, case 3) from
% "Implicit-explicit time-integration of nonlinar fractional differential
% equations", Applied Numerical Mathematics 156 (2020). The problem is
% defined as:
%
%   C_Du(t) = lambda*u(t) + f(u(t),t),
%
%   lambda = -1, f(u(t),t) = 0.01*u*(1-u^2) + 2*cos(2*pi*t),
%   with u(0) = 1
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
%% Parameters %%
t0 = 0;              % Initial time
T  = 50;             % Final time
dt = 2^-10;          % Time-step size
N = T/dt;            % Number of time steps
y0 = 1;              % Initial condition   
alpha = 0.7;         % Fractional order
epsil = 0.5*10^(-8); % Epsilon value for fast Topelitz inversion
lambdamrho = -1;     % linear term coefficient (lambda <= 0)
tol = 1e-6;          % Numerical tolerance for the Picard iteration (nonlinear cases)
imex_type = 1;       % Type of IMEX solver (1 = first order, 2 = second order)

%% Correction terms %%
M = 3;      % Number of terms
Mlu = M;    % linear term
sigma1 = [alpha, 2 * alpha, alpha+1, 4 * alpha, 5 * alpha];
Mu = M;     % history term
sigma2 = sigma1;    
Me = M;     % force term linearization
sigma3 = sigma1;
Mf = M;     % force term
sigma4 = sigma1;

%% IMEX Solver %%
if imex_type == 1
    [t, y11, ctime] = IMEX_I(@fun_f1,t0,T, y0,N,alpha,lambdamrho, ...
        Mlu,sigma1,Mf,sigma2,Mu,sigma3,Me,sigma4, tol,epsil);

elseif imex_type == 2
    [t, y11, ctime] = IMEX_II(@fun_f1,t0,T,y0, N,alpha,lambdamrho, ...
        Mlu,sigma1,Mf,sigma2,Mu,sigma3,Me,sigma4, tol,epsil);
end

%% Plotting the data %%
f1 = figure
set(gca,'FontSize',16);
hold on;
axis([0, T, -1, 2]);
xlabel('Time [s]');
ylabel('Solution u(t)');
p1 = plot(t, y11, '-b','LineWidth', 1);
% create a new pair of axes inside current figure
axes('position',[.5 .7 .3 .2])
box on % put box around new pair of axes
minplot = (T/2)/dt;
maxplot = (T/2+5)/dt;
plot(t(minplot:maxplot,1),y11(minplot:maxplot,1),'-b') % plot on new axes
axis tight

%% Right-hand-side definition
function f=fun_f1(t,y,alpha, lambda)
            
    f = 0.01*y * (1 - y^2) + 2*cos(2*pi*t);

end

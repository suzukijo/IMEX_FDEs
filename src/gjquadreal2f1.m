function [h]=gjquadreal2f1(a,b,c,zr,zi,n)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 2F1(a,b;c;z) for real    %
% a,b,c,z by applying Gauss-Jacobi quadrature on (4.8), as in Section 4.4.%

%-------------------------------------------------------------------------%
% Input:  a=Parameter value a                                             %
%         b=Parameter value b                                             %
%         c=Parameter value c                                             %
%         zr=Re(z)                                                        %
%         zi=Im(z)                                                        %
%         n=Number of mesh points, denoted as Nmesh in Section 3.6        %
% Output: h=Computed value of 1F1(a;b;z)                                  %
%-------------------------------------------------------------------------%

% Compute z in terms of zr,zi
z=zr+zi*1i;

% Apply code qrule.m [72] to integral (4.8) to find nodes and weights
[x,w]=qrule(n,7,c-b-1,b-1);

% Compute relevant Gamma functions
e1=gamma(b);e2=gamma(c);e3=gamma(c-b);

% Apply x,w,e1,e2,e3 to (4.8)
h=e2/e1/e3/(2^(c-1))*sum(w.*(1-z/2*(x+1)).^(-a));

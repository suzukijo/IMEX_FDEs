function [t, y, ctime]=IMEX_II(Rfun,t0,T,y0,N,alpha,lambda,Mlu,sigma1,Mf,sigma2,Mu,sigma3,Me,sigma4,tol,epsil)

% IMEX_II Second-order IMEX solver for nonlinear FDEs
%   IMEX_II solver a nonlinear FDE of the form:
%
%   {}^C_{0} D_t^\alpha u(t) = \lambda u(t) + f(t,u(t)), \alpha \in (0,1)
%
% Inputs:
%
%  Rfun: The right-hand side of the equation
%  t0: initial point
%  T: end point
%  y0: the value of the function at the initial point
%  N: number of step
%  alpha: the order of the fractional differential operator
%  lambda: a constant
%  Mlu: number of correction terms for \lambda u 
%  sigma1: singular powers for \lambda u
%  Mf: number of correction terms for the force 
%  sigma2: singular powers for the force
%  Mu: number of correction terms for the history
%  sigma3: singular powers for the history
%  Me: number of correction terms for f_(k+1) by extrapolation
%  sigma4: singular powers for f_(k+1)
%  tol: tolerate of the iteration algorithm
%  epsil: small number to construct the epsil-circulant matrix
%
% Outputs:
%
% t: the time grid where the solution is computed
% y: solution vector u(t)
% ctime: total computational time, including computation of correction
% terms.
%
% Revision history:
%
% Rev1 (September 2019) - Initial file
% Rev2 (November 2021)  - Added descriptions for fast inversion procedure
%

if nargin<17  
   epsil=0.5*10^(-8);
end

tstart = tic;

% Time-step size and allocating time/solution arrays
h=(T-t0)/N;
t=(t0:h:T)';
d=length(y0);
y=zeros(N+1,d);
y(1,:)=y0;

%%%%  Computing A and B coefficients for history load, which involves 
%%%% the computation of hypergeometric functions 
%%%%  a_{k,j}=\int_{t_k}^{t_k+1} (t_k+1 - s)^{alpha-1} (s - t_j)^{1-alpha} ds

% Gauss-Jacobi quadrature for hypergeometric functions
a(1,1)=1/alpha*gjquadreal2f1(alpha-1,1,alpha+1,-1,0,200);
a(1,2)=gamma(alpha)*gamma(2-alpha)/gamma(2);

% Evaluating coefficients and employing symmetry relationships
for k=2:3
    a(k,1)=k^(1-alpha)*1/alpha*gjquadreal2f1(alpha-1,1,alpha+1,-1/k,0,200);
    a(k,2)=a(k-1,1);
    a(k,3)=a(k-1,2);
end
a(3,4)=a(1,2);
for k=4:N
    a(k,1)=k^(1-alpha)*1/alpha*gjquadreal2f1(alpha-1,1,alpha+1,-1/k,0,200);
    a(k,2)=a(k-1,1);
    a(k,3)=a(k-1,2);
    a(k,4)=a(k-1,3);
    a(k,5)=a(k-1,4);
end
b(2,1)=1/alpha*gjquadreal2f1(alpha-2,1,alpha+1,-1,0,200);
b(2,2)=gamma(alpha)*gamma(3-alpha)/gamma(3);
b(3,1)=2^(2-alpha)*1/alpha*gjquadreal2f1(alpha-2,1,alpha+1,-1/2,0,200);
b(3,2)=b(2,1);
b(3,3)=b(2,2);
for k=4:N
    b(k,1)=(k-1)^(2-alpha)*1/alpha*gjquadreal2f1(alpha-2,1,alpha+1,1/(1-k),0,200);
    b(k,2)=b(k-1,1);
    b(k,3)=b(k-1,2);
    b(k,4)=b(k-1,3);
end

%% Checking for correction terms
if Mlu>0 || Mf>0 || Mu>0 || Me>0

   %%%% Computing the weights for lambda u
   V1 = zeros(Mlu);
   for i=1:Mlu
       for k=1:Mlu
           V1(i,k) = i^sigma1(k);
       end
   end
   rhs1 = zeros(N,Mlu);
   Flu=zeros(N,Mlu);  
   auxgamma1 = gamma(alpha+1);
   auxgamma2 = gamma(alpha+2);
   for k=1:Mlu
       Flu(1,k)=gamma(sigma1(k)+1)/gamma(sigma1(k)+alpha+1) - 1/auxgamma2;
       for j=2:N
           Flu(j,k)=(j-1)^sigma1(k)/auxgamma1 * gjquadreal2f1(-sigma1(k),1,alpha+1,-1/(j-1),0,200) ...
               - alpha*(j-1)^sigma1(k)/auxgamma2 - j^sigma1(k)/auxgamma2;
       end
       rhs1(:,k)=Flu(:,k);
   end
   Wlu = rhs1*inv(V1); % Obtaining correction weights for linear term
   
   %%%% Computing the weights for the force
   V2 = zeros(Mf);
   for i=1:Mf
       for k=1:Mf
           V2(i,k) = i^sigma2(k);
       end
   end
   rhs2 = zeros(N,Mf);
   Force=zeros(N,Mf);  
   for k=1:Mf
       sigma2k = sigma2(k);
       Force(1,k)=gamma(sigma2k+1)/gamma(sigma2k+alpha+1) - 1/auxgamma2;
       for j=2:N
           Force(j,k)=(j-1)^sigma2k/auxgamma1 * gjquadreal2f1(-sigma2k,1,alpha+1,-1/(j-1),0,200) ...
               - alpha*(j-1)^sigma2k/auxgamma2 - j^sigma2k/auxgamma2;
       end
       rhs2(:,k)=Force(:,k);
   end
   Wf = rhs2*inv(V2); % Getting correction weights
   
   %%%% Computing the weights for the history load
   V3 = zeros(Mu);
   for i=1:Mu
       for k=1:Mu
           V3(i,k) = i^sigma3(k);
       end
   end
   rhs3 = zeros(N,Mu);
   History=zeros(N,Mu);  
   r=zeros(1,N);
   r1=zeros(1,N);
   r1(2)=a(1,1)-a(1,2);
   r1(3)=a(2,1)-a(2,2)+2*a(2,3) - 2*(b(2,1) - b(2,2))/(2-alpha);
   for i=4:N
       r1(i)=(2*a(i-1,1)-2*a(i-1,2)+3*a(i-1,3)-a(i-1,4))/2 + (-2*b(i-1,1)+3*b(i-1,2)-b(i-1,3))/(2-alpha);
   end
   r(3)=(a(2,2)-3*a(2,3))/2 + (b(2,1)-b(2,2))/(2-alpha);
   r(4)=(a(3,2)-3*a(3,3)+4*a(3,4))/2 + (b(3,1)-3*b(3,2)+2*b(3,3))/(2-alpha);
   for i=5:N
       r(i)=(a(i-1,2)-3*a(i-1,3)+3*a(i-1,4)-a(i-1,5))/2 + (b(i-1,1)-3*b(i-1,2)+3*b(i-1,3)-b(i-1,4))/(2-alpha);
   end
   
   auxgamma = gamma(alpha)*gamma(2-alpha);

   for k=1:Mu
       
       sigma3k = sigma3(k);
       auxbeta = beta(sigma3k-alpha+1,alpha);
       
       % Performing fast history computation through FFT
       raux = [r(3:N)'; 0; zeros(N-3,1)]; % Circulant embedding (vector of unique coefficients)
       jaux = [((2:N-1).^sigma3k)'; zeros(N-2,1)];    
       res = [0; real(ifft( fft(raux).*fft(jaux) ))];
              
       auxgamma1 = gamma(alpha)*gamma(sigma3k-alpha+2);
       
       for j=2:N

           FF3=(r1(j)+res(j-1,1))/auxgamma;
                        
           History(j,k)=- gamma(1+sigma3k)*j^(alpha-1) ...
               * ((sigma3k-alpha+1)*(j)^(sigma3k-alpha+1)*betainc((j-1)/j,sigma3k-alpha+1,alpha)*auxbeta)/auxgamma1 ...
               + (j-1)^sigma3k - FF3;

       end

       rhs3(:,k)=History(:,k);

   end

   Wu = rhs3*inv(V3); % Correction weights for history term
    
   %%%% Computing the weights for f_(k+1) by extrapolation 
   V4 = zeros(Me);
   for i=1:Me
       for k=1:Me
           V4(i,k) = i^sigma4(k);
       end
   end
   rhs4 = zeros(N,Me);
   F4=zeros(N,Me); 
   for k=1:Me
       sigma4k = sigma4(k);
       for j=2:N
           F4(j,k)=j^sigma4k - 2*(j-1)^sigma4k + (j-2)^sigma4k;
       end
       rhs4(:,k)=F4(:,k);
   end
   We = rhs4*inv(V4);
    
   %%%%  computing y_1,y_2,...,y_N
   %%%%  (A-lambda*B)*Y=h^alpha*(B + C)*F + D (Eq 4.1 in Zhou et al.)

   auxgamma = gamma(alpha)*gamma(2-alpha);
   A1=zeros(N,1);
   A1(1)=1;
   A1(2)=-1+( (a(2,2)-3*a(2,3))/2 + (b(2,1)-b(2,2))/(2-alpha) )/auxgamma;
   A1(3)=( (a(3,2)-3*a(3,3)+4*a(3,4))/2 + (b(3,1)-3*b(3,2)+2*b(3,3))/(2-alpha) )/(gamma(alpha)*gamma(2-alpha));
   for i=4:N
       A1(i)=( (a(i,2)-3*a(i,3)+3*a(i,4)-a(i,5))/2 + (b(i,1)-3*b(i,2)+3*b(i,3)-b(i,4))/(2-alpha) )/auxgamma;
   end
   
   A2=zeros(N,1);
   A2(1)=1;
   A2(2)=-1+1/auxgamma*(a(1,1)-a(1,2));
   A2(3)=1/auxgamma*(a(2,1)-a(2,2)+2*a(2,3) - (2*b(2,1) - 2*b(2,2))/(2-alpha) );
   for i=4:N
       A2(i)=1/auxgamma*((2*a(i-1,1)-2*a(i-1,2)+3*a(i-1,3)-a(i-1,4))/2 + (-2*b(i-1,1)+3*b(i-1,2)-b(i-1,3))/(2-alpha));
   end

   A2=-A2 + A1;
   
   % B1 is stored as a sparse matrix
   B1col=h^alpha/gamma(alpha+2)*[1;alpha;zeros(N-2,1)];
   A1(1:2) = A1(1:2) - lambda * B1col(1:2);
   
   % Sparse Toeplitz allocation
   B1 = sptoeplitz(B1col, [B1col(1,1), zeros(1,N-1)]);
   
   % Similar procedure for B2
   B2col=h^alpha/gamma(alpha+2)*[-1;2;-1;zeros(N-3,1)];
   B2 = sptoeplitz(B2col, [B2col(1,1), zeros(1,N-1)]);
   B2(1,1)=0;
   
   % Vector with initial conditions and force
   C=zeros(N,d);
   C(1,:)=y0+alpha*h^alpha/gamma(alpha+2)*(lambda*y0 + feval(Rfun,t0,y0,alpha,lambda));
   C(2,:)=1/(auxgamma)*y0*(a(1,1)-a(1,2)) - h^alpha/gamma(alpha+2)*feval(Rfun,t0,y0,alpha,lambda);
   for i=3:N
       C(i,:)=1/(auxgamma)*y0*(a(i-1,1)-a(i-1,2)) ...
           -1/(auxgamma)*y0*(-(a(i-1,2)+a(i-1,3))/2 + (b(i-1,1)-b(i-1,2))/(2-alpha));
   end
   
   % Diagonal matrix for fast inverse approximation
   %%%% D_delt = diag(1,delt,...,delt^(N-1))
   delt=(epsil)^(1/N);
   D=zeros(N,1);
   for i=1:N
       D(i)=delt^(i-1);
   end

   % Computing Lambda = FN*D*K for fast Toeplitz inversion
   lmt=fft(D.*A1);
   
   %%%% Picard iteration algorithm
   auxgamma = gamma(alpha+2);
   e=ones(N,d);             % Residual vector
   F=zeros(N+1,d);          %%%% no need for circulant embedding here
   FC=zeros(N,d);           % Vector with correction terms
   y1=kron(ones(N,1),y(1,:));
   countpicard = 0;

   while (norm(e)>tol)

       Y=[y(1,:);y1];    %%%% initial iteration value set as initial condition

       for i=1:N+1
           F(i,:)=feval(Rfun,t(i),Y(i,:),alpha,lambda);
       end
       
       % Vector with correction terms
       FC(1:N,:)=lambda*h^alpha*Wlu(1:N,1:Mlu)*(Y(2:Mlu+1,:) - y0) ...
           + h^alpha/auxgamma*We(1:N,1:Me)*(F(2:Me+1,:)-F(1,:)) ...
           + h^alpha*Wf(1:N,1:Mf)*(F(2:Mf+1,:)-F(1,:)) ...
           - Wu(1:N,1:Mu)*(Y(2:Mu+1,:) - y0);

       % Right-hand-side vector with sparse Toeplitz Matvec, vector of
       % initial conditions, and vector with correction terms
       b2 = B1*F(2:end,:) + B2*F(2:end,:) + C + FC + A2.*Y(2,:);

       % FFT of rhs for fast inversion (Eq.4.8 from Zhou et al.)
       b1=fft(D.*b2);

       % diag(Lambda^-1)*(FN*D)*R
       b11=b1./lmt;

       % Solution for current Picard iteration
       y2=real(ifft(b11)./D);

       % Residual
       e=y2-y1;

       % Updating vectors
       y1=y2;

       countpicard = countpicard + 1;
   end

   % Attributing converged solution to output vector
   for i=2:N+1
       y(i,:)=y1(i-1,:);
   end

else

   %% No correction terms

   %%%%  computing y_1,y_2,...,y_N
   %%%%  (A-lambda*B)*Y=h^alpha*(B + C)*F (Eq 4.1 in Zhou et al. without D)
   auxgamma = gamma(alpha)*gamma(2-alpha);
   A1=zeros(N,1);
   A1(1)=1;
   A1(2)=-1+( (a(2,2)-3*a(2,3))/2 + (b(2,1)-b(2,2))/(2-alpha) )/auxgamma;
   A1(3)=( (a(3,2)-3*a(3,3)+4*a(3,4))/2 + (b(3,1)-3*b(3,2)+2*b(3,3))/(2-alpha) )/auxgamma;
   for i=4:N
       A1(i)=( (a(i,2)-3*a(i,3)+3*a(i,4)-a(i,5))/2 + (b(i,1)-3*b(i,2)+3*b(i,3)-b(i,4))/(2-alpha) )/auxgamma;
   end
   
   A2=zeros(N,1);
   A2(1)=1;
   A2(2)=-1+1/auxgamma*(a(1,1)-a(1,2));
   A2(3)=1/auxgamma*(a(2,1)-a(2,2)+2*a(2,3) - (2*b(2,1) - 2*b(2,2))/(2-alpha) );
   for i=4:N
       A2(i)=1/auxgamma*((2*a(i-1,1)-2*a(i-1,2)+3*a(i-1,3)-a(i-1,4))/2 + (-2*b(i-1,1)+3*b(i-1,2)-b(i-1,3))/(2-alpha));
   end

   A2=-A2 + A1;
  
   % B1 is stored as a sparse matrix
   B1col=h^alpha/gamma(alpha+2)*[1;alpha;zeros(N-2,1)];
   A1(1:2) = A1(1:2) - lambda * B1col(1:2);

   % Sparse Toeplitz allocation
   B1 = sptoeplitz(B1col, [B1col(1,1), zeros(1,N-1)]);
   
   % Similar procedure for B2
   B2col=h^alpha/gamma(alpha+2)*[-1;2;-1;zeros(N-3,1)];
   B2 = sptoeplitz(B2col, [B2col(1,1), zeros(1,N-1)]);
   B2(1,1)=0;
   
   % Vector with initial conditions and force
   C=zeros(N,d);
   C(1,:)=y0+alpha*h^alpha/gamma(alpha+2)*(lambda*y0 + feval(Rfun,t0,y0,alpha,lambda));
   C(2,:)=1/auxgamma*y0*(a(1,1)-a(1,2)) - h^alpha/gamma(alpha+2)*feval(Rfun,t0,y0,alpha,lambda);
   for i=3:N
       C(i,:)=1/auxgamma*y0*(a(i-1,1)-a(i-1,2)) ...
           -1/auxgamma*y0*(-(a(i-1,2)+a(i-1,3))/2 + (b(i-1,1)-b(i-1,2))/(2-alpha));
   end
   
   % Diagonal matrix for fast inverse approximation
   %%%% D_delt = diag(1,delt,...,delt^(N-1))
   delt=(epsil)^(1/N);
   D=zeros(N,1);
   for i=1:N
       D(i)=delt^(i-1);
   end

   % Computing Lambda = FN*D*K for fast Toeplitz inversion
   lmt=fft(D.*A1);
   
   %%%% Picard iteration algorithm
   e=ones(N,d);             % Residual vector
   F=zeros(N+1,d);          %%%% no need for circulant embedding here
   y1=kron(ones(N,1),y(1,:));
   countpicard = 0;

   while (norm(e)>tol)
       Y=[y(1,:);y1];    %%%% initial iteration value set as initial condition

       for i=1:N+1
           F(i,:)=feval(Rfun,t(i),Y(i,:),alpha,lambda);
       end

       % Right-hand-side vector with sparse Toeplitz Matvec and vector of
       % initial conditions
       b2 = B1*F(2:end,:) + B2*F(2:end,:) + C + A2.*Y(2,:);

       % FFT of rhs for fast inversion (Eq.4.8 from Zhou et al.)
       b1=fft(D.*b2);

       % diag(Lambda^-1)*(FN*D)*R
       b11=b1./lmt;

       % Solution for current Picard iteration
       y2=real(ifft(b11)./D);

       % Residual vector
       e=y2-y1;

       % Updating solution
       y1=y2;

       countpicard = countpicard + 1;
   end

   % Attributing converged solution to output vector
   for i=2:N+1
       y(i,:)=y1(i-1,:);
   end

end

ctime = toc(tstart);

function [a,b,d,e,f,x0]=dsr_oe(y,u,L,g,k,bmet,n)
% DSR_OE Deterministic system identification and Realization
%      [A,B,D,E,F,x0]=dsr_oe(Y,U,L)
%      [A,B,D,E,F,x0]=dsr_oe(Y,U,L,g)
%      PURPOSE:
%      Estimate the system order (n) and the matrices (A,B,D,E,F)
%      in the following discrete time output error model.
%
%      x_{t+1} = A x_t + B u_t,          x_{t=0}=x0
%      y_t     = D x_t + E u_t + e_t
%
%      where C     = CF*inv(F)       Kalman gain matrix, (C=inv(A) K))
%      Delta = E(e_t e_t^T) = F*F'   Innovations noise covariance matrix
%
%      ON INPUT:
%      Y         - Output time series matrix of size (N x m) 
%                  where N is the number of observations and
%                  m is the number of output variables.
%      U         - Input time series matrix of size (N x r)
%                  where r is the number of input variables.
%      L         - Number of block rows in extended observability matrix.
%                  Choose L .geq. 1. This means that one can estimate
%                  system order (n) bounded by, n .leq. L*m.
%      g         - g=1 default, g=0 force E to be the zero matrix.
%      ON OUTPUT:
%      A,B,D,E   - model system matrices
%      F         - Delta = F*F', The innovation noise covariance matrix
%      x0        - initial values for the state vector, x_t, i.e. state at t=0.
%
%                                       COPYRIGHT 1996, FANTOFT PROCESS
%                                       License belong to Terje Karstang
%                                       Product id: 10 0000
%------------------------------------------------------------------------

% Notes: 
% 1. Choose L as close to the observability index (n-d+1) as possible
%    L >= n-d+1, n >=d whenever the input is poor with frequencies
%    where d=rank(y_t) usually equal to m.
%    (however, this is not necessary)
% 2. In case that E=0, then chose parameter g=0.
% 3. k, bmet, n. Optional parameters for advanced use.
%    (k=L default, bmet=1
% 4. Algorithm: Di Ruscio (1996), A Method for ...,
%               In "Computer Aided Time series Modeling",
%               Ed: M. Aoki, Springer Verlag.
%------------------------------------------------------------------------

% low level functions: seye, sobsv, simpr, besolv, dread
if nargin == 6; n = 1; end
if nargin == 5; n = 1; bmet = 1; end
if nargin == 4; n = 1; bmet = 1; k = L; end
if nargin == 3; n = 1; bmet = 1; k = L; g = 1; end
if nargin == 2; n = 1; bmet = 1; k = 1; g = 1; L=1; end

if n <= 0
 n=1; 
 disp('DSR warning: choose system order n > 0')
end
if L <= 0
 L=1; k=1;
 disp('DSR warning: choose horizon L > 0')
end

[Ny,ny] = size(y);
[Nu,nu] = size(u);

k=0;

N  = min(Ny,Nu);
K  = N - L - k; 

% 1. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 7
 disp('Ordering the given input output data')
end
YL = zeros((L+1)*ny,K); UL = zeros((L+g)*nu,K);
for i=1:L+1
    YL(1+(i-1)*ny:i*ny,:) = y(i:K+i-1,:)';
end
for i=1:L+g
    UL(1+(i-1)*nu:i*nu,:) = u(i:K+i-1,:)';
end

%2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 7
 disp('QR decomposition')
end

R  = triu(qr([UL;YL]'))';

nr = (L+g)*nu+(L+1)*ny;
E  = diag(sign(diag(R(:,1:nr))));
R  = R(:,1:nr)*E/sqrt(K);

%3. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R11 = R(1:(L+g)*nu,1:(L+g)*nu);
R21 = R((L+g)*nu+1:(L+g)*nu+L*ny,1:(L+g)*nu);
R31 = R((L+g)*nu+L*ny+1:(L+g)*nu+(L+1)*ny,1:(L+g)*nu);
R31 = [R21(ny+1:L*ny,:);R31];

R22 = R((L+g)*nu+1:(L+g)*nu+L*ny,(L+g)*nu+1:(L+g)*nu+L*ny);
R32 = R((L+g)*nu+L*ny+1:(L+g)*nu+(L+1)*ny,(L+g)*nu+1:(L+g)*nu+L*ny);
R33 = R((L+g)*nu+L*ny+1:(L+g)*nu+(L+1)*ny,(L+g)*nu+L*ny+1:(L+g)*nu+(L+1)*ny);

f=R33;
% 4. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[U,S,V] = svd(R22);
nSi=L*ny;
sn = diag(S(1:nSi,1:nSi));
if nargin < 7
% sn'
 hold off;
 subplot(121);
 semilogy([1:L*ny],sn,'*')
 title('Singular Values')
 xlabel('System order')
 grid
 hold off
 subplot(122);
 semilogy([1:L*ny],sn(1) ./sn,'*')
 title('Condition Numbers')
 grid, subplot(111)
 xlabel('System order')
 n=0;
 for i=1:L*ny; if sn(i)/sn(1) > 1.0e-7; n=n+1; end; end 
 if L == 1 & ny == 1
  n = 1;
 else
  n = dread('System order ?',n);
 end
end
O  = U(:,1:n);
Sn = S(1:n,1:n); 
V1 = V(:,1:n);

a  = O'*[U(ny+1:L*ny,1:n);R32*V1*pinv(Sn)];
d  = O(1:ny,1:n);

%5. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bm=3;
if bm == 1
 p=2;
 On=O(1:p*ny,:);
 r11=R11(1:(p+g)*nu,1:(p+g)*nu);
 r31=R31(1:p*ny,1:(p+g)*nu);
 r41=R41(1:p*ny,1:(p+g)*nu);
 Atn=On*a*pinv(On'*On)*On';
 Z=r41-Atn*r31; U=r11;
 [b,e]=besolv(Z,U,Atn,a,d,p,g,nu);
elseif bm==2
 On=sobsv(a,d,k,1);
 At=On*a*pinv(On'*On)*On';
 Z=R21Y1-At*R21Y0; U=R21U0; % ok
 [b,e]=besolv(Z,U,At,a,d,k,g,nu);
elseif bm==3
 On=O(1:L*ny,:);
% On=sobsv(a,d,L,1);
 At=On*a*pinv(On'*On)*On';
 Z=R31-At*R21; U=R11;
 [b,e]=besolv(Z,U,At,a,d,L,g,nu);
end

%7. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hd=zeros(L*ny,L*nu);
for i=1:L
  Hd(:,(i-1)*nu+1:i*nu)=simpr(a,b,d,e,L,i);
end
x0=pinv(O)*(YL(1:L*ny,1)-Hd*UL(1:L*nu,1));
%
% END DSR

save r R22 R32
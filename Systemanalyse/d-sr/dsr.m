function [a,b,d,e,cf,f,x0,sn]=dsr(y,u,L,g,k,bmet,n)
% DSR  Deterministic and Stochastic system identification and Realization
%      [A,B,D,E,CF,F,x0]=dsr(Y,U,L)
%      [A,B,D,E,CF,F,x0]=dsr(Y,U,L,g)
%      [A,B,D,E,CF,F,x0]=dsr(Y,U,L,g,J,M,n)
%      PURPOSE:
%      Estimate the system order (n) and the matrices (A,B,D,E,CF,F)
%      in the following discrete time combined deterministic and 
%      stochastic dynamic model on innovations form
%
%      x_{t+1} = A x_t + B u_t + C e_t,    x_{t=0}=x0,
%      y_t     = D x_t + E u_t + e_t,
%
%      where C     = CF*inv(F),       (Kalman gain, K =inv(A)*C and C=A*K),
%      Delta = E(e_t e_t^T) = F*F',   (Innovations noise covariance matrix).
%
%      ON INPUT:
%      Y         - Output time series matrix of size (N x m) 
%                  where N is the number of observations and
%                  m is the number of output variables.
%      U         - Input time series matrix of size (N x r)
%                  where r is the number of input variables.
%      L         - Number of block rows in extended observability matrix.
%                  Choose L .geq. 1. This means that one can estimate
%                  system order (n) bounded by, 0 < n .leq. L*m.
%      OPTIONAL INPUT PARAMETERS:
%      g         - g=1 default, g=0 force E to be the zero matrix.
%      J         - Past horizon used to define instruments. Default, J=L.
%      M         - Default M=1. See tutorial.
%      n         - Optional specification of model order, 0 < n .leq. L m.
%      ON OUTPUT:
%      A,B,D,E   - Model system matrices.
%      CF, F     - C=CF*inv(F) is the Kalman filter gain matrix.
%      F         - Delta = F*F' is the innovation noise covariance matrix.
%      x0        - Initial values for the state vector, x_t, i.e. state at t=0.
%
%                                       COPYRIGHT 1996, 1999, DDIR
%                                       License belong to:
%                                       Product id: 10 000
%------------------------------------------------------------------------

% DATE: 9, january 2004. Better method for setting the system order.
%       1. november 1996
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
% 5. J must satisfy J > 0 in order to define the instruments.
%------------------------------------------------------------------------

% low level functions: seye, sobsv, simpr, besolv, dread
if nargin == 6; n = 1; end
if nargin == 5; n = 1; bmet = 1; end
if nargin == 4; n = 1; bmet = 1; k = L; end
if nargin == 3; n = 1; bmet = 1; k = L; g = 1; end
if nargin == 2; n = 1; bmet = 1; k = 1; g = 1; L=1; end

[Ny,ny] = size(y);
if isempty(u)==1; u=zeros(Ny,1); end                   % changed 17.01.2000

[Nu,nu] = size(u);

N  = min(Ny,Nu);
K  = N - L - k; 

% 1. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 7
 disp('Ordering the given input output data')
end
YL = zeros((L+k+1)*ny,K); UL = zeros((L+k+g)*nu,K);
for i=1:L+k+1
    YL(1+(i-1)*ny:i*ny,:) = y(i:K+i-1,:)';
end
for i=1:L+k+g
    UL(1+(i-1)*nu:i*nu,:) = u(i:K+i-1,:)';
end

%2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 7
 disp('QR decomposition')
end

R  = triu(qr([UL(k*nu+1:k*nu+(L+g)*nu,:);UL(1:k*nu,:);YL]'))';
ni = k*nu + k*ny;

nr = (L+k+g)*nu+(L+k+1)*ny;
E  = diag(sign(diag(R(:,1:nr))));
R  = R(:,1:nr)*E/sqrt(K);

%3. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R11 = R(1:(L+g)*nu          ,1:(L+g)*nu);
i21u=(L+g)*nu;
i21y=(L+g)*nu+k*nu;
R21U0=R(i21u+1:i21u+k*nu,1:(L+k+g)*nu);
R21U0=[R21U0;R(1:g*nu,1:(L+k+g)*nu)];

%R21Y0=R(i21y+1:i21y+k*nu,1:(L+k+g)*nu);                 % old
%R21Y1=R(i21y+1+ny:i21y+k*nu+ny,1:(L+k+g)*nu);           % old

R21Y0=R(i21y+1:i21y+k*ny,1:(L+k+g)*nu);                  % Changed 17/3-97
R21Y1=R(i21y+1+ny:i21y+k*ny+ny,1:(L+k+g)*nu);            % Changed 17/3-97

i31=(L+g)*nu+k*nu+k*ny;
R31=R(i31+1:i31+L*ny,1:(L+g)*nu);
R41=R(i31+1+ny:i31+(L+1)*ny,1:(L+g)*nu);

j32=(L+g)*nu;
R32 = R(i31+1:i31+L*ny,j32+1:j32+ni);
i41=(L+g)*nu+k*nu+k*ny+L*ny;
R42 = R(i41+1:i41+ny,j32+1:j32+ni);

R33 = R(i31+1:i31+L*ny,j32+ni+1:j32+ni+L*ny);
R43 = R(i41+1:i41+ny,j32+ni+1:j32+ni+L*ny);
R44 = R(i41+1:i41+ny,i41+1:i41+ny); 

f= R44;         % square root of noise inovation process

% 4. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[U,S,V] = svd(R32);
nSi=min(L*ny,k*(ny+nu));
sn = diag(S(1:nSi,1:nSi));
if nargin < 7
 %sn'
 hold off;
 subplot(121);
% semilogy([1:L*ny],sn,'*')                              % old
 semilogy([1:nSi],sn,'*')                                % Changed 14/5-1997
 title('Singular Values')
 xlabel('System order')
 grid
 hold off
 subplot(122);
% semilogy([1:L*ny],sn(1) ./sn,'*')                      % old
 semilogy([1:nSi],sn(1) ./sn,'*')                        % Changed 14/5-1997
 title('Condition Numbers')
 grid, subplot(111)
 xlabel('System order')
 n=0;
% for i=1:L*ny; if sn(i)/sn(1) > 1.0e-7; n=n+1; end; end % old
% Alternative search for default model order.
log_sn=log(sn);
n_def=min(L,max(find(log_sn>(max(log_sn)+min(log_sn))/2)));
%
for i=1:nSi; if sn(i)/sn(1) > 1.0e-7; n=n+1; end; end   % Changed 14/5-1997
 if L == 1 & ny == 1
  n = 1;
 else
  n = dread('System order ?',n_def);
 end
end
O  = U(:,1:n);
Sn = S(1:n,1:n); 
V1 = V(:,1:n);

a  = O'*[U(ny+1:L*ny,1:n);R42*V1*pinv(Sn)];
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
 Z=R41-At*R31; U=R11;
 [b,e]=besolv(Z,U,At,a,d,L,g,nu);
end

%6. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
At=O*a*pinv(O'*O)*O';
ic=bmet;
if ic == 1
 Hs=[R33 zeros(L*ny,ny);R43 R44];
 OCF=[R33(ny+1:L*ny,1:ny);R43(:,1:ny)];
 cf = pinv(O)*OCF;
else
 r43= R(i31+ny+1:i31+(L+1)*ny,j32+ni+1:j32+ni+L*ny);
 At = O*a*pinv(O'*O)*O';
 Ct = r43 - At*R33;
 Hs = zeros(L*ny,(L+1)*ny);
 Ei = zeros(ny*L,ny); Ei(ny*(L-1)+1:ny*L,:) = f;
 Hs(:,(L-0)*ny+1:(L+1)*ny) = Ei;
 E = zeros(ny,L*ny);
 for i=1:L
   E(1:ny,(i-1)*ny+1:i*ny) = Hs((L-i)*ny+1:(L+1-i)*ny,(L+1-i)*ny+1:(L+2-i)*ny);
   Hs(:,  (L-i)*ny+1:(L+1-i)*ny) = Ct(:,(L-i)*ny+1:(L+1-i)*ny) + At*Ei;
   Ei = Hs(:,  (L-i)*ny+1:(L+1-i)*ny);
   if i < L-1
     Ei(1:(L-i-1)*ny,1:ny) = zeros((L-i-1)*ny,ny);
   end
   Hs(:,(L-i)*ny+1:(L+1-i)*ny) = Ei;
 end
 cf = pinv(O)*Ei;
end

%7. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hd=zeros(L*ny,L*nu);
for i=1:L
  Hd(:,(i-1)*nu+1:i*nu)=simpr(a,b,d,e,L,i);
end
x0=pinv(O)*(YL(1:L*ny,1)-Hd*UL(1:L*nu,1));
%
% END DSR

function [a,d,cf,f,x0,sn]=dsr_s(y,L,k,bmet,n)
% DSR_S Stochastic system identification and Realization
%      [A,D,CF,F,x0]=dsr_s(Y,L)
%      [A,D,CF,F,x0]=dsr_s(Y,L,J)
%      [A,D,CF,F,x0]=dsr_s(Y,L,J,M,n)
%      PURPOSE:
%      Estimate the system order (n) and the matrices (A,D,CF,F)
%      in the following discrete time stochastic dynamic model 
%      on innovations form
%
%      x_{t+1} = A x_t + CF ep_t,  x_{t=0}=x0, E(ep_t ep_t^T)=I,
%      y_t     = D x_t + F ep_t,
%
%      which is equivalent with the standard innovations (prediction) form
%
%      x_{t+1} = A x_t + C e_t,  x_{t=0}=x0,
%      y_t     = D x_t + e_t,
%
%      where C     = CF*inv(F),       (Kalman gain, K =inv(A)*C and C=A*K),
%      Delta = E(e_t e_t^T) = F*F',   (Innovations noise covariance matrix).
%
%      ON INPUT:
%      Y         - Output time series matrix of size (N x m) 
%                  where N is the number of observations and
%                  m is the number of output variables.
%      L         - Number of block rows in extended observability matrix.
%                  Choose L .geq. 1. This means that one can estimate
%                  system order (n) bounded by, 0 < n .leq. L*m.
%      OPTIONAL INPUT PARAMETERS:
%      J         - Past horizon used to define instruments. Default, J=L.
%      M         - Default M=1. See tutorial.
%      n         - Optional specification of model order, 0 < n .leq. L m.
%      ON OUTPUT:
%      A,D       - Model system matrices.
%      CF, F     - C=CF*inv(F) is the Kalman filter gain matrix.
%      F         - Delta = F*F' is the innovation noise covariance matrix.
%      x0        - Initial values for the state vector, x_t, i.e. state at t=0.
%
%                                       COPYRIGHT 1996, 1999, DDIR
%                                       License belong to:
%                                       Product id: 10 000
%------------------------------------------------------------------------

% DATE: 1. november 1996
% Notes: 
% 4. Algorithm: Di Ruscio (1996), A Method for ...,
%               In "Computer Aided Time series Modeling",
%               Ed: M. Aoki, Springer Verlag.
% 5. J must satisfy J > 0 in order to define the instruments.
%------------------------------------------------------------------------

% low level functions: seye, sobsv, dread
if nargin == 4; n = 1; end
if nargin == 3; n = 1; bmet = 1; end
if nargin == 2; n = 1; bmet = 1; k = L; end
if nargin == 1; n = 1; bmet = 1; k = 1; L=1; end

[Ny,ny] = size(y);

N  = Ny;
K  = N - L - k; 

% 1. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 7
 disp('Ordering the given input output data')
end
YL = zeros((L+k+1)*ny,K);
for i=1:L+k+1
    YL(1+(i-1)*ny:i*ny,:) = y(i:K+i-1,:)';
end

%2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 7
 disp('QR decomposition')
end

R  = triu(qr(YL'))';
ni = k*ny;

nr = (L+k+1)*ny;
E  = diag(sign(diag(R(:,1:nr))));
R  = R(:,1:nr)*E/sqrt(K);

%3. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i31 = ni;
i41 = i31+L*ny;

R32 = R(i31+1:i31+L*ny,1:ni);
R42 = R(i41+1:i41+ny,1:ni);

R33 = R(i31+1:i31+L*ny,ni+1:ni+L*ny);
R43 = R(i41+1:i41+ny,ni+1:ni+L*ny);
R44 = R(i41+1:i41+ny,i41+1:i41+ny); 

f= R44;         % square root of noise innovation process

% 4. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[U,S,V] = svd(R32);
nSi=min(L*ny,k*ny);
sn = diag(S(1:nSi,1:nSi));
if nargin < 7
% sn'
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
% for i=1:nSi; if sn(i)/sn(1) > 1.0e-7; n=n+1; end; end   % Changed % % 14/5-1997 % old2
 % Alternative search for default model order.
 log_sn=log(sn);
 n_def=min(L,max(find(log_sn>(max(log_sn)+min(log_sn))/2)));
 
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
x0=pinv(O)*YL(1:L*ny,1);
%
% END DSR

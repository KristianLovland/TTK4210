function [a,b,d,e,cf,f,x0]=mdsr(y,u,L,g,J,bmet,n,Ni)
% [a,b,d,e,cf,f,x0]=mdsr(Y,U,L,g,J,M,n,NN)
% MDSR Multiple DSR. Deterministic and Stochastic system identification 
%      and Realization from possible multiple time series.
%      [A,B,D,E,CF,F,x0]=mdsr(Y,U,L)
%      [A,B,D,E,CF,F,x0]=mdsr(Y,U,L,g)
%      [A,B,D,E,CF,F,x0]=mdsr(Y,U,L,g,J,M,n,NN)
%      PURPOSE:
%      Estimate the system order (n) and the matrices (A,B,D,E,CF,F)
%      in the following discrete time combined deterministic and 
%      stochastic dynamic model on innovations form
%      x_{t+1} = A x_t + B u_t + C e_t,    x_{t=0}=x0,
%      y_t     = D x_t + E u_t + e_t,
%      where C     = CF*inv(F),       (Kalman gain matrix, (C=inv(A) K),
%      Delta = E(e_t e_t^T) = F*F',   (Innovations noise covariance matrix).
%
%      MULTIPLE TIME SERIES:
%      Assume that the following multiple time series are given.
%      Y_i of size (Ni x ny) and U_i of size (Ni x nu) for all i=1,...,Ne.
%      Specify parameters and data-matrices as follows:
%      NN=[N1 N2 ... Ne]
%      Y =[Y1;Y2;...;Y_Ne] an (sum(NN) x ny) data matrix.
%      U =[Y1;Y2;...;Y_Ne] an (sum(NN) x nu) data matrix.
%      Example of call: Given N1=100 samples in (Y1,U1) and N2=200 in (Y2,U2),
%      [a,b,d,e,cf,f,x0]=mdsr([Y1;Y2],[U1;U2],L,1,L,1,[],[100 200]);
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
%                  n=[] --> order identified.
%      NN        - In case of multiple time series. NN=[N1...Ni...Ne].
%      ON OUTPUT:
%      A,B,D,E   - Model system matrices.
%      CF, F     - C=CF*inv(F) is the Kalman filter gain matrix.
%      F         - Delta = F*F' is the innovation noise covariance matrix.
%      x0        - Initial values for the state vector, x_t, i.e. state at t=0.
%                  x0 is a matrix in case of multiple time series, i.e.
%                  when NN is specified as NN=[N1...Ne].
%------------------------------------------------------------------------

% low level functions: seye, sobsv, simpr, besolv, dread

% Written for project PR8-41464.01 by David Di Ruscio. 12/3-98.
% Based on an idea by Terje Karstang of making a 3-dimensional DSR algorithm.

%%%%%%%%%%%%%%%% SPECIFY DEFAULT PARAMETERS AND CHECK INPUTS %%%%%%%%%%%%%%%%%%
mc=0;
if nargin >= 7;  
 if isempty(n) == 1;         % Id./choose order if called with n=[] and L>1
  n=1; 
 else 
  mc=1;              % Order specified. Indicate Monte Carlo experiment,
 end;                % No display etc. on screen.
end

[Ny,ny]=size(y);
[Nu,nu]=size(u);
N=min(Ny,Nu);

if nargin == 7; Ni=N; end
if nargin == 6; Ni=N; n = 1; end
if nargin == 5; Ni=N; n = 1; bmet = 1; end
if nargin == 4; Ni=N; n = 1; bmet = 1; J = L; end
if nargin == 3; Ni=N; n = 1; bmet = 1; J = L; g = 1; end
if nargin == 2; Ni=N; n = 1; bmet = 1; J = 1; g = 1; L=1; end

%% START 0: %%%% COUNT NUMBER OF EXPERIMENTS OR PAIR OF yi, ui TIME SERIES %%%
n_e=length(Ni);                   % n_e equal to # of input and output series.
if n_e==1; Ni=N; end              
%% END 0: %%%%% THE NUMBER OF EXPERIMENTS IS NOW n_e=1,2, ...  %%%%%%%%%%%%%%%

%% START 1: %%%% MAKE CONCATENATED DATA MATRICES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 7  
 disp('Ordering the given input output data')
end

%%% Make space for data matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ks=0;
for ie=1:n_e
 K=Ni(ie)-L-J; Ks=Ks+K;
end
Yx0=zeros(L*ny,n_e); Ux0=zeros(L*nu,n_e);
%YL=zeros((J+L+1)*ny,Ks); UL=zeros((J+L+g)*nu,Ks);    

%%% Make large data matrices YL and UL by DYNAMIC ALLOCATION %%%%%%%%%%%%%%%%%
id=0; YL=[]; UL=[]; 
for ie = 1:n_e                        % make YL=[Y1 ... Yn_e], UL=[U1 ...Un_e] 
 K = Ni(ie) - L - J; 
 Yi=zeros(ny,K); Ui=zeros(nu,K);
 for i=1:L+J+1
     Yi(1+(i-1)*ny:i*ny,1:K) = y(id+i:id+K+i-1,:)';
 end
 for i=1:L+J+g
     Ui(1+(i-1)*nu:i*nu,1:K) = u(id+i:id+K+i-1,:)';
 end
 id=Ni(ie)+id;
 YL=[YL Yi];, UL=[UL Ui]; 
% Yx0(:,ie)=YL(1:L*ny,1+(ie-1)*K); Ux0(:,ie)=UL(1:L*nu,1+(ie-1)*K);
 Yx0(:,ie)=Yi(1:L*ny,1); Ux0(:,ie)=Ui(1:L*nu,1);
end
% END 1. %%%% Data ordered as YL =[Y1 Y2 ... Yn_e] and UL=[U1 U2 ... Un_e] %%%

k=J;                 % Past horizon for instruments.

%size(YL), size(UL), Ks
%save mdsrlog YL UL Ks Yx0 Ux0

%2 %%%%%%%%%%%% IDENTICAL TO ORDINARY DSR IN THE FOLLOWING %%%%%%%%%%%%%%%%%%%%
if nargin < 7
 disp('QR decomposition')
end

R  = triu(qr([UL(k*nu+1:k*nu+(L+g)*nu,:);UL(1:k*nu,:);YL]'))';
ni = k*nu + k*ny;

nr = (L+k+g)*nu+(L+k+1)*ny;
E  = diag(sign(diag(R(:,1:nr))));
R  = R(:,1:nr)*E/sqrt(Ks);

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
%if nargin < 7
if mc==0
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
 for i=1:nSi; if sn(i)/sn(1) > 1.0e-7; n=n+1; end; end   % Changed 14/5-1997
 if L == 1 & ny == 1
  n = 1;
 else
  n = dread('System order ?',n);
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
%x0=pinv(O)*(YL(1:L*ny,1)-Hd*UL(1:L*nu,1));
x0=pinv(O)*(Yx0-Hd*Ux0);
%
% END MDSR

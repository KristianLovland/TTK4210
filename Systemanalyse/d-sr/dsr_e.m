function [a,b,d,e,c,f,x0,E_J1]=dsr_e(y,u,L,g,J,n,im)
%DSR_e  Deterministic and Stochastic system identification and Realization
%       for open and closed loop systems.
%       [A,B,D,E,K,F,x0,Ef1]=dsr_e(Y,U,L,g,J,n)
%       PURPOSE:
%       Estimate the the matrices (A,B,D,E,K,F)
%       in the discrete time combined deterministic and 
%       stochastic dynamic model on innovations form 
%       and the initial state, x0.
%
% ON INPUT:
%       Y         - Output time series matrix of size (N x m) 
%                   where N is the number of observations and
%                   m is the number of output variables.
%       U         - Input time series matrix of size (N x r)
%                   where r is the number of input variables.
%       L         - Number of block rows in extended observability matrix.
%                   Choose L .geq. 1. This means that one can estimate
%                   system order (n) bounded by, 0 < n .leq. L*m.
%       g         - Chose g=0 for closed loop systems, i.e. E=0..
%       J         - Past horizon used to define instruments. Choosen such
%                   that (A-KD)J \approx 0.
%       n         - Model order, 0 < n .leq. L m.
%
%                                        COPYRIGHT 2004, DDIR
%                                        License belong to: unspecified
%                                        Product id: 10 000
% ------------------------------------------------------------------------

% DDIR, 030104.

if nargin==5; im=0; end
if nargin==6; im=0; end

[Ny,ny] = size(y);
[Nu,nu] = size(u);

N  = min(Ny,Nu);
K  = N - J; 
YL = zeros((J+1)*ny,K); UL = zeros((J+1)*nu,K);
for i=1:J+1
    YL(1+(i-1)*ny:i*ny,:) = y(i:K+i-1,:)';
    UL(1+(i-1)*nu:i*nu,:) = u(i:K+i-1,:)';
end

Up=UL(1:J*nu,:);
Uf=UL(J*nu+1:(J+1)*nu,:);
Yp=YL(1:J*ny,:);
Wp=[Up;Yp];
Yf=YL(J*ny+1:(J+1)*ny,:); 

y_J1=Yf(1:ny,:);
u_J1=Uf(1:nu,:);                              % Changed 23/2-05

ie=2;
if ie==1                                      % Using the definitions directly. 
    P1=y_J1*Wp'*pinv(Wp*Wp');
    DXJ=P1*Wp;                                % The signal part, Dx_J, in y_J=Dx_J+e_J.
    E_J1=y_J1-DXJ;                            % The innovations. 
   
elseif ie==2                                  % Using QR decomposition, a fast alternative to the above.
    [Q,R]=qr([Wp;y_J1]',0); 
    Q=Q'; R=R';
    R21=R(J*(nu+ny)+1:J*(nu+ny)+ny, 1:J*(nu+ny));
    R22=R(J*(nu+ny)+1:J*(nu+ny)+ny, J*(nu+ny)+1:J*(nu+ny)+ny);
    Q1=Q(1:J*(nu+ny),:);
    Q2=Q(J*(ny+nu)+1:J*(nu+ny)+ny,:);
    DXJ=R21*Q1;
    E_J1=R22*Q2;                              % May also use, E_J1=Q2;
elseif ie==3                                  % Using pem for testing purpose. poor alternative!
    th=pem([y,u],n,'nk',1);
    [a,b,d,e,k,x0]=th2ss(th);
    yp=dsropt(a,b,d,e,k,eye(ny),y,u,x0);
    E_J1=y(J+1:end,:)-yp(J+1:end,:); E_J1=E_J1';
    DXJ=yp(J+1:end,:)';
end
    

if im==1
    [a,b,d,e,f_oe,x0]=dsr_oe(y_J1',[u_J1' E_J1'],L,g,J,1,n);
    c=b(:,nu+1:nu+ny);                         % The Kalman gain matrix, K
    b=b(:,1:nu);
    f=e(:,nu+1:nu+ny);
    e=e(:,1:nu);
elseif im==0  % Works well 
    if nargin==5
       [a,b,d,e,f_oe,x0]=dsr_oe(DXJ',[u_J1' E_J1'],L,g,J,1);
    else
       [a,b,d,e,f_oe,x0]=dsr_oe(DXJ',[u_J1' E_J1'],L,g,J,1,n);
    end
    c=b(:,nu+1:nu+ny);                         % The Kalman gain matrix, K
    b=b(:,1:nu); 
    f=e(:,nu+1:nu+ny);
    e=e(:,1:nu);
elseif im==3
    ue=u_J1-u_J1*E_J1'*pinv(E_J1*E_J1')*E_J1;
    [a,b,d,e,f_oe,x0]=dsr_oe(DXJ',ue',L,0,J,1,n);
    c=b; f=e; 
elseif im==4 % gives bias
    ue=u_J1*Wp'*pinv(Wp*Wp')*Wp;
    [a,b,d,e,f_oe,x0]=dsr_oe(DXJ',ue',L,0,J,1,n); 
    c=b;
elseif im==2 % Unbiased but best???
    [a,b,d,e,f_oe,x0]=dsr_oe(y_J1'-E_J1',[u_J1' y_J1'],L,0,J,1,n);    % Id of Kalman filter
    c=b(:,nu+1:nu+ny);                         % The Kalman gain matrix, K
    b=b(:,1:nu);
    f=e(:,nu+1:nu+ny);
    e=e(:,1:nu);
    a=a+c*d;         
end
   
% Test
it=0;
if it==1
    Yp=dsropt(a,b,d,e,c,eye(2),y_J1',u_J1',x0);
    E2=y_J1'-Yp;
    [a,b,d,e,f_oe,x0]=dsr_oe(y_J1'-E2,[u_J1' E2],L,g,J,1,n);
%[a,b,d,e,f_oe,x0]=dsr_oe(y_J1',[u_J1' E2],L,g,J,1,n);
    c=b(:,nu+1:nu+ny);
    b=b(:,1:nu);
    f=e(:,nu+1:nu+ny);
    e=e(:,1:nu);
elseif it==2
    Yp=dsropt(a,b,d,e,c,eye(2),y,u,x0);
    E2=y-Yp;
    [a,b,d,e,f_oe,x0]=dsr_oe(Yp,[u E2],L,g,J,1,n);
    c=b(:,nu+1:nu+ny);
    b=b(:,1:nu);
    f=e(:,nu+1:nu+ny);
    e=e(:,1:nu);
elseif it==3  % Simulerer ny innovasjons serie.
    nit=5;
    for itt=1:nit
    ep=dsrsim(a-c*d,[b c],-d,[zeros(ny,nu) eye(ny,ny)],[u y]);       
%    [a,b,d,e,f_oe,x0]=dsr_oe(y,[u ep],L,g,J,1,n);
    [a,b,d,e,f_oe,x0]=dsr_oe(y_J1',[u_J1' ep(J*ny+1:end,:)],L,g,J,1,n);
    c=b(:,nu+1:nu+ny);
    b=b(:,1:nu);
    f=e(:,nu+1:nu+ny);
    e=e(:,1:nu);
    end
end


% The square root of the innovations covariance, the F matrix.
R = triu(qr(E_J1'))';
E = diag(sign(diag(R(:,1:ny))));
f = R(:,1:ny)*E/sqrt(K);
%
% END DSR_e

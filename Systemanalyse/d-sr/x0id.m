function x = x0id(y,u,a,b,d,e,c,L);
% X0ID Computes the initial state for a given model on
% innovations form (Kalman filter) and given data (Y and U.
% SYNTAX
% x1 = x0id(y,u,a,b,d,e,c,L);
% ON INPUT
% Y,U     - Output and input data matrices, respectively.
% a,b,d,e - Deterministic part of model
% c       - The kalman filter gain matrix.
%           Put c=0 for a pure deterministic model.
% L       - The state identification horizon.
% ON OUTPUT
% X1      - The initial state vector.

[Ny,m]=size(y); [Nu,r]=size(u); N=min(Ny,Nu);
[m,r]=size(e);
be=b-c*e;
ae=a-c*d;

Hd=zeros(L*m,L*r);
for i=1:L
  Hd((i-1)*m+1:m*i,(i-1)*r+1:r*i)=e;
  for j=i+1:L
    Hd((j-1)*m+1:m*j,(i-1)*r+1:r*i)=d*ae^(j-1-i)*be;
  end
end

Hs=zeros(L*m,L*m);
for i=1:L
  for j=i+1:L
    Hs((j-1)*m+1:m*j,(i-1)*m+1:m*i)=d*ae^(j-1-i)*c;
  end
end

for i=1:L
  oe((i-1)*m+1:i*m,:)=d*ae^(i-1);
end

K=1; % K=N-L;
for i=1:L
    YL(1+(i-1)*m:i*m,:)=y(i:K+i-1,:)';
    UL(1+(i-1)*r:i*r,:)=u(i:K+i-1,:)';
end

x = pinv(oe)*(YL-Hd*UL-Hs*YL);
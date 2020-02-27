function yp_M=dsr_prd(A,B,D,E,C,F,y,u,x0,M)
% DSR_PRD Computes the M-step ahead prediction of y_t.
% SYNTAX
% Yp_M = dsr_prd(A,B,D,E,CF,F,Y,U,x0,M)
% DESCRIPTION
% Old outputs up to time t-M and all relevant inputs 
% are used to predict the output at time t. 
% ON INPUT
% A B D E CF F x0 - Innovations form model matrices and initial state, x0,
%                   as computed by DSR, i.e., the outputs from DSR. 
% Y, U            - The output and input data matrices, respectively, i.e., 
%                   Y is an (N x m) data matrix and U is an (N x r) data matrix,
%                   where N is the number of observations, m is the number of outputs
%                   and r is the number of input variables.
% M               - The prediction horizon, i.e., an integer M>=1.
% ON OUTPUT
% Yp_M            - Matrix with the M-step ahead predictions, i.e. a (N x m) matrix.
%
% FUNCTIONS CALLED: none
% Copyright 2000, Dr. ing. David Di Ruscio. 

% Written 3/2-00

[N,m]=size(y);
yp_M=zeros(N,m);

K=C*pinv(F);
A_e=A-K*D;
B_e=B-K*E;

x=x0;

for t=1:M-1
   yp=D*x+E*u(t,:)';
   yp_M(t,:)=yp';
   x = A*x + B*u(t,:)';
end
yp=D*x+E*u(M,:)';
yp_M(M,:)=yp';

for t=M:N-1
   x=x0;
   x = A_e*x+B_e*u(t-M+1,:)'+K*y(t-M+1,:)';
   x0=x;  
   for i=1:M-1
      x = A*x + B*u(t-M+i+1,:)';
   end
   yp=D*x+E*u(t+1,:)';
   yp_M(t+1,:)=yp';   
end   
% END DSR_PRD   
      
   

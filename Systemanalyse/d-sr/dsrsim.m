
function [y,x] = dsrsim(a,b,d,e,u,x0);
%DSRSIM Simulation of discrete-time linear systems.
%       Y=DSRSIM(A,B,D,E,U)
%       Y=DSRSIM(A,B,D,E,U,x0)
%    PURPOSE:
%    Simulate the time response of the discrete system:
%
%               x_{t+1} = Ax_t + Bu_t
%               y_t     = Dx_t + Eu_t
%
%    ON INPUT:
%    A,B,D,E - Discrete dynamic model matrices.
%    U         Matrix U must have as many columns as there are inputs, u. 
%              Each row of U corresponds to a new time point, a (N x r) matrix
%              where r is the number of inputs.
%    x0      - Optional initial state vector. x0 is a (n x 1) vector where
%              n is the number of states. 
%    ON OUTPUT:
%    Y       - Matrix Y with system outputs, i.e. a (N x m) matrix where m is
%              the number of output variables.
%              When invoked with left hand arguments,
%              Y = DSRSIM(A,B,D,E,U) or Y=DSRSIM(A,B,D,E,U,x0)
%    X       - Matrix X with system states, i.e. a (N x n) matrix
%              When invoked with left hand arguments,
%              [Y,X] = DSRSIM(A,B,D,E,U) or [Y,X]=DSRSIM(A,B,D,E,U,x0)
%              returns the output and state time history in the 
%              matrices Y and X.
%---------------------------------------------------------------------------

% WRITTEN TO BE USED AS SUPPLEMENT TO THE DSR IDENTIFICATION ALGORITHM
% DATE: 26. august 1996

%1.
[N,r]=size(u); [m,n]=size(d);

%2. Define the initial state vector
if nargin == 5
 xm=zeros(n,1);
else
 xm=x0;
end

%3. Initialize output arrays
y=zeros(N,m); x=zeros(N,n);

%4. Simulation loop for evaluation of the states
for i=1:N
% store xm in array x for the states
  x(i,:)=xm';

% update (integrate) the state vector
  xm = a*xm + b*u(i,:)';

end

%5.  Compute the outputs outside the loop for increased speed
y=x*d'+u*e';
%
% END DSRSIM
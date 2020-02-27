function [W,P,C,T,U,E,F]=pls(X,Y,a);
%PLS The Partial Least Squares (PLS) algorithm.
%    [W,P,C,T,U,E,F]=pls(X,Y,a);
%    PURPOSE:
%    Given data matrices X and Y. This algorithm compute matrices P, C, T,
%    E and F which defines the decomposition X = T P' + E and Y = T C' + F.
%    The PLS estimate of the regression matrix B, in the bilinear equation 
%    Y = X B + E is given by B=W*inv(P'*W)*C', where B is an (r x m) matrix
%    of regression coefficients.
%
%    ON INPUT:
%    X     - An (N x r) data matrix with input variables.
%    Y     - An (N x m) data matrix with output variables.
%    a     - Number of components (factors) to estimate, 1 < a <= r.
%
%    ON OUTPUT:
%    T,P,C - Matrices in the PLS decomposition X = T P' + E and Y = T C' + F.
%    W     - Orthonormal weighting matrix.
%    P     - Matrix with lin. independent columns, loading matrix for X.
%    C     - Weight matrix for Y.
%    T     - Matrix of orthogonal score vectors for X.
%    U     - Matrix of orthogonal score vectors for Y.
%    E, F  - The residuals in the PLS decomposition.
%    REMARK:
%    The matrix of regression coefficients is B=W inv(P' W) C'.
%    See also D-SR functions PLS2 and PLS3.
%-----------------------------------------------------------------------------
%                                       COPYRIGHT 1996, FANTOFT PROCESS
%                                       License belong to:
%                                       Product id: 10 0000
%-----------------------------------------------------------------------------

% Last update: Sun Nov 24 20:23:46 GMT 1996, David Di Ruscio.

if nargin == 2                   
 a=r; a=fchan('Number of components:',a);
end

[Nx,r]=size(X); [Ny,m]=size(Y); N=min(Nx,Ny);
W=zeros(r,a); P=zeros(r,a); C=zeros(m,a); 
T=zeros(N,a); U=zeros(N,a);

for i=1:a                         %%%%% START MAIN ITERATION LOOP %%%%%%%%%%%%
 xy= Y'*X;
 [uw,sw,vw]=svd(xy'); w=uw(:,1);  % Solve eigenvector problem, X'YY'X*w=s*w,
                                  % as the SVD of the matrix X'Y for accuracy.
 t=X*w/(w'*w);                    % Score vector for X.
 c=Y'*t/(t'*t);                   % Weight vector for Y=TC'+F.
 u=Y*c/(c'*c);                    % Score vector for Y.
 p=X'*t/(t'*t);                   % loading vector for X=TP'+E.

 X=X-t*p';                        % Updating X and Y.          
 Y=Y-t*c';                 

 W(:,i)=w; P(:,i)=p; C(:,i)=c;
 T(:,i)=t; U(:,i)=u;
end                               %%%%% END MAIN ITERATION LOOP %%%%%%%%%%%%%%
E=X; F=Y;                         % The residuals.
%
% END D-SR TOOLBOX FUNCTION PLS
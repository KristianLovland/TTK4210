function [b,e]=besolv(Z,U,At,A,D,L,g,r);
% BESOLV Low level DSR function.
%        [B,E]=besolv(Z,U,At,A,D,L,g,r);
%        Optimal least squares metod to solve
%        Z = tilde(B) U for system matrices B and E
%        ALGORITHM
%        Solve equivalent system cs([B;E])=N\cs(Z)

% Version: 20. october 1996
[m,n]=size(D);

if g==1                                      % form the matrix N in LS problem
 E0=[sobsv(A,D,L,1) zeros(L*m,m)];
 EL=[zeros(L*m,n) seye(eye(m,m),L,L)];
else
 E0=[sobsv(A,D,L,1)];
 EL=[zeros(L*m,n)];
end
N=kron( U((L+g-1)*r+1:(L+g)*r,:)',EL);

for i=1:L
  if g==1
   E1=[sobsv(A,D,L,i+1) seye(eye(m,m),L,i)];
  else
   E1=[sobsv(A,D,L,i+1)];
  end
  N=N+kron( U((i-1)*r+1:i*r,:)',E0-At*E1);
  E0=E1;
end

theta=N\Z(:);   % cs([B;E])=pinv(N)*cs(Z)    % solve LS problem for theta

b=zeros(n,r);
e=zeros(m,r);
if g== 1                                     % form B and E from theta
 nt=r*(n+m);
else
 nt=n*r;
end
j=0;
if g==1
 for i=1:r
   b(:,i)=theta(1+j:n+j);
   e(:,i)=theta(1+j+n:m+j+n);
   j=j+n+m;
 end
else
 for i=1:r
   b(:,i)=theta(1+j:n+j);
   j=j+n;
 end
end
%
% END BESOLV
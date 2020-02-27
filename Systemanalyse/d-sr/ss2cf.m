function [ac,bc,dc,t]=ss2cf(a,b,d);
% [ac,bc,dc,t]=ss2cf(a,b,d);
% Purpose
% Transform a state space model (a,b,d) into observable canonical form
% (ac,bc,dc).
% xh = a x + b u     -->          xch = ac xc + bc u
% y  = d x + e u     -->          y   = dc xc + e  u
%
% Note:
% The matrix (e) is not influenced by state coordinate transformation.

[n,m]=size(b);           % check input parameters
o=d;                     % form observability matrix (ill conditioned step !)
for i=1:n-1
o=[o;d*a^i];
end
t=o(1:n,:);              % transformation matrix, ill-cond. if non-singular
%
oc=o*inv(t);             % obsv. matrix for obsv. can. form state space model
dc=d*inv(t);             % can. form (ac,bc,dc)
ac=pinv(oc)*o*a*inv(t);
bc=pinv(oc)*o*b;

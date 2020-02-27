function r=dread(text,f)
% DREAD Low level DSR function.
%       r=dread(text,f)
%       Example:  temp=dread('Temp',temp)
%                 var =dread('Variables',[3 4 6])
%                 temp=dread('Temp')
l=length(text);
if (nargin<2), f=0;, end
n=length(f);

s=[text ' ...'];
for i=l:33-3*n
s=[s '.'];

end
s=[s '. ('];
for i=1:n
s=[s sprintf(' %g ',f(i)) ];
end
s=[s ') =? '];
r=input(s);
if length(r) == 0, r=f; end
%
% END DREAD


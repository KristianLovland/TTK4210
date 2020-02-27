function I=seye(e,L,l);
% SEYE  Low level DSR function.
%       I=seye(e,L,l)
%       Return the extended matrix I_{m x r} with E in row position l

[m,r]=size(e);
I=zeros(L*m,r);
I((l-1)*m+1:l*m,:)=e;
%
% END SEYE
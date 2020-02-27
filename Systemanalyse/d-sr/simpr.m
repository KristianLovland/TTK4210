function h=simpr(a,b,d,e,L,l);
% SIMPR Low level DSR function.
%       Hd = simpr(a,b,d,e,L,l)
%       PURPOSE
%       Return the extended matrix H^d_L(:,1:r) with impulse response 
%       matrices shifted l rows down in the Lm times r array Hd.
%       INTERNAL ROUTINE

[m,n]=size(d); [n,r]=size(b);
h=zeros(L*m,r); h((l-1)*m+1:l*m,:)=e;
for i=l+1:L
  h((i-1)*m+1:i*m,:)=d*a^(i-l-1)*b;
end
%
% END SIMPR
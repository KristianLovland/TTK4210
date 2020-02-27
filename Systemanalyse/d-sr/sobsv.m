function o=sobsv(a,d,L,l);
% SOBSV  Low level DSR function.
%        O = sobsv(a,d,L,l)
%        Return the extended observability matrix O_L 
%        shifted l rows down in the Lm times n array.

[m,n]=size(d);
o=zeros(L*m,n);
for i=l:L
  o((i-1)*m+1:i*m,:)=d*a^(i-l);
end
%
% END SOBSV
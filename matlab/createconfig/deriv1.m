function [res]=deriv1(f,x, dim)
nel=size(f);
nel1=size(x);

f1=reshape(f,[nel1 1]);
num=nel(dim);

% if nel(1)>nel(2)
%     num=nel(1);
% else
%     num=nel(2);
% end

% if (nel ne nel1) then begin
%  print,'Inconsistant input, stop.'
%  stop
% endif

res=zeros(num,1);
for i=3:num-2 
    res(i)=(1./12.0./(x(i+1)-x(i)))*(8.0*f1(i+1)-8.0*f1(i-1)-f1(i+2)+f1(i-2));
end
%for i=1,nel-2 do res(i)=(1.d0/2.d0/(x(i+1)-x(i)))*(f(i+1)-f(i-1))
res(1)=res(3);
res(2)=res(3);
res(num)=res(num-2);
res(num-1)=res(num-2);

end




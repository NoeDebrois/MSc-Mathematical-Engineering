function [dudx] = df(x,u)
%df Numerical derivative of the function u(x) defined on the vector xi
n= length(x);
   xd = diff([x(3),x,x(n-2)]);  
   ud = diff([u(3),u,u(n-2)]);  
   dudx = (ud(1:end-1)./xd(1:end-1).*xd(2:end) ...
          + ud(2:end)./xd(2:end).*xd(1:end-1)) ...
          ./ (xd(2:end)+xd(1:end-1));
end
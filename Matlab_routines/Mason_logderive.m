function [fs,df,ddf] = Mason_logderive(x,f,width)

%   A special purpose routine for finding the first and second logarithmic 
%   derivatives of function f(x) with a low frequency behavior and high 
%   frequency noise.  It returns 'fs'-- a smoother version of 'f' and the 
%   first and second log derivative of fs.

df = zeros(size(f));
ddf = zeros(size(f));
fs = zeros(size(f));
x = repmat(x,1,size(f,2));
lx = log(x);
ly = log(f);

for i=1:size(f,1)
    %Introduce 'Gaussian' centered on lx(i) and width 'width'
    w = exp( -(lx-lx(i,:)).^2 / (2*width.^2) );
    for j = 1:size(f,2)
        % Truncate the Gaussian for faster computation
        ww = find(w(:,j) > 0.03);                          
        res = polyfitw(lx(ww,j),ly(ww,j),w(ww,j),2);
        fs(i,j) = exp(res(3) + res(2)*lx(i,j) + res(1)*(lx(i,j)^2));
        df(i,j) = res(2)+(2*res(1)*lx(i,j));
        ddf(i,j) = 2*res(1);
    end
end

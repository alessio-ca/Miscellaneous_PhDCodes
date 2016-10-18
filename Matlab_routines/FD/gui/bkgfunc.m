function y = bkgfunc(a,bkg,x,y)
y=  sum(sum((bkg-(a(1)+a(2)*x+a(3)*y+a(4)*(x.*x)+a(5)*(y.*y)+a(6)*(x.*y))).^2));

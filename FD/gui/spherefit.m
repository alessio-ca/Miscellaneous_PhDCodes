function [x,y,z,r]=spherefit(xc,yc,zc)
%usage: [x,y,z,r]=spherefit(xc,yc,zc)
% Returns the coordinates and radius of a sphere fitted to the supplied
% coordinates
% IH 24/11/02

TolX = 1.e-8; TolFun = 1.e-8; MaxFunEvals = 1.e+8;
Options = optimset('TolX',TolX,'TolFun',TolFun,'MaxFunEvals',MaxFunEvals);
startx=mean(xc);
starty=mean(yc);
startz=mean(zc);
startr=2;
astart=[startx starty startz startr];
afinish=fminsearch('sphfunc',astart,Options,xc,yc,zc);

x=afinish(1);
y=afinish(2);
z=afinish(3);
r=afinish(4);


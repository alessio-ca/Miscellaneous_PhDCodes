pdem=createpde();

%% Geometry (N.B! Units are scaled in um)
r1 = [3 4 -2 2 2 -2  -2 -2 2 2]';
c1 = [1 0 0 .4]';
% Append extra zeros to the circles
% so they have the same number of rows as the rectangle
c1=[c1;zeros(length(r1)-length(c1),1)];
gdm = [r1,c1];
g = decsg(gdm, 'r1-c1', char('r1','c1')');
geometryFromEdges(pdem,g);
figure(1) 
pdegplot(pdem, 'EdgeLabels', 'on','SubdomainLabels','on')
xlim([-3,3])
axis equal

%% Problem init (N.B! Inputs of the interp functions must be in SI units, non scaled)
load('Force.mat')
load('divF.mat')

spl_x=csapi({-2e-6:.5e-7:2e-6,-2e-6:.5e-7:2e-6},F.Vx);
spl_y=csapi({-2e-6:.5e-7:2e-6,-2e-6:.5e-7:2e-6},F.Vy);
spl_divF=csapi({-2e-6:.5e-8:2e-6,-2e-6:.5e-8:2e-6},divF);


[phi,r] = meshgrid((2*pi)/100:(2*pi)/100:2*pi,.4*ones(1,100));
[xx,yy] = Transform.Pol2Car(phi,r);


fun_xp = @(region,state) 1e+12*fnval(spl_x,{-2e-6,1e-6*region.y}); %Edge 3 (force expressed in pN)
fun_xm = @(region,state) -1e+12*fnval(spl_x,{2e-6,1e-6*region.y}); %Edge 1 (force expressed in pN)
fun_yp = @(region,state) 1e+12*(fnval(spl_y,{1e-6*region.x,-2e-6}))'; %Edge 4 (force expressed in pN)
fun_ym = @(region,state) -1e+12*(fnval(spl_y,{1e-6*region.x,2e-6}))'; %Edge 2 (force expressed in pN)
fun_r = @(region,state) -1e+12*(cos(atan2(1e-6*region.y,1e-6*region.x)).*fnval(spl_x,[1e-6*region.x;1e-6*region.y]) + sin(atan2(1e-6*region.y,1e-6*region.x)).*fnval(spl_y,[1e-6*region.x;1e-6*region.y]));

fun_divF = @(p,t,u,time) 1e+6*fnval(spl_divF,[1e-6*(p(1,t(1,:))+p(1,t(2,:))+p(1,t(3,:)))/3;1e-6*(p(2,t(1,:))+p(2,t(2,:))+p(2,t(3,:)))/3]);

applyBoundaryCondition(pdem,'Edge',1,'g',fun_xm,'Vectorized','on'); %Neumann, n . (c \nabla u) + qu = g on edge
applyBoundaryCondition(pdem,'Edge',3,'g',fun_xp,'Vectorized','on'); %Neumann, n . (c \nabla u) + qu = g on edge
applyBoundaryCondition(pdem,'Edge',2,'g',fun_ym,'Vectorized','on'); %Neumann, n . (c \nabla u) + qu = g on edge
applyBoundaryCondition(pdem,'Edge',4,'g',fun_yp,'Vectorized','on'); %Neumann, n . (c \nabla u) + qu = g on edge
applyBoundaryCondition(pdem,'Edge',5:8,'u',0,'Vectorized','on'); 
applyBoundaryCondition(pdem,'Edge',5:8,'g',fun_r,'Vectorized','on'); %Neumann, n . (c \nabla u) + qu = g on edge
c = 1;
a = 0;
f = fun_divF;

%% Init mesh
[p,e,t] = initmesh(g,'Hmax',.1);
fprintf('Mesh generated. Number of nodes: %d\n',size(p,2));
figure(2)
pdemesh(p,e,t);
axis equal

%%
u = assempde(pdem,p,e,t,c,a,f);
figure(3)
pdesurf(p,t,u);

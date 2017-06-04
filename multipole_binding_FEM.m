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
load('Efield.mat')

spl_x=csapi({-2e-6:.5e-7:2e-6,-2e-6:.5e-7:2e-6},F.Vx');
spl_y=csapi({-2e-6:.5e-7:2e-6,-2e-6:.5e-7:2e-6},F.Vy');
spl_divF=csapi({-2e-6:.5e-8:2e-6,-2e-6:.5e-8:2e-6},divF');
spl_E=csapi({-2e-6:.5e-7:2e-6,-2e-6:.5e-7:2e-6},(norm(Et).^2)');

fun_li = @(t) (1e+12*fnval(spl_x,[.4e-6*cos(t);.4e-6*sin(t)]).*sin(t)*0.4 - ...
               1e+12*fnval(spl_y,[.4e-6*cos(t);.4e-6*sin(t)]).*cos(t)*0.4);

fun_xp = @(region,state) 1e+12*fnval(spl_x,{-2e-6,1e-6*region.y}); %Edge 3 (force expressed in pN)
fun_xm = @(region,state) -1e+12*fnval(spl_x,{2e-6,1e-6*region.y}); %Edge 1 (force expressed in pN)
fun_yp = @(region,state) 1e+12*(fnval(spl_y,{1e-6*region.x,-2e-6}))'; %Edge 4 (force expressed in pN)
fun_ym = @(region,state) -1e+12*(fnval(spl_y,{1e-6*region.x,2e-6}))'; %Edge 2 (force expressed in pN)

fun_r = @(region,state) -1e+12*(cos(atan2(1e-6*region.y,1e-6*region.x)).*fnval(spl_x,[1e-6*region.x;1e-6*region.y]) + sin(atan2(1e-6*region.y,1e-6*region.x)).*fnval(spl_y,[1e-6*region.x;1e-6*region.y])); %Edge 5:8 (force expressed in pN)
fun_E = @(region,state) 1e-12*(fnval(spl_E,[1e-6*region.x;1e-6*region.y])); %(E expressed in V/um)
fun_lib = @(region,state) integral(fun_li,0,atan2(region.y,region.x));

fun_divF = @(p,t,u,time) 1e+6*fnval(spl_divF,[1e-6*(p(1,t(1,:))+p(1,t(2,:))+p(1,t(3,:)))/3;1e-6*(p(2,t(1,:))+p(2,t(2,:))+p(2,t(3,:)))/3]);

applyBoundaryCondition(pdem,'Edge',1,'g',fun_xm,'Vectorized','on'); %Neumann, n . (c \nabla u) + qu = g on edge
applyBoundaryCondition(pdem,'Edge',3,'g',fun_xp,'Vectorized','on'); %Neumann, n . (c \nabla u) + qu = g on edge
applyBoundaryCondition(pdem,'Edge',2,'g',fun_ym,'Vectorized','on'); %Neumann, n . (c \nabla u) + qu = g on edge
applyBoundaryCondition(pdem,'Edge',4,'g',fun_yp,'Vectorized','on'); %Neumann, n . (c \nabla u) + qu = g on edge
applyBoundaryCondition(pdem,'Edge',5:8,'u',fun_lib,'Vectorized','off'); 
applyBoundaryCondition(pdem,'Edge',5:8,'g',fun_r,'Vectorized','on'); %Neumann, n . (c \nabla u) + qu = g on edge
c = 1;
a = 0;
f = fun_divF;

%% Init mesh
generateMesh(pdem,'Hmax',.05);
[p,e,t] = meshToPet(pdem.Mesh);
fprintf('Mesh generated. Number of nodes: %d\n',size(p,2));
figure(2)
subplot(1,2,1)
pdemesh(p,e,t);
title('Mesh')
xlim([-3,3])
axis equal

%%
u = assempde(pdem,c,a,f);
subplot(1,2,2)
pdesurf(p,t,u);
xlabel('x (um)')
ylabel('y (um)')
zlabel('U (10^{-18} J)')
title('Dipole binding potential')

%% Check on solution
[ux,uy] = pdegrad(p,t,u);
ugrad = [ux;uy];
ugrad = pdeprtni(p,t,ugrad);
ugrad_n = sqrt(ugrad(:,1).^2 + ugrad(:,2).^2);
figure(3)
subplot(1,2,1)
pdesurf(p,t,ugrad_n);
title('Calculated Force')
xlabel('x (um)')
ylabel('y (um)')
subplot(1,2,2)
surf(norm(F),'EdgeColor','none')
title('Original force')
xlabel('x (um)')
ylabel('y (um)')

%%
result = createPDEResults(pdem,u);
[x,y]=meshgrid(-2:.005:2,-2:.005:2);
uintrp = interpolateSolution(result,x,y);
uintrp = reshape(uintrp,size(x));
spl_u=csapi({-2:.005:2,-2:.005:2},uintrp');
fun_u = @(X) fnval(spl_u,X);
%%
options=optimset('PlotFcns',@optimplotfval);
minPt=fminsearch(fun_u,[0;-1.1],options);
minPt=[minPt,fminsearch(fun_u,[0;1.1],options)];
minPt=[minPt,fminsearch(fun_u,[-1.1;0],options)];
minPt=[minPt,fminsearch(fun_u,[1.1;0],options)];

    
%%
figure(4)
surf(-2:.005:2,-2:.005:2,fun_u({-2:.005:2,-2:.005:2})','edgealpha',0)
axis equal
view(2)
hold on
plot3(minPt(1,:),minPt(2,:),10*ones(1,length(minPt)),'ro')
hold off







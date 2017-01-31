%% DEFINITION OF BEAM

w0 = 5e-3;
Ex0 = 1;
Ey0 = 0;
R = 20e-3;
Nphi = 16;
Nr = 10;
bg = BeamGauss(Ex0,Ey0,w0,R,Nphi,Nr);


f = 1.075*R;
ef = EFieldFocus(bg,f);
%% DEFINITION OF INDUCEDDIPOLES
ep = 2.25;
a = 200e-9;
lambda0 = 633e-9;
alpharc = InducedDipole.polarizability('corrected',a,ep,'lambda0',lambda0);
x1 = 0;
y1 = -1e-6;
z1 = 0;
r1=Point(x1,y1,z1);

x2=0;
y2= +1e-6;
z2=0;
r2=Point(x2,y2,z2);

id1 = InducedDipole(alpharc,'lambda0',lambda0,'rd',r1);
id2 = InducedDipole(alpharc,'lambda0',lambda0,'rd',r2);

%% DIPOLE 1&2, ELECTRIC FIELDS
[y,z,x] = meshgrid(-3e-6:1e-8:3e-6,-4e-6:1e-8:4e-6,0);
r = Point(x,y,z);
E=ef.E(r);
I = 0.5/PhysConst.Z0*sqrt(ef.er/ef.mr)*(E.*conj(E));

E1 = id1.E(r,E);
E2 = id2.E(r,E);
E1.Vx(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;
E1.Vy(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;
E1.Vz(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;
E2.Vx(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;
E2.Vy(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;
E2.Vz(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;


figure(1)
g=subplot(2,2,1);
surf(y*1e+9,z*1e+9,norm(real(E1+E2)),'edgealpha',0)
title('Perpendicular, zy plane')
axis equal
view(2)
xlabel('y [nm]')
ylabel('z [nm]')
color_scaling=g.CLim;
%%
E1 = id1.E(r,E.zrotation(pi/2));
E2 = id2.E(r,E.zrotation(pi/2));
E1.Vx(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;
E1.Vy(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;
E1.Vz(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;
E2.Vx(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;
E2.Vy(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;
E2.Vz(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;


subplot(2,2,3)
surf(y*1e+9,z*1e+9,norm(real(E1+E2)),'edgealpha',0)
title('Parallel, zy plane')
axis equal
view(2)
xlabel('y [nm]')
ylabel('z [nm]')
caxis(color_scaling)

[x,y,z] = meshgrid(-4e-6:1e-8:4e-6,-4e-6:1e-8:4e-6,0);
r = Point(x,y,z);
E=ef.E(r);
E1 = id1.E(r,E);
E2 = id2.E(r,E);
E1.Vx(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;
E1.Vy(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;
E1.Vz(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;
E2.Vx(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;
E2.Vy(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;
E2.Vz(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;

subplot(2,2,2)
surf(x*1e+9,y*1e+9,norm(real(E1+E2)),'edgealpha',0)
title('Perpendicular, xy plane')
axis equal
view(2)
xlabel('x [nm]')
ylabel('y [nm]')
caxis(color_scaling)

E1 = id1.E(r,E.zrotation(pi/2));
E2 = id2.E(r,E.zrotation(pi/2));
E1.Vx(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;
E1.Vy(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;
E1.Vz(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;
E2.Vx(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;
E2.Vy(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;
E2.Vz(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;

subplot(2,2,4)
surf(x*1e+9,y*1e+9,norm(real(E1+E2)),'edgealpha',0)
title('Parallel, xy plane')
axis equal
view(2)
xlabel('x [nm]')
ylabel('y [nm]')
caxis(color_scaling)

drawnow()


%% Plot Field + Beam
[y,z,x] = meshgrid(-3e-6:1e-8:3e-6,-4e-6:1e-8:4e-6,0);
r = Point(x,y,z);
E=ef.E(r);
I = 0.5/PhysConst.Z0*sqrt(ef.er/ef.mr)*(E.*conj(E));

E1 = id1.E(r,E);
E2 = id2.E(r,E);
E1.Vx(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;
E1.Vy(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;
E1.Vz(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;
E2.Vx(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;
E2.Vy(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;
E2.Vz(x.^2+(y-y2).^2+z.^2<(0.4e-6).^2 | x.^2+(y+y2).^2+z.^2<(0.4e-6).^2)=NaN;

figure(2)

cmapX=([ones(100,1), [1:-0.01:0.01]', [1:-0.01:0.01]']); %red colormap
cmapY=parula(64); %standard colormap
cmap=[parula(15);cmapY];
%colormap extrema
cmin = 0;
cmax = 1;
m1=15;
m2=64;
% CData for pcolor(1)
C1 = min(m1,round((m1-1)*I./max(max(I))+1)); 
% CData for pcolor(2)
C2 = m1+1+min((m2-1),round(m2*norm(real(E1+E2))./max(max(norm(real(E1+E2))))));

13
h(1)=pcolor(y*1e+6,z*1e+6,I./max(max(I)));
shading interp
set(h(1),'FaceColor','flat','CData',C1)
hold on

h(2)=pcolor(y*1e+6,z*1e+6,norm(real(E1+E2))./max(max(norm(real(E1+E2)))));
shading interp
set(h(2),'FaceAlpha',0.75,'CData',C2)

hold off

colormap(cmap)
axis equal
axis tight
xlabel('y [\mum]')
ylabel('z [\mum]')

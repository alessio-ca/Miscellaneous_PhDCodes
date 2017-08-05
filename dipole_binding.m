%% INITIALIZATION
ep = 2.25;
a = .1e-6;
lambda0 = 800e-9;
alpharc = InducedDipole.polarizability('corrected',a,ep,'lambda0',lambda0);
Iel=1e12;
El=sqrt(Iel/(PhysConst.c0*PhysConst.e0));
Ei = ComplexVector(0,0,0,El/2,0,0);
E_inc1 = EFieldPlaneWave(El/2,Vector(0,0,0,0,0,1),Vector(0,0,0,1,0,0));
E_inc2 = EFieldPlaneWave(0,Vector(0,0,0,0,0,-1),Vector(0,0,0,1,0,0));

x2=0;
y2=1e-6;
z2=0;
r2=Point(x2,y2,z2);

id1 = InducedDipole(alpharc,'lambda0',lambda0);
id2 = InducedDipole(alpharc,'lambda0',lambda0,'rd',r2);


%% DIPOLE 1&2, ELECTRIC FIELDS
id1 = InducedDipole(alpharc,'lambda0',lambda0);
[x,y,z] = meshgrid(-4e-6:1e-8:4e-6,-4e-6:1e-8:4e-6,0);
r = Point(x,y,z);
E1 = id1.E(r,Ei);
E2 = id2.E(r,Ei);
B2 = id2.B(r,Ei);
%%


E1.Vx(x.^2+(y-y2).^2+z.^2<(2*a).^2 | x.^2+(y+y2).^2+z.^2<(2*a).^2)=NaN;
E1.Vy(x.^2+(y-y2).^2+z.^2<(2*a).^2 | x.^2+(y+y2).^2+z.^2<(2*a).^2)=NaN;
E1.Vz(x.^2+(y-y2).^2+z.^2<(2*a).^2 | x.^2+(y+y2).^2+z.^2<(2*a).^2)=NaN;
E2.Vx(x.^2+(y-y2).^2+z.^2<(2*a).^2 | x.^2+(y+y2).^2+z.^2<(2*a).^2)=NaN;
E2.Vy(x.^2+(y-y2).^2+z.^2<(2*a).^2 | x.^2+(y+y2).^2+z.^2<(2*a).^2)=NaN;
E2.Vz(x.^2+(y-y2).^2+z.^2<(2*a).^2 | x.^2+(y+y2).^2+z.^2<(2*a).^2)=NaN;


figure(1)
g=subplot(2,2,1);
surf(y*1e+9,z*1e+9,norm(real(E1+E2)),'edgealpha',0)
title('Perpendicular, zy plane')
axis equal
view(2)
xlabel('y [nm]')
ylabel('z [nm]')
color_scaling=g.CLim;

E1 = id1.E(r,Ei.zrotation(pi/2));
E2 = id2.E(r,Ei.zrotation(pi/2));
E1.Vx(x.^2+(y-y2).^2+z.^2<(2*a).^2 | x.^2+(y+y2).^2+z.^2<(2*a).^2)=NaN;
E1.Vy(x.^2+(y-y2).^2+z.^2<(2*a).^2 | x.^2+(y+y2).^2+z.^2<(2*a).^2)=NaN;
E1.Vz(x.^2+(y-y2).^2+z.^2<(2*a).^2 | x.^2+(y+y2).^2+z.^2<(2*a).^2)=NaN;
E2.Vx(x.^2+(y-y2).^2+z.^2<(2*a).^2 | x.^2+(y+y2).^2+z.^2<(2*a).^2)=NaN;
E2.Vy(x.^2+(y-y2).^2+z.^2<(2*a).^2 | x.^2+(y+y2).^2+z.^2<(2*a).^2)=NaN;
E2.Vz(x.^2+(y-y2).^2+z.^2<(2*a).^2 | x.^2+(y+y2).^2+z.^2<(2*a).^2)=NaN;

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
E1 = id1.E(r,Ei);
E2 = id2.E(r,Ei);
E1.Vx(x.^2+(y-y2).^2+z.^2<(2*a).^2 | x.^2+(y+y2).^2+z.^2<(2*a).^2)=NaN;
E1.Vy(x.^2+(y-y2).^2+z.^2<(2*a).^2 | x.^2+(y+y2).^2+z.^2<(2*a).^2)=NaN;
E1.Vz(x.^2+(y-y2).^2+z.^2<(2*a).^2 | x.^2+(y+y2).^2+z.^2<(2*a).^2)=NaN;
E2.Vx(x.^2+(y-y2).^2+z.^2<(2*a).^2 | x.^2+(y+y2).^2+z.^2<(2*a).^2)=NaN;
E2.Vy(x.^2+(y-y2).^2+z.^2<(2*a).^2 | x.^2+(y+y2).^2+z.^2<(2*a).^2)=NaN;
E2.Vz(x.^2+(y-y2).^2+z.^2<(2*a).^2 | x.^2+(y+y2).^2+z.^2<(2*a).^2)=NaN;;

subplot(2,2,2)
surf(x*1e+9,y*1e+9,norm(real(E1+E2)),'edgealpha',0)
title('Perpendicular, xy plane')
axis equal
view(2)
xlabel('x [nm]')
ylabel('y [nm]')
caxis(color_scaling)

E1 = id1.E(r,Ei.zrotation(pi/2));
E2 = id2.E(r,Ei.zrotation(pi/2));
E1.Vx(x.^2+(y-y2).^2+z.^2<(2*a).^2 | x.^2+(y+y2).^2+z.^2<(2*a).^2)=NaN;
E1.Vy(x.^2+(y-y2).^2+z.^2<(2*a).^2 | x.^2+(y+y2).^2+z.^2<(2*a).^2)=NaN;
E1.Vz(x.^2+(y-y2).^2+z.^2<(2*a).^2 | x.^2+(y+y2).^2+z.^2<(2*a).^2)=NaN;
E2.Vx(x.^2+(y-y2).^2+z.^2<(2*a).^2 | x.^2+(y+y2).^2+z.^2<(2*a).^2)=NaN;
E2.Vy(x.^2+(y-y2).^2+z.^2<(2*a).^2 | x.^2+(y+y2).^2+z.^2<(2*a).^2)=NaN;
E2.Vz(x.^2+(y-y2).^2+z.^2<(2*a).^2 | x.^2+(y+y2).^2+z.^2<(2*a).^2)=NaN;

subplot(2,2,4)
surf(x*1e+9,y*1e+9,norm(real(E1+E2)),'edgealpha',0)
title('Parallel, xy plane')
axis equal
view(2)
xlabel('x [nm]')
ylabel('y [nm]')
caxis(color_scaling)

drawnow()
id1 = InducedDipole(alpharc,'lambda0',lambda0);
%% FORCE BY DIPOLE 1, 1D case
%Parallel polarization
[x,y,z] = meshgrid(2*a:.5e-8:5e-6,0,0);
r = Point(x,y,z);
E1 = id1.E(r,Ei);
B1 = id1.B(r,Ei);
Ei1 = E_inc1.E(r);
Bi1 = E_inc1.B(r);
Ei2 = E_inc2.E(r);
Bi2 = E_inc2.B(r);
[F,Fgrad,Fscat] = id2.force_general(r,E1+Ei1+Ei2,B1+Bi1+Bi2,1);

figure(2)
plot(F.X*1e+6,F.Vx*1e+12)
%Perpendicular polarization
Ei_rot = Ei.zrotation(pi/2);
E_inc1_rot = EFieldPlaneWave(El/2,Vector(0,0,0,0,0,1),Vector(0,0,0,0,1,0));
E_inc2_rot = EFieldPlaneWave(0,Vector(0,0,0,0,0,-1),Vector(0,0,0,0,1,0));
[x,y,z] = meshgrid(2*a:.5e-8:5e-6,0,0);
r = Point(x,y,z);
E1 = id1.E(r,Ei_rot);
B1 = id1.B(r,Ei_rot);
Ei1 = E_inc1_rot.E(r);
Bi1 = E_inc1_rot.B(r);
Ei2 = E_inc2_rot.E(r);
Bi2 = E_inc2_rot.B(r);
[F,~,~] = id2.force_general(r,E1+Ei1+Ei2,B1+Bi1+Bi2,1);

hold on
plot(F.X*1e+6,F.Vx*1e+12)
title('Optical binding force for two colloids (axial)')
view(2)
xlabel('\Deltax [\mum]')
ylabel('Force [pN]')
legend('Parallel','Perpendicular')
xlim = get(gca,'xlim');  %Get x range 
plot([xlim(1) xlim(2)],[0 0],'k')
hold off

drawnow()
%% FORCE BY DIPOLE 1, 2D case
[x,y,z] = meshgrid(-2e-6:.5e-8:2e-6,-2e-6:.5e-8:2e-6,0);
r = Point(x,y,z);
E1 = id1.E(r,Ei);
B1 = id1.B(r,Ei);
Ei1 = E_inc1.E(r);
Bi1 = E_inc1.B(r);
Ei2 = E_inc2.E(r);
Bi2 = E_inc2.B(r);
[F,Fgrad,Fscat] = id2.force_general(r,E1+Ei1+Ei2,B1+Bi1+Bi2,2);
F.Vx(x.^2+y.^2<(2*a).^2)=NaN;
F.Vy(x.^2+y.^2<(2*a).^2)=NaN;
E1.Vx(x.^2+y.^2<(2*a).^2)=NaN;
E1.Vy(x.^2+y.^2<(2*a).^2)=NaN;
E1.Vz(x.^2+y.^2<(2*a).^2)=NaN;
%%
zmin = min(min(min(F.Vy)));
zmax = max(max(max(F.Vy)));
zinc = (zmax - zmin) / 20;
zlevs = zmin:zinc:zmax;
figure(3)
g=subplot(1,2,1);
%surf(x*1e+9,y*1e+9,norm(F),'edgealpha',0)
contourf(x*1e+6,y*1e+6,F.Vy,zlevs);
hold on
h=streamslice(F.X(1:20:end,1:20:end)*1e+6,F.Y(1:20:end,1:20:end)*1e+6,F.Vx(1:20:end,1:20:end),F.Vy(1:20:end,1:20:end),4,'method','cubic');
%h=streamline(F.X(1:20:end,1:20:end)*1e+9,F.Y(1:20:end,1:20:end)*1e+9,F.Vx(1:20:end,1:20:end),F.Vy(1:20:end,1:20:end),F.X(1:30:end,1:30:end)*1e+9,F.Y(1:30:end,1:30:end)*1e+9);
%streakarrow(F.X(1:20:end,1:20:end)*1e+9,F.Y(1:20:end,1:20:end)*1e+9,F.Vx(1:20:end,1:20:end),F.Vy(1:20:end,1:20:end),1.5,1)
%quiver(F.X(1:20:end,1:20:end)*1e+9,F.Y(1:20:end,1:20:end)*1e+9,F.Vx(1:20:end,1:20:end),F.Vy(1:20:end,1:20:end),1)
for i=1:length(h)
    h(i).LineWidth=1;
end
set(h,'Color',[0.6 0.2 0.])
colormap(g,'parula');
hold off
title('F_y and streamlines, xy plane')
axis equal
view(2)
xlabel('x [\mum]')
ylabel('y [\mum]')

magnitude=norm(F);
subplot(1,2,2)
surf(x*1e+6,y*1e+6,magnitude,'edgealpha',0)
colormap jet
title('Magnitude of Force, xy plane')
axis equal
view(2)
xlabel('x [\mum]')
ylabel('y [\mum]')

% magnitude=norm(E1);
% figure(4)
% surf(x*1e+9,y*1e+9,magnitude,'edgealpha',0)
% title('xy plane')
% axis equal
% view(2)
% xlabel('x [nm]')
% ylabel('y [nm]')





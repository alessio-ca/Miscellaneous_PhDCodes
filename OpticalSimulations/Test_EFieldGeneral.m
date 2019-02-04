%% INITIALIZATION
ep = 2.25;
a = .2e-6;
lambda0 = 1064e-9;
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

%% FORCE BY DIPOLE 1, 1D case
%Parallel polarization
[x,y,z] = meshgrid(0.5e-6:.5e-8:5e-6,0,0);
r = Point(x,y,z);
E1 = id1.E(r,Ei);
B1 = id1.B(r,Ei);
Ei1 = E_inc1.E(r);
Bi1 = E_inc1.B(r);
Ei2 = E_inc2.E(r);
Bi2 = E_inc2.B(r);

Etot=E1+Ei1+Ei2;
Btot=B1+Bi1+Bi2;

[F,Fgrad,Fscat] = id2.force_general(r,Etot,Btot,1);

figure(1)
p1=plot(F.X*1e+6,F.Vx*1e+12);
ax=gca;
eTot=EFieldGeneral(E1+Ei1+Ei2,B1+Bi1+Bi2,1);
[Fgen,Fgradgen,Fscatgen] = id2.force(r,eTot);
hold on
plot(Fgen.X*1e+6,-Fgen.Vx*1e+12,'Color',ax.ColorOrder(1,:),'LineStyle','--')
%Check seems ok


%Perpendicular polarization
Ei_rot = Ei.zrotation(pi/2);
E_inc1_rot = EFieldPlaneWave(El/2,Vector(0,0,0,0,0,1),Vector(0,0,0,0,1,0));
E_inc2_rot = EFieldPlaneWave(0,Vector(0,0,0,0,0,-1),Vector(0,0,0,0,1,0));
[x,y,z] = meshgrid(0.5e-6:.5e-8:5e-6,0,0);
r = Point(x,y,z);
E1 = id1.E(r,Ei_rot);
B1 = id1.B(r,Ei_rot);
Ei1 = E_inc1_rot.E(r);
Bi1 = E_inc1_rot.B(r);
Ei2 = E_inc2_rot.E(r);
Bi2 = E_inc2_rot.B(r);
Etot=E1+Ei1+Ei2;
Btot=B1+Bi1+Bi2;
[F,Fgrad,Fscat] = id2.force_general(r,Etot,Btot,1);

p2=plot(F.X*1e+6,F.Vx*1e+12,'Color',ax.ColorOrder(2,:));
eTot=EFieldGeneral(E1+Ei1+Ei2,B1+Bi1+Bi2,1);
[Fgen,Fgradgen,Fscatgen] = id2.force(r,eTot);
plot(Fgen.X*1e+6,-Fgen.Vx*1e+12,'Color',ax.ColorOrder(2,:),'LineStyle','--')
title('Optical binding force for two colloids (axial), check on formulation')
view(2)
xlabel('\Deltax [\mum]')
ylabel('Force [pN]')
legend([p1 p2],'Parallel','Perpendicular')
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
Etot=E1+Ei1+Ei2;
Btot=B1+Bi1+Bi2;
[F,Fgrad,Fscat] = id2.force_general(r,Etot,Btot,2);
eTot=EFieldGeneral(E1+Ei1+Ei2,B1+Bi1+Bi2,2);
[Fgen,Fgradgen,Fscatgen] = id2.force(r,eTot);

F.Vx(x.^2+y.^2<(2*a).^2)=NaN;
F.Vy(x.^2+y.^2<(2*a).^2)=NaN;
Fgen.Vx(x.^2+y.^2<(2*a).^2)=NaN;
Fgen.Vy(x.^2+y.^2<(2*a).^2)=NaN;
%%
zmin = min(min(min(F.Vy)));
zmax = max(max(max(F.Vy)));
zinc = (zmax - zmin) / 20;
zlevs = zmin:zinc:zmax;
figure(3)
g=subplot(1,2,1);
contourf(x*1e+6,y*1e+6,F.Vy,zlevs);
title('F_y, xy plane')
axis equal
view(2)
xlabel('x [\mum]')
ylabel('y [\mum]')

subplot(1,2,2)
contourf(x*1e+6,y*1e+6,Fgen.Vy,zlevs);
title('F_y, EFieldClass, xy plane')
axis equal
view(2)
xlabel('x [\mum]')
ylabel('y [\mum]')

%Check seems ok




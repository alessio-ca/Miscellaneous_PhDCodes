%% DEFINITION OF PARTICLESPHERICAL & RAYS
c = Point(1,0,0);
r = 1;
nm = 1;
np = 1.5;
bead = ParticleSpherical(c,r,nm,np); %Particle in c, radius r and index np (against medium index nm)
close all

figure(1)
title('Initial configuration')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
bead.plot();


v = Vector(-1,-0.5,0,2,0,0);
P = 1;
pol = Vector(0,0,0,0,1,0); 
r = Ray(v,P,pol);
r.plotnopolnorm(bead.sp,'color','b','ShowArrowHead','off');

%% SCATTERING
r_vec = bead.scattering(r);
rr = r_vec(1).r;
rt = r_vec(1).t;
rr.plotnopol('color','r');
rt.plotnopolnorm(bead.sp,'color','b','ShowArrowHead','off');
rr = r_vec(2).r;
rt = r_vec(2).t;
rr.plotnopol('color','r');
rt.plotnopolnorm(bead.sp,'color','b');
%%
F = bead.force(r)*1e+15;
F.plot('Scale',[1 1/sqrt(F.Vx.^2 + F.Vy.^2)],'Color','k','LineWidth',2);



%% DEFINITION OF PARTICLESPHERICAL & RAYS
c = Point(1,0,0);
r = 1;
nm = 1;
np = 1.5;
bead = ParticleSpherical(c,r,nm,np); %Particle in c, radius r and index np (against medium index nm)
close all

figure(1)
title('Initial configuration')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
bead.plot();

n=1;
m=2;
v = Vector([-1 -1],[-2 2],zeros(n,m),repmat(1,n,m),[1.25 -1.25],zeros(n,m));
P = ones(1,2);
pol = Vector(zeros(n,m),zeros(n,m),zeros(n,m),zeros(n,m),ones(n,m),zeros(n,m)); 
r = Ray(v,P,pol);
r.plotnopolnorm(bead.sp,'color','b','ShowArrowHead','off');

%% SCATTERING
r_vec = bead.scattering(r);
rr = r_vec(1).r;
rt = r_vec(1).t;
rr.plotnopol('color','r');
rt.plotnopolnorm(bead.sp,'color','b','ShowArrowHead','off');
rr = r_vec(2).r;
rt = r_vec(2).t;
rr.plotnopol('color','r');
rt.plotnopolnorm(bead.sp,'color','b');
%%
F = bead.force(r)*1e+15;
F1=Vector(F.X(1),F.Y(1),F.Z(1),F.Vx(1),F.Vy(1),F.Vz(1));
F2=Vector(F.X(2),F.Y(2),F.Z(2),F.Vx(2),F.Vy(2),F.Vz(2));
F=plus(F1,F2);
F.plot('Scale',[1 1.25/sqrt(F.Vx.^2 + F.Vy.^2)],'Color','k','LineWidth',2);
ax = gca;



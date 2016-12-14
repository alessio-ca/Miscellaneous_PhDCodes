%% DEFINITION OF BEAMGAUSS

w0 = 5e-3;
Ex0 = 0.1;
Ey0 = 1i;
R = 1e-3;
Nphi = 16;
Nr = 50;
bg = BeamGauss(Ex0,Ey0,w0,R,Nphi,Nr);

res=Ray.beam2focused(bg,1.2e-3);
res.v=2*res.v;
CellRay={res.v.X,res.v.Y,res.v.Z,res.v.Vx,res.v.Vy,res.v.Vz,res.P,res.pol.X,res.pol.Y,res.pol.Z,res.pol.Vx,res.pol.Vy,res.pol.Vz};
selec=@(X) X(abs(res.v.Vx)<1e-6);
CellRay=cellfun(selec,CellRay,'UniformOutput', false);
res=Ray(Vector(CellRay{1},CellRay{2},CellRay{3},CellRay{4},CellRay{5},CellRay{6}),CellRay{7},Vector(CellRay{8},CellRay{9},CellRay{10},CellRay{11},CellRay{12},CellRay{13}));
res=res.yrotation(pi/2);
dP=topoint(-res.v./2);
res = translate(res,dP);
close all
figure(1)
hold on
axis equal
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
res.plotnopol();
view(3)
grid on
%% DEFINITION OF PARTICLESPHERICAL
c = Point(0.3,0.0,0);
r = 0.5;
nm = 1;
np = 1.5;
bead = ParticleSpherical(c,r,nm,np); %Particle in c, radius r and index np (against medium index nm)
bead.plot();
%%
F = bead.force(res)*1e+15;
Fx=sum(F.Vx(~isnan(F.Vx)));
Fy=sum(F.Vy(~isnan(F.Vx)));
Fz=sum(F.Vz(~isnan(F.Vx)));
Ftot=Vector(F.X(1),F.Y(1),F.Z(1),Fx,Fy,Fz);
Ftot.plot('Scale',[1 1./F0],'Color','k','LineWidth',2);
ax = gca;
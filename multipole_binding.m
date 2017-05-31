%% INITIALIZATION
ep = 2.25;
a = .2e-6;
lambda0 = 1064e-9;
alpharc = InducedDipole.polarizability('corrected',a,ep,'lambda0',lambda0);
Iel=1e12;
El=sqrt(Iel/(PhysConst.c0*PhysConst.e0));
Ei = ComplexVector(0,0,0,El/2,0,0);
E_inc = EFieldPlaneWave(El/2,Vector(0,0,0,0,0,1),Vector(0,0,0,1,0,0));

x2=0;
y2=1e-6;
z2=0;

r1=Point(0,0,0);
r2=Point(x2,y2,z2);

rP = [r1,r2];
idP = InducedDipole(alpharc,'lambda0',lambda0,'rd',rP(1));

for i=2:length(rP)
    idP = [idP,InducedDipole(alpharc,'lambda0',lambda0,'rd',rP(i))];
end

numDip=2;


%% FORCE BY DIPOLE 1, 2D case
% Interacting dipole
[x,y,z] = meshgrid(-2e-6:.5e-8:2e-6,-2e-6:.5e-8:2e-6,0);
r = Point(x,y,z);
tol = 1;
cc = 0;
Ei_n=[E_inc.E(rP(1)),E_inc.E(rP(2))];
chk = zeros(1,numDip);

while tol > 0.001
    Ei_o = Ei_n;
    cc = cc + 1;
    disp(['Iteration ',num2str(cc)]);
    
    Ei_n=[E_inc.E(rP(1)),E_inc.E(rP(2))];
    Bi_n=[E_inc.B(rP(1)),E_inc.B(rP(2))];
    for i=1:length(rP)
        for j=1:length(rP)
            if i~=j
                Ei_n(i) = Ei_n(i)+idP(j).E(rP(i),Ei_o(j));
                Bi_n(i) = Bi_n(i)+idP(j).B(rP(i),Ei_o(j));
            end
        end
    end
    
    for i=1:length(rP)
        chk(i) = norm(Ei_n(i) - Ei_o(i))./norm(Ei_o(i));
    end
    
    tol = max([real(chk),imag(chk)]);
end
%%
Et1 = E_inc.E(r) + idP(1).E(r,Ei_n(1));
Bt1 = E_inc.B(r) + idP(1).B(r,Ei_n(1));
Et2 = E_inc.E(r) + idP(2).E(r,Ei_n(2));
Bt2 = E_inc.B(r) + idP(2).B(r,Ei_n(2));

%%

[F,Fgrad,Fscat] = idP(2).force_general(r,Et1,Bt1,2);
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
figure(1)
g=subplot(1,2,1);
contourf(x*1e+6,y*1e+6,F.Vy,zlevs);
hold on
h=streamslice(F.X(1:20:end,1:20:end)*1e+6,F.Y(1:20:end,1:20:end)*1e+6,F.Vx(1:20:end,1:20:end),F.Vy(1:20:end,1:20:end),4,'method','cubic');

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





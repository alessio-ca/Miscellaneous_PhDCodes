%% INITIALIZATION
ep = 2.25;
a = .527e-6/2;
lambda0 = 1064e-9;
alpharc = InducedDipole.polarizability('corrected',a,ep,'lambda0',lambda0);
Iel=1e12;
El=sqrt(Iel/(PhysConst.c0*PhysConst.e0));
Ei = ComplexVector(0,0,0,El/2,0,0);
E_inc = EFieldPlaneWave(El/2,Vector(0,0,0,0,0,1),Vector(0,0,0,1,0,0));



%% 3 INTERACTING DIPOLES: ITERATIVE SOLVER

%Fixed positions
x1=0;
y1=-0.9769e-6/2;
z1=0;
r1=Point(x1,y1,z1);
x2=0;
y2=0.9769e-6/2;
z2=0;
r2=Point(x2,y2,z2);

meshStep = .01e-6;


% Variable position: Coarse mesh (with 3 mesh steps space)
[x,y,z] = meshgrid(-2e-6 - 3*meshStep:meshStep:2e-6 + 3*meshStep,-2e-6 - 3*meshStep:meshStep:2e-6 + 3*meshStep,0);
r = Point(x,y,z);

Ex=zeros(size(r,1),size(r,2));
Ey=zeros(size(r,1),size(r,2));
Ez=zeros(size(r,1),size(r,2));
Bx=zeros(size(r,1),size(r,2));
By=zeros(size(r,1),size(r,2));
Bz=zeros(size(r,1),size(r,2));
Fx=zeros(size(r,1),size(r,2));
Fy=zeros(size(r,1),size(r,2));
Fz=zeros(size(r,1),size(r,2));
rx=r.X(:);
ry=r.Y(:);
rz=r.Z(:);

%%
% %Enable progress bar for parallel pool
try
    parpool;
catch ME
    if ~strcmp(ME.identifier,'parallel:convenience:ConnectionOpen')
        rethrow(ME)
    end
end

targetWorkCount = size(r,1)*size(r,2);
barWidth= int32( 30 );
p =  TimedProgressBar( targetWorkCount, barWidth, ...
    'Computing, wait for ', ', completed ', 'Concluded in ' );

parfor k = 1:targetWorkCount
    r3=Point(rx(k),ry(k),rz(k));
    rP = [r1,r2,r3];
    idP = InducedDipole(alpharc,'lambda0',lambda0,'rd',rP(1));
    
    for i=2:length(rP)
        idP = [idP,InducedDipole(alpharc,'lambda0',lambda0,'rd',rP(i))];
    end
    
    numDip=3;
    
    if (norm(Point(rx(k),ry(k),rz(k)) - r1)<(2*a - 6*meshStep) || ...
             norm(Point(rx(k),ry(k),rz(k)) - r2)<(2*a - 6*meshStep)) %For force calculation, add 6 grid points inside
        Ex(k)=0;
        Ey(k)=0;
        Ez(k)=0;
        Bx(k)=0;
        By(k)=0;
        Bz(k)=0;
        Fx(k)=0;
        Fy(k)=0;
        Fz(k)=0;
    else
        tol = 1;
        cc = 0;
        Ei_n=[E_inc.E(rP(1)),E_inc.E(rP(2)),E_inc.E(rP(3))];
        Bi_n=[E_inc.B(rP(1)),E_inc.B(rP(2)),E_inc.B(rP(3))];
        chk = zeros(1,numDip);
        
        while tol > 0.001 
            Ei_o = Ei_n;
            cc = cc + 1;        
            for i=1:length(rP)
                Ei_n(i) = E_inc.E(rP(i));
                Bi_n(i) = E_inc.B(rP(i));
            end
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
        
        %Field calculation
        Ex(k)=E_inc.E(r3).Vx;
        Ey(k)=E_inc.E(r3).Vy;
        Ez(k)=E_inc.E(r3).Vz;
        Bx(k)=E_inc.B(r3).Vx;
        By(k)=E_inc.B(r3).Vy;
        Bz(k)=E_inc.B(r3).Vz;
        
        for i=1:length(rP)
            if i~=3
                Ex(k)=Ex(k)+idP(i).E(r3,Ei_n(i)).Vx;
                Ey(k)=Ey(k)+idP(i).E(r3,Ei_n(i)).Vy;
                Ez(k)=Ez(k)+idP(i).E(r3,Ei_n(i)).Vz;
                Bx(k)=Bx(k)+idP(i).B(r3,Ei_n(i)).Vx;
                By(k)=By(k)+idP(i).B(r3,Ei_n(i)).Vy;
                Bz(k)=Bz(k)+idP(i).B(r3,Ei_n(i)).Vz;
            end
        end

        %Force calculation, 5 pts interpol.
        [x_temp,y_temp,z_temp] = meshgrid(r3.X-2*meshStep:meshStep:r3.X+2*meshStep,r3.Y-2*meshStep:meshStep:r3.Y+2*meshStep,0);
        r_temp = Point(x_temp,y_temp,z_temp);
        
        
        E_temp = E_inc.E(r_temp);
        B_temp = E_inc.B(r_temp);
        for i=1:length(rP)
            if i~=3
                E_temp=E_temp + idP(i).E(r_temp,Ei_n(i));
                B_temp=B_temp + idP(i).B(r_temp,Ei_n(i));
            end
        end
        
        [F_temp,~,~] = idP(3).force_general(r3,E_temp,B_temp,2);
        
        Fx(k)=F_temp.Vx;
        Fy(k)=F_temp.Vy;
        Fz(k)=F_temp.Vz;
        
        if (norm(Point(rx(k),ry(k),rz(k)) - r1)<(2*a-3*meshStep) || ...
            norm(Point(rx(k),ry(k),rz(k)) - r2)<(2*a-3*meshStep))%For pre-divF force expression, add 3 grid point inside the particle (discard inner 3 pts)
            Ex(k)=0;
            Ey(k)=0;
            Ez(k)=0;
            Bx(k)=0;
            By(k)=0;
            Bz(k)=0;
            Fx(k)=0;
            Fy(k)=0;
            Fz(k)=0;
        end
        
    end
    p.progress;

end
p.stop;

Et = ComplexVector(r.X,r.Y,r.Z,Ex,Ey,Ez);
Bt = ComplexVector(r.X,r.Y,r.Z,Bx,By,Bz);
F = ComplexVector(r.X,r.Y,r.Z,Fx,Fy,Fz);

%% Plotting
figure(1)
subplot(2,2,1)
surf(r.X,r.Y,norm(Et),'edgealpha',0)
subplot(2,2,2)
surf(norm(Bt),'edgealpha',0)
subplot(2,2,3)
surf(norm(F),'edgealpha',0)

%% Optional: export of the calculated fields
save('Efield_3.mat','Et')
save('Force_3.mat','F')

%% Optional: load the calculated fields
load('Force_3.mat')

%% Spline interpolant

spl_x=csapi({-2e-6-3*meshStep:meshStep:2e-6+3*meshStep,-2e-6-3*meshStep:meshStep:2e-6+3*meshStep},F.Vx');
spl_y=csapi({-2e-6-3*meshStep:meshStep:2e-6+3*meshStep,-2e-6-3*meshStep:meshStep:2e-6+3*meshStep},F.Vy');

fun_x = @(X) fnval(spl_x,X);
fun_y = @(X) fnval(spl_y,X);

[xx,yy,zz] = meshgrid(-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,0);
fineFx = fun_x({-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep})';
fineFy = fun_y({-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep})';

fineFy(yy==0)=0;

figure(2)
subplot(1,2,1)
%surf(fineFy,'edgealpha',0)
plot(-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,fineFy(yy==0))
hold on
plot(-2e-6-3*meshStep:meshStep:2e-6+3*meshStep,F.Vy(F.Y==0))
hold off
subplot(1,2,2)
%surf(F.Vy,'edgealpha',0)
plot(-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,1e8*fineFy(:,1821),'x')
hold on
plot(-2e-6-3*meshStep:meshStep:2e-6+3*meshStep,1e8*F.Vy(:,183),'o')
divF = divergence(xx,yy,fineFx,fineFy);
plot(-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,divF(:,1821))
%plot(-2e-6-2*meshStep:meshStep:2e-6+2*meshStep,F.Vx(F.Y==0))
hold off

divF((xx-x1).^2 + (yy-y1).^2 < (2*a-3*meshStep/10)^2 | (xx-x2).^2 + (yy-y2).^2 < (2*a-3*meshStep/10)^2)=0; %Cut at particle diameter with an internal grid point
%% Optional: export of the calculated fields
save('divF_3.mat','divF')

%%
spl_divF=csapi({-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep},divF');
fun_divF = @(X) fnval(spl_divF,X);
figure(1)
subplot(2,2,4)
surf(fun_divF({-2e-6-3*meshStep:meshStep:2e-6+3*meshStep,-2e-6-3*meshStep:meshStep:2e-6+3*meshStep})','edgealpha',0)

%%
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


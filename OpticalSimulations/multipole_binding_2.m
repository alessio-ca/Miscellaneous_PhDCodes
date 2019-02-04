clear all
close all
%% INITIALIZATION
ep = 1.60^2;
a = .527e-6/2;
lambda0 = 1064e-9;
nm = 1.33;
alpharc = InducedDipole.polarizability('corrected',a,ep,'lambda0',lambda0,'em',nm^2);
%alpharc = InducedDipole.polarizability('clausius-mossotti',a,ep,'lambda0',lambda0,'em',nm^2);
Iel=1e12;
El=sqrt(Iel/(PhysConst.c0*PhysConst.e0));
Ei = ComplexVector(0,0,0,El/2,0,0);
E_inc = EFieldPlaneWave(El/2,Vector(0,0,0,0,0,1),Vector(0,0,0,1,0,0),'er',nm^2,'lambda0',lambda0);

Lx = 2000e-9;
Ly = 2000e-9;
Lz = 2000e-9;


%% 2 INTERACTING DIPOLES: ITERATIVE SOLVER

%Fixed position
x1=0;
y1=0;
z1=0;
r1=Point(x1,y1,z1);

meshStep = .01e-6;


% Variable position: Coarse mesh
[x,y,z] = meshgrid(-2e-6 - 3*meshStep:meshStep:2e-6 + 3*meshStep,-2e-6 - 3*meshStep:meshStep:2e-6 + 3*meshStep,0);
%[x,y,z] = meshgrid(-2e-6 - 3*meshStep:meshStep:2e-6 + 3*meshStep,0,0);
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
    r2=Point(rx(k),ry(k),rz(k));
    rP = [r1,r2];
    idP = InducedDipole(alpharc,'lambda0',lambda0,'rd',rP(1),'er',nm^2);
    
    for i=2:length(rP)
        idP = [idP,InducedDipole(alpharc,'lambda0',lambda0,'rd',rP(i),'er',nm^2)];
    end
    
    numDip=2;
    
    if norm(Point(rx(k),ry(k),rz(k)) - r1)<(2*a - 6*meshStep) %For force calculation, add 6 grid points inside
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
        Ei_n=[E_inc.E(rP(1)),E_inc.E(rP(2))];
        chk = zeros(1,numDip);
        
        while tol > 1e-3 
            Ei_o = Ei_n;
            cc = cc + 1;
            
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
        
        %Field calculation
        Ex(k)=E_inc.E(r2).Vx + idP(1).E(r2,Ei_n(1)).Vx;
        Ey(k)=E_inc.E(r2).Vy + idP(1).E(r2,Ei_n(1)).Vy;
        Ez(k)=E_inc.E(r2).Vz + idP(1).E(r2,Ei_n(1)).Vz;
        Bx(k)=E_inc.B(r2).Vx + idP(1).B(r2,Ei_n(1)).Vx;
        By(k)=E_inc.B(r2).Vy + idP(1).B(r2,Ei_n(1)).Vy;
        Bz(k)=E_inc.B(r2).Vz + idP(1).B(r2,Ei_n(1)).Vz;
        
        %Force calculation, 5 pts interpol.
        [x_temp,y_temp,z_temp] = meshgrid(r2.X-2*meshStep:meshStep:r2.X+2*meshStep,r2.Y-2*meshStep:meshStep:r2.Y+2*meshStep,0);
        %[x_temp,y_temp,z_temp] = meshgrid(r2.X-2*meshStep:meshStep:r2.X+2*meshStep,0,0);
        r_temp = Point(x_temp,y_temp,z_temp);


        E_temp = E_inc.E(r_temp) + idP(1).E(r_temp,Ei_n(1));
        B_temp = E_inc.B(r_temp) + idP(1).B(r_temp,Ei_n(1));
        [F_temp,~,~] = idP(2).force_general(r2,E_temp,B_temp,2);
        %[F_temp,~,~] = idP(2).force_general(r2,E_temp,B_temp,1);
        
        Fx(k)=F_temp.Vx;
        Fy(k)=F_temp.Vy;
        Fz(k)=F_temp.Vz;
        
        if norm(Point(rx(k),ry(k),rz(k)) - r1)<(2*a-3*meshStep) %For pre-divF force expression, add 3 grid point inside the particle
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
%%
clear i
k_wave = (2*pi*nm/lambda0);
Gxx = @(x) exp(i*k_wave*abs(x))./(4*pi*PhysConst.e0*nm^2*abs(x).^3).*(-2*i*k_wave*abs(x) + 2);
Gyy = @(x) exp(i*k_wave*abs(x))./(4*pi*PhysConst.e0*nm^2*abs(x).^3).*(k_wave^2*abs(x).^2 + i*k_wave*abs(x) - 1);
EyInc = E_inc.E(Point(0,0,0)).Vy;
EyDip = @(x) alpharc*Gyy(x).*EyInc; 
EyTot = @(x) (EyInc + EyDip(x))./(1 - alpharc^2 * Gyy(x).^2);

Fx_perp = @(x) abs(alpharc).^2 * (El/2)^2 ./ ((8*pi*PhysConst.e0*nm^2*abs(x).^4)) .* ...
    (3*cos(k_wave*abs(x)) + 3*k_wave.*abs(x).*sin(k_wave*abs(x)) - ...
    2*k_wave^2*abs(x).^2.*cos(k_wave*abs(x)) - k_wave^3*abs(x).^3.*sin(k_wave*abs(x)));
Fx_paral = @(x) abs(alpharc).^2 * (El/2)^2 ./ ((4*pi*PhysConst.e0*nm^2*abs(x).^4)) .* ...
    (-3*cos(k_wave*abs(x)) - 3*k_wave.*abs(x).*sin(k_wave*abs(x)) + k_wave^2*abs(x).^2.*cos(k_wave*abs(x)));
%%
plot(x(x>2*a),Fx_paral(x(x>2*a)))
hold on
plot(F.X(F.X>2*a),F.Vx(F.X>2*a))
hold off
%% Plotting
figure(1)
subplot(2,2,1)
surf(r.X,r.Y,norm(Et))
subplot(2,2,2)
surf(norm(Bt))
subplot(2,2,3)
surf(norm(F))

%% Optional: export of the calculated fields
save('Efield.mat','Et')
save('Force.mat','F')

%% Optional: load the calculated fields
load('Force.mat')

%% Spline interpolant

F.Vx((F.X-x1).^2 + (F.Y-y1).^2 < (2*a-3*meshStep)^2)=NaN;
F.Vy((F.X-x1).^2 + (F.Y-y1).^2 < (2*a-3*meshStep)^2)=NaN;

idxgood=~(isnan(F.Vx) | isnan(F.Vy));
FxG = griddata(x(idxgood),y(idxgood),F.Vx(idxgood),x,y,'natural');
FyG = griddata(x(idxgood),y(idxgood),F.Vy(idxgood),x,y,'natural');

spl_x=csapi({-2e-6-3*meshStep:meshStep:2e-6+3*meshStep,-2e-6-3*meshStep:meshStep:2e-6+3*meshStep},FxG');
spl_y=csapi({-2e-6-3*meshStep:meshStep:2e-6+3*meshStep,-2e-6-3*meshStep:meshStep:2e-6+3*meshStep},FyG');

fun_x = @(X) fnval(spl_x,X);
fun_y = @(X) fnval(spl_y,X);

[xx,yy,zz] = meshgrid(-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,0);
fineFx = fun_x({-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep})';
fineFy = fun_y({-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep})';
divF = divergence(xx,yy,fineFx,fineFy);
spl_divF=csapi({-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep},divF');
fun_divF = @(X) fnval(spl_divF,X);
%% Testing fineFx & fineFy on y=0 and x=-0.1979 (intersection point)
figure(2)
subplot(2,2,1)
plot(-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,fineFx(yy==0))
hold on
plot(-2e-6-3*meshStep:meshStep:2e-6+3*meshStep,F.Vx(F.Y==0),'o')
hold off

subplot(2,2,2)
plot(-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,fineFx(xx==0))
hold on
plot(-2e-6-3*meshStep:meshStep:2e-6+3*meshStep,F.Vx(F.X==0),'o')
hold off

subplot(2,2,3)
plot(-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,fineFy(yy==0))
hold on
plot(-2e-6-3*meshStep:meshStep:2e-6+3*meshStep,F.Vy(F.Y==0),'o')
hold off

subplot(2,2,4)
plot(-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,fineFy(xx==0))
hold on
plot(-2e-6-3*meshStep:meshStep:2e-6+3*meshStep,F.Vy(F.X==0),'o')
hold off
%% Divergence calculation & fitting
divF = divergence(xx,yy,fineFx,fineFy);
divF(xx.^2 + yy.^2 < (2*a - 3*0.05e-7)^2)=0; %Cut at particle diameter with three internal fine grid points
spl_divF=csapi({-2.1e-6:.5e-8:2.1e-6,-2.1e-6:.5e-8:2.1e-6},divF');
fun_divF = @(X) fnval(spl_divF,X);

%% Optional: export of the calculated fields
save('ForceGx.mat','FxG')
save('ForceGy.mat','FyG')
save('divF.mat','divF')
%% Test on divergence smoothness

figure(1)
subplot(2,2,4)
surf(fun_divF({-2e-6:.5e-7:2e-6,-2e-6:.5e-7:2e-6})')

%Test on line x=2e-7
figure(3) 
plot(-2.1e-6:.5e-8:2.1e-6,divF(:,381))
hold on
plot(-2.1e-6:.5e-9:2.1e-6,fun_divF([2e-7*ones(1,length(-2.1e-6:.5e-9:2.1e-6));-2.1e-6:.5e-9:2.1e-6]))
hold off
%%
zmin = min(min(min(F.Vy)));
zmax = max(max(max(F.Vy)));
zinc = (zmax - zmin) / 20;
zlevs = zmin:zinc:zmax;
figure(3)
g=subplot(1,2,1);
%surf(x*1e+9,y*1e+9,norm(F),'edgealpha',0)
%contourf(x*1e+6,y*1e+6,F.Vy,zlevs);
%hold on
h=streamslice(F.X(1:20:end,1:20:end)*1e+6,F.Y(1:20:end,1:20:end)*1e+6,F.Vx(1:20:end,1:20:end),F.Vy(1:20:end,1:20:end),4,'method','cubic');
for i=1:length(h)
    h(i).LineWidth=1;
end
set(h,'Color',[0.6 0.2 0.])
colormap(g,'parula');
%hold off
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
%%
[xprobe,yprobe] = meshgrid(-Lx+3*meshStep:10*meshStep:Lx-3*meshStep,-Ly+3*meshStep:10*meshStep:Ly-3*meshStep);
points=[xprobe(:),yprobe(:)];
points=points(xprobe.^2 + yprobe.^2 > (2*a)^2,:);
points_out=zeros(size(points));
traj_x = zeros(200,size(points,1));
traj_y = zeros(200,size(points,1));

for i = 1:length(points)
    x_in = [points(i,1),points(i,2)];
    [x_out,traj] = graddesc(x_in,fun_x,fun_y,[x1,y1],2*a,[Lx,Ly],[-Lx,-Ly],'PointOut',1);
    points_out(i,:) = x_out;
    traj_x(:,i)=traj(:,1);
    traj_y(:,i)=traj(:,2);
end
scatter(points_out(:,1),points_out(:,2))





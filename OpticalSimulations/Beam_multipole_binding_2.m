%% INITIALIZATION
clear all;
close all;
clc;

%Domain
Lx = 2000e-9;
Ly = 2000e-9;
Lz = 2000e-9;

L = 32;

ep = 1.50^2;
a = .527e-6/2;
lambda0 = 1064e-9;
nm = 1.33;
alpharc = InducedDipole.polarizability('corrected',a,ep,'lambda0',lambda0,'em',nm^2);

%% Setup in POM 0.26
% Beam waist = 10.7e-3 m
% f = 200 mm (tube length for infinity corrected) / 60 (obj magnification)
% Pmax = 60/100 (obj transmission) * 64/100 (AOD transmission) * 2 W

% Beam & objective

w0 = 10.7e-3;
Ex0 = 1;
Ey0 = 0;
NA = 1.20;
f = 200e-3/60;
R = f*NA/nm;
f = nm*R/NA;

% Focused field
Nphi = 16;
Nr = 10;
bg = BeamGauss(Ex0,Ey0,w0,R,Nphi,Nr,'lambda0',lambda0,'er',nm^2);
bg = bg.normalize(60*64*2/10000);
E_inc = EFieldFocus(bg,f,'lambda0',lambda0,'er',nm^2);

% Field
meshStep = .01e-6;
[x,y,z] = meshgrid(-Lx - 3*meshStep:meshStep:Lx + 3*meshStep,-Ly - 3*meshStep:meshStep:Ly + 3*meshStep,0);
ef = E_inc.E(Point(x,y,z),'er',nm.^2);
I = .5*(PhysConst.c0/nm)*(nm^2*PhysConst.e0)*norm(ef).^2;

figure(1)
surf(1e6*x,1e6*y,I)
shading flat
axis auto
xlabel('x [um]')
ylabel('y [um]')
zlabel('E [V/m]')



%% 2 INTERACTING DIPOLES: ITERATIVE SOLVER

%Fixed positions
x1=0;
y1=0;
z1=0;
r1=Point(x1,y1,z1);

phase_r1 = focus_phase(bg,f,x1,y1,z1,'er',nm^2);


% Variable position: Coarse mesh (with 3 mesh steps space)
r = Point(x,y,z);

Ex=zeros(size(r));
Ey=zeros(size(r));
Ez=zeros(size(r));
Bx=zeros(size(r));
By=zeros(size(r));
Bz=zeros(size(r));
Fx=zeros(size(r));
Fy=zeros(size(r));
Fz=zeros(size(r));

rx=r.X(:);
ry=r.Y(:);
rz=r.Z(:);

%%
%Enable progress bar for parallel pool
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
    phase_r2 = focus_phase(bg,f,r2.X,r2.Y,r2.Z,'er',nm^2);

    rP = [r1,r2];
    idP = InducedDipole(alpharc,'lambda0',lambda0,'er',nm^2,'rd',rP(1));
    
    for i=2:length(rP)
        idP = [idP,InducedDipole(alpharc,'lambda0',lambda0,'er',nm^2,'rd',rP(i))];
    end
    
    numDip=2;
    
    if (norm(Point(rx(k),ry(k),rz(k)) - r1)<(2*a - 6*meshStep)) 
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
        E_focus_1 = focus_fast_E(bg,f,phase_r1,'er',nm.^2);
        E_focus_2 = focus_fast_E(bg,f,phase_r2,'er',nm.^2);
        Ei_n_base = [E_focus_1(1,:);E_focus_2(1,:)];
        Bi_n_base=[focus_fast_B(E_inc,E_focus_1);focus_fast_B(E_inc,E_focus_2)];
        
        Ei_n=Ei_n_base;
        Bi_n=Bi_n_base;
        chk = zeros(1,numDip);

        
        while tol > 0.001 && cc < 1000
            Ei_o = Ei_n;
            cc = cc + 1;        
            for i=1:length(rP)
                Ei_n(i,:) = Ei_n_base(i,:);
                Bi_n(i,:) = Bi_n_base(i,:);
            end
            for i=1:length(rP)
                for j=1:length(rP)
                    if i~=j
                        eTemp = focus_fast_E_dip(idP(j),rP(i),Ei_o(j,:));
                        bTemp = focus_fast_B(E_inc,eTemp);
                        Ei_n(i,:) = Ei_n(i,:)+[eTemp(1,1),eTemp(1,2),eTemp(1,3)];
                        Bi_n(i,:) = Bi_n(i,:)+[bTemp(1),bTemp(2),bTemp(3)];
                    end
                end
            end
            
            for i=1:length(rP)
                chk(i) = norm(Ei_n(i,:) - Ei_o(i,:))./norm(Ei_o(i,:));
            end
            
            tol = max([real(chk),imag(chk)]);
        end
        
%         if cc == 1000
%             warning('Max iteration reached for r=(%d,%d,%d)',rx(k),ry(k),rz(k))
%         end
        %Field calculation
        Ex(k)=Ei_n_base(2,1);
        Ey(k)=Ei_n_base(2,2);
        Ez(k)=Ei_n_base(2,3);
        Bx(k)=Bi_n_base(2,1);
        By(k)=Bi_n_base(2,2);
        Bz(k)=Bi_n_base(2,3);
        
        for i=1:length(rP)-1
                eTemp = focus_fast_E_dip(idP(i),r2,Ei_n(i,:));
                bTemp = focus_fast_B(E_inc,eTemp);
                Ex(k) =  Ex(k)+eTemp(1,1);
                Ey(k) =  Ey(k)+eTemp(1,2);
                Ez(k) =  Ez(k)+eTemp(1,3);
                Bx(k) =  Bx(k)+bTemp(1);
                By(k) =  By(k)+bTemp(2);
                Bz(k) =  Bz(k)+bTemp(3);
        end

        %Force calculation, 5 pts interpol.
        [x_temp,y_temp,z_temp] = meshgrid(r2.X-2*meshStep:meshStep:r2.X+2*meshStep,r2.Y-2*meshStep:meshStep:r2.Y+2*meshStep,0);
        r_temp = Point(x_temp,y_temp,z_temp);

        E_temp = zeros(numel(r_temp),3);
        B_temp = zeros(numel(r_temp),3);
        for i=1:numel(r_temp)
            phaseTemp = focus_phase(bg,f,r_temp.X(i),r_temp.Y(i),r_temp.Z(i),'er',nm^2);
            efTemp = focus_fast_E(bg,f,phaseTemp,'er',nm.^2);
            E_temp(i,:) = efTemp(1,:);
            B_temp(i,:) = focus_fast_B(E_inc,efTemp);
        end
        E_temp = ComplexVector(r_temp.X,r_temp.Y,r_temp.Z,...
            reshape(E_temp(:,1),size(r_temp,1),size(r_temp,2)),...
            reshape(E_temp(:,2),size(r_temp,1),size(r_temp,2)),...
            reshape(E_temp(:,3),size(r_temp,1),size(r_temp,2)) ...
        );
        B_temp = ComplexVector(r_temp.X,r_temp.Y,r_temp.Z,...
            reshape(B_temp(:,1),size(r_temp,1),size(r_temp,2)),...
            reshape(B_temp(:,2),size(r_temp,1),size(r_temp,2)),...
            reshape(B_temp(:,3),size(r_temp,1),size(r_temp,2)) ...
        );
    
        for i=1:length(rP)-1
            eTemp = ComplexVector(rP(i).X,rP(i).Y,rP(i).Z,Ei_n(i,1),Ei_n(i,2),Ei_n(i,3));
            E_temp=E_temp + idP(i).E(r_temp,eTemp);
            B_temp=B_temp + idP(i).B(r_temp,eTemp);
        end

        [F_temp,~,~] = idP(2).force_general(r2,E_temp,B_temp,2);
        
        Fx(k)=F_temp.Vx;
        Fy(k)=F_temp.Vy;
        Fz(k)=F_temp.Vz;
        
        if (norm(Point(rx(k),ry(k),rz(k)) - r1)<(2*a-3*meshStep))
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
subplot(1,2,1)
surf(r.X,r.Y,norm(Et),'edgealpha',0)
subplot(1,2,2)
surf(norm(F),'edgealpha',0)

%% Optional: export of the calculated fields
save('Force_2_beam_reduced.mat','F')

%% Optional: load the calculated fields
load('Force_2_beam_reduced.mat')

%% Spline interpolant

F.Vx((F.X-x1).^2 + (F.Y-y1).^2 < (2*a-3*meshStep)^2)=NaN;
F.Vy((F.X-x1).^2 + (F.Y-y1).^2 < (2*a-3*meshStep)^2)=NaN;

idxgood=~(isnan(F.Vx) | isnan(F.Vy));
FxG = griddata(x(idxgood),y(idxgood),F.Vx(idxgood),x,y,'natural');
FyG = griddata(x(idxgood),y(idxgood),F.Vy(idxgood),x,y,'natural');

spl_x=csapi({-Lx-3*meshStep:meshStep:Lx+3*meshStep,-Ly-3*meshStep:meshStep:Ly+3*meshStep},FxG');
spl_y=csapi({-Lx-3*meshStep:meshStep:Lx+3*meshStep,-Ly-3*meshStep:meshStep:Ly+3*meshStep},FyG');

fun_x = @(X) fnval(spl_x,X);
fun_y = @(X) fnval(spl_y,X);

[xx,yy,zz] = meshgrid(-Lx-3*meshStep:meshStep/10:Lx+3*meshStep,-Ly-3*meshStep:meshStep/10:Ly+3*meshStep,0);
fineFx = fun_x({-Lx-3*meshStep:meshStep/10:Lx+3*meshStep,-Ly-3*meshStep:meshStep/10:Ly+3*meshStep})';
fineFy = fun_y({-Lx-3*meshStep:meshStep/10:Lx+3*meshStep,-Ly-3*meshStep:meshStep/10:Ly+3*meshStep})';
divF = divergence(xx,yy,fineFx,fineFy);
spl_divF=csapi({-Lx-3*meshStep:meshStep/10:Lx+3*meshStep,-Ly-3*meshStep:meshStep/10:Ly+3*meshStep},divF');
fun_divF = @(X) fnval(spl_divF,X);
%% Optional: export of the calculated fields
save('ForceGx_2_beam_reduced.mat','FxG')
save('ForceGy_2_beam_reduced.mat','FyG')
save('divF_2_beam_reduced.mat','divF')
%% Testing fineFx & fineFy on y=0 and x=0 
figure(2)
subplot(2,2,1)
plot(-Lx-3*meshStep:meshStep/10:Lx+3*meshStep,fineFx(yy==0))
hold on
plot(-Lx-3*meshStep:meshStep:Lx+3*meshStep,F.Vx(F.Y==0),'o')
hold off
title('Fx, y=0')

subplot(2,2,2)
plot(-Ly-3*meshStep:meshStep/10:Ly+3*meshStep,fineFx(xx==0))
hold on
plot(-Ly-3*meshStep:meshStep:Ly+3*meshStep,F.Vx(F.X==0),'o')
hold off
title('Fy, y=0')


subplot(2,2,3)
plot(-Lx-3*meshStep:meshStep/10:Lx+3*meshStep,fineFy(yy==0))
hold on
plot(-Lx-3*meshStep:meshStep:Lx+3*meshStep,F.Vy(F.Y==0),'o')
hold off
title('Fx, x=0')


subplot(2,2,4)
plot(-Ly-3*meshStep:meshStep/10:Ly+3*meshStep,fineFy(xx==0))
hold on
plot(-Ly-3*meshStep:meshStep:Ly+3*meshStep,F.Vy(F.X==0),'o')
hold off
title('Fy, x=0')


%% Testing divF on y=0 and x=-0.1979 (intersection point)
divF((xx-x1).^2 + (yy-y1).^2 < (2*a)^2)=NaN; %Cut at particle diameter with an internal grid point
edge=-Lx-3*meshStep:meshStep/10:Lx+3*meshStep;
subplot(1,2,1)
plot(-Lx-3*meshStep:meshStep/10:Lx+3*meshStep,divF(yy==0))
hold on
plot(edge,1e-4+fun_divF([edge;zeros(1,length(edge))]))
hold off
subplot(1,2,2)
plot(-Lx-3*meshStep:meshStep/10:Lx+3*meshStep,divF(xx==0))
hold on
plot(edge,1e-4+fun_divF([zeros(1,length(edge));edge]))
hold off

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
%%
[xprobe,yprobe] = meshgrid(-Lx+3*meshStep:20*meshStep:Lx-3*meshStep,-Ly+3*meshStep:20*meshStep:Ly-3*meshStep);
points=[xprobe(:),yprobe(:)];
points=points(xprobe.^2 + yprobe.^2 > (2*a)^2,:);
points_out=zeros(size(points,1),3);
traj_x = zeros(200,size(points,1));
traj_y = zeros(200,size(points,1));

parfor i = 1:length(points)
    x_in = [points(i,1),points(i,2)];
    [x_out,traj] = graddesc(x_in,fun_x,fun_y,[x1,y1],2*a,[Lx,Ly],[-Lx,-Ly],'PointOut',1);
    points_out(i,:) = x_out;
    traj_x(:,i)=traj(:,1);
    traj_y(:,i)=traj(:,2);
end
points_out=points_out(points_out(:,1)<Lx & points_out(:,1)>-Lx & points_out(:,2)<Ly & points_out(:,2)>-Ly,:);
points_out=points_out(points_out(:,1).^2 + points_out(:,2).^2 > 4*a^2,:);
figure(1)
scatter(points_out(:,1),points_out(:,2))
% hold on
% for ii=1:length(points)
%     sx = traj_x(:,ii);
%     sy = traj_y(:,ii);
%     plot(sx,sy,'r')
%     sx = traj_x(1:5:end,ii);
%     sy = traj_y(1:5:end,ii);
%     ux = gradient(sx);
%     uy = gradient(sy);
%     quiver(sx,sy,ux,uy,'b');
% end
% hold off

%% Load the energy field from other source
load('U_2_beam_reduced.mat')
%%
[xxx,yyy,zzz] = meshgrid(-Lx:meshStep:Lx,-Ly:meshStep:Ly,0);
zmin = min(min(min(uintrp)));
zmax = max(max(max(uintrp)));
zinc = (zmax - zmin) / 30;
zlevs = zmin:zinc:zmax;
figure(3)
g=subplot(1,2,1);
contourf(xxx*1e+6,yyy*1e+6,uintrp,zlevs);
%surf(xxx*1e+6,yyy*1e+6,uintrp,'EdgeColor','none');
hold on
% h=streamslice(F.X(1:1:end,1:1:end)*1e+6,F.Y(1:1:end,1:1:end)*1e+6,F.Vx(1:1:end,1:1:end),F.Vy(1:1:end,1:1:end),4,'method','cubic');
% for i=1:length(h)
%     h(i).LineWidth=1;
% end
% set(h,'Color',[0.6 0.2 0.])
colormap(g,'parula');
hold off
title('Energy Map & Force streamlines')
axis equal
view(2)
xlabel('x [\mum]')
ylabel('y [\mum]')

magnitude=norm(F);
subplot(1,2,2)
surf(x*1e+6,y*1e+6,magnitude,'edgealpha',0)
colormap jet
title('Force field (magnitude)')
axis equal
view(2)
xlabel('x [\mum]')
ylabel('y [\mum]')

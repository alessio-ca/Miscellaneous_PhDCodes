load('Force.mat')
load('divF.mat')
load('Efield.mat')
F.Vx(F.X.^2 + F.Y.^2 < (0.527e-6)^2)=0;
spl_x=csapi({-2e-6:.5e-7:2e-6,-2e-6:.5e-7:2e-6},F.Vx');
spl_y=csapi({-2e-6:.5e-7:2e-6,-2e-6:.5e-7:2e-6},F.Vy');
spl_divF=csapi({-2e-6:.5e-8:2e-6,-2e-6:.5e-8:2e-6},divF');
spl_E=csapi({-2e-6:.5e-7:2e-6,-2e-6:.5e-7:2e-6},(norm(Et).^2)');


fun_xp = @(x,y) 1e+12*fnval(spl_x,[x;y]); %Edge 3 (force expressed in pN)
fun_xm = @(x,y) -1e+12*fnval(spl_x,[x;y]); %Edge 1 (force expressed in pN)
fun_yp = @(x,y) 1e+12*fnval(spl_y,[x;y]); %Edge 4 (force expressed in pN)
fun_ym = @(x,y) -1e+12*fnval(spl_y,[x;y]); %Edge 2 (force expressed in pN)

y=-2e-6:.5e-7:2e-6;
x=2e-6*ones(1,length(y));
subplot(1,2,1)
plot(y,F.Vx(F.X==2e-6))
subplot(1,2,2)
plot(y,fun_xp(x,y))
%%
fun_xp = @(X) 1e+12*fnval(spl_x,X); %Edge 3 (force expressed in pN)
fun_yp = @(X) 1e+12*fnval(spl_y,X); %Edge 3 (force expressed in pN)

[x,y]=meshgrid(-2e-6:.5e-7:2e-6,-2e-6:.5e-7:2e-6);
subplot(1,2,1)
surf(x,y,F.Vx)
subplot(1,2,2)
surf(x,y,fun_xp({-2e-6:.5e-7:2e-6,-2e-6:.5e-7:2e-6})')

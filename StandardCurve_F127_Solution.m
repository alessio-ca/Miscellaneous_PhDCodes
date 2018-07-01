filename = 'F:\AC_2018_01_23\StandardCurveF127_H2O.txt';
series = dlmread(filename,'\t',4,0);
V1mL=1.012;
m_f_Salt = [0;series(:,1)]./([0;series(:,1)] + [V1mL;series(:,2)]);
density = [V1mL;series(:,3)]/V1mL;
x = [100*m_f_Salt,(100*m_f_Salt).^2]; %Constrain the first point
b = x\(density-1);
m_f_Salt_f = linspace(min(m_f_Salt),max(m_f_Salt))';
xf = [100*m_f_Salt_f,(100*m_f_Salt_f).^2]; %Constrain the first point
density_calc = xf*b + 1; %Add f(0) = 1

figure(1)
plot(100*m_f_Salt,density,'o')
hold on
plot(100*m_f_Salt_f,density_calc)
hold off

m_Salt_curve = linspace(0,100*36/(100+36))';
density_curve = [m_Salt_curve,m_Salt_curve.^2]*b + 1; %Constrain the first point & Add f(0) = 1
figure(2)
plot(m_Salt_curve,m_Salt_curve.*density_curve)
hold on
plot(m_Salt_curve,100*(1./(1-m_Salt_curve/100) - 1))
hold off

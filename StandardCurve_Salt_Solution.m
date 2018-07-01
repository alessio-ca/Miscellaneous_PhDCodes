filename = 'F:\AC_2017_12_21\StandardCurveNaCl_H2O.txt';
series = dlmread(filename,'\t',4,0);
V1mL=1.012;
m_f_Salt = [0;series(:,1)]./([0;series(:,1)] + [V1mL;series(:,2)]);
density = [V1mL;series(:,3)]/V1mL;
x = [ones(length(m_f_Salt),1),100*m_f_Salt];
b = x\density;
density_calc = x*b;

figure(1)
plot(100*m_f_Salt,density,'o')
hold on
plot(100*m_f_Salt,density_calc)
hold off

m_Salt_curve = linspace(0,100*36/(100+36))';
density_curve = [ones(length(m_Salt_curve),1),m_Salt_curve]*b;
figure(2)
plot(m_Salt_curve,m_Salt_curve.*density_curve)
hold on
plot(m_Salt_curve,100*(1./(1-m_Salt_curve/100) - 1))
hold off

%% Test on single exponential
dt = 1e-6;
t = logspace(log10(dt),0,200)';
tau_rel = 1e-3;
ACF_exp = exp(-t./tau_rel);
[tau,MSD,ACF,MSD_vec,ACF_vec]=DLS_Analysis_Javi(t,ACF_exp);

close all
semilogx(t,ACF_exp,'*',tau,ACF)

%% Test on double exponential
dt = 1e-6;
t = logspace(log10(dt),0,200)';
tau_rel_1 = 1e-3;
tau_rel_2 = 1e-1;
ACF_exp = 0.5*(exp(-t./tau_rel_1) + exp(-t./tau_rel_2));
semilogx(t,ACF_exp)
[tau,MSD,ACF,MSD_vec,ACF_vec]=DLS_Analysis_Javi(t,ACF_exp);

close all
semilogx(t,ACF_exp,'*',tau,ACF)
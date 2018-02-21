t = Time;
data = [x1,x2,x3];
[D,I] = DLS_DistroFit(t,data);
[Z,PDI]=DLS_CumFit(t,data);
%%
%K2 parameter increase
[~,~,~,~,Beta_vec]=DLS_CumFit(t,data);
epsilon = 0.9;
data_new = Beta_vec.^(1-epsilon).*data.^epsilon;
[D_new,I_new] = DLS_DistroFit(t,data_new);
[Z_new,PDI_new]=DLS_CumFit(t,data_new);
%%
%Interpolation
I_spline_old=csape(D*1e3',I');
I_interp_old = fnval(Size',I_spline_old)';
I_spline=csape(D_new*1e3',I_new');
I_interp = fnval(Size',I_spline)';
%%
semilogx(Size,I_Malv,Size,I_interp_old./sum(I_interp_old)*1e2,Size,I_interp./sum(I_interp)*1e2)
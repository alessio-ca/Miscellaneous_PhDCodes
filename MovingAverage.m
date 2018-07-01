function [acf_MA,tau_MA] = MovingAverage(fs,f_c,data)

ind = 1:length(data);
time = (1/fs)*ind';

%Window size
WinSize=round(fs/f_c);


% Centers of bins
ind_min=ind(ind-floor(WinSize)>0);
ind_min=ind_min(1);

ind_max=ind(ind+ceil(WinSize)-length(time)<0);
ind_max=ind_max(end);



binscenters = [1,(ind_min:WinSize:ind_max),length(time)];

acf_MA = movmean(data,WinSize);
tau_MA = movmean(time,WinSize);

acf_MA = acf_MA(binscenters,:);
tau_MA = tau_MA(binscenters);
end

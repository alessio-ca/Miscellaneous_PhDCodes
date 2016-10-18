folder1='/Users/AleCac/Documents/Python_stuff/';
filename=[folder1,'RheoData.csv'];
data=dlmread(filename,',');
[fs,spectrum]=Alessio_FFT(data(:,2),1000);
semilogy(fs,spectrum)
hold on
[fs,spectrum]=Alessio_FFT(data(:,3),1000);
semilogy(fs,spectrum)
hold off

data_noise(:,1)=data(:,1);
data_noise(:,2)=data(:,2)+5e-2*sin(2*pi*100*data(:,1)/1000);
plot(data_noise(:,1),data_noise(:,2))
hold on
plot(data(:,1),data(:,2))
hold off
[fs,spectrum]=Alessio_FFT(data_noise(:,2),1000);
semilogy(fs,spectrum)

data_noise(:,3)=data(:,3)+5e-3*sin(2*pi*100*data(:,1)/1000);
data_noise(:,3)=data_noise(:,3)+5e-2*sin(2*pi*120*data(:,1)/1000);
[fs,spectrum]=Alessio_FFT(data_noise(:,3),1000);
semilogy(fs,spectrum)
dlmwrite([folder1,'RheoData_noise.csv'],data_noise, 'delimiter', ',', 'precision', 9);

%%
Q=1000;
fps=1000;
W0=100/(fps/2);
BW=W0/Q;
[b,a] = iirnotch(W0,BW);
data_filtered(:,2)=filter(b,a,data_noise(:,2));
data_filtered(:,3)=filter(b,a,data_noise(:,3));
W0=120/(fps/2);
BW=W0/Q;
[b,a] = iirnotch(W0,BW);
[fs,spectrum]=Alessio_FFT(data_filtered(:,2),1000);
semilogy(fs,spectrum)
%%
data_filtered(:,3)=filter(b,a,data_filtered(:,3));
[fs,spectrum]=Alessio_FFT(data_filtered(:,3),1000);
semilogy(fs,spectrum)
data_filtered(:,1)=data_noise(:,1);
dlmwrite([folder1,'RheoData_filtered.csv'],data_filtered, 'delimiter', ',', 'precision', 9);
close all

%Folder, file name and acquisition fps
folder='/Users/AleCac/Documents/Python_stuff/';
filename=[folder,'RheoData.csv'];
fps=1000;

%Pure data, visualize FFT
data=dlmread(filename,',');
figure(1)
[fs,spectrum]=Alessio_FFT(data(:,2),fps);
semilogy(fs,spectrum)
hold on
[fs,spectrum]=Alessio_FFT(data(:,3),fps);
semilogy(fs,spectrum)
hold off
xlabel('Frequency (Hz)')
ylabel('|P(f)|')

%%
%Noisy data, construct and visualize FFT
%X traj: one peak
noise_f = 100;
data_noise(:,1)=data(:,1);
data_noise(:,2)=data(:,2)+5e-2*sin(2*pi*noise_f*data(:,1)/fps);
figure(1)
plot(data_noise(1:1000,1),data_noise(1:1000,2))
hold on
plot(data(1:1000,1),data(1:1000,2))
hold off
legend('X, Noisy','X, Pure')
xlabel('Time (ts)')
ylabel('Coordinate (um)')
figure(2)
[fs,spectrum]=Alessio_FFT(data_noise(:,2),1000);
semilogy(fs,spectrum)
xlabel('Frequency (Hz)')
ylabel('|P(f)|')
title('X, noisy spectrum')
%%
%Y traj: two peaks
noise_f1=100;
noise_f2=120;
data_noise(:,3)=data(:,3)+5e-3*sin(2*pi*noise_f1*data(:,1)/1000);
data_noise(:,3)=data_noise(:,3)+5e-2*sin(2*pi*noise_f2*data(:,1)/1000);
[fs,spectrum]=Alessio_FFT(data_noise(:,3),1000);
semilogy(fs,spectrum)
xlabel('Frequency (Hz)')
ylabel('|P(f)|')
title('Y, noisy spectrum')
%dlmwrite([folder,'RheoData_noise.csv'],data_noise, 'delimiter', ',', 'precision', 9);

%%
%Filter design
Q=1000; %High quality factor notch
notch_freq1=100;
notch_freq2=120;
W0=notch_freq1/(fps/2);
BW=W0/Q;
[b,a] = iirnotch(W0,BW);
%Filter once for the X (one peak)
data_filtered(:,2)=filter(b,a,data_noise(:,2));
data_filtered(:,3)=filter(b,a,data_noise(:,3));

W0=notch_freq2/(fps/2);
BW=W0/Q;
[b,a] = iirnotch(W0,BW);
%Filter twice for the Y (two peaks)
data_filtered(:,3)=filter(b,a,data_filtered(:,3));

%%
%Visualize
[fs,spectrum]=Alessio_FFT(data_filtered(:,2),1000);
semilogy(fs,spectrum)
hold on
[fs,spectrum]=Alessio_FFT(data_filtered(:,3),1000);
semilogy(fs,spectrum)
hold off
xlabel('Frequency (Hz)')
ylabel('|P(f)|')
legend('X, corrected','Y, corrected')
title('Corrected spectra')
%%
%Print 
data_filtered(:,1)=data_noise(:,1);
%dlmwrite([folder,'RheoData_filtered.csv'],data_filtered, 'delimiter', ',', 'precision', 9);
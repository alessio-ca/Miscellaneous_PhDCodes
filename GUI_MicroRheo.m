%Input: enter mother directory name (it must end with a "\" sign)
folder1='C:\Users\ac2014\Documents\MATLAB\Git-MATLAB\Data_Microrheology\';

%Input: enter the last stamp of the track_output files (a list if more than
%one is present)

cell_folder_list={'60x_1.20_Si02_1000FPS_1.5um_Tweezer_G_0.05.13Aug2016_15.27.31'};
%1.2 is found to be best for log-sampling param
%cell_folder_list={'Simulation_1000FPS_Kel_1e-6_D_2e-6'};
%1.25 is found to be best for log-sampling param

listlen=length(cell_folder_list);

%Input: enter frame rate
fps=1000;

%Input: FFT calculation
FFT=1;

%Input: microrheo param
Dbead=1.5e-6;
%Dbead=2e-6;

kB = 1.3806503 * 1E-23;
Temp = 298;


close all
%%
for frame=1:listlen
    filename=['track_particle_output',char(cell_folder_list(frame)),'.csv'];
    try traj=dlmread([folder1,'De-drifted_',filename],',');
        disp('De-drifted trajectory file found.')
    catch ME
        if strcmp(ME.identifier,'MATLAB:dlmread:FileNotOpened')
            disp('De-drifted trajectory file not found. Proceeding with un-drifted...')
            try traj=dlmread([folder1,'Linked_',filename],',');
                traj(:,4)=[];
                traj=circshift(traj,[0 1]);
                disp('Linked trajectory file found.')
            catch ME
                if strcmp(ME.identifier,'MATLAB:dlmread:FileNotOpened')
                    error('No linked trajectory file found! Please run the DataAnalysis_Standard routine.')
                else
                    rethrow(ME)
                end
            end
            const=traj(1,2:3);
            traj(:,2:3)=traj(:,2:3)-repmat(const,length(traj),1);
        else
            rethrow(ME)
        end
    end
    rheofactor=1e12*kB*Temp./[mean(traj(:,2).^2),mean(traj(:,3).^2)];
    rheofactor=rheofactor./(3*pi*Dbead);
    [acfx,~,~,~] = acf(1/fps,traj(:,2));
    [acfy,tau,~,~] = acf(1/fps,traj(:,3));
    
    idx=find(acfx<0, 1 );
    idy=find(acfy<0, 1 );
    xx=smooth(acfx(idx:end),length(acfx)/100);
    yy=smooth(acfy(idy:end),length(acfy)/100);
    acfxfiltered=[acfx(1:idx-1);xx];
    acfyfiltered=[acfy(1:idy-1);yy];

    

    
    maxlag=max(length(tau)-50000,50);
    maxlag=min(floor(length(tau)/4),maxlag);
    %Log sampling
    interval=1;
    res=1;
    while (interval < maxlag/4)
        res=res+1;
        interval=ceil(1.2^res);
    end
    res=res-1;
    interval=ceil(1.2.^(0:res));
    
    plotinterval=unique(floor(logspace(0,log10(maxlag),1000)));
    %Cut-off (if present)
    taucutoff=2;
    interval=interval(tau(interval)<taucutoff);
    tau_downsample=tau(interval);
    acfx_downsample=acfxfiltered(interval);
    acfy_downsample=acfyfiltered(interval);
    acfx_spline=csape(tau_downsample,acfx_downsample);
    acfy_spline=csape(tau_downsample,acfy_downsample);
    
    %Plotting section
    figure(1)
    ax1=subplot(2,2,1);
    semilogx(tau(plotinterval),acfx(plotinterval),'o')
    hold on
    semilogx(tau(plotinterval),acfxfiltered(plotinterval),'o')
    semilogx(tau_downsample,acfx_downsample,'s','MarkerEdgeColor','k','MarkerFaceColor',[0 0 1],'MarkerSize',10)
    semilogx(tau_downsample,fnval(tau_downsample,acfx_spline),'--k');
    xlabel('Time (s)')
    ylabel('A(t)')
    hold off
    
    ax2=subplot(2,2,3);
    semilogx(tau(plotinterval,:),acfy(plotinterval,:),'o');
    hold on
    semilogx(tau(plotinterval),acfyfiltered(plotinterval),'o')
    semilogx(tau_downsample,acfy_downsample,'s','MarkerEdgeColor','k','MarkerFaceColor',[0 0 1],'MarkerSize',10)
    semilogx(tau_downsample,fnval(tau_downsample,acfy_spline),'--k');
    xlabel('Time (s)')
    ylabel('A(t)')
    hold off
    
    linkaxes([ax1,ax2],'x')

    

    %%
    if FFT==0
        try acf_fft=dlmread([folder1,'FFT-ACF_',filename],',');
            disp('FFT-ACF file found.')
        catch ME
            if strcmp(ME.identifier,'MATLAB:dlmread:FileNotOpened')
                error('FFT-ACF file not found. Re-run with the FFT calculation enabled if you wish to continue with the microrheology analysis.');
            else
                rethrow(ME)
            end
        end
    elseif FFT==1
        
        Beta=5000;
        resfreq=fps/length(tau);
        minfreq=1/(resfreq*tau_downsample(end));
        maxfreq=fps/(2*resfreq);
        omega=resfreq*logspace(log10(minfreq),log10(maxfreq));
        %oversample_tau=(1/(Beta*fps):1/(Beta*fps):tau_downsample(end));
        oversample_tau=unique(logspace(log10(1/(Beta*fps)),log10(tau_downsample(end)),Beta*fps*tau_downsample(end)));
       %%

        %Enable progress bar for parallel pool
        try
            parpool;
        catch ME
            if ~strcmp(ME.identifier,'parallel:convenience:ConnectionOpen')
                rethrow(ME)
            end
        end
        warning('off','MATLAB:Java:DuplicateClass')
        pctRunOnAll javaaddpath java
        progressStep=ceil(length(omega)/100);
        
        %X analysis
        disp('Analysis of X direction in progress...')
        oversample_acf=fnval(oversample_tau,acfx_spline);
        
        fft_acf = 1i*2*pi*omega;
        fft_acf = fft_acf + (1-exp(-1i*2*pi*omega*oversample_tau(1)))*(oversample_acf(1)-1)/oversample_tau(1);
        ppm = ParforProgMon('Progress: ', length(omega), progressStep, 300, 80);
        
        parfor i=1:length(omega)
            fft_acf(i)=fft_acf(i)-sum(diff(exp(-1i*2*pi*omega(i)*oversample_tau)).*diff(oversample_acf)./diff(oversample_tau));
            if mod(i,progressStep)==0
                ppm.increment();
            end
        end
        ppm.delete()
        clear oversample_acf
        acf_fft = [omega',real(fft_acf)',imag(fft_acf)'];
        acf_fft(:,2) = acf_fft(:,2)./(-4*pi*pi*acf_fft(:,1).^2);
        acf_fft(:,3) = acf_fft(:,3)./(-4*pi*pi*acf_fft(:,1).^2);
        clear fft_acf

        %Y analysis
        disp('Analysis of Y direction in progress...')
        oversample_acf=fnval(oversample_tau,acfy_spline);
        
        fft_acf = 1i*omega*2*pi;
        fft_acf = fft_acf + (1-exp(-1i*2*pi*omega*oversample_tau(1)))*(oversample_acf(1)-1)/oversample_tau(1);
        ppm = ParforProgMon('Progress: ', length(omega), progressStep, 300, 80);
        
        parfor i=1:length(omega)
            fft_acf(i) = fft_acf(i)-sum(diff(exp(-1i*2*pi*omega(i)*oversample_tau)).*diff(oversample_acf)./diff(oversample_tau));
            if mod(i,progressStep)==0
                ppm.increment();
            end
        end
        ppm.delete()
        clear oversample_tau 
        clear oversample_acf
        acf_fft(:,4:5) = [real(fft_acf)',imag(fft_acf)'];
        acf_fft(:,4)=acf_fft(:,4)./(-4*pi*pi*acf_fft(:,1).^2);
        acf_fft(:,5)=acf_fft(:,5)./(-4*pi*pi*acf_fft(:,1).^2);
        clear fft_acf
        
        dlmwrite([folder1,'FFT-ACF_',filename],acf_fft, 'delimiter', ',', 'precision', 9);
    else
        error('Wrong FFT parameter!')
    end
    acf_fft(:,[2,4])=-acf_fft(:,[2,4]);
    acf_fft(:,[3,5])=-(acf_fft(:,[3,5]) + 1./(2*pi*acf_fft(:,[1,1])));
    
    ampl=acf_fft(:,[2,4]).^2 + acf_fft(:,[3,5]).^2;
    omega=acf_fft(:,1);
    Gpx = -rheofactor(1)*(acf_fft(:,3)./(2*pi*acf_fft(:,1).*ampl(:,1))); %Uncorrected for trap for visualization purposes
    Gppx = -rheofactor(1)*acf_fft(:,2)./(2*pi*acf_fft(:,1).*ampl(:,1));
    Gpy = -rheofactor(2)*(acf_fft(:,5)./(2*pi*acf_fft(:,1).*ampl(:,2))); %Uncorrected for trap for visualization purposes
    Gppy = -rheofactor(2)*acf_fft(:,4)./(2*pi*acf_fft(:,1).*ampl(:,2));
    clear acf_fft
    %%
    ax3=subplot(2,2,2);
    loglog(omega,Gpx,'*')
    hold on
    loglog(omega,Gppx,'x')
    loglog(omega,rheofactor(1)*ones(length(omega),1),'--k')
    hold off
    xlabel('Frequency(Hz)')
    ylabel('G'',G'''' (Pa)')
    legend('G''','G''''',['Trap: ',num2str(rheofactor(1),3),' Pa'],'Location','northwest')
    
    ax4=subplot(2,2,4);
    loglog(omega,Gpy,'*')
    hold on
    loglog(omega,Gppy,'x') 
    loglog(omega,rheofactor(2)*ones(length(omega),1),'--k')
    hold off
    xlabel('Frequency (Hz)')
    ylabel('G'',G'''' (Pa)')
    legend('G''','G''''',['Trap: ',num2str(rheofactor(2),3),' Pa'],'Location','northwest')
    linkaxes([ax3,ax4],'x')
    
    dlmwrite([folder1,'Gmoduli_',filename],[omega,Gpx,Gpy,Gppx,Gppy], 'delimiter', ',', 'precision', 9);


end
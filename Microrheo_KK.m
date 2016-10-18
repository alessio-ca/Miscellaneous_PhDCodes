%Input: enter mother directory name (it must end with a "\" sign)
folder1='/Users/AleCac/Documents/Python_stuff/';%folder = 'F:\CGMuRheoFile\';

%Input: enter the last stamp of the track_output files (a list if more than
%one is present)

cell_folder_list={'60x_1.20_Si02_1000FPS_1.5um_Tweezer_G_0.05.13Aug2016_15.27.31'};
cell_folder_list={'Simulation_1000FPS_Kel_1e-6_D_2e-6'};
listlen=length(cell_folder_list);

%Input: enter frame rate and pixel-to-nm conversion factor
fps=1000;

%Input: KK calculation
KK=0;

%Input: enter parameters for PSD characterization
kB = 1.3806503 * 1E-23;
Temp = 298;
Dbead = 2e-6;
eta = 1e-3;

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
    
    [psdx,~,~,~] = psd(1/fps,traj(:,2));
    [psdy,freq,~,~] = psd(1/fps,traj(:,3));
    
    
    gauss_points=21;
    h=gausswin(gauss_points);
    h=h/sum(h);
    
    interval=ceil(gauss_points/2):gauss_points:(gauss_points*(floor(length(psdx)/gauss_points-1))+ceil(gauss_points/2));
    freqsmooth=freq(interval);
    smoothpsdx=conv(psdx,h,'same');
    smoothpsdx=smoothpsdx(interval);
    smoothpsdy=conv(psdy,h,'same');
    smoothpsdy=smoothpsdy(interval);
    
    
    f_c=@(x) sqrt((4.14e-21/(3*pi*eta*Dbead))./(pi*pi*x)); 
    k_c=@(x) 2*pi*(3*pi*eta*Dbead)*x; 
    
    %Plotting section
    figure(1)
    ax1=subplot(2,2,1);
    figure(1)
    loglog(freqsmooth(2:end),smoothpsdx(2:end),'o')
    hold on
    
    fitobj_inx=fit(freqsmooth(2:35),smoothpsdx(2:35),'poly1');
    fc=1e6*f_c(fitobj_inx.p2); %In Hz
    kcx=k_c(fc); %In SI units
    loglog(freqsmooth(2:350),fitobj_inx.p2*ones(size(freqsmooth(2:350),1),1),'--r','LineWidth',2)
    fitobj_finx=fit(log(freqsmooth(347:end)),log(smoothpsdx(347:end)),'poly1');
    loglog(freqsmooth(85:end),exp(log(freqsmooth(85:end))*fitobj_finx.p1+fitobj_finx.p2*ones(size(freqsmooth(85:end),1),1)),'--g','LineWidth',2)
    hold off
    xlabel('Frequency (Hz)')
    ylabel('PSD (\mum^2/Hz)')
    legend('Data','D/(2\pi^2) {f_c}^{-2}','D/(2\pi^2) {f}^{ -2}')
    
    ax2=subplot(2,2,3);
    loglog(freqsmooth(2:end),smoothpsdy(2:end),'o')
    hold on
    fitobj_iny=fit(freqsmooth(2:35),smoothpsdy(2:35),'poly1');
    fc=1e6*f_c(fitobj_iny.p2); %In Hz
    kcy=k_c(fc); %In SI units
    loglog(freqsmooth(2:350),fitobj_iny.p2*ones(size(freqsmooth(2:350),1),1),'--r','LineWidth',2)
    fitobj_finy=fit(log(freqsmooth(347:end)),log(smoothpsdy(347:end)),'poly1');
    loglog(freqsmooth(85:end),exp(log(freqsmooth(85:end))*fitobj_finy.p1+fitobj_finy.p2*ones(size(freqsmooth(85:end),1),1)),'--g','LineWidth',2)
    hold off
    xlabel('Frequency (Hz)')
    ylabel('PSD (\mum^2/Hz)')
    legend('Data','D/(2\pi^2) {f_c}^{-2}','D/(2\pi^2) {f}^{ -2}')
    linkaxes([ax1,ax2],'x')

    
    if KK==0
        try KK=dlmread([folder1,'KK_',filename],',');
            disp('KK file found.')
        catch ME
            if strcmp(ME.identifier,'MATLAB:dlmread:FileNotOpened')
                error('KK file not found. Re-run with the KK calculation enabled if you wish to continue with the microrheological analysis.');
            else
                rethrow(ME)
            end
        end
    elseif KK==1
        
        disp('Kramers-Kronig calculation in progress...')
        
        interval=(1:length(freqsmooth));
        freqsample=freqsmooth(interval);
        psdxsample=smoothpsdx(interval);
        psdysample=smoothpsdy(interval);
        
        AImx = 2*pi*freqsample.*psdxsample/(2*kB*1e12*Temp); %In SI units
        ARex = kkrebook2(2*pi*freqsample,AImx,0);
        AImy = 2*pi*freqsample.*psdysample/(2*kB*1e12*Temp); %In SI units
        ARey = kkrebook2(2*pi*freqsample,AImy,0);
        
        KK=[freqsample,AImx,ARex',AImy,ARey'];
        clear AImx
        clear ARex
        clear AImy
        clear ARey
        
        dlmwrite([folder1,'KK_',filename],KK, 'delimiter', ',', 'precision', 9);
    else
        error('Wrong KK parameter!')
    end
    
    freqsample=KK(:,1);
    GmREALx=(KK(:,3)./(KK(:,3).^2 + KK(:,2).^2)) / (3*pi*Dbead);
    GmIMAGx=(KK(:,2)./(KK(:,3).^2 + KK(:,2).^2)) / (3*pi*Dbead);
    GmREALy=(KK(:,5)./(KK(:,5).^2 + KK(:,4).^2)) / (3*pi*Dbead);
    GmIMAGy=(KK(:,4)./(KK(:,5).^2 + KK(:,4).^2)) / (3*pi*Dbead);
    
    gauss_points=5;
    h=gausswin(gauss_points);
    h=h/sum(h);
    interval=ceil(gauss_points/2):gauss_points:(gauss_points*floor(length(freqsample)/gauss_points-1)+ceil(gauss_points/2));
    omega=freqsample(interval);
    Gpx=conv(GmREALx,h,'same');
    Gpx=Gpx(interval);
    Gpy=conv(GmREALy,h,'same');
    Gpy=Gpy(interval);
    Gppx=conv(GmIMAGx,h,'same');
    Gppx=Gppx(interval);
    Gppy=conv(GmIMAGy,h,'same');
    Gppy=Gppy(interval);
    
    ax3=subplot(2,2,2);
    loglog(omega,Gpx,'*')
    hold on
    loglog(omega,Gppx,'x')
    loglog(omega,(kcx/(3*pi*Dbead))*ones(length(omega),1),'--k')
    hold off
    xlabel('Frequency(Hz)')
    ylabel('G'',G'''' (Pa)')
    legend('G''','G''''',['Trap: ',num2str(kcx/(3*pi*Dbead),3),' Pa'],'Location','northwest')
    
    ax4=subplot(2,2,4);
    loglog(omega,Gpy,'*')
    hold on
    loglog(omega,Gppy,'x')
    loglog(omega,(kcy/(3*pi*Dbead))*ones(length(omega),1),'--k')
    hold off
    xlabel('Frequency(Hz)')
    ylabel('G'',G'''' (Pa)')
    legend('G''','G''''',['Trap: ',num2str(kcy/(3*pi*Dbead),3),' Pa'],'Location','northwest')
    linkaxes([ax3,ax4],'x')
    dlmwrite([folder1,'KKGmoduli_',filename],[omega,Gpx,Gpy,Gppx,Gppy], 'delimiter', ',', 'precision', 9);

end


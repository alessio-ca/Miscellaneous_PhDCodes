%This routine executes a single or multi-particle tracking for bright/dark field
%images in which the particle has a defined blob appearance (for instance,
%FLUORESCENCE MICROSCOPY pictures or BRIGHT FIELD MICROSCOPY pictures on
%small particles). It applies a standard centroid location algorithm with
%sub-pixel resolution (see Track_Fin.m documentation for more details).

%INPUT: requires knowledge of the approximate DIAMETER of the blobs.
%REQUIRES (STRONGLY SUGGESTED) pretracking routines SignalToNoise.m and
%ImageFiltering.m to get a good value of the blob's diameter in order to
%maximize the signal-to-noise ratio.


%Input: enter mother directory name (it must end with a "\" sign)
folder1='F:\AC_2016_07_18\5.0mmSDS\Acqui4\';
%Input: enter child directory name
folder2='60x_1.20_PPB_OD_FPS20_500nmPSR_g(r)*';
%Input: enter the diameter of the blobs (should be the same as chosen
%during the pretracking routines and shouldn't change during the series)
featdiam=4;

list = dir([folder1,folder2,'*']);
for i=1:length(list)
    folder2=list(i).name;
    stamp=[folder2,'_'];
    
    listlen = length(dir([folder1,folder2,'\','*.tiff']));
    %Input: by default, the program process all the pictures. Modify the 2nd
    %and 3rd field if necessary
    featsize = floor(featdiam/2);
    if mod(featsize,2) == 0
        featsize = featsize + 1;
    end
    poslist=Ftrack_Fin([folder1,folder2,'\',stamp],1,listlen-1,featsize);
    
    %Input: by default, the program will create an output file in the mother
    %directory. Modify if necessary
    fullFileName=fullfile(folder1,['track_particle_output',folder2,'.csv']);
    if exist(fullFileName,'file')==2
        %Input to specify how to deal with pre-existing file
        loop='A';
        while (loop ~= 'Y') && (loop ~= 'N')
            loop=input('Output track file already exists. Do you want to overwrite it? [Y/N] ','s');
            if loop=='Y'
                delete(fullFileName);
            elseif loop=='N'
                while exist(fullFileName,'file')==2
                    fullFileName = [fullFileName(1:end-4),'bis.csv'];
                end
            else
                disp('Incorrect input from user.')
                loop = 'A';
            end
        end
    end
    
    dlmwrite(fullFileName, poslist, 'delimiter', ',', 'precision', 9); 
    
%     param =struct('mem',0,'good',100,'dim',2,'quiet',0);
%     pos=poslist(:,[2:3,1]);
%     linked =track(pos,10,param);
    
    
    if (i<length(list))
        loop='A';
        while (loop ~= 'Y') && (loop ~= 'N')
            loop=input('Inspect the output file. Continue? [Y/N] ','s');
            if loop=='Y'
                continue
            elseif loop=='N'
                
            else
                disp('Loop interrupted.')
                return
            end
        end
    else
        disp('Everything is done!')
    end
end

%This routine executes a single-particle tracking for bright field images
%in which the particle has a doughnut-like appearance. It applies an
%2D autocorrelation filter with a custom filter input by the user and
%paraboloid fitting to get sub-pixel resolution (see AnalyseTrack_Parallel.m
%documentation for more details).
%It DOESN'T require pre-tracking procedures.
%It supports parallelization.
%It DOESN'T adapt the filter (required for parallelization). Frames which
%can't be tracked are lost (reported at the end of the execution).


%Filter syntax:
%rMax,rMid and rMin are mandatory and must be positive
%int1,int2 and int3 are optional
%Background intensity is always 0
%Outer rim is int1
%Mid rim is int2
%Inner rim is int3
%
%Default is: 1 (outer rim = 1)
%            -2 (mid rim = -2)
%            4 (inner rim = 4)
%



close all
%Input: enter directory name
folder_name=uigetdir('/Users/AleCac/Desktop/');

%Input: enter filename stamp
[fullpath,stamp,ext]=fileparts(folder_name);
stamp=[stamp,ext];
%Format this if you series doesn't start with 0 or has a different number of digits
img=imread([folder_name,'/',stamp,'_000000.tiff']); 
imshow(mat2gray(img));
[filter{1,1},filter{1,2},r_max]=CreateFilter(img);


listlen = length(dir([folder_name,'/','*.tiff']));

%Input: by default, the program process all the pictures. Modify the 2nd
%and 3rd field if necessary
[poslist,histostep]=AnalyseTrack_Parallel([folder_name,'/',stamp,'_'],1,listlen-1,filter,r_max*2); %The listlen-1 accounts for one frame being 00000
figure(1)
histogram(histostep);


%Input: by default, the program will create an output file in the mother
%directory. Modify if necessary
fullFileName=fullfile(fullpath,['track_particle_output',stamp,'.csv']);
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

%Plot section
figure(2)
plot(poslist(:,1),poslist(:,2))
hold on
plot(poslist(:,1),poslist(:,3))
hold off
title('Coordinates over time')
xlabel('Time step')
ylabel('Pixel')
legend('X coordinate','Y coordinate')


folder = 'F:\AC_2018_01_11\ME_F108_SiOil_Frames_OOF';
listFrames = dir(folder);
listFrames = listFrames([listFrames.isdir]==0);
diami = [];

for i=1:length(listFrames)
    disp(i);
    img  = imread([folder,'\',listFrames(i).name]);
    diami = [diami;DropletDetector(img,10,100,'canny',[0.7,1.],'dilDiam',3,'Rad',1.25)];
end
%%
%Conversion factors:
% -    10x 0.45 NA: 169.33 px to 100 um
% -    20x 0.75 NA: 342.67 px to 100 um
% -    40x 0.95 NA: 689.65 px to 100 um
% -    60x 1.20 NA: 1000.0 px to 100 um

convfac = 169.33;
h=histogram(diami*(100/convfac),'Normalization','probability','BinWidth',3);
disp(['Mean size: ',num2str(mean(diami*(100/convfac)))])
disp(['PDI: ',num2str(std(diami*(100/convfac))/mean(diami*(100/convfac)))])

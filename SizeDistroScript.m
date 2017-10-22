folder = 'J:\ME_F108_SiOil_Frames';
listFrames = dir(folder);
listFrames = listFrames([listFrames.isdir]==0);
diami = [];

for i=1:length(listFrames)
    disp(i);
    img  = imread([folder,'\',listFrames(i).name]);
    diami = [diami;DropletDetector(img,50,400,'canny',[0.6,7.],'dilDiam',5,'Rad',1.3)];
end
%%
%Conversion factors:
% -    10x 0.45 NA: 169.33
% -    20x 0.75 NA: 342.67
% -    40x 0.95 NA: 689.65 
convfac = 169.33;
h=histogram(diami*(100/convfac));
disp(['Mean size: ',num2str(mean(diami*(100/convfac)))])
disp(['PDI: ',num2str(std(diami*(100/convfac))/mean(diami*(100/convfac)))])

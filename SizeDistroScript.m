folder = 'F:\AC_2017_05_19\ME_F108_SiOil_Frames';
listFrames = dir(folder);
listFrames = listFrames([listFrames.isdir]==0);
diami = [];

for i=1:length(listFrames)
    disp(i);
    img  = imread([folder,'\',listFrames(i).name]);
    diami = [diami;DropletDetector(img,50,400,'canny',[0.4,5.])];
end
%%
histogram(diami*(100/169.33))
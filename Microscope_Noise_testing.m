[filename,fullpath]=uigetfile('*.movie','Select a file','C:');
[fileframes,filesize,filedata,fileheader]=movie2info([fullpath,filename]);
disp(['File is ',Bytes2str(filesize)]);
disp(['Number of frames: ',num2str(fileframes)]);
offset = 10000*filedata + fileheader;
[video,offset,nfr]=movie2frame([fullpath,filename],offset);
res=zeros(1,size(video,1));
x=zeros(1,size(video,1));
y=zeros(1,size(video,1));
for i=1:size(video,1)
    img=squeeze(video(i,:,:));
    res(i)=mean2(mat2gray(video(i,:,:)));
    [Y,iy]=max(img);
    [x(i),ix]=max(Y);
    y(i)=iy(ix);
end
%%
[f,spectrum]=Alessio_FFT(res,1000);
semilogy(f,spectrum)
%%
[f,spectrum]=Alessio_FFT(video(:,46,48),1000);
semilogy(f,spectrum)
%%
clear spectrum
for i=1:size(img,1)
    for j=1:size(img,2)
        [f,spectrum]=Alessio_FFT(mat2gray(video(:,i,j)),1000);
        noisepeak_100(i,j)=mean(spectrum(10006:10016));
        noisepeak_120(i,j)=mean(spectrum(12000:12010));
    end
end
%%
figure(1)
surf(noisepeak_100)
figure(2)
surf(noisepeak_120)
figure(3)
surf(img)
%%
[f,spectrum]=Alessio_FFT(y,1000);
semilogy(f,spectrum)
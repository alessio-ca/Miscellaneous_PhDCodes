%------------------------------ 1. Input and some starting basic calculations. ------------------------------

tic

tmp=load(fname);
names=fieldnames(tmp);                           % The curly brackets in the names{1,1} construct are necessary because its 
A=getfield(tmp,names{1});                        % a cell array; Matlab 6.5 requires(?) A=tmp.(names{1,1});
clear tmp;
clear names;
nPix=size(A,1);
nFrames=size(A,3);                               % simulation reads droplets from txt file, generates spheres  
startRow = 1;                                    % droplets and the data_file defining the simulated droplets.
randn('state', 0); rand('state', 0);             % sample centre: (nPix+1)/2, (nPix+1)/2
disp('boxX   boxY   boxZ '), disp([boxX boxY boxZ])
disp('nPix   nFrames '), disp([nPix nFrames])

    for jj =1:maxjj, R(jj) = rFirst + stepR * (jj - 1); end

disp('       Rmin    Rnum     ThreshMaxI'), disp([R(maxjj) maxjj  ThreshMaxI])

midZ = double(uint8(nFrames/2));            % middle of the data box
midY = round(nPix/2);
ddx = boxX/nPix
ddy = boxY/nPix
ddz = boxZ/nFrames

vertPlane = zeros(nFrames,nPix);
endresult=zeros(maxRegionsPerRadius, 4, maxjj);
conn3d = 26;


%------------------------------2. Drawing horizontal, vertical planes of real and fourier transform pics for middle planes of whole stack.-----


figure, imagesc(A(:,:,midZ),[0 255]); title ('input spheres'); colormap(gray);colorbar;axis('image');
drawnow

for i = 1:nFrames
    for j = 1:nPix
        vertPlane(i,j) = A(midY,j,i);
    end
end
figure, imagesc(vertPlane,[0 255]); title ('input spheres Z'); colormap(gray);colorbar;axis('image');
daspect([ddz ddx 1])                                    % for vertical plane only
drawnow

[mg sg qg] = ndgrid((1:nPix)-0.5,(1:nPix)-0.5,(1:nFrames)-0.5);
ipDist = sqrt(((mg-(nPix)/2)*ddx).^2 + ((sg-(nPix)/2)*ddy).^2 + ((qg-(nFrames)/2)*ddz).^2); % IJ: looks distances from center.
save ipdist.mat ipDist;             
clear ipDist;
clear mg; clear sg; clear qg;
fii = fftshift(fftn(fftshift(A)));                      % fourier transform of the box with all spheres in ii. First does shift of A, than fft, and then again shift.
    clear A;

figure, imagesc(abs(fii(:,:,midZ))); title ('fii'); colorbar;axis('image');
drawnow
for i = 1:nFrames
    for j = 1:nPix
        vertPlane(i,j) = abs(fii(midY,j,i));
    end
end
figure, imagesc(vertPlane); title ('fii Z'); colorbar;axis('image');
daspect([ddz ddx 1])                                    % when xz used
drawnow
save fii.mat fii '-v7.3';
clear fii;



%-------3. Make one sample sphere in 3D in the center of stack by making intensity 200 plus normal noise. Draw horiz. and vert plane of it.--
   
for jj = 1:maxjj
        disp('Test Radius'), disp(R(jj));
    ip = randn(nPix, nPix, nFrames);
    ip = abs(ip) * maxNoise*0.01;
    load ipdist.mat;
    iii = find(ipDist <= R(jj));                        % iii = find(ipDist <= 1 & ipDist > 0.9); %ring
    clear ipDist;
    ip(iii) = 100 + ip(iii);                            % IJ: add to sphere 200 and normally generated noise. Was 100 before.
    clear iii;
    if jj == 1
        figure, imagesc(ip(:,:,round(nFrames/2)),[0 255]); colormap(gray);colorbar;axis('image'); 
        drawnow
        for i = 1:nFrames
            for j = 1:nPix
                vertPlane(i,j) = ip(round(nPix/2),j,i);
            end
        end
        figure, imagesc(vertPlane,[0 255]); colormap(gray);colorbar;axis('image');
        daspect([ddz ddx 1]) % when xz used
        drawnow
        clear vertPlane;
    end
    
    
%------------------------- 4. Find fft of the sample droplet and do Wiener filtering of data for particular radius. --------------
    
    fip = fftshift(fftn(fftshift(ip)));   %fourier transform of the single sphere in 3D
    clear ip;
   
    Rf = fip.*conj(fip);                    %filter high frequency due to sharp interfaces. if you don't do it you would get whole drops instead of only edges in the figure below.
    Rf = Rf./(Rf+RfConst);  % use to be fil1 = ...   Wiener filter multiplier of a single sample sphere.
    
    nZfip = find (fip == 0);        % IJ: I think it doesn't really happen.
    if nZfip 
        for i=1:nPix
            disp ('Zeros in fip!!!!!')
            for j=1:nPix
                for k=1:nFrames
                    if abs(fip(i,j,k)) > 0 
                        Rf(i,j,k) = (fii(i,j,k) / fip(i,j,k)) * Rf(i,j,k);  % IJ: Was nii = ... I changed to Rf = ..., same below.
                    else
                        Rf(i,j,k) = (fii(i,j,k) / 0.1 ) * Rf(i,j,k);   % IJ changed, was: nii = fii(i,j,k) * Rf(i,j,k);
                    end
               end
           end
        end
    else
        Rf = Rf./fip;
        clear fip
        load fii.mat;
        Rf = fii.*Rf;    % IJ: this is Wiener filter H.
        clear fii
          
    end
   
    Rf = abs(fftshift(ifftn(fftshift(Rf))));
    nii = Rf;
    clear Rf;
    if jj == 1
         
        figure, imagesc(nii(:,:,midZ)); title ('nii'); colorbar;axis('image');  % What you get after Wiener filtering.
        drawnow

    end
    
    
    
    %------ 5. For given radius you request that number of accepted pixels to be not more than maxnumPixels and same for maxnumobjects. I removed it and left thresholding only by intensity cause i don't know what should be the number of objects.
    
    bw = nii(:,:,1:nFrames);
    c = max(bw,[],2);
    c = max(c,[],1);
    c = max(c,[],3);                      %finding the max intensity of fft in 3d
    cc(jj) = c;             % Brightest pixel for this radius.
    disp('maxI')
    disp(c)
    i = 1;
   % numpixels = maxnumPixels + 1;     IJ: I simplified code below, need check.
     numobjects = maxnumObjects + 1; 
   % while numpixels > maxnumPixels 
   %     bw = nii(:,:,1:nFrames) > c * Thresholding*(0.9 + i * 0.1); %?? bwlabeln works on binary ??
   %     numpixels = nnz(bw);
   %     disp(' numPixels #loop');
   %     disp([numpixels i]);
   %     i = i +1;
   % end 
   % i = i - 1;
    while numobjects > maxnumObjects 
        bw = nii(:,:,1:nFrames) > c * Thresholding*(0.9 + i * 0.1); %?? bwlabeln works on binary ??
        bnii = bw.*nii(:,:,1:nFrames);
        numpixels = nnz(bw);
        disp(' numPixels #loop');
       disp([numpixels i]);
        clear bw;
        [L,numobjects] = bwlabeln (bnii(:,:,:), 26);
        L = uint16(L);
        disp(' numObjects #loop');
        disp([numobjects i]);
        i = i +1;
    end 
    clear nii;
    
    
    
    %------------------------------ 6. Calculate CM of each found and filtered in section above object.--------------------------
    
    stats = regionprops(L,'Area','PixelList');  % centre of mass
    clear L;
    xcmplot=zeros(numobjects,3);
    indexM = 1;
    for i=1:numobjects                          % for every object whose maxI > c*ThreshMaxI compute CM
        ball=[stats(i).PixelList];
        ar=[stats(i).Area];
        intensity = 0;
        sort_temp = zeros(1,ar);                % form intensity vector                   
        for ss = 1:ar
            sort_temp(ss) = bnii(ball(ss,2),ball(ss,1),ball(ss,3));
        end
        [sorted, idx] = sort(sort_temp);        % sort intensity (maximum is last)
        if sorted(ar) > c * ThreshMaxI      
            xcm=zeros(1,3);
            for j=1:ar                          % for every voxel
                intensity = intensity + bnii(ball(j,2), ball(j,1), ball(j,3));
            end
            for k=1:3
                for j=1:ar
                    xcm(1,k) = xcm(1,k) + bnii(ball(j,2),ball(j,1), ball(j,3))*ball(j,k);  % Calc CM
                end     
            end
            for k=1:3
                xcm(1,k) = xcm(1,k)/intensity;      % Normalize.
                xcmplot(i,k) = xcm(1,k);        % gives centre of mass weighted by intensity
            end
            tmpresult(indexM,1)=xcmplot(i,1) - .5;    % Shifting CM pixels to their centers.
            tmpresult(indexM,2)=xcmplot(i,2) - .5;
            tmpresult(indexM,3)=xcmplot(i,3) - .5;
            tmpresult(indexM,4)=sorted(ar);         % IJ: Intensity of the brightest pixel for this object not normalized by the brightest pix for this radius.
            indexM = indexM + 1;                % store only if maxI > ThreshMxI * c
        end
   end
   clear bnii;
   
%--------- 7. Removes all the deltas by intensities that are beyond maxRegionperRadius. I think ThreshMax and Thresholding should be set well enough already to avoid this but seems like it makes no harm, left. 
   
   sorted = sortrows(tmpresult,4);              % sorts tmpresult
   endresult(maxRegionsPerRadius, 1, jj) = 0;
   if indexM >= maxRegionsPerRadius
       %disp('too many filtered objects!')
       for i = 1:maxRegionsPerRadius-1
           endresult(i,:,jj) = sorted(indexM-maxRegionsPerRadius+i,:);    % Save only brightest ones.
       end
       disp(endresult(maxRegionsPerRadius-5:maxRegionsPerRadius-1,:,jj))    % show 5 only  
       %disp(endresult(:,:,jj))    % show 5 only  
   else
       for i = 1:indexM-1
           endresult(i,:,jj) = sorted(i,:);
       end
       if indexM > 5
           disp(endresult(indexM-5:indexM-1,:,jj)) 
       else
           disp(endresult(1:indexM-1,:,jj))
       end
   end
   clear tmpresult;
end

%------------------------ 8. Saving results all together for all radii into grandresult. Combining close center graphs. --------------

clear vertPlane;
grandresult = zeros(1,5);
k = 1;
for i = 1:maxjj
    for j = 1:maxRegionsPerRadius
        if endresult(j,1,i)
            for ik = 1:3, grandresult(k, ik) = endresult(j,ik,i); end
            grandresult(k,4) = R(i);
            grandresult(k,5) = endresult(j,4,i);
            k = k + 1;
        end
    end
end
clear endresult
sortedX = sortrows(grandresult,1);          %sort in X
totD = size(sortedX,1);
gres = zeros(1, 5);
combined = zeros(1, totD);
deleted = zeros(1, totD);
k = 1;
for i = 1 : totD                                            % combine spheres having close centres  
    if ~deleted(i)
        for j = i + 1:totD
            if ~deleted(j)
                 %if (ddx*abs(sortedX(j,1) - sortedX(i,1)) + ddy*abs(sortedX(j,2) - sortedX(i,2)) + ddz*abs(sortedX(j,3) - sortedX(i,3))) > (min(sortedX(i,4),sortedX(j,4)))  % IJ: rstricted condition so that as many as possible spheres get rejected. 
                 %   break;              % IJ: Not continue because first j comes right after i and supposed to be the only possible same one as i.
                 %else
                    dsq = (sortedX(i,1) - sortedX(j,1))*(sortedX(i,1) - sortedX(j,1))*ddx*ddx;
                    dsq = dsq + (sortedX(i,2) - sortedX(j,2))*(sortedX(i,2) - sortedX(j,2))*ddy*ddy;
                    dsq = sqrt(dsq + (sortedX(i,3) - sortedX(j,3))*(sortedX(i,3) - sortedX(j,3))*ddz*ddz);   % IJ: added sqrt outside, was without.
                    temp = min(sortedX(i,4),sortedX(j,4));
                    if dsq < temp
                        if sortedX(i,5) < sortedX(j,5)
                            deleted(i) = 1;
                            combined(j) = combined(j) + combined(i) + 1;
                            combined(i) = 0;            
                            break;
                        else
                            combined(i) = combined(i) + combined(j) + 1;
                            deleted(j) = 1;
                            combined(j) = 0;    % IJ: I added this.
                        end
                    end
                %end
            end
        end
        if combined(i) >= minRadiiNum - 1   %remove spheres belonging to one radius only. -1 to avoid itself.
            gres(k,:) = sortedX(i,:);
            k = k + 1;
        end
    end
end
multS = k - 1;
grandresult = zeros(1,5);
grandresult = sortrows(gres,5);
disp('grandresult (pixels)')
disp(grandresult)                                         % print results   
disp('#spheres')
disp(multS)

toc

save(FDres)

save(out_test, 'grandresult', '-ascii')

 disp('Hit any key to proceed to choosing drops by hand')
 pause

 
 
%------------------------------ 9. UI accept/reject visualization.---------------------------------
                                 % read / construct 'A' again for visualisation
    tmp=load(fname);
    names=fieldnames(tmp);
    A=getfield(tmp,names{1});   %Matlab 6.5 requires(?) A=tmp.(names{1,1});

j = 1;
acc = 0; tmp = 1;
i = multS;


while(i >= 1)              % for all spheres found
    disp([i acc])
    
    x = double(uint16(grandresult(i,2)+0.5)); if x == 0, x = 1; end         %draw cross
    y = double(uint16(grandresult(i,1)+0.5)); if y == 0, y = 1; end
    z = double(uint16(grandresult(i,3)+0.5)); if z == 0, z = 1; end
    deltaPix = double(uint16(grandresult(i,4)/ddx+0.5))-1;
    if z > startRow - 1 && z < startRow + nFrames 
        zd = z;
    elseif z < startRow 
        zd = startRow;
    elseif z >= startRow + nFrames
        zd = startRow + nFrames - 1;
    end
    for n = -deltaPix:deltaPix
        if  x+n > 0 && x+n < nPix
            xt(n+deltaPix+1) = A(x+n, y, zd);
            A(x+n, y, zd) = 0;                %black cross
        end
        if y+n > 0 && y+n < nPix
            yt(n+deltaPix+1) = A(x, y+n, zd);
            A(x, y+n, zd) = 0;                %black cross
        end
    end
    figure, imshow(A(:,:, zd),[0 255]); title(['x= ', num2str(grandresult(i,2)),' y= ', num2str(grandresult(i,1)),' z= ', num2str(grandresult(i,3)-startRow+1),' r= ', num2str(grandresult(i,4)/ddx),' i= ', num2str(grandresult(i,5))]);  % IJ: removed colorbar, gave errors don't know why.
    tmp = input ('accept? (ent/1/4/2)');                                    % accept, reject, reject all
    close
    
    
    if tmp == 4    % Accepts the current sphere and inputed # of spheres after the current one. No white crosses on these spheres.
        tmp1 = input('Please enter the number of sphere in the raw u would like to accept besides the current one (should be >=0). \n');
        selectedSpheres(j:(j+tmp1),:) = grandresult(i:-1:(i-tmp1),:);
        
        for k=i:-1:(i-tmp1)         
            x = double(uint16(grandresult(k,2)+0.5)); if x == 0, x = 1; end         %draw cross
            y = double(uint16(grandresult(k,1)+0.5)); if y == 0, y = 1; end
            z = double(uint16(grandresult(k,3)+0.5)); if z == 0, z = 1; end
            deltaPix = double(uint16(grandresult(k,4)/ddx+0.5))-1;
            if z > startRow - 1 && z < startRow + nFrames 
                zd = z;
            elseif z < startRow 
                zd = startRow;
            elseif z >= startRow + nFrames
                zd = startRow + nFrames - 1;
            end
            
            for n = -deltaPix:deltaPix 
                if x+n > 0 && x+n < nPix,  A(x+n, y, zd) = 255; end    % ? cross
                if y+n > 0 && y+n < nPix,  A(x, y+n, zd) = 255; end    % ? cross
            end
            
        end
        clear x y z zd n deltaPix
        
        j   = (j+tmp1) + 1;
        acc = acc + tmp1 + 1;
        i   = i - (tmp1+1);
                
    
    elseif tmp == 3    % Accepts the sphere if tmp=3
        selectedSpheres(j,:) = grandresult(i,:);
        for n = -deltaPix:deltaPix 
            if x+n > 0 && x+n < nPix,  A(x+n, y, zd) = 255; end    % ? cross
            if y+n > 0 && y+n < nPix,  A(x, y+n, zd) = 255; end    % ? cross
        end
        j   = j + 1;
        acc = acc + 1;
        i   = i - 1;
    elseif tmp == 1   % Rejectes the sphere if tmp=1
        for n = -deltaPix:deltaPix
            if x+n > 0 && x+n < nPix,  A(x+n, y, zd) = xt(n+deltaPix+1); end
            if y+n > 0 && y+n < nPix,  A(x, y+n, zd) = yt(n+deltaPix+1); end
        end
        i   = i - 1;
    elseif tmp == 2  % Rejects this and all the following spheres if tmp =2.
        break
    else             % Accepts the sphere if tmp is not equal to one of the stated above.         
        selectedSpheres(j,:) = grandresult(i,:);
        for n = -deltaPix:deltaPix 
            if x+n > 0 && x+n < nPix,  A(x+n, y, zd) = 255; end    % ? cross
            if y+n > 0 && y+n < nPix,  A(x, y+n, zd) = 255; end    % ? cross
        end
        j   = j + 1;
        acc = acc + 1;
        i   = i - 1;
    end
end
disp('Selected Spheres (pixels)')
for i = 1:j -1      
    selectedSpheres(i,3) = (selectedSpheres(i,3) - startRow + 1);
end
save(out_file, 'selectedSpheres', '-ascii')
disp(selectedSpheres)                                         % print results   
for i = 1:j -1      %transform into world coordinates
    worldresult(i,1) = (selectedSpheres(i,1) - 0.5) * ddx;
    worldresult(i,2) = (selectedSpheres(i,2) - 0.5) * ddy;
    worldresult(i,3) = (selectedSpheres(i,3) - 0.5) * ddz;
    worldresult(i,4) = selectedSpheres(i,4);          % R(index_c(maxjj-nSpR(r0,i,1)+1));
end
disp('world result (world coord.)')
disp(worldresult)                                         % print results 
disp('#selected')
disp(j - 1)

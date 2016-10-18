%Currently has projection of crosses through z-slices turned off

titl = 'Series001'
fname = [titl '.mat']      % Output of bigstackprocJB.
out_file = [titl '.txt']   % Output of FD.
out_test = [titl '_test.txt']  % Outputs all drops with garbage, but can check in gui if they at least cover mosto of the spheres.
FDres    = [titl '_FDres.mat']; % Think this is useless
nSn      = 1.403/1.474 ; %1.405/1.45; Snells law correction
boxX=50; % Box in xy
boxY=boxX;
boxZ=nSn*50.14;            % Each step is boxZ/(nslices - 1). I here divided since the original nSlices was 213.


%------------------------------ 1. Input and some starting basic calculations. ------------------------------

tmp=load(fname);
names=fieldnames(tmp);
A=getfield(tmp,names{1});   %Matlab 6.5 requires(?) A=tmp.(names{1,1});

grandresult=load(out_test);

nPix=size(A,1);
nFrames=size(A,3);
startRow = 1;

ddx = boxX/nPix
ddy = boxY/nPix
ddz = boxZ/nFrames

%------------------------------ 9. UI accept/reject visualization.---------------------------------
                                 % read / construct 'A' again for visualisation

j = 1;
acc = 0; tmp = 1;
i = size(grandresult(:,1),1);


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
    tmp = input ('Is this a droplet which I see before me? (ent/1/4/2)');                                    % accept, reject, reject all
    close
    
    
    if tmp == 4    % Accepts the current sphere and inputed # of spheres after the current one. No white crosses on these spheres.
        tmp1 = input('Please enter the number of sphere in the raw u would like to accept besides the current one (should be >=0). \n');
        selectedSpheres(j:(j+tmp1),:) = grandresult(i:-1:(i-tmp1),:);
        
        for k=i:-1:(i-tmp1)         
            x = double(uint16(grandresult(k,2)+0.5)); if x == 0, x = 1; end         %draw cross
            y = double(uint16(grandresult(k,1)+0.5)); if y == 0, y = 1; end
            z = double(uint16(grandresult(k,3)+0.5)); if z == 0, z = 1; end
            deltaZ = double(uint16(grandresult(k,4)/ddz+.5))-1;
            if z > startRow - 1 && z < startRow + nFrames 
                zd = z;
            elseif z < startRow 
                zd = startRow;
            elseif z >= startRow + nFrames
                zd = startRow + nFrames - 1;
            end
            for m = 0:0
                if zd + m > 0 && zd + m < nFrames
                   deltaPix = double(uint16(sqrt(double(uint16(grandresult(k,4)/ddx+.5)-1)^2 - double(uint16(m*(ddz/ddx)+.5)-1)^2)));
                   for n = -deltaPix:deltaPix 
                      if x+n > 0 && x+n < nPix,  A(x+n, y, zd+m) = 255; end    % ? cross
                      if y+n > 0 && y+n < nPix,  A(x, y+n, zd+m) = 255; end    % ? cross
                   end
                end
            end     
        end
        clear x y z zd n deltaPix
        
        j   = (j+tmp1) + 1;
        acc = acc + tmp1 + 1;
        i   = i - (tmp1+1);
                
    
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
        deltaZ = double(uint16(grandresult(i,4)/ddz+.5))-1;
        for m = 0:0
            if zd + m > 0 && zd + m < nFrames
               deltaPix = double(uint16(sqrt(double(uint16(grandresult(i,4)/ddx+.5)-1)^2 - double(uint16(m*(ddz/ddx)+.5)-1)^2)));
                   for n = -deltaPix:deltaPix 
                      if x+n > 0 && x+n < nPix,  A(x+n, y, zd+m) = 255; end    % ? cross
                      if y+n > 0 && y+n < nPix,  A(x, y+n, zd+m) = 255; end    % ? cross
                   end
            end
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

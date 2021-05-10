close all, clear all

% documents
% 
% inverse integrate gradient : https://fr.mathworks.com/matlabcentral/fileexchange/9734-inverse-integrated-gradient
% 
% PIVlab
% https://blog.espci.fr/mecaflu/files/2019/02/PIVlab_tuto.pdf


%% launch PIVlab_GUI
%% read results
%% run fhat = intgrad2(U,V)
 %%
cd('C:\Users\darcy\Desktop\TP_Image')
load('set02_movie01_all_frames.mat')
figure
U =  u_original{40, 1} ;
V =  v_original{40, 1} ;
quiver(U,V)
%%

figure, box on, hold on
set(gcf,'position',[680         203        1193         775])
view(-14,43)
for it = 1 : 1 : size(u_original,1)
    if exist('hs')
        delete(hs)
    end
    U =  u_original{it, 1} ;
    V =  v_original{it, 1} ;
    
    % remove the NaNs
    U(find(isnan(U))) = 0;     
    V(find(isnan(V))) = 0;
    
    fhat = intgrad2(U,V);
    h0 = 40; % (mm)
    h = sqrt(h0^2 - ( 2/0.25 ) * fhat );
    hs = surf(h,'edgecolor','none');
    title(sprintf('time: %0.0f',it))
    pause(.2)
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OLD STUFF 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% generate image sequence for PIVlab


iFolder = 'C:\Users\darcy\Desktop\TP_Image\set02\movie01';
oFolder = 'C:\Users\darcy\Desktop\TP_Image\set02\movie014PIVlab';
mkdir(oFolder)
cd(iFolder)
listIm = dir('*.tif')
cd(oFolder)
for i = 1 : length(listIm) 
    
end

%%

% define folder
cd('C:\Users\Lenovo\Jottacloud\ENSEIGNEMENT\2020_2021\L3_TP_images\ondesDeSurface\imseqBeadDrop_cutrot')

% load images
listImages = dir('*.png');
ImRef = imread(listImages(62).name);

%%
% figure
% hi = imagesc(imDiff);
hsurf = figure; hold on
 view(10,15)
clear dr imTL drTime
wtmplt = 10;
wRef = 2*wtmplt;
tic
for it =  24:26%1 : 5 : 66
    cd('C:\Users\Lenovo\Jottacloud\ENSEIGNEMENT\2020_2021\L3_TP_images\ondesDeSurface\imseqBeadDrop_cutrot')

    fprintf('current images being analysed: %0.0f \n',it)
    clear dr
im       = imread(listImages(it).name);
imDiff   = imsubtract(ImRef,im);

iix = 0;
iiy = 0;
% compute dr
for ix =   200 : 1 : 700  %200 : 10 : 700
    iix = iix + 1;
    iiy = 0;
    for iy = 250 : 260 %350%  250 : 1 : 500  %250 : 10 : 1000
        iiy = iiy + 1;
        [dx,dy,c] = fastCrossCorr(im,ImRef,ix,iy,wtmplt,wRef);
        dr(iiy,iix).x  = ix;
        dr(iiy,iix).y  = iy;
        dr(iiy,iix).dx = dx;
        dr(iiy,iix).dy = dy;
        dr(iiy,iix).dr = sqrt((dr(iiy,iix).dx)^2 + (dr(iiy,iix).dy)^2);
    end
end

%imTL(it,1:length(dr)) = [dr.dr];

% hi.CData = imDiff;
% pause(.2)

drTime(it).dr = dr;

 clear fhat fx fy
 for iy = 1 : size(dr,1)
     for ix = 1 : size(dr,2)
         fx(iy,ix) = dr(iy,ix).dx;
         fy(iy,ix) = dr(iy,ix).dy;
     end
 end
 cd('C:\Users\Lenovo\Jottacloud\ENSEIGNEMENT\2020_2021\L3_TP_images\ondesDeSurface\weFunction\integratedgradient')
fhat = intgrad2(fx,fy);
drTime(it).fhat = fhat;

% if it == 1
figure(hsurf)
 hs = surf(fhat,'edgecolor','none');
% else
%     hs.CData = fhat;
% end

title(sprintf('time is : %0.0f',it))
pause(.1)
end
toc
fprintf('correlation finished \n')


%% integrate the gradient
cd('C:\Users\Lenovo\Jottacloud\ENSEIGNEMENT\2020_2021\L3_TP_images\ondesDeSurface\weFunction\integratedgradient')

% fhat = intgrad2(fx,fy)

%  xp = 0:.1:1;
%  yp = [0 .1 .2 .4 .8 1];
%  [x,y]=meshgrid(xp,yp);
%  f = exp(x+y) + sin((x-2*y)*3);
%  [fx,fy]=gradient(f,.1,yp);
%  tic,fhat = intgrad2(fx,fy,.1,yp,1);toc

% Example usage 2: Large grid, 101x101
 xp = 0:.01:1;
 yp = 0:.01:1;
 [x,y]=meshgrid(xp,yp);
 f = exp(x+y) + sin((x-2*y)*3);
 [fx,fy]=gradient(f,.01);
 tic,fhat = intgrad2(fx,fy,.01,.01,1);toc


%%

%% browse results

figure, box on, hold on
set(gcf,'position',[680         203        1193         775])
view(-14,43)
for it = 1 : size(u_original,1)
    while(1)
    if exist('hs')
        delete(hs)
    end
    U =  u_original{it, 1} ;
    V =  v_original{it, 1} ;
    U(find(isnan(U))) = 0;     V(find(isnan(V))) = 0;
    fhat = intgrad2(U,V);
    hs = surf(fhat,'edgecolor','none');
    title(sprintf('time: %0.0f',it))

    end
end


him = figure('defaultAxesFontSize',20); hold on, box on


while contains('azesdxcr',get(gcf,'currentchar'))  % which gets changed when key is pressed
   figure(him)
   [x,y] = ginput(1);
   if get(gcf,'currentchar')=='z'
       iz = iz -  1;
   elseif get(gcf,'currentchar')=='e'
       iz = iz +  1;
   end
   
   clear xCH yCH
   xCH = [];
   yCH = [];
   for ich = 1 : length(channels)
       idx = find(iz == channels(ich).z);
       if idx
           xCH(ich) = channels(ich).x(idx);
           yCH(ich) = channels(ich).y(idx);
       end
   end
   
   
   figure(him), hold on
   im2D = im3D(:,:,iz);
   him2D.CData = im2D;
   axis(gca,[xAxis(1) xAxis(2) yAxis(1) yAxis(2)])
   hCH.XData = xCH;
   hCH.YData = yCH;
   %set(gcf,'position',[100 100 1000 1000])
   title(sprintf('z: %0.0f',iz))
   figure(hminiMap)
   hp.Vertices(:,3) = iz * [1,1,1,1,1];

end
 %%
 clear hs
cd('C:\Users\Lenovo\Jottacloud\ENSEIGNEMENT\2020_2021\L3_TP_images\ondesDeSurface\weFunction\integratedgradient')
hsurf = figure; hold on, box on
view(-11,23)
for it =  1 : 1 : 66
    if exist('hs')
        hs.FaceAlpha = .4;
    end
   figure(hsurf)
   hs = surf(drTime(it).fhat,'edgecolor','none');
   pause(.4)
end

%%
figure
%plot([dr.dr])
imagesc(imTL)
%%
imDR = zeros(size(dr,1),size(dr,2));
imDX = zeros(size(dr,1),size(dr,2));
imDY = zeros(size(dr,1),size(dr,2));
for ix = 1 : size(dr,1)
    for iy = 1 : size(dr,2)
        imDR(ix,iy) = dr(ix,iy).dr;
        imDX(ix,iy) = dr(ix,iy).dx;
        imDY(ix,iy) = dr(ix,iy).dy;
    end
end

clims = [0 6];
figure
imagesc(imDR,clims)

clims = [-6 6];
figure
imagesc(imDX,clims)

clims = [-10 10];
figure
imagesc(imDY,clims)

figure
plot(imDX(50,:))
%%
% build sample image
tmplt    = imcrop(im,    [ix-wtmplt iy-wtmplt 2*wtmplt-1 2*wtmplt-1]);
ImRefLoc = imcrop(ImRef, [ix-1*wRef iy-1*wRef 2*wRef 2*wRef]);
% do correlation
c = normxcorr2(tmplt,ImRefLoc);
% Gaussian interpolation
Nwidth = 1;
Ip = double(c(ypeak-Nwidth:ypeak+Nwidth,xpeak-Nwidth:xpeak+Nwidth));

% determine dx and dy
[ypeak,xpeak] = find(c==max(c(:)));
dx = xpeak-(wRef+wtmplt)+0.5*log(Ip(2,3)/Ip(2,1))/(log((Ip(2,2)*Ip(2,2))/(Ip(2,1)*Ip(2,3))));
dy = ypeak-(wRef+wtmplt)+0.5*log(Ip(3,2)/Ip(1,2))/(log((Ip(2,2)*Ip(2,2))/(Ip(1,2)*Ip(3,2))));

figure, hold on
h = imagesc(c)
%axis([200-50 200+50 250-50 250+50])

plot(xpeak,ypeak,'or')

% working on Gaussian interpolation
Nwidth = 1;
% if 1%(out(j,2)-Nwidth >0)&&(out(j,1)-Nwidth>0)&&(out(j,2)+Nwidth<Ny)&&(out(j,1)+Nwidth<Nx)
% cnt = cnt+1;

Ip = double(c(ypeak-Nwidth:ypeak+Nwidth,xpeak-Nwidth:xpeak+Nwidth));

xxx = xpeak + 0.5*log(Ip(2,3)/Ip(2,1))/(log((Ip(2,2)*Ip(2,2))/(Ip(2,1)*Ip(2,3))));
yyy = ypeak + 0.5*log(Ip(3,2)/Ip(1,2))/(log((Ip(2,2)*Ip(2,2))/(Ip(1,2)*Ip(3,2))));

figure
imagesc(Ip,[-0.5996    0.9874])


%%
clear X Y U V
ii = 0;
for ix = 1 : size(dr,1)
    for iy = 1 : size(dr,2)
        ii = ii +1;
        
        X(ii) = dr(ix,iy).x;
        Y(ii) = dr(ix,iy).y;
        U(ii) = dr(ix,iy).dx;
        V(ii) = dr(ix,iy).dy;
    end
end

figure
imagesc(imDiff), colormap gray, hold on
quiver(X,Y,U,V)


%%
irun = 0;
while(1)
    
[a,b] = max([dr.dr]);
if b<20 
    break
end
irun = irun + 1;
[row,col] = ind2sub(size(dr),b);


a1 = dr(row,col).dx;
a2 = dr(row,col).dy;
dr(row,col).dx = 0;
dr(row,col).dy = 0;
dr(row,col).dr = 0;
a5 = dr(row,col).dx;
a6 = dr(row,col).dy;

fprintf('run: %0.0f - before dx: %0.0f dy: %0.0f , after dx: %0.0f dy: %0.0f \n',irun,a1,a2,a5,a6)
end
%%
x = 482;
y = 334;

figure
title('origin figure')
imagesc(ImRef)
axis([x-20 x+20 y-20 y+20])
hold on
plot(x,y,'or')

figure
title('moved figure')
imagesc(im)
axis([x-20 x+20 y-20 y+20])
hold on
plot(x,y,'or')

tmplt = imcrop(im,[ix iy 2*wtmplt 2*wtmplt]);
% do correlation
c = normxcorr2(tmplt,ImRef);
% determine dx and dy
[ypeak,xpeak] = find(c==max(c(:)));
yoffSet = ypeak-size(tmplt,1);
xoffSet = xpeak-size(tmplt,2);
dx = xoffSet-ix+wtmplt;
dy = yoffSet-iy+wtmplt;
dx
dy


%%


%%


%%

%%
function [dx,dy,c] = fastCrossCorr(im,ImRef,ix,iy,wtmplt,wRef)
% 
%
% im,ImRef,ix,iy,wtmplt

% build sample image
tmplt    = imcrop(im,    [ix-wtmplt iy-wtmplt 2*wtmplt-1 2*wtmplt-1]);
ImRefLoc = imcrop(ImRef, [ix-1*wRef iy-1*wRef 2*wRef 2*wRef]);
% do correlation
c = normxcorr2(tmplt,ImRefLoc);
[ypeak,xpeak] = find(c==max(c(:)));
% Gaussian interpolation
Nwidth = 1;
Ip = double(c(ypeak-Nwidth:ypeak+Nwidth,xpeak-Nwidth:xpeak+Nwidth));
% determine dx and dy
dx = xpeak-(wRef+wtmplt)+0.5*log(Ip(2,3)/Ip(2,1))/(log((Ip(2,2)*Ip(2,2))/(Ip(2,1)*Ip(2,3))));
dy = ypeak-(wRef+wtmplt)+0.5*log(Ip(3,2)/Ip(1,2))/(log((Ip(2,2)*Ip(2,2))/(Ip(1,2)*Ip(3,2))));


end

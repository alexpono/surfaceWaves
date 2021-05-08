close all, clear all

%%


cd('C:\Users\Lenovo\Jottacloud\ENSEIGNEMENT\2020_2021\L3_TP_images\ondesDeSurface\imseqBeadDrop_cutrot')
listImages = dir('*.png');

ImRef = imread(listImages(62).name);

%%
% figure
% hi = imagesc(imDiff);
for it = 26 %1 : 62
    clear dr
im       = imread(listImages(it).name);
imDiff   = imsubtract(ImRef,im);

iix = 0;
iiy = 0;
% compute dr
for ix = 200 : 5 : 700  %200 : 10 : 700
    ix
    iix = iix + 1;
    iiy = 0;
    for iy = 250 : 5 : 350% 1000 %250 : 10 : 1000
        iiy = iiy + 1;
        % build sample image
        tmplt = imcrop(im,[ix-10 iy-10 19 19]);
        c = normxcorr2(tmplt,ImRef);
        [ypeak,xpeak] = find(c==max(c(:)));
        yoffSet = ypeak-size(tmplt,1);
        xoffSet = xpeak-size(tmplt,2);
        dr(iix,iiy).x  = ix;
        dr(iix,iiy).y  = iy;
        dr(iix,iiy).dx = xoffSet-ix+10;
        dr(iix,iiy).dy = yoffSet-iy+10;
    end
end


% hi.CData = imDiff;
% pause(.2)

end

fprintf('correlation finished \n')
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
for ix = 1 : size(dr,1)
    for iy = 1 : size(dr,2)
        ii = ii +1;
        dr(ix,iy).dr = sqrt((dr(ix,iy).dx)^2 + (dr(ix,iy).dy)^2);
    end
end
%%
[a,b] = max([dr.dr]);
a
b
[row,col] = ind2sub(size(dr),b)
%%
a1 = dr(row,col).dx;
a2 = dr(row,col).dy;
dr(row,col).dx = 0;
dr(row,col).dy = 0;
dr(row,col).dr = 0;
a5 = dr(row,col).dx;
a6 = dr(row,col).dy;

fprintf('before dx: %0.0f dy: %0.0f , after dx: %0.0f dy: %0.0f \n',a1,a2,a5,a6)
%%
function [ix,iy,xoffSet,yoffSet] = fastCrossCorr(im,ImRef,ix,iy)

end

%%
close all, clear all

%%


cd('C:\Users\Lenovo\Jottacloud\ENSEIGNEMENT\2020_2021\L3_TP_images\ondesDeSurface\imseqBeadDrop_cutrot')
listImages = dir('*.png');

ImRef = imread(listImages(62).name);

%%
% figure
% hi = imagesc(imDiff);
for it = 40 %1 : 62
    clear dr
im       = imread(listImages(it).name);
imDiff   = imsubtract(ImRef,im);

iix = 0;
iiy = 0;
% compute dr
for ix = 400 : 10 : 500  %200 : 10 : 700
    ix
    iix = iix + 1;
    iiy = 0;
    for iy = 500 : 10 : 600 %250 : 10 : 1000
        iiy = iiy + 1;
        % build sample image
        tmplt = imcrop(im,[ix-10 iy-10 19 19]);
        c = normxcorr2(tmplt,ImRef);
        [ypeak,xpeak] = find(c==max(c(:)));
        yoffSet = ypeak-size(tmplt,1);
        xoffSet = xpeak-size(tmplt,2);
        dr(iix,iiy).x  = ix;
        dr(iix,iiy).y  = iy;
        dr(iix,iiy).dx = xoffSet;
        dr(iix,iiy).dy = yoffSet;
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
quiver(X,Y,U,V)



%%




%%
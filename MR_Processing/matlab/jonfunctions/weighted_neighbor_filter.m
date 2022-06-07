%Messing around with a weighted-neighbor Gaussian filter
%clear all; close all;

%%
G3 = 1/16*[1 2 1; 2 4 2; 1 2 1];
U3 = [0 0 0; 0 1 0; 0 0 0];

%matrix to interpolate, e.g. displacement
%I = (sin(peaks))^2;
I = ux;%angle(ph_x);
%matrix used for weights, e.g. the anatomical or phase magnitude
%wij = peaks; 
wij = ph_mag;
wij = abs(wij)/max(abs(wij(:)));
t=0.1;
wij(wij>t)=t; wij = wij/t;

padsize = [floor(size(G3,1)/2) floor(size(G3,2)/2)];
wij = padarray(wij,padsize,'symmetric');
I = padarray(I,padsize,'symmetric');
W3 = cell(size(I));
I1 = zeros(size(I));

for i=(padsize(1)+1):(size(I,1)-padsize(1))
    for j=(padsize(2)+1):(size(I,2)-padsize(2))
%         if wij(i,j)==1
%             a=0;
%         end
        Iij = [I(i-1,j-1) I(i,j-1) I(i+1,j-1);...
            I(i-1,j) I(i,j) I(i+1,j);...
            I(i-1,j+1) I(i,j+1) I(i+1,j+1);];
        wm = [wij(i-1,j-1) wij(i,j-1) wij(i+1,j-1);...
            wij(i-1,j) wij(i,j) wij(i+1,j);...
            wij(i-1,j+1) wij(i,j+1) wij(i+1,j+1);];
        wm = wm/sum(wm(:));
        W3{i,j} = wij(i,j)*U3 + (1-wij(i,j))*G3.*wm;
        I1(i,j) = sum(sum(Iij'.*W3{i,j}));
    end
end

szI1 = size(I,1); szI2 = size(I,2);
I = I((padsize(1)+1):(szI1-padsize(1)),(padsize(2)+1):(szI2-padsize(2)));
I1 = I1((padsize(1)+1):(szI1-padsize(1)),(padsize(2)+1):(szI2-padsize(2)));
wij = wij((padsize(1)+1):(szI1-padsize(1)),(padsize(2)+1):(szI2-padsize(2)));


figure(11);
imagesc(wij)
figure(12)
imagesc(I);

figure(13);
imagesc(I1); colormap magma; axis image;
set(gca,'YDir','normal');

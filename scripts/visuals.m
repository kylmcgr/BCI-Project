% figure
% plot(0:80,accuracy)
% title('Square window sizes around 80,70 vs accuracy')
% xlabel('Radius of square (single pixel to whole image)')
% ylabel('Accuracy')

% figure
% hold on
% for pix=1:100
%     plot(fib,squeeze(accuracy(1,pix,:)))
% end
% title('Square window sizes vs accuracy')
% xlabel('Radius of square (single pixel to whole image)')
% ylabel('Accuracy')
% 
% figure
% shadedErrorBar(fib,squeeze(accuracy(1,:,:)),{@mean,@std}); 
% title('Square window sizes vs accuracy')
% xlabel('Radius of square (single pixel to whole image)')
% ylabel('Accuracy')

% load('Y:\Raw Data\20180308_S2R1\S2R1_searchlight_results.mat');
load('S25R1_searchlight_results.mat');
mask = zeros(size(pvalue_mask,1),size(pvalue_mask,2));
cz = 70; 
cx = 80;
% cz = 78; 
% cx = 68;
window_size = 20;
[x, z] = meshgrid(1:size(pvalue_mask,2), 1:size(pvalue_mask,1));
mask(pvalue_mask) = 5;
% mask((x - cx).^2 + (z - cz).^2 <= window_size.^2 & (x - cx).^2 + (z - cz).^2 > (window_size-1).^2) = 1;
mask(((x<cx-(window_size-1))|(x>cx+(window_size-1))|(z<cz-(window_size-1))|(z>cz+(window_size-1)))&((x>=cx-window_size)&(x<=cx+window_size)&(z>=cz-window_size)&(z<=cz+window_size))) = 3;
% window_size = 30;
% mask(((x<cx-(window_size-1))|(x>cx+(window_size-1))|(z<cz-(window_size-1))|(z>cz+(window_size-1)))&((x>=cx-window_size)&(x<=cx+window_size)&(z>=cz-window_size)&(z<=cz+window_size))) = 3;
cz2 = 2; 
cx2 = 26;
window_size = 20;
% mask((x - cx2).^2 + (z - cz2).^2 <= window_size.^2 & (x - cx2).^2 + (z - cz2).^2 > (window_size-1).^2) = 1;
mask(((x<cx2-(window_size-1))|(x>cx2+(window_size-1))|(z<cz2-(window_size-1))|(z>cz2+(window_size-1)))&((x>=cx2-window_size)&(x<=cx2+window_size)&(z>=cz2-window_size)&(z<=cz2+window_size))) = 2;
% window_size = 30;
% mask(((x<cx2-(window_size-1))|(x>cx2+(window_size-1))|(z<cz2-(window_size-1))|(z>cz2+(window_size-1)))&((x>=cx2-window_size)&(x<=cx2+window_size)&(z>=cz2-window_size)&(z<=cz2+window_size))) = 2;
mask(cz2-1:cz2+1,cx2-1:cx2+1) = 4;
mask(cz-1:cz+1,cx-1:cx+1) = 4;
mask(cz,cx) = 6;
mask(cz2-1,cx2) = 6;
imagesc(mask);
% 
% load('Z:\Acq_161502\UF');
% pixelsize = 0.1;
% X_img_mm = pixelsize/2 + (0:128-1)*pixelsize; % So that the center of the image is at true halfway point, not shifted.
% Z_img_mm = pixelsize/2 + (0:103-1)*pixelsize + UF.Depth(1);
% cmap = plotDuplexImage(X_img_mm, Z_img_mm, mask, zeros(103,128));


% load('S25R1_searchlight_results.mat');
% mask = zeros(103,128);
% cz = 70; 
% cx = 80;
% window_size = 20;
% [x, z] = meshgrid(1:128, 1:103);
% % mask(pvalue_mask) = 1;
% mask(((x<max(cx-(window_size-1),2))|(x>min(cx+(window_size-1),127))|(z<max(cz-(window_size-1),2))|(z>min(cz+(window_size-1),102)))&((x>=max(cx-window_size,1))&(x<=min(cx+window_size,128))&(z>=max(cz-window_size,1))&(z<=min(cz+window_size,103)))) = 3;
% window_size = 50;
% mask(((x<max(cx-(window_size-1),2))|(x>min(cx+(window_size-1),127))|(z<max(cz-(window_size-1),2))|(z>min(cz+(window_size-1),102)))&((x>=max(cx-window_size,1))&(x<=min(cx+window_size,128))&(z>=max(cz-window_size,1))&(z<=min(cz+window_size,103)))) = 3;
% window_size = 80;
% mask(((x<max(cx-(window_size-1),2))|(x>min(cx+(window_size-1),127))|(z<max(cz-(window_size-1),2))|(z>min(cz+(window_size-1),102)))&((x>=max(cx-window_size,1))&(x<=min(cx+window_size,128))&(z>=max(cz-window_size,1))&(z<=min(cz+window_size,103)))) = 3;
% mask(cz,cx) = 3;
% imagesc(mask);

% mask = zeros(103,128,100);
% cz = 70; 
% cx = 80;
% window_size = 15;
% [x, z,~] = meshgrid(1:128, 1:103,1:100);
% mask((x - cx).^2 + (z - cz).^2 <= window_size.^2 & (x - cx).^2 + (z - cz).^2 > (window_size-1).^2) = 1;

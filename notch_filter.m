clc;
clear;
close all;

img = imread('1.jpeg');
img = double(img);

[M,N,K] = size(img);

t = 0:pi/20:2*pi;
Xc = (M-N/2)/2;
Yc = (N+N/2)/2;
r_min = N/4;
r_max = N/2;

Xr_min = r_min * cos(t) + Xc;
Yr_min =  r_min * sin(t) + Yc;

Xr_max = r_max * cos(t) + Xc;
Yr_max =  r_max * sin(t) + Yc;

internal_mask = poly2mask(Xr_min,Yr_min,M,N);
filter_mask = poly2mask(Xr_max,Yr_max,M,N);
filter_mask(internal_mask) = 0;

res = zeros(size(img));

for ch = 1:size(img,3)
    color_channel = img(:,:,ch);
    
    Channel_FFT = fft2(color_channel);
    Channel_FFT = fftshift(Channel_FFT);
    
    figure(1);
    subplot(K,2,2*ch-1); imshow(log(1+abs(Channel_FFT)),[]); title(sprintf('channel %d',ch));
    
    Channel_FFT(filter_mask) = 0;
        
    subplot(K,2,2*ch);   imshow(log(1+abs(Channel_FFT)),[]); title(sprintf('channel %d',ch));
    
    res(:,:,ch) = ifft2(ifftshift(Channel_FFT));
end
res = uint8(abs(res));

figure,
subplot(1,2,1); imshow(uint8(img),[])
subplot(1,2,2); imshow(res,[])
imwrite(res,'res.bmp')


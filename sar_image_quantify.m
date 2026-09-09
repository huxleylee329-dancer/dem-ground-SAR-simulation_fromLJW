function [G_filter] = sar_image_quantify(im,win_size,sigma,displayflag)
% copyright @ youyanan 20160114
% input: im(complex sar image)
%        win_size(filter window size)
%        sigma( standard deviation sigma (positive),the default value for sigma is 0.5)
%        displayflag(display figure flag, no:0,yes:1 or nonzero number)
% output: G_ml_filter(image after multi-look and filter processing)
% -----------------------------------------------------------------------%
im_magnitude     = abs(im); %·ù¶ÈÍ¼Ïñ
% preprossing for filtering by three sigma
im_magnitude_std = std(im_magnitude(:));
im_magnitude(im_magnitude>=3*im_magnitude_std) = 3*im_magnitude_std;
G = (im_magnitude-min(im_magnitude(:)))./(max(im_magnitude(:))-min(im_magnitude(:))).*255;
% filter processing by using gaussian mode
J = mat2gray(G);
h = fspecial('gaussian',[win_size win_size],sigma);
G_filter = imfilter(J,h);
% output figure
if displayflag ~=0
    figure;imshow(G_filter)
end
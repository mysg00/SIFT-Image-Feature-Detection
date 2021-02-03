
clear; 
clc;

%read in the image
orig_img = imread("synthetic.jpg");
%resize the image
input_img = imresize(orig_img, 0.5);
%set any possible color to gray
input_img = rgb2gray(input_img);

%initialize count value for labeling
count = 1;
testcount = 0;

for j = 1: +1: 16
    %calculate the sigma values for Gaussian scale space layers j and j+1
    oldsigma = 2.^((j-1)/4);
    sigma = 2.^(j/4);
    %apply the gaussian filters to the input image
    filter1 = fspecial('gaussian', round(oldsigma * 3), oldsigma);
    filteredI1 = imfilter(input_img, filter1);
    filter2 = fspecial('gaussian', round(sigma * 3), sigma);
    filteredI2 = imfilter(input_img, filter2);
    %calculate the difference between the images
    filteredImage = filteredI2 - filteredI1;
    imageToGrey = mat2gray(filteredImage);
    %plot the difference of gaussians scale space
    subplot(4,4,count), imshow(imageToGrey), title("Scale space: " + count, "FontSize", 8);
    count = count + 1; 

end

% with the colorbar, the bottom last image will not properly appear

%{
%set up the colorbar by using manual position
top = get(subplot(4,4,4),'Position');
bottom = get(subplot(4,4,16),'Position');
colorbar('Position', [top(1)+top(3)+0.01 bottom(2) 0.05 0.8], 'LineWidth', 0.1)
caxis([0 256])
%}

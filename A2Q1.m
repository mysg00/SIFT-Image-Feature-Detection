clear; 
clc;

%read the image
original_img = imread("bird.jpg");
%resize the image to fit without necessary boundaries
input_img = imresize(original_img, 0.5);
%set the image to gray
input_img = rgb2gray(input_img);

%initialize values necessary for finding the 1st Gaussian scale space slice
sigma = 1;
m = 4;
hsize = round(sigma * 3);

%initialize count which will be used for labeling the subplots
count = 1;

%apply gaussian filter to the 1st scale space slice
filter = fspecial('gaussian', hsize, sigma);
filteredImage = imfilter(input_img, filter);

%loop through the 16 additional scale spaces
for j = 1: +1: 16
    %set the sigma and size of the filter
    sigma = 2.^(j/4);   
    hsize = round(sigma * 3);
    %apply the filter to the image
    filter = fspecial('gaussian', hsize, sigma);
    filteredImage = imfilter(input_img, filter);
    %plot the slices from 2 to 17
    subplot(4,4,count), imshow(filteredImage), title("Slice: " + (count+1) + ", Sigma " + sigma, "FontSize", 8);
    count = count + 1; 
end 

%set up the colorbar by using manual position
top = get(subplot(4,4,4),'Position');
bottom = get(subplot(4,4,16),'Position');
colorbar('Position', [top(1)+top(3)+0.01 bottom(2) 0.05 0.8], 'LineWidth', 0.1)







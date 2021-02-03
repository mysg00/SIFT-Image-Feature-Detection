clc;
clear;

%read the image
orig_img = imread("synthetic.jpg");
%resize the image for easier computation
input_img = imresize(orig_img, 0.5);
%set to gray
input_img = rgb2gray(input_img);

xFilter = [1 0 -1];
yFilter = [1; 0; -1];
%calculate the horizontal and vertical image gradients
Gx = conv2(input_img, xFilter, 'same');
Gy = conv2(input_img, yFilter, 'same');
[rows, cols] = size(input_img);

%Todo: confirm image gradient and second moment matrix calculations
%calculate matrices to be used in the second moment matrix
Ix = Gx.^2;
Iy = Gy.^2;
Ixy = Gx.*Gy;

%initialize count variable for labeling
count = 0;

%for loop for the slices the Gaussian scale space
for k = 0: +1: 16
   
    %compute and filter the gaussian scale space
    sigma = 2.^(k/4); 
    filter = fspecial('gaussian', round(3 * sigma), sigma);
    filteredImage = imfilter(input_img, filter);

    %compute and apply the gaussian filter to the second moment matrix
    filter2Matrix = fspecial('gaussian', round(6 * sigma), 2 * sigma);
    Ix = imfilter(Ix, filter2Matrix);
    Iy = imfilter(Iy, filter2Matrix);
    Ixy = imfilter(Ixy, filter2Matrix);

    %initialize a matrix to hold the harris-stevens operators
    hs = zeros(rows, cols);
    for i = 1:rows
        for j = 1:cols
            %compute the second moment matrix
            M = [sum(Ix(i,j), 'all') sum(Ixy(i,j), 'all');
                sum(Ixy(i,j), 'all') sum(Iy(i,j), 'all')];
            %calculate the harris-stevens operator
            hs(i,j) = det(M)- 0.01*trace(M);
        end
    end

    %Todo: choose better harris-stevens neighborhood operator
    %select the maximum harris-stevens operators within the neighborhood
    result = zeros(rows,cols);
    count1 = 0;
    for i = 2:rows-1
        for j = 2:cols-1
            %compares harris-stevens values
            if hs(i,j) > hs(i-1,j-1) && hs(i,j) > hs(i-1,j) && hs(i,j) > hs(i-1,j+1) && hs(i,j) > hs(i,j-1) && hs(i,j) > hs(i,j+1) && hs(i,j) > hs(i+1,j-1) && hs(i,j) > hs(i+1,j) && hs(i,j) > hs(i+1,j+1)
            result(i,j) = 1;
            count1 = count1 + 1;
            end
        end
    end

    %determine position of the maximum harris-stevens value
    hsV = zeros(count1,1);
    [rowPix,colPix] = find(result == 1);
    %place the position into hsV
    for i = 1:count1
        hsV(i) = hs(rowPix(i),colPix(i));
    end
    %sort the hsV in downwards direction
    [hsV,sortIndex] = sort(hsV, 'descend');
    rowPix = rowPix(sortIndex);
    colPix = colPix(sortIndex);
    %select the top 10 maximum harris-stevens operator
    top10 = round(0.1 * count1);
    pos = [colPix rowPix];
    top10Pos = pos(1:top10,:);
    %from the 2nd to the 17th gaussian scale space slices, insert a marker
    %for the top 10 harris-stevens values
    if (k > 0) 
        RGB = insertMarker(filteredImage, top10Pos, 'Size', 10);
        subplot(4,4,count), imshow(RGB), title("Slice: " + (count+1) + ", Sigma " + sigma, "FontSize", 8);
    end
    count = count + 1; 
    
end 

%set up the colorbar by using manual position
colormap('gray')
top = get(subplot(4,4,4),'Position');
bottom = get(subplot(4,4,16),'Position');
colorbar('Position', [top(1)+top(3)+0.01 bottom(2) 0.05 0.8], 'LineWidth', 0.1)
caxis([0 256])






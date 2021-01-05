clear; 
clc;

%read in the image and set to gray
input_img = imread("synthetic.jpg");
input_img = rgb2gray(input_img);
[rows,cols] = size(input_img);

%have a base array that will be used for concatenation
arrayHold = [0 0 0];

%initialize count variable
countTotal = 0;

%loop through scale spaces
for i = 2: +1: 15
    
    %compute the difference of gaussians scale space for the center layer
    oldsigmaCenter = 2.^((i-1)/4);
    sigmaCenter = 2.^(i/4);
    filter1C = fspecial('gaussian', round(oldsigmaCenter * 3), oldsigmaCenter);
    filteredI1C = imfilter(input_img, filter1C);
    filter2C = fspecial('gaussian', round(sigmaCenter * 3), sigmaCenter);
    filteredI2C = imfilter(input_img, filter2C);
    centerSlice = filteredI2C - filteredI1C;
    %calculate threshold boundary so that only pixels with 26 neighbors can
    %be analyzed
    thresholdImageBoundary = 2 * round(sigmaCenter);
	%rows and cols start at thresholdImageBoundary
    rowsEnd = rows - thresholdImageBoundary;
    colsEnd = cols - thresholdImageBoundary;

    %compute the difference of gaussians scale space for the bottom layer
    oldsigmaBot = 2.^((i-2)/4);
    sigmaBot = 2.^((i-1)/4);
    filter1B = fspecial('gaussian', round(oldsigmaBot * 3), oldsigmaBot);
    filteredI1B = imfilter(input_img, filter1B);
    filter2B = fspecial('gaussian', round(sigmaBot * 3), sigmaBot);
    filteredI2B = imfilter(input_img, filter2B);
    bottomSlice = filteredI2B - filteredI1B;

    %compute the difference of gaussians scale space for the top layer
    oldsigmaTop = 2.^((i)/4);
    sigmaTop = 2.^((i+1)/4);
    filter1T = fspecial('gaussian', round(oldsigmaTop * 3), oldsigmaTop);
    filteredI1T = imfilter(input_img, filter1T);
    filter2T = fspecial('gaussian', round(sigmaTop * 3), sigmaTop);
    filteredI2T = imfilter(input_img, filter2T);
    topSlice = filteredI2T - filteredI1T;
    
    %initialize an array to hold all potential candidate keypoints
    isCandidate = zeros(rows,cols);
    countCandidate = 0;
    thresholdDOG = 0;
    
    %loop through the image pixels within the image boundaries
    for j = thresholdImageBoundary:+1:rowsEnd
        for k = thresholdImageBoundary:+1:colsEnd
            
            %value is the center pixel point
            value = centerSlice(j,k);
            
            %the following are the areas that encompass the neighborhood
            areaCenter = centerSlice(j-1:j+1, k-1:k+1);
            areaTop = topSlice(j-1:j+1, k-1:k+1);
            areaBot = bottomSlice(j-1:j+1,k-1:k+1);
            
            %hessian matrix computation
            dxx = centerSlice(j, k+1) + centerSlice(j, k-1) - 2 * centerSlice(j,k);
            dyy = centerSlice(j+1,k) + centerSlice(j-1,k) - 2 * centerSlice(j,k);
            dxy = (centerSlice(j+1,k+1) + centerSlice(j-1,k-1) - centerSlice(j+1,k-1) - centerSlice(j-1,k+1))/4;
            trace = dxx + dyy;
            deter = dxx * dyy - dxy * dxy;
            hessThreshold = 11;
           
            %tests if the center pixel is the maximum among its neighbors
            %and if it fits within the hessian threshold
            if (value < min(areaTop(:)) && value < min(areaBot(:)) && value > (thresholdDOG + max(areaCenter(1,:))) && value > (thresholdDOG + max(areaCenter(3,:))) && value > (thresholdDOG + max(areaCenter(:,1))) && value > (thresholdDOG + max(areaCenter(:,3)))&& ((trace^2/deter) < (hessThreshold)) && deter > 0)
                %if the center pixel satisfies the condition, then it
                %becomes a candidate point
                isCandidate(j,k) = 1;
                countCandidate = countCandidate + 1;
            %tests if the center pixel is the minimum among its neighbors
            %and if it fits within the hessian threshold
            elseif (value < min(areaTop(:)) && value < min(areaBot(:)) && value < (min(areaCenter(1,:)) - thresholdDOG) && value < (min(areaCenter(3,:)) - thresholdDOG) && value < (min(areaCenter(:,1)) - thresholdDOG) && value < (min(areaCenter(:,3)) - thresholdDOG) && ((trace^2/deter) < (hessThreshold)) && deter > 0)
                %if the center pixel satisfies the condition, then it
                %becomes a candidate point
                isCandidate(j,k) = 1;
                countCandidate = countCandidate + 1;
            else
                isCandidate(j,k) = 0;
            end
        end
    end
    %scaleValue keeps track of the sigma value corresponding to each
    %keypoint
    scaleValue = sigmaCenter .* (ones(countCandidate,1));
    %rowPix and colPix keep track of the position value corresponding to
    %each keypoint
    [rowPix,colPix] = find(isCandidate == 1);
    %arrayHold contains the position and scale value for each candidate
    %point
    arrayHold = cat(1,arrayHold,[rowPix colPix scaleValue]);
    countTotal = countTotal + countCandidate;
end

%separating the arrayHold into its columns: row of the pixel, column of the
%pixel and sigma value 
arrayHold(1,:) = [];
arrayRow = arrayHold(:,1);
arrayCol = arrayHold(:,2);
arrayRad = arrayHold(:,3);

%show the image with the formed circle keypoints
imshow(input_img); 
hold on;
colorbar;
viscircles([arrayCol arrayRow], arrayRad ,'Color','b');
hold off;







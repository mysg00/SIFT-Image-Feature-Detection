clear; 
clc;

%read in the image
input_img = imread('synthetic.jpg');
%resize is used for question 7
%input_img = imresize(input_img, 1.2);
%set image to gray
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
            if (value > max(areaTop(:)) && value > max(areaBot(:)) && value > (thresholdDOG + max(areaCenter(1,:))) && value > (thresholdDOG + max(areaCenter(3,:))) && value > (thresholdDOG + max(areaCenter(:,1))) && value > (thresholdDOG + max(areaCenter(:,3))) && ((trace^2/deter) < (hessThreshold)) && deter > 0)
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
[numRows, numCols] = size(arrayHold);

%horizontal and vertical filters
xFilter = [1.0 0.0 -1.0];
yFilter = [1.0; 0.0; -1.0];

%number of bins
num_bins = 36;

%size of each bin
hist_step = 2*pi/num_bins;
%
hist_orient=[0:hist_step:2*pi-hist_step];
%initialize base array for concatenation
arrNewhold = [0,0,0,0];

input_img = double(input_img); 
for i = 1:numRows
    %apply the gaussian filter to the image based on the sigma value
    %corresponding to the keypoint
    sigma = arrayHold(i,3);
    filter = fspecial('gaussian', round(sigma * 3), sigma);
    filteredImage = imfilter(input_img, filter);
    
    %identify the image horizontal and vertical gradients
    Gx = conv2(filteredImage, xFilter, 'same');
    Gy = conv2(filteredImage, yFilter, 'same');
    
    %create the gaussian filter for applying to the magnitude
    filterGauss = fspecial('gaussian', round(sigma * 1.5), 1.5);
    %compute gradient magnitude
    gradientMagnitude = sqrt(Gy.^2 + Gx.^2);
    %compute gradient orientation
    gradOrientation = atan2(Gy,Gx);
    gradientMagnitude = conv2(gradientMagnitude, filterGauss, 'same');
    %compute shift for the area of the neighborhood
    shift = floor(length(filter)/2);
    %set r to the rows of the keypoints
    r = arrayHold(i,1);
    %set c to the columns of the keypoints
    c = arrayHold(i,2);
    %find the window that will be used for the weighting
    wght = gradientMagnitude((r-shift):(r+shift),(c -shift):(c+shift));
    %find the window that will be used for the gradient orientation
    %identifications
     grad_window = gradOrientation((r-shift):(r+shift),(c -shift):(c+shift));
     %initialize array to hold the orientation
     orient_hist=zeros(length(hist_orient),1);
     %loop for all 36 bins
     for bin=1:length(hist_orient)
        %deals with the negative values obtained from the tan values
        diff = mod( grad_window - hist_orient(bin), 2*pi );
        %add the weighted gradient magnitudes to the respective bins
        orient_hist(bin) = orient_hist(bin)+sum(sum(wght.*max(1 - abs(diff)/hist_step,0)));
     end
     %identify the max value in the histogram
     maxValue = max(orient_hist);
     %identify the bin id of all bins with peaks that are above 80% of the
     %max value
     [idbin,a]=find(orient_hist>0.8*maxValue);
     %concatenate the values to the arrNewhold
     arrNewhold = cat(1,arrNewhold,[r*(ones(length(idbin),1)) c*(ones(length(idbin),1)) sigma*(ones(length(idbin),1)) idbin]);
     
        %used for creating the histogram
        %{
        if(i == 18) 
           disp(r)
           disp(c)
           xAll = 1:36;
           xValues = xAll';
           figure;
           plot(xValues, orient_hist, '-x');
            xlabel('Bins');
            ylabel('Sum of weighted gradient magnitudes');
            title('Orientation Histogram');
        end
        %}
     
end


%remove the initialize row of zeros
arrNewhold(1,:) = [];
%arrNewhold contains rows, cols, sigma, idbin of the keypoints in
%respective order
for i=1:size(arrNewhold,1)
    r = arrNewhold(i,1);
    c = arrNewhold(i,2);
    sigma = arrNewhold(i,3);
    idbin = arrNewhold(i,4); 
end

%add the circles of the keypoints to the image
colormap('gray');
imagesc(input_img); 
hold on;
colorbar;
viscircles([arrayCol arrayRow], arrayRad ,'Color','b');

%plot the lines for all keypoints
sizeNewArray = size(arrNewhold);
for j = 1:sizeNewArray
    x = arrNewhold(j,2);
    y = arrNewhold(j,1);
    length = arrNewhold(j,3);
    bin = arrNewhold(j,4);
    %as the position plot is done with 2 values, calculate the x and y
    %coordinates of the second points
    y2=y+((2*length)*sin((bin/36)*2*pi-0.175));
    x2=x+((2*length)*cos((bin/36)*2*pi-0.175));
    %plot the line
    plot([x x2],[y y2], 'Color', 'r', 'LineWidth', 2)
end
hold off;



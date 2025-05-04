close all; clear all;
global mask mask2 incidence1 elevation n m cells nbors data melt nsmelt prop frac means incidence2 incidence3 incidence4 incidence5

greenwhite(1,1:3) = [0 0.2 0];
greenwhite(2,1:3) = [1 1 1];

%% 
% read data and set up two masks for the domain extent
mask = imread('Binary_Mask.tif');  
mask2= mask;
mask2(mask2<0) = NaN;   %substitute -3.4028235e+38 with NaN
mask(mask<0) = 0;       %substitute -3.4028235e+38 with 0

cells = sum(sum(mask)); %number of pixels in the region of interest (ROI)
boundarylength=(sum(sum(abs(diff(mask))))+sum(sum(abs(diff(mask'))))); %length of domain boundary
[m n] = size(mask); %rows and columns of mask

%calculate the number of neighbouring (domain) patches for each (domain) patch 
nbors = zeros(m,n);
nbors(2:m-1,2:n-1)=mask(1:m-2,2:n-1)+mask(3:m,2:n-1)+mask(2:m-1,1:n-2)+mask(2:m-1,3:n);
nbors(1,2:n-1)=mask(2,2:n-1)+mask(1,1:n-2)+mask(1,3:n);
nbors(m,2:n-1)=mask(m-1,2:n-1)+mask(m,1:n-2)+mask(m,3:n);
nbors(2:m-1,1)=mask(1:m-2,1)+mask(3:m,1)+mask(2:m-1,2);
nbors(2:m-1,n)=mask(1:m-2,n)+mask(3:m,n)+mask(2:m-1,n-1);
nbors(1,1)=mask(1,2)+mask(2,1);
nbors(1,n)=mask(1,n-1)+mask(2,n);
nbors(m,1)=mask(m-1,1)+mask(m,2);
nbors(m,n)=mask(m-1,n)+mask(m,n-1);


%-----------------------------Adding SAT as controlling variables-----%
I1 = imread('SAT_4D_20220401_Pad.tif');
hmax1 = max(max(I1));
hmin1 = min(min(I1));
%[r c val] = find(isnan(I5)==false); 
incidence1 = -1*ones(m,n); 
for i = 1:m
    for j = 1:n
        if (isnan(I1(i,j)) == false)     
            incidence1(i,j) = (hmax1-I1(i,j))/(hmax1-hmin1);
        end
    end
end

I2 = imread('SAT_4D_20220502_Pad.tif');
hmax2 = max(max(I2));
hmin2 = min(min(I2));
%[r c val] = find(isnan(I5)==false); 
incidence2 = -1*ones(m,n); 
for i = 1:m
    for j = 1:n
        if (isnan(I2(i,j)) == false)     
            incidence2(i,j) = (hmax2-I2(i,j))/(hmax2-hmin2);
        end
    end
end

I3 = imread('SAT_4D_20220526_Pad.tif');
hmax3 = max(max(I3));
hmin3 = min(min(I3));
%[r c val] = find(isnan(I5)==false); 
incidence3 = -1*ones(m,n); 
for i = 1:m
    for j = 1:n
        if (isnan(I3(i,j)) == false)     
            incidence3(i,j) = (hmax3-I3(i,j))/(hmax3-hmin3);
        end
    end
end
I4 = imread('SAT_4D_20220628_Pad.tif');
hmax4 = max(max(I4));
hmin4 = min(min(I4));
%[r c val] = find(isnan(I5)==false); 
incidence4 = -1*ones(m,n); 
for i = 1:m
    for j = 1:n
        if (isnan(I4(i,j)) == false)     
            incidence4(i,j) = (hmax4-I4(i,j))/(hmax4-hmin4);
        end
    end
end

I5 = imread('SAT_4D_20220722_Pad.tif');
hmax5 = max(max(I5));
hmin5 = min(min(I5));
%[r c val] = find(isnan(I5)==false); 
incidence5 = -1*ones(m,n); 
for i = 1:m
    for j = 1:n
        if (isnan(I5(i,j)) == false)     
            incidence5(i,j) = (hmax5-I5(i,j))/(hmax5-hmin5);
        end
    end
end






% read and covert elevation data
b = imread('DEM_Pad.tif'); 
b(b<0) = NaN;
hmax = max(max(b)); %elevation raster maximum
hmin = min(min(b)); %elevation raster minimum
elevation = -1*ones(m,n);
for i = 1:m
    for j = 1:n
        if (isnan(b(i,j)) == false)     
            elevation(i,j) = (hmax-b(i,j))/(hmax-hmin);
        end
    end
end

% calculate domain mean values for elevation and incidence angle
esum = 0;
[ii jj] = find(incidence1>=0);
for i = 1:length(ii)
    esum = esum+incidence1(ii(i),jj(i));
end
means(1) = esum/length(ii);

%second mean values of incidence
esum = 0;
[ii jj] = find(incidence2>=0);
for i = 1:length(ii)
    esum = esum+incidence2(ii(i),jj(i));
end
means(2) = esum/length(ii);


%third mean values of incidence
esum = 0;
[ii jj] = find(incidence3>=0);
for i = 1:length(ii)
    esum = esum+incidence3(ii(i),jj(i));
end
means(3) = esum/length(ii);

%fourth mean values of incidence
esum = 0;
[ii jj] = find(incidence4>=0);
for i = 1:length(ii)
    esum = esum+incidence4(ii(i),jj(i));
end
means(4) = esum/length(ii);

%fifth mean values of incidence
esum = 0;
[ii jj] = find(incidence5>=0);
for i = 1:length(ii)
    esum = esum+incidence5(ii(i),jj(i));
end
means(5) = esum/length(ii);


hsum = 0;
[ii jj] = find(elevation>=0);
for i = 1:length(ii)
    hsum = hsum+elevation(ii(i),jj(i));
end
means(6) = hsum/length(ii);

% clear all temporary variables
clear a b c ii jj i j r read val hmax hmin hsum esum

%% 
% load in the snowcover masks for fitting
data(:,:,1)= imread('Snow_3_Class_20220401_Pad.tif'); 
data(:,:,2)= imread('Snow_3_Class_20220502_Pad.tif');
data(:,:,3)= imread('Snow_3_Class_20220526_Pad.tif'); 
data(:,:,4)= imread('Snow_3_Class_20220629_Pad.tif'); 
data(:,:,5)= imread('Snow_3_Class_20220704_Pad.tif');




dimdata = size(data);  
if (length(dimdata)<3) 
    dimdata = 1; 
else
    dimdata = dimdata(:,3); % n0layers
end

% measure fraction of snow cover for each snow mask.
for i = 1:dimdata
    prop(i)= sum(sum(data(:,:,i)));    
    l1 = find(abs(diff(data(:,:,i).*mask2))==1);  
    l2 = find(abs(diff((data(:,:,i).*mask2)'))==1);    
    frac(i)=(length(l1)+length(l2))/cells;
end
%%%%%%%% can impose equation based on the NDSI Values %%%%%%%
%% 
% [rho alpha beta gamms pp qq rr]
lb = [2 0 0 0 0 0 0]; 
ub = [10 9 9 9 3 3 3];

% LHS : Latin Hypercube Sampling
nsamples = 10;
nruns = 1;
samplepars = lhsdesign(nsamples,length(lb));
drawon = true;

xin = (lb+(ub-lb).*samplepars(8,:));
% % If the user wants to run automatically the un-comment the line below, else
% % run the script of Stochastic_proxy.m seperately


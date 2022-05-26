% Code for PCA
close all;
nframes = 10;
[video, ~] = mmread('video2.mp4',1:nframes);
% mmread function used for reading the video file
for r=1:nframes
    video1(:,:,r) = rgb2gray(video.frames(r).cdata);
end
[m,n,~] = size(video1);
I = video1;
nI = video1;
I1 = video1;
% I1 is the original video
nI = imnoise(nI,'gaussian',0,0.01);
nI = imnoise(nI,'poisson');
nI = imnoise(nI,'salt & pepper');
% nI is the noisy video
for i=1:nframes
    I(:,:,i) = adapmedfilt(nI(:,:,i), 11);
end
% I is the median filtered output
%% Code for generating patches
m1 = m/4-1;
n1 = n/4-1;
patches = zeros(64,m1,n1,nframes);
% patches is the matrix of overlapping patches for each image and each
% frame
noisypixels = (I~=nI);
% missing pixels that had different values due to noise
% 4*4 sample interval
for r = 1:nframes
    for i=1:m1
        for j=1:n1
            patches(:,i,j,r) = reshape(I(1+(i-1)*4:4*(i+1),1+4*(j-1):4*(j+1),r),[64, 1]);
        end
    end
end
%% Code for patchwise denoising using pacthmatcher and LRMC
Omega = zeros(64,m1,n1);
% Omega is the set of know values i.e. 1 where pixel value is know and 0
% else
for i=1:m1
    for j=1:n1
        Omega(:,i,j) = reshape(noisypixels(1+(i-1)*4:4*(i+1),1+4*(j-1):4*(j+1),1),[64, 1]);
    end
end
Omega = ~Omega;
num_match = 5;
% number of matching patches per frame
denoised = zeros(64,num_match*nframes,m1,n1);
for i=1:m1
    for j=1:n1
        index = [i j 1];
        % denoising frame 1
        [mindicesi,mindicesj,mindicesf] = patchmatcher(patches,index,num_match);
        mindicesi = reshape(mindicesi,[num_match*nframes 1]);
        mindicesj = reshape(mindicesj,[num_match*nframes 1]);
        mindicesf = reshape(mindicesf,[num_match*nframes 1]);
%         i,j and f matching indices for a given patch with index = [i j 1]
        matchedpatches = zeros(64,num_match*nframes);
        for k=1:nframes*5
            matchedpatches(:,k) = reshape(patches(:,mindicesi(k),mindicesj(k),mindicesf(k)),[64 1]); 
        end
        P = matchedpatches;
        C=6;
        [Q] = pca1(P,C);
%         Completing matrix P using PCA algoritm
        denoised(:,:,i,j) = Q;
    end
end
%%
avgmatchedpatches = zeros(64,m1,n1);
% matrix for average matched patches
for i=1:m1
    for j=1:n1
        for k=1:num_match*nframes
            avgmatchedpatches(:,i,j) = denoised(:,k,i,j) + avgmatchedpatches(:,i,j);
        end
        avgmatchedpatches(:,i,j) = avgmatchedpatches(:,i,j)/(num_match*nframes); 
    end
end
% averaging matched patches
denoisedimage = zeros(m,n);
nums = zeros(m,n);
for i=1:m1
    for j=1:n1
        denoisedimage(1+(i-1)*4:4*(i+1),1+4*(j-1):4*(j+1)) = denoisedimage(1+(i-1)*4:4*(i+1),1+4*(j-1):4*(j+1)) + reshape(avgmatchedpatches(:,i,j),[8,8]);
        nums(1+(i-1)*4:4*(i+1),1+4*(j-1):4*(j+1)) = nums(1+(i-1)*4:4*(i+1),1+4*(j-1):4*(j+1)) + 1;
    end
end
% reshaping matched patches to form images
denoisedimage = denoisedimage./nums;
% final denoised average image
%% Results
figure(1);
imshow(I1(:,:,1));
title('Original Frame');
figure(2);
imshow(nI(:,:,1));
title('Noisy Frame');
figure(3);
imshow(I(:,:,1));
title('Median-Filtered Frame');
figure(4);
imshow(denoisedimage/255);
title('Reconstructed Frame');
rmse = immse(double(I1(:,:,1)),denoisedimage)/(sum(I1(:,:,1).^2,"all"));
disp(rmse);
% rmse between reconstructed frame and original frame
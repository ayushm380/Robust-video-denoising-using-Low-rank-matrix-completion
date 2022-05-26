function [mindicesi,mindicesj,mindicesf] = patchmatcher(patches,index,num_match)
% num_match -> number of matching patches per frame
% patches [8*8 m1 n1 nframes] patches
% index 3*1-> vector of the required patch index
[~, m1, n1, nframes] = size(patches);
mindicesi = [];
mindicesj = [];
mindicesf = [];
% matching i,j,frameno indices
reqpatch = patches(:,index(1),index(2),index(3));
% the patch to be matched
for f=1:nframes
    MAD = sum(abs(reqpatch-patches(:,:,:,f)));
%     L1 norm
    MAD = reshape(MAD,[m1*n1 1]);
    [~,matchindices] = mink(MAD,num_match);
%     first num_match locations per frame based on L1 norm
    [i,j] = ind2sub([m1 n1],matchindices);
    i = reshape(i,[1 num_match]);
    j = reshape(j,[1 num_match]);
    f1 = repelem(f,num_match);
    mindicesi = [mindicesi i];
    mindicesj = [mindicesj j];
    mindicesf = [mindicesf f1];
end
end
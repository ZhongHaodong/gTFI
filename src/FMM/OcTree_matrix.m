function ocTree = OcTree_matrix(mask, centroid, voxel_size, depth)
% OcTree_matrix:
% Generate an Octree based on the data matrix and assign corresponding centroids of triangular faces.
% Inputs:
%   mask       - Logical mask defining the region of interest.
%   centroid   - Points (centroids) within the volume.
%   voxel_size - Size of each voxel.
%   depth      - Depth of the octree.
% Output:
%   ocTree     - Structure representing the octree.
%
%   Created by Sven Holcombe.
%   1.0     - 2013-03 Initial release
%   1.1     - 2013-03 Added shrinking bins and allocate/deallocate space
%
%  Modified by Haodong Zhong on 2024.11.26.

% Initialize the structure
ocTree = struct();
ocTree.Centroid = centroid;  % Points renamed to Centroid
ocTree.PointBins = ones(size(centroid, 1), 1, 'single');
ocTree.BinCount = 1;
ocTree.BinIndex_Range = single([[1, 1, 1] size(mask)]);
ocTree.BinDepths = single(0);
ocTree.Mask = logical(mask); % Renamed to Mask
ocTree.BinParents = zeros(0, 1);
ocTree.voxel_size = voxel_size; % Renamed to VoxelSize
ocTree.Depth = depth; % Renamed to Depth

% Initialize bin parents
ocTree.BinParents(1) = 0;

% Start dividing bins
ocTree = preallocateSpace(ocTree);
ocTree = divide(ocTree, 1);

% Nested helper functions

    function ocTree = preallocateSpace(ocTree)
        ocTree.BinDepths(:) = single(0);
        ocTree.BinParents(:) = single(0);
    end

    function ocTree = divide(ocTree, startingBins)
        % Divide bins recursively
        for i = 1:length(startingBins)
            binNo = startingBins(i);
            % Prevent dividing beyond a minimum size
            thisBinIndex_Range = ocTree.BinIndex_Range(binNo, :);
            binIndexSize = diff(thisBinIndex_Range([1:3; 4:6])) + 1;
            maxbinIndexSize = max(binIndexSize);
            this_mask = ocTree.Mask(...
                thisBinIndex_Range(1):thisBinIndex_Range(4), ...
                thisBinIndex_Range(2):thisBinIndex_Range(5), ...
                thisBinIndex_Range(3):thisBinIndex_Range(6));

            oldCount = ocTree.BinCount;
            if maxbinIndexSize > max(size(ocTree.Mask)) / 2^ocTree.Depth && ...
                    (nnz(ocTree.PointBins == binNo) || any(this_mask, 'all'))
                ocTree = divideBin(ocTree, binNo);
                ocTree = divide(ocTree, oldCount + 1:ocTree.BinCount);
            end
        end
    end

    function ocTree = divideBin(ocTree, binNo)
        % Subdivide a bin into smaller bins
        binPtMask = ocTree.PointBins == binNo;
        thisBinsPoints = ocTree.Centroid(binPtMask, :);

        % Calculate new bin boundaries
        oldMin = ocTree.BinIndex_Range(binNo, 1:3);
        oldMax = ocTree.BinIndex_Range(binNo, 4:6);
        newDiv = mean([oldMin; oldMax], 1);

        % Build the new boundaries of 8 subdivisions
        minMidMax = [oldMin floor(newDiv) ceil(newDiv) oldMax];
        newBounds = minMidMax([...
            1 2 3 4 5 6;
            1 2 9 4 5 12;
            1 8 3 4 11 6;
            1 8 9 4 11 12;
            7 2 3 10 5 6;
            7 2 9 10 5 12;
            7 8 3 10 11 6;
            7 8 9 10 11 12]);

        % Assign points to new bins
        binMap = cat(3, [0 0 0], [0 0 1], [0 1 0], [0 1 1], ...
            [1 0 0], [1 0 1], [1 1 0], [1 1 1]);
        gtMask = bsxfun(@gt, thisBinsPoints, ...
            (newDiv - size(ocTree.Mask) ./ 2 - 1) .* ocTree.voxel_size');
        [~, binAssignment] = max(all(bsxfun(@eq, gtMask, binMap), 2), [], 3);

        % Update structure with new bins
        newBinInds = ocTree.BinCount + 1:ocTree.BinCount + 8;
        ocTree.BinIndex_Range(newBinInds, :) = single(newBounds);
        ocTree.BinDepths(newBinInds) = single(ocTree.BinDepths(binNo) + 1);
        ocTree.BinParents(newBinInds) = single(binNo);
        ocTree.PointBins(binPtMask) = single(newBinInds(binAssignment));
        ocTree.BinCount = single(ocTree.BinCount + 8);
    end

end

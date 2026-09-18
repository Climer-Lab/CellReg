function sp=footprints_to_sparse(spatial_footprints)
% This function converts the spatial footprints of one session into a sparse
% representation that allows corr2-equivalent spatial correlations between
% cells to be computed from a single sparse dot product (see sparse_corr2).
% Since each footprint covers only a small fraction of the FOV, this is
% orders of magnitude cheaper than corr2 on the full frames.

% Inputs:
% 1. spatial_footprints - either a matrix of size (number_of_cells x
%    y_size x x_size), or a path / matrix accepted by get_spatial_footprints

% Outputs:
% 1. sp - struct with fields:
%    S    - sparse (number_of_pixels x number_of_cells) matrix, one footprint per column
%    P    - number of pixels in the FOV
%    n    - number of cells
%    sum1 - (number_of_cells x 1) sum of each footprint
%    sum2 - (number_of_cells x 1) sum of squares of each footprint

if ~isnumeric(spatial_footprints)
    footprint_info=get_spatial_footprints(spatial_footprints);
    footprint_info=footprint_info.load_footprints;
    spatial_footprints=footprint_info.footprints;
end

number_of_cells=size(spatial_footprints,1);
% reshape to (cells x pixels) is free because cells are the first dimension;
% transposing gives one column per cell so that column slicing is cheap.
S=sparse(reshape(spatial_footprints,number_of_cells,[])).';

sp=struct;
sp.S=S;
sp.P=size(S,1);
sp.n=number_of_cells;
sp.sum1=full(sum(S,1)).';
sp.sum2=full(sum(S.^2,1)).';

end

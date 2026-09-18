function [patch,rows,cols]=transform_footprint_in_disc(image,image_row_offset,image_col_offset,N,M,map_coordinates,disc_center,centroid_location,maximal_radius)
% This function evaluates a spatial transformation (given as an inverse
% mapping from output pixels to input coordinates) for the output pixels
% that map to within maximal_radius of a cell's centroid, and returns only
% the bounding box of that disc. It is the shared vectorized core of
% rotate_spatial_footprint, translate_spatial_footprint and rotate_cell:
% pixels outside the disc are zero by definition, so nothing else has to be
% computed.

% Inputs:
% 1. image - the (whole or cropped) input image
% 2. image_row_offset - row index in the full frame of image(1,1) minus 1 (0 for a whole image)
% 3. image_col_offset - column index in the full frame of image(1,1) minus 1
% 4. N - number of rows of the full frame
% 5. M - number of columns of the full frame
% 6. map_coordinates - function handle mapping (2 x K) output pixel
%    coordinates to (2 x K) input coordinates
% 7. disc_center - (2 x 1) output-frame coordinates of the centroid's image
%    (used only to bound the disc; the mask itself is exact)
% 8. centroid_location - (1 x 2) centroid in input coordinates
% 9. maximal_radius - in pixels

% Outputs:
% 1. patch - transformed values for rows x cols of the output frame
% 2. rows - row indexes (in the full frame) covered by patch
% 3. cols - column indexes covered by patch

if any(~isfinite(disc_center))
    % a NaN centroid (empty footprint) never fails the radius test in the
    % original per-pixel code, so the whole frame is transformed
    rows=1:N;
    cols=1:M;
else
    rows=max(1,floor(disc_center(1)-maximal_radius)-1):min(N,ceil(disc_center(1)+maximal_radius)+1);
    cols=max(1,floor(disc_center(2)-maximal_radius)-1):min(M,ceil(disc_center(2)+maximal_radius)+1);
end
patch=zeros(numel(rows),numel(cols));
if isempty(rows) || isempty(cols)
    return
end

[p,q]=ndgrid(rows,cols);
wanted_coords=map_coordinates([p(:)' ; q(:)']);
distance_from_centroid=sqrt(sum((wanted_coords-centroid_location').^2,1));
inside_disc=~(distance_from_centroid>maximal_radius);
if ~any(inside_disc)
    return
end
wanted_coords=wanted_coords(:,inside_disc);
wanted_coords(wanted_coords<0)=0;
wanted_coords(1,wanted_coords(1,:)>N+1)=N+1;
patch(inside_disc)=interp_zero_outside(image,wanted_coords(1,:)-image_row_offset,wanted_coords(2,:)-image_col_offset);

end

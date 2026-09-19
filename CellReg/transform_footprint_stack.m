function [transformed_spatial_footprints]=transform_footprint_stack(spatial_footprints,centroid_locations,theta,translations,center_of_FOV,microns_per_pixel,show_progress)
% This function rotates (by theta, about center_of_FOV) and translates every
% spatial footprint of a session, producing exactly what a per-cell loop of
% rotate_spatial_footprint (theta ~= 0) or translate_spatial_footprint
% (theta == 0) applied to the transposed footprints would produce, but it
% only reads and writes the small window around each cell instead of the
% full frame.

% Inputs:
% 1. spatial_footprints - (number_of_cells x y_size x x_size)
% 2. centroid_locations - (number_of_cells x 2), (x,y) per cell
% 3. theta - rotation in degrees (0 = translation only)
% 4. translations - [y_translation x_translation] as passed by align_images
% 5. center_of_FOV - axis of rotation, (y,x) as in align_images
% 6. microns_per_pixel
% 7. show_progress - true to display a progress bar

% Outputs:
% 1. transformed_spatial_footprints - same size as spatial_footprints

if nargin<7
    show_progress=false;
end

[number_of_cells,y_size,x_size]=size(spatial_footprints);
transformed_spatial_footprints=zeros(number_of_cells,y_size,x_size);

% the per-cell functions work on the transposed footprint (x_size x y_size),
% so in their frame N=x_size and M=y_size and rows are x coordinates:
N=x_size;
M=y_size;
a=translations(2);
b=translations(1);
if theta~=0
    maximal_cell_radius=30; % in microns (rotate_spatial_footprint)
    theta_rad=theta*pi/180;
    transformation=[cos(theta_rad) -sin(theta_rad) ; sin(theta_rad) cos(theta_rad)]';
    trans_inv=transformation^-1;
    center=[center_of_FOV(1) ; center_of_FOV(2)];
    map_coordinates=@(coords) trans_inv*(coords-center+[a ; b])+center;
else
    maximal_cell_radius=25; % in microns (translate_spatial_footprint)
    map_coordinates=@(coords) coords+[a ; b];
end
normalized_maximal_radius=maximal_cell_radius/microns_per_pixel;
crop_margin=ceil(normalized_maximal_radius)+2; % bilinear neighbours + rounding

for k=1:number_of_cells
    if show_progress
        progress_tick(k,number_of_cells)
    end
    centroid_location=centroid_locations(k,:);
    if theta~=0
        disc_center=transformation*(centroid_location(:)-center)+center-[a ; b];
    else
        disc_center=centroid_location(:)-[a ; b];
    end
    if all(isfinite(centroid_location))
        % input window (in the transposed frame: rows = x, cols = y)
        in_rows=max(1,floor(centroid_location(1))-crop_margin):min(N,ceil(centroid_location(1))+crop_margin);
        in_cols=max(1,floor(centroid_location(2))-crop_margin):min(M,ceil(centroid_location(2))+crop_margin);
    else
        in_rows=1:N;
        in_cols=1:M;
    end
    % transposed crop: I(p,q) = spatial_footprints(k,q,p)
    crop=reshape(spatial_footprints(k,in_cols,in_rows),numel(in_cols),numel(in_rows))';
    [patch,rows,cols]=transform_footprint_in_disc(crop,in_rows(1)-1,in_cols(1)-1,N,M,map_coordinates,disc_center,centroid_location,normalized_maximal_radius);
    if ~isempty(patch)
        transformed_spatial_footprints(k,cols,rows)=reshape(patch',1,numel(cols),numel(rows));
    end
end

end

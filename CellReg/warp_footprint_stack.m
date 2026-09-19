function [warped_spatial_footprints]=warp_footprint_stack(spatial_footprints,displacement_field,show_progress)
% This function applies a non-rigid displacement field (from imregdemons) to
% every spatial footprint of a session. It is equivalent to calling
% imwarp(footprint,displacement_field) on each full-frame footprint, but each
% cell is only warped inside a window around its support (the displacement
% field is bounded, so everything further away stays zero), which avoids a
% full-frame warp and a strided full-frame slice per cell.

% Inputs:
% 1. spatial_footprints - (number_of_cells x y_size x x_size)
% 2. displacement_field - (y_size x x_size x 2) as returned by imregdemons
% 3. show_progress - true to display a progress bar

% Outputs:
% 1. warped_spatial_footprints - same size as spatial_footprints

if nargin<3
    show_progress=false;
end

[number_of_cells,y_size,x_size]=size(spatial_footprints);
warped_spatial_footprints=zeros(number_of_cells,y_size,x_size);

% pixels can move at most this far, so a cell's warped support lies within
% its original support expanded by the margin (and the sources of those
% output pixels within twice the margin):
maximal_displacement=ceil(max(abs(displacement_field(:))))+1;
window_margin=2*maximal_displacement+2;

% support of every cell along y and x (two passes over the data):
rows_with_signal=any(spatial_footprints~=0,3);
cols_with_signal=any(spatial_footprints~=0,2);

for k=1:number_of_cells
    if show_progress
        progress_tick(k,number_of_cells)
    end
    rows=find(rows_with_signal(k,:));
    cols=find(cols_with_signal(k,1,:));
    if isempty(rows) || isempty(cols)
        continue % empty footprint stays empty
    end
    rows=max(1,rows(1)-window_margin):min(y_size,rows(end)+window_margin);
    cols=max(1,cols(1)-window_margin):min(x_size,cols(end)+window_margin);
    crop=reshape(spatial_footprints(k,rows,cols),numel(rows),numel(cols));
    warped_crop=imwarp(crop,displacement_field(rows,cols,:));
    warped_spatial_footprints(k,rows,cols)=reshape(warped_crop,1,numel(rows),numel(cols));
end

end

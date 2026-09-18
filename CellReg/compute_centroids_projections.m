function [centroid_projections]=compute_centroids_projections(centroid_locations,spatial_footprints)
% This function projects the centroid locations of all the cells onto the FOV

% Inputs:
% 1. centroid_locations
% 2. spatial_footprints

% Outputs:
% 1. centroid_projections

number_of_sessions=size(centroid_locations,2);

centroid_projections=cell(1,number_of_sessions);
for n=1:number_of_sessions
    this_session_centroids=centroid_locations{n};
    this_session_footprint_info = get_spatial_footprints(spatial_footprints{n});
    footprints_size=this_session_footprint_info.size;
    centroid_projections{n}=project_centroids(this_session_centroids,footprints_size(2),footprints_size(3));
end

end

function projection=project_centroids(centroids,y_size,x_size)
% Each cell contributes a 3x3 weighted stamp (or a single pixel near the
% border) around its rounded centroid. The stamps are accumulated with a
% sparse constructor instead of allocating a full (cells x y x x) array and
% summing over cells; the weights are dyadic so the sum is exact.
number_of_cells=size(centroids,1);
rows=zeros(9*number_of_cells,1);
cols=zeros(9*number_of_cells,1);
vals=zeros(9*number_of_cells,1);
count=0;
stamp_dy=[-1 -1 -1 0 0 0 1 1 1];
stamp_dx=[-1 0 1 -1 0 1 -1 0 1];
stamp_w=[1/4 1/2 1/4 1/2 1 1/2 1/4 1/2 1/4];
for k=1:number_of_cells
    y=round(centroids(k,2));
    x=round(centroids(k,1));
    if y>1.5 && x>1.5 && y<y_size-1 && x<x_size-1
        rows(count+1:count+9)=y+stamp_dy;
        cols(count+1:count+9)=x+stamp_dx;
        vals(count+1:count+9)=stamp_w;
        count=count+9;
    elseif y>0 && x>0 && y<y_size && x<x_size
        count=count+1;
        rows(count)=y;
        cols(count)=x;
        vals(count)=1;
    end
end
projection=full(sparse(rows(1:count),cols(1:count),vals(1:count),y_size,x_size));
end

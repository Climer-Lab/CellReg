function [rotated_spatial_footprint]=rotate_spatial_footprint(original_spatial_footprint,theta,translations,center_of_FOV,centroid_location,microns_per_pixel)
% This function rotates a spatial footprint according to the
% angle theta, and translates the footprint in (x,y) according to translations.

% Inputs:
% 1. original_spatial_footprint
% 2. theta - rotation angle in degrees
% 3. translations
% 4. center_of_FOV - axis of rotation
% 5. centroid_location - is used to rotate only values in a certain radius
% 6. microns_per_pixel

% Outputs:
% 1. rotated_spatial_footprint

maximal_cell_radius=30; % in microns

theta=theta*pi/180;
N=size(original_spatial_footprint,1);
M=size(original_spatial_footprint,2);
a=translations(2);
b=translations(1);
transformation=[cos(theta) -sin(theta) ; sin(theta) cos(theta)]';
trans_inv=transformation^-1;
normalized_maximal_radius=maximal_cell_radius/microns_per_pixel;

% only pixels within the maximal radius of the centroid are non-zero, so
% only that disc is evaluated (vectorized version of the per-pixel loop).
% The disc of output pixels is centered on the forward image of the centroid:
center=[center_of_FOV(1) ; center_of_FOV(2)];
map_coordinates=@(coords) trans_inv*(coords-center+[a ; b])+center;
disc_center=transformation*(centroid_location(:)-center)+center-[a ; b];
rotated_spatial_footprint=zeros(N,M);
[patch,rows,cols]=transform_footprint_in_disc(original_spatial_footprint,0,0,N,M,map_coordinates,disc_center,centroid_location,normalized_maximal_radius);
rotated_spatial_footprint(rows,cols)=patch;

end

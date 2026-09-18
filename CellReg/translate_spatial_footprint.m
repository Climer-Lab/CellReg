function [translated_spatial_footprint]=translate_spatial_footprint(original_spatial_footprint,translations,centroid_location,microns_per_pixel)
% This function translates a spatial footprint according to the provided
% translations.

% Inputs:
% 1. original_spatial_footprint
% 2. translations
% 3. centroid_location - is used to translate only values in a certain radius
% 4. microns_per_pixel

% Outputs:
% 1.translated_spatial_footprint

maximal_cell_radius=25; % in microns

N=size(original_spatial_footprint,1);
M=size(original_spatial_footprint,2);
a=translations(2);
b=translations(1);
normalized_maximal_radius=maximal_cell_radius/microns_per_pixel;

% only pixels within the maximal radius of the centroid are non-zero, so
% only that disc is evaluated (vectorized version of the per-pixel loop):
map_coordinates=@(coords) coords+[a ; b];
disc_center=centroid_location(:)-[a ; b];
translated_spatial_footprint=zeros(N,M);
[patch,rows,cols]=transform_footprint_in_disc(original_spatial_footprint,0,0,N,M,map_coordinates,disc_center,centroid_location,normalized_maximal_radius);
translated_spatial_footprint(rows,cols)=patch;

end

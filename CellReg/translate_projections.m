function [translated_image]=translate_projections(original_image,translations)
% This function translates an image in (x,y) according to provided translations.

% Inputs:
% 1. original_image
% 2. translations

% Outputs:
% 1. translated_image

N=size(original_image,1);
M=size(original_image,2);
a=translations(2);
b=translations(1);

% inverse-mapping every output pixel at once (vectorized version of the
% original per-pixel loop):
[p,q]=ndgrid(1:N,1:M);
wanted_coords=[p(:)'+a ; q(:)'+b];
wanted_coords(wanted_coords<0)=0;
wanted_coords(1,wanted_coords(1,:)>N+1)=N+1;
translated_image=reshape(interp_zero_outside(original_image,wanted_coords(1,:),wanted_coords(2,:)),N,M);

end

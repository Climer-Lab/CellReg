function values=interp_zero_outside(image,rows,cols)
% This function computes bilinearly interpolated pixel values of an image at
% (possibly non-integer) coordinates, treating everything outside the image
% as zero. It is a vectorized equivalent of interpolate_pixel_value: the
% image is padded with one ring of zeros so that interp2 reproduces the same
% partial-weight behaviour at the borders (e.g. coordinates between 0 and 1,
% or between N and N+1), and the extrapolation value 0 covers everything
% further out.

% Inputs:
% 1. image
% 2. rows - row coordinates of the requested pixels (any shape)
% 3. cols - column coordinates of the requested pixels (same shape as rows)

% Outputs:
% 1. values - interpolated values, same shape as rows

[N,M]=size(image);
padded_image=zeros(N+2,M+2,'like',image);
padded_image(2:N+1,2:M+1)=image;
values=interp2(padded_image,cols+1,rows+1,'linear',0);

end

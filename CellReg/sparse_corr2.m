function correlations=sparse_corr2(sp_a,cell_a,sp_b,cells_b)
% This function computes the spatial correlation (equivalent to corr2 on the
% full-frame footprints) between one cell of session a and several cells of
% session b, using the sparse representation from footprints_to_sparse.
%
% corr2(A,B) = sum((A-mean(A)).*(B-mean(B))) / sqrt(sum((A-mean(A)).^2)*sum((B-mean(B)).^2))
%            = (A.B - P*mA*mB) / sqrt((A.A - P*mA^2)*(B.B - P*mB^2))
% where P is the number of pixels, so only the dot products over the
% (sparse) footprint supports are needed.

% Inputs:
% 1. sp_a - sparse footprints of session a (footprints_to_sparse)
% 2. cell_a - index of the cell in session a
% 3. sp_b - sparse footprints of session b
% 4. cells_b - indexes of the cells in session b

% Outputs:
% 1. correlations - row vector, one correlation per cell in cells_b

P=sp_a.P;
cells_b=cells_b(:);
dot_products=full((sp_a.S(:,cell_a).'*sp_b.S(:,cells_b)));
mean_a=sp_a.sum1(cell_a)/P;
mean_b=(sp_b.sum1(cells_b)/P).';
numerator=dot_products-P*mean_a*mean_b;
variance_a=sp_a.sum2(cell_a)-P*mean_a^2;
variance_b=(sp_b.sum2(cells_b)).'-P*mean_b.^2;
correlations=numerator./sqrt(variance_a*variance_b);

end

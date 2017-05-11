%{
Getter functions which return the transmission function of a correction lens.
 This is just a phase shift as the light propagates through the lens.
---------------------------------------------------------------------------
Input:

k = complex wavenumber.

r_matrix  = matrix with all the distances from origin to a (x,y) point in
the current plane.

obj = An object which contains all the information of the lens, such as the
focal length and the power.
---------------------------------------------------------------------------
Output:

TF = the transmission function of a correction lens which can be put just
before the cornea to adjust a near-sight or long-sight eye.

%}
function TF = TF_correct_lens(k,r_matrix)
f_correction = -0.47; % a positive value means farsightedness
TF = exp(-1i*k*r_matrix.^2/(2*f_correction));
end
%{
Getter functions which return the transmission function of a lens. This is
just a phase shift as the light propagates through the lens.
---------------------------------------------------------------------------
Input:

k = complex wavenumber.

r_matrix  = matrix with all the distances from origin to a (x,y) point in
the current plane.

obj = An object which contains all the information of the lens, such as the
focal length and the power.
---------------------------------------------------------------------------
Output:

TF = the transmission function of a lens.

%}
function TF = get_TF(k,r_matrix,obj)
            TF=exp(-1i*k*r_matrix.^2/(2*obj.focal));
end
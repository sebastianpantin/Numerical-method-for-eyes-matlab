%{
Getter functions which return the transmission function of a cylindric lens
that simulates astigmatism. This lens will get a rough edge by just varying
it in one directiion, in this case the y-direction.
---------------------------------------------------------------------------
Input:

ymat = A matrix with all the y sample points.
---------------------------------------------------------------------------
Output:

TF = the transmission function of a cylindric lens
that simulates astigmatism. which can be put just
before the cornea.
 --------------------------------------------------------------------------
Properties:

f_cyl = focal length of the cylindric lens. This changes the degree of
astigmatism.

k = wavenumber.
%}
function [ TF ] = TF_astigatism(ymat)
f_cyl=1/2;
k = (2*pi*1.0003)/(550e-9);
TF = exp(-1i*k*ymat.^2/(2*f_cyl));

end
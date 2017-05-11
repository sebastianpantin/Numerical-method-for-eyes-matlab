%{
BPM is a function based on the beam propagation method which calculates the
electromagnetic field for small changes in positition with the help of
angular spectrum propagation in the eye. For the first part we multiply the
corneas transmission function with the field right before the cornea to get
the resulting field after the cornea. The second part does the same thing
but for the lens and the last part propagates the wave through the aqueous
and vitreous humor. 

---------------Input---------------------
E1 = The electromagnetic field right before the next propagation,
delta_z = the length of propagation in meters,
delta_image = Sample distance at the image plane,
wavelength = the wavelength of the wave,
length_position = the current position of the total length of propagation
Lens = A class that constructs the cornea and the lens and contains focal
length etc,
T_aperatur = Transmission function of the pupil.

--------------Output---------------------
E2 = the resulting E-field after the cornea, the lens and at the retina.
%}

function E2=BPM_Gullstrand(E1,delta_z,delta_image,wavelength, length_position, Lens, TF_pupil)
if length_position == 0
    E2 = Lens(1).TF.*E1;

elseif length_position > 0 && length_position < 0.50e-3

    E2 = E1;
    
elseif length_position >= 0.50e-3 && length_position < 3.6e-3
    
    E2 = Angular_propagation(E1,delta_z,delta_image,wavelength,1.3374);
    
elseif length_position == 3.6e-3
    
    E2 = Lens(2).TF.*E1.*TF_pupil;
    
elseif length_position  > 3.6e-3 && length_position <= 4.1e-3
   
    E2 = E1;
    
elseif length_position == 4.1e-3
    
    E2 = Lens(3).TF.*E1;
    
elseif length_position  > 4.1e-3 && length_position <= 6.6e-3
   
     E2 = E1;

elseif length_position == 6.6e-3
    
    E2 = Lens(4).TF.*E1;
    
elseif length_position  > 6.6e-3 && length_position <= 7.2e-3
   
    E2 = E1;
    
else
    
    E2 = Angular_propagation(E1,delta_z,delta_image,wavelength,1.336);
    
end

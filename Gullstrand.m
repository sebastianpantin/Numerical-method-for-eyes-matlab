%%
close all
clc

global N
load Intensity_matrix

N=2048;     % Matrix size
delta_object=100e-2/N; % Sampling distance in the object plane

wavelength=550e-9;
refractive_index_air=1;

% Different method, uses radius of curvature.
%Lens(1) = Lenses(7.259e-3, -5.585e-3, 1, 1.376, 0);
%Lens(2) = Lenses(8.672e-3, -6.328e-3, 1.336, 1.406, 4.979e-3);

% Current method, uses optical power.
Lens(1) = Lenses(48.8, -5.88, 1, 1.376, 0.50e-3);
Lens(2) = Lenses(5.00, 2.53, 1.336, 1.386, 0.5e-3);
Lens(3) = Lenses(2.53, 3.47, 1.386, 1.406, 2.5e-3);
Lens(4) = Lenses(3.47, 8.33, 1.406, 1.386, 0.6e-3);

% Focal length
Lens(1).focal = 1/get_power(Lens(1));
Lens(2).focal = 1/get_power(Lens(2));
Lens(3).focal = 1/get_power(Lens(3));
Lens(4).focal = 1/get_power(Lens(4));

p(1)=-L; p(2)=3.1e-3;
for i=1:length(p)
    if i==1
        p(i)=p(i);
    else
        p(i)=p(i)-q(i-1);
    end
    q(i)=1./(1./Lens(i).focal-1./p(i));
    m(i)=(-1).*q(i)/p(i);
end

M=-1*prod(m);

delta_image=M.*delta_object;

x_vector=-N/2*delta_image:delta_image:(N/2-1)*delta_image;
y_vector=x_vector;
[x_matrix,y_matrix]=meshgrid(x_vector,y_vector); 
r_matrix=sqrt(x_matrix.^2+y_matrix.^2);

% Wave number
k0=2*pi*refractive_index_air/wavelength;
r_k=sqrt(L^2-(x_matrix).^2-(y_matrix).^2);
E_in=exp(1i*k0*r_k)./r_k;
f_cyl=1/0.75;

% Constants for the lenses/parts of the eye
Pupill_diameter=4e-3;
T_apertur=r_matrix<(Pupill_diameter/2);
k(1)=(2*pi*1.376)/(632.8e-9/1.376);
k(2)=2*pi*1.386/(632.8e-9/1.386);
k(3)=2*pi*1.406/(632.8e-9/1.406);
k(4)=2*pi*1.386/(632.8e-9/1.386);
for i=1:length(Lens)
    Lens(i).TF=get_TF(k(i),r_matrix,Lens(i));
end

% Calculation of the electrical field
f_korrektion =-0.47; % a positive value means farsightedness
% Focal length help
T_lins_korrektion=exp(-1i*k0*r_matrix.^2/(2*f_korrektion));
L=24e-3;
delta_z=0.1e-3;
Lvekt=0:delta_z:L;
E1=E_in;
I_norm=zeros(N,length(Lvekt));
steg_nummer=0;
for akt_L=Lvekt
    steg_nummer=steg_nummer+1;
    
    E2=BPM_Gullstrand(E1,delta_z,delta_image,wavelength, akt_L, Lens, T_apertur);
    
    I2=abs(E2).^2; 
    
    E2_y=E2(:,N/2+1);
    I2_y=abs(E2_y).^2;
    I2_y_norm=I2_y/max(I2_y);
    I_norm(:,steg_nummer)=I2_y_norm;
        
    E1=E2;
end
%%
figure(11)
imagesc(Lvekt*1e3,y_vector*1e3,I_norm/max(max(I_norm))*64)
title('       Intensity distribution along the propagation direction', 'Fontsize',14)
set(gca,'FontSize',14)
yticks([-2 -1 0 1 2])
xticks([0 6 12 18 24])
xlabel('z [mm]')
ylabel('y [mm]')
colorbar
colormap('jet');
%pause

%%
PSF=abs(E2).^2;
PSF1=fft2c(B).*fft2c(PSF);


I=ifft2c(fft2c(B).*fft2c(PSF));

figure('Name','Plane 2, Retina','NumberTitle','off')
%image(xvekt*1e3*delta_b,yvekt*1e3*delta_b,PSF/max(max(PSF))*64)
image(x_vector*1e3*delta_image,y_vector*1e3*delta_image,I/max(max(I))*64)
colormap(gray)
figure('Name','Plane 1, Eye Chart','NumberTitle','off')
image(x_vector*1e3*delta_image,y_vector*1e3*delta_image,B/max(max(B))*64)
colormap(gray)

%% plot intensity, phase and mag
figure(1)
plot(x_vector,PSF(N/2+1,:)/max(PSF(N/2+1,:)));
axis([-1e-3 1e-3 0 1])
set(gca,'FontSize',14)
yticks([0 0.5 1])
xticks([-1e-3 0 1e-3])
xlabel('x (m)'); ylabel('Intensity');
title(['Intensity distribution in the image plane']);

figure(2)
%plot obs field mag
plot(x_vector,abs(E2(N/2+1,:)));
xlabel('x (m)'); ylabel('Magnitude');
title(['Magnitude in the image plane']);

figure(3)
%plot obs field phase
plot(x_vector,unwrap(angle(E2(N/2+1,:))));
xlabel('x (m)'); ylabel('Fas (rad)');
title(['E-field phase distribution in the image plane']);

%% Calculate a number of different errors
err1 = immse(I/max(max(I))*64, B/max(max(B))*64);
fprintf('\n The mean-squared error is %0.4f\n', err1);
K = imabsdiff(I/max(max(I))*64, B/max(max(B))*64);
[ssimval, ssimmap] = ssim(I/max(max(I))*64,B/max(max(B))*64);
fprintf('\n The SSIM value is %0.4f.\n',ssimval);
[peaksnr, snr] = psnr(I/max(max(I))*64,B/max(max(B))*64);
fprintf('\n The Peak-SNR value is %0.4f \n', peaksnr);
fprintf('\n The SNR value is %0.4f \n', snr);
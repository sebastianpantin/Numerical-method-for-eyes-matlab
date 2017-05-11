%%
clc; clear all;
format short
global N
load Intensity_matrix


N=2048;     %matrisstorlek
delta_o=30e-2/N; %samplinsavst�nd i objektplan

x_k=0;  %PSF fr�n origo!
y_k=0;  %OBS 100mm = 0.1 !!!
L=5;    %propargationsstr�cka 5m
lambda_noll=550e-9;
n_medium=1;
%vanlig lins

%Konstanter f�r hornhinnan om den ses som en "converging menicus"
%Kallar p� klassen Lenses
%Lens(1) = Lenses(7.259e-3, -5.585e-3, 1, 1.376, 0);
Lens(1) = Lenses(48.35, -6.11, 1, 1.3771, 0.5e-3);
%Konstanter f�r linsen, tar medelv�rde f�r brytningsindex f�r linsen
%Lens(2) = Lenses(8.672e-3, -6.328e-3, 1.336, 1.406, 4.979e-3);
Lens(2) = Lenses(8.10, 14.00, 1.336, 1.42, 3.6e-3);
%R�knar ut fokall�ngd
Lens(1).focal = 1/get_power(Lens(1))+0.01;
Lens(2).focal = 1/get_power(Lens(2));

%Lens(1).f = Fok_lens(Lens(1));
%Lens(2).f = Fok_lens(Lens(2));
%Avst�nd i systemet
p(1)=-5; p(2)=3.1e-3;
for i=1:2
    if i==1
        p(i)=p(i);
    else
        p(i)=p(i)-q(i-1);
    end
    q(i)=1./(1./Lens(i).focal-1./p(i));
    m(i)=(-1).*q(i)/p(i);
end
M=-1*m(1).*m(2);
delta_b=M.*delta_o;  %samplingsavstånd i bildplan

xvekt=-N/2*delta_b:delta_b:(N/2-1)*delta_b; %objektplanet samma som Ein planet
yvekt=xvekt;
[xmat,ymat]=meshgrid(xvekt,yvekt); 
rmat=sqrt(xmat.^2+ymat.^2);

%V�gvektorn k
k0=2*pi*n_medium/lambda_noll;
r_k=sqrt(L^2+(x_k-xmat).^2+(y_k-ymat).^2);
E_in=exp(1i*k0*r_k)./r_k;
f_cyl=1/0.75;

%Konstanter f�r linser/�gats delar
D_apertur=4e-3;     %pupilldiameter
T_apertur=rmat<(D_apertur/2);
k(1)=(2*pi*1.3771)/(632.8e-9/1.3771);
k(2)=2*pi*1.42/(632.8e-9/1.42);
for i=1:length(Lens)
    Lens(i).TF=get_TF(k(i),rmat,Lens(i));
end

%Ber�kning av elf�lt
f_korrektion =-0.47; %+ betyder långsynt
%fokallangd  hjalpmedel
T_lins_korrektion=exp(-1i*k0*rmat.^2/(2*f_korrektion));
f_cyl=0.5;
E_ut1=Lens(1).TF.*E_in.*TF_astigatism(xmat,f_cyl,k0);%.*T_lins_korrektion; %.*TF_astigatism(xmat,f_cyl,k);
E_ut2=Angular_propagation(E_ut1,3.1e-3,delta_b,632.8e-9/1.336,1.336);
E_ut3=Lens(2).TF.*E_ut2.*T_apertur;

%Ber�kna E_f�ltet efter propagation till näthinnan
E_ut4=Angular_propagation(E_ut3,16.7e-3,delta_b,632.8e-9/1.336,1.336);
PSF=abs(E_ut4).^2;
PSF1=fft2c(B).*fft2c(PSF);


I=ifft2c(fft2c(B).*fft2c(PSF));

figure('Name','Plan 2, n�thinnan','NumberTitle','off')
%image(xvekt*1e3*delta_b,yvekt*1e3*delta_b,PSF/max(max(PSF))*64)
image(xvekt*1e3*delta_b,yvekt*1e3*delta_b,I/max(max(I))*64)
colormap(gray) 
figure('Name','Plan 1, syntavlan','NumberTitle','off')
image(xvekt*1e3*delta_b,yvekt*1e3*delta_b,B/max(max(B))*64)
colormap(gray)

%% plot intensity, phase and mag
figure(1)
plot(xvekt,PSF(N/2+1,:));
xlabel('x (m)'); ylabel('Intensitet');
title(['Propagationssträcka= ',num2str(L),' m']);

figure(2)
%plot obs field mag
plot(xvekt,abs(E_ut4(N/2+1,:)));
xlabel('x (m)'); ylabel('Magnitud');
title(['Propagationssträcka= ',num2str(L),' m']);

figure(3)
%plot obs field phase
plot(xvekt,unwrap(angle(E_ut4(N/2+1,:))));
xlabel('x (m)'); ylabel('Fas (rad)');
title(['Propagationssträcka= ',num2str(L),' m']);

%% get errors 
err1 = immse(I/max(max(I))*64, B/max(max(B))*64);
fprintf('\n The mean-squared error is %0.4f\n', err1);
K = imabsdiff(I/max(max(I))*64, B/max(max(B))*64);
[ssimval, ssimmap] = ssim(I/max(max(I))*64,B/max(max(B))*64);
fprintf('\n The SSIM value is %0.4f.\n',ssimval);
[peaksnr, snr] = psnr(I/max(max(I))*64,B/max(max(B))*64);
fprintf('\n The Peak-SNR value is %0.4f \n', peaksnr);
fprintf('\n The SNR value is %0.4f \n', snr);
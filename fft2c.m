function output=fft2c(x);
output=fftshift(fft2(fftshift(x)));
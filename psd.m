%*************************************************************************
%                COMUNICACIONES II(CURSO 2012/2013)
%*************************************************************************
%            Carlos Medina Rodr�guez & Alberto Olivares Vicente
%                       {cmedina,aolivares}@.ugr.es
%*************************************************************************
% function psd(y,Fs,Nfft)
% y    : Vector de muestras de se�al
% Fs   : Frecuencia de muestreo (Hz) 
% Nfft : Longitud de ventana para el c�lculo de la FFT (Nfft<=length(y))
%*************************************************************************
function [Pyy,f]=psd(y,Fs,Nfft)

% Esta funci�n permite obtener la densidad de potencia espectral 
% haciendo uso de la funci�n 'pwelch' 
[Pyy,f]=pwelch(y,hamming(Nfft),[],Nfft,Fs,'onesided');

% En este caso 'pwelch' calcula la parte positiva de frecuencias multiplicada 
% por 2, por tanto hay que tener en cuenta ese factor en los valores de psd
% devueltos:
Pyy=Pyy*0.5; %[W/Hz]

end


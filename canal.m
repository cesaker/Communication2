%*************************************************************************
%                COMUNICACIONES II (CURSO 2014/2015)
%*************************************************************************
%          Jos� Carlos Segura, Carlos Medina, �ngel de la Torre
%                   {segura, cmedina, atv}@.ugr.es
%*************************************************************************
% function [y]=canal(y,B,Fs,SNR)
%
% Argumentos de entrada:
% - y   : se�al a la entrada del canal
% - B   : Tasa de bits por segundo
% - Fs  : Frecuencia de muestreo seleccionada
% - SNR : SNR de se�al por muestra 
%*************************************************************************
% Salida:
% - yfr : se�al a la salida del canal
%*************************************************************************
function [yfr]=canal(y,B,Fs,SNR)

%% Filtro paso-baja de fase cero (corrige el desfase de se�al) 
N=100;                 % Orden del filtro FIR
fnorm=B./(Fs/2);       % Frecuencia de corte normalizada
c=fir1(N,fnorm);       % Vector de coeficientes del filtro FIR

% Modificaci�n para permitir filtfilt
ly=length(y);
if(ly<=3*length(c))
    y(3*length(c))=0;
end

yf=filtfilt(c,1,y);   
yf=yf(1:ly); 

%% A�adimos ruido:
yfr=awgn(yf,SNR,'measured');
%yfr=yf+randn(size(yf)).*sigma;

end
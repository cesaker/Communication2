%*************************************************************************
%                COMUNICACIONES II(CURSO 2012/2013)
%*************************************************************************
%          Jos� Carlos Segura, Carlos Medina, �ngel de la Torre
%                   {segura, cmedina, atv}@.ugr.es
%*************************************************************************
%% Par�metros de la funci�n:
%*****************************************************************
% s      : muestras de se�al a la entrada del filtro
% Fs     : Frecuencia de muestreo seleccionada (Hz)
% Rb     : Tasa binaria (bits/s)
% Nb     : N�mero de bits a decodificar
%*****************************************************************
%% Salidas:
% s      : se�al a la entrada del filtro adaptado (truncada a Nb bits)
% sbo    : Secuencia binaria recuperada
%*****************************************************************
function [s,sbo]=filtro_adaptado_plantilla_final(s,Fs,Rb,Nb)

% N�mero de muestras por s�mbolo
Nss=Fs/Rb;  

% Enventanado de la se�al para decodificar s�lo los Nb primeros bits
s=s(1:Nss*Nb);

% Definici�n del vector de tiempos
dt=1/Fs;
tv=[1:length(s)]*dt;

%% Definici�n del pulso para los s�mbolos "1" y "0" seg�n el c�digo de l�nea empleado.
%% C�digo polar_nrz (pulso rectangular): 
% (---Completar)
p = ones(length(Nss)); % (omitiendo signo de amplitud de los s�mbolos)
s1=p;
s2=-p;

% Dise�o del filtro adaptado:
% (---Completar)
filtro_adp=flipud(p);

% Filtrado de la se�al usando el filtro adaptado:
% (---Completar)
fo=filter(filtro_adp,1,s);

% Representaci�n
figure;
subplot(2,1,1);
plot(tv,s);
title('Se�al de entrada al filtro');
xlabel('t(s)');
grid on;

subplot(2,1,2);
plot(tv,fo);
title('Salida del filtro adaptado');
xlabel('t(s)');

%% Determinaci�n de los instantes �ptimos de muestreo 
% Considere la descripci�n del algoritmo propuesto en el gui�n 
% de la pr�ctica. Obtenga un vector "Ts" con los instantes muestrales
% que permita la correcta detecci�n de los s�mbolos.
% ---(Completar)

iinicial = Nss;
n=1;
for i=iinicial:1:iinicial+Nss-1
       j=1;
   k=1;
   v0=0;
   v1=0;
   
   
   for aux=i:1:length(s)/Nss;
       
      if s(aux)<0
          v1(j)=s(aux);
          j=j+1;
      else
          v0(k)=s(aux);
          k=k+1;
      end
   
   end
   
   media0=mean(v0);
   media1=mean(v1);
   Deltan(n)=media1-media0;
   n=n+1;
end

[valor ioptimo]=max(Deltan);
Ts=[ioptimo:Nss:length(s)];
hold on
plot(tv(Ts),fo(Ts),'or');
hold off
grid on

if nargin==6
    saveas(gcf,ruta_imagen);
    saveas(gcf,ruta_imagen2);
end

% Representamos con c�rculos rojos los instantes de muestreo �ptimos
hold on
plot(tv(Ts),fo(Ts),'or');
hold off
grid on

%% Recuperaci�n de la secuencia binaria a partir de los s�mbolos. Tenga en cuenta el 
%% tipo de c�digo de l�nea utilizado para fijar el umbral (en este caso, polar-NRZ) 
% ---(Completar)
sbo= int64(s(Ts)>0); 


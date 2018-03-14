%*************************************************************************
%                COMUNICACIONES II(CURSO 2012/2013)
%*************************************************************************
%          José Carlos Segura, Carlos Medina, Ángel de la Torre
%                   {segura, cmedina, atv}@.ugr.es
%*************************************************************************
%% Parámetros de la función:
%*****************************************************************
% s      : muestras de señal a la entrada del filtro
% Fs     : Frecuencia de muestreo seleccionada (Hz)
% Rb     : Tasa binaria (bits/s)
% Nb     : Número de bits a decodificar
%*****************************************************************
%% Salidas:
% s      : señal a la entrada del filtro adaptado (truncada a Nb bits)
% sbo    : Secuencia binaria recuperada
%*****************************************************************
function [s,sbo]=filtro_adaptado_plantilla_final(s,Fs,Rb,Nb)

% Número de muestras por símbolo
Nss=Fs/Rb;  

% Enventanado de la señal para decodificar sólo los Nb primeros bits
s=s(1:Nss*Nb);

% Definición del vector de tiempos
dt=1/Fs;
tv=[1:length(s)]*dt;

%% Definición del pulso para los símbolos "1" y "0" según el código de línea empleado.
%% Código polar_nrz (pulso rectangular): 
% (---Completar)
p = ones(length(Nss)); % (omitiendo signo de amplitud de los símbolos)
s1=p;
s2=-p;

% Diseño del filtro adaptado:
% (---Completar)
filtro_adp=flipud(p);

% Filtrado de la señal usando el filtro adaptado:
% (---Completar)
fo=filter(filtro_adp,1,s);

% Representación
figure;
subplot(2,1,1);
plot(tv,s);
title('Señal de entrada al filtro');
xlabel('t(s)');
grid on;

subplot(2,1,2);
plot(tv,fo);
title('Salida del filtro adaptado');
xlabel('t(s)');

%% Determinación de los instantes óptimos de muestreo 
% Considere la descripción del algoritmo propuesto en el guión 
% de la práctica. Obtenga un vector "Ts" con los instantes muestrales
% que permita la correcta detección de los símbolos.
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

% Representamos con círculos rojos los instantes de muestreo óptimos
hold on
plot(tv(Ts),fo(Ts),'or');
hold off
grid on

%% Recuperación de la secuencia binaria a partir de los símbolos. Tenga en cuenta el 
%% tipo de código de línea utilizado para fijar el umbral (en este caso, polar-NRZ) 
% ---(Completar)
sbo= int64(s(Ts)>0); 


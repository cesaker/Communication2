%**************************************************************************
%                   COMUNICACIONES II(CURSO 2014/2015)
%          MACRO: REPRESENTACIÓN DEL DIAGRAMA DE OJO DE UNA SEÑAL
%**************************************************************************
%          José Carlos Segura, Carlos Medina, Ángel de la Torre
%                   {segura, cmedina, atv}@.ugr.es
%%*************************************************************************
% function diaojo(y,Nss,Fs)
%
% Argumentos de entrada:
% - y    : Señal de entrada
% - Nss  : Número de muestras por bit
% - Fs   : Frecuencia de muestreo [en Hz] 
%*************************************************************************
function diaojo(y,Nss,Fs)

Nst=length(y); % Número de muestras de señal
xmax=max(y);   % Amplitud máxima de señal
xmin=min(y);   % Amplitud mínima de señal

np=3;          % Trazamos tres periodos de señal
n=1;
%tv=[0:(np*Nss)-1]/Fs;
tv=[1:(np*Nss)]/Fs;

while(n<=(Nst-(np*Nss)))
  yy=y(n:n+(np*Nss)-1);
  plot(tv,yy);
  hold on;
  n=n+(np*Nss);
end

hold off;
grid on;
axis([min(tv) max(tv) xmin*1.1 xmax*1.1]);

end

%**************************************************************************
%                   COMUNICACIONES II(CURSO 2014/2015)
%          MACRO: REPRESENTACI�N DEL DIAGRAMA DE OJO DE UNA SE�AL
%**************************************************************************
%          Jos� Carlos Segura, Carlos Medina, �ngel de la Torre
%                   {segura, cmedina, atv}@.ugr.es
%%*************************************************************************
% function diaojo(y,Nss,Fs)
%
% Argumentos de entrada:
% - y    : Se�al de entrada
% - Nss  : N�mero de muestras por bit
% - Fs   : Frecuencia de muestreo [en Hz] 
%*************************************************************************
function diaojo(y,Nss,Fs)

Nst=length(y); % N�mero de muestras de se�al
xmax=max(y);   % Amplitud m�xima de se�al
xmin=min(y);   % Amplitud m�nima de se�al

np=3;          % Trazamos tres periodos de se�al
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

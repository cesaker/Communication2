%*************************************************************************
%                COMUNICACIONES II (CURSO 2014/2015)
%*************************************************************************
% -------------------------------------------------------------------------
%          José Carlos Segura, Carlos Medina, Ángel de la Torre
%                   {segura, cmedina, atv}@.ugr.es
% -------------------------------------------------------------------------
function [xR]=demod_fsk_nocoherente(sr,f1,f2,Fs,Nss,mostrar,Nmostrar)

% FUNCTION DEMOD_FSK_NOCOHERENTE demodula una señal modulada en FSK binaria
% de forma no coherente.
%
% Parámetros de entrada:
% |_ sr (vector): Señal modulada.
% |_ f1:  Frecuencia de la primera portadora generada localmente.
% |_ f2:  Frecuencia de la segunda portadora generada localmente.
% |_ Fs:  Frecuencia de muestreo.
% |_ Nss: Número de muestras por símbolo.
% |_ mostrar: Flag para visualizar gráficas. 
% |_ Nmostrar: Numero de periodos a mostrar en las gráficas.
%
% Parámetros de salida:
% |_ xR (vector): Señal binaria recuperada del canal 1.
% -------------------------------------------------------------------------

% *************************************************************************
% 1) Filtramos cada uno de los canales a la frecuencia de la portadora
% *************************************************************************

% ****************************************
% Filtro Paso Banda: CANAL 1
% ****************************************
N=300;                             % Orden del filtro FIR
B=[f1-1 f1+1];                     % Ancho de banda del filtro
Bnorm=B./(Fs/2);                   % Frecuencias de ancho de banda normalizadas  
c=fir1(N,Bnorm);                   % Vector de coeficientes del filtro FIR
sr_filt1=filtfilt(c,1,sr);         % Filtrado de señal usando un filtro de fase cero filtfilt de matlab

% ****************************************
% Filtro Paso Banda: CANAL 2
% ****************************************
N=300;                             % Orden del filtro FIR
B=[f2-1 f2+1];                     % Ancho de banda del filtro
Bnorm=B./(Fs/2);                   % Frecuencias de ancho de banda normalizadas  
c=fir1(N,Bnorm);                   % Vector de coeficientes del filtro FIR
sr_filt2=filtfilt(c,1,sr);         % Filtrado de señal usando un filtro de fase cero filtfilt de matlab

if mostrar
    figure
    subplot(3,1,1)
    plot(sr(1:Nmostrar*Nss))
    title('Señal modulada (antes de filtrar)','FontSize',12)
    xlabel('Nº de muestra','FontSize',12)
    ylabel('Amplitud (V)','FontSize',12)
    grid on
    subplot(3,1,2)
    plot(sr_filt1(1:Nmostrar*Nss))
    title('Señal a la salida del filtro (portadora 1)','FontSize',12)
    xlabel('Nº de muestra','FontSize',12)
    ylabel('Amplitud (V)','FontSize',12)
    grid on
    subplot(3,1,3)
    plot(sr_filt2(1:Nmostrar*Nss))
    title('Señal a la salida del filtro (portadora 2)','FontSize',12)
    xlabel('Nº de muestra','FontSize',12)
    ylabel('Amplitud (V)','FontSize',12)
    grid on
end

% *************************************************************************
% 2) Pasamos las señales por el detector de envolvente
% *************************************************************************
% Tomamos el valor absoluto de los dos canales. Usamos un filtro adaptado
% a la envolvente de la señal(es rectangular, de duración igual al periodo de 
% símbolo)
h=ones(Nss,1)/Nss;
sr_envelope1=filter(h,1,abs(sr_filt1));   
sr_envelope2=filter(h,1,abs(sr_filt2)); 

if mostrar
    figure
    subplot(2,1,1)
    hold on
    plot(sr_filt1(1:Nmostrar*Nss));
    plot(sr_envelope1(1:Nmostrar*Nss),'r');
    title('Señal a la salida del canal 1 (portadora 1)','FontSize',12)
    legend('Señal filtrada','Envolvente de señal')
    xlabel('Nº de muestra','FontSize',12)
    ylabel('Amplitud (V)','FontSize',12)
    grid on
    subplot(2,1,2)
    hold on
    plot(sr_filt2(1:Nmostrar*Nss))
    plot(sr_envelope2(1:Nmostrar*Nss),'r')
    title('Señal a la salida del canal 2 (portadora 2)','FontSize',12)
    legend('Señal filtrada','Envolvente de señal')
    xlabel('Nº de muestra','FontSize',12)
    ylabel('Amplitud (V)','FontSize',12)
    grid on
end

% *************************************************************************
% 3) Muestreamos las señales
% *************************************************************************
zr1=sr_envelope1(Nss:Nss:end);
zr2=sr_envelope2(Nss:Nss:end);

% *************************************************************************
% 4) Pasamos las señales por el decisor
% *************************************************************************
% xR=(zr1>zr2); El valor decodificado es '1' cuando (zr1>zr2) y '0' cuando (zr1<zr2)
% También se puede calcular considerando la diferencia entre zr1 y zr2: zd=(zr1-zr2)
% El valor decodificado es '1' cuando zd>0 (que implica zr1>zr2) y '0' cuando zd<0 
% (que implica zr1<zr2)

zd=(zr1-zr2);
xR=(zd>0);

if mostrar
    figure
    subplot(3,1,1)
    plot(zr1(1:Nmostrar),'-o');
    title('Señal recuperada a la salida del canal 1','FontSize',12)
    xlabel('Nº de muestra','FontSize',12)
    ylabel('Amplitud (V)','FontSize',12)
    grid on
    subplot(3,1,2)
    plot(zr2(1:Nmostrar),'-o');
    title('Señal recuperada a la salida del canal 2','FontSize',12)
    xlabel('Nº de muestra','FontSize',12)
    ylabel('Amplitud (V)','FontSize',12)
    grid on
    subplot(3,1,3)
    plot(zr1(1:Nmostrar)-zr2(1:Nmostrar),'-o');
    hold on
    stem(xR(1:Nmostrar),'k','linewidth',2)
    plot(0*ones(1,Nmostrar),'r'); % Trazamos el umbral de decisión (está en 0)
    legend('Señal muestreada','Señal binaria recuperada','Umbral de decisión')
    title('Símbolos recuperados a la salida del decisor','FontSize',12)
    xlabel('Nº de muestra','FontSize',12)
    ylabel('Amplitud (V)','FontSize',12)
    grid on
end

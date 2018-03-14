%*************************************************************************
%                COMUNICACIONES II (CURSO 2014/2015)
%*************************************************************************
%          Jos� Carlos Segura, Carlos Medina, �ngel de la Torre
%                   {segura, cmedina, atv}@.ugr.es
%*************************************************************************
% function [yfr, pulso]=gen_code(Nbits,Vp,code,Rb,kb,B,SNR)
%
% Argumentos de entrada:
% - Nbits: Longitud de la secuencia binaria 
% - Vp   : Amplitud de los bits (V)
% - code : C�digo de l�nea:  ['unipolar_nrz'] ['polar_nrz'] ['unipolar_rz']  
%                            ['bipolar_rz']   ['manchester_nrz'] 
% - Rb   : Tasa de bits por segundo 
% - kb   : factor de sobremuestreo [Fs=Rb*kb]
% - B    : Ancho de banda del canal (Hz) (paso-baja) 
% - SNR  : SNR por muestra (en dB)
%*************************************************************************
% Salida:
% - yfr  : Se�al a la salida del canal
% - pulso: Pulso usado por el transmisor
%
% (*) Esta versi�n tambi�n traza la salida del filtro adaptado
%*************************************************************************
function [yfr, pulso]=gen_code(Nbits,Vp,code,Rb,kb,B,SNR)

% Numero de muestras que se desplaza el muestreo a la salida del demodulador 
%(debe estar entre 0 y kb)
dsinc=0;

%% Generamos una secuencia aleatoria de N bits
x=randn(Nbits,1)>0;

%% Comentar lo de arriba para cargar una secuencia binaria concreta (p.e)
%x=[1 0 1 0 1 1 0 0 1 1]';

%% Variables:
Fs=Rb*kb;              % Frecuencia de muestreo
Ts=1/Fs;		       % Periodo de muestreo
Tb=1/Rb;		       % Periodo de bit/s�mbolo 
Nss=(Tb/Ts);           % N�mero de muestras por bit/s�mbolo
Nst=Nbits*(Tb/Ts);     % N�mero de muestras totales

if dsinc<0 || dsinc>Nss
    fprintf('Error: dsinc debe estar entre 0 y %d\n',Nss);
    yfr=[]; pulso=[];
    return;
end

%% Conformaci�n de pulso (p(t)) usando se�al cuadrada
switch code
    case 'unipolar_nrz'
        str=sprintf('%s','Codificaci�n UNIPOLAR-NRZ');
        pulso=ones(Nss,1); 
        cb=x(:);

    case 'polar_nrz'
        str=sprintf('%s','Codificaci�n POLAR-NRZ');
        pulso=ones(Nss,1); 
        cb=2*x(:)-1;

    case 'unipolar_rz'
        str=sprintf('%s','Codificaci�n UNIPOLAR-RZ');
        pulso=[ones(Nss/2,1);zeros(Nss/2,1)]; 
        cb=x(:);
    
    case 'bipolar_rz'
        str=sprintf('%s','Codificaci�n BIPOLAR-RZ');
        pulso=[ones(Nss/2,1);zeros(Nss/2,1)];  
        cb=2*x(:)-1;
      
        m=1;
        for q=1:size(cb,1)
            if(cb(q)==1)
                cb(q)=m;
                m=-m;
            else
                cb(q)=0;
            end
        end
    
    case 'manchester_nrz'
        str=sprintf('%s','C�dificaci�n MANCHESTER-NRZ)');
        pulso=[ones(Nss/2,1);-ones(Nss/2,1)];
        cb=2*x(:)-1;

    otherwise
        error('c�digo de l�nea desconocido');
end

%% Generamos un tren de impulsos con una separaci�n de Nss muestras (Tb) 
xn=zeros(Nst,1);
xn(1:Nss:end)=cb;

%% Convoluci�n con el filtro de transmisi�n h(t)=p(t):
y=filter(pulso,1,xn);

%% Aplicamos el factor de amplitud (Vp):
y=Vp*y;

%% Representaci�n de la se�al binaria x(t) codificada (los 10 primeros bits �nicamente)
% Si hay m�s de 10 bits
nb=floor(length(y)/kb);
if(nb>10)
    nb=10;
end
tl=[0:(nb*Nss)-1]*Ts;
yl=y(1:nb*Nss);

figure(1);
plot(tl,yl,'b','LineWidth',2.5);
xlabel('Tiempo (sec)'); ylabel('Amplitud (V)'); title(str);
%axis([min(tl) max(tl) min(yl-1) max(yl+1)]);
axis auto
grid on;

%%***********************************************************
%% Efectos de canal: Perturbaciones de se�al (ruido + ISI)
%***********************************************************
[yfr]=canal(y,B,Fs,SNR);

hold on;
yfrl=yfr(1:nb*Nss);
plot(tl,yfrl,'r','LineWidth',0.5);
hold off;

%%*************************************************************
%% Recuperaci�n del sincronismo 
%%*************************************************************
% Por defecto se fija el instante �ptimo de muestreo a la primera muestra 
% de cada periodo de bit.

%//////////////////////////////////////////////////////////////////////
%Instante �ptimo de muestreo determinado experimentalmente a partir del 
%diagrama de ojo (caso particular para "unipolar_nrz") 
%//////////////////////////////////////////////////////////////////////
%if(strcmp(code,'unipolar_nrz'))
   
    %% Usando informaci�n del diagrama de ojo
    figure(3);
    diaojo(yfr,Nss,Fs);
    xlabel('Tiempo (seg)');
    ylabel('Amplitud (V)');
    title(sprintf('Diagrama de ojo a salida del canal (Vp=%1.f V, Rb=%d bps,SNR=%1.f dB)',Vp,Rb,SNR));
    
    %% A�adido para hacer bien los diagramas de ojo: a la salida del
    %% demodulador!!!!   ATV OCT-2013
    norma=(sum(pulso.*pulso));
    filtro_adp=flipud(pulso)/norma;
    zfr = filter(filtro_adp,1,yfr);
    figure(4);
    diaojo(zfr,Nss,Fs);
    xlabel('Tiempo (seg)');
    ylabel('Amplitud (V)');
    title(sprintf('Diagrama de ojo a la salida del filtro adaptado (Vp=%1.f V, Rb=%d bps,SNR=%1.f dB)',Vp,Rb,SNR));
    
    % Instante �ptimo de muestreo seleccionado
    tm = Nss-dsinc;
    
    %% Representaci�n de los instantes de muestreo de se�al
    zm=zfr(tm:Nss:end);
    figure(1);
    hold on;
    tt=[1:nb*Nss]*Ts;
    plot(tt,zfr(1:nb*Nss),'k','LineWidth',1.5);
    zml=zm(1:nb);
    t=[tm:Nss:Nst]*Ts;
    tl=t(1:nb);
    stem(tl,zml,'k','MarkerFaceColor','g');
    hold off;  
    legend('s.transmitida','s.recibida','s.filtro adap.');
%end

%% Interpretamos la amplitud del pulso para saber si es un 1 � un 0
dAmp=Vp/2; % Fijamos el umbral de decisi�n a la mitad de amplitud del pulso

switch code
    case 'unipolar_nrz'
       xr=(zm>dAmp);
        
    case 'polar_nrz'
       xr=(zm>0);

    case 'unipolar_rz'
       xr=(zm>dAmp);
    
    case 'bipolar_rz'
       xr=(abs(zm)>dAmp);
    
    case 'manchester_nrz'
       xr=(zm>0);

    otherwise
       error('c�digo de l�nea desconocido');
end

%% N�mero bits erroneos (probabilidad de error)
Neb=sum(xr~=x);
Peb=(Neb/Nbits)+1e-6;

fprintf('\n**** Estad�stica errores: %s ****\n',code);
fprintf('N� bits transmitidos  = %d\n',Nbits);
fprintf('N� bits err�neos      = %d\n',Neb);
fprintf('Peb (experimental)    = %.5f\n',Peb);

end



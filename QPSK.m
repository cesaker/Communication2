%*************************************************************************
%                COMUNICACIONES II (CURSO 2014/2015)
%*************************************************************************
%          José Carlos Segura, Carlos Medina, Ángel de la Torre
%                   {segura, cmedina, atv}@.ugr.es
%*************************************************************************
function [Peb,Pes,Pebt,Pest,EbNodB,EsNodB,sr]=QPSK(Nbits,sigma,desfase,flag)

%% Argumentos de entrada:
% - Nbits  : Número de bits a transmitir
% - sigma  : Valor de desviación típica del ruido
% - desfase: desfase de la constelación (radianes)
% - flag   : (0) oculta gráficas y resultados
%            (1) muestra gráficas y resultados
%
%% Argumentos de salida:
% - Peb   : Probabilidad de error de bit (experimental)
% - Pes   : Probabilidad de error de símbolo (experimental)
% - Pebt  : Probabilidad de error de bit (teórica)
% - Pest  : Probabilidad de error de símbolo (teórica)
% - EbNodB: Energía de bit/Densidad Pot. Ruido (Eb/No) (en dB)
% - EsNodB: Energía de símbolo/Densidad Pot. Ruido (Es/No) (en dB)
% - sr    : Señal a la salida del canal (sin demodular)

%% Objetivo:
% Esta macro ilustra el uso del modulador en cuadratura para la generación y detección 
% de modulaciones en cuadratura de fase. En concreto, el código que sigue simula un 
% sistema QPSK síncrono.

%% Parámetros:
T=1e-3;      % Periodo de símbolo T=1ms
dt=0.025e-3; % Periodo de muestreo dt=0.025ms
Fs=1/dt;     % Fs=1/dt = 40KHz
Nss=T/dt;    % Número de muestras por periodo de símbolo
Fc=2e3;      % Frecuencia de portadora Fc=2KHz

% Trazado de las señales base del espacio portadora en fase 'ci(t)=cos(2*pi*Fc*t)'
% y en cuadratura 'cq(t)=sin(2*pi*Fc*t)':
t=[0:Nss-1]'*dt;   % Vector de tiempos
ci=cos(2*pi*Fc*t); % señal en fase
cq=sin(2*pi*Fc*t); % señal en cuadratura

if(flag)
    figure(1);
    subplot(2,1,1);
    plot(t,ci);
    xlabel('t(s)');
    title('Señal portadora en fase: ci(t)=cos(2*pi*Fc*t)');
    grid on;
    subplot(2,1,2);
    plot(t,cq);
    xlabel('t(s)');
    title('Señal portadora en cuadratura: cq(t)=sin(2*pi*Fc*t)');
    grid on;
end

% Generación de la cadena de bits (datos):
Bits=(randn(1,Nbits)>0);

% Valores de los símbolos (valor entero m):
Simbs = 2*Bits(1:2:end) + Bits(2:2:end);
Nbs   = 2;         % Número de bits por símbolo (QPSK)
Nsimb = Nbits/Nbs; % Número de símbolos

%% CONSTELACIÓN DE SEÑAL
% Tabla de valores para cada pareja de bits de entrada b1,b2 (pareja de bits
% que definen el símbolo (m), con m={0,1,2,3}, dado que M=4). En este caso 
% seleccionamos una energía de símbolo Es=1, por lo que las proyecciones de
% las componentes en fase (ai) y cuadratura (aq) tomarán el valor 1/sqrt(2).
% Es=(ai^2)+(aj^2)

% b1 b2   m      ai         aq
% -- --  --- ----------  ----------
% 0   0   0   +1/sqrt(2) +1/sqrt(2)
% 0   1   1   -1/sqrt(2) +1/sqrt(2)
% 1   0   2   +1/sqrt(2) -1/sqrt(2)
% 1   1   3   -1/sqrt(2) -1/sqrt(2)

% (NOTA): La tabla es construída considerando la constelación de
% símbolos mostrada en el guión de la práctica.

% Almacenamos las componentes en una tabla bidimensional:
% Amp(m,1) = amplitud en fase del símbolo m
% Amp(m,2) = amplitud en cuadratura del símbolo m

Amp=[ +1/sqrt(2) +1/sqrt(2)
      -1/sqrt(2) +1/sqrt(2)
      +1/sqrt(2) -1/sqrt(2)
      -1/sqrt(2) -1/sqrt(2)
    ];

%% Señales en banda base (en fase y en cuadratura)
% Generamos un tren de impulsos con separación de Nss muestras
% de amplitudes correpondientes a la constelación de señal
xil=[];
xql=[];
ximp=[1 ; zeros(Nss-1,1)];
for n=1:Nsimb
    m = Simbs(n); % Símbolo actual
    % NOTA: Se indexa m+1 porque Matlab usa índices comenzando en 1
    xil = [xil ; Amp(m+1,1)*ximp];
    xql = [xql ; Amp(m+1,2)*ximp];
end

% Filtro del modulador (conformación de pulsos):
p=ones(Nss,1); % Forma de onda NRZ-L (igual en fase y cuadratura)
p=p/norm(p);   % Para normalizar en energía el pulso
xi=filter(p,1,xil);
xq=filter(p,1,xql);

% Composición de la señal modulada:
t=[0:length(xi)-1]'*dt; % Vector de tiempos
si = xi.*cos(2*pi*Fc*t);
sq = xq.*sin(2*pi*Fc*t);
s  = si - sq;

% Trazado de los datos (primeros símbolos):
Ns=16;
N=Ns*Nss;
if(flag)
    figure(2);
    subplot(3,1,1); plot(t(1:N),si(1:N)); title('Componente de señal: fase (si)');grid on;
    subplot(3,1,2); plot(t(1:N),sq(1:N)); title('Componente de señal: cuadratura (sq)');grid on;
    subplot(3,1,3); plot(t(1:N),s(1:N));  title('Señal modulada: s=si-sq');grid on; 
end

% Ruido del canal:
sr = s + randn(size(s))*sigma;

if(flag)
    figure(3);
    subplot(2,1,1); plot(t(1:N),si(1:N));
    plot(t(1:N),s(1:N));
    title('Señal a la entrada del canal: s');
    xlabel('t(s)');
    grid on;

    subplot(2,1,2); plot(t(1:N),sq(1:N));
    plot(t(1:N),sr(1:N));
    title('Señal a la salida del canal: sr=s+noise');
    xlabel('t(s)');
    grid on;
end

%************************************
%% Demodulación QPSK
%************************************
% Producto por las funciones base (demodulación coherente):
% Nota: cos(a)cos(b)=1/2[cos(a+b)+cos(a-b)].
% Multiplicamos por un factor 2 para compensar ese escalado 
sri =  sr.*2.*cos(2*pi*Fc*t+desfase); 
srq = -sr.*2.*sin(2*pi*Fc*t+desfase);

% Integración (filtro adaptado a la forma de onda en banda-base):
pa=flipud(p);
zri = filter(pa,1,sri);
zrq = filter(pa,1,srq);

% Muestreo al final del periodo de símbolo:
Zri = zri(Nss:Nss:end);
Zrq = zrq(Nss:Nss:end);

if(flag)
    % Trazado de los datos (primeros N símbolos)
    figure(4);
    subplot(2,1,1);
    plot(t(1:N),zri(1:N));
    hold on
    plot(t(Nss:Nss:N),Zri(1:Ns),'or');
    hold off
    title('Salida del filtro adaptado (fase): zri=fadap(sri)');
    xlabel('t(s)');
    grid on
    axis([-Inf Inf -1.5 1.5]);

    subplot(2,1,2);
    plot(t(1:N),zrq(1:N));
    hold on
    plot(t(Nss:Nss:N),Zrq(1:Ns),'or');
    hold off
    title('Salida del filtro adaptado (cuadratura): zrq=fadap(srq)');
    xlabel('t(s)');
    grid on
    axis([-Inf Inf -1.5 1.5]);

    % Constelación de salida (todos los datos)
    figure(5);
    plot(Zri,Zrq,'.');
    title('Constelación de la señal QPSK')
    xlabel('Componente en fase (zri)');
    ylabel('Componente en cuadratura (zrq)');
    axis([-2 2 -2 2]);
    axis equal
    grid on
    
    for p=1:size(Amp,1)
       text(Amp(p,1),Amp(p,2)+0.2,sprintf(dec2bin(p-1,2)));
    end 
end

%% Detección de los bits (mínima distancia):
M=size(Amp,1); % Número de símbolos en el alfabeto
Rsimbs=[];     % Símbolos recibidos
Rbits=[];      % Bits recibidos

for n=1:length(Simbs)
    x=Zri(n);
    y=Zrq(n);
    dmin=1e9;
    mmin=0;
    for m=1:M
        d=sqrt((Amp(m,1)-x).^2+(Amp(m,2)-y).^2);
        if(d<dmin)
            dmin=d;
            mmin=m;
        end
    end
    simbolo = mmin-1; % Nota: consideramos la indexación en MATLAB
    Rsimbs=[Rsimbs simbolo];
   
    % Separación de los bits
    b1 = floor(simbolo/2); % El bit impar
    b2 = simbolo-2*b1;     % El bit par
    Rbits=[Rbits b1 b2];
end

% Número de errores de símbolo y de bit:
Nes=sum(Simbs ~= Rsimbs);
Neb=sum(Bits ~= Rbits);
Pes=(Nes/Nsimb);
Peb=(Neb/Nbits);

if(flag)
    fprintf('Número de símbolos=%d\n',Nsimb);
    fprintf('Número de bits=%d\n',Nbits);
    fprintf('Número de símbolos erróneos=%d\t Pes= %f\n',Nes,Nes/Nsimb);
    fprintf('Número de bits erróneos=%d\t\t Peb = %f\n',Neb,Neb/Nbits);
end

%*************************************************************************
%% Estimación teórica de Eb/No (denotado como EbNo):
%*************************************************************************
% 1º Solución:
%        N=(No/2)*Fs=sigma.^2 => No=(2/Fs)*(sigma.^2)        [1.1]
%        Eb=(A.^2)*Nss/(Nbs*Fs) con Nbs=2 (QPSK,M=2.^Nbs=4)  [1.2]
%        EbNo=(Eb/No)=(A.^2)*Nss/(4*sigma.^2) 
%
% 2º Solución:
%        S/N=(A.^2)/sigma.^2;                                        [1.3]
%        Tb*W= Tb*(Fs/2)=(Nsb/Fs)*(Fs/2)=Nsb/2=Nss/4 con Nsb=Nss/2   [1.4] 
%        EbNo=(S/N)*(Tb*W)=[((A.^2)/sigma.^2)*(Nss/4)]=(A.^2)*Nss/(4*sigma.^2)
%*************************************************************************
% Estimamos la amplitud de la señal modulada (A) para calcular la energía 
% de bit (Eb). Hay que considerar el factor 1/sqrt(2) que introduce la 
% modulación de señal.
A=max(s)/sqrt(2);    
EbNo=(A.^2)*Nss/(4*sigma.^2);
EbNodB=10*log10(EbNo);
EsNodB=10*log10(2*EbNo); % Es=Eb*log2(M)=2*Eb

%*************************************************************************
%% Estimación experimental de Eb/No 
%*************************************************************************
% Usamos el valor de Peb para despejar el argumento de la función Q y con 
% él el valor de Eb/No
%************************************************************************
% xb=erfcinv(2*Peb)*sqrt(2);
% EbNo=(xb.^2)/2;
% EbNodB=10*log10(EbNo);
% EsNodB=10*log10(2*EbNo); % Es=Eb*log2(M)=2*Eb

%% Peb=Q(sqrt(2*Eb/No)): 
xb=sqrt(2*EbNo);
Pebt=0.5*erfc(xb/sqrt(2));

%% Pes=2*Q(sqrt(Es/No)); Es=Eb*log2(M)=2*Eb; (aprox.):
xs=sqrt(2*EbNo);
Pest=erfc(xs/sqrt(2));      

if(flag)
    fprintf('------------------------------------------\n');
    fprintf('Comprobación (Teoria <-> Experimental)\n');
    fprintf('------------------------------------------\n');
    fprintf('Datos:\t sigma=%1.3f\t Nbits=%d\n',sigma,Nbits);
    fprintf('      \t EbNodB=%1.3f\t EsNodB=%1.3f\n',EbNodB,EsNodB); 
    fprintf('Prob:\t Peb=%1.4f\t\t Pebt=%1.4f\n',Peb,Pebt); 
    fprintf('     \t Pes=%1.4f\t\t Pest=%1.4f\n',Pes,Pest);
    fprintf('-------------------------------------------\n');
end

end


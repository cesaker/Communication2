%*************************************************************************
%                COMUNICACIONES II (CURSO 2014/2015)
%*************************************************************************
%          Jos� Carlos Segura, Carlos Medina, �ngel de la Torre
%                   {segura, cmedina, atv}@.ugr.es
%*************************************************************************
function [Peb,Pes,Pebt,Pest,EbNodB,EsNodB,sr]=QPSK(Nbits,sigma,desfase,flag)

%% Argumentos de entrada:
% - Nbits  : N�mero de bits a transmitir
% - sigma  : Valor de desviaci�n t�pica del ruido
% - desfase: desfase de la constelaci�n (radianes)
% - flag   : (0) oculta gr�ficas y resultados
%            (1) muestra gr�ficas y resultados
%
%% Argumentos de salida:
% - Peb   : Probabilidad de error de bit (experimental)
% - Pes   : Probabilidad de error de s�mbolo (experimental)
% - Pebt  : Probabilidad de error de bit (te�rica)
% - Pest  : Probabilidad de error de s�mbolo (te�rica)
% - EbNodB: Energ�a de bit/Densidad Pot. Ruido (Eb/No) (en dB)
% - EsNodB: Energ�a de s�mbolo/Densidad Pot. Ruido (Es/No) (en dB)
% - sr    : Se�al a la salida del canal (sin demodular)

%% Objetivo:
% Esta macro ilustra el uso del modulador en cuadratura para la generaci�n y detecci�n 
% de modulaciones en cuadratura de fase. En concreto, el c�digo que sigue simula un 
% sistema QPSK s�ncrono.

%% Par�metros:
T=1e-3;      % Periodo de s�mbolo T=1ms
dt=0.025e-3; % Periodo de muestreo dt=0.025ms
Fs=1/dt;     % Fs=1/dt = 40KHz
Nss=T/dt;    % N�mero de muestras por periodo de s�mbolo
Fc=2e3;      % Frecuencia de portadora Fc=2KHz

% Trazado de las se�ales base del espacio portadora en fase 'ci(t)=cos(2*pi*Fc*t)'
% y en cuadratura 'cq(t)=sin(2*pi*Fc*t)':
t=[0:Nss-1]'*dt;   % Vector de tiempos
ci=cos(2*pi*Fc*t); % se�al en fase
cq=sin(2*pi*Fc*t); % se�al en cuadratura

if(flag)
    figure(1);
    subplot(2,1,1);
    plot(t,ci);
    xlabel('t(s)');
    title('Se�al portadora en fase: ci(t)=cos(2*pi*Fc*t)');
    grid on;
    subplot(2,1,2);
    plot(t,cq);
    xlabel('t(s)');
    title('Se�al portadora en cuadratura: cq(t)=sin(2*pi*Fc*t)');
    grid on;
end

% Generaci�n de la cadena de bits (datos):
Bits=(randn(1,Nbits)>0);

% Valores de los s�mbolos (valor entero m):
Simbs = 2*Bits(1:2:end) + Bits(2:2:end);
Nbs   = 2;         % N�mero de bits por s�mbolo (QPSK)
Nsimb = Nbits/Nbs; % N�mero de s�mbolos

%% CONSTELACI�N DE SE�AL
% Tabla de valores para cada pareja de bits de entrada b1,b2 (pareja de bits
% que definen el s�mbolo (m), con m={0,1,2,3}, dado que M=4). En este caso 
% seleccionamos una energ�a de s�mbolo Es=1, por lo que las proyecciones de
% las componentes en fase (ai) y cuadratura (aq) tomar�n el valor 1/sqrt(2).
% Es=(ai^2)+(aj^2)

% b1 b2   m      ai         aq
% -- --  --- ----------  ----------
% 0   0   0   +1/sqrt(2) +1/sqrt(2)
% 0   1   1   -1/sqrt(2) +1/sqrt(2)
% 1   0   2   +1/sqrt(2) -1/sqrt(2)
% 1   1   3   -1/sqrt(2) -1/sqrt(2)

% (NOTA): La tabla es constru�da considerando la constelaci�n de
% s�mbolos mostrada en el gui�n de la pr�ctica.

% Almacenamos las componentes en una tabla bidimensional:
% Amp(m,1) = amplitud en fase del s�mbolo m
% Amp(m,2) = amplitud en cuadratura del s�mbolo m

Amp=[ +1/sqrt(2) +1/sqrt(2)
      -1/sqrt(2) +1/sqrt(2)
      +1/sqrt(2) -1/sqrt(2)
      -1/sqrt(2) -1/sqrt(2)
    ];

%% Se�ales en banda base (en fase y en cuadratura)
% Generamos un tren de impulsos con separaci�n de Nss muestras
% de amplitudes correpondientes a la constelaci�n de se�al
xil=[];
xql=[];
ximp=[1 ; zeros(Nss-1,1)];
for n=1:Nsimb
    m = Simbs(n); % S�mbolo actual
    % NOTA: Se indexa m+1 porque Matlab usa �ndices comenzando en 1
    xil = [xil ; Amp(m+1,1)*ximp];
    xql = [xql ; Amp(m+1,2)*ximp];
end

% Filtro del modulador (conformaci�n de pulsos):
p=ones(Nss,1); % Forma de onda NRZ-L (igual en fase y cuadratura)
p=p/norm(p);   % Para normalizar en energ�a el pulso
xi=filter(p,1,xil);
xq=filter(p,1,xql);

% Composici�n de la se�al modulada:
t=[0:length(xi)-1]'*dt; % Vector de tiempos
si = xi.*cos(2*pi*Fc*t);
sq = xq.*sin(2*pi*Fc*t);
s  = si - sq;

% Trazado de los datos (primeros s�mbolos):
Ns=16;
N=Ns*Nss;
if(flag)
    figure(2);
    subplot(3,1,1); plot(t(1:N),si(1:N)); title('Componente de se�al: fase (si)');grid on;
    subplot(3,1,2); plot(t(1:N),sq(1:N)); title('Componente de se�al: cuadratura (sq)');grid on;
    subplot(3,1,3); plot(t(1:N),s(1:N));  title('Se�al modulada: s=si-sq');grid on; 
end

% Ruido del canal:
sr = s + randn(size(s))*sigma;

if(flag)
    figure(3);
    subplot(2,1,1); plot(t(1:N),si(1:N));
    plot(t(1:N),s(1:N));
    title('Se�al a la entrada del canal: s');
    xlabel('t(s)');
    grid on;

    subplot(2,1,2); plot(t(1:N),sq(1:N));
    plot(t(1:N),sr(1:N));
    title('Se�al a la salida del canal: sr=s+noise');
    xlabel('t(s)');
    grid on;
end

%************************************
%% Demodulaci�n QPSK
%************************************
% Producto por las funciones base (demodulaci�n coherente):
% Nota: cos(a)cos(b)=1/2[cos(a+b)+cos(a-b)].
% Multiplicamos por un factor 2 para compensar ese escalado 
sri =  sr.*2.*cos(2*pi*Fc*t+desfase); 
srq = -sr.*2.*sin(2*pi*Fc*t+desfase);

% Integraci�n (filtro adaptado a la forma de onda en banda-base):
pa=flipud(p);
zri = filter(pa,1,sri);
zrq = filter(pa,1,srq);

% Muestreo al final del periodo de s�mbolo:
Zri = zri(Nss:Nss:end);
Zrq = zrq(Nss:Nss:end);

if(flag)
    % Trazado de los datos (primeros N s�mbolos)
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

    % Constelaci�n de salida (todos los datos)
    figure(5);
    plot(Zri,Zrq,'.');
    title('Constelaci�n de la se�al QPSK')
    xlabel('Componente en fase (zri)');
    ylabel('Componente en cuadratura (zrq)');
    axis([-2 2 -2 2]);
    axis equal
    grid on
    
    for p=1:size(Amp,1)
       text(Amp(p,1),Amp(p,2)+0.2,sprintf(dec2bin(p-1,2)));
    end 
end

%% Detecci�n de los bits (m�nima distancia):
M=size(Amp,1); % N�mero de s�mbolos en el alfabeto
Rsimbs=[];     % S�mbolos recibidos
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
    simbolo = mmin-1; % Nota: consideramos la indexaci�n en MATLAB
    Rsimbs=[Rsimbs simbolo];
   
    % Separaci�n de los bits
    b1 = floor(simbolo/2); % El bit impar
    b2 = simbolo-2*b1;     % El bit par
    Rbits=[Rbits b1 b2];
end

% N�mero de errores de s�mbolo y de bit:
Nes=sum(Simbs ~= Rsimbs);
Neb=sum(Bits ~= Rbits);
Pes=(Nes/Nsimb);
Peb=(Neb/Nbits);

if(flag)
    fprintf('N�mero de s�mbolos=%d\n',Nsimb);
    fprintf('N�mero de bits=%d\n',Nbits);
    fprintf('N�mero de s�mbolos err�neos=%d\t Pes= %f\n',Nes,Nes/Nsimb);
    fprintf('N�mero de bits err�neos=%d\t\t Peb = %f\n',Neb,Neb/Nbits);
end

%*************************************************************************
%% Estimaci�n te�rica de Eb/No (denotado como EbNo):
%*************************************************************************
% 1� Soluci�n:
%        N=(No/2)*Fs=sigma.^2 => No=(2/Fs)*(sigma.^2)        [1.1]
%        Eb=(A.^2)*Nss/(Nbs*Fs) con Nbs=2 (QPSK,M=2.^Nbs=4)  [1.2]
%        EbNo=(Eb/No)=(A.^2)*Nss/(4*sigma.^2) 
%
% 2� Soluci�n:
%        S/N=(A.^2)/sigma.^2;                                        [1.3]
%        Tb*W= Tb*(Fs/2)=(Nsb/Fs)*(Fs/2)=Nsb/2=Nss/4 con Nsb=Nss/2   [1.4] 
%        EbNo=(S/N)*(Tb*W)=[((A.^2)/sigma.^2)*(Nss/4)]=(A.^2)*Nss/(4*sigma.^2)
%*************************************************************************
% Estimamos la amplitud de la se�al modulada (A) para calcular la energ�a 
% de bit (Eb). Hay que considerar el factor 1/sqrt(2) que introduce la 
% modulaci�n de se�al.
A=max(s)/sqrt(2);    
EbNo=(A.^2)*Nss/(4*sigma.^2);
EbNodB=10*log10(EbNo);
EsNodB=10*log10(2*EbNo); % Es=Eb*log2(M)=2*Eb

%*************************************************************************
%% Estimaci�n experimental de Eb/No 
%*************************************************************************
% Usamos el valor de Peb para despejar el argumento de la funci�n Q y con 
% �l el valor de Eb/No
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
    fprintf('Comprobaci�n (Teoria <-> Experimental)\n');
    fprintf('------------------------------------------\n');
    fprintf('Datos:\t sigma=%1.3f\t Nbits=%d\n',sigma,Nbits);
    fprintf('      \t EbNodB=%1.3f\t EsNodB=%1.3f\n',EbNodB,EsNodB); 
    fprintf('Prob:\t Peb=%1.4f\t\t Pebt=%1.4f\n',Peb,Pebt); 
    fprintf('     \t Pes=%1.4f\t\t Pest=%1.4f\n',Pes,Pest);
    fprintf('-------------------------------------------\n');
end

end


%*************************************************************************
%                COMUNICACIONES II (CURSO 2014/2015)
%*************************************************************************
%          José Carlos Segura, Carlos Medina, Ángel de la Torre
%                   {segura, cmedina, atv}@.ugr.es
%*************************************************************************
function [Pes,Peb,Pebt,Pest,EbNodB,EsNodB,sr]=QAM16(Nbits,sigma,desfase,flag)

%% Argumentos de entrada:
% - Nbits  : Número de bits a transmitir
% - sigma  : Valor de desviación típica del ruido
% - desfase: desfase de la constelación (radianes)
% - flag   : (0) oculta gráficas y resultados
%            (1) muestra gráficas y resultados
%
%% Objetivo:
% Esta macro simula un sistema un sistema QAM16 síncrono.

% Parámetros:
T=1e-3;      % Periodo de símbolo T=1ms
dt=0.025e-3; % Periodo de muestreo dt=0.025ms
Fs=1/dt;     % Fs=1/dt = 40KHz
Nss=T/dt;    % Número de muestras por periodo de símbolo
Fc=2e3;      % Frecuencia de portadora Fc=2KHz

% Señales base en fase 'ci(t)=cos(2*pi*Fc*t)' y en cuadratura
% 'cq(t)=sin(2*pi*Fc*t)':
t=[0:Nss-1]'*dt;       % Vector de tiempos
si = cos(2*pi*Fc*t);   % señal en fase 
sq = sin(2*pi*Fc*t);   % señal en cuadratura

if(flag)
    figure(1);
    subplot(2,1,1);
    plot(t,si);
    xlabel('t(s)');
    title('Señal portadora en fase: ci(t)=cos(2*pi*Fc*t)');
    subplot(2,1,2);
    plot(t,sq);
    xlabel('t(s)');
    title('Señal portadora en cuadratura: cq(t)=sin(2*pi*Fc*t)');
end

% Generación de la cadena de bits (datos):
Bits=(randn(1,Nbits)>0);


%% Modulación QAM16
% Tabla de símbolos para la constelación QAM16. Cada símbolo está formado
% por un conjunto de 4 bits {b1,b2,b3,b4}. Su representación en el espacio 
% de señal se puede expresar a través de las componentes en fase (ai) y
% cuadratura (aq).
% NOTA: El factor 1/sqrt(2) se aplica después para no complicar la tabla

% b1 b2 b3 b4  m  ai aq
% -- -- -- -- --- -- --
%  0  0  0  0   0 -3 -3
%  0  0  0  1   1 -3 -1
%  0  0  1  0   2 -3  3
%  0  0  1  1   3 -3  1
%  0  1  0  0   4 -1 -3
%  0  1  0  1   5 -1 -1
%  0  1  1  0   6 -1  3
%  0  1  1  1   7 -1  1
%  1  0  0  0   8  3 -3
%  1  0  0  1   9  3 -1
%  1  0  1  0  10  3  3
%  1  0  1  1  11  3  1
%  1  1  0  0  12  1 -3
%  1  1  0  1  13  1 -1
%  1  1  1  0  14  1  3
%  1  1  1  1  15  1  1

% Rellenamos una matriz (Amp) con los valores de ai y aq de la tabla.
% La primera columna corresponde a los valores de "ai" y la segunda a 
% los valores de "aq" (normalizados ya con 1/sqrt(2)):
Amp=[-3 -3; -3 -1; -3 3; -3 1; -1 -3 ; -1 -1; -1 3; -1 1; 3 -3; 3 -1; 3 3; 3 1; 1 -3; 1 -1; 1 3; 1 1]./sqrt(2); 

% Recuperamos el vector de símbolos (valor entero m) a partir de la 
% secuencia binaria generada:
Nbs   = 4;            % Número de bits por símbolo 
Nsimb = Nbits/Nbs;    % Número de símbolos
Simbs = 8*Bits(1:4:end) + 4*Bits(2:4:end) + 2*Bits(3:4:end) + Bits(4:4:end);

%% Señales en banda base (en fase y en cuadratura)
% Generamos un tren de impulsos con separación de Nss muestras y
% amplitudes correpondientes a la constelación:
xil=[];
xql=[];
ximp=[1 ; zeros(Nss-1,1)];

for n=1:Nsimb
    m = Simbs(n); % Símbolo actual
    % NOTA: Se indexa con k+1 porque Matlab usa índices comenzando en 1
    xil = [xil ; Amp(m+1,1)*ximp];
    xql = [xql ; Amp(m+1,2)*ximp];
end

% Filtro del modulador (conformación de pulsos):
p=ones(Nss,1); % Forma de onda NRZ-L (igual en fase y cuadratura)
p=p/norm(p);   % Para normalizar en energía el pulso
xi=filter(p,1,xil);
xq=filter(p,1,xql);

% Composición de la señal modulada:
t=[0:length(xi)-1]'*dt;
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

%% Ruido del canal:
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

%% Demodulador QAM
% Proyecciones en fase y cuadratura

% Producto por las funciones base (demodulación coherente):
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
    axis([-Inf Inf -3.5 3.5]);

    subplot(2,1,2);
    plot(t(1:N),zrq(1:N));
    hold on
    plot(t(Nss:Nss:N),Zrq(1:Ns),'or');
    hold off
    title('Salida del filtro adaptado (cuadratura): zrq=fadap(srq)');
    xlabel('t(s)');
    axis([-Inf Inf -3.5 3.5]);
   
    % Constelación de salida (todos los datos):
    figure(5);
    plot(Zri,Zrq,'ob','MarkerFaceColor','r','LineWidth',1.0);
    title('Constelación QAM16')
    xlabel('Componente en fase (zri)');
    ylabel('Componente en cuadratura (zrq)');
    axis([-4.5 4.5 -4.5 4.5]);
    axis equal
    grid on
    
    %Etiquetado de los símbolos:
    for p=1:size(Amp,1)
      text(Amp(p,1)-0.3,Amp(p,2)+0.3,sprintf(dec2bin(p-1,4)));
    end 
end

%% Deteción de los símbolos y los bits: 
% Nos basamos en el cálculo de mínima distancia entre el símbolo recibido 
% (representado por sus componentes en fase y cuadratura) y los símbolos 
% de la constelación dados por la matriz "Amp".

M=size(Amp,1);   % Número de símbolos en el alfabeto
Rsimbs=[];       % Símbolos recibidos
Rbits=[];        % Bits recibidos

for n=1:length(Simbs)
    x=Zri(n);
    y=Zrq(n);
    dmin=1e9;
    mmin=0;
    
    for m=1:M
      d=(x-Amp(m,1)).^2+(y-Amp(m,2)).^2;
      if(d<dmin)
          dmin=d;
          mmin=m;
      end
    end
    
    simbolo = mmin-1; % Nota: consideramos la indexación en MATLAB
    Rsimbs=[Rsimbs simbolo];
   
    % Recuperamos la secuencia de bits dado el símbolo: 
    bb=dec2bin(simbolo,4);
    b1 = bb(1)-'0'; 
    b2 = bb(2)-'0';
    b3 = bb(3)-'0';
    b4 = bb(4)-'0';
    Rbits=[Rbits b1 b2 b3 b4];
end

% Número de errores de símbolo y de bit:
Nes=sum(Simbs~=Rsimbs);
Neb=sum(Bits~=Rbits);
Pes=Nes/length(Rsimbs);
Peb=Neb/length(Rbits);

if(flag)
    fprintf('Número de símbolos=%d\n',Nsimb);
    fprintf('Número de bits=%d\n',Nbits);
    fprintf('Número de símbolos erróneos = %d Probabilidad = %f\n',Nes,Pes);
    fprintf('Número de bits erróneos = %d Probabilidad = %f\n',Neb,Peb);
end

%***************************************************
%% Estimación experimental de la relación Eb/No
%***************************************************
% Usamos el valor de Peb para despejar el argumento de la función
% Q y con él valor de Eb/No.
M=16;
L=sqrt(M);
xb=sqrt(2)*erfcinv(log2(L)*Peb/(1-1/L));
EbNo=0.5*(xb^2*(L^2-1)/(3*log2(L)));
EsNo=4*EbNo;

% Hacemos uso de la función Q(x)=0.5*erfc(x/sqrt(2)) y el vector EbNo
% estimado para obtener los valores de probabilidad teóricos:
L=sqrt(M);
xb=sqrt(2*EbNo*(3*log2(L))/(L^2-1));
Pebt=(2*(1-1/L)/log2(L))*0.5*erfc(xb/sqrt(2));
Pest=log2(M)*Pebt;
EbNodB=10*log10(EbNo);
EsNodB=10*log10(EsNo);

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


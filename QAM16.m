%*************************************************************************
%                COMUNICACIONES II (CURSO 2014/2015)
%*************************************************************************
%          Jos� Carlos Segura, Carlos Medina, �ngel de la Torre
%                   {segura, cmedina, atv}@.ugr.es
%*************************************************************************
function [Pes,Peb,Pebt,Pest,EbNodB,EsNodB,sr]=QAM16(Nbits,sigma,desfase,flag)

%% Argumentos de entrada:
% - Nbits  : N�mero de bits a transmitir
% - sigma  : Valor de desviaci�n t�pica del ruido
% - desfase: desfase de la constelaci�n (radianes)
% - flag   : (0) oculta gr�ficas y resultados
%            (1) muestra gr�ficas y resultados
%
%% Objetivo:
% Esta macro simula un sistema un sistema QAM16 s�ncrono.

% Par�metros:
T=1e-3;      % Periodo de s�mbolo T=1ms
dt=0.025e-3; % Periodo de muestreo dt=0.025ms
Fs=1/dt;     % Fs=1/dt = 40KHz
Nss=T/dt;    % N�mero de muestras por periodo de s�mbolo
Fc=2e3;      % Frecuencia de portadora Fc=2KHz

% Se�ales base en fase 'ci(t)=cos(2*pi*Fc*t)' y en cuadratura
% 'cq(t)=sin(2*pi*Fc*t)':
t=[0:Nss-1]'*dt;       % Vector de tiempos
si = cos(2*pi*Fc*t);   % se�al en fase 
sq = sin(2*pi*Fc*t);   % se�al en cuadratura

if(flag)
    figure(1);
    subplot(2,1,1);
    plot(t,si);
    xlabel('t(s)');
    title('Se�al portadora en fase: ci(t)=cos(2*pi*Fc*t)');
    subplot(2,1,2);
    plot(t,sq);
    xlabel('t(s)');
    title('Se�al portadora en cuadratura: cq(t)=sin(2*pi*Fc*t)');
end

% Generaci�n de la cadena de bits (datos):
Bits=(randn(1,Nbits)>0);


%% Modulaci�n QAM16
% Tabla de s�mbolos para la constelaci�n QAM16. Cada s�mbolo est� formado
% por un conjunto de 4 bits {b1,b2,b3,b4}. Su representaci�n en el espacio 
% de se�al se puede expresar a trav�s de las componentes en fase (ai) y
% cuadratura (aq).
% NOTA: El factor 1/sqrt(2) se aplica despu�s para no complicar la tabla

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

% Recuperamos el vector de s�mbolos (valor entero m) a partir de la 
% secuencia binaria generada:
Nbs   = 4;            % N�mero de bits por s�mbolo 
Nsimb = Nbits/Nbs;    % N�mero de s�mbolos
Simbs = 8*Bits(1:4:end) + 4*Bits(2:4:end) + 2*Bits(3:4:end) + Bits(4:4:end);

%% Se�ales en banda base (en fase y en cuadratura)
% Generamos un tren de impulsos con separaci�n de Nss muestras y
% amplitudes correpondientes a la constelaci�n:
xil=[];
xql=[];
ximp=[1 ; zeros(Nss-1,1)];

for n=1:Nsimb
    m = Simbs(n); % S�mbolo actual
    % NOTA: Se indexa con k+1 porque Matlab usa �ndices comenzando en 1
    xil = [xil ; Amp(m+1,1)*ximp];
    xql = [xql ; Amp(m+1,2)*ximp];
end

% Filtro del modulador (conformaci�n de pulsos):
p=ones(Nss,1); % Forma de onda NRZ-L (igual en fase y cuadratura)
p=p/norm(p);   % Para normalizar en energ�a el pulso
xi=filter(p,1,xil);
xq=filter(p,1,xql);

% Composici�n de la se�al modulada:
t=[0:length(xi)-1]'*dt;
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

%% Ruido del canal:
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

%% Demodulador QAM
% Proyecciones en fase y cuadratura

% Producto por las funciones base (demodulaci�n coherente):
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
    axis([-Inf Inf -3.5 3.5]);

    subplot(2,1,2);
    plot(t(1:N),zrq(1:N));
    hold on
    plot(t(Nss:Nss:N),Zrq(1:Ns),'or');
    hold off
    title('Salida del filtro adaptado (cuadratura): zrq=fadap(srq)');
    xlabel('t(s)');
    axis([-Inf Inf -3.5 3.5]);
   
    % Constelaci�n de salida (todos los datos):
    figure(5);
    plot(Zri,Zrq,'ob','MarkerFaceColor','r','LineWidth',1.0);
    title('Constelaci�n QAM16')
    xlabel('Componente en fase (zri)');
    ylabel('Componente en cuadratura (zrq)');
    axis([-4.5 4.5 -4.5 4.5]);
    axis equal
    grid on
    
    %Etiquetado de los s�mbolos:
    for p=1:size(Amp,1)
      text(Amp(p,1)-0.3,Amp(p,2)+0.3,sprintf(dec2bin(p-1,4)));
    end 
end

%% Deteci�n de los s�mbolos y los bits: 
% Nos basamos en el c�lculo de m�nima distancia entre el s�mbolo recibido 
% (representado por sus componentes en fase y cuadratura) y los s�mbolos 
% de la constelaci�n dados por la matriz "Amp".

M=size(Amp,1);   % N�mero de s�mbolos en el alfabeto
Rsimbs=[];       % S�mbolos recibidos
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
    
    simbolo = mmin-1; % Nota: consideramos la indexaci�n en MATLAB
    Rsimbs=[Rsimbs simbolo];
   
    % Recuperamos la secuencia de bits dado el s�mbolo: 
    bb=dec2bin(simbolo,4);
    b1 = bb(1)-'0'; 
    b2 = bb(2)-'0';
    b3 = bb(3)-'0';
    b4 = bb(4)-'0';
    Rbits=[Rbits b1 b2 b3 b4];
end

% N�mero de errores de s�mbolo y de bit:
Nes=sum(Simbs~=Rsimbs);
Neb=sum(Bits~=Rbits);
Pes=Nes/length(Rsimbs);
Peb=Neb/length(Rbits);

if(flag)
    fprintf('N�mero de s�mbolos=%d\n',Nsimb);
    fprintf('N�mero de bits=%d\n',Nbits);
    fprintf('N�mero de s�mbolos err�neos = %d Probabilidad = %f\n',Nes,Pes);
    fprintf('N�mero de bits err�neos = %d Probabilidad = %f\n',Neb,Peb);
end

%***************************************************
%% Estimaci�n experimental de la relaci�n Eb/No
%***************************************************
% Usamos el valor de Peb para despejar el argumento de la funci�n
% Q y con �l valor de Eb/No.
M=16;
L=sqrt(M);
xb=sqrt(2)*erfcinv(log2(L)*Peb/(1-1/L));
EbNo=0.5*(xb^2*(L^2-1)/(3*log2(L)));
EsNo=4*EbNo;

% Hacemos uso de la funci�n Q(x)=0.5*erfc(x/sqrt(2)) y el vector EbNo
% estimado para obtener los valores de probabilidad te�ricos:
L=sqrt(M);
xb=sqrt(2*EbNo*(3*log2(L))/(L^2-1));
Pebt=(2*(1-1/L)/log2(L))*0.5*erfc(xb/sqrt(2));
Pest=log2(M)*Pebt;
EbNodB=10*log10(EbNo);
EsNodB=10*log10(EsNo);

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


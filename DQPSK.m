%*************************************************************************
%                COMUNICACIONES II (CURSO 2014/2015)
%*************************************************************************
%          Jos� Carlos Segura, Carlos Medina, �ngel de la Torre
%                   {segura, cmedina, atv}@.ugr.es
%*************************************************************************
function [Pes,Peb,Pebt,Pest,EbNodB,EsNodB,sr]=DQPSK(Nbits,sigma,desfase,flag)
%% Argumentos de entrada:
% - Nbits  : N�mero de bits a transmitir
% - sigma  : Valor de desviaci�n t�pica del ruido
% - desfase: desfase de la constelaci�n (radianes)
% - flag   : (0) oculta gr�ficas y resultados
%            (1) muestra gr�ficas y resultados
%
%% Objetivo:
% Esta macro simula un sistema DQPSK as�ncrono.

% Par�metros:
T=1e-3;      % Periodo de simbolo T=1ms
dt=0.025e-3; % Periodo de muestreo dt=0.025ms
Fs=1/dt;     % Fs=1/dt = 40KHz
Nss=T/dt;    % N�mero de muestras por periodo de s�mbolo
Fc=2e3;      % Frecuencia de portadora Fc=2KHz

% Generaci�n de la cadena de bits (Bits): Nbits=6e3
Bits=(randn(1,Nbits)>0);
 
% Valores de los s�mbolos (valor entero m): 
Nbs   = 2;                               % N�mero de bits por s�mbolo 
Nsimb = Nbits/Nbs;                       % N�mero de s�mbolos
Simbs = 2*Bits(1:2:end) + Bits(2:2:end); % valores m={0,1,2,3}

%************************************
%% Modulaci�n DQPSK
%************************************
% Aplicamos codificaci�n de Gray para representar los s�mbolos en 
% el espacio de se�al (usamos la funci�n "bin2gray" de matlab)
M=4;
[x_gray,gray_map]=bin2gray(Simbs,'dpsk',M); % gray_map=(0,1,3,2) => nuevo orden de los s�mbolos  
[tf,index]=ismember(Simbs,gray_map);        % ordenamos los s�mbolos (el 3 cambia por el 2)
x=index-1;                                  % restamos 1 para obtener un �ndice de s�mbolo entre 0 y 3

%% DIFERENCIA DE FASE ASIGNADA A LOS S�MBOLOS EN DQPSK (Sin codificaci�n de Gray)
% NOTA: Mantenemos el orden de s�mbolos sin codificaci�n de Gray, dado que 
% dicha codificaci�n se considera al aplicar bin2gray a la secuencia
% de s�mbolos. Por tanto, las diferencias de fase entre s�mbolos son:
% b1 b2   m   dif_Fase       
% -- --  --- ----------  
% 0   0   0       0  
% 0   1   1     +pi/2 
% 1   0   2     +pi  
% 1   1   3     +3pi/2 

dif_Fase=(x==0)*0+(x==1)*pi/2+(x==2)*pi+(x==3)*3*pi/2;

% Calculamos el vector de fases de la se�al modulada en DQPSK. Comenzamos
% con una fase inicial igual a cero:
yPhase=[];
fase=0;
for k=1:length(dif_Fase)
    fase=fase+dif_Fase(k);
    fase=mod(fase,2*pi);
    yPhase=[yPhase fase];
end

%% Efectos del canal: desfase
yPhase=yPhase+desfase;

% Calculamos la se�al banda-base (consideramos s�lo informaci�n de fase):
s=exp(j*yPhase);

% Nos aseguramos que la salida del modulador de se�al es compleja: 
if isreal(s)
    s=complex(s,0);
end

%% Efectos del canal: a�adimos ruido complejo
sr=s+randn(size(s))*sigma+j*randn(size(s))*sigma; 

%% Representaci�n del diagrama de constelaci�n de se�al y las
%% transiciones entre s�mbolos.
if (flag)
    figure(1);
    if(desfase==0)
      des=sprintf('0');  
    elseif(desfase==pi/4)
      des=sprintf('pi/4');  
    else
      des=sprintf('fase incorrecta');
    end

    plot(sr,'.');
    grid on;
    xlabel('Componente en fase (ai)');
    ylabel('Componente en cuadratura (aq)');
    title(sprintf('Constelaci�n de se�al (modulaci�n DQPSK): desfase:%s',des));
end

%************************************
%% Demodulaci�n DQPSK
%************************************
% Interpretamos el valor de desfase adicional de se�al para dejarlo
% en el intervalo [-pi,+pi]
desfase=mod(desfase,2*pi);
if(desfase>pi)
   desfase=desfase-2*pi;
end

% Recuperamos el valor de fase de los s�mbolos. Tened en cuenta que los 
% s�mbolos son n�meros complejos. El �ngulo puede ser determinado usando 
% la funci�n "angle". Para recuperar la fase de cada s�mbolo hay que
% deshacer el acumulado realizando la diferencia. Para ello tenemos en
% cuenta que la fase inicial es cero:

z=diff(unwrap([0 angle(sr)]));

% Convertimos z al dominio lineal y redondeamos al valor entero m�s
% pr�ximo. Esto recupera el s�mbolo (con codificaci�n de Gray):

Norm_Factor=M/(2*pi);      % Factor de nomalizaci�n para pasar al dominio lineal 
z=ceil(z*Norm_Factor-0.5); % Aplicamos el criterio de redondeo con umbral 0.5

% Recuperamos el valor de los s�mbolos:
z(z<0)=M + z(z<0);
z=z';

% Aplicamos la decodificaci�n de Gray (usamos la funci�n "gray2bin"): 
[z_degray,gray_map]=gray2bin(z,'dpsk',M); 

% Secuencia de s�mbolos recuperada:  
z=gray_map(z+1); 

%******************************************************************
%% Conversi�n de la secuencia de simbolos a una secuencia binaria
%******************************************************************
Rsimbs=[]; % S�mbolos recibidos
Rbits=[];  % Bits recibidos

% Convertimos los s�mbolos (n�meros decimales) a binario y concatenamos 
% la secuencia de bits en un �nico vector:
%% Correspondencia de s�mbolos a bits:
% 0 => 0 0
% 1 => 0 1
% 2 => 1 0
% 3 => 1 1

for n=1:length(z)
    simb=z(n);
    Rsimbs=[Rsimbs simb];
    if(simb==0) bits=[0 0];
    elseif(simb==1) bits=[0 1];
    elseif(simb==2) bits=[1 0];
    else bits=[1 1];
    end
    Rbits=[Rbits bits];
end

%************************************************
%% Estimaci�n de las probabilidades de error
%************************************************
% Comparamos los bits y los s�mbolos recuperados en el receptor con los
% transmistidos. Ello nos permite determinar el n�mero de errores durante 
% la transmisi�n y por tanto la probabilidad de error.
Nes=sum(Simbs~=Rsimbs);
Neb=sum(Bits~=Rbits);

Pes=Nes/length(Simbs);
Peb=Neb/length(Bits);

if(flag)
    fprintf('N�mero de s�mbolos=%d\n',Nsimb);
    fprintf('N�mero de bits=%d\n',Nbits);
    fprintf('N�mero de s�mbolos err�neos=%d\t Pes= %f\n',Nes,Pes);
    fprintf('N�mero de bits err�neos=%d\t\t Peb = %f\n',Neb,Peb);
end


%% Estimaci�n de Eb/No 
% Estimaci�n experimental: Usamos el valor de Pes para despejar el
% argumento de la funci�n Q y con �l el valor de Eb/No
xs=sqrt(2)*erfcinv(Pes);
EsNo=0.5*(xs/sin(pi/(sqrt(2)*4))).^2;
EbNo=EsNo/2;

% Expresiones de la probabilidad de error de bit (Pebt) y de s�mbolo (Pest) te�ricas. 
% En ellas sustituimos el vector EbNo obtenido de forma experimental a partir
% de Peb y Pes experimentales.

%% Pes=2*Q(sqrt(2*Es/No)*sin(pi/(sqrt(2)*M))); con M=4(QPSK) Es=Eb*log2(M)=2*Eb;
xs=sqrt(4*EbNo)*sin(pi/(sqrt(2)*4));
Pest=erfc(xs/sqrt(2));      

%% Peb=Pes/log2(M)=Pes/2;
Pebt=Pest/2;

% Los valores de EsNo y EbNo en decibelios son:
EsNodB=10*log10(EsNo);
EbNodB=10*log10(EbNo);

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


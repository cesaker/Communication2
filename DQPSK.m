%*************************************************************************
%                COMUNICACIONES II (CURSO 2014/2015)
%*************************************************************************
%          José Carlos Segura, Carlos Medina, Ángel de la Torre
%                   {segura, cmedina, atv}@.ugr.es
%*************************************************************************
function [Pes,Peb,Pebt,Pest,EbNodB,EsNodB,sr]=DQPSK(Nbits,sigma,desfase,flag)
%% Argumentos de entrada:
% - Nbits  : Número de bits a transmitir
% - sigma  : Valor de desviación típica del ruido
% - desfase: desfase de la constelación (radianes)
% - flag   : (0) oculta gráficas y resultados
%            (1) muestra gráficas y resultados
%
%% Objetivo:
% Esta macro simula un sistema DQPSK asíncrono.

% Parámetros:
T=1e-3;      % Periodo de simbolo T=1ms
dt=0.025e-3; % Periodo de muestreo dt=0.025ms
Fs=1/dt;     % Fs=1/dt = 40KHz
Nss=T/dt;    % Número de muestras por periodo de símbolo
Fc=2e3;      % Frecuencia de portadora Fc=2KHz

% Generación de la cadena de bits (Bits): Nbits=6e3
Bits=(randn(1,Nbits)>0);
 
% Valores de los símbolos (valor entero m): 
Nbs   = 2;                               % Número de bits por símbolo 
Nsimb = Nbits/Nbs;                       % Número de símbolos
Simbs = 2*Bits(1:2:end) + Bits(2:2:end); % valores m={0,1,2,3}

%************************************
%% Modulación DQPSK
%************************************
% Aplicamos codificación de Gray para representar los símbolos en 
% el espacio de señal (usamos la función "bin2gray" de matlab)
M=4;
[x_gray,gray_map]=bin2gray(Simbs,'dpsk',M); % gray_map=(0,1,3,2) => nuevo orden de los símbolos  
[tf,index]=ismember(Simbs,gray_map);        % ordenamos los símbolos (el 3 cambia por el 2)
x=index-1;                                  % restamos 1 para obtener un índice de símbolo entre 0 y 3

%% DIFERENCIA DE FASE ASIGNADA A LOS SÍMBOLOS EN DQPSK (Sin codificación de Gray)
% NOTA: Mantenemos el orden de símbolos sin codificación de Gray, dado que 
% dicha codificación se considera al aplicar bin2gray a la secuencia
% de símbolos. Por tanto, las diferencias de fase entre símbolos son:
% b1 b2   m   dif_Fase       
% -- --  --- ----------  
% 0   0   0       0  
% 0   1   1     +pi/2 
% 1   0   2     +pi  
% 1   1   3     +3pi/2 

dif_Fase=(x==0)*0+(x==1)*pi/2+(x==2)*pi+(x==3)*3*pi/2;

% Calculamos el vector de fases de la señal modulada en DQPSK. Comenzamos
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

% Calculamos la señal banda-base (consideramos sólo información de fase):
s=exp(j*yPhase);

% Nos aseguramos que la salida del modulador de señal es compleja: 
if isreal(s)
    s=complex(s,0);
end

%% Efectos del canal: añadimos ruido complejo
sr=s+randn(size(s))*sigma+j*randn(size(s))*sigma; 

%% Representación del diagrama de constelación de señal y las
%% transiciones entre símbolos.
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
    title(sprintf('Constelación de señal (modulación DQPSK): desfase:%s',des));
end

%************************************
%% Demodulación DQPSK
%************************************
% Interpretamos el valor de desfase adicional de señal para dejarlo
% en el intervalo [-pi,+pi]
desfase=mod(desfase,2*pi);
if(desfase>pi)
   desfase=desfase-2*pi;
end

% Recuperamos el valor de fase de los símbolos. Tened en cuenta que los 
% símbolos son números complejos. El ángulo puede ser determinado usando 
% la función "angle". Para recuperar la fase de cada símbolo hay que
% deshacer el acumulado realizando la diferencia. Para ello tenemos en
% cuenta que la fase inicial es cero:

z=diff(unwrap([0 angle(sr)]));

% Convertimos z al dominio lineal y redondeamos al valor entero más
% próximo. Esto recupera el símbolo (con codificación de Gray):

Norm_Factor=M/(2*pi);      % Factor de nomalización para pasar al dominio lineal 
z=ceil(z*Norm_Factor-0.5); % Aplicamos el criterio de redondeo con umbral 0.5

% Recuperamos el valor de los símbolos:
z(z<0)=M + z(z<0);
z=z';

% Aplicamos la decodificación de Gray (usamos la función "gray2bin"): 
[z_degray,gray_map]=gray2bin(z,'dpsk',M); 

% Secuencia de símbolos recuperada:  
z=gray_map(z+1); 

%******************************************************************
%% Conversión de la secuencia de simbolos a una secuencia binaria
%******************************************************************
Rsimbs=[]; % Símbolos recibidos
Rbits=[];  % Bits recibidos

% Convertimos los símbolos (números decimales) a binario y concatenamos 
% la secuencia de bits en un único vector:
%% Correspondencia de símbolos a bits:
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
%% Estimación de las probabilidades de error
%************************************************
% Comparamos los bits y los símbolos recuperados en el receptor con los
% transmistidos. Ello nos permite determinar el número de errores durante 
% la transmisión y por tanto la probabilidad de error.
Nes=sum(Simbs~=Rsimbs);
Neb=sum(Bits~=Rbits);

Pes=Nes/length(Simbs);
Peb=Neb/length(Bits);

if(flag)
    fprintf('Número de símbolos=%d\n',Nsimb);
    fprintf('Número de bits=%d\n',Nbits);
    fprintf('Número de símbolos erróneos=%d\t Pes= %f\n',Nes,Pes);
    fprintf('Número de bits erróneos=%d\t\t Peb = %f\n',Neb,Peb);
end


%% Estimación de Eb/No 
% Estimación experimental: Usamos el valor de Pes para despejar el
% argumento de la función Q y con él el valor de Eb/No
xs=sqrt(2)*erfcinv(Pes);
EsNo=0.5*(xs/sin(pi/(sqrt(2)*4))).^2;
EbNo=EsNo/2;

% Expresiones de la probabilidad de error de bit (Pebt) y de símbolo (Pest) teóricas. 
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
    fprintf('Comprobación (Teoria <-> Experimental)\n');
    fprintf('------------------------------------------\n');
    fprintf('Datos:\t sigma=%1.3f\t Nbits=%d\n',sigma,Nbits);
    fprintf('      \t EbNodB=%1.3f\t EsNodB=%1.3f\n',EbNodB,EsNodB); 
    fprintf('Prob:\t Peb=%1.4f\t\t Pebt=%1.4f\n',Peb,Pebt); 
    fprintf('     \t Pes=%1.4f\t\t Pest=%1.4f\n',Pes,Pest);
    fprintf('-------------------------------------------\n');
end

end


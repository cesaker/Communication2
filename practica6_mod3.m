%*************************************************************************
%                COMUNICACIONES II (CURSO 2014/2015)
%*************************************************************************
% -------------------------------------------------------------------------
%          Jos� Carlos Segura, Carlos Medina, �ngel de la Torre
%                   {segura, cmedina, atv}@.ugr.es
% -------------------------------------------------------------------------
clear all;
close all;
clc;
% *************************************************************************
% 1) Definimos los par�metros de la simulaci�n ----------------------------
% *************************************************************************
SNRdB=1:0.5:9;                          % Vector de SNR en dB
SNR=10.^(SNRdB./10);                     % Vector de SNR en escala lineal  
Npalabras=10000;                          % N�mero de palabras de informaci�n  
k=4;                                     % N�mero de bits del mensaje  
n=7;                                     % N�mero de bits de la palabra de c�digo
BER_soft=zeros(length(SNR),1);           % Iniciamos el vector de BER para decodificaci�n SOFT 
BER_hard=zeros(length(SNR),1);           % Iniciamos el vector de BER para decodificaci�n HARD 
BER_bit=zeros(length(SNR),1);            % Iniciamos el vector de BER sin codificaci�n 
mensajes=floor(2*rand(k,Npalabras));     % Generaci�n de los mensajes

% *************************************************************************
% 2) Definimos los bits de paridad ----------------------------------------
% *************************************************************************
b0=xor(mensajes(1,:),xor(mensajes(3,:),mensajes(4,:)));   % Primer bit de paridad
b1=xor(mensajes(1,:),xor(mensajes(2,:),mensajes(3,:)));   % Segundo bit de paridad
b2=xor(mensajes(2,:),xor(mensajes(3,:),mensajes(4,:)));   % Tercer bit de paridad

% *************************************************************************
% 3) Codificamos las palabras de informaci�n  -----------------------------
% *************************************************************************
palabras_codigo=[mensajes;b0;b1;b2];      % A�adimos los bits de paridad
pcodbin=palabras_codigo;                  % Guardamos las palabras originales
palabras_codigo(palabras_codigo==0)=-1;   % Aplicamos codificaci�n bipolar (+1,-1)

% *************************************************************************
% 3) Definimos los par�metros de la decodificaci�n  -----------------------
% *************************************************************************
bits_decodificados_hard=zeros(n,Npalabras);     % Inicializamos el vector de bits decodificados (HARD).  
bits_decodificados_soft=zeros(n,Npalabras);     % Inicializamos el vector de bits decodificados (SOFT).  
H=[1 0 1 1 1 0 0;1 1 1 0 0 1 0;0 1 1 1 0 0 1]'; % Construimos la matriz de paridad.

C=de2bi((0:2^(k)-1),4);                     % Construimos la tabla de 2^k mensajes posibles.
C(1:16,5)=xor(C(:,1),xor(C(:,3),C(:,4)));   % A�adimos el primer bit de paridad para completar la tabla.
C(1:16,6)=xor(C(:,1),xor(C(:,2),C(:,3)));   % A�adimos el segundo bit de paridad para completar la tabla.
C(1:16,7)=xor(C(:,2),xor(C(:,3),C(:,4)));   % A�adimos el tercer bit de paridad para completar la tabla.

Cp=2*C-1; % Mensajes bipolares (-1,+1)

distancia=zeros(1,2^k);   % Inicializamos el vector de distancias (decodificaci�n SOFT). 

% *************************************************************************
% 4) Iniciamos la decodificaci�n iterando los valores de SNR  -------------
% *************************************************************************
for i=1:length(SNR)
    
    y=palabras_codigo+randn(n,Npalabras)/sqrt(SNR(i));     % A�adimos ruido para simular la transmisi�n por el canal.
    
    % 4.1) Decodificaci�n HARD --------------------------------------------
    bits_decodificados_hard(y>0)=1;   % Todos los bits positivos recibidos son convertidos a 1
    bits_decodificados_hard(y<0)=0;   % Todos los bits negativos recibidos son convertidos a 0
    
    % 4.1.0) BER de los bits sin codificar ----------------------------------
    BER_bit(i)=length(find(bits_decodificados_hard~=pcodbin));
    
    % 4.1.1) Detecci�n y correcci�n de errores HARD -----------------------
    for l=1:Npalabras
        % Calculamos el s�ndrome:
        sindrome=mod(bits_decodificados_hard(:,l)'*H,2);  
        
        % Buscamos el s�ndrome calculado en la matriz de comprobaci�n de paridad 
        %(que es equivalente a buscar en la matriz de s�ndromes posibles):
        for j=1:n
            % Si encontramos el s�ndrome, corregimos el bit correspondiente (neg�ndolo):
            if (sindrome==H(j,:))
                bits_decodificados_hard(j,l)=~bits_decodificados_hard(j,l);
            end
        end
    end
    % Calculamos el BER:
    BER_hard(i)=length(find(bits_decodificados_hard(1:4,:)~=mensajes));  
    
    % 4.2) Decodificaci�n SOFT --------------------------------------------
    for l=1:Npalabras
        % Para cada una de las palabras recibidas calculamos su distancia
        % respecto a cada una de las palabras en la tabla de mensajes
        % posibles (C):
        for m=1:(2^k)           
            distancia(m)=sqrt(sum((y(:,l)-Cp(m,:)').^2));
        end
        % Buscamos la m�nima distancia y su �ndice asociado:
        [minval,minind]=min(distancia);  
        
        % Escogemos como palabra decodificada la palabra cuya distancia es
        % m�nima respecto a la palabra recibida:
        bits_decodificados_soft(:,l)=C(minind,:);      
    end
    % Calculamos el BER:
    BER_soft(i)=length(find(bits_decodificados_soft(1:4,:)~=mensajes));   
    
end

% *************************************************************************
% 4) Representamos gr�ficamente los BER  ----------------------------------
% *************************************************************************
BER_soft=BER_soft/(k*Npalabras); % S�lo bits decodificados
BER_hard=BER_hard/(k*Npalabras); % S�lo bits decodificados
BER_bit=BER_bit/(n*Npalabras);   % Todos los bits (no hay codificaci�n)

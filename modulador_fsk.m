%*************************************************************************
%                COMUNICACIONES II (CURSO 2014/2015)
%*************************************************************************
% -------------------------------------------------------------------------
%          José Carlos Segura, Carlos Medina, Ángel de la Torre
%                   {segura, cmedina, atv}@.ugr.es
% -------------------------------------------------------------------------

function [s_mod]=modulador_fsk(x,amplitud,f1,f2,Rb,Nss,Es,mostrar,Nmostrar)

% FUNCTION MODULADOR_FSK modula una señal binaria utilizando modulación
% FSK.
% Parámetros de entrada:
% |_ x (vector): Señal binaria.
% |_ amplitud (escalar): amplitud de la codificación.
% |_ f1:  Frecuencia de la primera portadora.
% |_ f2:  Frecuencia de la segunda portadora.
% |_ Rb:  Tasa de transmisión (bits/s).
% |_ Nss: Número de muestras por símbolo.
% |_ Es:  Energía de símbolo.
% |_ mostrar: Flag para visualizar gráficas. 
% |_ Nmostrar: Numero de periodos a mostrar en las gráficas.
%
% Parámetros de salida:
% |_ s_mod (vector): Señal modulada.

% *************************************************************************
% 1) Codificamos la señal binaria utilizando codificación unipolar NRZ
% *************************************************************************

% Separamos la secuencia de símbolos "1": 
x_cod1=linecodeunipolar(x,amplitud,Nss);
x_cod1(1+length(x_cod1))=x_cod1(length(x_cod1));

% Separamos la secuencia de símbolos "0":
% Para ello calculamos el negativo del vector, es decir, ponemos un "0" 
% donde hay 1's y un "1" donde hay 0's  
xnot=not(x); 
x_cod2=linecodeunipolar(xnot,amplitud,Nss);
x_cod2(1+length(x_cod2))=x_cod2(length(x_cod2));

if mostrar
    figure
    subplot(2,2,1)
    stem(x(1:10),'r')
    xlabel('Nº de símbolo','FontSize',12)
    ylabel('Amplitud (V)','FontSize',12)
    legend('Secuencia original')
    grid on
    subplot(2,2,2)
    stem(xnot(1:10),'r')
    xlabel('Nº de símbolo','FontSize',12)
    ylabel('Amplitud (V)','FontSize',12)
    legend('Secuencia original negada')
    grid on
    subplot(2,2,3)
    plot(x_cod1(1:Nmostrar*Nss),'LineWidth',2)
    xlabel('Nº de muestra','FontSize',12)
    ylabel('Amplitud (V)','FontSize',12)
    legend('Sec. codificada')
    grid on
    subplot(2,2,4)
    plot(x_cod2(1:Nmostrar*Nss),'LineWidth',2)
    xlabel('Nº de muestra','FontSize',12)
    ylabel('Amplitud (V)','FontSize',12)
    legend('Sec. codificada negada')
    grid on
end

% *************************************************************************
% 2) Generamos las portadoras 
% *************************************************************************
t=0:2/Nss:2;
T_simbolo=Nss/Rb;
portadora1=sqrt(2*Es/T_simbolo)*sin(pi*f1*t);
portadora2=sqrt(2*Es/T_simbolo)*sin(pi*f2*t);

if mostrar
    figure
    subplot(2,1,1)
    plot(portadora1,'LineWidth',2)
    legend('Portadora 1')
    xlabel('Nº de muestra','FontSize',12)
    ylabel('Amplitud (V)','FontSize',12)
    axis([0 100 min(portadora1) max(portadora1)]);
    grid on
    subplot(2,1,2)
    plot(portadora2,'b','LineWidth',2)
    legend('Portadora 2')
    xlabel('Nº de muestra','FontSize',12)
    ylabel('Amplitud (V)','FontSize',12)
    axis([0 100 min(portadora2) max(portadora2)]);
    grid on
end

% *************************************************************************
% 3) Modulamos la señal
% *************************************************************************
j=1;
k=Nss+1;
for n=1:length(x)
    fsk1(j:k)=x_cod1(j:k).*portadora1;
    j=j+Nss;
    k=k+Nss;
end
j=1;
k=Nss+1;
for n=1:length(xnot)
    fsk2(j:k)=x_cod2(j:k).*portadora2;
    j=j+Nss;
    k=k+Nss;
end
s_mod=fsk1+fsk2;

if mostrar
    figure
    subplot(3,1,1)
    plot(fsk1(1:Nmostrar*Nss),'LineWidth',2)
    hold on
    plot(x_cod1(1:Nmostrar*Nss),'r','LineWidth',2)
    legend('Señal modulada (portadora 1)','Sec. cod.')
    xlabel('Nº de muestra','FontSize',12)
    ylabel('Amplitud (V)','FontSize',12)
    title('Modulación de los símbolos 1')
    grid on
    subplot(3,1,2)
    plot(fsk2(1:Nmostrar*Nss),'LineWidth',2)
    hold on
    plot(x_cod2(1:Nmostrar*Nss),'r','LineWidth',2)
    legend('Señal modulada (portadora 2)','Sec. cod. negada')
    xlabel('Nº de muestra','FontSize',12)
    ylabel('Amplitud (V)','FontSize',12)
    title('Modulación de los símbolos 0')
    grid on
    subplot(3,1,3)
    plot(s_mod(1:Nmostrar*Nss),'LineWidth',2)  
    legend('Señal modulada (FSK)')
    xlabel('Nº de muestra','FontSize',12)
    ylabel('Amplitud (V)','FontSize',12)
    title('Señal modulada con FSK')
    grid on
end



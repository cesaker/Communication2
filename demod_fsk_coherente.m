%*************************************************************************
%                COMUNICACIONES II (CURSO 2014/2015)
%*************************************************************************
% -------------------------------------------------------------------------
%          Jos� Carlos Segura, Carlos Medina, �ngel de la Torre
%                   {segura, cmedina, atv}@.ugr.es
% -------------------------------------------------------------------------
function [xR]=demod_fsk_coherente(sr,f1,f2,Rb,Nss,Es,desfase,mostrar,Nmostrar)

% FUNCTION DEMOD_FSK_COHERENTE demodula una se�al modulada en FSK binaria
% de forma coherente.
%
% Par�metros de entrada:
% |_ sr (vector): Se�al modulada.
% |_ f1:  Frecuencia de la primera portadora generada localmente.
% |_ f2:  Frecuencia de la segunda portadora generada localmente.
% |_ Rb:  Tasa de transmisi�n (bits/s).
% |_ Nss: N�mero de muestras por s�mbolo.
% |_ Es:  Energ�a de s�mbolo.
% |_ desfase: Desfase de las portadoras locales (idealmente es nulo).
% |_ mostrar: Flag para visualizar gr�ficas. 
% |_ Nmostrar: Numero de periodos a mostrar en las gr�ficas.
%
% Par�metros de salida:
% |_ xR (vector): Se�al binaria recuperada.

% *************************************************************************
% 1) Generamos las portadoras locales
% *************************************************************************
t=0:2/Nss:2;
T_simbolo=Nss/Rb;
portadora1=(sqrt(2*Es/T_simbolo))*sin(pi*f1*t+desfase);
portadora2=(sqrt(2*Es/T_simbolo))*sin(pi*f2*t+desfase);

% *************************************************************************
% 2) Multiplicamos la se�al recibida por las portadoras locales.
% *************************************************************************
o=1;
p=Nss+1;
for n=1:round(length(sr)/Nss)
sr1(o:p)=portadora1.*sr(o:p);
sr2(o:p)=portadora2.*sr(o:p);
o=o+Nss;
p=p+Nss;
end

if mostrar
    figure
    subplot(2,1,1)
    plot(sr1(1:Nmostrar*Nss))
    xlabel('N� de muestra','FontSize',12)
    ylabel('Amplitud (V)','FontSize',12)
    legend('Se�al demodulada (canal de los 1s)')
    title('Demodulaci�n FSK coherente')
    grid on
    subplot(2,1,2)
    plot(sr2(1:Nmostrar*Nss))
    xlabel('N� de muestra','FontSize',12)
    ylabel('Amplitud (V)','FontSize',12)
    legend('Se�al demodulada (canal de los 0s)')
    grid on
end

% *************************************************************************
% 3) Muestreamos la se�al y la pasamos por el decisor.
% *************************************************************************
u=1;
v=Nss+1;
xR=[];
for s=1:round(length(sr)/Nss)
    sum_sr1=sum(sr1(u:v));
    sum_sr2=sum(sr2(u:v));
    sumco=sum_sr1-sum_sr2;
    s_demod(u:v)=sumco;
    if sumco>0;
        xR(s)=1;
    else
        xR(s)=0;
    end
    u=u+Nss;
    v=v+Nss;
end

if mostrar
    figure
    subplot(2,1,1)
    plot(s_demod(1:Nmostrar*Nss));
    hold on;
    plot(zeros(1,Nmostrar*Nss),'r','LineWidth',2);
    xlabel('N� de muestra','FontSize',12)
    ylabel('Amplitud (V)','FontSize',12)
    legend('Se�al codificada recuperada');
    title('Codificaci�n de se�al y salida del decisor')
    grid on
    subplot(2,1,2);
    stem(xR(1:Nmostrar),'b','fill');
    xlabel('N� de muestra','FontSize',12)
    ylabel('Amplitud (V)','FontSize',12)
    legend('Se�al binaria recuperada');
    grid on
end




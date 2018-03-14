%*************************************************************************
%                COMUNICACIONES II (CURSO 2014/2015)
%*************************************************************************
% -------------------------------------------------------------------------
%          José Carlos Segura, Carlos Medina, Ángel de la Torre
%                   {segura, cmedina, atv}@.ugr.es
% -------------------------------------------------------------------------

% *************************************************************************
% 1) Definimos los parámetros de modulación -------------------------------
% *************************************************************************
M=2;            % Número de símbolos.
amplitud=1;     % Amplitud de la codificación unipolar.
f1=3;           % Frecuencia de la primera portadora.
f2=7;           % Frecuencia de la segunda portadora.
Rb=1;           % Tasa de bits/s.
Fs=1e2;         % Frecuencia de muestreo. 
Nss=Fs/Rb;      % Número de muestras por símbolo.
Es=10;          % Energía de símbolo.
mostrar=0;      % Flag para mostrar gráficas (poner a cero para no mostrar).
Nmostrar=10;    % Número de símbolos/bits a mostrar en las gráficas.

% *************************************************************************
% 2) Generamos la cadena de bits que componen el mensaje a transmitir -----
% *************************************************************************
Nbits=3e4;
x = randi([0 M-1],Nbits,1);

% *************************************************************************
% 3) Modulamos la señal usando FSK ----------------------------------------
s_mod = modulador_fsk(x,amplitud,f1,f2,Rb,Nss,Es,mostrar,Nmostrar);
% *************************************************************************

% *************************************************************************
% 4) Simulamos la transmisión de la señal por un canal AWGN ---------------
% *************************************************************************
i=1;
for sigma=0.6:0.2:3.0
    
sr_mod = s_mod + randn(size(s_mod))*sigma;
if mostrar
    figure
    subplot(2,1,1)
    plot(s_mod(1:Nmostrar*Nss))
    legend('Señal modulada (FSK) a la entrada del canal')
    xlabel('Nº de muestra','FontSize',12)
    ylabel('Amplitud (V)','FontSize',12)
    subplot(2,1,2)
    plot(sr_mod(1:Nmostrar*Nss))
    legend('Señal modulada (FSK) a la salida del canal')
    xlabel('Nº de muestra','FontSize',12)
    ylabel('Amplitud (V)','FontSize',12)
end

% *************************************************************************
% 5) Demodulamos la señal y calculamos el BER -----------------------------
% *************************************************************************
% 5.1) Demodulación coherente 
desfase=pi/2;  % Desfase de la portadora local.
xR_coh=demod_fsk_coherente(sr_mod,f1,f2,Rb,Nss,Es,desfase,mostrar,Nmostrar);

fprintf('Desfase = %.2f x pi\n',desfase/pi);
fprintf('Sigma   = %.2f\n',sigma);

error=sum(xor(x,xR_coh'));
BER_coh(i)=error/length(xR_coh);
fprintf('BER   coh:%0.3e\n',BER_coh);

% 5.2) Demodulación no coherente
xR_nocoh = demod_fsk_nocoherente(sr_mod,f1,f2,Fs,Nss,mostrar,Nmostrar);

error=sum(xor(x,xR_nocoh'));
BER_nocoh(i)=error/length(xR_nocoh);
fprintf('BER nocoh:%0.3e\n',BER_nocoh);
i=i+1;
end
sigma=[0.6:0.2:3.0];
semilogy(sigma,BER_coh, 'r');
hold on
semilogy(sigma,BER_nocoh,'b');
legend('COH', 'NOCOH');
hold off



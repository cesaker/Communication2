%*************************************************************************
%                COMUNICACIONES II (CURSO 2014/2015)
%*************************************************************************
% -------------------------------------------------------------------------
%          Jos� Carlos Segura, Carlos Medina, �ngel de la Torre
%                   {segura, cmedina, atv}@.ugr.es
% -------------------------------------------------------------------------

% *************************************************************************
% 1) Definimos los par�metros de modulaci�n -------------------------------
% *************************************************************************
M=2;            % N�mero de s�mbolos.
amplitud=1;     % Amplitud de la codificaci�n unipolar.
f1=3;           % Frecuencia de la primera portadora.
f2=7;           % Frecuencia de la segunda portadora.
Rb=1;           % Tasa de bits/s.
Fs=1e2;         % Frecuencia de muestreo. 
Nss=Fs/Rb;      % N�mero de muestras por s�mbolo.
Es=10;          % Energ�a de s�mbolo.
mostrar=0;      % Flag para mostrar gr�ficas (poner a cero para no mostrar).
Nmostrar=10;    % N�mero de s�mbolos/bits a mostrar en las gr�ficas.

% *************************************************************************
% 2) Generamos la cadena de bits que componen el mensaje a transmitir -----
% *************************************************************************
Nbits=3e4;
x = randi([0 M-1],Nbits,1);

% *************************************************************************
% 3) Modulamos la se�al usando FSK ----------------------------------------
s_mod = modulador_fsk(x,amplitud,f1,f2,Rb,Nss,Es,mostrar,Nmostrar);
% *************************************************************************

% *************************************************************************
% 4) Simulamos la transmisi�n de la se�al por un canal AWGN ---------------
% *************************************************************************
i=1;
for sigma=0.6:0.2:3.0
    
sr_mod = s_mod + randn(size(s_mod))*sigma;
if mostrar
    figure
    subplot(2,1,1)
    plot(s_mod(1:Nmostrar*Nss))
    legend('Se�al modulada (FSK) a la entrada del canal')
    xlabel('N� de muestra','FontSize',12)
    ylabel('Amplitud (V)','FontSize',12)
    subplot(2,1,2)
    plot(sr_mod(1:Nmostrar*Nss))
    legend('Se�al modulada (FSK) a la salida del canal')
    xlabel('N� de muestra','FontSize',12)
    ylabel('Amplitud (V)','FontSize',12)
end

% *************************************************************************
% 5) Demodulamos la se�al y calculamos el BER -----------------------------
% *************************************************************************
% 5.1) Demodulaci�n coherente 
desfase=pi/2;  % Desfase de la portadora local.
xR_coh=demod_fsk_coherente(sr_mod,f1,f2,Rb,Nss,Es,desfase,mostrar,Nmostrar);

fprintf('Desfase = %.2f x pi\n',desfase/pi);
fprintf('Sigma   = %.2f\n',sigma);

error=sum(xor(x,xR_coh'));
BER_coh(i)=error/length(xR_coh);
fprintf('BER   coh:%0.3e\n',BER_coh);

% 5.2) Demodulaci�n no coherente
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



[Pes,Peb,Pebt,Pest,EbNodB,EsNodB,sr]=QPSK(10e3,0,0,0);
dt=0.025e-3; % Periodo de muestreo dt=0.025ms
Fs=1/dt; 
[Pyy,f]=psd(sr,Fs,length(sr)/4)
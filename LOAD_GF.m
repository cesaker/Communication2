%*************************************************************************
%                  COMUNICACIONES II (CURSO 2014/2015)
%           MACRO: CARGAR FORMA DE ONDA ARBITRARIA EN GF-857
%*************************************************************************
%         José Carlos Segura, Carlos Medina, Ángel de la Torre
%                   {segura, cmedina, atv}@.ugr.es
%*************************************************************************
% function LOAD_GF(x,Fs)
%
% Argumentos de entrada:
% - x  : Vector de muestras de señal que queremos cargar en el generador
% - Fs : Frecuencia de muestreo (Fs=Rb*kb) 
%*************************************************************************
function LOAD_GF(x,Fs)

% Variables:
Nblock=4000;
 
% Configuración del puerto serie:
s1=serial('COM1','BaudRate',19200,'Parity','none','StopBits',1,'ByteOrder','littleEndian');
set(s1,'Timeout',5);
set(s1,'Terminator','CR');
set(s1,'Flowcontrol','none');
set(s1,'OutputBufferSize',2*Nblock);
ter=char(10);
fopen(s1);

% Escalado de señal al rango [-1:1]:
xmax=max(abs(x));
if(xmax>1)
    x=x/xmax;
    fprintf('Valor x re-escalado %f\n',1/xmax);
end

% Conversión 
d=ceil(2047*x);
d(d<0)=32768-d(d<0);
b=[];

for k=1:length(d)
    vv=d(k);
    hb=floor(vv/256);
    lb=vv-256*hb;
    hb=uint8(hb);
    lb=uint8(2*floor(lb/2)+1);
    b=[b ; lb ; hb];
end
b=uint8(b);

% Configura el GF a modo arbitrario (ARB):
s='SOUR:FUNC:ARB';
fwrite(s1,[s ter],'sync');
pause(2);

% Establece la frecuencia de onda sintetizada (Fs):
s=sprintf('SOUR:FREQ:SYNT %.2f',Fs);
fwrite(s1,[s ter]);
pause(2);

% Comprueba la frecuencia de la onda sintetizada:
s=sprintf('SOUR:FREQ:SYNT?');
fwrite(s1,[s ter]);
v=fscanf(s1);
if(length(v)<1)
    fprintf('(2) No puedo contactar con el generador\n');
    fclose(s1);
    return
end

if(v(1)==1)
    Fsr=sscanf(v(2:end),'%f');
    fprintf('Frecuencia fijada (Hz): %10.2f\n',Fsr);
end

% Inicia la carga de datos:
Nw=length(b);
s=sprintf('SOUR:FUNC:LDWF %d',Nw/2);
fwrite(s1,[s ter],'sync');
pause(1);

% Comprueba el estado del equipo antes de escribir los datos:
n=0;
while(s1.bytesavailable<3)
    pause(0.1);
    n=n+1;
    if(n>8)
        fclose(s1);
        fprintf('(1) No puedo cargar los datos\n');
        return
    end
end
v=fread(s1,3,'char');

if(isempty(v))
    fprintf('(2) No puedo enviar los datos\n');
    fclose(s1);
    return
end

if(v(1)~=1)
    fprintf('(3) No puedo enviar los datos\n');
    fclose(s1);
    return
end

% Escribe los datos usando bloques de N bytes:
set(s1,'Timeout',5);
kfin=0;
Nblock=4000;

while(1)
    kini=kfin+1;
    kfin=min(kini+Nblock-1,Nw);
    nw=kfin-kini+1;
    if(nw<1),break;end
    fwrite(s1,b(kini:kfin),'sync');
    pause(0.1);
end
fwrite(s1,ter);
pause(0.5);

% % Pone nuevamente el GF en modo arbitrario(ARB):
s='SOUR:FUNC:ARB';
fwrite(s1,[s ter]);

% Cierra el puerto y termina:
fclose(s1);


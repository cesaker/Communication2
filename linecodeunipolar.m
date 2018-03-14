%*************************************************************************
%                COMUNICACIONES II (CURSO 2014/2015)
%*************************************************************************
% -------------------------------------------------------------------------
%          José Carlos Segura, Carlos Medina, Ángel de la Torre
%                   {segura, cmedina, atv}@.ugr.es
% -------------------------------------------------------------------------
% PARÁMETROS DE ENTRADA:
% |_ binary: Señal binaria (x).
% |_ A     : amplitud de la codificación (amplitud).
% |_ nos   : Número de muestras por símbolo (Nss).
%
% PARÁMETROS DE SALIDA:
% |_ bb: señal codificada usando código de línea UNIPOLAR NRZ y forma de
%        de pulso rectangular. 

function [bb]=linecodeunipolar(binary,A,nos)

l=1;
k=nos;
bb=zeros(1,nos*length(binary));

for n=1:length(binary)
    bit=binary(n);
    if bit==1
        b=A*ones(1,nos);
        bb(l:k)=b;
    else
        b=zeros(1,nos);
        bb(l:k)=b;
    end
    l=l+nos;
    k=k+nos;
end

end
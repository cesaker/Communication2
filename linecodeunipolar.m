%*************************************************************************
%                COMUNICACIONES II (CURSO 2014/2015)
%*************************************************************************
% -------------------------------------------------------------------------
%          Jos� Carlos Segura, Carlos Medina, �ngel de la Torre
%                   {segura, cmedina, atv}@.ugr.es
% -------------------------------------------------------------------------
% PAR�METROS DE ENTRADA:
% |_ binary: Se�al binaria (x).
% |_ A     : amplitud de la codificaci�n (amplitud).
% |_ nos   : N�mero de muestras por s�mbolo (Nss).
%
% PAR�METROS DE SALIDA:
% |_ bb: se�al codificada usando c�digo de l�nea UNIPOLAR NRZ y forma de
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
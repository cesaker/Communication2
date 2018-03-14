%*************************************************************************
%                COMUNICACIONES II (CURSO 2014/2015)
%*************************************************************************
%          José Carlos Segura, Carlos Medina, Ángel de la Torre
%                   {segura, cmedina, atv}@.ugr.es
%*************************************************************************
% Esta función recibe como argmento el nombre del código de línea
%*************************************************************************
% function line_code_mod(code)
%
% Argumentos de entrada:
% code : ['unipolar_nrz'] ['polar_nrz'] ['unipolar_rz']  
%        ['bipolar_rz']   ['manchester_nrz'] 
%*************************************************************************
function line_code(code)

%% Generamos una secuencia aleatoria de 10 bits:
%%Nbits=10;
%%x=randn(Nbits,1)>0;

%% Podemos introducir también una secuencia binaria concreta (p.e):
x=[0 0 0 1 1 1 1 0 0 1];

switch code
    case 'unipolar_nrz'
        hold off;
        n=1;
        x(11)=1;
        while n<=10;
            t=n-1:0.001:n;
            if x(n) == 0
                if x(n+1)==0  
                    y=(t>n);
                else
                    y=(t==n);
                end
            else
                if x(n+1)==0
                    y=(t<n)-0*(t==n);
                else
                    y=(t<n)+1*(t==n);
                end
            end
            figure(1);
            d=plot(t,y);title('Código UNIPOLAR-NRZ');grid on;
            set(d,'LineWidth',2.5);
            xlabel('Tiempo (ms)');
            ylabel('Amplitud (V)');
            hold on;
            axis([0 10 -1.5 1.5]);
            n=n+1;
        end
    case 'polar_nrz'
        hold off;
        n=1;
        x(11)=1;
        while n<=10;
            t=n-1:0.001:n;
            if x(n) == 0
                if x(n+1)==0  
                    y=-(t<n)-(t==n);
                else
                    y=-(t<n)+(t==n);
                end
            else
                if x(n+1)==0
                    y=(t<n)-1*(t==n);
                else
                    y=(t<n)+1*(t==n);
                end
            end
            figure(1);
            d=plot(t,y);title('Código POLAR-NRZ');grid on
            set(d,'LineWidth',2.5); 
            xlabel('Tiempo (ms)');
            ylabel('Amplitud (V)');
            hold on;
            axis([0 10 -1.5 1.5]);
            n=n+1;
        end
        
    case 'unipolar_rz'
        hold off;
        n=1;
        x(11)=1;
        while n<=10;
            t=n-1:0.001:n;
            
            %Graficación de los CEROS (0)
            if x(n) == 0
                if x(n+1)==0  
                    y=(t>n);
                else
                    y=(t==n);
                end
          
            %Graficación de los UNOS (1)
            else
                if x(n+1)==0
                    y=(t<n-0.5);
                else
                    y=(t<n-0.5)+1*(t==n);
                end
            end
            figure(1);
            d=plot(t,y);title('Código UNIPOLAR-RZ');grid on
            set(d,'LineWidth',2.5);
            xlabel('Tiempo (ms)');
            ylabel('Amplitud (V)');
            hold on;
            axis([0 10 -1.5 1.5]);
            n=n+1;
        end
        
    case 'bipolar_rz'
        hold off;
        n=1;
        m=1;
        x(11)=1;
     
        while n<=10;
            t=n-1:0.001:n;
            
            %Graficación de los CEROS (0)
            if x(n)==0
                if x(n+1)==0  
                    y=(t>n);
                else
                    y=(m)*(t==n);
                end
          
            %Graficación de los UNOS (1)
            else
                if x(n+1)==0
                    y=(m)*(t<n-0.5);
                else
                    y=(m)*(t<n-0.5)+(-m)*(t==n);
                end
                m=-m;
            end
            
            figure(1);
            d=plot(t,y);title('Código BIPOLAR-RZ');grid on
            set(d,'LineWidth',2.5);
            xlabel('Tiempo (ms)');
            ylabel('Amplitud (V)');
            hold on;
            axis([0 10 -1.5 1.5]);
            n=n+1;
        end

    case 'manchester_nrz'
        hold off;
        x=~x;
        n=1;
        x(11)=1;
        
        while n<=10;
            t=n-1:0.001:n;
            if x(n)==0
                if x(n+1)==0  
                    y=-(t<n)+2*(t<n-0.5)+1*(t==n);
                else
                    y=-(t<n)+2*(t<n-0.5)-1*(t==n);
                end
            else
                if x(n+1)==0
                    y=(t<n)-2*(t<n-0.5)+1*(t==n);
                else
                    y=(t<n)-2*(t<n-0.5)-1*(t==n);
                end
            end
            figure(1);
            d=plot(t,y);title('Código MANCHESTER-NRZ');grid on
            set(d,'LineWidth',2.5);
            xlabel('Tiempo (ms)');
            ylabel('Amplitud (V)');
            hold on;
            axis([0 10 -1.5 1.5]);
            n=n+1;
        end
end

% Mostramos en la parte superior de la figura la secuencia binaria codificada:
figure(1);
str=num2str(x(1:end-1)); 
text(0.5,1.25,sprintf('Secuencia binaria: %s',str),'FontSize',13);

%*****************************************************************************
% Representa la PSD teórica del código de línea seleccionado
%*****************************************************************************

switch code
    
    case 'unipolar_nrz'
        hold off;
        A=sqrt(2);
        Tb=1.0e-3;
        R=1/Tb;
        L=2*R;
        f=0:L/100:L;
        del=0;
        P=(A.^2*Tb)/4*(sinc(f*Tb)).^2*(1+(1/Tb)*del);
        figure(2);
        g=plot(f,P);
        title('PSD UNIPOLAR-NRZ');
        hold on;xlabel('Frecuencia (Hz)');ylabel('Densidad espectral de potencia (W/Hz)');
        axis([0 L 0 1.1*Tb]);set(g,'LineWidth',2.5);
        stem(0,(A.^2*Tb)/2,'LineWidth',2.5);hold off;
        axis([0 L 0 1.04*Tb]);set(g,'LineWidth',2.5);
        set(gca,'XTickMode','manual','XTick',[0,0.5*R,R,1.5*R,2*R]);grid on;
        set(gca,'YTickMode','manual','YTick',[0,0.25*Tb,0.5*Tb,0.75*Tb,Tb]);
        set(gca,'XTickLabel',{['0'];['0.5R'];['R'];['1.5R'];['2R']});
        set(gca,'YTickLabel',{['0'];['0.25Tb'];['0.5Tb'];['0.75Tb'];['Tb']});
        
    case 'polar_nrz'
        hold off;
        A=1;
        Tb=1.0e-3;
        R=1/Tb;
        L=2*R;
        f=0:L/100:L;
        del=0;
        P=(A.^2*Tb)*(sinc(f*Tb)).^2;
        figure(2);
        g=plot(f,P);
        title('PSD POLAR-NRZ');
        hold on;xlabel('Frecuencia (Hz)');ylabel('Densidad espectral de potencia (W/Hz)');
        axis([0 L 0 1.05*Tb]);set(g,'LineWidth',2.5);
        set(gca,'XTickMode','manual','XTick',[0,0.5*R,R,1.5*R,2*R]);grid on;
        set(gca,'YTickMode','manual','YTick',[0,0.25*Tb,0.5*Tb,0.75*Tb,Tb]);
        set(gca,'XTickLabel',{['0'];['0.5R'];['R'];['1.5R'];['2R']});
        set(gca,'YTickLabel',{['0'];['0.25Tb'];['0.5Tb'];['0.75Tb'];['Tb']});
        
    case 'unipolar_rz'    
        hold off;
        A=2;
        Tb=1.0e-3;
        R=1/Tb;
        L=2*R;
        f=0:L/100:L;
        P=(A.^2*Tb)/16*(sinc(f*Tb/2)).^2;
        figure(2);
        g=plot(f,P);
        title('PSD UNIPOLAR-RZ');
        hold on;xlabel('Frecuencia (Hz)');ylabel('Densidad espectral de potencia (W/Hz)');
        axis([0 L 0 1.05*Tb]);set(g,'LineWidth',2.5);
        stem([0 R],[(A*Tb)/4 (A*Tb)/10],'LineWidth',2.5);hold off;
        axis([0 L 0 1.04*Tb]);set(g,'LineWidth',2.5);
        set(gca,'XTickMode','manual','XTick',[0,0.5*R,R,1.5*R,2*R]);grid on;
        set(gca,'YTickMode','manual','YTick',[0,0.25*Tb,0.5*Tb,0.75*Tb,Tb]);
        set(gca,'XTickLabel',{['0'];['0.5R'];['R'];['1.5R'];['2R']});
        set(gca,'YTickLabel',{['0'];['0.25Tb'];['0.5Tb'];['0.75Tb'];['Tb']});
        
     case 'bipolar_rz'
        hold off;
        A=2;
        Tb=1.0e-3;
        R=1/Tb;
        L=2*R;
        f=0:L/100:L;
        P=(A.^2*Tb)/8*(sinc(f*Tb/2)).^2.*(1-cos(2*pi*f*Tb));
        figure(2);
        g=plot(f,P);
        title('PSD BIPOLAR-RZ');
        hold on;xlabel('Frecuencia (Hz)');ylabel('Densidad espectral de potencia (W/Hz)');
        axis([0 L 0 1.05*Tb]);set(g,'LineWidth',2.5);
        set(gca,'XTickMode','manual','XTick',[0,0.5*R,R,1.5*R,2*R]);grid on;
        set(gca,'YTickMode','manual','YTick',[0,0.25*Tb,0.5*Tb,0.75*Tb,Tb]);
        set(gca,'XTickLabel',{['0'];['0.5R'];['R'];['1.5R'];['2R']});
        set(gca,'YTickLabel',{['0'];['0.25Tb'];['0.5Tb'];['0.75Tb'];['Tb']});
               
     case 'manchester_nrz'
        hold off;
        A=1;
        Tb=1.0e-3;
        R=1/Tb;
        L=2*R;
        f=0:L/100:L;
        P=(A.^2*Tb)*(sinc(f*Tb/2)).^2.*(sin(pi*f*Tb/2)).^2;
        figure(2);
        g=plot(f,P);
        title('PSD MANCHESTER-NRZ');
        hold on;xlabel('Frecuencia (Hz)');ylabel('Densidad espectral de potencia (W/Hz)');
        axis([0 L 0 1.05*Tb]);set(g,'LineWidth',2.5);
        set(gca,'XTickMode','manual','XTick',[0,0.5*R,R,1.5*R,2*R]);grid on;
        set(gca,'YTickMode','manual','YTick',[0,0.25*Tb,0.5*Tb,0.75*Tb,Tb]);
        set(gca,'XTickLabel',{['0'];['0.5R'];['R'];['1.5R'];['2R']});
        set(gca,'YTickLabel',{['0'];['0.25Tb'];['0.5Tb'];['0.75Tb'];['Tb']});
end

end

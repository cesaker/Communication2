sigma=[0.05:0.05:1.0];
for i=1:length(sigma)
[Peb(i),Pes(i),Pebt(i),Pest(i),EbNodB(i),EsNodB(i),sr]=QPSK(10e3,sigma(i),0,0);
end

hold on
figure(1)
plot(EbNodB,Peb)
plot(EbNodB,Pebt,'r')
title('Peb');
hold off
hold on
figure(2)
plot(EbNodB,Pes)
plot(EbNodB,Pest,'r')
title('Pes');
hold off
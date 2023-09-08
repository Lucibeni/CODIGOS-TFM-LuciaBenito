nT = [303:50:953];
azules = [0 0 0; 0 0 0.15; 0 0 0.3; 0 0 0.6; 0 0 0.8; 0 0 1; 0 0.2 1; 0 0.5 1; 0 0.8 1; 0 1 1; 0 1 0.8; ;0 1 0.5 ;0 1 0.2; 0 1 0];

for ii=1:length(nT)
    ii
    evolucion_tumorales=distintasT(nT(ii));

    tumorales_finales_sigue(ii) = evolucion_tumorales(end); %Guarda las células tumorales que hay en el último paso de tiempo
    tumorales_finales10(ii) = evolucion_tumorales(432); %Guarda las células tumorales que hay en el paso de tiempo en el que han muerto todas las células T 

    minimos10(ii) = min(evolucion_tumorales); %Calcula cual es el número de células tumorales más pequeño que se alcanza
    times = find(evolucion_tumorales==min(evolucion_tumorales)); 
    tiempo_min(ii) = times(1); %Calcula el tiempo para el cual se alcanza ese mínimo

    if tumorales_finales_sigue(ii) >= 903 %Con este 'if' lo que hacemos es ver para que paso de tiempo el tumor vuelve a tener su tamaño inicial
        tiempos = find(evolucion_tumorales >= 903);
        k=1;
        while tiempos(k) < tiempo_min(ii)
            k = k+1;            
        end
        reset10(ii) = tiempos(k);
    else
        reset10(ii) = 0;
    end

    %FIGURA: Gráfica de la evolución de las células tumorales a lo largo
    %del tiempo para distintas dosis de células T
% % % %     figure(2)
% % % %     P1(ii)=plot(1:435+520, [evolucion_tumorales], 'Color', azules(ii,1:3),'LineWidth',1.3, 'DisplayName',['T_0 = ' num2str(nT(ii))]);
% % % %     hold on
% % % %     plot(1:reset(ii), (903)*ones(1,reset(ii)), ':r','LineWidth',1.2)
% % % %     plot(reset(ii)*ones(1,903-500+1),500:903, ':r','LineWidth',1.2)
end


sdfghjkl

xline(432)
xlabel('Tiempo en días') 
ylabel('Número de células tumorales') 
title('Evolución del número de células tumorales a lo largo de 35 días')
legend(P1)
set(gca, 'xtick',[0,120,240,360,480,600,720,840,960], 'xticklabels', {'0','5', '10', '15', '20', '25', '30', '35', '40'}) %Cambio el eje a días en vez de en horas

%legend('nT_0 = 653','nT_0 = 678','nT_0 = 703','nT_0 = 728','nT_0 = 753','nT_0 = 778','nT_0 = 803','nT_0 = 828','nT_0 = 853','nT_0 = 878','nT_0 = 903','nT_0 = 928','nT_0 = 953')


%FIGURA: Número de células tumorales finales en función del número inicial de células T

%tumorales_finales_media = (tumorales_finales1 + tumorales_finales2 + tumorales_finales3 + tumorales_finales4 + tumorales_finales5 + tumorales_finales6 + tumorales_finales7 + tumorales_finales8 + tumorales_finales9 + tumorales_finales10)/10;

figure(3)
plot(nT,tumorales_finales,'Marker','*', "Color", azules(7,1:3), 'LineWidth',2.8)
hold on
xlabel('Número inicial de células T') 
ylabel('Número de células tumorales finales') 
title('Número de células tumorales finales en función del número inicial de células T')
yline(903, ':b','LineWidth',2)
legend('Número de células tumorales finales', 'Número inicial de células tumorales')
polinomio = polyfit(nT,tumorales_finales,2);
plot(nT, polinomio(1)*nT.^2+polinomio(2)*nT+polinomio(3),':r', 'LineWidth',1.8)


%FIGURA: Número mínimo de células tumorales alcanzado en función del número inicial de células T
figure(4)
plot(nT,minimos,'Marker','*', "Color", azules(7,1:3), 'LineWidth',2.8)
xlabel('Número inicial de células T') 
ylabel('Número mínimo de células tumorales que se alcanza durante el tratamiento') 
title('Número mínimo de células tumorales alcanzado en función del número inicial de células T')
hold on
polinomio_min = polyfit(nT,minimos,2);
plot(nT, polinomio_min(1)*nT.^2+polinomio_min(2)*nT+polinomio_min(3),':r', 'LineWidth',1.8)


%FIGURA: Número de horas transcurridas cuando se alcanza el mínimo de células tumorales en función del número inicial de células T
figure(5)
plot(nT,tiempo_min,'Marker','*', "Color", azules(7,1:3), 'LineWidth',2.8)
xlabel('Número inicial de células T') 
ylabel('Número de horas transcurridas cuando se alcanza mínimo de células tumorales') 
title('Número de horas transcurridas cuando se alcanza el mínimo de células tumorales en función del número inicial de células T')

% figure(6)
% plot(nT,reset,'Marker','*', "Color", azules(12,1:3), 'LineWidth',2.8)
% xlabel('Número inicial de células T') 
% ylabel('Número de horas transcurridas cuando se alcanza el número inicial de células tumorales') 
% title('Número de horas transcurridas cuando se alcanza el número inicial de células tumorales en función del número inicial de células T')

hold off
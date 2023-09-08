aporte = [300:50:1600];
azules = [0 0 0; 0 0 0.2; 0 0 0.3; 0 0 0.5; 0 0 0.6; 0 0 0.7; 0 0 0.8; 0 0 0.9; 0 0.1 1; 0 0.2 1; 0 0.3 1; 0 0.4 1; 0 0.5 1; 0 0.6 1; 0 0.7 1; 0 0.8 1; 0 0.9 1; 0 1 1;0 1 0.9;0 1 0.8;0 1 0.7;0 1 0.6;0 1 0.5;0 1 0.4;0 1 0.3;0 1 0.2;0 1 0.1;0 1 0];

for ii=1:length(aporte)
    ii
    evolucion_tumorales=aporteT(aporte(ii)); 

    tumorales_finales(ii) = evolucion_tumorales(end);  %Guarda las células tumorales que hay en el último paso de tiempo (pasado unos 3 meses)

    minimos(ii) = min(evolucion_tumorales); %Calcula cual es el número de células tumorales más pequeño que se alcanza
    times = find(evolucion_tumorales==min(evolucion_tumorales));
    tiempo_min(ii) = times(1); %Calcula el tiempo para el cual se alcanza ese mínimo

    if tumorales_finales(ii) >= 903 %Con este 'if' lo que hacemos es ver para que paso de tiempo el tumor vuelve a tener su tamaño inicial
        tiempos = find(evolucion_tumorales >= 903);
        k=1;
        while tiempos(k) < tiempo_min(ii)
            k = k+1;            
        end
        reset(ii) = tiempos(k);
    else
        reset(ii) = 0;
    end

    %FIGURA: Gráfica de la evolución de las células tumorales a lo largo
    %del tiempo para distintas dosis de células T
    figure(2)
    plot(1:2200, [evolucion_tumorales], 'Color', azules(ii,1:3),'LineWidth',1.3)
    hold on
end
xlabel('Tiempo en horas') 
ylabel('Número de células tumorales') 
title('Evolución del número de células tumorales a lo largo de meses')
legend('nT_0 = 700','nT_0 = 725','nT_0 = 750','nT_0 = 775','nT_0 = 800','nT_0 = 825','nT_0 = 850','nT_0 = 875','nT_0 = 900','nT_0 = 925','nT_0 = 950','nT_0 = 975','nT_0 = 1000','nT_0 = 1025','nT_0 = 1050','nT_0 = 1075','nT_0 = 1100','nT_0 = 1125','nT_0 = 1150','nT_0 = 1175','nT_0 = 1200','nT_0 = 1225','nT_0 = 1250','nT_0 = 1275','nT_0 = 1300','nT_0 = 1325','nT_0 = 1350','nT_0 = 1375','nT_0 = 1400','nT_0 = 1425','nT_0 = 1450')
set(gca, 'xtick',[0,120,240,360,480,600,720,840,960,1080,1200,1320,1440,1560,1680,1800,1920,2040,2160,2280], 'xticklabels', {'0','5', '10', '15', '20', '25', '30', '35', '40', '45','50', '55','60','65','70','75','80','85','90','95'}) %Cambio el eje a días en vez de en horas

%FIGURA: Número de células tumorales finales en función del número de
%células T que se van aportando
figure(3)
plot(aporte,tumorales_finales,'Marker','*', "Color", azules(13,1:3), 'LineWidth',2.8)
hold on
xlabel('Aporte de células T') 
ylabel('Número de células tumorales finales') 
title('Número de células tumorales finales en función del número de células T')
yline(903, ':b','LineWidth',2)
legend('Número de células tumorales finales', 'Número inicial de células tumorales')
% polinomio = polyfit(aporte,tumorales_finales,2);
% plot(aporte, polinomio(1)*aporte.^2+polinomio(2)*aporte+polinomio(3),':r', 'LineWidth',1.8)
% ajuste = fit(aporte',tumorales_finales','exp1');
% plot(aporte, ajuste.a*exp(ajuste.b*aporte),':m', 'LineWidth',1.8)

ajuste2 = fit(aporte',tumorales_finales','exp2');
plot(aporte, ajuste2.a*exp(ajuste2.b*aporte) + ajuste2.c*exp(ajuste2.d*aporte),':r', 'LineWidth',1.8)


%FIGURA: Número mínimo de células tumorales alcanzado en función del aporte 
%de células T
figure(4)
plot(aporte,minimos,'Marker','*', "Color", azules(13,1:3), 'LineWidth',2.8)
xlabel('Aporte de células T') 
ylabel('Número mínimo de células tumorales que se alcanza durante el tratamiento') 
title('Número mínimo de células tumorales alcanzado en función del número de células T')


%FIGURA: Número de horas transcurridas cuando se alcanza el mínimo de células
%tumorales en función del aporte de células T
figure(5)
plot(aporte,tiempo_min,'Marker','*', "Color", azules(13,1:3), 'LineWidth',2.8)
xlabel('Aporte de células T') 
ylabel('Número de horas transcurridas cuando se alcanza mínimo de células tumorales') 
title('Número de horas transcurridas cuando se alcanza el mínimo de células tumorales en función del número de células T')

% figure(6)
% plot(nT,reset,'Marker','*', "Color", azules(12,1:3), 'LineWidth',2.8)
% xlabel('Número inicial de células T') 
% ylabel('Número de horas transcurridas cuando se alcanza el número inicial de células tumorales') 
% title('Número de horas transcurridas cuando se alcanza el número inicial de células tumorales en función del número inicial de células T')

hold off
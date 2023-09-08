aporte = [400:5:500];
azules = [0 0 0; 0 0 0.4;  0 0 0.8;  0 0 1;  0 0.2 1; 0 0.4 1;  0 0.6 1;  0 0.8 1; 0 1 1;0 1 0.8;0 1 0.6;0 1 0.4;0 1 0.2;0 1 0];

for ii=1:length(aporte)
    ii
    evolucion_tumorales=aporteT(aporte(ii));

    %tumorales_finales(ii) = evolucion_tumorales(end); %Guarda las células tumorales que hay en el último paso de tiempo (pasado unos 6 meses)
    pendiente(ii) = (evolucion_tumorales(end)-evolucion_tumorales(2200))/4400; %Calcula la pendiente de la evolución de las células tumorales desde el día 90 al 180, es decir, aprox los últimos 3 meses


% % %     figure(2)
% % %     plot(1:2200, [evolucion_tumorales], 'Color', azules(ii,1:3),'LineWidth',1.3)
% % %     hold on
% % %     plot(1:2200, 903+pendiente(ii)*[1:2200],':', 'Color', azules(ii,1:3),'LineWidth',1)
end
% % % xlabel('Tiempo en horas') 
% % % ylabel('Número de células tumorales') 
% % % title('Evolución del número de células tumorales a lo largo de meses')
% % % %legend('nT_0 = 700','nT_0 = 725','nT_0 = 750','nT_0 = 775','nT_0 = 800','nT_0 = 825','nT_0 = 850','nT_0 = 875','nT_0 = 900','nT_0 = 925','nT_0 = 950','nT_0 = 975','nT_0 = 1000','nT_0 = 1025','nT_0 = 1050','nT_0 = 1075','nT_0 = 1100','nT_0 = 1125','nT_0 = 1150','nT_0 = 1175','nT_0 = 1200','nT_0 = 1225','nT_0 = 1250','nT_0 = 1275','nT_0 = 1300','nT_0 = 1325','nT_0 = 1350','nT_0 = 1375','nT_0 = 1400','nT_0 = 1425','nT_0 = 1450')

%FIGURA: En esta gráfica se representa la pendiente de la evolución de las
%células tumorales en los 3 últimos meses de simulación (del dia 90 al 180)
%en función del aporte de células T
figure(3)
plot(aporte, pendiente, 'Marker','*', 'Color',azules(5,1:3),'LineWidth',2.8)
xlabel('Aporte de células T') 
ylabel('Pendiente de crecimiento') 
title('Dependencia de la pendiente y el aporte de células T')
hold on
yline(0)


hold off
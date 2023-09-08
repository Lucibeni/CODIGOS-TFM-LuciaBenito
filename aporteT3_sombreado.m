semillas = [10, 100, 500, 888, 1002, 1500, 1999, 2000, 2067, 2110];

%%%%%%%%%% 300 rojo-grante
for jj=1:length(semillas)
        evolucion_tumorales_300(jj,:)=aporteT3(300,semillas(jj)); 
end
    figure(2)
    media_evolucion_tumorales_300 = sum(evolucion_tumorales_300)./jj;
    maximo_300 = max(evolucion_tumorales_300);
    minimo_300 = min(evolucion_tumorales_300);
    fill([1:4400, fliplr(1:4400)], [maximo_300, fliplr(minimo_300)], [255 148 162]./255, 'FaceAlpha',.3)
    hold on
    plot(1:4400, media_evolucion_tumorales_300, 'Color', [0.6350 0.0780 0.1840],'LineWidth',4)
    plot(1:4400, maximo_300, 'Color', [255 148 162]./255,'LineWidth',4)
    plot(1:4400, minimo_300, 'Color', [255 148 162]./255,'LineWidth',4)
  
%%%%%%%%%% 525 azules
for jj=1:length(semillas)
        evolucion_tumorales_525(jj,:)=aporteT3(525,semillas(jj)); 
end
figure(1)
    media_evolucion_tumorales_525 = sum(evolucion_tumorales_525)./jj;
    maximo_525 = max(evolucion_tumorales_525);
    minimo_525 = min(evolucion_tumorales_525);
    fill([1:4400, fliplr(1:4400)], [maximo_525, fliplr(minimo_525)], [0.3010 0.7450 0.9330], 'FaceAlpha',.3)
    plot(1:4400, media_evolucion_tumorales_500, 'Color', 'b','LineWidth',4)
    plot(1:4400, maximo_525, 'Color', [0.3010 0.7450 0.9330],'LineWidth',4)
    plot(1:4400, minimo_525, 'Color', [0.3010 0.7450 0.9330],'LineWidth',4)



%%%%%%%%%% 1600 verdes
for jj=1:length(semillas)
        evolucion_tumorales_1600(jj,:)=aporteT3(1600,semillas(jj));
end
figure(2)
    media_evolucion_tumorales_1600 = sum(evolucion_tumorales_1600)./jj;
    maximo_1600 = max(evolucion_tumorales_1600);
    minimo_1600 = min(evolucion_tumorales_1600);
    fill([1:4400, fliplr(1:4400)], [maximo_1600, fliplr(minimo_1600)], [0.4660 0.6740 0.1880], 'FaceAlpha',.3)
    plot(1:4400, media_evolucion_tumorales_1600, 'Color', 'g','LineWidth',4)
    plot(1:4400, maximo_1600, 'Color', [0.4660 0.6740 0.1880],'LineWidth',4)
    plot(1:4400, minimo_1600, 'Color', [0.4660 0.6740 0.1880],'LineWidth',4)


xlabel('Tiempo en días') 
ylabel('Número de células tumorales') 
title('Evolución del número de células tumorales a lo largo de 6 meses')
set(gca, 'xtick',[0,360,720,1080,1440, 1800, 2160, 2520, 2880, 3240, 3600, 3960, 4320], 'xticklabels', {'0','15', '30', '45', '60', '75', '90', '105', '120', '135', '150', '165', '180'})

hold off
n = 100; % nxn tamaño del autómata (en ocasiones le haremos bordes a las matrices para evitar coger posiciones fuera de rango)
ntum = 900; %nº de celulas infectadas
A = zeros(n+1,n+1); %Matriz que guarda las posiciones de las cosas y genera la representación gráfica
B = zeros(n+1,n+1); %Matriz que guarda la carga citotóxica que tiene una célula T para disparar
D = zeros(n+1,n+1); %Matriz auxiliar para poder guardar donde estan las células tumorales que crecen por debajo del vaso sanguíneo pero poder pintar el vaso sanguíneo en esa posición
MT = zeros(n+1,n+1); %Matriz que guarda el tiempo de vida que le queda a una célula T para desaparecer del autómata
k=0; %contador de celulas infectadas

rng('default');
rng(2110); %el numerito es la semilla


%%%%%% Colocar vasos sanguíneos %%%%%%%%

% Parto la forma en la que dibujo los
% vasos sanguineos para poder meter en una lista las posibles posiciones de
% las celulas T en los bordes de los vasos sanguíneos
% En los índices se encuentran ciertas operaciones que me ayudaron a
% controlar la situación que tenían las distintas ramificaciones de los
% vasos sanguíneos
A(13,36:100) = 3*ones(1,65);
A(15,36:100) = 3*ones(1,65);
for i = 15:20
    A(i-8:i-7,i+24-10:i+1+24-10) = 3*ones(2,2);
end
for i = 9:20
    A(-i+35:-i+36,i+24-10:i+1+24-10) = 3*ones(2,2);
end

A(36,30:100) = 3*ones(1,71);
A(38,30:100) = 3*ones(1,71);

A(63,47:100) = 3*ones(1,54);
A(65,47:100) = 3*ones(1,54);
for i = 15:26
    A(i-8+44:i-7+44,i+24-5:i+1+24-5) = 3*ones(2,2);
end
for i = 8:26
    A(-i+35+56:-i+36+56,i+24-5:i+1+24-5) = 3*ones(2,2);
end
for i = 8:15
    A(i-8+68:i-7+68,i+24-5:i+1+24-5) = 3*ones(2,2);
end

A(90,33:100) = 3*ones(1,68);
A(92,33:100) = 3*ones(1,68);
for i = 18:21
    A(-i+35+77:-i+36+77,i+14-3:i+1+14-3) = 3*ones(2,2);
end

[fila_posT, col_posT] = find(A==3); %fila_postT y col_posT guardan todas las posiciones(fila,columna) de bordes de vaso sanguíneo

%Relleno las posiciones del centro del vaso sanguíneo
A(14,36:100) = 3*ones(1,65);
A(37,30:100) = 3*ones(1,71);
A(64,47:100) = 3*ones(1,54);
A(91,33:100) = 3*ones(1,68);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%% TUMOR LATERAL %%%%%%%%%%%%%%%%%%%

%En las matrices A y D se van colocando las células tumorales. Se va
%creando a la izquierda un tumor compacto, es decir, una célula contigua a
%otra.
%Con fila1 y col1 se guardan todas las posiciones(fila,columna) de la celda
%superior izquierda de cada célula tumoral.

A(1:end-1, 1:8+16) = 2;
D(1:end-1, 1:8+16) = 2;
col1(1:n/2) = ones(1,n/2);
fila1(1:n/2) = 2.*[1:n/2]-1;
col1(n/2+1:n) = 3*ones(1,n/2);
fila1(n/2+1:n) = 2.*[1:n/2]-1;
col1(n+1:n+n/2) = 5*ones(1,n/2);
fila1(n+1:n+n/2) = 2.*[1:n/2]-1;
col1(n+(n/2)+1:2*n) = 7*ones(1,n/2);
fila1(n+(n/2)+1:2*n) = 2.*[1:n/2]-1;

col1(2*n+1:2*n+n/2) = 9*ones(1,n/2);
fila1(2*n+1:2*n+n/2) = 2.*[1:n/2]-1;
col1(2*n+n/2+1:3*n) = 11*ones(1,n/2);
fila1(2*n+n/2+1:3*n) = 2.*[1:n/2]-1;
col1(3*n+1:3*n+n/2) = 13*ones(1,n/2);
fila1(3*n+1:3*n+n/2) = 2.*[1:n/2]-1;
col1(3*n+(n/2)+1:4*n) = 15*ones(1,n/2);
fila1(3*n+(n/2)+1:4*n) = 2.*[1:n/2]-1;

col1(4*n+1:4*n+n/2) = 17*ones(1,n/2);
fila1(4*n+1:4*n+n/2) = 2.*[1:n/2]-1;
col1(4*n+n/2+1:5*n) = 19*ones(1,n/2);
fila1(4*n+n/2+1:5*n) = 2.*[1:n/2]-1;
col1(5*n+1:5*n+n/2) = 21*ones(1,n/2);
fila1(5*n+1:5*n+n/2) = 2.*[1:n/2]-1;
col1(5*n+(n/2)+1:6*n) = 23*ones(1,n/2);
fila1(5*n+(n/2)+1:6*n) = 2.*[1:n/2]-1;

%Tras ello, se va creando una zona del tumor más infiltrativa, es decir,
%las células no son contiguas, sino que van proliferando aleatoriamente
%hasta insertar un número total de ntum=900 células tumorales

while k < ntum
    [A(1:n,1:n),B(1:n,1:n),D(1:n,1:n),MT(1:n,1:n),fila1,col1] = proliferar_tfm(fila1,col1,A(1:n,1:n),B(1:n,1:n),D(1:n,1:n),MT(1:n,1:n),0.7);
    k = length(col1); %Contar nº de células que llevamos insertadas
end


%%%%%%%%% Rellenar con células T %%%%%%%%%%%%%%%%%%
nT = 300; %Número de células T que se van a introducir inicialmente en la simulación
kT=0; %contador de celulas T
while kT < nT
    numerito = ceil(length(col_posT)*rand(1,1)); %Generamos una posicion de la lista de posiciones de los vasos sanguíneos para introducir en ella una célula T
    if A(fila_posT(numerito,1),col_posT(numerito,1)) == 3 % si el valor de la posicion es 3, añadimos la celula T (pues si en dicha posición ya hay una célula T no puede insertarse otra)
        A(fila_posT(numerito,1), col_posT(numerito,1)) = 1.25; %Pongo 1.25 para diferenciar las células T que están dentro de los vasos sanguíneos
        %La células T no tienen nada en B hasta que no salen del vaso sanguíneo, ahí es donde tomaran valor en B 
        vida = round(12 + 6*rand(1,1)); %Generamos aleatoriamente los días (entre 12 y 18) de vida que tendrá la célula T que estamos introduciendo. 
        MT(fila_posT(numerito,1), col_posT(numerito,1)) = 24*vida; %Guardamos este valor en horas (por ello *24) en la matriz MT para que este valor siga los movimientos de esa célula T
    end  
    kT = numel(A(A==1.25)); %Contar nº de células que llevamos insertadas
end

%Amplio las matrices A y B una casilla por arriba y por la izquierda
%Ya las teniamos ampliadas por abajo y por la derecha
A = [zeros(1,n+1); A];
A = [zeros(n+2,1) A];
B = [zeros(1,n+1); B];
B = [zeros(n+2,1) B];
D = [zeros(1,n+1); D];
D = [zeros(n+2,1) D];
MT = [zeros(1,n+1); MT];
MT = [zeros(n+2,1) MT];



tasa_intro = 1600/360; %Tasa de células T que entran por el vaso sanguíneo cada hora para que se mantenga un número constante de células T
entero = floor(tasa_intro); %Número entero de células T que entran cada hora
decimales = tasa_intro - entero; %Lo que se 'acumula' cada hora para saber cada cuántas horas se debe meter una célula T extra para corresponder con la tasa calculada
posiciones = [13, 14, 15, 36, 37, 38, 63, 64, 65, 90, 91, 92]; %Filas en las que puede aparecer la célula T dentro del vaso sanguíneo por la izquierda
acumulador = 0;


%%%% Graficamos los movimientos de las células T:
        % salir del vaso sanguineo,
        % disparar, 
        % matar, 
        % mover,
        % morir +o- 15 días (entre 12 y 18 días)
%%%% Y la proliferación de células tumorales
T =4400; % tiempo en horas de la simulación (equivale a unos 3 meses)

figure(1)

for tiempo = 1:T
      carga_toxica_ctumoral = zeros(1,length(fila1)); %Vector que guarda en cada posición la carga citotóxica que se le ha disparado a cada célula tumoral.(Se crea de nuevo en cada paso de tiempo)

      COLOR=[178/255 218/255 250/255; 147/255 112/255 219/255; 0.6350 0.0780 0.1840; 1 0 0]; %Código de colores de la representación gráfica del autómata
      clim([0 3]) %Para la representación gráfica del autómata
      colormap(COLOR) %Para la representación gráfica del autómata
      imagesc(A(2:end-1,2:end-1)) %Para la representación gráfica del autómata
      title(['La simulación va por la hora ' num2str(tiempo) ' de  ' num2str(T+200)]) %Para la representación gráfica del autómata

      [A,B,MT] = salir_tfm(A(2:end-1,2:end-1),B(2:end-1,2:end-1),MT(2:end-1,2:end-1)); %Función mediante la cuál se llevan a cabo los cambios relativos a que las células T salgan del vaso sanguíneo. %pierdo los bordes
      [A,MT] = moverDentro_tfm(A,MT);%Función mediante la cuál se llevan a cabo los cambios relativos al avance de las células T a través de los vasos sanguíneos
            
      %Introducimos celulas T en el lado derecho de los vasos sanguineos
      acumulador = acumulador + decimales; %Va acumulando la parte decimal de la tasa de entrada. Cuando esta llega a 1, se debe introcudir una célula T para corresponder a la tasa calculada
      if acumulador > 1
          acumulador = acumulador - 1; %Renovamos el acumulador para que sea menor que 1 pero sin perder lo que nos hemos pasado de 1
          num_posi = round(1+11*rand(1,1)); %Generamos una de las posiciones de la lista donde se guardan las posibles posiciones por donde pueden llegar las células T en los vasos sanguíneos
          while A(posiciones(num_posi), 100) == 1.25 %Si en la posición generada ya hay una célula T, se genera otra posición hasta que esta esté libre
              num_posi = round(1+11*rand(1,1));
          end
          A(posiciones(num_posi), 100) = 1.25; %Introducimos la célula T. Vale 1.25 para distinguir que está dentro del vaso
          vida = round(12 + 6*rand(1,1)); %Generamos aleatoriamente los días (entre 12 y 18) de vida que tendrá la célula T que estamos introduciendo. 
          MT(posiciones(num_posi), 100) = 24*vida; %Guardamos este valor en horas (por ello *24) en la matriz MT para que este valor siga los movimientos de esa célula T
      end
      for iintro = 1:entero %Ahora, del mismo modo que antes, introducimos tantas células T como indica la tasa en cada paso de tiempo
          num_posi = round(1+11*rand(1,1));
          while A(posiciones(num_posi), 100) == 1.25
              num_posi = round(1+11*rand(1,1));
          end
          A(posiciones(num_posi), 100) = 1.25;
          vida = round(12 + 6*rand(1,1));
          MT(posiciones(num_posi), 100) = 24*vida;
      end


      %Vuelvo a añadir los bordes:
      A = [zeros(1,n); A; zeros(1,n)];
      A = [zeros(n+2,1) A zeros(n+2,1)];
      B = [zeros(1,n); B; zeros(1,n)];
      B = [zeros(n+2,1) B zeros(n+2,1)];
      MT = [zeros(1,n); MT; zeros(1,n)];
      MT = [zeros(n+2,1) MT zeros(n+2,1)];

      [carga_toxica_ctumoral,B] = disparar_tfm(fila1, col1, carga_toxica_ctumoral, B);  %Función mediante la cuál se contabiliza la carga citotóxica que una célula tumoral recibe en total en un paso de tiempo 
      %%% Ahora, en función de esta carga citotóxica que ha recibido cada
      %%% célula tumoral, van a morirse con mayor o menor probabilidad

      %Lo primero es inicializar nuevos vectores que guarden la fila, la columna y 
      %la carga citotóxica de cada célula tumoral, para ir guardando la
      %información de las que han sobrevivido al ataque de las células T
      fila_new = [];
      col_new = [];
      carga_toxica_ctumoral_new = [];
      for i = 1:length(carga_toxica_ctumoral)
        matamos = morir_tfm(carga_toxica_ctumoral(i)); %Para cada célula se decide si muere o no en función de la carga citotóxica recibida 
        if matamos == 1 %Si la célula tumoral muere:
            D(fila1(i)+1, col1(i)+1) = 0; %Quitamos las 4 posiciones de la célula tumoral en la matriz auxiliar D
            D(fila1(i)+2, col1(i)+1) = 0;
            D(fila1(i)+1, col1(i)+2) = 0;
            D(fila1(i)+2, col1(i)+2) = 0;

            if A(fila1(i)+1, col1(i)+1) == 2 %Quitamos de la matriz A las posiciones de la célula tumoral que estaban representadas gráficamente
                A(fila1(i)+1, col1(i)+1) = 0;
            end
            if A(fila1(i)+2, col1(i)+1) == 2
                A(fila1(i)+2, col1(i)+1) = 0;
            end
            if A(fila1(i)+1, col1(i)+2) == 2
                A(fila1(i)+1, col1(i)+2) = 0;
            end
            if A(fila1(i)+2, col1(i)+2) == 2
                A(fila1(i)+2, col1(i)+2) = 0;
            end
        else %Si la célula tumoral NO muere, se guarda su fila, su columna y su carga citotóxica recibida en los nuevos vectores
            fila_new = [fila_new fila1(i)];            
            col_new = [col_new col1(i)];
            carga_toxica_ctumoral_new = [carga_toxica_ctumoral_new carga_toxica_ctumoral(i)];
        end
      end     
      %Ahora cambiamos los vectores que estabamos considerando por los nuevos vectores que se han ido creando
      fila1 = fila_new;
      col1 = col_new;
      carga_toxica_ctumoral = carga_toxica_ctumoral_new;

      %Con 'mover' pierdo los bordes:
      [A,B,MT] = mover_tfm(A(2:end-1,2:end-1),B(2:end-1,2:end-1),MT(2:end-1,2:end-1)); %Esta función lleva a cabo los cambios relativos a los desplazamientos de las células T fuera de los vasos sanguíneos y a la recarga de su carga citotóxica
      [A,B,D,MT,fila1,col1] = proliferar_tfm(fila1,col1,A,B,D(2:end-1,2:end-1),MT,0.01); %Esta función lleva a cabo los cambios relativos a la proliferación de las células tumorales
      
      [fila_adios,col_adios] = find(MT==1); %fila_adios y col_adios guardan las posiciones de aquellas células T las cuales estan en su última hora de vida
      for dis = 1:length(col_adios) %Cada célula T que está en su última hora de vida se elimina del autómata
          if A(fila_adios(dis), col_adios(dis)) == 1 %Si está fuera del vaso sanguíneo:
            A(fila_adios(dis), col_adios(dis)) = 0; %La celda pasa a ser de espacio intercelular
          else
            A(fila_adios(dis), col_adios(dis)) = 3; %Sino, la celda se convierte en vaso sanguíneo
          end
          B(fila_adios(dis), col_adios(dis)) = 0; %Su carga citotóxica disponible para disparar desaparece
      end

      MT = MT - (MT>0); %En cada paso temporal se le resta una hora de vida 

      %Vuelvo a añadir los bordes:
      A = [zeros(1,n); A; zeros(1,n)];
      A = [zeros(n+2,1) A zeros(n+2,1)];
      B = [zeros(1,n); B; zeros(1,n)];
      B = [zeros(n+2,1) B zeros(n+2,1)];
      D = [zeros(1,n); D; zeros(1,n)];
      D = [zeros(n+2,1) D zeros(n+2,1)];
      MT = [zeros(1,n); MT; zeros(1,n)];
      MT = [zeros(n+2,1) MT zeros(n+2,1)];

      drawnow %Para la representación gráfica del autómata
      %pause(0.05)
      evolucion(tiempo) = length(fila1); %Hace el recuento de células tumorales que hay en cada paso de tiempo y lo guarda en un vector
      evolucionT(tiempo) = numel(A(A==1))+numel(A(A==1.25)); %Hace el recuento de células T que hay en cada paso de tiempo y lo guarda en un vector
end


      %FIGURAS DE PRUEBA: Representaciones de la evolución del número de células
      %tumorales y células T totales, respectivamente, a lo largo del tiempo
      figure(2)
      plot(1:T, [evolucion],'m')
      figure(3)
      plot(1:T, [evolucionT],'c')
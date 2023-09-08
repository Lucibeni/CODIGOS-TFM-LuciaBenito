%% CóDIGO tfm_introducirT HECHO FUNCION

function evolucion = aporteT3(flujoT,semilla)
n = 100; % nxn tamaño de la matriz
ntum = 900; %nº de celulas infectadas
A = zeros(n+1,n+1);
B = zeros(n+1,n+1);
D = zeros(n+1,n+1);
MT = zeros(n+1,n+1);
k=0; %contador de celulas infectadas

rng('default');
rng(semilla); %el numerito es la semilla


%%%%%% Colocar vasos sanguíneos %%%%%%%%

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

[fila_posT, col_posT] = find(A==3);

A(14,36:100) = 3*ones(1,65);
A(37,30:100) = 3*ones(1,71);
A(64,47:100) = 3*ones(1,54);
A(91,33:100) = 3*ones(1,68);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%% TUMOR LATERAL %%%%%%%%%%%%%%%%%%%

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

while k < ntum
    [A(1:n,1:n),B(1:n,1:n),D(1:n,1:n),MT(1:n,1:n),fila1,col1] = proliferar_tfm(fila1,col1,A(1:n,1:n),B(1:n,1:n),D(1:n,1:n),MT(1:n,1:n),0.7);
    k = length(col1); %Contar nº de células que llevamos insertadas
end


%%%%%%%%% Rellenar con células T %%%%%%%%%%%%%%%%%%
nT = 300;
kT=0; %contador de celulas T
while kT < nT
    numerito = ceil(length(col_posT)*rand(1,1)); %Generamos un posicion de la lista de posiciones de las celulas T dentro de los vasos sanguíneos
    if A(fila_posT(numerito,1),col_posT(numerito,1)) == 3 % si el valor de la posicion es 3, añadimos la celula T
        A(fila_posT(numerito,1), col_posT(numerito,1)) = 1.25; %Pongo 1.25 para diferenciar las células T que están dentro de los vasos sanguíneos
        %La células T no tienen nada en B hasta que no salen del vaso sanguíneo, ahí es donde tomaran valor en B 
        %B(fila_posT(numerito,1),col_posT(numerito,1)) = 9;
        vida = round(12 + 6*rand(1,1));
        MT(fila_posT(numerito,1), col_posT(numerito,1)) = 24*vida;
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

%%%% Graficamos los movimientos de las células T:
% salir del vaso sanguineo,
% disparar, 
% matar, 
% mover,
% morir +o- 15 días (entre 12 y 18 días)
%%%% Y la proliferación de células tumorales
T =4400; % tiempo en horas

tasa_intro = flujoT/360;
entero = floor(tasa_intro);
decimales = tasa_intro - entero;
posiciones = [13, 14, 15, 36, 37, 38, 63, 63, 65, 90, 91, 92];
acumulador = 0;

%figure(1)

for tiempo = 1:T
      carga_toxica_ctumoral = zeros(1,length(fila1));
      COLOR=[178/255 218/255 250/255; 147/255 112/255 219/255; 0.6350 0.0780 0.1840; 1 0 0];
      clim([0 3])
      colormap(COLOR)
      %imagesc(A(2:end-1,2:end-1))
      %title(['hora ' num2str(tiempo) ' de  ' num2str(T)])
      [A,B,MT] = salir_tfm(A(2:end-1,2:end-1),B(2:end-1,2:end-1),MT(2:end-1,2:end-1)); %pierdo doble borde
      [A,MT] = moverDentro_tfm(A,MT);
      
      %Introducimos celulas T en el lado derecho de los vasos sanguineo
      acumulador = acumulador + decimales;
      if acumulador > 1
          acumulador = acumulador - 1;
          num_posi = round(1+11*rand(1,1));
          while A(posiciones(num_posi), 100) == 1.25
              num_posi = round(1+11*rand(1,1));
          end
          A(posiciones(num_posi), 100) = 1.25;
          vida = round(12 + 6*rand(1,1));
          MT(posiciones(num_posi), 100) = 24*vida;
      end
      for iintro = 1:entero
          num_posi = round(1+11*rand(1,1));
          while A(posiciones(num_posi), 100) == 1.25
              num_posi = round(1+11*rand(1,1));
          end
          A(posiciones(num_posi), 100) = 1.25;
          vida = round(12 + 6*rand(1,1));
          MT(posiciones(num_posi), 100) = 24*vida;
      end


      %Vuelvo a añadir el doble borde:
      A = [zeros(1,n); A; zeros(1,n)];
      A = [zeros(n+2,1) A zeros(n+2,1)];
      B = [zeros(1,n); B; zeros(1,n)];
      B = [zeros(n+2,1) B zeros(n+2,1)];
      MT = [zeros(1,n); MT; zeros(1,n)];
      MT = [zeros(n+2,1) MT zeros(n+2,1)];
      [carga_toxica_ctumoral,B] = disparar_tfm(fila1, col1, carga_toxica_ctumoral, B);
      fila_new = [];
      col_new = [];
      carga_toxica_ctumoral_new = [];
      for i = 1:length(carga_toxica_ctumoral)
        matamos = morir_tfm(carga_toxica_ctumoral(i));
        if matamos == 1
            D(fila1(i)+1, col1(i)+1) = 0;
            D(fila1(i)+2, col1(i)+1) = 0;
            D(fila1(i)+1, col1(i)+2) = 0;
            D(fila1(i)+2, col1(i)+2) = 0;

            if A(fila1(i)+1, col1(i)+1) == 2
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
        else
            fila_new = [fila_new fila1(i)];            
            col_new = [col_new col1(i)];
            carga_toxica_ctumoral_new = [carga_toxica_ctumoral_new carga_toxica_ctumoral(i)];
        end
      end      
      fila1 = fila_new;
      col1 = col_new;
      carga_toxica_ctumoral = carga_toxica_ctumoral_new;
      
      %Con 'mover' pierdo el doble borde:
      [A,B,MT] = mover_tfm(A(2:end-1,2:end-1),B(2:end-1,2:end-1),MT(2:end-1,2:end-1));
      [A,B,D,MT,fila1,col1] = proliferar_tfm(fila1,col1,A,B,D(2:end-1,2:end-1),MT,0.01);
      
      %%Las posiciones donde haya 1, esa posicion en A y en B desaparecen 
      [fila_adios,col_adios] = find(MT==1);
      for dis = 1:length(col_adios)
          if A(fila_adios(dis), col_adios(dis)) == 1
            A(fila_adios(dis), col_adios(dis)) = 0;
          else
            A(fila_adios(dis), col_adios(dis)) = 3;
          end
          B(fila_adios(dis), col_adios(dis)) = 0;
      end

      MT = MT - (MT>0); %Resto un día de vida
      %Vuelvo a añadir el doble borde:
      A = [zeros(1,n); A; zeros(1,n)];
      A = [zeros(n+2,1) A zeros(n+2,1)];
      B = [zeros(1,n); B; zeros(1,n)];
      B = [zeros(n+2,1) B zeros(n+2,1)];
      D = [zeros(1,n); D; zeros(1,n)];
      D = [zeros(n+2,1) D zeros(n+2,1)];
      MT = [zeros(1,n); MT; zeros(1,n)];
      MT = [zeros(n+2,1) MT zeros(n+2,1)];
      %drawnow
      %pause(0.05)
      evolucion(tiempo) = length(fila1);
end
end
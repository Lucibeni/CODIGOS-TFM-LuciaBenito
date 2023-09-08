function [AA,BB,DD,MMTT,FILA,COLUMNA] = proliferar_tfm(fila, columna, A, B, D, MT, rho)
    %ES LA FUNCIÓN ENCARGADA DE LA PROLIFERACIÓN DE LAS CÉLULAS TUMORALES
    
    % rho es el parámetro de la tasa de proliferación
    % A, B, D, MT, fila y columna son las ya definidas en el código
    % principal (haciendo referencia fila y columna a fila1 y col1)
    % Meter las matrices sin los bordes en el principal 
    
    % Ampliamos de nuevo los bordes dentro de esta funcion:
    n = length(A);
    A = [2*ones(1,n); 2*ones(1,n); A; 2*ones(1,n); 2*ones(1,n)];
    A = [2*ones(n+4,1) 2*ones(n+4,1) A 2*ones(n+4,1) 2*ones(n+4,1)];
    D = [2*ones(1,n); 2*ones(1,n); D; 2*ones(1,n); 2*ones(1,n)];
    D = [2*ones(n+4,1) 2*ones(n+4,1) D 2*ones(n+4,1) 2*ones(n+4,1)];
    B = [2*ones(1,n); 2*ones(1,n); B; 2*ones(1,n); 2*ones(1,n)];
    B = [2*ones(n+4,1) 2*ones(n+4,1) B 2*ones(n+4,1) 2*ones(n+4,1)];
    MT = [2*ones(1,n); 2*ones(1,n); MT; 2*ones(1,n); 2*ones(1,n)];
    MT = [2*ones(n+4,1) 2*ones(n+4,1) MT 2*ones(n+4,1) 2*ones(n+4,1)];
    %(Se escoge esta ampliación de doble borde de 2s para que si una célula 
    % tumoral quiere proliferar fuera del espacio, lo detecte como si ya hubiese
    % una célula tumoral en esa posición y no escoja esa posición para proliferar)

    % Como hemos modificado el tamaño de A, hay que modificar tambien donde
    % estan las células tumorales, que es en dos posiciones más de la
    % original
    for a=1:length(fila)
        fila(a)=fila(a)+2;
        columna(a)=columna(a)+2;
    end

    longitud=length(fila); %Da el número de células tumorales en el sistema
    contador_vector=longitud; %Se utilizará para indicar la posición del vector en la que insertaremos la posición de una nueva célula tumoral
    for j=1:longitud %El proceso se sigue para todas las células tumorales
        prolifera = Bernu_tfm(rho); %Decidimos si la célula prolifera o no utilizando una Bernoulli
        if prolifera == 1 %Si la célula prolifera, debemos de elegir dónde, teniendo como opción un máximo de 16 posiciones
            F = zeros(1,2); %Este vector guardará la fila de la celda superior izquierda de cada posible posición en la que pueda surgir una nueva célula tumoral
            C = zeros(1,2); %Este vector guardará la columna de la celda superior izquierda de cada posible posición en la que pueda surgir una nueva célula tumoral
            count=0; %Contador de en cuántas posibles posiciones puede proliferar una célula tumoral
            
            %En los siguientes 'if' lo que se hace es comprobar si la
            %posición contigua correspondiente está completamente libre. 
            % Si lo está, el contador suma 1 y los vectores F y C guardan
            % la fila y columna de la celda superior izquierda de la
            % posible posición donde puede surgir una nueva célula tumoral

            %1%%%%%%%%%%%%%%%%%%%%%%%%
            if D(fila(j)-2,columna(j)-2) ~=2 & D(fila(j)-2,columna(j)-1) ~=2 & D(fila(j)-1,columna(j)-2) ~=2 & D(fila(j)-1,columna(j)-1) ~=2
                    count = count + 1;
                    F(count) = fila(j)-2;
                    C(count) = columna(j)-2;
            end

            %2%%%%%%%%%%%%%%%
            if D(fila(j)-2,columna(j)-1) ~=2 & D(fila(j)-2,columna(j)) ~=2 & D(fila(j)-1,columna(j)-1) ~=2 & D(fila(j)-1,columna(j)) ~=2
                    count = count + 1;
                    F(count) = fila(j)-2;
                    C(count) = columna(j)-1;
            end
            
            %3%%%%%%%%%%%%%%%%%%
            if D(fila(j)-2,columna(j)) ~=2 & D(fila(j)-2,columna(j)+1) ~=2 & D(fila(j)-1,columna(j)) ~=2 & D(fila(j)-1,columna(j)+1) ~=2
                    count = count + 1;
                    F(count) = fila(j)-2;
                    C(count) = columna(j);
            end
            
            %4%%%%%%%%%%%%%%%%%%%%%
            if D(fila(j)-2,columna(j)+1) ~=2 & D(fila(j)-2,columna(j)+2) ~=2 & D(fila(j)-1,columna(j)+1) ~=2 & D(fila(j)-1,columna(j)+2) ~=2
                    count = count + 1;
                    F(count) = fila(j)-2;
                    C(count) = columna(j)+1;
            end
            
            %5%%%%%%%%%%%%%%%%%%
            if D(fila(j)-2,columna(j)+2) ~=2 & D(fila(j)-2,columna(j)+3) ~=2 & D(fila(j)-1,columna(j)+2) ~=2 & D(fila(j)-1,columna(j)+3) ~=2
                    count = count + 1;
                    F(count) = fila(j)-2;
                    C(count) = columna(j)+2;
            end
            
            %6%%%%%%%%%%%%%%%
             if D(fila(j)-1,columna(j)+2) ~=2 & D(fila(j)-1,columna(j)+3) ~=2 & D(fila(j),columna(j)+2) ~=2 & D(fila(j),columna(j)+3) ~=2
                    count = count + 1;
                    F(count) = fila(j)-1;
                    C(count) = columna(j)+2;
             end
            
             %7%%%%%%%%%%%%%%%
             if D(fila(j),columna(j)+2) ~=2 & D(fila(j),columna(j)+3) ~=2 & D(fila(j)+1,columna(j)+2) ~=2 & D(fila(j)+1,columna(j)+3) ~=2
                    count = count + 1;
                    F(count) = fila(j);
                    C(count) = columna(j)+2;
             end
             
             %8%%%%%%%%%%%%%%%%%%%
             if D(fila(j)+1,columna(j)+2) ~=2 & D(fila(j)+1,columna(j)+3) ~=2 & D(fila(j)+2,columna(j)+2) ~=2 & D(fila(j)+2,columna(j)+3) ~=2
                    count = count + 1;
                    F(count) = fila(j)+1;
                    C(count) = columna(j)+2;
             end
             
             %9%%%%%%%%%%%%%%%%%%%
             if D(fila(j)+2,columna(j)+2) ~=2 & D(fila(j)+2,columna(j)+3) ~=2 & D(fila(j)+3,columna(j)+2) ~=2 & D(fila(j)+3,columna(j)+3) ~=2
                    count = count + 1;
                    F(count) = fila(j)+2;
                    C(count) = columna(j)+2;
             end
             
             %10%%%%%%%%%%%%%%%%%%%
             if D(fila(j)+2,columna(j)+1) ~=2 & D(fila(j)+2,columna(j)+2) ~=2 & D(fila(j)+3,columna(j)+1) ~=2 & D(fila(j)+3,columna(j)+2) ~=2
                    count = count + 1;
                    F(count) = fila(j)+2;
                    C(count) = columna(j)+1;
             end
             
             %11%%%%%%%%%%%%%%%%%%%
             if D(fila(j)+2,columna(j)) ~=2 & D(fila(j)+2,columna(j)+1) ~=2 & D(fila(j)+3,columna(j)) ~=2 & D(fila(j)+3,columna(j)+1) ~=2
                    count = count + 1;
                    F(count) = fila(j)+2;
                    C(count) = columna(j);
             end
             
             %12%%%%%%%%%%%%%%%%%%%
             if D(fila(j)+2,columna(j)-1) ~=2 & D(fila(j)+2,columna(j)) ~=2 & D(fila(j)+3,columna(j)-1) ~=2 & D(fila(j)+3,columna(j)) ~=2
                    count = count + 1;
                    F(count) = fila(j)+2;
                    C(count) = columna(j)-1;
             end
             
             %13%%%%%%%%%%%%%%%%%%%
             if D(fila(j)+2,columna(j)-2) ~=2 & D(fila(j)+2,columna(j)-1) ~=2 & D(fila(j)+3,columna(j)-2) ~=2 & D(fila(j)+3,columna(j)-1) ~=2
                    count = count + 1;
                    F(count) = fila(j)+2;
                    C(count) = columna(j)-2;
             end
             
             %14%%%%%%%%%%%%%%%%%%%
             if D(fila(j)+1,columna(j)-2) ~=2 & D(fila(j)+1,columna(j)-1) ~=2 & D(fila(j)+2,columna(j)-2) ~=2 & D(fila(j)+2,columna(j)-1) ~=2
                    count = count + 1;
                    F(count) = fila(j)+1;
                    C(count) = columna(j)-2;
             end
             
             %15%%%%%%%%%%%%%%%%%%%
             if D(fila(j),columna(j)-2) ~=2 & D(fila(j),columna(j)-1) ~=2 & D(fila(j)+1,columna(j)-2) ~=2 & D(fila(j)+1,columna(j)-1) ~=2
                    count = count + 1;
                    F(count) = fila(j);
                    C(count) = columna(j)-2;
             end
             
             %16%%%%%%%%%%%%%%%%%%%
             if D(fila(j)-1,columna(j)-2) ~=2 & D(fila(j)-1,columna(j)-1) ~=2 & D(fila(j),columna(j)-2) ~=2 & D(fila(j),columna(j)-1) ~=2
                    count = count + 1;
                    F(count) = fila(j)-1;
                    C(count) = columna(j)-2;
             end
             
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             
             if count ~= 0 %Si hay alguna posible posición vacía:
                contador_vector = contador_vector + 1; %Este contador suma uno para ser la posición del vector donde se van a guardar la fila y columna de la nueva célula
                k = round(1+(count-1)*rand(1,1)); %Se genera una posición de los vectores donde hemos guardado las posibles posiciones que puede tomar la nueva célula
                fila(contador_vector)=F(k); %Añadimos al vector la fila de la posición de la celda superior izquierda de la nueva célula
                columna(contador_vector)=C(k); %Añadimos al vector la columna de la posición de la celda superior izquierda de la nueva célula
                
                if k > 0 %Nos aseguramos que la posición del vector generada no es el 0
                    %A continuación, vamos a mover las células T que
                    %hubiese en las posiciónes donde ha surgido la nueva
                    %célula tumoral

                    %1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    cont_moverT=0; %Contador del número de posibles posiciones a las que se puede mover la célula T
                    FF=zeros(1,2); %Guarda las filas de las posibles posiciones a  las que se puede mover la célula T
                    CC=zeros(1,2); %Guarda las columnas de las posibles posiciones a  las que se puede mover la célula T
                    if A(F(k),C(k)) == 1 %Si en esta posición donde ha surgido la célula tumoral había una célula T:
                    %En los siguientes 'if' lo que se hace es comprobar si la
                    %posición contigua correspondiente está libre. 
                    % Si lo está, el contador suma 1 y los vectores F y C guardan
                    % la fila y columna la posible posición a donde se puede 
                    % mover la célula T

                        %1%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k)-1;
                        end
                        
                        %2%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k);
                        end
                        
                        %3%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)+1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k)+1;
                        end
                        
                        %4%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k)+2;
                        end
                        
                        %5%%%%%%%%%%%%%%%
                        if A(F(k),C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k);
                            CC(cont_moverT)=C(k)+2;
                        end
                        
                        %6%%%%%%%%%%%%%%%
                        if A(F(k)+1,C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+1;
                            CC(cont_moverT)=C(k)+2;
                        end
                        
                        %7%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k)+2;
                        end
                        
                        %8%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)+1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k)+1;
                        end
                        
                        %9%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k);
                        end
                        
                        %10%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k)-1;
                        end
                        
                        %11%%%%%%%%%%%%%%%
                        if A(F(k)+1,C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+1;
                            CC(cont_moverT)=C(k)-1;
                        end
                        
                        %12%%%%%%%%%%%%%%%
                        if A(F(k),C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k);
                            CC(cont_moverT)=C(k)-1;
                        end
                        
                        if cont_moverT>0 %Si hay posibles posiciones para moverse:
                            kk=round(1+(cont_moverT-1)*rand(1,1)); %Se genera la posición del vector donde se guardan las posibles posiciones
                            A(FF(kk),CC(kk)) = 1; %Se inserta en ella la célula T
                            B(FF(kk),CC(kk)) = B(F(k),C(k)); %Este cambio también se produce en la matriz que guarda su carga citotóxica
                            B(F(k),C(k))=0; 
                            MT(FF(kk),CC(kk)) = MT(F(k),C(k)); %Este cambio también se produce en la matriz que guarda las horas de vida que le quedan                        
                        end
                        MT(F(k),C(k))=0;
                    end

                    %%%%% Se realiza el mismo proceso para las posiciones 
                    % de la célula tumoral nueva:

                    %2%%%%%%%% 
                    cont_moverT=0;
                    if A(F(k),C(k)+1) == 1
                        %1%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k)-1;
                        end
                        %2%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k);
                        end
                        %3%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)+1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k)+1;
                        end
                        %4%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k)+2;
                        end
                        %5%%%%%%%%%%%%%%%
                        if A(F(k),C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k);
                            CC(cont_moverT)=C(k)+2;
                        end
                        %6%%%%%%%%%%%%%%%
                        if A(F(k)+1,C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+1;
                            CC(cont_moverT)=C(k)+2;
                        end
                        %7%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k)+2;
                        end
                        %8%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)+1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k)+1;
                        end
                        %9%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k);
                        end
                        %10%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k)-1;
                        end
                        %11%%%%%%%%%%%%%%%
                        if A(F(k)+1,C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+1;
                            CC(cont_moverT)=C(k)-1;
                        end
                        %12%%%%%%%%%%%%%%%
                        if A(F(k),C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k);
                            CC(cont_moverT)=C(k)-1;
                        end
                        if cont_moverT>0
                            kk=round((cont_moverT-1)*rand(1,1))+1;
                            A(FF(kk),CC(kk)) = 1;
                            B(FF(kk),CC(kk)) = B(F(k),C(k)+1);
                            B(F(k),C(k)+1)=0;
                            MT(FF(kk),CC(kk)) = MT(F(k),C(k)+1);                            
                        end
                        MT(F(k),C(k)+1)=0;
                    end

                    %3%%%%%%%%
                    cont_moverT=0;
                    if A(F(k)+1,C(k)) == 1
                        %1%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k)-1;
                        end
                        %2%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k);
                        end
                        %3%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)+1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k)+1;
                        end
                        %4%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k)+2;
                        end
                        %5%%%%%%%%%%%%%%%
                        if A(F(k),C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k);
                            CC(cont_moverT)=C(k)+2;
                        end
                        %6%%%%%%%%%%%%%%%
                        if A(F(k)+1,C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+1;
                            CC(cont_moverT)=C(k)+2;
                        end
                        %7%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k)+2;
                        end
                        %8%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)+1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k)+1;
                        end
                        %9%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k);
                        end
                        %10%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k)-1;
                        end
                        %11%%%%%%%%%%%%%%%
                        if A(F(k)+1,C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+1;
                            CC(cont_moverT)=C(k)-1;
                        end
                        %12%%%%%%%%%%%%%%%
                        if A(F(k),C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k);
                            CC(cont_moverT)=C(k)-1;
                        end
                        if cont_moverT>0
                            kk=round((cont_moverT-1)*rand(1,1))+1;
                            A(FF(kk),CC(kk)) = 1;
                            B(FF(kk),CC(kk)) = B(F(k)+1,C(k));
                            B(F(k)+1,C(k))=0;
                            MT(FF(kk),CC(kk)) = MT(F(k)+1,C(k));
                        end
                        MT(F(k)+1,C(k))=0;
                    end

                    %4%%%%%%%%
                    cont_moverT=0;
                    if A(F(k)+1,C(k)+1) == 1
                        %1%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k)-1;
                        end
                        %2%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k);
                        end
                        %3%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)+1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k)+1;
                        end
                        %4%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k)+2;
                        end
                        %5%%%%%%%%%%%%%%%
                        if A(F(k),C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k);
                            CC(cont_moverT)=C(k)+2;
                        end
                        %6%%%%%%%%%%%%%%%
                        if A(F(k)+1,C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+1;
                            CC(cont_moverT)=C(k)+2;
                        end
                        %7%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k)+2;
                        end
                        %8%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)+1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k)+1;
                        end
                        %9%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k);
                        end
                        %10%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k)-1;
                        end
                        %11%%%%%%%%%%%%%%%
                        if A(F(k)+1,C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+1;
                            CC(cont_moverT)=C(k)-1;
                        end
                        %12%%%%%%%%%%%%%%%
                        if A(F(k),C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k);
                            CC(cont_moverT)=C(k)-1;
                        end
                        if cont_moverT>0
                            kk=round(1+(cont_moverT-1)*rand(1,1));
                            A(FF(kk),CC(kk)) = 1;
                            B(FF(kk),CC(kk)) = B(F(k)+1,C(k)+1);
                            B(F(k)+1,C(k)+1)=0;
                            MT(FF(kk),CC(kk)) = MT(F(k)+1,C(k)+1);
                        end
                        MT(F(k)+1,C(k)+1)=0;
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %El siguiente paso es insertar la célula tumoral. En la
                    %matriz D se hará siempre en todas las posiciones. Sin 
                    % embargo, en la matriz A solo se cambiaran las celdas 
                    % que no correspondan a un vaso sanguíneo (ya sea con 
                    % célula T o el vaso vacío):

                    if A(F(k),C(k)) ~= 3 & A(F(k),C(k)) ~= 1.25
                        A(F(k),C(k)) = 2;
                    end
                    if A(F(k),C(k)+1) ~= 3 & A(F(k),C(k)+1) ~= 1.25
                        A(F(k),C(k)+1) = 2;
                    end
                    if A(F(k)+1,C(k)) ~= 3 & A(F(k)+1,C(k)) ~= 1.25
                        A(F(k)+1,C(k)) = 2;
                    end
                    if A(F(k)+1,C(k)+1) ~= 3 & A(F(k)+1,C(k)+1) ~= 1.25
                        A(F(k)+1,C(k)+1) = 2;
                    end
                    D(F(k),C(k)) = 2;
                    D(F(k),C(k)+1) = 2;
                    D(F(k)+1,C(k)) = 2;
                    D(F(k)+1,C(k)+1) = 2;                    
                end                 
             end
        end
    end
    %Volvemos a cambiar los siguientes vectores para considerar las
    %posiciones sin los dobles bordes artificiales
    for a=1:length(fila)
        fila(a)=fila(a)-2;
        columna(a)=columna(a)-2;
    end
    %Las matrices se devuelven los dobles bordes
    AA = A(3:end-2,3:end-2);
    BB = B(3:end-2,3:end-2);
    DD = D(3:end-2,3:end-2);
    MMTT = MT(3:end-2,3:end-2);
    FILA = fila;
    COLUMNA = columna;
end
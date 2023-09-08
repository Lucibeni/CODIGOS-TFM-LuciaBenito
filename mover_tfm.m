function [nuevoA, nuevoB, nuevoMT] = mover(A,B,MT)
    %ES LA FUNCIÓN ENCARGADA DE LOS MOVIMIENTOS DE LAS CÉLULAS T POR EL
    %ESPACIO INTERCELULAR
    
        n = length(A); %La matriz entra con bordes
        [filaM,columnaM] = find(A==1); %Los vectores filaM y columnaM guardan las filas y columnas de las celdas donde se encuentran las células T que ya han salido de los vasos sanguíneos
        for j = 1:length(filaM)
            F = zeros(1,2); %Este vector guardará la fila de cada posible posición a la que pueda mover la célula T
            C = zeros(1,2); %Este vector guardará la columna de cada posible posición a la que pueda mover la célula T
            count=0; %Contador de a cuántas posibles posiciones puede moverse una célula T (máximo 8)
            
            %Considerando a parte los casos especiales para los cuales la 
            %célula T podría salirse del autómata, en los siguientes 'if' lo
            % que se va a hacer es comprobar si las celdas contiguas a la 
            % célula T que son parte del autómata están libres. 
            % Por cada celda libre, el contador suma 1 y los vectores F y C 
            % guardan las filas y columnas de las posibles posiciones a las
            % que puede moverse la célula T

            % Si la célula no se en cuentra en los bordes del sistema:
            if filaM(j) > 1 & columnaM(j) > 1 & filaM(j) < n & columnaM(j) < n
                if A(filaM(j)-1,columnaM(j)-1) == 0
                    count = count + 1;
                    F(count) = filaM(j)-1;
                    C(count) = columnaM(j) -1;
                end
                if A(filaM(j)-1,columnaM(j)) == 0
                    count = count + 1;
                    F(count) = filaM(j)-1;
                    C(count) = columnaM(j);
                end
                if A(filaM(j)-1,columnaM(j)+1) == 0
                    count = count + 1;
                    F(count) = filaM(j)-1;
                    C(count) = columnaM(j) +1;
                end
                if A(filaM(j),columnaM(j)-1) == 0
                    count = count + 1;
                    F(count) = filaM(j);
                    C(count) = columnaM(j) -1;
                end
                if A(filaM(j),columnaM(j)+1) == 0
                    count = count + 1;
                    F(count) = filaM(j);
                    C(count) = columnaM(j) +1;
                end
                if A(filaM(j)+1,columnaM(j)-1) == 0
                    count = count + 1;
                    F(count) = filaM(j)+1;
                    C(count) = columnaM(j) -1;
                end
                if A(filaM(j)+1,columnaM(j)) == 0
                    count = count + 1;
                    F(count) = filaM(j)+1;
                    C(count) = columnaM(j) ;
                end
                if A(filaM(j)+1,columnaM(j)+1) == 0
                    count = count + 1;
                    F(count) = filaM(j)+1;
                    C(count) = columnaM(j) + 1;
                end
                if count > 0 %Si hay alguna posible posición vacía:
                    k = round((count-1)*rand(1,1))+1; %Se genera una posición de los vectores donde hemos guardado las posibles posiciones a las que puede moverse la célula T
                    A(filaM(j), columnaM(j)) = 0; %La posición en la que estaba queda libre
                    A(F(k),C(k)) = 1; %La célula se desplaza a una nueva posición
                    if B(filaM(j), columnaM(j)) < 9 %Si su carga citotóxica no es la dosis máxima:
                        B(F(k),C(k)) = B(filaM(j), columnaM(j))+1; %Su carga citotóxica se recarga sumando 1 y realiza el movimiento en la matriz B también
                    else
                        B(F(k),C(k)) = B(filaM(j), columnaM(j)); %Sino, se deja el mismo valor pero realiza el movimiento en la matriz B también
                    end
                    B(filaM(j), columnaM(j)) = 0;
                    MT(F(k),C(k)) = MT(filaM(j), columnaM(j)); %El movimiento también se realiza para su correspondiente esperanza de vida en la matriz MT
                    MT(filaM(j), columnaM(j)) = 0;
                end
            end
            
            %Para los siguientes casos, se procede del mismo modo:

            %Si la célula T se encuentra en la esquina superior izquierda:
            count = 0;
            if filaM(j) == 1 & columnaM(j) == 1
                if A(filaM(j),columnaM(j)+1) == 0
                    count = count + 1;
                    F(count) = filaM(j);
                    C(count) = columnaM(j) +1;
                end
                if A(filaM(j)+1,columnaM(j)) == 0
                    count = count + 1;
                    F(count) = filaM(j)+1;
                    C(count) = columnaM(j);
                end
                if A(filaM(j)+1,columnaM(j)+1) == 0
                    count = count + 1;
                    F(count) = filaM(j)+1;
                    C(count) = columnaM(j) +1;
                end
                if count > 0
                    k = round((count-1)*rand(1,1))+1;
                    A(filaM(j), columnaM(j)) = 0;
                    A(F(k),C(k)) = 1;
                    if B(filaM(j), columnaM(j)) < 9
                        B(F(k),C(k)) = B(filaM(j), columnaM(j))+1;
                    else
                        B(F(k),C(k)) = B(filaM(j), columnaM(j));
                    end
                    B(filaM(j), columnaM(j)) = 0;
                    MT(F(k),C(k)) = MT(filaM(j), columnaM(j));
                    MT(filaM(j), columnaM(j)) = 0;
                end
            end


            %Si la célula T se encuentra en la esquina superior derecha:
            count = 0;
            if filaM(j) == 1 & columnaM(j) == n
                if A(filaM(j),columnaM(j)-1) == 0
                    count = count + 1;
                    F(count) = filaM(j);
                    C(count) = columnaM(j) - 1;
                end
                if A(filaM(j)+1,columnaM(j)-1) == 0
                    count = count + 1;
                    F(count) = filaM(j)+1;
                    C(count) = columnaM(j)-1;
                end
                if A(filaM(j)+1,columnaM(j)) == 0
                    count = count + 1;
                    F(count) = filaM(j)+1;
                    C(count) = columnaM(j);
                end
                if count > 0
                    k = round((count-1)*rand(1,1))+1;
                    A(filaM(j), columnaM(j)) = 0;
                    A(F(k),C(k)) = 1;
                    if B(filaM(j), columnaM(j)) < 9
                        B(F(k),C(k)) = B(filaM(j), columnaM(j))+1;
                    else
                        B(F(k),C(k)) = B(filaM(j), columnaM(j));
                    end
                    B(filaM(j), columnaM(j)) = 0;
                    MT(F(k),C(k)) = MT(filaM(j), columnaM(j));
                    MT(filaM(j), columnaM(j)) = 0;
                end
            end


            %Si la célula T se encuentra en la esquina inferior izquierda:
            count = 0;
            if filaM(j) == n & columnaM(j) == 1
                if A(filaM(j)-1,columnaM(j)) == 0
                    count = count + 1;
                    F(count) = filaM(j)-1;
                    C(count) = columnaM(j);
                end
                if A(filaM(j)-1,columnaM(j)+1) == 0
                    count = count + 1;
                    F(count) = filaM(j)-1;
                    C(count) = columnaM(j)+1;
                end
                if A(filaM(j),columnaM(j)+1) == 0
                    count = count + 1;
                    F(count) = filaM(j);
                    C(count) = columnaM(j)+1;
                end
                if count > 0
                    k = round((count-1)*rand(1,1))+1;
                    A(filaM(j), columnaM(j)) = 0;
                    A(F(k),C(k)) = 1;
                    if B(filaM(j), columnaM(j)) < 9
                        B(F(k),C(k)) = B(filaM(j), columnaM(j))+1;
                    else
                        B(F(k),C(k)) = B(filaM(j), columnaM(j));
                    end
                    B(filaM(j), columnaM(j)) = 0;
                    MT(F(k),C(k)) = MT(filaM(j), columnaM(j));
                    MT(filaM(j), columnaM(j)) = 0;
                end
            end


            %Si la célula T se encuentra en la esquina inferior derecha:
            count = 0;
            if filaM(j) == n & columnaM(j) == n
                if A(filaM(j)-1,columnaM(j)-1) == 0
                    count = count + 1;
                    F(count) = filaM(j)-1;
                    C(count) = columnaM(j)-1;
                end
                if A(filaM(j)-1,columnaM(j)) == 0
                    count = count + 1;
                    F(count) = filaM(j)-1;
                    C(count) = columnaM(j);
                end
                if A(filaM(j),columnaM(j)-1) == 0
                    count = count + 1;
                    F(count) = filaM(j);
                    C(count) = columnaM(j)-1;
                end
                if count > 0
                    k = round((count-1)*rand(1,1))+1;
                    A(filaM(j), columnaM(j)) = 0;
                    A(F(k),C(k)) = 1;
                    if B(filaM(j), columnaM(j)) < 9
                        B(F(k),C(k)) = B(filaM(j), columnaM(j))+1;
                    else
                        B(F(k),C(k)) = B(filaM(j), columnaM(j));
                    end
                    B(filaM(j), columnaM(j)) = 0;
                    MT(F(k),C(k)) = MT(filaM(j), columnaM(j));
                    MT(filaM(j), columnaM(j)) = 0;
                end
            end


            %Si la célula T se encuentra en el borde superior:
            count = 0;
            if filaM(j) == 1 & columnaM(j) > 1 & columnaM(j) < n
                if A(filaM(j),columnaM(j)-1) == 0
                    count = count + 1;
                    F(count) = filaM(j);
                    C(count) = columnaM(j) -1;
                end
                if A(filaM(j),columnaM(j)+1) == 0
                    count = count + 1;
                    F(count) = filaM(j);
                    C(count) = columnaM(j) +1;
                end
                if A(filaM(j)+1,columnaM(j)-1) == 0
                    count = count + 1;
                    F(count) = filaM(j)+1;
                    C(count) = columnaM(j) -1;
                end
                if A(filaM(j)+1,columnaM(j)) == 0
                    count = count + 1;
                    F(count) = filaM(j)+1;
                    C(count) = columnaM(j) ;
                end
                if A(filaM(j)+1,columnaM(j)+1) == 0
                    count = count + 1;
                    F(count) = filaM(j)+1;
                    C(count) = columnaM(j) + 1;
                end
                if count > 0
                    k = round((count-1)*rand(1,1))+1;
                    A(filaM(j), columnaM(j)) = 0;
                    A(F(k),C(k)) = 1;
                    if B(filaM(j), columnaM(j)) < 9
                        B(F(k),C(k)) = B(filaM(j), columnaM(j))+1;
                    else
                        B(F(k),C(k)) = B(filaM(j), columnaM(j));
                    end
                    B(filaM(j), columnaM(j)) = 0;
                    MT(F(k),C(k)) = MT(filaM(j), columnaM(j));
                    MT(filaM(j), columnaM(j)) = 0;
                end
            end


            %Si la célula T se encuentra en el borde inferior:
            count = 0;
            if filaM(j) == n & columnaM(j) > 1 & columnaM(j) < n
                if A(filaM(j)-1,columnaM(j)-1) == 0
                    count = count + 1;
                    F(count) = filaM(j)-1;
                    C(count) = columnaM(j) -1;
                end
                if A(filaM(j)-1,columnaM(j)) == 0
                    count = count + 1;
                    F(count) = filaM(j)-1;
                    C(count) = columnaM(j);
                end
                if A(filaM(j)-1,columnaM(j)+1) == 0
                    count = count + 1;
                    F(count) = filaM(j)-1;
                    C(count) = columnaM(j) +1;
                end
                if A(filaM(j),columnaM(j)-1) == 0
                    count = count + 1;
                    F(count) = filaM(j);
                    C(count) = columnaM(j) -1;
                end
                if A(filaM(j),columnaM(j)+1) == 0
                    count = count + 1;
                    F(count) = filaM(j);
                    C(count) = columnaM(j) +1;
                end
                if count > 0
                    k = round((count-1)*rand(1,1))+1;
                    A(filaM(j), columnaM(j)) = 0;
                    A(F(k),C(k)) = 1;
                    if B(filaM(j), columnaM(j)) < 9
                        B(F(k),C(k)) = B(filaM(j), columnaM(j))+1;
                    else
                        B(F(k),C(k)) = B(filaM(j), columnaM(j));
                    end
                    B(filaM(j), columnaM(j)) = 0;
                    MT(F(k),C(k)) = MT(filaM(j), columnaM(j));
                    MT(filaM(j), columnaM(j)) = 0;
                end
            end


            %Si la célula T se encuentra en el borde izquierdo:
            count = 0;
            if columnaM(j) == 1 & filaM(j) < n & filaM(j) > 1
                if A(filaM(j)-1,columnaM(j)) == 0
                    count = count + 1;
                    F(count) = filaM(j)-1;
                    C(count) = columnaM(j);
                end
                if A(filaM(j)-1,columnaM(j)+1) == 0
                    count = count + 1;
                    F(count) = filaM(j)-1;
                    C(count) = columnaM(j) +1;
                end
                if A(filaM(j),columnaM(j)+1) == 0
                    count = count + 1;
                    F(count) = filaM(j);
                    C(count) = columnaM(j) +1;
                end
                if A(filaM(j)+1,columnaM(j)) == 0
                    count = count + 1;
                    F(count) = filaM(j)+1;
                    C(count) = columnaM(j) ;
                end
                if A(filaM(j)+1,columnaM(j)+1) == 0
                    count = count + 1;
                    F(count) = filaM(j) + 1;
                    C(count) = columnaM(j) + 1;
                end
                if count > 0
                    k = round((count-1)*rand(1,1))+1;
                    A(filaM(j), columnaM(j)) = 0;
                    A(F(k),C(k)) = 1;
                    if B(filaM(j), columnaM(j)) < 9
                        B(F(k),C(k)) = B(filaM(j), columnaM(j))+1;
                    else
                        B(F(k),C(k)) = B(filaM(j), columnaM(j));
                    end
                    B(filaM(j), columnaM(j)) = 0;
                    MT(F(k),C(k)) = MT(filaM(j), columnaM(j));
                    MT(filaM(j), columnaM(j)) = 0;
                end
            end


            %Si la célula T se encuentra en el borde derecho:
            count = 0;
            if columnaM(j) == n & filaM(j) < n & filaM(j) > 1
                if A(filaM(j)-1,columnaM(j)-1) == 0
                    count = count + 1;
                    F(count) = filaM(j)-1;
                    C(count) = columnaM(j) -1;
                end
                if A(filaM(j)-1,columnaM(j)) == 0
                    count = count + 1;
                    F(count) = filaM(j)-1;
                    C(count) = columnaM(j);
                end
                if A(filaM(j),columnaM(j)-1) == 0
                    count = count + 1;
                    F(count) = filaM(j);
                    C(count) = columnaM(j) -1;
                end
                if A(filaM(j)+1,columnaM(j)-1) == 0
                    count = count + 1;
                    F(count) = filaM(j)+1;
                    C(count) = columnaM(j) -1;
                end
                if A(filaM(j)+1,columnaM(j)) == 0
                    count = count + 1;
                    F(count) = filaM(j)+1;
                    C(count) = columnaM(j) ;
                end
                if count > 0
                    k = round((count-1)*rand(1,1))+1;
                    A(filaM(j), columnaM(j)) = 0;
                    A(F(k),C(k)) = 1;
                    if B(filaM(j), columnaM(j)) < 9
                        B(F(k),C(k)) = B(filaM(j), columnaM(j))+1;
                    else
                        B(F(k),C(k)) = B(filaM(j), columnaM(j));
                    end
                    B(filaM(j), columnaM(j)) = 0;
                    MT(F(k),C(k)) = MT(filaM(j), columnaM(j));
                    MT(filaM(j), columnaM(j)) = 0;
                end
            end
        end
    nuevoA = A;
    nuevoB = B;
    nuevoMT = MT;
end
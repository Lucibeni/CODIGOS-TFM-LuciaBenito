function [nuevoA, nuevoMT] = moverDentro_tfm(A,MT)
    %ES LA FUNCIÓN ENCARGADA DEL AVANCE DE LAS CÉLULAS T HACIA EL TUMOR
    %DENTRO DE LOS VASOS SANGUÍNEOS

        n = length(A); %Matriz sin bordes

        [filaM,columnaM] = find(A==1.25); %En los vectores filaM y columnaM se guardan las posiciones de las células T que se encuentran dentro del vaso sanguíneo
        for j = 1:length(filaM)
            F = zeros(1,2); %Este vector guardará la fila de cada posible posición a la que pueda desplazar la célula T
            C = zeros(1,2); %Este vector guardará la columna de cada posible posición a la que pueda desplazar la célula T
            count=0; %Contador de a cuántas posibles posiciones puede desplazarse una célula T
            
            % En los siguientes 'if' lo que se hace es comprobar si las
            % celdas del vaso sanguíneo contiguas a la célula T que no estén más alejadas
            % del tumor de lo que se encuentra la célula en ese momento están
            % libres. 
            % Por cada celda libre, el contador suma 1 y los vectores F y C guardan
            % las filas y columnas de las posibles posiciones a las que
            % puede avanzar la célula T

            if filaM(j) > 1 & columnaM(j) > 1 & filaM(j) < n & columnaM(j) < n 
                %1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if A(filaM(j)-1,columnaM(j)-1) == 3
                    count = count + 1;
                    F(count) = filaM(j)-1;
                    C(count) = columnaM(j) -1;
                end

                %2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if A(filaM(j)-1,columnaM(j)) == 3
                    count = count + 1;
                    F(count) = filaM(j)-1;
                    C(count) = columnaM(j);
                end
                
                %3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if A(filaM(j),columnaM(j)-1) == 3
                    count = count + 1;
                    F(count) = filaM(j);
                    C(count) = columnaM(j) -1;
                end

                %4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if A(filaM(j)+1,columnaM(j)-1) == 3
                    count = count + 1;
                    F(count) = filaM(j)+1;
                    C(count) = columnaM(j) -1;
                end

                %5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if A(filaM(j)+1,columnaM(j)) == 3
                    count = count + 1;
                    F(count) = filaM(j)+1;
                    C(count) = columnaM(j) ;
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if count > 0 %Si hay alguna posible posición vacía:
                    k = round((count-1)*rand(1,1))+1; %Se genera una posición de los vectores donde hemos guardado las posibles posiciones a las que puede avanzar la célula T
                    A(filaM(j), columnaM(j)) = 3; %La posición en la que estaba queda libre, se convierte en vaso sanguíneo
                    A(F(k),C(k)) = 1.25; %La célula se desplaza por el vaso sanguíneo a una nueva posición, sigue tomando el valor 1.25
                    MT(F(k),C(k)) = MT(filaM(j), columnaM(j)); %El movimiento también se aplica en la matriz MT en la que a cada célula T le acompaña su esperanza de vida
                    MT(filaM(j), columnaM(j)) = 0;
                end
            end
            
            %Borde derecho (Mismo proceso para las célula del borde derecho)
            count = 0;
            if columnaM(j) == n & filaM(j) < n & filaM(j) > 1
                %1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if A(filaM(j)-1,columnaM(j)-1) == 3
                    count = count + 1;
                    F(count) = filaM(j)-1;
                    C(count) = columnaM(j) -1;
                end

                %2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if A(filaM(j)-1,columnaM(j)) == 3
                    count = count + 1;
                    F(count) = filaM(j)-1;
                    C(count) = columnaM(j);
                end

                %3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if A(filaM(j),columnaM(j)-1) == 3
                    count = count + 1;
                    F(count) = filaM(j);
                    C(count) = columnaM(j) -1;
                end

                %4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if A(filaM(j)+1,columnaM(j)-1) == 3
                    count = count + 1;
                    F(count) = filaM(j)+1;
                    C(count) = columnaM(j) -1;
                end

                %5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if A(filaM(j)+1,columnaM(j)) == 3
                    count = count + 1;
                    F(count) = filaM(j)+1;
                    C(count) = columnaM(j) ;
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if count > 0
                    k = round((count-1)*rand(1,1))+1;
                    A(filaM(j), columnaM(j)) = 3;
                    A(F(k),C(k)) = 1.25;
                    MT(F(k),C(k)) = MT(filaM(j), columnaM(j));
                    MT(filaM(j), columnaM(j)) = 0;
                end
            end
        end
    nuevoA = A;
    nuevoMT = MT;
end
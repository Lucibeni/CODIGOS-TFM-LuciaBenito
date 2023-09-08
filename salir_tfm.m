function [nuevoA, nuevoB, nuevoMT] = salir_tfm(A,B,MT)
    %ES LA FUNCIÓN ENCARGADA DE LA SALIDA DE LAS CÉLULAS T FUERA DE LOS
    %VASOS SANGUÍNEOS

        n = length(A); %Estoy metiendo la matriz sin los bordes

        [filaTVS,colTVS] = find(A==1.25); %En los vectores filaTVS y colTVS se guardan las posiciones de las células T que se encuentran dentro del vaso sanguíneo
        for j = 1:length(filaTVS)
            sale = Bernu_tfm(0.07); %Primero decido si sale o no del vaso sanguíneo
            if sale == 1 %Si la célula T sale del vaso, hay que decidir hacia dónde
                F = zeros(1,2); %Este vector guardará la fila de cada posible posición a la que pueda salir la célula T
                C = zeros(1,2); %Este vector guardará la columna de cada posible posición a la que pueda salir la célula T
                count=0; %Contador de a cuántas posibles posiciones puede salir una célula T
            
                %En los siguientes 'if' lo que se hace es comprobar si cada
                %posición contigua al vaso sanguíneo está libre. 
                % Si lo está, el contador suma 1 y los vectores F y C guardan
                % la fila y columna de la posible posición a donde puede salir
                % la célula T

                if filaTVS(j) > 1 & colTVS(j) > 1 & filaTVS(j) < n & colTVS(j) < 65 %Tomamos que la columna donde se encuentra la célula T sea <65 porque necesitamos que se encuentre suficientemente cerca del tumor para activarse y salir del vaso
                    
                    %1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if A(filaTVS(j)-1,colTVS(j)-1) == 0
                        count = count + 1;
                        F(count) = filaTVS(j)-1;
                        C(count) = colTVS(j) -1;
                    end
                    
                    %2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if A(filaTVS(j)-1,colTVS(j)) == 0
                        count = count + 1;
                        F(count) = filaTVS(j)-1;
                        C(count) = colTVS(j);
                    end
                    
                    %3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if A(filaTVS(j)-1,colTVS(j)+1) == 0
                        count = count + 1;
                        F(count) = filaTVS(j)-1;
                        C(count) = colTVS(j) +1;
                    end
                    
                    %4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if A(filaTVS(j),colTVS(j)-1) == 0
                        count = count + 1;
                        F(count) = filaTVS(j);
                        C(count) = colTVS(j) -1;
                    end
                    
                    %5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if A(filaTVS(j),colTVS(j)+1) == 0
                        count = count + 1;
                        F(count) = filaTVS(j);
                        C(count) = colTVS(j) +1;
                    end
                    
                    %6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if A(filaTVS(j)+1,colTVS(j)-1) == 0
                        count = count + 1;
                        F(count) = filaTVS(j)+1;
                        C(count) = colTVS(j) -1;
                    end
                    
                    %7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if A(filaTVS(j)+1,colTVS(j)) == 0
                        count = count + 1;
                        F(count) = filaTVS(j)+1;
                        C(count) = colTVS(j) ;
                    end
                    
                    %8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if A(filaTVS(j)+1,colTVS(j)+1) == 0
                        count = count + 1;
                        F(count) = filaTVS(j)+1;
                        C(count) = colTVS(j) + 1;
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    if count > 0 %Si hay alguna posible posición vacía:
                        k = round((count-1)*rand(1,1))+1; %Se genera una posición de los vectores donde hemos guardado las posibles posiciones a las que puede salir la célula T
                        A(filaTVS(j), colTVS(j)) = 3; %La posición donde estaba la célula pasa a ser una posición de vaso sanguíneo
                        A(F(k),C(k)) = 1; %Ahora la célula T está fuera del vaso, por lo que se le asigna el valor 1
                        B(F(k),C(k)) = 9; %Al salir, se le asigna su carga citotóxica completa
                        MT(F(k),C(k)) = MT(filaTVS(j), colTVS(j)); %El cambio del movimiento de la célula T también afecta a la matriz MT, pues cada célula lleva asignada su esperanza de vida y este valor debe seguir los movimientos de su correspondiente célula
                        MT(filaTVS(j), colTVS(j)) = 0;
                    end
                end
            end
        end
    nuevoA = A;
    nuevoB = B;
    nuevoMT = MT;
end
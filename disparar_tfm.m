function [pium,BB] = disparar_tfm(fila, columna, cargatoxicactumoral, B)
    %ES LA FUNCIÓN ENCARGADA DE LOS DISPAROS DE CARGA CITOTÓXICA POR PARTE
    %DE LAS CÉLULAS T A LAS CÉLULAS TUMORALES

    %Hemos metido B con borde, entonces las posiciones de las c.tumorales cambian y
    %las de las células T también. Porque la fila1 y columna1 que metemos están
    %calculadas con el tamaño inicial de la matriz.
    %Por tanto, en los índices de la matriz B habrá un desfase de una
    %posición que hemos tenido que compensar

    %En el siguiente 'for', lo que se lleva a cabo es, para cada célula
    %tumoral, comprobar si tiene células T contiguas a ella. Si en una de
    %las posiciones contiguas (12 en total) hay una célula T, entonces su
    %carga citotóxica se suma como daño en la célula tumoral.

    %Podemos observar que a cada dosis citotóxica que se acumula, se le
    %resta 1, esto se debe a que las dosis citotóxicas se cuantifican de 1 
    %a 8, pero en la matriz B se asignan los valores de 2 a 9 y el valor 1
    %hace referencia a una célula T sin carga citotóxica debido a que ha 
    %disparado en el anterior paso de tiempo y aún no se ha recargado.

    %Las cargas citotóxicas de las células T que han disparado
    %pasan a valer 0 (antes de pasar al siguiente paso de tiempo, 
    %pasarán a valer 1, por ello de considerar de 2 a 9 y restar 1).

        for i = 1:length(fila)

            %1 %%%%%%%%%%%%%%%%%%%%%%%%
            if B(fila(i),columna(i)) >1 %Si la posición vale 1 es una célula T sin carga
                cargatoxicactumoral(i)= cargatoxicactumoral(i)+B(fila(i),columna(i))-1; 
                B(fila(i),columna(i)) = 0; %
            end

            %2 %%%%%%%%%%%%%%%%%%%%%%%%
            if B(fila(i),columna(i)+1) >1
                cargatoxicactumoral(i)= cargatoxicactumoral(i)+B(fila(i),columna(i)+1)-1;
                B(fila(i),columna(i)+1) = 0;
            end

            %3 %%%%%%%%%%%%%%%%%%%%%%%%
            if B(fila(i),columna(i)+2) >1
                cargatoxicactumoral(i)= cargatoxicactumoral(i)+B(fila(i),columna(i)+2)-1;
                B(fila(i),columna(i)+2) = 0;
            end

            %4 %%%%%%%%%%%%%%%%%%%%%%%%
            if B(fila(i),columna(i)+3) >1
                cargatoxicactumoral(i)= cargatoxicactumoral(i)+B(fila(i),columna(i)+3)-1;
                B(fila(i),columna(i)+3) = 0;
            end

            %5 %%%%%%%%%%%%%%%%%%%%%%%%
            if B(fila(i)+1,columna(i)) >1
                cargatoxicactumoral(i)= cargatoxicactumoral(i)+B(fila(i)+1,columna(i))-1;
                B(fila(i)+1,columna(i)) = 0;
            end

            %6 %%%%%%%%%%%%%%%%%%%%%%%%
            if B(fila(i)+1,columna(i)+3) >1
                cargatoxicactumoral(i)= cargatoxicactumoral(i)+B(fila(i)+1,columna(i)+3)-1;
                B(fila(i)+1,columna(i)+3) = 0;
            end
            
            %7 %%%%%%%%%%%%%%%%%%%%%%%%
            if B(fila(i)+2,columna(i)) >1
                cargatoxicactumoral(i)= cargatoxicactumoral(i)+B(fila(i)+2,columna(i))-1;
                B(fila(i)+2,columna(i)) = 0;
            end

            %8 %%%%%%%%%%%%%%%%%%%%%%%%
            if B(fila(i)+2,columna(i)+3) >1
                cargatoxicactumoral(i)= cargatoxicactumoral(i)+B(fila(i)+2,columna(i)+3)-1;
                B(fila(i)+2,columna(i)+3) = 0;
            end

            %9 %%%%%%%%%%%%%%%%%%%%%%%%
            if B(fila(i)+3,columna(i)) >1
                cargatoxicactumoral(i)= cargatoxicactumoral(i)+B(fila(i)+3,columna(i))-1;
                B(fila(i)+3,columna(i)) = 0;
            end

            %10 %%%%%%%%%%%%%%%%%%%%%%%%
            if B(fila(i)+3,columna(i)+1) >1
                cargatoxicactumoral(i)= cargatoxicactumoral(i)+B(fila(i)+3,columna(i)+1)-1;
                B(fila(i)+3,columna(i)+1) = 0;
            end

            %11 %%%%%%%%%%%%%%%%%%%%%%%%
            if B(fila(i)+3,columna(i)+2) >1
                cargatoxicactumoral(i)= cargatoxicactumoral(i)+B(fila(i)+3,columna(i)+2)-1;
                B(fila(i)+3,columna(i)+2) = 0;
            end

            %12 %%%%%%%%%%%%%%%%%%%%%%%%
            if B(fila(i)+3,columna(i)+3) >1
                cargatoxicactumoral(i)= cargatoxicactumoral(i)+B(fila(i)+3,columna(i)+3)-1;
                B(fila(i)+3,columna(i)+3) = 0;
            end
        end
        pium = cargatoxicactumoral;
        BB = B;
    end
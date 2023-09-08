function muerto = morir_tfm(cargatoxicactumoral)
    %ES LA FUNCIÓN ENCARGADA DE LA MUERTE DE LAS CÉLULAS TUMORALES EN
    %FUNCIÓN DE LA CARGA CITOTÓXICA RECIBIDA

    %Cuanto mayor sea la dosis citotóxica recibida por la célula tumoral,
    %mayor probabilidad va a tener de morir

    muerto = 0;
    %Si la carga citotóxica está entre 0 y 1 dosis completa de
    %citotoxicidad, morirá con un 0.05 de probabilidad.
    if cargatoxicactumoral> 0 & cargatoxicactumoral <9 
        muerto = Bernu_tfm(0.05);
    end

    %Si la carga citotóxica está entre 1 y 2 dosis completas de
    %citotoxicidad, morirá con un 0.12 de probabilidad.
    if cargatoxicactumoral> 8 & cargatoxicactumoral <17
        muerto = Bernu_tfm(0.12);
    end

    %Si la carga citotóxica está entre 2 y 4 dosis completas de
    %citotoxicidad, morirá con un 0.5 de probabilidad.
    if cargatoxicactumoral> 16 & cargatoxicactumoral < 33
        muerto = Bernu_tfm(0.5);
    end

    %Si la carga citotóxica está entre 4 y 6 dosis completas de
    %citotoxicidad, morirá con un 0.8 de probabilidad.
    if cargatoxicactumoral> 32 & cargatoxicactumoral < 49
        muerto = Bernu_tfm(0.8);
    end

    %Si la carga citotóxica es mayor que 6 dosis completas de
    %citotoxicidad, su muerte será practicamente segura.
    if cargatoxicactumoral> 48
        muerto = Bernu_tfm(0.99);
    end    
end



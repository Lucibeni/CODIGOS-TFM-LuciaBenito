function muerte = Bernu_tfm(p)
    %Genera muestras de una Bernoulli de parámetro p
    muerte = 0;
    num = rand(1,1);
    if num < p
        muerte = 1;
    end
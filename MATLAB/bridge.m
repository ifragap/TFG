function lgraph = bridge(inputSize,profundidad)
    
    % inputSize=[128,128,1];
    % profundidad=3;
    % Creamos el bridge y lo conectamos al encoder
    lgraph=encoder(inputSize,profundidad);
    res=resBlock(32*(2^(profundidad+1)),"Bridge");
    lgraph = addLayers(lgraph, res);
    lgraph = connectLayers(lgraph,  'Downsample'+string(profundidad), ...
                           "Bridge" + "_bn1");
    % Conectamos la salida de nuestro encoder a la capa adition de nuestro
    % bridge
    lgraph = connectLayers(lgraph, 'Downsample'+string(profundidad), ...
                           "Bridge_add1/in2");

end


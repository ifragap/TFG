function lgraph = decoder(lgraph, profundidad)
    
    % profundidad = 3 (por ejemplo)
    % lgraph=bridge_1;
    % profundidad=3;
    i=profundidad;
    % Hacemos el upsample y lo conectamos a nuestro grafo(Tambien hacemos la concatenacion pero aun no la conectamos)
    upsample = transposedConv2dLayer(2, 32*(2^(i)), 'Stride', 2, 'Cropping', 'same', 'Name', 'Upsample' + string(i));
    concat = depthConcatenationLayer(2, 'Name', 'Concat' + string(i));
    lgraph = addLayers(lgraph, upsample);
    lgraph = addLayers(lgraph, concat);  
    lgraph = connectLayers(lgraph, "Bridge_add1", "Upsample" + string(i));
    % Hacemos el bloque res y conectamos la parte de los grafos
    % correspondientes a la concatenación
    lgraph = connectLayers(lgraph, "profu" + string(i) + "_add1", "Concat" + string(i) + "/in2");
    lgraph = connectLayers(lgraph, "Upsample" + string(i), "Concat" + string(i) + "/in1");
    res = resBlock(32*(2^i), "profu_decoder" + string(i));
    lgraph = addLayers(lgraph, res);
    % Se reduce el numero de filtros de la salida de la concatenación
    convReduce = convolution2dLayer(1, (32*(2^i)), "Padding", "same", "Name", "reduceFilters"+string(i));
    lgraph = addLayers(lgraph, convReduce);
    lgraph = connectLayers(lgraph, "Concat" + string(i), "reduceFilters"+string(i));
    % Conectamos con el bloque RES
    lgraph = connectLayers(lgraph, "reduceFilters"+string(i), "profu_decoder" + string(i) + "_bn1");
    % Conectamos la capa de entrada con la capa de suma del RES
    lgraph = connectLayers(lgraph, "reduceFilters"+string(i), "profu_decoder" + string(i) + "_add1/in2");


    % conv = convolution2dLayer(3, 32*(2^(i)), "Padding", "same", "Name", "ReduceFilters"+ string(i));
    % lgraph = addLayers(lgraph, conv);
    % lgraph = addLayers(lgraph, concat);  
    % lgraph = connectLayers(lgraph, "Conv2d_bridge", "Upsample" + string(i));
    % lgraph = connectLayers(lgraph, "profu" + string(i) + "_add3", "Concat" + string(i) + "/in2");
    % lgraph = connectLayers(lgraph, "Upsample" + string(i), "Concat" + string(i) + "/in1");    
    % res = resBlock(32*(2^i), "profu_decoder" + string(i));
    % lgraph = addLayers(lgraph, res);
    % lgraph = connectLayers(lgraph, "Concat" + string(i), "profu_decoder" + string(i) + "_bn1");
    % lgraph = connectLayers(lgraph, "Concat" + string(i), "ReduceFilters"+ string(i));
    % lgraph = connectLayers(lgraph, "ReduceFilters"+ string(i), "profu_decoder" + string(i) + "_add1/in2");
    % lgraph = connectLayers(lgraph, 'profu_decoder' + string(i) + '_add1', "profu_decoder" + string(i) + "_add2/in2");
    % lgraph = connectLayers(lgraph, 'profu_decoder' + string(i) + '_add2', "profu_decoder" + string(i) + "_add3/in2");

    for i = profundidad-1:-1:1
        % Se hace el upsampling y la capa de concatenación  y lo conectamos con el bloque res anterior 
        upsample = transposedConv2dLayer(2, 32*(2^(i)), 'Stride', 2, 'Cropping', 'same', 'Name', 'Upsample' + string(i));
        concat = depthConcatenationLayer(2, 'Name', 'Concat' + string(i));
        lgraph = addLayers(lgraph, upsample);
        lgraph = addLayers(lgraph, concat);  
        lgraph = connectLayers(lgraph, "profu_decoder" + string(i+1) + "_add1", "Upsample" + string(i));
        lgraph = connectLayers(lgraph, "profu" + string(i) + "_add1", "Concat" + string(i) + "/in2");
        lgraph = connectLayers(lgraph, "Upsample" + string(i), "Concat" + string(i) + "/in1");
        % Reducimos el numero de filtros de la salida
        convReduce = convolution2dLayer(1, (32*(2^i)), "Padding", "same", "Name", "reduceFilters"+string(i));
        lgraph = addLayers(lgraph, convReduce);
        lgraph = connectLayers(lgraph, "Concat" + string(i), "reduceFilters"+string(i));
        % Creamos el bloque RES y lo conectamos a la entrada y a la capa de
        % suma
        res = resBlock(32*(2^i), "profu_decoder" + string(i));
        lgraph = addLayers(lgraph, res);
        lgraph = connectLayers(lgraph, "reduceFilters"+string(i), "profu_decoder" + string(i) + "_bn1");
        lgraph = connectLayers(lgraph, "reduceFilters"+string(i), "profu_decoder" + string(i) + "_add1/in2");
        
    
    % % Concatenate with the corresponding ResBlock output from the encoder
    %     concat = depthConcatenationLayer(2, 'Name', 'Concat' + string(i));
    %     lgraph = addLayers(lgraph, concat);
    %     lgraph = connectLayers(lgraph, 'profu_decoder' + string(i+1) + "_add3", 'Upsample' + string(i));
    %     lgraph = connectLayers(lgraph, "Upsample" + string(i), "Concat" + string(i) + "/in1");
    %     lgraph = connectLayers(lgraph, "profu" + string(i) + "_add3", "Concat" + string(i) + "/in2");
    %     % Residual Block
    %     res = resBlock(32*(2^i), "profu_decoder" + string(i));
    %     lgraph = addLayers(lgraph, res);
    %     lgraph = connectLayers(lgraph, "Concat" + string(i), "profu_decoder" + string(i) + "_bn1");
    %     conv = convolution2dLayer(3, 32*(2^(i)), "Padding", "same", "Name", "ReduceFilters"+ string(i));
    %     lgraph = addLayers(lgraph, conv);
    %     lgraph = connectLayers(lgraph, "Concat" + string(i), "ReduceFilters"+ string(i));
    %     lgraph = connectLayers(lgraph, "ReduceFilters"+ string(i), "profu_decoder" + string(i) + "_add1/in2");
    %     lgraph = connectLayers(lgraph, 'profu_decoder' + string(i) + '_add1', "profu_decoder" + string(i) + "_add2/in2");
    %     lgraph = connectLayers(lgraph, 'profu_decoder' + string(i) + '_add2', "profu_decoder" + string(i) + "_add3/in2");
    end

    % Capa de salida
    % outputLayer = convolution2dLayer(1, 1, 'Padding', 'same', 'Name', 'output');
    % lgraph = addLayers(lgraph, outputLayer);
    % lgraph = connectLayers(lgraph, "profu_decoder" + string(1) + "_add3", 'output');
end

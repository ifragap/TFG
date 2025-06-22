function lgraph = encoder(inputSize,profundidad)
    
    % inputSize=[544,720,1];
    % profundidad=8;
    lgraph = layerGraph();
    
    % Entrada
    lgraph = addLayers(lgraph, imageInputLayer(inputSize, "Name", "input"));
    conv1_1 = convolution2dLayer(3,64,'Padding','same','Name','conv1_1');
    capa_batch=batchNormalizationLayer('Name',"bn1_1");
    lgraph = addLayers(lgraph, conv1_1);
    lgraph = addLayers(lgraph, capa_batch);
    lgraph = connectLayers(lgraph, 'input', 'conv1_1');
    lgraph = connectLayers(lgraph, 'conv1_1', 'bn1_1');
    
    % Primer bloque del res
    res=initial_resBlock(64,"profu1");
    lgraph = addLayers(lgraph, res);
    lgraph = connectLayers(lgraph, 'bn1_1', 'profu1_conv1');
    % Realizar dowsampling y batch normalization
    dow_sample=convolution2dLayer(1, 64*2, "Stride", 2, "Padding", "same", "Name", "Downsample1");
    % capa_batch=batchNormalizationLayer('Name',"bn1_2");
    lgraph = addLayers(lgraph, dow_sample);
    % lgraph = addLayers(lgraph, capa_batch);
    lgraph = connectLayers(lgraph, 'profu1_add1', 'Downsample1');
    % lgraph = connectLayers(lgraph, 'Downsample1', 'bn1_2');
    % Conectar el input con la capa de suma del bloque de RES
    lgraph = connectLayers(lgraph, 'bn1_1', 'profu1_add1/in2');
    
    for i=2:1:profundidad
        % Creamos el siguiente bloque RES y lo conectamos con la anterior 
        % capa de profundidad
        res=resBlock(32*(2^i),"profu"+string(i));
        lgraph=addLayers(lgraph,res);
        lgraph = connectLayers(lgraph, 'Downsample'+string(i-1), "profu"+string(i)+"_bn1");
        lgraph = connectLayers(lgraph, 'Downsample'+string(i-1), "profu"+string(i)+"_add1/in2");
        % Hacemos el dowsampling
        dow_sample = convolution2dLayer(1, 32*(2^(i+1)), "Stride", 2, "Padding", "same", "Name", "Downsample" + string(i));
        % capa_batch=batchNormalizationLayer('Name',"bn"+string(i)+"_2");
        lgraph = addLayers(lgraph, dow_sample);
        % lgraph = addLayers(lgraph, capa_batch);
        lgraph = connectLayers(lgraph, 'profu'+string(i)+'_add1', 'Downsample'+string(i));
        % lgraph = connectLayers(lgraph, 'Downsample'+string(i), "bn"+string(i)+"_2");
    end

end

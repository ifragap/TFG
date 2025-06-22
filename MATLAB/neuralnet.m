function red_lista=neuralnet(dimensiones,profundidad,numero_clases)
    
    bridge_1=bridge(dimensiones,profundidad);
    encoder_bridge_decoder=decoder(bridge_1,profundidad);
    outputLayer = convolution2dLayer(1, numero_clases, 'Padding', 'same', 'Name', 'output');
    net_full = addLayers(encoder_bridge_decoder, outputLayer);
    net_full = connectLayers(net_full, "profu_decoder" + string(1) + "_add1", 'output');
    capasoftmax=softmaxLayer("Name", "Softmax");
    net_full = addLayers(net_full, capasoftmax);
    net_full=connectLayers(net_full,"output","Softmax");
    red_lista = dlnetwork(net_full);
    
end
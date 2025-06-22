function layers = initial_resBlock(numFilters,name)


    layers=[
            convolution2dLayer(3, numFilters, "Padding", "same","DilationFactor",1,"Name", name + "_conv1")
            batchNormalizationLayer('Name', name + "_bn1")
            reluLayer('Name', name + "_relu1")
            convolution2dLayer(3, numFilters, "Padding", "same","DilationFactor",3,"Name", name + "_conv2")
            additionLayer(2, 'Name', name + "_add1")
    ];
    
    % layers = [
    %     batchNormalizationLayer('Name', name + "_bn1")
    %     reluLayer('Name', name + "_relu1")
    %     convolution2dLayer(3, numFilters, "Padding", "same","Name", name + "_conv1")
    %     additionLayer(2, 'Name', name + "_add1")
    %     %
    %     batchNormalizationLayer('Name', name + "_bn2")
    %     reluLayer('Name', name + "_relu2")
    %     convolution2dLayer(3, numFilters, "Padding", "same","Name", name + "_conv2")
    %     additionLayer(2, 'Name', name + "_add2")
    %     %
    %     batchNormalizationLayer('Name', name + "_bn3")
    %     reluLayer('Name', name + "_relu3")
    %     convolution2dLayer(3, numFilters, "Padding", "same","Name", name + "_conv3")
    %     additionLayer(2, 'Name', name + "_add3")
    % ];
end


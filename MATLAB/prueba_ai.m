% % % data_etiquetado = load('imagenes_etiquetadas.mat');
%    imageFolder = 'D:\IMAGENES_TRABAJO\Validacion\20_02_2025\imagenes';
% % % % % % 
%    outputFolder = 'D:\IMAGENES_TRABAJO\Validacion\20_02_2025\imagenes_redimensionadas';
% % % % % % 
%    imds = imageDatastore(imageFolder);
% % % % % % 
% % % % % imagen=imread("D:\IMAGENES_TRABAJO\ETIQUETADO ANTERIOR_SEMANA_FINAL_ENERO\redimensionado_etiquetado\Basler_acA720-520um__40083930__20241024_104751790_1301.tiff");
% % % % % resizedImg = imresize(imagen, [544 720]);
% % % % % imwrite(resizedImg,"C:\Users\Diego Ñaña Mejía\Documents\CTB-CHAMBA\CODIGO_MATLAB\codigo_imagenes_v1.6\ETIQUETAR\redimensionado\Basler_acA720-520um__40083930__20241115_120613835_14002.tiff")
%  while hasdata(imds)
%      [img, info] = read(imds);
%      %img=img(97:end,:);
% %     img=imrotate(img,90);
% %     % Redimensionar la imagen
%       resizedImg = img(14:525,104:615);
% %     resizedImg = img(:,104:615);
%       %resizedImg = imresize(img, [512 512],"nearest");
% %     Reescalar imagen
%       %resizedImg = im2double(img);  
% %     % Guardar la imagen redimensionada
%       [~, name, ext] = fileparts(info.Filename);
%       imwrite(resizedImg, fullfile(outputFolder, [name, ext]));
% % 
%  end
% % 
%    outputLabelFolder = 'D:\IMAGENES_TRABAJO\Validacion\20_02_2025\etiquetas_redimensionadas';
% % % 
%    direccion_etiqueta=gTruthMed.LabelData;
% % % % 
% for i=1:21
%     variable=load(direccion_etiqueta(i,1));
%     imagen=variable.labels;
%     imagen = imagen(14:525,104:615);
%     %resizedImg = imresize(imagen, [544 720],"nearest");
%     [~, name, ext] = fileparts(direccion_etiqueta(i,1));
%     imwrite(imagen, outputLabelFolder+"\"+name+".bmp");
% end
% for i=1:1:52
%     variable=load(gTruthMed.LabelData(i));
%     imagen=variable.labels;
%     [~, name, ext] = fileparts(gTruthMed.LabelData(i));
%     imwrite(imagen, "D:\IMAGENES_TRABAJO\ETIQUETADO\etiquetado_mas"+"\"+name+".tiff");
%     % x=gTruthMed.LabelData(i);
%     % y=load(gTruthMed.LabelData(i));
%     % y=y.labels;
%     % save("D:\IMAGENES_TRABAJO\ETIQUETADO\etiquetado_mas\")
% end
% Configuración Inicial
% imageFolder = 'D:\IMAGENES_TRABAJO\total_imagenes\imagen';
% labelFolder = 'D:\IMAGENES_TRABAJO\total_imagenes\etiqueta';
% 
% 
% 
% 
% 
% %[encoder,~] = imagePretrainedNetwork("efficientnetb0");
% 
% % Crear un Datastore para las Imágenes y las Etiquetas
% imds = imageDatastore(imageFolder);
% 
% pxds = pixelLabelDatastore(labelFolder, ["fondo","canal","celula"], [0,1,2]);
% % pxds = pixelLabelDatastore(labelFolder, ["canal","celula"], [1,2]);
% 
% % Configuración de la U-Net
% inputSize = [512 512 1]; % Tamaño de entrada de las imágenes
% numClasses = 3; % Número de clases de segmentación
% profundidad=4;
% %lgraph = unet(inputSize, numClasses);
% %resunet_propia = unet([544 720 1],4);
% resunet_propia=neuralnet(inputSize,profundidad,numClasses);
% % Opciones de Entrenamiento
% options = trainingOptions("adam", ...
%     'InitialLearnRate', 1e-4, ...
%     'MaxEpochs',25, ...
%     'MiniBatchSize', 3, ...
%     'Shuffle', 'every-epoch', ...
%     'Plots', 'training-progress', ...
%     'Metrics','accuracy', ...
%     'Verbose', false, ...
%     'ExecutionEnvironment', 'gpu'); % Usar la GPU para el entrenamiento
% % Preparar los Datos de Entrenamiento
% ds = combine(imds, pxds);
% %ds = pixelLabelImageDatastore(imds,pxds);
% 
% 
% 
% 
% 
% 
% net_n=trainnet(ds,resunet_propia,"crossentropy",options);
% respaldo=net_n;
% %seguir entrenando 
% net_2=trainnet(ds,net_n,"crossentropy",options);
% net_3=trainnet(ds,net_2,"crossentropy",options);
% net_4=trainnet(ds,net_3,"crossentropy",options);
% net_5=trainnet(ds,net_4,"crossentropy",options);
% % dataOriginal = read(ds);
% % dataAugmented = read(dsAugmented);
% 
% % figure;
% % subplot(2,2,1); imshow(dataOriginal{1}); title('Imagen Original');
% % subplot(2,2,2); imshow(uint8(dataOriginal{2}),[0 3]); title('Segmentación Original');
% % subplot(2,2,3); imshow(dataAugmented{1}); title('Imagen Aumentada');
% % subplot(2,2,4); imshow(uint8(dataAugmented{2}),[0 3]); title('Segmentación Aumentada');
% 
% 
% 
% 
% 
% 

% Entrenar la Red
 
% Guardar el Modelo Entrenado
% save('resunet_propia_v7.mat', 'net_n');

load('resunet_propia_v7.mat')
probando_imagen=imread("D:\IMAGENES_TRABAJO\27_03_2025\1\test3\Basler_acA720-520um__40083930__20250327_113128502_3744.tiff");
% imagen_francia=imrotate(probando_imagen,90);
% imagen_francia=imagen_francia(:,104:615);
% resizedImg=zeros(512,512);
% resizedImg=im2uint8(resizedImg);
% resizedImg(:)=112;
% resizedImg(45:468,1:512)=imagen_francia;
%resizedImg = imresize(probando_imagen, [544 720]);

% probando_imagen=imrotate(probando_imagen,90);
% resizedImg = probando_imagen(:,104:615);
% %Estiramos la imagen 20.75% en el eje y 
% resizedImg = imresize(resizedImg, [512 512],"bicubic");

resizedImg = probando_imagen(14:525,104:615);
%resizedImg=probando_imagen;
imshow(resizedImg)
C = semanticseg(resizedImg, net_n);

% Mostrar la segmentación superpuesta a la imagen original
B = labeloverlay(resizedImg, C);
figure
imshow(B);
celula=C=="C3";

figure
imshow(celula)

[regiones,n]= bwlabel(celula,4);

propiedades=regionprops(regiones==1,"FilledArea");
areas=zeros(n,1);
areas(1)=propiedades.FilledArea;

for i=2:1:n
        region_n=regiones==i;
        propiedades=regionprops(region_n,"FilledArea");
        areas(i)=propiedades.FilledArea;
        %propiedades.Orientation
end

[valor,pos]=max(areas);

imshow(regiones==pos)



% Mostrar la segmentación en un color diferente
figure;
imshow(C, 'InitialMagnification', 'fit');
colormap('jet');
colorbar;
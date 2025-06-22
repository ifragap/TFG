function [pos,region,regiones] = etiquetado(imagen)
    
    % Vamos a binarizar la imagen para que el crecimiento de regiones sea
    % mas sencillo y no tengamos que meter un rango.
    
    bordes= edge(imagen,'canny');
    bordes=imclose(bordes,ones(5));
    bordes(1:end,1)=1;
    bordes(1:end,end)=1;
    figure, imshow(imcomplement(bordes))
    
    % Se introduce la semilla de manera interactiva
    
    [x,y]=ginput(1);

    % Etiquete la máscara
    [regiones,n]= bwlabel(imcomplement(bordes),4);
    
    
    % Cree una nueva máscara binaria de aquellos píxeles con la misma
    % etiqueta que el punto seleccionado
    v=regiones(round(y),round(x));
    region=regiones==v;
    
    coordenadas1=regionprops(region,"Orientation","FilledArea");
    propiedades=regionprops(regiones==1,"Orientation","FilledArea","Extrema");

    for i=2:1:n
        region_n=regiones==i;
        propiedades=[propiedades regionprops(region_n,"Orientation","FilledArea","Extrema")];
    end

    similitud=zeros(n,2);

    for i=1:1:n
        ratio=propiedades(1,i).FilledArea/coordenadas1.FilledArea;
        similitud(i,:)=[abs(ratio-1),i];

    end

    similitud=sortrows(similitud);

    pos=similitud;
    region=propiedades;


    % Represente la imagen original superponiendo la región calculada en
    % rojo. 
    % imagen=double(imagen);
    % imagen = (imagen - min(imagen(:)))/(max(imagen(:))-min(imagen(:))); 
    % r = imagen; 
    % r(region)=1;
    % g = imagen;
    % g(region)=0;
    % b = imagen;
    % b(region)=0;
    % 
    % 
    % 
    % rgb(:,:,1) = r;
    % rgb(:,:,2) = g;
    % rgb(:,:,3) = b;
    % 
    % 
    % imshow(rgb);
    % 
end
% function [limite_izquierdo,limite_derecho,limite_superior,limite_inferior,n_region,region,region2]= get_mascara(imagen)
function [posicion_tub]= find12channel(imagen_sucia,imagen_fondo)
    
    % Quitamos el fondo de la imagen 
    imagen_limpia=filtro1(imagen_sucia,imagen_fondo);
    
    % Conseguimos las regiones y posiciones de la imagen
    [pos,regi,~]=etiquetado(imagen_limpia);
    close all

    % Obtenemos solo las regiones de las paredes del canal
    regiones_features=regi(1,pos(1,2));

    for i=2:1:11
        regiones_features=[regiones_features regi(1,pos(i,2))];
    end
    
    % Ordenamos las posiciones de las paredes del canal
    posicion_y_top=zeros(11,3);
    posicion_y_top(1,:)=[pos(1,2) regiones_features(1,1).Extrema(8,2)  1];

    for i=2:1:11
    posicion_y_top(i,:)=[pos(i,2) regiones_features(1,i).Extrema(8,2) i];
    end

    posicion_y_top=sortrows(posicion_y_top,2);

    % Obtenemos las posiciones de los canales en si
    posicion_tub=zeros(12,2);
    posicion_tub(1,:)=regiones_features(1,posicion_y_top(1,3)).Extrema(8,:);
    posicion_tub(1,2)=posicion_tub(1,2)-23;
    for i=1:1:11
        posicion_tub(i+1,:)=regiones_features(1,posicion_y_top(i,3)).Extrema(7,:);
    end

    posicion_tub=int32(round(posicion_tub));
    
end
% function [limite_izquierdo,limite_derecho,limite_superior,limite_inferior,n_region,region,region2]= get_mascara(imagen)
function [posicion_tub]= find3channel(imagen_sucia,imagen_fondo)
    
    % Quitamos el fondo de la imagen 
    imagen_limpia=filtro1(imagen_sucia,imagen_fondo);
    
    % Conseguimos las regiones y posiciones de la imagen
    [pos,regi,~]=etiquetado(imagen_limpia);
    close all

    % Obtenemos solo las regiones de las paredes del canal
    regiones_features=regi(1,pos(1,2));
    regiones_features=[regiones_features regi(1,pos(2,2))];
   
    
    % Ordenamos las posiciones de las paredes del canal
    posicion_y_top=zeros(2,3);

    posicion_y_top(1,:)=[pos(1,2) regiones_features(1,1).Extrema(1,2)  1];
    posicion_y_top(2,:)=[pos(2,2) regiones_features(1,2).Extrema(1,2) 2];

    posicion_y_top=sortrows(posicion_y_top,2);

    % Obtenemos las posiciones de los canales en si
    posicion_tub=zeros(3,2);
    posicion_tub(1,:)=regiones_features(1,posicion_y_top(1,3)).Extrema(1,:);
    posicion_tub(1,2)=posicion_tub(1,2)-40;
    for i=1:1:2
        posicion_tub(i+1,:)=regiones_features(1,posicion_y_top(i,3)).Extrema(6,:);
    end

    posicion_tub=int32(round(posicion_tub));
    
    

end


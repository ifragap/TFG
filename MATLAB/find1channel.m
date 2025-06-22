function [posicion_tub,posicion_pend,posicion_cons,ancho_canal,ancho_cons,largo] = find1channel(imagen_sucia,imagen_fondo,hay_fondo)
    
    % Quitamos el fondo de la imagen 
    if(hay_fondo == true)
        imagen_limpia=filtro1(imagen_sucia,imagen_fondo);
        [posicion_tub,posicion_pend,posicion_cons,ancho_canal,ancho_cons,largo] = new_etiquetado(imagen_limpia);
        close all
    else
        [posicion_tub,posicion_pend,posicion_cons,ancho_canal,ancho_cons,largo] = new_etiquetado(imagen_sucia);
        close all
    end

end


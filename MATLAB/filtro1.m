function [filtrado] = filtro1(imagen_sucia,imagen_fondo)
    
    imagen_fondo=im2double(imagen_fondo);
    imagen_sucia=im2double(imagen_sucia);
    filtrado=imcomplement(-log(imagen_sucia./imagen_fondo));
    %filtrado=im2uint16(filtrado);ç

end


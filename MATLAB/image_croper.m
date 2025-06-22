function [coordenada_izquierda_superior,coordenada_derecha_superior,altura] = image_croper(imagen,izquierdo_superior,largo,ancho)
        
    [dimension_y,dimension_x]=size(imagen);
    espacio_izquierda=izquierdo_superior(1,1);
    espacio_derecha=dimension_x-(izquierdo_superior(1,1)+largo);
    espacio_superior=izquierdo_superior(1,2);
    espacio_inferior=dimension_y-(izquierdo_superior(1,2)+ancho);
    
    espacio_disponible_x=dimension_x-largo;
    espacio_necesario_x=512-largo;
    espacio_disponible_y=dimension_y-ancho;
    espacio_necesario_y=512-ancho;

    resize_x=false;
    resize_y=false;

    if espacio_necesario_x>espacio_disponible_x
        resize_x=true;
    end

    if espacio_necesario_y>espacio_disponible_y
        resize_y=true;
    end
    
    modulo_div=mod(espacio_necesario_x,2);
    
    if modulo_div==0
        corte_izquierdo_x=ceil(espacio_necesario_x/2.0);
        corte_derecha_x=ceil(espacio_necesario_x/2.0);
    else
        corte_izquierdo_x=ceil((espacio_necesario_x-1)/2.0);
        corte_derecha_x=ceil((espacio_necesario_x-1)/2.0);
        corte_izquierdo_x=corte_izquierdo_x+1;
    end

    
    if espacio_izquierda<corte_izquierdo_x
        corte_derecha_x=(corte_izquierdo_x-espacio_izquierda)+corte_derecha_x;
        corte_izquierdo_x=espacio_izquierda;
    end

    if espacio_derecha<corte_derecha_x

        corte_izquierdo_x=(corte_derecha_x-espacio_derecha)+corte_izquierdo_x;
        corte_derecha_x=espacio_derecha;
    end
    
    modulo_div=mod(espacio_necesario_y,2);

    if modulo_div==0
        corte_superior_y=ceil(espacio_necesario_y/2.0);
        corte_inferior_y=ceil(espacio_necesario_y/2.0);
    else
        corte_superior_y=ceil((espacio_necesario_y-1)/2.0);
        corte_inferior_y=ceil((espacio_necesario_y-1)/2.0);
        corte_superior_y=corte_superior_y+1;
    end

    if espacio_superior<corte_superior_y
        corte_inferior_y=(corte_superior_y-espacio_superior)+corte_inferior_y;
        corte_superior_y=espacio_superior;
    end

    if espacio_inferior<corte_inferior_y

        corte_superior_y=(corte_inferior_y-espacio_inferior)+corte_superior_y;
        corte_inferior_y=espacio_inferior;
    end

    coordenada_izquierda_superior=[1+(izquierdo_superior(1,1)-corte_izquierdo_x),1+(izquierdo_superior(1,2)-corte_superior_y)];
    coordenada_derecha_superior=[1+(izquierdo_superior(1,1)+largo+corte_derecha_x),1+(izquierdo_superior(1,2)-corte_superior_y)];
    
    altura=coordenada_izquierda_superior(1,2)-coordenada_izquierda_superior(1,2)+corte_superior_y+ancho+corte_inferior_y;

end
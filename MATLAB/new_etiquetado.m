function [posicion_tub,posicion_pend,posicion_cons,ancho_canal,ancho_cons,largo] = new_etiquetado(imagen_limpia)
 
    % Vamos a binarizar la imagen para que el crecimiento de regiones sea
    % mas sencillo y no tengamos que meter un rango.
    [bordes1]=edge(imagen_limpia,'prewitt',"horizontal");
    [bordes2]=edge(imagen_limpia,'Roberts',"horizontal");
    [bordes3]=edge(imagen_limpia,'Roberts',"vertical");

    bordes=bordes1+bordes2+bordes3;
    bordes=bordes>0;
    clear bordes1 bordes2 bordes3

    % Volvemos los bordes de la imagen a 1 para que no haya errores mas
    % adelante cuando se binariza
    bordes(1:end,1)=1;
    bordes(1:end,end)=1;
    figure, imshow(imcomplement(bordes))
    
    % Introduzca la semilla de manera interactiva usando el comando ginput
    [x,y]=ginput(1);

    % Etiquete la máscara
    [regiones,~]= bwlabel(imcomplement(bordes),4);
    
    % Cree una nueva máscara binaria de aquellos píxeles con la misma
    % etiqueta que el punto seleccionado
    v=regiones(round(y),round(x));

    % Se rellena los agujeros de la region binarizada
    region=regiones==v;
    region=imfill(region,"holes");
  
    % Obtenemos los perfiles a lo largo del eje x de la imagen, vamos de la
    % posicion 2 a la 719 para no tomar los bordes a 1
    tamano=size(imagen_limpia);
    x=zeros(1,tamano(1,2));
    for i=2:tamano(1,2)-1
        for j=1:1:tamano(1,1)
            valor=region(j,i);
            if (valor==1)
                x(i)=j;
                break
            end
        end
    end

    % Ya teniendo el perfil rellenamos los valores que pusimos a drede en
    % los bordes para que estos se correspondan con el valor de Y de 
    % nuestra region correspondientes
    x(1)=x(2);
    x(end)=x(end-1);

    % Suavizamos el perfil 
    m=10;
    ventana=ones(1,m)/m;
    x_movil=filtfilt(ventana,1,x);
    derivada=diff(x_movil);
    flip_derivada=flip(derivada);
    
    % A traves de la derivada obtenemos los puntos que nos interesa los
    % cuales son donde debemos cortar la imagen (punto_3), donde tenemos la
    % pendiente de nuestra constriccion (punto_1) y donde esta nuestra
    % constriccion(punto_2)
    [~,punto_1]=max(derivada);

    % El punto_1_next es el punto donde se tiene la pendiente de la
    % constriccion pero desde es el otro extremo de la imagen, lo mismo que
    % punto_2_next y punto_3_next
    [~,punto_1_next]=min(flip_derivada);
    
    for i=punto_1:1:length(derivada)
        if(derivada(i)<=0)
            punto_2=i;
            break
        end
    end
    for i=punto_1:-1:1
        if(derivada(i)<=0)
            punto_3=i;
            break
        end
    end
    
    for i=punto_1_next:-1:1
        if(flip_derivada(i)>=0 || i==1)
            punto_3_next=i;
            break
        end
    end

    punto_3_next=size(region,2)-punto_3_next;
    largo=punto_3_next-punto_3;

    % Obtenemos perfiles a lo largo del eje y correspondientes al punto 1 y
    % 3 de la propia imagen limpia para hacer correciones en la posición. 
    % El punto 2 tambien se hace lo mismo pero a traves de un eje a 45 
    % grados
    perfil_posicion_tub=flip(imagen_limpia(round(x_movil(punto_3))-10:round(x_movil(punto_3)),punto_3));
    perfil_posicion_constriccion=flip(imagen_limpia(round(x_movil(punto_2))-10:round(x_movil(punto_2)),punto_2));
    perfil_posicion_pendiente=zeros(11,1);
    perfil_posicion_pendiente(1,1)=imagen_limpia(round(x_movil(punto_1)),punto_1);
    for i=2:1:11
        perfil_posicion_pendiente(i,1)=imagen_limpia(round(x_movil(punto_1))-i+1,punto_1+i-1);
    end
   
    % Corregimos las posiciones y obtenemos los puntos anteriormente
    % mencionados
    [~,correccion_posicion_tub]=min(perfil_posicion_tub);
    [~,correccion_posicion_constriccion]=min(perfil_posicion_constriccion);
    [~,correccion_posicion_pendiente]=min(perfil_posicion_pendiente);
    correccion_posicion_tub=correccion_posicion_tub-1;
    correccion_posicion_constriccion=correccion_posicion_constriccion-1;
    correccion_posicion_pendiente=correccion_posicion_pendiente-1;

    posicion_tub=[punto_3,round(x_movil(punto_3))-(correccion_posicion_tub)];
    posicion_pend=[punto_1+correccion_posicion_pendiente,round(x_movil(punto_1))-correccion_posicion_pendiente];
    posicion_cons=[punto_2,round(x_movil(punto_2))-(correccion_posicion_constriccion)];
    
    % Conseguir ancho en pixeles del canal y constriccion
    tamano=size(imagen_limpia);
    
    ancho_canal = 0;
    contador = 0;
    for i=posicion_tub(1,2)+10:1:tamano(1,1)-1
         valor=region(i,posicion_tub(1,1));
         contador=contador+1;
         if(valor==0)
             ancho_canal=i-posicion_tub(1,2);
             break
         end    
    end
    
    ancho_cons = 0;
    contador = 0;
    for i=posicion_cons(1,2)+5:1:tamano(1,1)-1
        valor=region(i,posicion_cons(1,1));
        contador=contador+1;
        if(valor==0)
            ancho_cons=i-posicion_cons(1,2)-1;
            break
        end    
    end

    figure
    imshowpair(imagen_limpia, region, 'blend')
    hold on;

    plot(posicion_tub(1),posicion_tub(2),'ro','MarkerSize',8);
    plot(posicion_pend(1),posicion_pend(2),'go','MarkerSize',8);
    plot(posicion_cons(1),posicion_cons(2),'ko','MarkerSize',8);

    legend('Zona ancha', 'Pendiente', 'Constriccion');

    saveas(gcf, 'canal.png')
end
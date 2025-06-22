function [direcciones_imagenes,imagen_fondo,imagen_base] = search_adress_images(volteado,hay_fondo)
    
    direc_experimentos = uigetdir(path,"Selecciona la carpeta del experimento");
    imagen_fondo=[];
    if(hay_fondo==true)
        direc_fondos = uigetdir(path,"Selecciona la carpeta de los fondos");
    end
    direc_base=uigetdir(path,"Seleccione la carpeta base de los test");
    filePattern = fullfile(direc_experimentos, '*.tiff');
    imagenes=dir(filePattern);
    y=size(imagenes);
    if (y(1,1)==0)
        x=input("Que formato de imagen utilizas \n","s");
        filePattern = fullfile(direc_experimentos, "*."+x);
        imagenes=dir(filePattern);
        if(hay_fondo==true)
            imagen_fondo=mean_images(direc_fondos,x);
        end
        imagen_base=mean_images(direc_base,x);
    else
        if(hay_fondo==true)
            imagen_fondo=mean_images(direc_fondos,"tiff");
        end
    imagen_base=mean_images(direc_base,"tiff");
    end
    
    if(volteado==true)
        imagen_base = imrotate(imagen_base,90);
        if(hay_fondo==true)
            imagen_fondo= imrotate(imagen_fondo,90);
        end
    end

    imagenes=natsortfiles(imagenes);
    %direccion_fondos=dir(filePattern_2);
    %direccion_fondos=direc_fondos + "\" + direccion_fondos(1,1).name;

    direcciones_imagenes(1:length(imagenes))="";

    
    for i=1:1:length(imagenes)
        direcciones_imagenes(1,i) = direc_experimentos + "\" + imagenes(i,1).name;
    end
    
end


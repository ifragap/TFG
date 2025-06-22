classdef Dispositivo < handle

    properties
        id                  % Integer - ID of the device used in the test
        imagen_base         % Image - Base image
        imagen_base_ori     % Imagen base sin modificaciones
        imagen_fondo        % Image - Background image of the microscope
        imagen_fondo_ori    % Imagen fondo sin modificaciones
        hay_fondo           % Boolean - Indicates if a bakcground image was
                            % passed as input
        resize_x            % Boolean that indicates if a resize in x axis 
                            % was done
        resize_y            % Boolean that indicates if a resize in y axis 
                            % was done                    
        volteado            % Boolean - Indicates if the images are rotated
        direc_imagenes      % Array de String - Array con todas las 
                            % direcciones de las imagenes en memoria
        num_de_imagenes     % Integer - Numero de imagenes en el test
        tipo_dispositivo    % Integer - Tipo de dispositivo (1, 3 o 12) 
        tamano_canal        % Double - Ancho de la constriccion de los 
                            % canales (5, 7.5 o 10) 
        fs                  % Integer - Frecuencia de muestreo de imagen
        Q                   % Double - Flow of the experiment
        long                % Double - Ancho para el corte del canal
        height              % Double - Altura para el corte del canal
        t_start             % Double - Epoch Time of the start of the test
        date_start          % String - Date and time of the start of the 
                            % test
        Canales             % Array de objetos Canal - Array de canales con
                            % la ubicacion y los eventos de estos
        ElectricData        % ElectricData - Object with all electric data
    end


    methods

        function obj=Dispositivo(id,tipo_dispositivo,tamano_canal,fs,Q,hay_fondo,volteado)
            
            % In case there is no ID as input, size and number channels
            % uses the other inputs
            if(id == 0)
                obj.id = "Device previous to inventory system";
                obj.tipo_dispositivo=tipo_dispositivo;
                obj.tamano_canal=tamano_canal;
            else
                obj.id = id;
                
                var = num2str(obj.id);

                switch var(1:2)
                    case '05'
                        obj.tamano_canal = 5;
                    case '07'
                        obj.tamano_canal = 7.5;
                    case '10'
                        obj.tamano_canal = 10;
                end

                switch var(3:4)
                    case '01'
                        obj.tipo_dispositivo = 1;
                    case '03'
                        obj.tipo_dispositivo = 3;
                    case '12'
                        obj.tipo_dispositivo = 12;
                end
            end
       

            %NUEVA PARTE 24/04/2025
            % [obj.direc_imagenes,obj.imagen_fondo,obj.imagen_base]=search_adress_images(volteado,hay_fondo);
            [obj.direc_imagenes,imagen_fondo,imagen_base]=search_adress_images(volteado,hay_fondo);
            if hay_fondo==false
                imagen_fondo=imagen_base;
                imagen_fondo(:)=0;
            end



            dimension=size(imagen_base); %[tamaño en y, tamaño en x]
            obj.imagen_base_ori=imagen_base;
            obj.imagen_fondo_ori=imagen_fondo;
            obj.resize_y=false;
            obj.resize_x=false;
            imagen_base_r=imagen_base;
            imagen_fondo_r=imagen_fondo;
            if(dimension(1,1)<512)
                imagen_base_r = imresize(imagen_base, [512 dimension(1,2)],"bicubic");
                imagen_fondo_r = imresize(imagen_fondo, [512 dimension(1,2)],"bicubic");
                resize_x=true;
                obj.resize_x=resize_x;
            end

            if(dimension(1,2)<512)
                imagen_base_r = imresize(imagen_base, [dimension(1,1) 512],"bicubic");
                imagen_fondo_r = imresize(imagen_fondo, [dimension(1,1) 512],"bicubic");
                resize_y=true;
                obj.resize_y=resize_y;
            end

            obj.volteado=volteado;
            obj.hay_fondo=hay_fondo;
            
            obj.num_de_imagenes=length(obj.direc_imagenes);
            
            %obj.setLong_Heigth(tipo_dispositivo,tamano_canal);
            obj.fs=fs;
            obj.Q = Q;
            
            %Añadi esto para los metadatos de francia

            origen_imagen= input("Are these images from the UPM, 1=yes,0=no");
            if (origen_imagen==1)
                obj.t_start = epoch_time(obj.direc_imagenes(1));
                obj.date_start = string(datetime(obj.t_start/1000,'ConvertFrom','posixtime'));
            end
            if(origen_imagen==0)
               obj.t_start = epoch_time_france_images(obj.direc_imagenes(1));
               obj.date_start = string(datetime(obj.t_start/1000,'ConvertFrom','posixtime'));
            end

            switch tipo_dispositivo
                case 1
                    % [posicion_tub,posicion_pend,posicion_cons,ancho_canal,ancho_cons,largo] = find1channel(obj.imagen_base,obj.imagen_fondo,hay_fondo);
                    [posicion_tub,posicion_pend,posicion_cons,ancho_canal,ancho_cons,largo] = find1channel(imagen_base_r,imagen_fondo_r,hay_fondo);
                    %añadir linea de codigo aqui que calcule el corte 
                    [position_izquierda_superior,coordenada_derecha_superior,altura]=image_croper(imagen_base_r,posicion_tub,largo,ancho_canal+3);

                    obj.imagen_base=imagen_base_r(position_izquierda_superior(1,2):position_izquierda_superior(1,2)+altura-1,position_izquierda_superior(1,1):coordenada_derecha_superior(1,1)-1);
                    obj.imagen_fondo=imagen_fondo_r(position_izquierda_superior(1,2):position_izquierda_superior(1,2)+altura-1,position_izquierda_superior(1,1):coordenada_derecha_superior(1,1)-1);
                    posicion_tub=posicion_tub-position_izquierda_superior;
                    posicion_pend=posicion_pend-position_izquierda_superior;
                    posicion_cons=posicion_cons-position_izquierda_superior;
                    obj.height = ancho_canal + 3;
                    obj.long = largo;
                    obj.Canales = Canal(position_izquierda_superior,posicion_tub,posicion_pend,posicion_cons,ancho_cons,obj);
                    
                %hay que modificar el case 3 y case 12, ahora solo funciona case 1   
                case 3
                    x=find3channel(obj.imagen_base,obj.imagen_fondo);
                    canales=Canal(x(1,:),obj);
                    for i=1:1:2
                        canales=[canales Canal(x(i+1,:),obj)];
                    end
                    obj.Canales=canales;
                case 12
                    x=find12channel(obj.imagen_base,obj.imagen_fondo);
                    canales=Canal(x(1,:),obj);    
                    for i=1:1:11
                        canales=[canales Canal(x(i+1,:),obj)];
                    end
                    obj.Canales=canales;
            end
            
            % Obtaining both the directory where all CSVs are stored and
            % the TXT file with t start for the test and/or date
            electric_path = uigetdir(path,['Select the directory of ' ...
                            'the electric data...']);
            [file_txt, path_txt] = uigetfile('*.txt', ['Select the TXT' ...
                                   ' file with the t_start...']);
            
            
            if(isequal(electric_path,0) | isequal(file_txt,0))
                warning("No electric data available.")
                obj.ElectricData = [];

            else
                txt_file = fullfile(path_txt, file_txt);
                csv_path = process_csv(electric_path);
                obj.ElectricData = ElectricData(csv_path,txt_file,obj);
            end
        end
        

        function find_events(obj)
            tic
            for i=1:1:length(obj.Canales)
                obj.Canales(i).event_finder();
                if obj.Canales(i).num_frame_eventos~=0
                    obj.Canales(i).create_events();
                end
            end
            toc
            beep
        end


        function analyse_events(obj)
            for i=1:1:length(obj.Canales)
                tic
                for j=1:1:length(obj.Canales.Eventos)
                    obj.Canales(i).Eventos(1,j).is_cell_and_radious()
                end
                beep

                if(~isempty(obj.ElectricData))
                    obj.ElectricData.set_C1_R1(0);
                    obj.ElectricData.asign_electric_data();
                end
                beep
                
                for j=1:1:length(obj.Canales.Eventos)
                    if(obj.Canales(i).Eventos(1,j).is_cell == 1)

                        fprintf("\nAnalysing: Evento %s\n", num2str(j));
                        obj.Canales(i).Eventos(1,j).cell_front_back();
                        obj.Canales(i).Eventos(1,j).estimating_Rc;
                        obj.Canales(i).Eventos(1,j).calc_final();
                        
                        if(~isempty(obj.ElectricData))
                            if(~isempty(obj.Canales(i).Eventos(1,j).data_electric()))
                                obj.Canales(i).Eventos(1,j).set_C2_R2(0);
                                obj.Canales(i).Eventos(1,j).centering_electric_data();
                            end
                        end

                        % obj.Canales.Eventos(1,j).final_analysis();
                    end
                end

                toc
                beep
            end
        end


        %Funcion que ya no se usa por el momento
        function  setLong_Heigth(obj,tipo,tamano)
            switch tipo
                case 1
                    switch tamano
                        case 5
                            obj.long=205;
                            %obj.height=102;
                        case 7.5
                            obj.long=205;
                            %obj.height=105;
                        case 10
                            obj.long=205;
                            %obj.height=113;
                    end
                case 3
                    switch tamano
                        case 5
                            obj.long=205;
                            obj.height=110;
                        case 7.5
                            obj.long=205;
                            obj.height=110;
                        case 10
                            obj.long=205;
                            obj.height=110;
                    end
                case 12
                    switch tamano
                        case 5
                            obj.long=205;
                            obj.height=110;
                        case 7.5
                            obj.long=205;
                            obj.height=110;
                        case 10
                            obj.long=205;
                            obj.height=110;
                    end                   
            end

        end


        %FUNCION PENSADA PARA UN SOLO CANAL, FALTA CAMBIARLO PARA QUE
        %ACEPTE MAS CANALES
        function writeEventsTxt(obj)
            nombre=input("Que nombre quieres para el archivo \n","s");
            direc= uigetdir(path,"Selecciona la carpeta donde desea guardar el archivo");
            %frame_eventos=zeros(length(obj.Canales.Eventos),2);
            Frame_inicial=zeros(length(obj.Canales.Eventos),1);
            Frame_final=zeros(length(obj.Canales.Eventos),1);
            Numero_evento=zeros(length(obj.Canales.Eventos),1);
            Duracion=zeros(length(obj.Canales.Eventos),1);
            for i=1:1:length(obj.Canales)
                for j=1:1:length(obj.Canales.Eventos)
                    Frame_inicial(j,1)=obj.Canales(i).Eventos(j).frame_inicial;
                    Frame_final(j,1)=obj.Canales(i).Eventos(j).frame_final;
                    Duracion(j,1)=obj.Canales(i).Eventos(j).frame_final-obj.Canales(i).Eventos(j).frame_inicial + 1;
                    Numero_evento(j,1)=j;
                end
            end
            tabla=table(Numero_evento,Frame_inicial,Frame_final,Duracion);
            writetable(tabla,direc+"\"+nombre)
        end
        
        function writeRealEventsTxt(obj)
            nombre=input("Que nombre quieres para el archivo \n","s");
            direc= uigetdir(path,"Selecciona la carpeta donde desea guardar el archivo");
            eventos=obj.Canales(1).Eventos;
            i=1;
            es_celula=0;
            Numero_evento=0;
            numero_eventos_totales=length(eventos);

            while(true)
                es_celula=eventos(i).is_cell;

                if (es_celula==1)
                    evento_real=eventos(i);
                    Numero_evento=i;
                    break
                end
                i=i+1;
                 if(i>numero_eventos_totales)
                     break
                 end
            end
            if Numero_evento==0
                disp("No hay eventos reales")
                return
            end
            es_celula=0;
            i=i+1;
            for j=i:1:length(eventos)
                es_celula=eventos(j).is_cell;
                if (es_celula==1)
                    evento_real=[evento_real eventos(j)];
                    Numero_evento=[Numero_evento; j];
                end

            end








            %frame_eventos=zeros(length(obj.Canales.Eventos),2);
            Frame_inicial=zeros(length(evento_real),1);
            Frame_final=zeros(length(evento_real),1);
            %Numero_evento=zeros(length(evento_real),1);
            Duracion=zeros(length(evento_real),1);
            %for i=1:1:length(obj.Canales)
            for j=1:1:length(evento_real)
                Frame_inicial(j,1)=evento_real(j).frame_inicial;
                Frame_final(j,1)=evento_real(j).frame_final;
                Duracion(j,1)=evento_real(j).frame_final-evento_real(j).frame_inicial + 1;
                %Numero_evento(j,1)=j;
            end
            %end
            tabla=table(Numero_evento,Frame_inicial,Frame_final,Duracion);
            writetable(tabla,direc+"\"+nombre)
        end









    end

end
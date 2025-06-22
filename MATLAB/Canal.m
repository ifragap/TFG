% CLASS "Canal" - Defined to store all events with cells, channel's width,
% coordinates, among others.

classdef Canal < handle

    properties 
        position_izquierda_superior % Array (2 elements) Double - Position
                                    % [x,y] to cut images to 512x512
        posicion_dispositivo        % Array (2 elements) Double - Position
                                    % [x,y] of start of the slope
        posicion_pendiente          % Array (2 elements) Double - Position
                                    % [x,y] of the mid-slope of the channel
        posicion_constriccion       % Array (2 elements) Double - Position 
                                    % [x,y] of the constriction
        bloqueo                     % Array Logical - Stores in which of 
                                    % all frames the channel is blocked
        linea                       % Double - Vertical line to obtain AL
        wch                         % Double - Semi-width of the channel
        alpha                       % Double - Angle formed by parallel and 
                                    % tangent lines
        rho                         % Double - R of outer circumference
        Dispositivo                 % Dispositivo - References the device 
                                    % from which the channel comes from
        num_frame_eventos           % Double - Number of frames in which 
                                    % something is passing
        Eventos                     % Array Evento - Storage of Eventos

    end


    methods
        
        % Constructor of the class
        function obj = Canal(position_izquierda_superior,posicion_dispositivo, posicion_pendiente,posicion_constriccion,ancho_cons,Dispositivo)
            
            obj.position_izquierda_superior=position_izquierda_superior;
            obj.posicion_dispositivo = posicion_dispositivo;
            obj.posicion_pendiente = posicion_pendiente;
            obj.posicion_constriccion = posicion_constriccion;
            obj.wch = ancho_cons/2.0;
            obj.bloqueo = false(1,Dispositivo.num_de_imagenes);            
            obj.Dispositivo = Dispositivo;
            set_linea_alpha(obj);
            set_rho(obj);

        end

        % Function 
        function event_finder(obj)

            disp=obj.Dispositivo; %utilizamos el objeto padre referenciandolo(Apunta a la misma direccion de memoria)
            imagen_base=disp.imagen_base; %La imagen base ya esta cortada desde antes
            imagen_fondo=disp.imagen_fondo; %La imagen fondo ya esta cortada desde antes
            hay_fondo=disp.hay_fondo;
            if(hay_fondo==true)
                imagen_base=filtro1(imagen_base,imagen_fondo);
            end
            %OJO EN LA SIGUIENTE LINEA, Es valido cortar el largo del canal a la mitad? mejor seria verlo todo creo yo
            %imagen_base=imagen_base(obj.posicion_dispositivo(1,2):obj.posicion_dispositivo(1,2)+(disp.height-3),obj.posicion_dispositivo(1,1):obj.posicion_dispositivo(1,1)+round(disp.long/2));
            imagen_base=imagen_base(obj.posicion_dispositivo(1,2):obj.posicion_dispositivo(1,2)+(disp.height-3),obj.posicion_dispositivo(1,1):obj.posicion_dispositivo(1,1)+disp.long); 
            xi=[1 size(imagen_base,2)];
            yi=[ceil(size(imagen_base,1)/2.0) ceil(size(imagen_base,1)/2.0)];
            perfil_base=improfile(imagen_base,xi,yi);
            perfil_base=rescale(perfil_base,0,1);
            %autocorrelacion=xcorr(perfil_base,perfil_base,'normalized');
            %max_base=max(autocorrelacion);
            contador=0;
            
            a=obj.posicion_dispositivo(1,2);
            b=obj.posicion_dispositivo(1,1);
            direc=disp.direc_imagenes;
            height=disp.height;
            long=disp.long;
            bloqueo_2=obj.bloqueo;
            volteado=disp.volteado;
            corte_imagen=obj.position_izquierda_superior;
            resize_x=disp.resize_x;
            resize_y=disp.resize_y;
            dimension=size(disp.imagen_base_ori);
            parfor i=1:disp.num_de_imagenes
                imagen_test=imread(direc(1,i));
                if (volteado==true)
                    imagen_test=imrotate(imagen_test,90);
                end

                %primero tomaremos el caso sin resize
                if resize_x==true
                    imagen_test=imresize(imagen_test, [512 dimension(1,2)],"bicubic");
                end

                if resize_y==true
                    imagen_test=imresize(imagen_test, [dimension(1,1) 512],"bicubic");
                end
                imagen_test=imagen_test(corte_imagen(1,2):corte_imagen(1,2)+511,corte_imagen(1,1):corte_imagen(1,1)+511);
                
                if (hay_fondo==true)
                    imagen_test=filtro1(imagen_test,imagen_fondo);
                end 
                
                %LO MISMO AQUI, mejor ver todo el largo del canal en mi
                %opinion
                %imagen_test=imagen_test(a:a+height,b:b+round(long/2));
                imagen_test=imagen_test(a:a+height-3,b:b+long);
                perfil_test=improfile(imagen_test,xi,yi);
                perfil_test=rescale(perfil_test,0,1);
                correlacion=xcorr(perfil_base,perfil_test,'normalized');
                max_test=max(correlacion);
                dif_porcen=(1-max_test)*100;
                if dif_porcen>2.5
                    contador=contador+1;
                    bloqueo_2(1,i)=true;
                end

            end

            % toc
            % beep
            obj.bloqueo=bloqueo_2;
            obj.num_frame_eventos=contador;
            contador
            msgbox("Se ha terminado de encontrar los eventos del test", ...
                    "Éxito");
        end
        
        
        % Function to 
        function create_events(obj)
              
            frame_inicial = 0;
            frame_final = 0;
            contador = 1;
            sin_evento = true;
            
            while(sin_evento)
                if(contador>1)
                    if(obj.bloqueo(contador) == 1 && obj.bloqueo(contador-1) == 0)
                        frame_inicial = contador-1;
                    end

                    if(obj.bloqueo(contador) == 0 && obj.bloqueo(contador-1) == 1)
                        frame_final = contador-1;
                        obj.Eventos = Evento(frame_inicial,frame_final,obj);
                        sin_evento = false;
                    end
                end
                contador = contador + 1;
            end

            for i=contador:1:length(obj.bloqueo)
                if(obj.bloqueo(i) == 1 && obj.bloqueo(i-1) == 0)
                    frame_inicial = i-1;
                end

                if(obj.bloqueo(i) == 0 && obj.bloqueo(i-1) == 1)
                    frame_final = i-1;
                    obj.Eventos = [obj.Eventos Evento(frame_inicial,frame_final,obj)];
                end

            end
            
        end
        
        
        % Function 
        function set_linea_alpha(obj) 
            
            figure
            imagen_base=obj.Dispositivo.imagen_base;  %La imagen base ya se guardo cortada en Dispositivo
            imshow(imagen_base);
            
            x_45=[obj.posicion_pendiente(1,1)-100,obj.posicion_pendiente(1,1)+100];
            y_45=[obj.posicion_pendiente(1,2)-100,obj.posicion_pendiente(1,2)+100];
            x_p=[obj.posicion_constriccion(1,1)-200,obj.posicion_constriccion(1,1)+200];
            y_p=[obj.posicion_constriccion(1,2),obj.posicion_constriccion(1,2)];
        
            [x,~] = polyxpoly(x_p,y_p,x_45,y_45);
        
            obj.linea = round(x);

            u = [x_p(2)-x_p(1) y_p(2)-y_p(1)];
            v = [x_45(2)-x_45(1) y_45(2)-y_45(1)];

            obj.alpha = rad2deg(acos(dot(u,v)/(norm(u)*norm(v))));
            drawline("Position",[obj.posicion_pendiente(1,1)-50,obj.posicion_pendiente(1,2)-50;obj.posicion_pendiente(1,1)+50,obj.posicion_pendiente(1,2)+50])
            drawline("Position",[obj.posicion_constriccion(1,1)-200,obj.posicion_constriccion(1,2);obj.posicion_constriccion(1,1)+200,obj.posicion_constriccion(1,2)])
            drawline("Position",[obj.linea,1;obj.linea,500])
            drawline("Position",[obj.posicion_dispositivo(1,1),obj.posicion_dispositivo(1,2);obj.posicion_dispositivo(1,1),obj.posicion_dispositivo(1,2)+obj.Dispositivo.height])
        end
        

        % Function to define de radius of the outer circumference, it
        % depends on the size of the channel
        function set_rho(obj)
            
            % Obtaining the size of channel
            tamano_canal = obj.Dispositivo.tamano_canal;
            
            if(tamano_canal == 5)
                obj.rho = 8.85;
            elseif(tamano_canal == 7.5)
                obj.rho = 13.275;
            elseif(tamano_canal == 10)
                obj.rho = 17.7;
            else
                error("The size of the channel can only be: 5, 7.5 or 10.")
            end

        end

    end

end

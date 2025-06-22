% CLASS "Evento" - Defined to store all the information about a set of
% frames in which something is passing through the constriction

classdef Evento < handle

    properties 
        frame_inicial   % Unsigned Integer - Frame 1 inicial del evento
        frame_final     % Unsigned Integer - Frame N final del evento
        is_cell         % Booleano - Indica si o no una celula 
        data_xf_xb      % Array Double - 
        Rc_real         % Double - Radio real de la celula  
        estimated_Rc    % Array Double - Dimentionless estimated Rc
        Rc_previo       % Double - Radio de la celula a partir de area al 
                        % ser deformada
        Ac              % Double - Area de la celula siendo deformada |CAMBIAR A ARRAY PARA AREA EN CADA FRAME
        x0_plana        % Double - x0 para el caso 1 en que la celula es
                        % tangente a la zona plana
        x0_cilin        % Double - x0 para el caso 2 en que la celula es
                        % tangente a la zona cilindrica|
        vector_Al       % Array de Double - Array de posicion del frente de
                        % la celula a lo largo del evento
        vector_back     % Array de Double - Array de posicion del frente 
                        % trasero de la celula a lo largo del evento
        Al_pl           % Array de Double - Array de frentes de la celula 
                        % normalizado habiendo restado xo_plana
        Al_cl           % Array de Double - Array de frentes de la celula 
                        % normalizado habiendo restado xo_cilin
        data_pl         % Array de Double - Array con cada par de datos
                        % tiempo y Al-x0_plana/Wch
        data_cl         % Array de Double - Array con cada par de datos
                        % tiempo y Al-x0_cilin/Wch
        data_electric   % Array Double - Impedance, phase and time of the 
                        % event (empty if not cell)
        C2              % Double - Capacitance of the cell
        R2              % Double - Resistance of the cell
        R_leak          % Double - Leak resistance due to the cell passing   
                        % through the constriction    
        electric_pl     % Array Double - Electric data associated to data_pl
        electric_cl     % Array Double - Electric data associated to data_cl
        E_0
        g_infinite
        tau_C
        Canal           % Canal - Hace referencia al canal del evento
    end


    methods

        function obj = Evento(frame_inicial,frame_final,Canal)
            obj.frame_inicial = frame_inicial;
            obj.frame_final = frame_final;
            obj.data_xf_xb = [
                0.365005  -2.83499  -3.64513  -4.48023  -6.12355
                0.457524  -2.74713  -3.57887  -4.42076  -6.07289
                0.510148  -2.69954  -3.54118  -4.38693  -6.04408
                0.556687  -2.65840  -3.50785  -4.35701  -6.01860
                0.599775  -2.62100  -3.47699  -4.32932  -5.99501
                0.640978  -2.58579  -3.44748  -4.30283  -5.97245
                0.692357  -2.54348  -3.41069  -4.26980  -5.94432
                0.738109  -2.50686  -3.37792  -4.24039  -5.91927
                0.784486  -2.47018  -3.34470  -4.21058  -5.89388
                0.831229  -2.43371  -3.31123  -4.18053  -5.86828
                0.880206  -2.39578  -3.27615  -4.14905  -5.84147
                0.931670  -2.35416  -3.23929  -4.11597  -5.81329
                0.990808  -2.30850  -3.19694  -4.07795  -5.78091
                1.057000  -2.25767  -3.14953  -4.03540  -5.74467
                1.133260  -2.19941  -3.09492  -3.98638  -5.70292
                1.232490  -2.12359  -3.02384  -3.92259  -5.64859
                1.238330  -2.12052  -3.01967  -3.91884  -5.64539
                1.241510  -2.11758  -3.01739  -3.91680  -5.64365
                1.244670  -2.11537  -3.01513  -3.91477  -5.64192
                1.248100  -2.11277  -3.01267  -3.91256  -5.64004
                1.252640  -2.10933  -3.00942  -3.90964  -5.63756
                1.258810  -2.10462  -3.00500  -3.90567  -5.63418
                1.266770  -2.09849  -2.99929  -3.90056  -5.62982
                1.276430  -2.09104  -2.99238  -3.89435  -5.62453
                1.287450  -2.08250  -2.98448  -3.88726  -5.61849
                1.299400  -2.07327  -2.97593  -3.87959  -5.61196
                1.311800  -2.06374  -2.96704  -3.87161  -5.60516
                1.324320  -2.05411  -2.95808  -3.86356  -5.59831
                1.336850  -2.04445  -2.94910  -3.85551  -5.59145
                1.349250  -2.03488  -2.94023  -3.84754  -5.58466
                1.361540  -2.02541  -2.93142  -3.83964  -5.57793
                1.373850  -2.01588  -2.92260  -3.83172  -5.57119
                1.386300  -2.00631  -2.91369  -3.82373  -5.56438
                1.399200  -1.99639  -2.90445  -3.81543  -5.55731
                1.413070  -1.98571  -2.89452  -3.80652  -5.54972
                1.428480  -1.97381  -2.88348  -3.79661  -5.54128
                1.445950  -1.96032  -2.87097  -3.78538  -5.53172
                1.466400  -1.94441  -2.85632  -3.77224  -5.52052
                1.490350  -1.92587  -2.83917  -3.75684  -5.50741
                1.518710  -1.90385  -2.81886  -3.73861  -5.49188
                1.553080  -1.87701  -2.79424  -3.71651  -5.47306
                1.596030  -1.84350  -2.76348  -3.68890  -5.44954
                1.651010  -1.80033  -2.72410  -3.65356  -5.41944
            ];
            obj.Rc_real = 0;
            obj.Rc_previo = 0;
            obj.vector_Al = zeros(frame_final-frame_inicial+1,2);
            obj.vector_back = zeros(frame_final-frame_inicial+1,2);
            obj.Al_pl = zeros(frame_final-frame_inicial+1,1);
            obj.Al_cl = zeros(frame_final-frame_inicial+1,1);
            %convertiendo area en un array para cada frame del evento
            obj.Ac = zeros(frame_final-frame_inicial+1,1);
            obj.Canal = Canal;
            obj.E_0 = zeros(1,2);
            obj.g_infinite = zeros(1,2);
            obj.tau_C = zeros(1,2);
        end

           
        function calc_x0_plana(obj)
            obj.x0_plana = obj.Rc_real + (obj.Canal.wch/tan(obj.Canal.alpha)) - (obj.Rc_real/sin(obj.Canal.alpha));
        end


        function calc_x0_cilin(obj)

            %obj.x0_cilin = obj.Rc_real - (sqrt((obj.Rc_real + obj.Canal.rho)^2 - (obj.Canal.rho - obj.Canal.wch)^2) - obj.Canal.rho*tan(obj.Canal.alpha/2));
            obj.x0_cilin = obj.Rc_real - (sqrt((obj.Rc_real + obj.Canal.rho)^2 - (obj.Canal.rho + obj.Canal.wch)^2) - obj.Canal.rho*tan(obj.Canal.alpha/2));
        end

        
        function is_cell_and_radious(obj)
            %tomando en consideracion las imagenes de francia cuando
            %volteado==true

            volteado=obj.Canal.Dispositivo.volteado;
            resize_x=obj.Canal.Dispositivo.resize_x;
            resize_y=obj.Canal.Dispositivo.resize_y;
            corte_imagen=obj.Canal.position_izquierda_superior;
            dimension=size(obj.Canal.Dispositivo.imagen_base_ori);
            area_compensacion=1.0;
            valor = 0;

            for k = 1:5
                direccion_frame = obj.Canal.Dispositivo.direc_imagenes(obj.frame_inicial+k);
                frame_evento = imread(direccion_frame);
                
                if(volteado==true)
                    frame_evento = imrotate(frame_evento,90);
                end
                if resize_x==true
                    frame_evento=imresize(frame_evento, [512 dimension(1,2)],"bicubic");
                    compensar=512/dimension(1,1);
                    area_compensacion=area_compensacion*compensar;
                end
                if resize_y==true
                    frame_evento=imresize(frame_evento, [dimension(1,1) 512],"bicubic");
                    compensar=512/dimension(1,2);
                    area_compensacion=area_compensacion*compensar;
                end
    
                frame_evento=frame_evento(corte_imagen(1,2):corte_imagen(1,2)+511,corte_imagen(1,1):corte_imagen(1,1)+511);
   
                modelo = load("resunet_propia_v7.mat");
                modelo = modelo.net_n;
                C = semanticseg(frame_evento,modelo);
                celula=C=="C3";
                hay_celula=find(celula==1);
                hay_celula=size(hay_celula,1);
    
                if(hay_celula==0)
                    continue
                end
                
                [regiones,n]= bwlabel(celula,4);
                
                %HAY QUE COMPESAR EL AREA SI ES QUE SE HIZO UN RESIZE
                propiedades=regionprops(regiones==1,"FilledArea");
                areas=zeros(n,1);
                areas(1)=propiedades.FilledArea;
                areas(1)=areas(1)/area_compensacion;
                for i=2:1:n
                    region_n=regiones==i;
                    propiedades=regionprops(region_n,"FilledArea");
                    areas(i)=propiedades.FilledArea;
                    areas(i)=areas(i)/area_compensacion;
                end
    
                [valor,~] = max(areas);
    
                if ((~(valor < 1500 || valor > 3000) && ~volteado)||(~(valor < 2200 || valor > 3200) && volteado))
                    obj.is_cell=1;
                    obj.Ac(k) = valor;
                    obj.Rc_previo = round(sqrt(obj.Ac(k)/pi),3,"significant");
                    return
                end
            end    
         
            obj.Ac(5) = valor; %El area ya esta compensada de antes, no es necesario volver a dividir
           
            obj.is_cell=0;
          
        end 
        

        function cell_front_back(obj)

            model = load("resunet_propia_v7.mat");
            model = model.net_n;
            volteado=obj.Canal.Dispositivo.volteado;
            resize_x=obj.Canal.Dispositivo.resize_x;
            resize_y=obj.Canal.Dispositivo.resize_y;
            corte_imagen=obj.Canal.position_izquierda_superior;
            dimension=size(obj.Canal.Dispositivo.imagen_base_ori);
            area_compensacion=1.0;
            for i = obj.frame_inicial:obj.frame_final
                frame_path = obj.Canal.Dispositivo.direc_imagenes(i+1);
                frame = imread(frame_path);
                if(volteado==true)
                    frame = imrotate(frame,90);
                end
                if resize_x==true
                    frame=imresize(frame, [512 dimension(1,2)],"bicubic");
                    compensar=512/dimension(1,1);
                    area_compensacion=area_compensacion*compensar;

                end
                if resize_y==true
                    frame=imresize(frame, [dimension(1,1) 512],"bicubic");
                    compensar=512/dimension(1,2);
                    area_compensacion=area_compensacion*compensar;
                end

                frame=frame(corte_imagen(1,2):corte_imagen(1,2)+511,corte_imagen(1,1):corte_imagen(1,1)+511);

                C = semanticseg(frame,model);
                cell = C == "C3";

                [regions, ~] = bwlabel(cell, 4);
                properties = regionprops(regions, "FilledArea");

                areas = [properties.FilledArea]/area_compensacion;

                [value,pos] = max(areas);
                if (isempty(areas)==1)
                    continue
                end
                
                if (((value < 1500 || value > 3000) && volteado==false)||((value < 2200 || value > 3200) && volteado==true))
                    continue
                else
                    
                    
                    
                %Aqui se suman ciertas numeros a las localizaciones, un
                %poco raro, creo que seran los numeros de cuando haciamos
                %cortes arbitrarios
                    cell_frame = regions == pos;
                    % cons_start = obj.Canal.posicion_constriccion(1,2) + 3 - 13;
                    cons_start = obj.Canal.posicion_constriccion(1,2) + 3;
                    % cons_end = obj.Canal.posicion_constriccion(1,2) + obj.Canal.wch*2 - 3 - 13;
                    cons_end = obj.Canal.posicion_constriccion(1,2) + obj.Canal.wch*2 - 3;

                    cons = (cons_start:cons_end)';

                    x_interested = cell_frame(cons,:);

                    cell_front = zeros(length(cons),1);
                    cell_back = zeros(length(cons),1);

                    for k = 1:1:length(cons)
                        %Todo este for lo he cambiado con lo de iñaki
                        f = find(x_interested(k,:),1,"last");
                        b = find(x_interested(k,:),1,"first");

                        if ~isempty(f)
                            cell_front(k) = f;
                        end
                    
                        if ~isempty(b)
                            cell_back(k) = b;
                        end
                    end

                    [Al,y1] = parab(cell_front,cons,"front");
                    [back,y2] = parab(cell_back,cons,"back");

                    Al = Al - (obj.Canal.linea);
                    back = back - (obj.Canal.linea);

                    obj.vector_Al(i-obj.frame_inicial+1,:) = [Al,y1];
                    obj.vector_back(i-obj.frame_inicial+1,:) = [back,y2];
                    %metiendo las areas en el array
                    obj.Ac(i-obj.frame_inicial+1)=value;
                end
            end
        end


        % Estimating cell size from images in the constriction
        % April 2025, A. Abarca, I. Fraga, D. Ñaña, B. González, G.R. Plaza
        function estimating_Rc(obj)

            if(isempty(obj.vector_Al) || isempty(obj.vector_back))
                return
            end
            
            obj.estimated_Rc = zeros(size(obj.vector_Al,1),1);

            for i = 1:size(obj.vector_Al,1)
                
                xf = obj.vector_Al(i,1);
                xb = obj.vector_back(i,1);

                if(xf == 0 || xb == 0 || xf == xb || isfinite(xf) || isfinite(xb))
                    continue
                end

                xf = xf/obj.Canal.wch;
                xb = xb/obj.Canal.wch;

                if(xf > obj.data_xf_xb(end,1))
                    continue
                end

                % Evaluating interpolations
                xb_1_6 = interp1(obj.data_xf_xb(:,1), obj.data_xf_xb(:,2), xf, 'pchip', 'extrap');
                xb_2_0 = interp1(obj.data_xf_xb(:,1), obj.data_xf_xb(:,3), xf, 'pchip', 'extrap');
                xb_2_4 = interp1(obj.data_xf_xb(:,1), obj.data_xf_xb(:,4), xf, 'pchip', 'extrap');
                xb_3_2 = interp1(obj.data_xf_xb(:,1), obj.data_xf_xb(:,5), xf, 'pchip', 'extrap');
                
                vec_Rc = [1.6, 2.0, 2.4, 3.2];
                vec_xb = [xb_1_6, xb_2_0, xb_2_4, xb_3_2];
                
                obj.estimated_Rc(i) = interp1(vec_xb, vec_Rc, xb, 'pchip'); % Valor estimado de Rc/Wch
            end

            all_Rc = nonzeros(obj.estimated_Rc);
            
            if(~isnan(all_Rc))
                obj.Rc_real = mean(all_Rc)*obj.Canal.wch;
            end

        end


        function calc_final(obj)

            if(obj.Rc_real == 0)
                obj.Rc_real = obj.Rc_previo;
            end

            obj.calc_x0_plana();
            obj.calc_x0_cilin();
            
            t_0 = obj.frame_inicial * (1/obj.Canal.Dispositivo.fs);
            t = (obj.frame_inicial:obj.frame_final) * (1/obj.Canal.Dispositivo.fs) - t_0;

            if(obj.x0_plana > 0)
                obj.Al_pl = (obj.vector_Al(:,1)-obj.x0_plana)/obj.Canal.wch;
                idx = obj.Al_pl > 0;
                obj.data_pl = zeros(nnz(idx),2);
                obj.data_pl(:,1) = t(idx);
                obj.data_pl(:,2) = obj.Al_pl(idx);
            elseif(obj.x0_cilin > 0)
                obj.Al_cl = (obj.vector_Al(:,1)-obj.x0_cilin)/obj.Canal.wch;
                idx = obj.Al_cl > 0;
                obj.data_cl = zeros(nnz(idx),2);
                obj.data_cl(:,1) = t(idx);
                obj.data_cl(:,2) = obj.Al_cl(idx);
            end

        end 

        
        % Function to set C2, R2 and R_leak
        % Non-linear fitting using the Levemberg-Marquardt algorithm
        % Created December 2024, I. Fraga, D. Ñaña, B. González, G. Plaza
        % Last Edited April 2025 
        % Fitting of electrical parameters
        function [] = set_C2_R2(obj,num)

            % Getting the data where no cells are passing through (needs
            % modifying to obtain that data)
            data_freq = obj.Canal.Dispositivo.ElectricData.frequencies;

            if(isempty(obj.data_electric))
                warning("No impedance data for this event...")
                return
            end

            freq_Imp = data_freq;
            
            obj.C2 = zeros(size(obj.data_electric,1),2);
            obj.R2 = zeros(size(obj.data_electric,1),2);
            obj.R_leak = zeros(size(obj.data_electric,1),2);

            obj.C2(:,1) = obj.data_electric(:,1)/1000;
            obj.R2(:,1) = obj.data_electric(:,1)/1000;
            obj.R_leak(:,1) = obj.data_electric(:,1)/1000;

            for l = 1:size(obj.data_electric,1)
                argument_Imp = obj.data_electric(l,4:5); % argument of impedance in deg
                modulus_Imp = obj.data_electric(l,2:3); % modulus of impedance in Ohm
                
                if(num == 1)
                    plot(freq_Imp,modulus_Imp, '*');
                end
    
                % Extracting the vectors of interest from the matrix
                Real_Exp = modulus_Imp.*cos(argument_Imp*pi/180);
                Imag_Exp = modulus_Imp.*sin(argument_Imp*pi/180);
    
                num_parameters = 3;
                cal_R1 = obj.Canal.Dispositivo.ElectricData.R1;
                inv_C1 = 1/obj.Canal.Dispositivo.ElectricData.C1;
                Beta=[cal_R1/2 inv_C1/2 cal_R1/2]; % Beta(1) = R2; 
                                                   % Beta(2) = inv_C2; 
                                                   % Beta(3) = R_leak;
    
                num_datos = length(freq_Imp);
                num_iter = 5000;
                levenberg = 0.1; % Levenberg-Marquardt parameter
                All_Exp = [Real_Exp(:); Imag_Exp(:)];
                [Real_Teo, Imag_Teo] = obj.Impedance_theor_event(cal_R1,inv_C1,Beta(1),Beta(2),Beta(3),freq_Imp);
                All_Teo = [Real_Teo(:); Imag_Teo(:)];
                DifferenceAll = All_Teo - All_Exp;
                % new_error = dot(DifferenceAll,DifferenceAll)/dot(All_Exp,All_Exp);
    
                for j = 1:num_iter
    
                    epsilon_beta = Beta/100;
                    J = zeros(2*num_datos,num_parameters);
    
                    for k = 1:num_parameters
	                    for i = 1:num_datos
                            [Real_2, Imag_2] = obj.Impedance_theor_event(cal_R1,inv_C1, Beta(1)+obj.delta_event(1,k)*epsilon_beta(1), ...
                                                               Beta(2)+obj.delta_event(2,k)*epsilon_beta(2), ...
                                                               Beta(3)+obj.delta_event(3,k)*epsilon_beta(3),freq_Imp(i));
                            [Real_1, Imag_1] = obj.Impedance_theor_event(cal_R1,inv_C1,Beta(1),Beta(2),Beta(3),freq_Imp(i));
		                    J(i,k) = (Real_2 - Real_1)/epsilon_beta(k);
                            J(i + num_datos,k) = (Imag_2 - Imag_1)/epsilon_beta(k);
	                    end
                    end
    
                    deltaBeta = transpose(-inv(transpose(J)*J + levenberg*eye(num_parameters))*(transpose(J)*DifferenceAll));
                    Beta = Beta + deltaBeta;
    
                    [Real_Teo, Imag_Teo] = obj.Impedance_theor_event(cal_R1,inv_C1,Beta(1),Beta(2),Beta(3),freq_Imp);
                    All_Teo = [Real_Teo(:); Imag_Teo(:)];
                    DifferenceAll = All_Teo - All_Exp;
                    new_error = dot(DifferenceAll,DifferenceAll)/dot(All_Exp,All_Exp);
    
                    if(num == 1)
                        fprintf('\nIteración %d - Error: %f\n', j, new_error);
                        fprintf('\nBeta(1): %f, Beta(2): %f\n', Beta(1), Beta(2));
                    end
                end
                
                if(num == 1)
                    % Plot the results
                    freq_Simula = min(freq_Imp):20:max(freq_Imp);
                    [Real_Simula, Imag_Simula] = obj.Impedance_theor_event(R1,invC1, Beta(1),Beta(2),Beta(3),freq_Simula);
                    figure;
                    plot(freq_Simula,Real_Simula,'-');
                    hold on;
                    plot(freq_Imp,Real_Exp,'*');
                    figure;
                    plot(freq_Simula,Imag_Simula,'-');
                    hold on;
                    plot(freq_Imp,Imag_Exp,'*');
                end
    
                obj.R2(l,2) = Beta(1);
                obj.C2(l,2) = 1/Beta(2);
                obj.R_leak(l,2) = Beta(3);
            end

        end  


        % Impedance as a function of the five parameters
        function [Real_Imp, Imag_Imp] = Impedance_theor_event(~, R1, invC1, R2, invC2, R3, freq_Imp)
           % Freq_Imp: frequency in Hz
           % R1, R2, R3: resistences in Ohm
           % C1, C2: capacitances in F
           C1 = 1/invC1;
           calc_C2 = 1/invC2;
           w = 2*pi*freq_Imp;
           Real_Imp = R1+(R3+R2*R3*(R2+R3)*((w*calc_C2).^2))./(1+(((R2+R3)*w*calc_C2).^2));
           Imag_Imp = -1./(w*C1)-((R3)^2)*w*calc_C2./(1+(((R2+R3)*w*calc_C2).^2));
        end 


        % Delta function to be used in obtainment of C1 and R1
        function g = delta_event(~, a, b)
            if (a == b)
                g = 1;
            else
                g = 0;
            end
        end

        
        % Function to centering electric parameters to the mechanical ones
        function [] = centering_electric_data(obj)
            
            if(isempty(obj.data_pl) && isempty(obj.data_cl))
                warning("No mechanical data to center the electric parameters...");
                return
            end

            t_0 = obj.frame_inicial * (1/obj.Canal.Dispositivo.fs);
            t = (obj.frame_inicial:obj.frame_final) * (1/obj.Canal.Dispositivo.fs) - t_0;

            electric_times = obj.data_electric(:,1)/1000;
            
            for i = 1:2
                
                if(i == 1 && ~isempty(obj.Al_pl) && ~isempty(obj.data_pl))
                    Al_vector = obj.Al_pl;
                    t_data = obj.data_pl(:,1);
                elseif(i == 2 && ~isempty(obj.Al_cl) && ~isempty(obj.data_cl))
                    Al_vector = obj.Al_cl;
                    t_data = obj.data_cl(:,1);
                else
                    continue
                end
                
                idx = find(Al_vector > 0,1,"first");
                t_start = t(idx);
                
                t_electric = electric_times - t_start;

                idx = (t_electric > t_data(1) & t_electric < t_data(end));
                
                if(~any(idx))
                    continue
                end
                
                times = electric_times(idx);

                valid_C2 = obj.C2(idx,2);
                valid_R2 = obj.R2(idx,2);
                valid_R_leak = obj.R_leak(idx,2);

                if(isscalar(times)) 

                    v = ones(length(t_data),1);

                    if(i == 1)
                        obj.electric_pl = [t_data v*valid_C2 ...
                                           v*valid_R2 v*valid_R_leak];
                    else
                        obj.electric_cl = [t_data v*valid_C2 ...
                                           v*valid_R2 v*valid_R_leak];
                    end

                else
                
                    indexes = obj.assign_by_time_threshold(t_data,times');

                    data = ones(length(t_data),3);

                    for j = 1:length(times)
                        data(logical(indexes(:,j)),1) = valid_C2(j);
                        data(logical(indexes(:,j)),2) = valid_R2(j);
                        data(logical(indexes(:,j)),3) = valid_R_leak(j);
                    end
                    
                    if(i == 1)
                        obj.electric_pl = [t_data data];
                    else
                        obj.electric_cl = [t_data data];
                    end
                end

            end
            
        end 

        
        % Function to asign the electric time to mechanical time axis
        function assign_indexes = assign_by_time_threshold(~,data_mec,data_elec)
        
            N = length(data_mec);
            M = length(data_elec);
        
            assign_indexes = zeros(N, M);
            thrs = [-Inf (data_elec(1:end-1)+data_elec(2:end))/2 Inf];
        
            for i = 1:N
                t = data_mec(i);
                for j = 1:M
                    if t > thrs(j) && t <= thrs(j+1)
                        assign_indexes(i,j) = 1;
                        break
                    end
                end
            end
        end


        % Analyses of aspiration curves for rectangular constrictions
        % Aldo, Blanca, Cristina, María, Gustavo. 2022

        % Adapted to MATLAB and modified May 2025, 
        % I. Fraga, D. Ñaña, B. González, G. Plaza 

        % Function to include all Aldo's Python analysis and saves 
        function [] = final_analysis(obj,blch)

            cwd = pwd;
            folder = cwd + "\results\";
            folderqs = cwd + "\sim\";

            rce = obj.Rc_real;
            wch = obj.Canal.wch;

            InputQSFilepreName = folderqs + 'rc_'; % csv format, with ",". First 
                                       % column: lambda=AL/Wc, second 
                                       % column: Pi=DP/E0. No headings
            List_of_QS_curves = ["1.6", "2.0", "2.4", "3.2"]; % List of values of 
                                                  % Rc/Wch of the curves 
                                                  % obtained by numerical 
                                                  % analyses

            OutputFileName = folder + 'Outf_';
            nb_points_QS_curves = 1000; % The numerical curves from Abaqus are 
                                        % discretized in this number of points

            Rc_to_Wch = rce/wch;

            num_cons = obj.Canal.Dispositivo.tipo_dispositivo;

            if(num_cons == 1)
                
                P_exp = obj.delta_pressure_0();

            elseif(num_cons == 3)

            elseif(num_cons == 12)

                % Relaxation function: gA + gB*math.exp(-t/tau) with gA + gB = 1.

                % Relation of P_exp to blocked channels, using Q=10ml/hr
                % rblqch = [[0,7413.9], [2, 8014.5], [4, 8972.3], [6, 10610.3], ...
                % [8, 14038.7], [10, 18644.0]]
    
                % SI ES LINEAL
                Q_ori = 10; % ml/h
                Q_exp = obj.Canal.Dispositivo.Q; 
                rel = Q_ori/Q_exp;
                P_exp = obj.delta_pressure_12(blch,rel);
            end

            if(P_exp <= 0)
                warning("P_exp = %s",num2str(P_exp))
                return
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Import and compute numerical quasistatic-curve data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            List_of_QS_curves_float = double(List_of_QS_curves);
            nb_integers = double(0:nb_points_QS_curves-1);
            lambda_values = zeros(nb_points_QS_curves, length(List_of_QS_curves));
            Pi_values = zeros(nb_points_QS_curves, length(List_of_QS_curves));

            for i = 1:length(List_of_QS_curves) 
                InputQSFileName = InputQSFilepreName + List_of_QS_curves(i) + '.txt';

                % Primera columna el tiempo, segunda la longitud aspirada normalizada
                % con el radio de la celula
                A = readmatrix(InputQSFileName, 'Delimiter', '\t');

                A_lambda_corr = A(1,1);
                A_Pi_corr = A(1,1);
                k = 1;

                for j = 2:size(A,1)
                    if A(j,1) > A_lambda_corr(k) && A(j,2) > A_Pi_corr(k)
                        A_lambda_corr(end+1,1) = A(j,1);
                        A_Pi_corr(end+1,1) = A(j,2);
                        k = k + 1;
                    end
                end

                Pi_vs_lambda_el = @(x) interpd1(A_lambda_corr, A_Pi_corr, x);
                max_lambda_el = max(A_lambda_corr);

                lambda_values(:,i) = nb_integers*max_lambda_el/ ...
                                     (nb_points_QS_curves - 1);
                Pi_values(:,i) = Pi_vs_lambda_el(lambda_values(:,i));

            end

            lambda_Pi_matrix = zeros(nb_points_QS_curves, 2);

            for i = 1:nb_points_QS_curves

                interp_lambda = @(x) interpd1(List_of_QS_curves_float, ...
                                lambda_values(i,:), x);
                lambda_Pi_matrix(i,1) = interp_lambda(Rc_to_Wch);

                interp_Pi = @(x) interpd1(List_of_QS_curves_float, Pi_values(i,:), ...
                            x);
                lambda_Pi_matrix(i,2) = interp_Pi(Rc_to_Wch);
            end

            Pi_vs_lambda_el = @(x) interpd1(lambda_Pi_matrix(:,1), ...
                              lambda_Pi_matrix(:,2), x);
            lambda_vs_Pi_el = @(x) interpd1(lambda_Pi_matrix(:,2), ...
                              lambda_Pi_matrix(:,1), x);

            max_Pi_el = max(lambda_Pi_matrix(:,2));
            max_lambda_el = max(lambda_Pi_matrix(:,1));

            figure; clf; hold on;

            for i = 1:length(List_of_QS_curves)
                plot(lambda_values(:,i), Pi_values(:,i), 'DisplayName', ...
                     sprintf('Numerical curve - R_c^*: %.2f', List_of_QS_curves(i)));
            end

            ylabel('Delta P / E_0');
            xlabel('A_L / W_ch');

            plot(lambda_Pi_matrix(:,1), lambda_Pi_matrix(:,2), '--', ...
                'DisplayName', sprintf('Extrapolated curve - R_c^*: %.2f', Rc_to_Wch));

            legend('Location', 'northwest');
            xlim([0, 3.0]);
            ylim([0, 0.8]);
            grid on;
            set(gca, 'GridAlpha', 0.5);
            legend('Box', 'on');

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Import experiment data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            for i = 1:2
                if(i == 1 && ~isempty(obj.data_pl))
                    B = obj.data_pl;
                    str_end = "_pl.csv";
                elseif(i == 2 && ~isempty(obj.data_cl))
                    B = obj.data_cl;
                    str_end = "_cl.csv";
                else
                    continue
                end

                lambda_exp_vec = B(:,2);
                t_exp_vec = B(:,1);
                max_t_exp = max(B(:,1));
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Run of the fitting process
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                E0_ini = P_exp/ ...
                (Pi_vs_lambda_el(mean(lambda_exp_vec(1:4)))); % We consider instant 
                % deformation equal to the average of the _ just the first measurement
                gA_ini = 0.05;
                tau_ini = max_t_exp*1;

                disp(E0_ini)
    
                [E0, gA, tau, lambda_vs_t_fit, best_error] = obj.fitting_nlls(t_exp_vec, ...
                                                             lambda_exp_vec, P_exp, ...
                                                             E0_ini, gA_ini, tau_ini, ...
                                                             lambda_vs_Pi_el, max_Pi_el,...
                                                             max_lambda_el);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Export the results
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                fid = fopen((OutputFileName + str_end), 'w');

                FirstLine = sprintf('E0 = %.6f; gA = %.6f; tau = %.6f\n', E0, gA, tau);
                fprintf(fid, FirstLine);

                fprintf(fid, 'time (s);AL/Wc exp;AL/Wc fit\n');

                for j = 1:length(lambda_exp_vec)
                    NewLine = sprintf('%.6f;%.6f;%.6f\n', t_exp_vec(j), ...
                                      lambda_exp_vec(j), lambda_vs_t_fit(t_exp_vec(j)));
                    fprintf(fid, NewLine);
                end

                fclose(fid);

                disp([E0, gA, tau, best_error]);
                disp(E0_ini);
                disp(P_exp);

                figure; clf; hold on;

                plot(t_exp_vec, lambda_exp_vec, '-b', 'DisplayName', ...
                     'Experimental curve', 'LineWidth', 1.2, 'Color', [0, 0, 1, 0.5]);
                scatter(t_exp_vec, lambda_exp_vec, 20, 'filled', 'MarkerFaceAlpha', ...
                        0.5, 'MarkerEdgeAlpha', 0.5);
                plot(t_exp_vec, lambda_vs_t_fit(t_exp_vec), '-r', 'LineWidth', 1.5, ...
                     'DisplayName', sprintf(['Fitted curve, E_0: %.2f [Pa], g_∞: ' ...
                     '%.3f, τ_C: %.3f [s], Error: %.3f %%'], E0, gA, tau, best_error));
                xlabel('Time (s)');
                ylabel('A_L/W_ch');

                xlim([0, max(t_exp_vec)*1.2]);
                ylim([0, max(lambda_exp_vec)*1.2]);

                grid on;
                set(gca, 'GridAlpha', 0.5);
                legend('Location', 'best', 'Box', 'on');

                obj.E_0(1,i) = E0;
                obj.g_infinite(1,i) = gA;
                obj.tau_C(1,i) = tau;

            end
            
        end    


        % Function for eQLV simulating the viscous entry process
        function [lambda_vec, lambda_vs_t_teo, t_vs_lambda_teo, ...
                 t_max_teo] = theoret_aspiration(~, lambda_vs_Pi_el, ...
                 max_Pi_el, P_exp, E0, gA, tau_fit)

            gB = 1-gA;
            tau = 1;

            % Relaxation function: gA + gB*math.exp(-t/tau)
            Pi_max = P_exp/E0; % Maximun value (constant during the creep 
                               % process) of the non-dimensional differential
                               % pressure

            if(Pi_max >= max_Pi_el)
                disp(['E0 is too low in the computations, since the ' ...
                      'non-dimensional aspiration pressure exceeds the ' ...
                      'maximum value in que quasistatic curve'])
                input("Press any key to continue")
            end

            function val = creep_f(t)
                val = (1+gB/gA)-(gB/gA)*exp(-t/(tau*(1+gB/gA)));
            end

            T_max = -tau*(1+gB/gA)*log((gA/gB)*((gB/gA+1)-max_Pi_el/Pi_max));
            nb_points_creep = 20;
            T_increment = T_max/(nb_points_creep-1);

            lambda_vec = zeros(1,nb_points_creep); % Vector of the non-dimensional 
                                                   % advance of the cell lambda = 
                                                   % AL/Wc
            Pi_el_vec = zeros(1,nb_points_creep); % Vector of the non-dimensional 
                                                  % differential pressure Pi=DP/E0
            time_T_vec = zeros(1,nb_points_creep); % Vector of non-dimensional time
                                                   % T

            for k = 1:nb_points_creep
                T_value = (k-1)*T_increment;
                time_T_vec(k) = T_value;
                Pi_el_vec(k) = Pi_max*creep_f(T_value);
                lambda_vec(k) = lambda_vs_Pi_el(Pi_el_vec(k));
            end

            time_t_vec = tau_fit*time_T_vec;

            t_vs_lambda_teo = @(lambda) interpd1(lambda_vec, ...
                                                time_t_vec, ...
                                                lambda);
            lambda_vs_t_teo = @(t) interpd1(time_t_vec, ...
                                           lambda_vec, t);

            t_max_teo = max(time_t_vec);

        end  


        % Function for computing mechanical parameters by fitting non-linear 
        % least-squares method based on a script by G. Guinea
        function [bestBeta_0, bestBeta_1, bestBeta_2, bestlambda_vs_t_teo2, ...
                 best_error] = fitting_nlls(obj, t_exp_vec, lambda_exp_vec, P_exp, ...
                 E0_ini, gA_ini, tau_ini, lambda_vs_Pi_el, max_Pi_el, ...
                 max_lambda_el)

            % This gA_ini not used, because a set of different values are 
            % tried as initial value...

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INICIALIZACION DE DATOS PARA QUE NO EXPLOTE EL CODIGO
            bestBeta_0 = NaN;
            bestBeta_1 = NaN;
            bestBeta_2 = NaN;
            bestlambda_vs_t_teo2 = [];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            y = t_exp_vec;
            x = lambda_exp_vec;
            num_data = length(x);
            best_error = 10000000;
            iter = 50;

            num_param = 3;

            for gA_index = 0:9

                gA_ini = gA_index/10;

                fprintf("\ngA_ini %s\n",num2str(gA_ini))

                beta = [E0_ini, gA_ini, tau_ini];
                epsilon_beta = 1e-5*beta;

                levenberg = 100;
                levenberg_factor_mult = 2;
                levenberg_factor_div = 2;

                if(beta(1) <= P_exp/max_Pi_el)
                    beta(1) = 1.0001*P_exp/max_Pi_el;
                end
                max_gA_possible = (P_exp/(beta(1)*max_Pi_el))*0.5;

                if(beta(2) >= max_gA_possible)
                    beta(2) = max_gA_possible*0.95;
                end

                if(beta(2) < 1.e-2)
                    beta(2) = 1.e-2;
                end

                if(beta(3) < 1.e-3)
                    beta(3) = 1.e-3;
                end

                [~, ~, t_vs_lambda_teo, ~] = ...
                obj.theoret_aspiration(lambda_vs_Pi_el, max_Pi_el, P_exp, ...
                beta(1), beta(2), beta(3));     

                R_2 = 0;
                y_2 = 0;
                R = zeros(num_data,1);

                for i = 1:num_data
                    R(i) = t_vs_lambda_teo(x(i))-y(i);
                    R_2 = R_2 + R(i)*R(i);
                    y_2 = y_2 + y(i)*y(i);
                end

                prev_error = 100*R_2/y_2;

                fprintf("\nFitting error: %s; Beta: %s\n", ...
                        num2str(prev_error),num2str(beta));

                J = zeros(num_data, num_param);

                for j = 1:iter

                    for k = 1:num_param
                        epsilon_beta = beta*1.e-3;
                        E01 = beta(1) + obj.delta_event(k,1)*epsilon_beta(1);
                        gA1 = beta(2) + obj.delta_event(k,2)*epsilon_beta(2);
                        tau1 = beta(3) + obj.delta_event(k,3)*epsilon_beta(3);

                        [~, ~, t_vs_lambda_teo1, ~] = ...
                        obj.theoret_aspiration(lambda_vs_Pi_el, max_Pi_el, ...
                        P_exp, E01, gA1, tau1);

                        E02 = beta(1);
                        gA2 = beta(2);
                        tau2 = beta(3);

                        [~, ~, t_vs_lambda_teo2, ~] = ...
                        obj.theoret_aspiration(lambda_vs_Pi_el, max_Pi_el, ...
                        P_exp, E02, gA2, tau2);

                        for i = 1:num_data
                            J(i,k) = (t_vs_lambda_teo1(x(i))- ...
                                      t_vs_lambda_teo2(x(i)))/epsilon_beta(k);
                        end
                    end

                    delta_beta = transpose(-inv(transpose(J)*J + levenberg* ...
                                 eye(num_param))*(transpose(J)*R));

                    for k = 1:num_param
                        beta(k) = beta(k) + delta_beta(k);
                    end

                    if(beta(1) < P_exp/max_Pi_el)
                        beta(1) = 1.001*P_exp/max_Pi_el;
                    end

                    max_gA_possible = P_exp/(beta(1)*max_Pi_el);

                    if(beta(2) > max_gA_possible)
                        beta(2) = max_gA_possible*0.95;
                    end

                    if(beta(2) < 1.e-8)
                        beta(2) = gA1;
                    end

                    if(beta(3) < 1.e-8)
                        beta(3) = tau1;
                    end

                    R_2 = 0;
                    y_2 = 0;

                    E02 = beta(1);
                    gA2 = beta(2);
                    tau2 = beta(3);

                    [~, lambda_vs_t_teo2, t_vs_lambda_teo2, ~] = ...
                    obj.theoret_aspiration(lambda_vs_Pi_el, max_Pi_el, P_exp, ...
                    E02, gA2, tau2);

                    for i = 1:num_data
                        R(i) = t_vs_lambda_teo2(x(i))-y(i);
                        R_2 = R_2 + R(i)*R(i);
                        y_2 = y_2 + y(i)*y(i);
                    end

                    now_error = 100*R_2/y_2;

                    disp(now_error);

                    if(now_error > prev_error && levenberg < 1.e7)
                        levenberg = levenberg * levenberg_factor_mult;
                        disp('now_error > prev_error');
                    else
                        levenberg = levenberg/levenberg_factor_div;
                        if(now_error < best_error)
                            bestBeta_0 = beta(1);
                            bestBeta_1 = beta(2);
                            bestBeta_2 = beta(3);
                            bestlambda_vs_t_teo2 = lambda_vs_t_teo2;
                            best_error = now_error;
                        end
                    end

                    prev_error = now_error;

                    fprintf("\nlevenberg: %s; Fitting error: %s; Beta: %s\n", ...
                            num2str(levenberg), num2str(now_error), ...
                            num2str(beta));

                end
            end
        end


        % Function to obtain the relation of pressure to number of blocked
        % channel in a 12 constriction device
        function delta_P = delta_pressure_12(~,x,rel)
            delta_P = (0.5056*x^6 - 15.468*x^5 + 184.86*x^4 - 1076.9*x^3 + ...
                      3160.9*x^2 - 4062.7*x + 3519.2)*rel;
        end

        
        % Function to obtain the differential pressure of a 1 constriction
        % device. Basic and provisional 
        function [delta_P] = delta_pressure_0(obj)

            % Obtaining all parameters needed
            h = obj.Canal.wch*2/2.95;    % Height of the channel
            w = 20;   % Width of the channel
            Q = obj.Canal.Dispositivo.Q;    % Flow of the experiment
            L = 40;   % Length of the constriction
            eta = 0.001;    % Viscosity

            h_w = h/w;

            Rh = 0;

            if(h_w < 0.2)
                Rh = obj.two_plates_calc(eta,L,h,w);
            elseif(h_w < 1)
                Rh = obj.rectangular_calc(eta,L,h,w);
            elseif(h_w == 1)
                Rh = obj.square_calc(eta,L,h,w);
            end
            
            % Conversion of Pa*um^-3/s to Pa*ml/h
            Rh = Rh*1000000000000/3600;
            
            % If Rh is greater than 0, the differential pressure is
            % calculated
            
            if(Rh > 0)
                delta_P = Rh*Q;
            else
                delta_P = 0;
            end
 
        end


        % Function to obtain the differential pressure of a 1 constriction
        % device
        function [delta_P] = delta_pressure_1(obj)
            
            % Obtaining all parameters needed
            h = obj.Canal.wch*2/2.925;    % Height of the channel
            w = 20;   % Width of the channel
            Q = obj.Canal.Dispositivo.Q;    % Flow of the experiment
            L = 40;   % Length of the constriction
            eta = 0.001;    % Viscosity
            D_cell = obj.Rc_real*2/2.925; % Cell diameter
            
            % Obtaining the width not covered by the cell
            w1 = w-D_cell;

            if(w1 <= 0)
                warning("Value w1 is <= 0.")
                return
            end
            
            % Obtaining all non-dimensionless Als to take them from the
            % start of the constriction (have to see if it is okay)
            Als = nonzeros(obj.vector_Al(:,1));
            Ls = Als + obj.Canal.linea - obj.Canal.posicion_constriccion(1);

            idx = find(diff(Ls) < 0,1);

            if(~isempty(idx))
                Ls = Ls(1:idx);
            end
            
            % Calculating the mean of the Ls that the cell ocuppies in the
            % channel to obtain L1 (for cell) and L2 (rest of the 
            % constriction)
            L0 = Ls(1);
            L1 = mean(Ls);
            L2 = (L - L1);
            
            % Relation of h and w of the cell 
            h_w_1 = h/w1;
            Rh = 0;
            
            % Calculation of the hydraulic resistance depending on the
            % condition
            if(h_w_1 < 0.2)
                Rh = obj.two_plates_calc(eta,L1,h,w1);
            elseif(h_w_1 < 1)
                Rh = obj.rectangular_calc(eta,L1,h,w1);
            elseif(h_w_1 == 1)
                Rh = obj.square_calc(eta,L1,h,w1);
            end
            
            % Relation of h and w of the rest of the channel
            h_w_2 = h/w;
            
            % Again, calculation of the hydraulic resistance depending on 
            % the condition
            if(h_w_2 < 0.2)
                Rh = Rh + obj.two_plates_calc(eta,L2,h,w);
            elseif(h_w_2 < 1)
                Rh = Rh + obj.rectangular_calc(eta,L2,h,w);
            elseif(h_w_2 == 1)
                Rh = Rh + obj.square_calc(eta,L2,h,w);
            end
            
            % Conversion of Pa*um^-3/s to Pa*ml/h
            Rh = Rh*1000000000000/3600;
            
            % If Rh is greater than 0, the differential pressure is
            % calculated
            
            if(Rh > 0)
                delta_P = Rh*Q;
            else
                delta_P = 0;
            end

        end

        
        % Function to calculate the hydraulic resistance if two plates
        % condition
        function Rh = two_plates_calc(~,eta,L,h,w)
            Rh = (12*eta*L)/(w*h^3);
        end   


        % Function to calculate the hydraulic resistance if rectangular
        % condition
        function Rh = rectangular_calc(~,eta,L,h,w)
            Rh = (8*eta*L*(1+h/w)^4)/(h^4*pi);
        end   


        % Function to calculate the hydraulic resistance if square
        % condition
        function Rh = square_calc(~,eta,L,h)
            Rh = (28.4*eta*L)/(h^4);
        end  

    end

end


% Function to do the exact Python interpolation function "interpd1"
% with parameter "fill_value = 'extrapolate'"
function fq = interpd1(x,y,xq)

    x = x(:);
    y = y(:);
    xq = xq(:);

    valid_idx = isfinite(x) & isfinite(y);
    x_valid = x(valid_idx);
    y_valid = y(valid_idx);

    if(numel(x_valid) < 2)
        fq = NaN(size(xq));
        return
    end

    [x_sorted, sort_idx] = sort(x_valid);
    y_sorted = y_valid(sort_idx);

    [x_unique, ~, ic] = unique(x_sorted, 'stable');
    y_grouped = accumarray(ic, y_sorted, [], @mean);

    fq = interp1(x_unique, y_grouped, xq, "linear", "extrap");

    if(~isreal(xq))
        fq = fq + 1i*imag(xq);
    end

    fq(~isfinite(xq)) = NaN;
    
end
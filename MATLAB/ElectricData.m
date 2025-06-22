% CLASS "ElectricData" - Defined to store all electric data related to a
% test, frequencies, times, impedance, among others.

classdef ElectricData < handle

    properties 
        csv_path            % String - Unified CSV file path
        txt_path            % String - TXT file with time and/or date
        t_start             % Double - Time EPOCH when the test was started
        date_start          % Datetime - Date when the test was started
        frequencies         % Double - Frequencies used in the test
        times               % Double - Time of each measure
        impedance           % Array Double - Impedance measures, each 
                            % column corresponds to a different frequency
        phase               % Array Double - Phase measures, each column
                            % corresponds to a different frequency
        medium_imp          % Array Double - Medium impedance mean and std
        medium_phase        % Array Double - Medium phase mean and std
        C1                  % Double - Cell medium capacitance
        R1                  % Double - Cell medium resistance
        Device              % Dispositivo - References the device from 
                            % which the electric data comes from
    end


    methods
        
        % Constructor of the class
        function obj = ElectricData(csv_path,t_start_path,Device)
            
            obj.csv_path = csv_path;
            obj.txt_path = t_start_path;
            [new_t_start, new_date_start] = time_and_date(obj);
            obj.t_start = new_t_start;
            obj.date_start = new_date_start;
            [obj.frequencies,obj.times,obj.impedance,obj.phase] = ...
                load_csv(obj);
            obj.Device = Device;
           
            obj.times = obj.times + obj.t_start - obj.Device.t_start;

            % Sum of a hour due to the time change in Spain
            obj.times = obj.times + 3600*1000;

            obj.C1 = 0;
            obj.R1 = 0;
            
        end
        

        % Function to load data from the t_start TXT file with the time
        % and/or date
        function [t_start,date_start] = time_and_date(obj)
            
            txt_file = obj.txt_path;

            % Checking the file exists
            if ~isfile(txt_file)
                error("The file does not exist.");
            end

            % Getting the lines in the file and counting how many there are
            lines = readlines(txt_file);
            num_lines = numel(lines);

            % Checking is the file has at least one line
            if num_lines < 1
                error("The file does not contain data.")
            end
            
            % Asigning the start time if it is not a NaN
            t_start = str2double(lines(1));
            if isnan(t_start)
                error("The variable t_start is NaN.")
            end

            % If there are not 2 lines, then the date have not been saved,
            % so it is not loaded
            if(num_lines == 2)
                date_start = lines(2);
            else
                date_start = "Data previous to 'impedance_analyzer_v4.js'";
            end

        end  


        % Function to load all data from the unified CSV file
        function [frequencies,times,impedance,phase] = load_csv(obj)

            csv_file = obj.csv_path;
            t_file = obj.t_start;
    
            % Checking if the file could not be opened
            fid = fopen(csv_file, 'r');
            if fid == -1
                error("Couldn't open the file: %s", csv_file);
            end
            
            % Reading new lines till no '#' appears
            line = fgetl(fid);
            while startsWith(line, '#')
                line = fgetl(fid);
            end
            
            % Gets the headers of all data inside the CSV file (divided by
            % ',')
            headers = strsplit(line, ',');
            
            % Reads the rest of the data, saves it and closes the file
            data = textscan(fid, repmat('%f', 1, length(headers)), ...
                            'Delimiter', ',', 'HeaderLines', 0);                            
            fclose(fid);
            
            % Substracting the starting time from the time of all measures 
            % and adding the time difference between clicks
            if(data{1,3}(1) > t_file)
                data{1,3} = data{1,3} - t_file;
            end
            
            % Counting the number of frequencies the test has
            num_freq = length(unique(data{1,4}));
            
            % All frequencies in an array plus names for future figures
            freq = zeros(1,num_freq);
            for i = 1:num_freq
                freq(i) = round(data{1,4}(i));
            end
            
            time_data = zeros((length(data{1,1})/num_freq),1);
            imp_data = zeros((length(data{1,1})/num_freq),2);
            phase_data = zeros((length(data{1,1})/num_freq),2);
            
            % Saving the data from each frequency to a structure with names 
            % like 'data_(frequency)', then representing for each frequency
            for i = 1:num_freq
                    
                counter = 1;
                
                % Saving for each frequency the Test, T_start, T_end, Z and
                % Phase
                for j = i:num_freq:length(data{1,1})
                    time_data(counter) = data{1,3}(j);
                    imp_data(counter,i) = data{1,6}(j);
                    phase_data(counter,i) = data{1,5}(j);
            
                    counter = counter + 1;
                end
            end
            
            % Asigning the data to output
            frequencies = freq;
            times = time_data;
            impedance = imp_data;
            phase = phase_data;
        end


        % Function to obtain all electric data not corresponding to events,
        % corresponding to cells, and event that are not cells
        function [medium_data,no_cell_data,cell_data] = event_data(obj)
            
            % Preparing the time axis and logical array for electric data
            time_axis = obj.times;
            no_event = ones(length(time_axis),1);
            
            % Preparing the time axis and logical array for images
            num_img = length(obj.Device.direc_imagenes);
            time_res = 1/obj.Device.fs*1000;
            time_img = (0:time_res:(num_img-1)*10);
            cell_event = zeros(length(time_img),1);
            no_cell_event = zeros(length(time_img),1);
            
            % Saving the events of the device
            events = obj.Device.Canales.Eventos;
            num_events = length(events);
            
            % For every event we safe the frames and if it is a cell its
            % value in the plot will be 1 and if not it will be 0.5
            for i = 1:num_events
                
                event = events(i);
            
                starting_frame = event.frame_inicial/obj.Device.fs*1000;
                final_frame = event.frame_final/obj.Device.fs*1000;
                
                idx_electric = ((time_axis >= starting_frame) & ...
                               (time_axis <= final_frame));
                idx_images = ((time_img >= starting_frame) & ...
                             (time_img <= final_frame));
            
                if event.is_cell
                    cell_event(idx_images) = 1;
                else
                    no_cell_event(idx_images) = 1;
                end
                
                no_event(idx_electric) = 0;
            end
            
            % The output are all indexes for each type of data
            medium_data = logical(no_event);
            no_cell_data = logical(no_cell_event);
            cell_data = logical(cell_event);
        end
        

        % Function to get mean of impedance and phase for cell medium
        function [mean_std_imp, mean_std_phase] = get_medium_means(obj)
            
            % Getting all data where an event happens
            [medium_data,~,~] = obj.event_data();
            
            % Declaring the output variables
            mean_std_imp = zeros(numel(obj.frequencies),2);
            mean_std_phase = zeros(numel(obj.frequencies),2);
            
            % For each frequency a figure is created
            for i = 1:numel(obj.frequencies)
                
                % Obtaining mean and std of medium data to know where the
                % peaks are
                imp_medium = obj.impedance(:,i);
                phase_medium = obj.phase(:,i);

                imp_medium = imp_medium(medium_data);
                phase_medium = phase_medium(medium_data);

                imp_mean = mean(imp_medium);
                phase_mean = mean(phase_medium);
                imp_std = std(imp_medium);
                phase_std = std(phase_medium);

                % Obtaining the idx for data that exceeds thresholds
                imp_idx = ((obj.impedance(:,i) > (imp_mean - ...
                          3*imp_std)) & (obj.impedance(:,i) < ...
                          (imp_mean + 3*imp_std)));

                phase_idx = ((obj.phase(:,i) > (phase_mean - ...
                            3*phase_std)) & (obj.phase(:,i) < ...
                            (phase_mean + 3*phase_std)));
    
                imp_no_peaks = obj.impedance(:,i);
                phase_no_peaks = obj.phase(:,i);
                
                imp_no_peaks = imp_no_peaks(imp_idx);
                phase_no_peaks = phase_no_peaks(phase_idx);

                imp_mean = mean(imp_no_peaks);
                phase_mean = mean(phase_no_peaks);
                imp_std = std(imp_no_peaks);
                phase_std = std(phase_no_peaks);
                   
                mean_std_imp(i,1) = imp_mean;
                mean_std_imp(i,2) = imp_std;
                mean_std_phase(i,1) = phase_mean;
                mean_std_phase(i,2) = phase_std;
            end

            obj.medium_imp = mean_std_imp;
            obj.medium_phase = mean_std_phase;
        end    

        
        % Function to set C1 and R1
        % Non-linear fitting using the Levemberg-Marquardt algorithm
        % Created Decemeber 2024, I. Fraga, D. Ñaña, B. González, G. Plaza
        % Last Edited April 2025 
        % Fitting of electrical parameters
        function [] = set_C1_R1(obj,num)

            % Getting the data where no cells are passing through (needs
            % modifying to obtain that data)
            data_freq = obj.frequencies;
            [data_imp, data_arg] = get_medium_means(obj);

            freq_Imp = data_freq;
            argument_Imp = data_arg(:,1); % argument of impedance in deg
            modulus_Imp = data_imp(:,1); % modulus of impedance in Ohm
            
            if(num == 1)
                plot(freq_Imp,modulus_Imp, '*');
            end

            % Extracting the vectors of interest from the matrix
            Real_Exp = modulus_Imp.*cos(argument_Imp*pi/180);
            Imag_Exp = modulus_Imp.*sin(argument_Imp*pi/180);

            num_parameters = 2;
            cal_R1 = 1.e5;
            inv_C1 = 1.e8;
            Beta = [cal_R1/2 inv_C1/2]; % Beta(1) = R1; Beta(2) = 1/C1; 

            num_datos = length(freq_Imp);
            num_iter = 5000;
            levenberg = 0.1; % Levenberg-Marquardt parameter
            All_Exp = [Real_Exp(:); Imag_Exp(:)];
            [Real_Teo, Imag_Teo] = obj.Impedance_theor(Beta(1),Beta(2), ...
                                   freq_Imp);
            All_Teo = [Real_Teo(:); Imag_Teo(:)];
            DifferenceAll = All_Teo - All_Exp;
            % new_error = dot(DifferenceAll,DifferenceAll)/dot(All_Exp,All_Exp);

            for j = 1:num_iter

                epsilon_beta = Beta/100;
                J = zeros(2*num_datos,num_parameters);

                for k = 1:num_parameters
	                for i = 1:num_datos
                        [Real_2, Imag_2] = obj.Impedance_theor(Beta(1) + obj.delta(1,k)*epsilon_beta(1), ...
                                                               Beta(2) + obj.delta(2,k)*epsilon_beta(2), ...
                                                               freq_Imp(i));
                        [Real_1, Imag_1] = obj.Impedance_theor(Beta(1), ...
                                                               Beta(2), ...
                                                               freq_Imp(i));
		                J(i,k) = (Real_2 - Real_1)/epsilon_beta(k);
                        J(i+num_datos,k) = (Imag_2-Imag_1)/epsilon_beta(k);
	                end
                end

                deltaBeta = transpose(-inv(transpose(J)*J + ...
                            levenberg*eye(num_parameters))* ...
                            (transpose(J)*DifferenceAll));

                Beta = Beta + deltaBeta;

                [Real_Teo, Imag_Teo] = obj.Impedance_theor(Beta(1), ...
                                       Beta(2),freq_Imp);
                All_Teo = [Real_Teo(:); Imag_Teo(:)];
                DifferenceAll = All_Teo - All_Exp;
                new_error = dot(DifferenceAll,DifferenceAll)/ ...
                            dot(All_Exp,All_Exp);

                if(num == 1)
                    fprintf('Iteración %d - Error: %f\n', j, new_error);
                    fprintf('Beta(1): %f, Beta(2): %f\n', Beta(1), Beta(2));
                end
            end
            
            if(num == 1)
                % Plot of the results
                freq_Simula = min(freq_Imp):20:max(freq_Imp);
                [Real_Simula, Imag_Simula] = Impedance_theor(Beta(1),Beta(2),freq_Simula);
                figure;
                plot(freq_Simula,Real_Simula,'-');
                hold on;
                plot(freq_Imp,Real_Exp,'*');
                figure;
                plot(freq_Simula,Imag_Simula,'-');
                hold on;
                plot(freq_Imp,Imag_Exp,'*');
            end

            obj.R1 = Beta(1);
            obj.C1 = 1/Beta(2);
        end  


        % Function to get impedance as a function of the five parameters
        function [Real_Imp, Imag_Imp] = Impedance_theor(~, R1, inv_C1, ...
                                        freq_Imp)
           % Freq_Imp: frequency in Hz
           % R1, resistence in Ohm
           % C1 capacitance in F
           cal_C1 = 1/inv_C1;
           w = 2*pi*freq_Imp;
           Real_Imp = R1*ones(size(freq_Imp));
           Imag_Imp = -1./(w*cal_C1);
        end 


        % Delta function to be used in obtainment of C1 and R1
        function g = delta(~, a, b)
            if (a == b)
                g = 1;
            else
                g = 0;
            end
        end


        % Function to asign electric data to each event with is_cell == 1
        function [] = asign_electric_data(obj)
            
            % Preparing the time axis for electric data
            time_axis = obj.times;
            
            % Saving the events of the device
            events = obj.Device.Canales.Eventos;
            
            % Obtaining the events with cells
            idx = find([events.is_cell] == 1);
            
            % If no event is a cell, a warning is displayed and function is
            % cancelled
            if(isempty(idx))
                warning("No cell events found...")
                return
            end
            
            % Only the events with cells are analyzed
            num_events = length(idx);
            
            % For every event the electric data between the starting and
            % ending time of the event is asigned to the event in the
            % properties: "imp" and "arg"
            for i = 1:num_events
                
                event = events(idx(i));
            
                starting_frame = event.frame_inicial/obj.Device.fs*1000;
                final_frame = event.frame_final/obj.Device.fs*1000;

                idx_electric = ((time_axis >= starting_frame) & (time_axis <= final_frame));
                
                times_data = time_axis(idx_electric,:) - starting_frame;
                impedance_data = obj.impedance(idx_electric,:);
                phase_data = obj.phase(idx_electric,:);
                
                is_extrema = zeros(size(impedance_data,2),1);

                for j = 1:size(impedance_data,2)
                    
                    medium_imp_mean = obj.medium_imp(j,1);
                    medium_imp_std = obj.medium_imp(j,2);

                    medium_phase_mean = obj.medium_phase(j,1);
                    medium_phase_std = obj.medium_phase(j,2);
                    
                    imp = impedance_data(:,j);
                    arg = phase_data(:,j);

                    is_imp_extrema = any(abs(imp-medium_imp_mean) > 3*medium_imp_std);
                    is_arg_extrema = any(abs(arg-medium_phase_mean) > 3*medium_phase_std);

                    if(any(is_imp_extrema == 1) || any(is_arg_extrema == 1))
                        is_extrema(j,1) = 1;
                    end

                end

                if(~any(is_extrema == 1))
                    warning("All electric data for this event is " + ...
                            "similar to medium data, passing to " + ...
                            "the next one...")
                    continue
                end

                obj.Device.Canales.Eventos(idx(i)).data_electric = ...
                [times_data, impedance_data, phase_data];
            end
        end    


        % Function to display only the raw electric data
        function [] = display_electric_data(obj)

            % For each frequency a figure is created
            for i = 1:numel(obj.frequencies)
                
                % Creation of the time vector of seconds
                time = obj.times/1000;               
    
                % Creating a figure with Z and Phase for each frequency
                frequency = obj.frequencies(i);
        
                % Plots are displayed
                figure(NumberTitle="off")
                
                % Subplot for impedance modulus
                subplot(2,1,1)
                plot(time,obj.impedance(:,i)/1000000);
                title("Impedance Modulus at " + frequency + " Hz"); 
                xlabel("Time, t (s)"); 
                ylabel("Impedance Modulus, |Z| (M\Omega)"); 
                xlim([(min(time)-1),(max(time)+1)])
                ylim([(min(obj.impedance(:,i)/1000000)-0.1), ...
                    (max(obj.impedance(:,i)/1000000)+0.1)]);
                
                % Subplot for phase
                subplot(2,1,2)
                plot(time,obj.phase(:,i)); 
                title("Phase at " + frequency + " Hz")
                xlabel("Time, t (s)"); 
                ylabel("Phase, \Phi (°)"); 
                xlim([(min(time)-1),(max(time)+1)])
                ylim([(min(obj.phase(:,i))-0.2), ...
                    (max(obj.phase(:,i))+0.2)]);
            end

        end    
        

        % Function to display all information (electrical and events) in a 
        % time axis to visualize which are synchronized, the medium is also
        % displayed with mean and std
        function [] = all_display(obj)
            
            % Remove previous figures
            close all;

            % Preparing the time axis and logical array for electric data
            time_axis = obj.times;           
            
            % Preparing the time axis and logical array for images
            num_img = length(obj.Device.direc_imagenes);
            time_res = 1/obj.Device.fs*1000;
            time_img = (0:time_res:(num_img-1)*10)/1000;
            cell_events = zeros(length(time_img),1);
            no_cell_events = zeros(length(time_img),1);
            
            % All indexes are obtained for data
            [medium_data,no_cell_data,cell_data] = obj.event_data();
            
            % All event data is binarized
            cell_events(cell_data) = 1;
            no_cell_events(no_cell_data) = 1;

            % For each frequency a figure is created
            for i = 1:numel(obj.frequencies)
                
                % Obtaining mean and std of medium data to know where the
                % peaks are
                imp_medium = obj.impedance(:,i);
                phase_medium = obj.phase(:,i);

                imp_medium = imp_medium(medium_data);
                phase_medium = phase_medium(medium_data);

                imp_mean = mean(imp_medium);
                phase_mean = mean(phase_medium);
                imp_std = std(imp_medium);
                phase_std = std(phase_medium);

                % Obtaining the idx for data that exceeds thresholds
                imp_idx = ((obj.impedance(:,i) > (imp_mean - ...
                          3*imp_std)) & (obj.impedance(:,i) < ...
                          (imp_mean + 3*imp_std)));

                phase_idx = ((obj.phase(:,i) > (phase_mean - ...
                            3*phase_std)) & (obj.phase(:,i) < ...
                            (phase_mean + 3*phase_std)));

                imp_peaks = zeros(length(time_axis),1);
                phase_peaks = zeros(length(time_axis),1);

                imp_peaks(imp_idx) = 1;
                phase_peaks(phase_idx) = 1;

                % Creating a figure with Z and Phase for each frequency
                frequency = obj.frequencies(i);

                % Plots are displayed
                figure(NumberTitle="off")
                
                % Subplot for impedance modulus
                subplot(2,1,1)                
                plot(time_axis/1000,imp_peaks,'b'); 
                hold on
                plot(time_img,no_cell_events,'r');

                if any(cell_events)
                    plot(time_img,cell_events,'g');
                    legend(['Impedance Peaks ($|Z|_i>\overline{|Z|}+' ...
                          '3\sigma_{|Z|}$)'], 'No cell events', ...
                          'Cell events', 'Interpreter', 'latex')
                else
                    legend(['Impedance Peaks ($|Z|_i>\overline{|Z|}+' ...
                          '3\sigma_{|Z|}$)'], 'No cell events', ...
                          'Interpreter', 'latex')
                end

                hold off
                title("Comparison with impedance at " + frequency + " Hz"); 
                xlabel("Time, t (s)");                 
                xlim([(min(time_axis/1000)-1),(max(time_axis/1000)+1)])
                ylim([0, 1.5]);           
                
                % Subplot for phase
                subplot(2,1,2)
                plot(time_axis/1000,phase_peaks,'b'); 
                hold on
                plot(time_img,no_cell_events,'r');

                if any(cell_events)
                    plot(time_img,cell_events,'g');
                    legend(['Phase Peaks ($\phi_i>\overline{\phi}+' ...
                           '3\sigma_{\phi}$)'], 'No cell events', ...
                           'Cell events', 'Interpreter', 'latex')
                else
                    legend(['Phase Peaks ($\phi_i>\overline{\phi}+' ...
                           '3\sigma_{\phi}$)'], 'No cell events', ...
                           'Interpreter', 'latex')
                end    

                hold off
                title("Comparison with phase at " + frequency + " Hz")
                xlabel("Time, t (s)"); 
                xlim([(min(time_axis/1000)-1),(max(time_axis/1000)+1)])
                ylim([0, 1.5]);

                fprintf("\nThreshold at %d Hz for impedance = %.4f + %.4f \x03A9\n",frequency,imp_mean,3*imp_std)
                fprintf("Threshold at %d Hz for phase = %.4f + %.4f°\n",frequency,phase_mean,3*phase_std)

                imp_no_peaks = obj.impedance(:,i);
                phase_no_peaks = obj.phase(:,i);
                
                imp_no_peaks = imp_no_peaks(imp_idx);
                phase_no_peaks = phase_no_peaks(phase_idx);

                imp_mean = mean(imp_no_peaks);
                phase_mean = mean(phase_no_peaks);
                imp_std = std(imp_no_peaks);
                phase_std = std(phase_no_peaks);
            
                % Plots for mean are displayed
                figure(NumberTitle="off")

                % Subplot for impedance modulus
                subplot(2,1,1)                
                plot(time_axis(imp_idx)/1000,imp_no_peaks/1000000,'b'); 
                hold on
                yline(imp_mean/1000000,"Color","r");
                yline((imp_mean + imp_std)/1000000,"Color","k");
                yline((imp_mean - imp_std)/1000000,"Color","k");
                hold off
                title("Only medium data of impedance at " + frequency + " Hz"); 
                xlabel("Time, t (s)");    
                ylabel("Impedance Modulus, |Z| (M\Omega)")
                xlim([(min(time_axis(imp_idx)/1000)-1),(max(time_axis(imp_idx)/1000)+1)])
                legend(["Medium Impedance", "Mean of Medium Impedance", ...
                        "\pmSTD of Medium Impedance"])

                % Subplot for phase
                subplot(2,1,2)
                plot(time_axis(phase_idx)/1000,phase_no_peaks,'b');
                hold on
                yline(phase_mean,"Color","r");
                yline((phase_mean + phase_std),"Color","k");
                yline((phase_mean - phase_std),"Color","k");
                hold off
                title("Only medium data of phase at " + frequency + " Hz")
                xlabel("Time, t (s)"); 
                ylabel("Phase, \Phi (°)")
                xlim([(min(time_axis(phase_idx)/1000)-1),(max(time_axis(phase_idx)/1000)+1)])
                legend(["Medium Phase", "Mean of Medium Phase", ...
                        "\pmSTD of Medium Phase"])

                fprintf("\nMean and std at %d Hz for medium impedance:\nMean = %.4f \nSTD = %.4f \x03A9\n",frequency,imp_mean,imp_std)
                fprintf("\nMean and std at %d Hz for medium phase:\nMean = %.4f \nSTD = %.4f°\n",frequency,phase_mean,phase_std)
            end
        end

    end
   
end
% Function to process all CSV files in a folder, these CSV files follow the
% format "A_B_C.csv", where A is the number of measurement, B is the start
% time of said measurement (ms) and C is the end time (ms)

function [file_created] = process_csv(electric_path)
    
    % % The user selects the directory of the CSV files
    % electric_path = uigetdir(path,['Select the directory of the ' ...
    %                 'electric data...']);
    
    % Loading all the files in a list
    files = dir(fullfile(electric_path, '*.csv'));
    
    % If the list is empty the function ends because there are no CSV files
    % in the folder selected
    if isempty(files)
        disp('No CSV files found.');
        return;
    end
    
    % Saving the name of the folder to eliminate any file that is not the ones
    % to analyze
    [~, folder_name] = fileparts(electric_path);
    file_del = startsWith({files.name}, folder_name);
    files(file_del) = [];
    
    % Saving every file name into a new variable 
    files_names = {files.name};
    
    % Now, a function is applied to the files where:
    % 1. regexp() searchs for a pattern in the files names provided where:
    %   - "^" to start from the start of the name
    %   - "\d+" to search for consecutive digits
    %   - "match" to return the patterns
    %   - "once" to only return the first for each file name, A, the number
    %     of the measure
    % 2. str2double() converts the first part of each file name to a number
    % so that it can later be sorted
    % 3. cellfun(@x),...,files_names) applies the previous two steps to
    % every element in "files_names, saving all the numbers 
    numbers = cellfun(@(x) str2double(regexp(x, '^\d+', 'match', 'once')), files_names);
    
    % The numbers are sorted and their index is saved to then sort the CSV
    % files
    [~, idx] = sort(numbers);
    sorted_files = files_names(idx);
    
    % Once files have been sorted, reads the information that comes in the 
    % first of sorted CSV files (lines with "#" start that are common in 
    % all files)
    first_file = fullfile(electric_path, sorted_files{1});
    fid = fopen(first_file, 'r');
    info = {};
    
    % Reads till the end of the file
    while ~feof(fid)
        line = fgetl(fid);
        if startsWith(line, '#')
            info{end+1} = line;
        else
            break;
        end
    end
    
    fclose(fid);
    
    % The output file name is created following the name "folder_name" +  
    % "_" + "date", where folder_name will be test1, test2, etc. and date 
    % is the date in which data is being processed (CAMBIAR A LA FECHA DEL TXT)
    date = datetime("now", "TimeZone", "local", "Format", "d-MMM-y");
    new_file = folder_name + "_" + string(date) + ".csv";
    output_file = fullfile(electric_path, new_file);
    
    % The output file is created and the info data is poured and headers of
    % each magnitude too
    file_out = fopen(output_file,"w");
    
    for i = 1:length(info)
        fprintf(file_out, '%s\n', info{i});
    end
    
    headers = {'Test', 'T_start (ms)', 'T_end (ms)', 'Frequency (Hz)', ...
               'Trace th (deg)', 'Trace |Z| (Ohm)'};
    fprintf(file_out, '%s\n', strjoin(headers, ','));
    
    % The extraction of data is done for each sorted file
    for i = 1:length(sorted_files)
        
        % From each file with regexp A, B and C are extracted
        file_path = fullfile(electric_path, sorted_files{i});
        name_parts = regexp(sorted_files{i}, '(\d+)_([\d]+)_([\d]+)', ...
                     'tokens');
        
        % If the "A_B_C.csv" format is not followed, the function stops
        if isempty(name_parts)
            disp(['File ', num2str(i), ' does not follow the format', ...
                 ' "A_B_C.csv".']);
            fclose(file_out);
            return;
        end
    
        % For each, A, B and C, are saved in double type
        name_parts = name_parts{1};
        test = str2double(name_parts{1});
        t_start = str2double(name_parts{2});
        t_end = str2double(name_parts{3});
        
        % Information in the file is read as a table
        data = readtable(file_path, 'FileType', 'text', 'Delimiter', ...
               ',', 'ReadVariableNames', false);
        
        % All files that starts with # or is part from the headers is not
        % analised
        % data = data(~startsWith(data.Var1, '#') & ~contains(lower(data.Var1), 'frequency'),:);
        
        % All useful information is saved in the output file
        for j = 1:height(data)
            fprintf(file_out, '%d,%d,%d,%.10f,%.10f,%.10f\n', test, ...
                t_start, t_end, data.Var1(j), data.Var2(j), data.Var3(j));
        end  
    end  
    
    % Closing of the output file
    fclose("all");
    fprintf("Unified file created: %s\n", output_file);

    file_created = output_file;
end
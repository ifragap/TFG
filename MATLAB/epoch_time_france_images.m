function [epoch_time] = epoch_time_france_images(directorio)
    
    % directorio="D:\IMAGENES_TRABAJO\NUEVOS_TEST_FRANCIA\test2\_1\Default\img_channel000_position000_time000000000_z000.tif";
    [path,~,~]=fileparts(directorio);
    direc_experimentos = uigetfile(path+"\*.txt","Selecciona el txt de la metadata");
    direc_metadata=path+"\"+string(direc_experimentos);
    fid = fopen(direc_metadata);
    raw = fread(fid,inf);
    str = char(raw');
    fclose(fid);
    val = jsondecode(str);
    val_cell=struct2cell(val);
    metadata_primera_imagen=val_cell{3};
    tiempo_string=metadata_primera_imagen.UserData.TimeReceivedByCore.scalar;
    tiempo_string=tiempo_string(1:23);

    
    % With the pattern we create a final date to convert to date time
    dt = datetime(tiempo_string,'Format','yyyy-MM-dd HH:mm:ss.SSS');
    
    % Date time data is used to convert to Epoch time in seconds
    epoch_sec = posixtime(dt);
    
    % Conversion from seconds to miliseconds
    epoch_time = epoch_sec*1000;
   
end


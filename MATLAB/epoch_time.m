function [epoch_time] = epoch_time(image_path)
    
    % The path is splitted by "\" character to then only take the end, the
    % name of the image file
    paths = split(image_path,"\");
    image_name = paths(end);
    
    % The image name file is splitted by "_" character to then take the
    % end-2 part (which is the date) and end-1 part (which is the time)
    name_parts = split(image_name,"_");
    
    new_date = char(name_parts(end-2));
    new_time = char(name_parts(end-1));
    
    % Every part of the date is turned into a variable to then follow a
    % pattern
    year = new_date(1:4);
    month = new_date(5:6);
    day = new_date(7:8);    
    hour = new_time(1:2);
    min = new_time(3:4);
    sec = new_time(5:6);
    msec = new_time(7:end);
    
    % With the pattern we create a final date to convert to date time
    final_date = sprintf('%s-%s-%s %s:%s:%s.%s',year,month,day,hour,min,sec,msec);
    dt = datetime(final_date,'Format','yyyy-MM-dd HH:mm:ss.SSS');
    
    % Date time data is used to convert to Epoch time in seconds
    epoch_sec = posixtime(dt);
    
    % Conversion from seconds to miliseconds
    epoch_time = epoch_sec*1000;
   
end


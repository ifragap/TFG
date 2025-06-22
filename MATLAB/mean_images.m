% Function that with the path of the folder where images are results in a
% single image that is the mean of all of them.

function [mean_image] = mean_images(images_path,x)
    
    % Loading of images as Strings of their path in a variable
    filePattern = fullfile(images_path, "*."+x);
    images = dir(filePattern);
    %images = natsortfiles(images);
    images_paths(1:length(images)) = "";
    
    % Saving their path in the variable
    for i=1:1:length(images)
        images_paths(1,i) = images_path + "\" + images(i,1).name;
    end
    
    % Making sure the folder is not empty
    if (isempty(images))
        error("No images in the folder.");
    end
    
    % Loading of the first image to start the summatory
    first_image = imread(images_paths(1));
    image_sum = im2double(first_image);
    
    % Continuing the processing of the rest of the images
    for i = 2:length(images_paths)
        
        % Reading the images
        current_image = imread(images_paths(i));
        
        % Checking the images have the size they should
        % if (isequal(size(current_image),[540 720]))
        %     problematic_image = replace(dir_imagenes(1),(imageFolder + ...
        %                                 "\"),"");
        %     error("The image " + problematic_image + " does not have" + ...
        %           " the size [540 720].");
        % end
        
        % Acumulation of images values for the summatory
        image_sum = image_sum + im2double(current_image);
    end
    
    % Calculation the mean image by dividing the number of images to the
    % summatory
    mean_image = image_sum/length(images_paths);
    
    % Changing format to uint8
    mean_image = im2uint8(mean_image);

end
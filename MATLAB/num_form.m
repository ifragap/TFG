function [] = num_form(data,file_name)
    
    file_name = "res_exp_" + file_name + ".txt";

    file_ID = fopen(file_name,"w");

    for i = 1:size(data,1)
        
        if(i == size(data,1))
            salto = "";
        else
            salto = "\n";
        end    

        fprintf(file_ID,("%0.18e %0.18e"+salto),data(i,1),data(i,2));
    end    

    fclose(file_ID);
end


function savefield(filename, format,format_num,bfield,title,debug_flag)
    [fid,message] = fopen(filename,'w');
    if( fid == -1) 
        disp(message);
        stop;
    end
    string = format;num1 = format_num;
    if(title(1) ~= 0)
        for i = 1:4
            fprintf(fid,'%15.8e\n',title(i)); 
        end
        for i = 5:6
            fprintf(fid,'%d\n',title(i));
        end
    end
    m = length(bfield(1,:));n = length(bfield(:,1));
    m1 = floor(m/num1);num2 = mod(m,num1);
    string1 = string;string2 =string;
    for i = 2 : num1
       string1 =[string1,string];
    end
    string1 = [string1,'\n'];
    for i = 2 : num2
       string2 =[string2,string];
    end
    string2 = [string2,'\n'];
    
    if(debug_flag == 1)
        num1
        num2
        string1
        string2
    end
    for i = 1 : n
        k = 0;
        for j = 1 : m1 
            fprintf(fid,string1,bfield(i,k+1:k+num1)); 
            k = k + num1;
        end
        if( num2 > 0)
            fprintf(fid,string2,bfield(i,k+1:k+num2));
        end
    end
    fclose(fid);    
end
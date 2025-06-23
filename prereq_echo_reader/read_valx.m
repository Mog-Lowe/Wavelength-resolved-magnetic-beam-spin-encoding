function [ valx ] = read_valx( linex,par )

if  regexp(linex, par)>0
        vl=strsplit(linex,':');
        if length(vl)==1
            vl=strsplit(linex,' ');
        end
        numstr = vl{2};
        valx = str2num(numstr);
else
    valx=[];
end

end


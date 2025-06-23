function [v] = getvar(filename,var)
ok=1;
if exist(filename,'file')~=2
    display('file not found');
    v=0;
    return
end
fid = fopen(filename,'r');
if fid==-1 
    v='no file';
end
while (1)
    tline = fgetl(fid);   
    if tline==-1
        break
    end
    if regexp(tline, var)>0
        ok=0;
        break
    end
end
if ok==1 
    v='9e99';
    fclose(fid);
    return
    
end
a=tline;
a1=strsplit(a,':');
if length(a1)>1
    v=a1{2};
else
    a1=strsplit(a,' ');
    v=a1{2};
end
fclose(fid);


end

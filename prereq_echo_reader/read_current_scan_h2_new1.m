function [Iscan,counts,datavecX,lam,spec_lamX,err] = read_current_scan_h2_new1(filename,NoPlot,mirror,Energyplot,FixOn) 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% A note to future self that's getting annoyed because this isn't working -
% you probably need to delete the hack on line 435 where the number of
% columns in the file has been hard coded as 6 :-)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

global sht
global numsteps
ext_sig=0;
countsd=[];
% Tnz=100
% mirror ==1 No Mirror  of the signal   mirror =0 mirror the signal
% NoPlot=1 = No plotting
% Energyplot=1 No energy plot (becuase sometime it's meaningless)
% FixOn=1 drift Fix On FixOn=0 Fix Off (No fix)
if nargin<2
    Energyplot=0;
    NoPlot=0;
    mirror=0;
    FixOn=0;
end

if nargin<3
    mirror=0;
    Energyplot=0;
    FixOn=0;
end

if nargin<4
    Energyplot=0;
    FixOn=0;
end

if nargin<5
    FixOn=0;
end

% for H2 NO mirror ir allowed
mirror=1;

BGON=0;

fid = fopen(filename,'r');
if fid==-1
    display('File Not Found')
    ek=0;spec_pos=0;Iscan=0;counts=0;datavecX=0;lam=0;spec_lamX=0;
    return
end


bgg=getvar(filename,'# BG');
if regexp(bgg,'ON')>0
    BGON=1;
end
BGON = 0;
sampletemp=getvar(filename,'crystal temperature');



if FixOn==1
    BGON=0;
end

fl=0;
Only2=0;
clear datavec

callibI1=0;
wI=0;
sht=0;
numpoints=1;
while (1)
    tline = fgetl(fid);   
    if tline==-1
        break
    end
    
    pyy=read_valx(tline, '# solenoid');
    if ~isempty(pyy) ;scan_psu=round(pyy);end    
    
    pyy=read_valx(tline, 'solenoid:');
    if ~isempty(pyy) ;scan_psu=round(pyy);end    
    
    pyy=read_valx(tline, 'manipulator desired alpha:');
    if ~isempty(pyy) ;alpha=pyy;end    
	
    if regexp(tline,'Calibrated current I1')>0
        callibI1=1;
    end
    
    pyy=read_valx(tline, 'I1-I2-A');
    if ~isempty(pyy) ;a1=pyy;end

    pyy=read_valx(tline, 'I1-I2-B');
    if ~isempty(pyy) ;b1=pyy;end
    
    pyy=read_valx(tline, 'other_solenoid_current');
    if ~isempty(pyy) ;const_curr=pyy;end
    
    if regexp(tline,'SHT')
        sht=1;
    end
    pyy=read_valx(tline, 'nozzle temperature');
    if ~isempty(pyy) ;Tnz=pyy;end
    
    
    pyy=read_valx(tline, 'numloops');
    if ~isempty(pyy) ;nloops=pyy;end

    pyy=read_valx(tline, 'numsteps');
    if ~isempty(pyy) ;numsteps=pyy;end
    
    pyy=read_valx(tline, 'min_value');
    if ~isempty(pyy) ;min_value=pyy;end
    
    pyy=read_valx(tline, 'max_value');
    if ~isempty(pyy) ;max_value=pyy;end
    % Identify the start of the data block
    if regexp(tline,'# loop number')>0
        fl=1;
        break
    end
    if regexp(tline,'# angular data')>0
        fl=1;
        break
    end
end
if fl==0;display('no data in file');return;end
if (exist ('Tnz')==0),    Tnz=100,end;
Tnz=round(Tnz);
ndata=numsteps*numpoints;
%datavec=zeros(ndata,2);
counts=zeros(ndata,1);
k=1;
ok=0;
ttll=['I',num2str(3-(scan_psu)),'=',num2str(const_curr),'  \alpha=',num2str(alpha), '   T_s=',num2str(sampletemp)];
[datavec,nloops]=read_file(fid,NoPlot);

%%%%%%%%%%%%%%%%%%%%%% FIX  %%%%%%%%%%%%%%%%%%%%%
k=nloops;
if k==0
    k=1;
end

datavecX=datavec(1:k,:,:);
if FixOn>0
    q=datavec(1:k,:,1);
    vq=reshape(q',[size(q,1)*size(q,2),1]);
    X=1:length(vq);
    p=polyfit(X,vq',3);
    Y=polyval(p,X);
%     if NoPlot==0
%         figure
%         plot(X,Y,X,vq)
%         title('Fix')
%         drawnow
%     end
    q1=reshape(vq./Y',[size(q,2),size(q,1)]);
    datavec(1:k,:,1)=q1';
    
end

if k==1
    counts=datavec(:,:,1);
    counts0=datavec(:,:,1);
    if BGON==1
    bgc=datavec(:,:,5);
    
%         countsm=counts-bgc;
%         countsd=counts./bgc;
    end
     
    [Iscan,counts,err]=getavg(datavec,k,2,1);

        if BGON==1
            datavec(:,:,1)=counts0-bgc;
            [Iscan,countsm]=getavg(datavec,k,2,1);
            datavec(:,:,1)=counts0./bgc;
            [Iscan,countsd]=getavg(datavec,k,2,1);
        end
    if ext_sig==1 
        ext_counts=datavec(1,:,5);
        ext_Iscan=datavec(1,:,3);
    end
    if Only2==1
        Iscan=datavec(1,:,1);
    else
%         Iscan=datavec(1,:,3);
%         [Iscan,counts]=getavg(datavec,k,2,1);
    end
else 
    [Iscan,counts,err]=getavg(datavec,k,2,1);
    if ext_sig==1 ,[ext_Iscan,ext_counts]=getavg(datavec,k,2,5);end
        if BGON==1
            counts0=datavec(1:k,:,1);
            bgc=datavec(1:k,:,5);
            datavec(1:k,:,1)=counts0-bgc;
            [Iscan,countsm]=getavg(datavec,k,2,1);
            datavec(1:k,:,1)=counts0./bgc;
            [Iscan,countsd]=getavg(datavec,k,2,1);
        end
        
    
end
% datavec=datavec(1:k,:,:);

if scan_psu==1
    imm=1/a1*const_curr-b1/a1;
else
    imm=a1*const_curr+b1;
end



if BGON==1
    
    qsig=counts0(1:k,:);
    qbgold=datavec(1:k,:,5);
    datavec(1:k,:,5)=fix_bg(datavec(1:k,:,5));
    
    qbg=datavec(1:k,:,5);
    
    
    qsig1=reshape(qsig',[size(qsig,1)*size(qsig,2),1]);
    qbg1=reshape(qbg',[size(qbg,1)*size(qbg,2),1]);
    
    qbgold1=reshape(qbgold',[size(qbgold,1)*size(qbgold,2),1]);
    
    if NoPlot==0
        figure
        plot(qsig1)
        title(' Long sig')
        figure
        plot(qsig1./qbg1)
        title(' Long sig/bg')
        figure
        plot(qbg1)
        hold on
        plot(qbgold1)
        
        title(' Long bg')
        
    end
end
%  counts=counts-Y;

%%%%%%%%%%%%%%%%%%%%%% FIX  %%%%%%%%%%%%%%%%%%%%%

% fclose(fid);
if NoPlot==0 
%     figure
%     
%     
%     plot(Iscan,counts,'o-')
%     title ('measure')
%     drawnow
    if sht==1
        figure
%         errorbar(Iscan,counts,err,'-o')
        plot(Iscan,counts/mean(counts),'-o')
        title ('measure with error bars')
    end

    if ext_sig==1 
        figure
        plot(ext_Iscan,ext_counts,'o-')
    end
end

%return
[I_full,data_full] = inter_n_mirror(Iscan,counts,imm,1,mirror);
if ext_sig==1 ,[ext_I_full,ext_data_full] = inter_n_mirror(ext_Iscan,ext_counts,imm,1,mirror);end
if BGON==1
    [I_fullm,data_fullm] = inter_n_mirror(Iscan,countsm,imm,1,mirror);
    [I_fulld,data_fulld] = inter_n_mirror(Iscan,countsd,imm,1,mirror);
    if NoPlot==0
    figure
    plot(Iscan,countsm,'o-')
    title('sig-bg')
    figure
    plot(Iscan,countsd,'o-')
    title('sig/bg')
    end
end
%%%%%%%%%%%%%%%%%%%% End FIX %%%%%%%%%%%%%%%%%%%%

[ek,spec_pos,lam,spec_lam]=getspec(I_full,data_full,scan_psu,Tnz);
if sht==1
    

end

if BGON==1
    [ekm,spec_posm,lamm,spec_lamm]=getspec(I_fullm,data_fullm,scan_psu,Tnz);
    [ekd,spec_posd,lamd,spec_lamd]=getspec(I_fulld,data_fulld,scan_psu,Tnz);
end
if BGON==1
    spec_lamX=spec_lamm;
else
    spec_lamX=spec_lam;
end
    

if NoPlot==1
    return
end
if ~isempty(countsd)
    counts=countsd;
end

if ext_sig==0
    figure
    plot(lam,spec_lam,'.-')
    if FixOn==1
        ttll=[ttll,' Fix'];
    end
    title(ttll)
    
        
    if BGON==1
        hold on
        plot(lamm,spec_lamm,'.-')
        drawnow
    end
    
else
    plot(lam,spec_lam/1+0*max(spec_lam),'.-')
    legend('detectoe','Extrel')
    title(ttll)
    drawnow
end
%axis([0,3e-10,-inf,inf])

if Energyplot==0
    if length(ek)>1
        figure
        if ext_sig==0
            plot(ek,spec_pos)
        end
    end
end

figure(2)
plot(Iscan,counts/mean(counts))

return


end
% end




function [I_full,data_full,okk] = inter_n_mirror(I,data,imm,mpoint,mirror)
% if mpoint =0 mirror at point 1

% resample the data with interpolation 

okk=0;
I_full=I-mean(I);
data_full=data;
return


end

function [ek,spec_pos,lambda_pos,spec_pos_lam]= getspec(I_full,data_full,scan_psu,Tnz)



%  display('PARAMETES OF H2')
beff1=0.010957903656743;
beff2=0.010911076157630;

if  scan_psu==1
    beff=beff1;
elseif scan_psu==2
    beff=beff2;
else
    fprint(1,'wrong sol number..');
    return
end    
[lambda_pos,spec_pos_lam]=calc_spectrum(I_full,beff,data_full,'H2',Tnz);
ek=0;
spec_pos=0;
% uncomment this to aloow plots

return
end

function [newx,newy,err]= getavg(dv,ldv,cx,cy)
global sht
if size(dv,3)<3; cx=1;cy=2; end
err=0;

% when addidng two noisey measurements 1 A with  std=sa
% and the other B with std=sb then the resulting signal is A+B 
% and std(A+b)=sqrt(sa^2 + sb^2)
if sht==1
    lps=size(dv,1);
    col_std_sig=7;
	col_std_bg=8;
    if lps>1
        
        bgstd=norm(dv(:,:,col_std_bg),1)/sqrt(lps);
        sigstd=norm(dv(:,:,col_std_sig),1)/sqrt(lps);
        loopstd=std(dv(:,:,1))/sqrt(lps);
        err=sqrt(bgstd.^2+sigstd.^2+loopstd.^2);
    else
        bgstd=dv(1,:,col_std_bg);
        sigstd=dv(1,:,col_std_sig);
        err=sqrt(bgstd.^2+sigstd.^2);
    end
end
newx=dv(1,:,cx);
if size(dv,1)==1
    newy=dv(1,:,cy);
else
    newy=mean(dv(:,:,cy));
end
%newy=newy(nxx);

end


function [datavec,nloops]=read_file(fid,NoPlot)
global numsteps
A = fread(fid,'*char')';
fclose(fid);
b1=regexpi(A,'loop');
b2=regexpi(A,char(10));
ncol=length(strsplit(A(b2(2)+1:b2(3)-1)));
%
%
%
%**************************************************************************
% ncol = 5;
%**************************************************************************
%
%
%
b3=regexpi(A,'finish');
nloops=length(b1);

if isempty(b3)
    nloops=nloops-1;
    if NoPlot==0
        display(['there are ',num2str(nloops), ' loops last loop is ignored'])
    end
else
    b1=[b1,b3];
    if NoPlot==0
        display(['there are ',num2str(nloops), ' loops'])
    end
end
lpposition=[];
k=nloops;

for i=1:nloops
    c3=b2(b2>b1(i));
    lpposition(i)=c3(1);
    txtloop=A(lpposition(i)+1:b1(i+1)-3);
    C = textscan(txtloop,'%f');
    if length(C{1})<numsteps*ncol
        nloops=nloops-1;
        break
    end
    datavec(i,:,:)=reshape(C{1},ncol,length(C{1})/ncol)';
end
if nloops==0
    c3=b2(b2>b1);
    txtloop=A(c3(1)+1:end);
    C = textscan(txtloop,'%f');
    C=C{1};
%     length(C{1})
    nrows=floor(length(C)/ncol);
    C=C(1:ncol*nrows);
    datavec(1,:,:)=reshape(C,ncol,nrows)';
    k=1;
end
% datavec=datavec(:,1:250,:);
end


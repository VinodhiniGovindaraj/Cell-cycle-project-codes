% Paper: Transcriptional fluctuations govern the serum dependent cell cycle duration heterogeneities in Mammalian cells
% Author: Vinodhini Govindaraj, Subrot Sarma, Atharva Karulkar, Rahul Purwar and Sandip Kar
% e-mail about the code: vinodhinigovindaraj@gmail.com,sandipkar@iitb.ac.in
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;

% no of cell lineages and generations
    nk=1;
    nrows=3;
    ku=23;
    
% CV for transcription rates
    cvp=0.2;
    sig=1.0;
    sigd=0.05;
    vf=6;
    flagnoise=1;

    dat=[];
    seed=100;%100,300,600
    rng(seed);  %start from the same seed each time


for jk=1:nk
    
%	initial time,final time and initialization of number of reactions
    td=[];
    zd=[];
    akc=[];
    

    
    z=zeros(45,1);
    idiv=1;
    ij=0;

%	initial values of mrna and proteins
	
	is=1;
	%isc=1;

	%[GFi(idiv),z0(1)]=deal(dGF);
    
	[inMV,zd(1,1)]=deal(5769*is);
	[inMycm,zd(2,1)]=deal(138*is); 
	[inMycp,zd(3,1)]=deal(2*is);
	[inE2Fm,zd(4,1)]=deal(62*is);
	[inY,zd(5,1)]=deal(9392*is);
	[inDE,zd(6,1)]=deal(6101*is);
	[inRbm,zd(7,1)]=deal(480*is);
	[inRbt,zd(8,1)]=deal(10206*is);
	[iniDE,zd(9,1)]=deal(1599*is);
	[iniDEP,zd(10,1)]=deal(1*is);
	[iniRb1,zd(11,1)]=deal(8*is);
	[iniRb2,zd(12,1)]=deal(8538*is);
	[inCycDm,zd(13,1)]=deal(171*is);
	[inCycD,zd(14,1)]=deal(829*is);
	[inCycEm,zd(15,1)]=deal(308*is);
	[inCycE,zd(16,1)]=deal(8882*is);
	[inIm,zd(17,1)]=deal(204*is);
	[inI,zd(18,1)]=deal(524*is);
	[inIp,zd(19,1)]=deal(548*is);
	[inXI,zd(20,1)]=deal(155*is);
	[inSkp2m,zd(21,1)]=deal(205*is);
	[inSkp2T,zd(22,1)]=deal(13615*is);
	[inSkp2a,zd(23,1)]=deal(13507*is);
	[inCycAm,zd(24,1)]=deal(154*is);
	[inCycA,zd(25,1)]=deal(5567*is);
    
	[inC20m,zd(26,1)]=deal(10*is);
	[inCdc20T,zd(27,1)]=deal(622*is);
	[inCdc20a,zd(28,1)]=deal(65*is);
	[inCdh1m,zd(29,1)]=deal(199*is);
	[inCdh1T,zd(30,1)]=deal(10000*is);
	[inCdh1a,zd(31,1)]=deal(27*is);	
	[inIEP,zd(32,1)]=deal(1105*is);	
	[inCdt1m,zd(33,1)]=deal(214*is);
	[inGemm,zd(34,1)]=deal(214*is);
	[inGem,zd(35,1)]=deal(6699*is);
	[inCdt1,zd(36,1)]=deal(8*is);
    
    [inCycBm,zd(37,1)]=deal(445*is);
	[inCycBT,zd(38,1)]=deal(4259*is);
	[inCycBa,zd(39,1)]=deal(174*is);
    
	[inC25m,zd(40,1)]=deal(102*is);
	[inC25T,zd(41,1)]=deal(5007*is);
	[inC25a,zd(42,1)]=deal(231*is); 
    
    [inWee1m,zd(43,1)]=deal(159*is);
	[inWee1T,zd(44,1)]=deal(7999*is);
	[inWee1a,zd(45,1)]=deal(7618*is);
    
  %{  
    [inSkp2i(idiv),z0(1)]=inSkp2T-inSkp2a(idiv);
    [inCdc20i(idiv),z0(1)]=inCdc20T-inCdc20a(idiv);
    [inCdh1i(idiv),z0(1)]=inCdh1T-inCdh1a(idiv);
    [inIEPi(idiv),z0(1)]=(10000*is)-inIEP(idiv);
	[inE2Fp(idiv),z0(1)]=inY-inDE(idiv)-iniDE(idiv)-iniDEP(idiv);	
	[ipRb,z0(1)]=iniRb1(idiv)+iniRb2(idiv);
	[inRb(idiv,z0(1)inRbT-ipRb-iniDE(idiv)-iniDEP(idiv);
	tid(idiv)=0;
%}
%The rate constants 
    
%	initial rates for varying rate constants(without scaling factor)
%	#mycm	
	[bmycn,ark(1)] = deal(0.000033);
	[dkMcm1n,ark(2)] = deal(0.033);
	[dkMcm2n,ark(3)] = deal(0.000015);
%	#E2FM		
	[be2fn,ark(4)] = deal(0.0000033);
	[dkEm1n,ark(5)] = deal(0.005);
	[dkEm2n,ark(6)] = deal(0.0005);
%	#Rbm
	[dkrbmn,ark(7)] = deal(0.05);
%	#CycDm
	[dkCDm1n,ark(8)] = deal(0.000001);
	[dkCDm2n,ark(9)] = deal(0.01);
%	#CycEm
	[bcycemn,ark(10)] = deal(0.000002);
	[dksemn,ark(11)] = deal(0.04);
%	#p27m
	[dksimn,ark(12)] = deal(0.0142);
%	#Skp2m
	[ak17mn,ark(13)] = deal(0.0036);
%	#CycAm
	[bsamn,ark(14)] = deal(0.000001);
	[dksamn,ark(15)] = deal(0.02);
%	#Cdh1m
	[ak3mn,ark(16)] = deal(0.06);
%	#Cdc20m
	[ak5amn,ark(17)] = deal(0.0001);
	[ak5bmn,ark(18)] = deal(0.004);
%	#Cdt
	[ak21mn,ark(19)] = deal(0.0075);
%	#Gem
	[ak20mn,ark(20)] = deal(0.0075);
%	#CycBm
	[dkscbmn,ark(21)] = deal(0.02);
%	#C25m
	[dksc25mn,ark(22)] = deal(0.02);
%	#Wee1m
	[dkswmn,ark(23)] = deal(0.016);
    %{
    akm=ark;
    cf=1;
    param=parameter(akm,cf);
    filename='param.csv';
    xlswrite(filename,param');
    %}
    
%	calculating log values for transcription rates

    tcv=zeros(ku);
    rsigma=zeros(ku);
    rmu=zeros(ku);
    
%	#CV for the transcription rates
    for it=1:ku
		tcv(it)=cvp;
    end

%	putting transcription rates in an array

%	picking the new rates for each daughter pair
    for ir=1:ku

		rsigma(ir)=sqrt(log(1.0+(tcv(ir)*tcv(ir))));
		rmu(ir)=log(ark(ir))-(rsigma(ir)*rsigma(ir))/2.0;
    end


    for jr=1:ku
     	rk1 = rand();
        rk2 = rand();
        api = 4.0*atan(1.0);

        if (rk1<=0.0)
            rk1 = 0.00001;
        end
               
        ba  = -2.0*sig*sig*log(rk1);
        ba1 = cos(2.0*api*rk2);
        sf  = sqrt(ba)*ba1;
        akc(jr,1)=exp(rmu(jr)+rsigma(jr)*sf);
   end   


td(1) = 0.0;
tfinal = 72.0; %set to some very large number
tspan  = [td(1) tfinal];

% using events -- very helpful in computing precisely when cell division occurs.
%options = odeset('RelTol',1e-5,'AbsTol',1e-5);

options = odeset('Events',@events_cellcycle,'RelTol',1e-5,'AbsTol',1e-5);

cc=[];
  
for ii=0:nrows
   
    	if (ii==0)
			cf=1;
		else
			cf=cf+0.4;
	    end 

    for jj=1:(2^ii)
  
    time     = [];
    zsol     = [];
    
    tevent   = [];
    zevent   = [];
    ievent   = [];

    
    ij=ij+1;
    
    

    tstart=td(ij);
    z0=zd(:,ij);
    akm=akc(:,ij);
    param=parameter(akm,cf);
    

    ie=0;
  % skip cycles where simulation time ended  
    if(tstart>=tfinal)
        for it=1:2   
             idiv=idiv+1;
             td=[td;tstart];
             zd(:,idiv)=z0;
             akc(:,idiv)=akm;
        end
     break;
    end
   icount=0;
 % ode solver for each cell
    while tstart<tfinal
        
        % event loop starts
         while (ie ~= 1) & (tstart < tfinal)
            
        
 [t,z,xe,ye,ie]=ode15s(@cellcycle_rhs,[tstart tfinal],z0,options,param);
 
            time     = [time;t];
            zsol     = [zsol;z];
           
            tevent   = [tevent;xe];
            zevent   = [zevent;ye];
            ievent   = [ievent;ie];
           
   
            tstart=t(end);
            z0=zsol(end,:);
               
            
                if (isempty(ie))
                    ie=4;
        
                elseif (ie(end)==1)
                           if (ii>0)
                             
                             tcc =round(t(end)-td(ij),2);

                             tg1=round(tcmax-td(ij),2);
                             tsg2m=tcc-tg1;
                             dat=[dat;[tcc tg1 tsg2m jk ii jj]]; 
                             
                           end
    
                elseif (ie(end)==2)
                            akn=trvariation(akm,ku,vf);
                            akm=akn;
                            param=parameter(akm,cf);
                            [mCdt1, idx]=max(z(:,36));
                            tcmax=t(idx);
                            icount=icount+1;
                            if icount==2
                                t(end)=72;
                                break;
                            end                      
             
               elseif (ie(end)==3)  
                    
                            akn=trvariation(akm,ku,vf);
                            akm=akn;
                            param=parameter(akm,cf);

               end

          end
             

        %event loop ends    
        
% initial conditions of daughter pairs after cell division or end of
% simulation
for id=1:2
    
    
    if (ie(end)==1)
          
         idiv=idiv+1;
    
         td=[td;t(end)];
         akc(:,idiv)=akm;
         
            if flagnoise==0

                f=0.5;
                zd(:,idiv)=round(zsol(end,:)*f);
         
            elseif flagnoise==1 
                if id==1
                f= randn(1)* sigd + 0.5;
                
                zd(:,idiv)=round(zsol(end,:)*f);

                else
                zd(:,idiv)=round(zsol(end,:)*(1-f)); 

                end
                
            end
        
    else
            
          idiv=idiv+1;
          td=[td;t(end)];  
          zd(:,idiv)=zsol(end,:);
         
          akc(:,idiv)=akm;
          
    end
end



 cc =[cc;round(t(end)-td(ij),2)];

 break;
 
    end
  % end of while loop for the cell     
%{ 

  if jk==1 & ii>=0 & jj<=1
    
%plot(time,zsol(:,39),'LineWidth',2,'Color','r');
%hold on;
plot(time,zsol(:,9),'LineWidth',2,'Color','g');
hold on;
ax=plot(time,zsol(:,6),'LineWidth',2,'Color','r');
hold on;
xlabel('Time (h)');
ylabel('Protein level');
%title('active Wee1')
set(gca,'FontSize',35,'FontName','Times New Roman');
set(findobj(gcf,'type','axes'),'XColor','k','YColor','k','LineWidth',3);
figurepalette('show') ;
pbaspect([2 1.5 1])
  end
%}
    [t(end) jk ii jj]

    end
end

D=flip(cc);
B=[1 2;3 4;5 6;7 8;9 10;11 12;13 14];
%B=[1,2]%
%plot(phytree(B,D));



end

fid=fopen('timings.txt','w');

dat


fprintf(fid,'%f %f %f %d %d %d\n',dat');
fclose(fid);



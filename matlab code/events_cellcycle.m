
%     Paper: Exploring the Roles of Noise in the Eukaryotic Cell Cycle
%     (Submitted to PNAS)
%     Author: Sandip Kar, William Baumann, Mark R. Paul and John J. Tyson
%     e-mail about the code: skar@vt.edu
%     DATE: 1st October, 2008
%
function [value,isterminal,direction] = events_cellcycle(t,z,param)
    
    
    %	initialization for period calculation
	maxc=0;
	maxg=0;
	tc=0;
	ncell=1;
	eg(ncell)=maxc;
	teg(ncell)=tc;
	mcount=1;
	ncount=0;
	  isc=10;      
    
	MV=z(1);
	Mycm=z(2); 
    Mycp=z(3);
	E2Fm=z(4);
	Y=z(5);
	DE=z(6);
	Rbm=z(7);
	RbT=z(8);
	iDE=z(9);
	iDEP=z(10);
	iRb1=z(11);
	iRb2=z(12);
	CycDm=z(13);
	CycD=z(14);
	CycEm=z(15);
	CycE=z(16);
	Im=z(17);
	I=z(18);
	Ip=z(19);
	EI=z(20);
	Skp2m=z(21);
	Skp2T=z(22);
	Skp2=z(23);
	CycAm=z(24);
	CycA=z(25);
	Cdc20m=z(26);
	Cdc20T=z(27);
	Cdc20A=z(28);
	Cdh1m=z(29);
	Cdh1T=z(30);
	Cdh1=z(31);
	IEP=z(32);
	Cdt1m=z(33);
	Gemm=z(34);
	Gem=z(35);
	Cdt1=z(36);
     %
    CycBm=z(37);
	CycBT=z(38);
	CycBa=z(39);
    %
	C25m=z(40);
	C25T=z(41);
    C25a=z(42); 
    %
    Wee1m=z(43);
	Wee1T=z(44);
	Wee1a=z(45);
    

    
               value=[CycBa-1000;CycA-2000;CycBa-4000];
               isterminal=[1;1;1];
               direction=[-1;1;1]; 

end





% Paper: Transcriptional fluctuations govern the serum dependent cell cycle duration heterogeneities in Mammalian cells
% Author: Vinodhini Govindaraj, Subrot Sarma, Atharva Karulkar, Rahul Purwar and Sandip Kar
% e-mail about the code: vinodhinigovindaraj@gmail.com,sandipkar@iitb.ac.in
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zp=cellcycle_rhs(t,z,param)
%
% all the variables in the cell cycle problem
%	initial values of mrna and proteins
	

	GF=10.0;
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
	%Skp2i=inSkp2T-inSkp2a(idiv);
	CycAm=z(24);
	CycA=z(25);
	Cdc20m=z(26);
	Cdc20T=z(27);
	Cdc20A=z(28);
	%Cdc20i=inCdc20T-inCdc20a(idiv);
	Cdh1m=z(29);
	Cdh1T=z(30);
	Cdh1=z(31);
	%Cdh1i=inCdh1T-inCdh1a(idiv);
	IEP=z(32);
	%IEPi=(10000*is)-inIEP(idiv);
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


% all the parameters in the cell cycle problem
	%parameter values

	%MV
	kmv=param(1);
    smv=param(2);
    gmv=param(3);
	kdmv=param(4);
	
	%mycm	
	bmyc=param(5);
	kMcm1=param(6);
	kMcm2=param(7);	
	kmm=param(8);
	n=param(9);
	kdMcm=param(10);
	%Mycp
	kMcp=param(11);
	kdMcp=param(12);
	kdMcp1=param(13);
	%E2FM
	be2f=param(14);
	kEm1=param(15);
	kEm2=param(16);
	kdEm=param(17);
	%E2Fp
	kEp1=param(18);
	kdEp=param(19);
	kdep1=param(20);
	Jde=param(21);
	%DE
	kDE1=param(22);
	kDE2=param(23);
	Dpt=param(24);
	kdep2=param(25);
	%Rbm
	krbm=param(26);
	kdrbm=param(27);
	%Rb
	krbp=param(28);
	kdrbp=param(29);
	%nRbt=1;
	%iDE
	kRD1=param(30);
	kRD2=param(31);
	kRpD1=param(32);
	kRpD2=param(33);
	%Rip
	kiR1=param(34);
	JRD=param(35);
	kiR2=param(36);
	JR1=param(37);
	kiR3=param(38);
	JRE=param(39);
	kiR4=param(40);
	JR2=param(41);
	%CycDm
	kCDm1=param(42);
	kCDm2=param(43);
	kdCDm=param(44);
	%cyc D
	kCDp1=param(45);
	kdCDp=param(46);
	kdcp1=param(47);
	%CycEm
	bcycem=param(48);
	ksem=param(49);
	kdcem=param(50);
	%cyc E
	kse=param(51);
	kde=param(52);
	kd=param(53);
	kdea=param(54);
	%p27m
	k11m=param(55);
	k11dm=param(56);
	%p27
	ksi=param(57);
	kdi=param(58);
	%EI
	kaei=param(59);
	kdei=param(60);
	%p27p
	kip=param(61);
	kip1=param(62);
	Jce=param(63);
	kdip1=param(64);
	kdip=param(65);	
	%Skp2m
	k17m=param(66);
	k17dm=param(67);
	%Skp2
	k17=param(68);
    k18a=param(69);
	k18b=param(70);
	J18=param(71);
	k19a=param(72);
	k19b=param(73);
	k18c=param(74);
	k18d=param(75);
	J18b=param(76);
	%CycAm
	bsam=param(77);
	ksam=param(78);
	kdam=param(79);
	%Cycad
	ksa=param(80);
	kda=param(81);
	kda1=param(82);
	%Cdh1m
	k3m=param(83);
	k3dm=param(84);
	%Cdh1
	k3=param(85);
	k3d=param(86);
	k3a=param(87);
	k3b=param(88);
	J3=param(89);
	k4=param(90);
	k4b=param(91);
	J4=param(92);
	%Cdc20m
	k5am=param(93);
	k5bm=param(94);
	J5=param(95);
	m=param(96);
	k5dm=param(97);
	%Cdc20
	k5a=param(98);
	k6=param(99);
	k7=param(100);
	J7=param(101);
	k8=param(102);
	J8=param(103);
	%IEP
	k9=param(104);
	k10=param(105);
	Mad=param(106);
	IEPT=param(107);
	%Cdt
	k21m=param(108);
	k21dm=param(109);
	k21=param(110);
	k22a=param(111);
	k22b=param(112);
	J22=param(113);
	%Gem
	k20m=param(114);
	k20dm=param(115);
	k19=param(116);
	k20a=param(117);
	k20b=param(118);
	J20=param(119);
        
%	#CycBm
	kscbm=param(120);
	kdcbm=param(121);
%	#CycB
	kscb=param(122);
	kdcb=param(123);
	kdcb1=param(124);
	kcba=param(125);
	kcbi=param(126);

%	#cdc25m
    ksc25m=param(127);
	kdc25m=param(128);

%	#cdc25
	ksc25=param(129);
	kdc25=param(130);
	kc25a=param(131);
	kc25i=param(132);
	kmc1=param(133);
	kmc2=param(134);

%	#wee1m
	kswm=param(135);
	sgw=param(136);
	kgw=param(137);
	kdwm=param(138);

%	Wee1
	ksw=param(139);
	kdw=param(140);
	kwa=param(141);
	kwi=param(142);
	kmw1=param(143);
	kmw2=param(144);
    
    Jms=param(145);
    Jcd=param(146);
    

 
% all the differential equations in the population space
%Volume
%dV/dt=kve*V*GF/(Vmv+GF)
%dV/dt=kv*GF/(Vmv+GF)

%s=10000
Dp1p=((Dpt)-DE-iDEP);
E2Fp=(Y-DE-iDE-iDEP);
Rb=RbT-iRb1-iRb2-iDE-iDEP;
CycBi=CycBT-CycBa;
C25i=C25T-C25a;
Wee1i=Wee1T-Wee1a;




zp(1)=kmv*GF/(smv+(gmv*GF))-kdmv*MV;
%Myc
zp(2)=((bmyc)+kMcm1*(MV)+kMcm2*(DE^n)/(((kmm)^n)...
    +(DE^n))-kdMcm*Mycm);
zp(3)=(kMcp*(Mycm)-kdMcp*Mycp-(kdMcp1)*Mycp*Skp2/(Jms+Mycp));
%E2F
zp(4)=((be2f)+kEm1*(Mycp)+kEm2*(DE^n)/(((kmm)^n)+...
    (DE^n))-kdEm*E2Fm);
zp(5)=(kEp1*(E2Fm)-kdEp*E2Fp-kdEp*iDE-kdEp*DE-...
    ((kdep1)*CycA*DE/((Jde)+DE))-((kdep1)*CycA*E2Fp/((Jde)+E2Fp)));
zp(6)=((kDE1)*Dp1p*E2Fp-kDE2*DE-(kRD1)*DE*Rb+kRD2*iDE-...
    (kdEp*DE)+(kiR3)*iDEP*CycE/((JRE)+iDEP)...
    -(kRpD1)*DE*iRb1+kRpD2*iDEP-((kdep1)*CycA*DE/((Jde)+DE)));
%RB
zp(7)=krbm-kdrbm*Rbm;
zp(8)=krbp*(Rbm)-kdrbp*(Rb+iRb1+iRb2)-kdEp*(iDE);
%dRb/dt=krbp*Rbm-kdrbp*Rb-kdep*iDE-(kRD1)*DE*Rb+kRD2*iDE-(kiR1)*Rb*CycD/((JRD)+Rb)+kiR2*iRb1/((s*JR1)+iRb1)
%
zp(9)=(kRD1)*DE*Rb-kRD2*iDE-kdEp*iDE-(kiR1)*iDE*CycD/((JRD)+iDE)...
    +kiR2*iDEP/((JR1)+iDEP);
zp(10)=(kiR1)*iDE*CycD/((JRD)+iDE)-kiR2*iDEP/((JR1)+iDEP)-...
    (kiR3)*iDEP*CycE/((JRE)+iDEP)+(kRpD1)*DE*iRb1-kRpD2*iDEP;
zp(11)=(kiR1)*Rb*CycD/((JRD)+Rb)-kiR2*iRb1/((JR1)+iRb1)-...
    (kiR3)*iRb1*CycE/((JRE)+iRb1)+kiR4*iRb2/((JR2)+iRb2)-...
    (kRpD1)*DE*iRb1+kRpD2*iDEP-kdrbp*iRb1;
zp(12)=(kiR3)*iRb1*CycE/((JRE)+iRb1)-kiR4*iRb2/((JR2)+iRb2)+...
    (kiR3)*iDEP*CycE/((JRE)+iDEP)-kdrbp*iRb2;


%
%CyclinD
zp(13)=(kCDm1*(MV))+kCDm2*(Mycp)-kdCDm*CycDm;
zp(14)=kCDp1*(CycDm)-kdCDp*CycD-(kdcp1)*CycD*CycA/(Jcd+CycD);
%Cyclin E
zp(15)=bcycem+ksem*(DE)-kdcem*CycEm;
zp(16)=kse*(CycEm)-kde*CycE-(kd)*CycE*Cdh1-(kaei)*CycE*I+ ...
kdei*EI-(kdea)*(Skp2)*CycE+(kdip)*Skp2*EI;
%CKI
zp(17)=(k11m-k11dm*Im);
zp(18)=ksi*(Im)-kip*(CycE)*I/((Jce)+I)-(kaei)*CycE*I+...
    (kdei)*EI-kdi*(I);
zp(19)=kip*(CycE)*I/((Jce)+I)-kdip1*Ip-(kdip)*Skp2*Ip;
zp(20)=(kaei)*CycE*I-kdei*EI-kde*EI-kdi*EI-(kdip)*Skp2*EI;
%Skp2
zp(21)=(k17m-k17dm*Skp2m);
zp(22)=(k17*(Skp2m)-k18a*(Skp2T-Skp2)-...
    k18b*(Skp2T-Skp2)*Cdh1/((J18)+(Skp2T-Skp2))-...
    k18c*Skp2-k18d*Skp2*Cdh1/((J18b)+Skp2));
zp(23)=((k19a)*CycE*(Skp2T-Skp2)-(k19b*Skp2)-k18c*Skp2-...
    k18d*Skp2*Cdh1/((J18b)+Skp2));

%Cyclin A
zp(24)=(bsam)+ksam*(DE)-kdam*CycAm;
zp(25)=(ksa*(CycAm)-kda*CycA-(kda1)*CycA*Cdh1);
%Cdc20
zp(26)=((k5am+k5bm*((CycBa/(J5))^(m))/(1+(CycBa/(J5))^(m)))-...
    k5dm*Cdc20m);
zp(27)=(k5a*(Cdc20m)-k6*Cdc20T);
zp(28)=((k7)*IEP*(Cdc20T-Cdc20A)/((J7)+Cdc20T-Cdc20A)-...
    k8*Mad*Cdc20A/((J8)+Cdc20A)-k6*Cdc20A);
%Cdh1;
zp(29)=k3m-k3dm*Cdh1m;
zp(30)=(k3*(Cdh1m)-k3d*Cdh1T);
zp(31)=(((k3a)+k3b*Cdc20A)*(Cdh1T-Cdh1)/((J3)+Cdh1T-Cdh1)-...
    ((k4*CycBa))*Cdh1/((J4)+Cdh1)-((k4b*CycA))*Cdh1/((J4)+Cdh1)-...
    k3d*Cdh1);
%IEP
zp(32)=((k9)*CycBa*(IEPT-IEP)-k10*IEP);

%Cdt1 and Geminin
zp(33)=(k21m-k21dm*Cdt1m);
zp(34)=(k20m-k20dm*Gemm);
zp(35)=(k19*(Gemm)-k20a*Gem-k20b*Gem*Cdh1/((J20)+Gem));
zp(36)=(k21*(Cdt1m)-k22a*Cdt1-k22b*Cdt1*Skp2/((J22)+Cdt1));

%CycB
zp(37)=kscbm*(CycA)-kdcbm*CycBm;
zp(38)=(kscb*(CycBm)-kdcb*CycBT-(kdcb1)*CycBT*Cdh1);
zp(39)=(kscb*(CycBm)+(kcba)*CycBi*C25a-(kcbi)*CycBa*Wee1a-...
    kdcb*CycBa-(kdcb1)*CycBa*Cdh1);


%Cdc25
zp(40)=(ksc25m*(CycA))-kdc25m*C25m;
zp(41)=(ksc25*(C25m)-kdc25*C25T);
zp(42)=(kc25a*CycBa*(C25i/(C25i+(kmc1)))-...
    kc25i*(C25a/(C25a+(kmc2)))-kdc25*C25a);
%Wee1
zp(43)=(kswm*GF/(sgw+(kgw*GF)))-kdwm*Wee1m;
zp(44)=(ksw*(Wee1m)-kdw*Wee1T);
zp(45)=(ksw*(Wee1m)+kwa*(Wee1i/(Wee1i+(kmw1)))-...
    kwi*CycBa*(Wee1a/(Wee1a+(kmw2)))-kdw*Wee1a);



zp=zp(:);



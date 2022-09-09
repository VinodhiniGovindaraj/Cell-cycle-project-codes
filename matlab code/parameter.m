function param = parameter(akm,cf)
       %	parameter values

	s1=10000;
	s2=10000;
    trm=1.0;
%	constant
	dh=60.0;
	ds=17.0;
	d=dh/ds;
   

%	#MV
	[dkmv,param(1)]=deal(3.0*s2/cf);
    [smv,param(2)]=deal(2.2);
    [gmv,param(3)]=deal(0.2);
	[dkdmv,param(4)]=deal(4.0);
	

%	#mycm	
	[bmyc,param(5)]=deal(akm(1)*trm*s1);
	[dkMcm1,param(6)]=deal(akm(2)*trm*(s1/s2));
	[dkMcm2,param(7)]=deal(akm(3)*trm*s1);
	[dkmm,param(8)]=deal(0.33*s2);
	[an,param(9)]=deal(2);
	[dkdMcm,param(10)]=deal(1.38);

%	#Mycp
	[dkMcp,param(11)]=deal(40*(s2/s1));
	[dkdMcp,param(12)]=deal(0.7);
	[dkdMcp1,param(13)]=deal(2.0);
  

%	#E2FM
	[be2f,param(14)]=deal(akm(4)*trm*s1);
	[dkEm1,param(15)]=deal(akm(5)*trm*(s1/s2));
	[dkEm2,param(16)]=deal(akm(6)*trm*s1);
	[dkdEm,param(17)]=deal(0.25);

    
%	#E2Fp
	[dkEp1,param(18)]=deal(50.0*(s2/s1));
	[dkdEp,param(19)]=deal(0.25);
	[dkdEp1,param(20)]=deal(1);
	[djde,param(21)]=deal(1*s2);
    
%	#DE
	[dkDE1,param(22)]=deal(871.4/s2);
	[dkDE2,param(23)]=deal(55.0);
	[nDp1t,param(24)]=deal(1*s2);
	[dkdep2,param(25)]=deal(0/s2);
        
%	#Rbm
	[dkrbm,param(26)]=deal(akm(7)*trm*s1);
	[dkdrbm,param(27)]=deal(1.04);

%	#Rb
	[dkrbp,param(28)]=deal(5*(s2/s1));
	[dkdrbp,param(29)]=deal(0.231);
    
%	nRbt=1*s2
%	#iDE
	[dkRD1,param(30)]=deal(100.0/s2);
	[dkRD2,param(31)]=deal(0.5);
	[dkRpD1,param(32)]=deal(0.1/s2);
	[dkRpD2,param(33)]=deal(0.05);
    
%	#Rip
	[dkiR1,param(34)]=deal(0.8);
	[aJRD,param(35)]=deal(0.01*s2);
	[dkiR2,param(36)]=deal(1.0*s2);
	[aJR1,param(37)]=deal(0.01*s2);
	[dkiR3,param(38)]=deal(60.0);
	[aJRE,param(39)]=deal(0.001*s2);
	[dkiR4,param(40)]=deal(25.0*s2);
	[aJR2,param(41)]=deal(0.01*s2);
    
%	#CycDm
	[dkCDm1,param(42)]=deal(akm(8)*trm*(s1/s2));
	[dkCDm2,param(43)]=deal(akm(9)*trm*(s1/s2));
	[dkdCDm,param(44)]=deal(0.173);
  
%	#cyc D
	[dkCDp1,param(45)]=deal(56.67*(s2/s1));
	[dkdCDp,param(46)]=deal(1.386);
	[dkdCp1,param(47)]=deal(2); 
    
%	#CycEm
	[bcycem,param(48)]=deal(akm(10)*trm*s1);
	[dksem,param(49)]=deal(akm(11)*trm*(s1/s2));
	[dkdcem,param(50)]=deal(1.0);
    
%	#cyc E
	[dkse,param(51)]=deal(100.0*(s2/s1));
	[dkde,param(52)]=deal(1);
	[dkd,param(53)]=deal(0.0/s2);
	[dkdea,param(54)]=deal(2/s2);
    
%	#p27m
	[dksim,param(55)]=deal(akm(12)*trm*s1);
	[dkdim,param(56)]=deal(0.693);
    
%	#p27
	[dksi,param(57)]=deal(100.0*(s2/s1));
	[dkdi,param(58)]=deal(1.0);
    
%	#EI
	[dkaei,param(59)]=deal(10.0/s2);
	[dkdei,param(60)]=deal(1.0);
    
%	#p27p
	[dkip,param(61)]=deal(5.0);
	[dkip1,param(62)]=deal(0.0);
	[aJce,param(63)]=deal(0.1*s2);
	[dkdip1,param(64)]=deal(1.0);
	[dkdip,param(65)]=deal(20.0/s2);
%	#Skp2m
	[ak17m,param(66)]=deal(akm(13)*trm*s1);
	[ak17dm,param(67)]=deal(0.17);
   
%	#Skp2
	[ak17,param(68)]=deal(7.5*d*(s2/s1));
	[ak18a,param(69)]=deal(5.5*d);
	[ak18b,param(70)]=deal(15.0*d);
	[aJ18,param(71)]=deal(0.001*s2);
	[ak19a,param(72)]=deal(20.0*d/s2);
	[ak19b,param(73)]=deal(0.1*d);
	[ak18c,param(74)]=deal(0.008*d);
	[ak18d,param(75)]=deal(0.1*d);
	[aJ18b,param(76)]=deal(0.01*s2);
%	#CycAm
	[bsam,param(77)]=deal(akm(14)*trm*s1);
	[dksam,param(78)]=deal(akm(15)*trm*(s1/s2));
	[dkdam,param(79)]=deal(1.0);
    
%	#Cycad
	[dksa,param(80)]=deal(2.0*d*(s2/s1));
	[dkda,param(81)]=deal(0.04*d);
	[dkda1,param(82)]=deal(1.0*d/s2);
    
%	#Cdh1m
	[ak3m,param(83)]=deal(akm(16)*trm*s1);
	[ak3dm,param(84)]=deal(3.0);
    
%	#Cdh1
	[ak3,param(85)]=deal(50*d*(s2/s1));
	[ak3d,param(86)]=deal(1*d);
	[ak3a,param(87)]=deal(1.0*d*s2);
	[ak3b,param(88)]=deal(50.0*d);
	[aJ3,param(89)]=deal(0.04*s2);
	[ak4,param(90)]=deal(35*d);
	[ak4b,param(91)]=deal(35*d);
	[aJ4,param(92)]=deal(0.04*s2);
%	#Cdc20m
	[ak5am,param(93)]=deal(akm(17)*trm*s1*d);
	[ak5bm,param(94)]=deal(akm(18)*trm*d*s1);
	[aJ5,param(95)]=deal(0.3*s2);
	[am,param(96)]=deal(4);
	[ak5dm,param(97)]=deal(0.1*d);
%	#Cdc20
	[ak5a,param(98)]=deal(5*d*(s2/s1));
	[ak6,param(99)]=deal(0.1*d);
	[ak7,param(100)]=deal(1.0*d);
	[aJ7,param(101)]=deal(0.05*s2);
	[ak8,param(102)]=deal(0.5*d);
	[aJ8,param(103)]=deal(0.05*s2);
%	#IEP
	[ak9,param(104)]=deal(0.1*d/s2);
	[ak10,param(105)]=deal(0.02*d);
	[nMad,param(106)]=deal(1*s2);
	[nIEPT,param(107)]=deal(1*s2);
%	#Cdt
	[ak21m,param(108)]=deal(akm(19)*trm*s1);
	[ak21dm,param(109)]=deal(0.35);
	[ak21,param(110)]=deal(2.5*d*(s2/s1));
	[ak22a,param(111)]=deal(0.06*d);
	[ak22b,param(112)]=deal(2*d);
	[aJ22,param(113)]=deal(0.04*s2);
%	#Gem
	[ak20m,param(114)]=deal(akm(20)*trm*s1);
	[ak20dm,param(115)]=deal(0.35);
	[ak19,param(116)]=deal(2.5*d*(s2/s1));
	[ak20a,param(117)]=deal(0.06*d);
	[ak20b,param(118)]=deal(1*d);
	[aJ20,param(119)]=deal(0.5*s2);
    
%	#CycBm
	[dkscbm,param(120)]=deal(akm(21)*trm*(s1/s2));
	[dkdcbm,param(121)]=deal(0.1);
%	#CycB
	[dkscb,param(122)]=deal(1*d*(s2/s1));
	[dkdcb,param(123)]=deal(0.04*d);
	[dkdcb1,param(124)]=deal(1.0*d/s2);
	[dkcba,param(125)]=deal(10.0*d/s2);
	[dkcbi,param(126)]=deal(10.0*d/s2);

%	#cdc25m
	[dksc25m,param(127)]=deal(akm(22)*trm*(s1/s2));
	[dkdc25m,param(128)]=deal(1.0);

%	#cdc25
	[dksc25,param(129)]=deal(50.0*d*(s2/s1));
	[dkdc25,param(130)]=deal(1.0*d);
	[dkc25a,param(131)]=deal(10.0*d);
	[dkc25i,param(132)]=deal(0.4*s2*d);
	[dkmc1,param(133)]=deal(0.05*s2);
	[dkmc2,param(134)]=deal(0.05*s2);

%	#wee1m
	[dkswm,param(135)]=deal(akm(23)*trm*s1);
	[dsgw,param(136)]=deal(1);
	[dkgw,param(137)]=deal(0.5);
	[dkdwm,param(138)]=deal(1.0);

%	Wee1
	[dksw,param(139)]=deal(50.0*d*(s2/s1));
	[dkdw,param(140)]=deal(1.0*d);
	[dkwa,param(141)]=deal(4.0*s2*d);
	[dkwi,param(142)]=deal(50.0*d);
	[dkmw1,param(143)]=deal(0.2*s2);
	[dkmw2,param(144)]=deal(0.2*s2);
    
    [Jms,param(145)]=deal(0.001*s2);
    [Jcd,param(146)]=deal(0.01*s2); 
end


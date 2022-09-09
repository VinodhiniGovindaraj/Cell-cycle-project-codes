c	Paper: Transcriptional fluctuations govern the serum dependent cell cycle duration heterogeneities in Mammalian cells
c	Author: Vinodhini Govindaraj, Subrot Sarma, Atharva Karulkar, Rahul Purwar and Sandip Kar
c	e-mail about the code: vinodhinigovindaraj@gmail.com,sandipkar@iitb.ac.in

c	10% CV for asymmetric cell division
c	Transcription rate variation


	implicit real*8(a-h,q-z)

	parameter(nu=129,nr=100000000,ncc=2000,nc=2000)
	dimension a(nu),c(nu)
	dimension nG(nr),tG(nr)
	dimension teg(ncc),eg(ncc),tmm(ncc),mm(ncc),tem(ncc),em(ncc)
	dimension tcc(ncc),g1(ncc),sg2m(ncc)


c	array to store  no. of molecules after division and division time for each cell
c	This used will be used as initial conditions for daughter cells

	dimension inMV(nc),inMycm(nc),inMycp(nc),inE2Fm(nc),inE2Fp(nc)
	dimension inDE(nc),iniDE(nc),iniDEP(nc),inRbm(nc),inRb(nc)
	dimension iniRb1(nc),iniRb2(nc),inCycDm(nc),inCycD(nc)
	dimension inCycEm(nc),inCycE(nc),inI(nc),inIp(nc),inXI(nc)
	dimension inSkp2a(nc),inSkp2i(nc),inCdh1a(nc),inCdh1i(nc)
	dimension inCdc20a(nc),inCdc20i(nc),inIEPi(nc),inIEP(nc)
	dimension inIm(nc),inSkp2m(nc),inCdt1m(nc),inGemm(nc),inCdt1(nc)
	dimension inGem(nc),inCycAm(nc),inCycA(nc),inCdh1m(nc),inC20m(nc)
	dimension inCycBm(nc),inCycBa(nc),inCycBi(nc)
	dimension inWee1m(nc),inWee1i(nc),inWee1a(nc)
	dimension inC25m(nc),inC25i(nc),inC25a(nc)

	dimension GFi(nc),tid(nc)

c	array for varying rate constants
	parameter(ku=23)
	dimension akm(ku),akn(ku),akc(ku,ncc)
	dimension ark(ku),rsigma(ku),rmu(ku),tcv(ku)


	
	interface
		function rates(akm,ku) result(akn)

			implicit real*8 (a-h,q-z)
			dimension akm(ku),akn(ku)
			end
	end interface

c	no. of generations to be followed
	nk=100
	nrows=1
c	cell or generation for which output file is created
	kl=1
	ki=0
	nj=1

c	For calculating log mean value of TRs
	sig=1.0
	
c	transcription rate CV for the initial population
	cvp=0.2

c	for calculating size at division
	sigma=0.05

c	some initial parameters
	l100=500
	imaximum=1000
	nrxn=0




c	open(400,file='g1phase.out',status='unknown')
c	open(401,file='sg2mphase.out',status='unknown')
c	open(402,file='total.out',status='unknown')
c	open(403,file='tcg.out',status='unknown')
c	
	open(301,file='Volume.out',status='unknown')
	open(302,file='GF.out',status='unknown')

	open(201,file='CycA.out',status='unknown')
c	open(202,file='Cdh1i.out',status='unknown')
	open(203,file='Cdh1a.out',status='unknown')
c	open(204,file='Cdc20i.out',status='unknown')
	open(205,file='Cdc20a.out',status='unknown')
c	open(206,file='IEPi.out',status='unknown')
c	open(207,file='IEPa.out',status='unknown')
c	open(208,file='Rb.out',status='unknown')
	open(209,file=' Cdh1m.out',status='unknown')
c	open(210,file=' Cdc20m.out',status='unknown')
c
c	open(211,file=' Mycm.out',status='unknown')
c	open(212,file='Mycp.out',status='unknown')
c	open(213,file=' E2Fm.out',status='unknown')
c	open(214,file='E2FT.out',status='unknown')
c	open(215,file='DE.out',status='unknown')
c	open(216,file='IDE.out',status='unknown')	
c	open(217,file=' CycDm.out',status='unknown')
c	open(218,file='CycDp.out',status='unknown')
c	open(219,file=' CycEm.out',status='unknown')
c	open(220,file='CycEp.out',status='unknown')
c	open(221,file='Rbp.out',status='unknown')
c	open(222,file='Rbpp.out',status='unknown')
c
c	open(223,file='Skp2a.out',status='unknown')
c	open(224,file='Skp2i.out',status='unknown')
	open(225,file='Cdt1.out',status='unknown')
	open(226,file='Geminin.out',status='unknown')
c	open(227,file=' Cdt1m.out',status='unknown')
c	open(228,file=' Gemm.out',status='unknown')
c	open(229,file=' Skp2m.out',status='unknown')
c	open(230,file=' CKIm.out',status='unknown')
c	open(231,file='CKI.out',status='unknown')
c	open(232,file='CKIp.out',status='unknown')
c	open(233,file='EI.out',status='unknown')
c	open(234,file=' CycAm.out',status='unknown')
c	open(235,file=' Rbm.out',status='unknown')
c	open(236,file='RbT.out',status='unknown')
c	open(237,file='MV.out',status='unknown')
c	open(238,file='IDEP.out',status='unknown')

	open(239,file=' CycBm.out',status='unknown')
	open(240,file='CycBa.out',status='unknown')
c	open(241,file=' C25m.out',status='unknown')
	open(242,file='C25a.out',status='unknown')
c	open(243,file=' Wee1m.out',status='unknown')
	open(244,file='Wee1a.out',status='unknown')
c	open(245,file=' CycDm.out',status='unknown')
c	open(246,file='CycDp.out',status='unknown')


	
c	open(245,file='bmyc_kv_full.out',status='unknown')
c	open(247,file='aatimings_raw.txt',status='unknown')
	open(248,file='aatimings.txt',status='unknown')
c	open(249,file='d1.txt',status='unknown')
c	open(250,file='d2.txt',status='unknown')

c	for writing transcripton rates
	open(601,file='trates_birth.out',status='unknown')
	open(602,file='trates_division.out',status='unknown')
	open(603,file='trates_mitosis.out',status='unknown')
	open(604,file='trates_sphase.out',status='unknown')
	open(605,file='trates_mphase.out',status='unknown')

	
c	write(601,*)'bmyc','kMcm1','kMcm2','be2f','kEm1','kEm2'
c	write(602,*)'kRbm','kCDm1','kCDm2','bCycEm','kCycEm','kCKIm'
c	write(603,*)'kSkp2m','bCycAm','kCycAm','kCdh1m','bC20m','kC20m'
c	write(604,*)'kCdt1m','kGemm','kMV'
	
	do 300 jk=1,nk
c
c	call random seed
		l100=l100+100
		iseed=l100
		useed=dustar(iseed)	


c	cycle counts
	ij=0
	itr=0
	itrt=0
	itb=1
	itrn=1
	idiv=1
	idiva=0
	idivb=1


c	initial time,final time and initialization of number of reactions
	ti=0.0 
	tf=72.0
c	Growth factor and volume related values 

	dGF=2.0
	dkdGF=0.0



c	initial values of mrna and proteins
	


	GFi(idiv)=dGF

	inMV(idiv)=5769
	inMycm(idiv)=138
	inMycp(idiv)=2
	inE2Fm(idiv)=62
	inY=9392
	inDE(idiv)=6101
	inRbm(idiv)=480
	inRbt=10206
	iniDE(idiv)=1599
	iniDEP(idiv)=1
	iniRb1(idiv)=8
	iniRb2(idiv)=8538
c
	inCycDm(idiv)=171
	inCycD(idiv)=829
	inCycEm(idiv)=308
	inCycE(idiv)=8882
	inIm(idiv)=204
	inI(idiv)=524
	inIp(idiv)=548
	inXI(idiv)=155
	inSkp2m(idiv)=205
	inSkp2T=13615
	inSkp2a(idiv)=13507
	inSkp2i(idiv)=inSkp2T-inSkp2a(idiv)
c
	inCycAm(idiv)=154
	inCycA(idiv)=5567


    	inCycBm(idiv)=445
	inCycBa(idiv)=174
    	inCycBT=4259
    	
    	inCycBi(idiv)=inCycBT-inCycBa(idiv)

    	inC25m(idiv)=102
    	inC25T=5007
    	inC25a(idiv)=231
    	inC25i(idiv)=inC25T-inC25a(idiv)

    	inWee1m(idiv)=159
    	inWee1T=7999
    	inWee1a(idiv)=7618
    	inWee1i(idiv)=inWee1T-inWee1a(idiv)

	inCdh1m(idiv)=199
	inCdh1T=10000
	inCdh1a(idiv)=27
	inCdh1i(idiv)=inCdh1T-inCdh1a(idiv)

	inC20m(idiv)=10
	inCdc20T=622
	inCdc20a(idiv)=65
	inCdc20i(idiv)=inCdc20T-inCdc20a(idiv)

	inIEP(idiv)=1105
	inIEPi(idiv)=(10000)-inIEP(idiv)

	inCdt1m(idiv)=214
	inGemm(idiv)=214
	inGem(idiv)=6699
	inCdt1(idiv)=8

	inE2Fp(idiv)=inY-inDE(idiv)-iniDE(idiv)-iniDEP(idiv)	
	ipRb=iniRb1(idiv)+iniRb2(idiv)
	inRb(idiv)=inRbT-ipRb-iniDE(idiv)-iniDEP(idiv)

	tid(idiv)=0


c	initial rates for varying rate constants(without scaling factor)
c	#mycm	
	bmycn=0.000033 
	dkMcm1n=0.033 
	dkMcm2n=0.000015 
c	#E2FM		
	be2fn=0.0000033 
	dkEm1n=0.005 
	dkEm2n=0.0005 
c	#Rbm
	dkrbmn=0.05 
c	#CycDm
	dkCDm1n=0.000001
	dkCDm2n=0.01
c	#CycEm
	bcycemn=0.000002
	dksemn=0.04
c	#p27m
	dksimn=0.0142
c	#Skp2m
	ak17mn=0.0036
c	#CycAm
	bsamn=0.000001
	dksamn=0.02
c	#Cdh1m
	ak3mn=0.06
c	#Cdc20m
	ak5amn=0.0001
	ak5bmn=0.004
c	#Cdt
	ak21mn=0.0075
c	#Gem
	ak20mn=0.0075
c	#CycBm
	dkscbmn=0.02
c	#C25m
	dksc25mn=0.02
c	#Wee1m
	dkswmn=0.016


c	putting transcription rates in an array
	
c	#mycm
	ark(1)=bmycn
	ark(2)=dkMcm1n
	ark(3)=dkMcm2n	
c	#E2FM
	ark(4)=be2fn
	ark(5)=dkEm1n
	ark(6)=dkEm2n
c	#Rbm
	ark(7)=dkrbmn
c	#CycDm
	ark(8)=dkCDm1n
	ark(9)=dkCDm2n
c	#CycEm
	ark(10)=bcycemn
	ark(11)=dksemn
c	#p27m
	ark(12)=dksimn
c	#Skp2m
	ark(13)=ak17mn
c	#CycAm
	ark(14)=bsamn
	ark(15)=dksamn
c	#Cdh1m
	ark(16)=ak3mn
c	#Cdc20m
	ark(17)=ak5amn
	ark(18)=ak5bmn
c	#Cdt
	ark(19)=ak21mn
c	#Gem
	ark(20)=ak20mn
c	#CycBm
	ark(21)=dkscbmn
c	#C25m
	ark(22)=dksc25mn
c	#Wee1m
	ark(23)=dkswmn

c	calculating log values for transcription rates


c	#CV for the transcription rates
	do 39 it=1,ku
		tcv(it)=cvp
 39	continue

c	putting transcription rates in an array

c	picking the new rates for each daughter pair

	do  43 ir=1,ku

		rsigma(ir)=sqrt(log(1.0+(tcv(ir)*tcv(ir))))
		rmu(ir)=log(ark(ir))-(rsigma(ir)*rsigma(ir))/2.0

 43	continue


 	 do 63 jr=1,ku
     	       rk1 = duni()
               rk2 = duni()
               api = 4.0*atan(1.0)

               if(rk1.le.0.0)then
                  rk1 = 0.00001
               endif
               
                ba  = -2.0*sig*sig*log(rk1)
                ba1 = dcos(2.0*api*rk2)
                sf  = sqrt(ba)*ba1
                akc(jr,1)=exp(rmu(jr)+rsigma(jr)*sf)


 63	continue

c	write(*,*)akc(1,1),rmu(1),sf,rsigma(1),tcv(1)
	
c	following the generations

	do 90 ii=0,nrows

c	confluence factor

		if(ii.lt.1)then
			cf=1
		else
			cf=cf+0.4
		end if


	do 80 jj=1,(2**ii)

c	initialization for period calculation
	maxc=0
	maxg=0

	tc=0
	ncell=1
	eg(ncell)=maxc
	teg(ncell)=tc
	mcount=1
	ncount=0
	tcount=0


c	To update initial no. of molecules and initial time for the daughter cells
c	ij is updated for every odd jj count,to get the initial conditions from corresponding mother cell


c 		if (mod(jj,2).ne.0) then
				ij=ij+1
c			endif

c	initial values for the proteins	
			nMV=inMV(ij)
			nMycm=inMycm(ij)
			nMycp=inMycp(ij)
			nE2Fm=inE2Fm(ij)
			nE2Fp=inE2Fp(ij)
			nDE=inDE(ij)
			niDE=iniDE(ij)
			niDEP=iniDEP(ij)
			nRbm=inRbm(ij)
			nRb=inRb(ij)
			niRb1=iniRb1(ij)
			niRb2=iniRb2(ij)
			nCycDm=inCycDm(ij)
			nCycD=inCycD(ij)
			nCycEm=inCycEm(ij)
			nCycE=inCycE(ij)
			nI=inI(ij)
			nIp=inIp(ij)
			nXI=inXI(ij)
			nSkp2a=inSkp2a(ij)
			nSkp2i=inSkp2i(ij)
			nCdc20a=inCdc20a(ij)
			nCdc20i=inCdc20i(ij)
			nCdh1a=inCdh1a(ij)
			nCdh1i=inCdh1i(ij)
c			nIEPi=inIEPi(ij)
			nIEP=inIEP(ij)
			nIm=inIm(ij)
			nSkp2m=inSkp2m(ij)
			nCdt1m=inCdt1m(ij)
			nGemm=inGemm(ij)
			nCdt1=inCdt1(ij)
			nGem=inGem(ij)
			nCycAm=inCycAm(ij)
			nCycA=inCycA(ij)
			nCdh1m=inCdh1m(ij)
			nC20m=inC20m(ij)


			nCycBm=inCycBm(ij)
			nCycBa=inCycBa(ij)
			nCycBi=inCycBi(ij)

			nC25m=inC25m(ij)
			nC25a=inC25a(ij)
			nC25i=inC25i(ij)

			nWee1m=inWee1m(ij)
			nWee1i=inWee1i(ij)
			nWee1a=inWee1a(ij)


c	initial time point for the cell
			ti=tid(ij)
			tmit=tid(ij)
			

c			GF at birth

			GFb=GFi(ij)

c   copying mother transcription rates
	

		do 26 jr=1,ku

			akm(jr)=akc(jr,ij)		
        
 26		continue

		write(601,*)jk,ii,jj
		write(601,*)(akm(jr),jr=1,ku)





c	to stop following a lineage after end of time 
	if(ti.ge.tf) goto 65
	
	mstart=0
	mend=1
	ecount=0



 10	dGF=GFb*exp(-dkdGF*(ti-tb))


	if((ti.lt.tmit).and.(mstart.eq.1))then

		trm=1
		mend=1
		mstart=0

	write(603,*)jk,ii,jj,ti,tmit,trm,bmycn*trm,'3'

		
	else if((ti.ge.tmit).and.(mend.eq.1)) then
		trm=1.0
		mstart=1
		mend=0
		ecount=ecount+1

	  if(ecount.eq.2)then


		akn=rates(akm,ku)

		do  jr=1,ku

			akm(jr)=akn(jr)
			akc(jr,ij)=akn(jr)
		

		end do

		write(605,*)jk,ii,jj
		write(605,*)(akm(jr),jr=1,ku)

		ecount=0
 
		write(603,*)jk,ii,jj,ti,tmit,trm,akc(1,ij)*trm,'4'

	  end if				
		
	end if

c	parameter values

	s=10000
	


c	#mycm
	bmycn=akc(1,ij)
	dkMcm1n=akc(2,ij)
	dkMcm2n=akc(3,ij)	
c	#E2FM
	be2fn=akc(4,ij)
	dkEm1n=akc(5,ij)
	dkEm2n=akc(6,ij)
c	#Rbm
	dkrbmn=akc(7,ij)
c	#CycDm
	dkCDm1n=akc(8,ij)
	dkCDm2n=akc(9,ij)
c	#CycEm
	bcycemn=akc(10,ij)
	dksemn=akc(11,ij)
c	#p27m
	dksimn=akc(12,ij)
c	#Skp2m
	ak17mn=akc(13,ij)
c	#CycAm
	bsamn=akc(14,ij)
	dksamn=akc(15,ij)
c	#Cdh1m
	ak3mn=akc(16,ij)
c	#Cdc20m
	ak5amn=akc(17,ij)
	ak5bmn=akc(18,ij)
c	#Cdt
	ak21mn=akc(19,ij)
c	#Gem
	ak20mn=akc(20,ij)
c	#CycBm
	dkscbmn=akc(21,ij)
c	#C25m
	dksc25mn=akc(22,ij)
c	#Wee1m
	dkswmn=akc(23,ij)



c	constant
	dh=60.0
	ds=17.0
	d=dh/ds
c	write(*,*)nC25a,nC25i,nCycBa,nCycBi


c	#MV
	dkmv=3*s/cf
	dsmv=2.2
	dgmv=0.2
	dkdmv=4.0
c	#mycm	
	bmyc=bmycn*trm*s
	dkMcm1=dkMcm1n*trm
	dkMcm2=dkMcm2n*trm*s	
	dkdMcm=1.38
	dkmm=0.33*s
	an=2
c	#Mycp
	dkMcp=40
	dkdMcp=0.7
	dkdMcp1=2.0
	dJms=0.001*s
c	#E2FM
	be2f=be2fn*trm*s
	dkEm1=dkEm1n*trm
	dkEm2=dkEm2n*trm*s
	dkdEm=0.25
c	#E2Fp
	dkEp1=50.0
	dkdEp=0.25
	dkdEp1=1
	djde=1*s
c	#DE
	dkDE1=871.4/s
	dkDE2=55.0
	nDp1t=1*s
c	#Rbm
	dkrbm=dkrbmn*trm*s
	dkdrbm=1.04
c	#Rb
	dkrbp=5
	dkdrbp=0.231
c	#iDE
	dkRD1=100.0/s
	dkRD2=0.5
	dkRpD1=0.1/s
	dkRpD2=0.05
c	#Rip
	dkiR1=0.8
	aJRD=0.01*s
	dkiR2=1*s
	aJR1=0.01*s
	dkiR3=60.0
	aJRE=0.001*s
	dkiR4=25.0*s
	aJR2=0.01*s
c	#CycDm
	dkCDm1=dkCDm1n*trm
	dkCDm2=dkCDm2n*trm
	dkdCDm=0.173
c	#cyc D
	dkCDp1=56.67
	dkdCDp=1.386
	dkdCp1=2
	djcd=0.01*s
c	#CycEm
	bcycem=bcycemn*trm*s
	dksem=dksemn*trm
	dkdcem=1.0
c	#cyc E
	dkse=100.0
	dkde=1
	dkdea=2/s
c	#p27m
	dksim=dksimn*trm*s
	dkdim=0.693
c	#p27
	dksi=100.0
	dkdi=1.0
c	#EI
	dkaei=10.0/s
	dkdei=1.0
c	#p27p
	dkip=5
	aJce=0.1*s
	dkdip=20.0/s
	dkdip1=1.0	
c	#Skp2m
	ak17m=ak17mn*trm*s
	ak17dm=0.17
c	#Skp2
	ak17=7.5*d
	ak18a=5.5*d
	ak18b=15.0*d
	aJ18=0.001*s
	ak19a=20.0*d/s
	ak19b=0.1*d
	ak18c=0.008*d
	ak18d=0.1*d
	aJ18b=0.01*s
c	#CycAm
	bsam=bsamn*trm*s
	dksam=dksamn*trm
	dkdam=1
c	#Cycad
	dksa=2*d
	dkda=0.04*d
	dkda1=1*d/s
c	#CycBm
	dkscbm=dkscbmn*trm
	dkdcbm=0.1
c	#CycB
	dkscb=1*d
	dkdcb=0.04*d
	dkdcb1=1.0*d/s
	dkcba=10.0*d/s
	dkcbi=10.0*d/s

c	#cdc25m
	dksc25m=dksc25mn*trm
	dkdc25m=1.0

C	#cdc25
	dksc25=50.0*d
	dkdc25=1.0*d
	dkc25a=10.0*d
	dkc25i=0.4*s*d
	dkmc1=0.05*s
	dkmc2=0.05*s

c	#wee1m
	dkswm=dkswmn*trm*s
	dsgw=1
	dkgw=0.5
	dkdwm=1.0

c	Wee1
	dksw=50.0*d
	dkdw=1.0*d
	dkwa=4.0*s*d
	dkwi=50.0*d
	dkmw1=0.2*s
	dkmw2=0.2*s

c	#Cdh1m
	ak3m=ak3mn*trm*s
	ak3dm=3
c	#Cdh1
	ak3=50*d
	ak3d=1*d
	ak3a=1.0*d*s
	ak3b=50.0*d
	aJ3=0.04*s
	ak4=35*d
	ak4b=35*d
	aJ4=0.04*s
c	#Cdc20m
	ak5am=ak5amn*trm*s*d
	ak5bm=ak5bmn*trm*d*s
	aJ5=0.3*s
	am=4
	ak5dm=0.1*d
c	#Cdc20
	ak5a=5*d
	ak6=0.1*d
	ak7=1.0*d
	aJ7=0.05*s
	ak8=0.5*d
	aJ8=0.05*s
c	#IEP
	ak9=0.1*d/s
	ak10=0.02*d
	nMad=1*s
	nIEPT=1*s
c	#Cdt
	ak21m=ak21mn*trm*s
	ak21dm=0.35
	ak21=2.5*d
	ak22a=0.06*d
	ak22b=2*d
	aJ22=0.04*s
c	#Gem
	ak20m=ak20mn*trm*s
	ak20dm=0.35	
	ak19=2.5*d
	ak20a=0.06*d
	ak20b=1*d
	aJ20=0.5*s


c	propensity functions	
c	some total proteins 
	nY=nE2Fp+nDE+niDE+niDEP	
	nDp1p=nDp1t-nDE-niDE-niDEP
	nRbt=nRb+niRb1+niRb2+niDE+niDEP
	nIEPi=nIEPT-nIEP
	nSkp2T=nSkp2a+nSkp2i
	nCdc20T=nCdc20a+nCdc20i
	nCdh1T=nCdh1i+nCdh1a
	
	nCycBT=nCycBi+nCycBa
	nWee1T=nWee1a+nWee1i
	nC25T=nC25a+nC25i

c	#MV
	a(1)=dkmv*dGF/(dsmv+(dgmv*dGF))
	a(2)=dkdmv*nMV
c	#mycm
	a(3)=bmyc
	a(4)=dkMcm1*nMV
	a(5)=(dkMcm2*(nDE**an))/((dkmm**an)+(nDE**an))
	a(6)=dkdMcm*nMycm
c	#Mycp
	a(7)=dkMcp*nMycm
	a(8)=dkdMcp*nMycp
	a(9)=dkdMcp1*nMycp*nSkp2a/(dJms+nMycp)
c	#E2Fm
	a(10)=be2f
	a(11)=dkEm1*nMycp
	a(12)=(dkEm2*(nDE**an))/((dkmm**an)+(nDE**an))
	a(13)=dkdEm*nE2Fm
c	#E2FP
	a(14)=dkEp1*nE2Fm
	a(15)=dkdEp*nE2Fp
	a(16)=dkdEp1*nCycA*nE2Fp/(djde+nE2Fp)
c	#DE
	a(17)=dkDE1*nDp1p*nE2Fp
	a(18)=dkDE2*nDE
	a(19)=dkdep*nDE
	a(20)=dkdep1*nDE*nCycA/(djde+nDE)
c	#rbm
	a(21)=dkrbm
	a(22)=dkdrbm*nRbm
	a(23)=dkrbp*nRbm
	a(24)=dkdrbp*nRb
	a(25)=dkdrbp*niRb1
	a(26)=dkdrbp*niRb2
c	#iDE	
	a(27)=dkRD1*nDE*nRb
	a(28)=dkRD2*niDE
	a(29)=dkdep*niDE
c	#IDEP
	a(30)=dkRpD1*nDE*niRb1
	a(31)=dkRpD2*niDEp
	a(32)=dkiR1*niDE*nCycD/(aJRD+niDE)
	a(33)=dkiR2*niDEP/(aJR1+niDEP)
	a(34)=dkiR3*niDEP*nCycE/(aJRE+niDEP)
c	#Rb
	a(35)=dkiR1*nRb*nCycD/(aJRD+nRb)
	a(36)=dkiR2*niRb1/(aJR1+niRb1)
	a(37)=dkiR3*niRb1*nCycE/(aJRE+niRb1)
	a(38)=dkiR4*niRb2/(aJR2+niRb2)
c	#CycDm
	a(39)=dkCDm1*nMV
	a(40)=dkCDm2*nMycp
	a(41)=dkdCDm*nCycDm
c	#CycD
	a(42)=dkCDp1*nCycDm
	a(43)=dkdCDp*nCycD
	a(44)=dkdCp1*nCycD*nCycA/(djcd+nCycD)
c	#CycEm 
	a(45)=bcycem
	a(46)=dksem*nDE
c	#CycE
	a(47)=dkdcem*nCycEm
	a(48)=dkse*nCycEm
	a(49)=dkde*nCycE	
	a(50)=dkdea*nSkp2a*nCycE
c	#p27m
	a(51)=dksim
	a(52)=dkdim*nIm
c	#p27
	a(53)=dksi*nIm
	a(54)=dkdi*nI
c	#p27p
	a(55)=dkip*nCycE*nI/(aJce+nI)
	a(56)=dkdip1*nIp
	a(57)=dkdip*nIp*nSkp2a
c	#EI
	a(58)=dkaei*nCycE*nI
	a(59)=dkdei*nXI
	a(60)=dkdi*nXI
	a(61)=dkde*nXI
	a(62)=dkdip*nXI*nSkp2a
c	#Skp2m
	a(63)=ak17m
	a(64)=ak17dm*nSkp2m
c	#Skp2p
	a(65)=ak17*nSkp2m
	a(66)=ak18a*nSkp2i
	a(67)=ak18b*nSkp2i*nCdh1a/(aJ18+nSkp2i)
	a(68)=ak19a*nCycE*nSkp2i
	a(69)=ak19b*nSkp2a
	a(70)=ak18c*nSkp2a
	a(71)=ak18d*nSkp2a*nCdh1a/(aJ18b+nSkp2a) 
c	#CycAm
	a(72)=bsam
	a(73)=dksam*nDE
	a(74)=dkdam*nCycAm
c	#CycA
	a(75)=dksa*nCycAm
	a(76)=dkda*nCycA
	a(77)=dkda1*nCycA*nCdh1a
c	#CycBm
	a(78)=dkscbm*nCycA
	a(79)=dkdcbm*nCycBm
c	#CycB
	a(80)=dkscb*nCycBm
	a(81)=dkdcb*nCycBa
	a(82)=dkdcb1*nCycBa*nCdh1a
	a(83)=dkcbi*nCycBa*nWee1a
	a(84)=dkcba*nCycBi*nC25a
	a(85)=dkdcb*nCycBi
	a(86)=dkdcb1*nCycBi*nCdh1a
c	#Cdc25
	a(87)=dksc25m*nCycA
	a(88)=dkdc25m*nC25m
	a(89)=dksc25*nC25m
	a(90)=dkdc25*nC25i
	a(91)=dkdc25*nC25a
	a(92)=dkc25a*nCycBa*(nC25i/(nC25i+dkmc1))
	a(93)=dkc25i*(nC25a/(nC25a+dkmc2))
c	#Wee1
	a(94)=dkswm*dGF/(dsgw+(dkgw*dGF))
	a(95)=dkdwm*nWee1m
	a(96)=dksw*nWee1m
	a(97)=dkdw*nWee1a
	a(98)=dkwa*(nWee1i/(nWee1i+dkmw1))
	a(99)=dkwi*nCycBa*(nWee1a/(nWee1a+dkmw2))
	a(100)=dkdw*nWee1i
c	#Cdh1m
	a(101)=ak3m
	a(102)=ak3dm*nCdh1m
c	#Cdh1
	a(103)=ak3*nCdh1m
	a(104)=ak3d*nCdh1i
	a(105)=ak3a*nCdh1i/(aJ3+nCdh1i)
	a(106)=ak3b*nCdh1i*nCdc20a/(aJ3+nCdh1i)
	a(107)=ak4*nCycBa*nCdh1a/(aJ4+nCdh1a)
	a(108)=ak4b*nCycA*nCdh1a/(aJ4+nCdh1a)
	a(109)=ak3d*nCdh1a
c	#Cdc20m
	a(110)=ak5am
	a(111)=ak5bm*((nCycBa/aJ5)**am)/(1+((nCycBa/aJ5)**am))
	a(112)=ak5dm*nC20m
c	#Cdc20
	a(113)=ak5a*nC20m
	a(114)=ak6*nCdc20i
	a(115)=ak7*nCdc20i*nIEP/(aJ7+nCdc20i)
	a(116)=ak8*nCdc20a*nMad/(aJ8+nCdc20a)
	a(117)=ak6*nCdc20a
c	#IEP
	a(118)=ak9*nIEPi*nCycBa
	a(119)=ak10*nIEP
c	#Cdt1
	a(120)=ak21m
	a(121)=ak21dm*nCdt1m
	a(122)=ak21*nCdt1m
	a(123)=ak22a*nCdt1
	a(124)=ak22b*nCdt1*nSkp2a/(aJ22+nCdt1)
c	#Gem
	a(125)=ak20m
	a(126)=ak20dm*nGemm
	a(127)=ak19*nGemm
	a(128)=ak20a*nGem
	a(129)=ak20b*nGem*nCdh1a/(aJ20+nGem)


c	Updating the number of molecules for each reaction
c
	sum1=0.0
	do 45 i=1,nu
		sum1=sum1+a(i)
	if((jk.eq.kl).and.(ii.eq.ki).and.(jj.eq.nj))then
c	write(*,*)sum1,i,a(i)
	end if
 45	continue
	at=sum1
c
c	Calculation of random number for tau
	r1=duni()
        if(r1.eq.0.0)then
         r1=0.0001
         endif
	tau=log(1/r1)/at
c	write(*,*)at,tau
c
c	Calculation of random number for next reaction
	r2=duni()
	r2at=r2*at
	sum2=0.0
	do 30 k=1,nu
	sum2=sum2+a(k)
	mu=k
	if(sum2.ge.r2at)then
c	write(*,*)mu
	goto 100
	end if
 30	continue
c
c
c	Updating the number of molecules for each reaction
 100      goto(101,102,103,104,105,106,107,108,109,110,111,112,113,
     $         114,115,116,117,118,119,120,121,122,123,124,125,
     $         126,127,128,129,130,131,132,133,134,135,136,137,
     $         138,139,140,141,142,143,144,145,146,147,148,149,
     $         150,151,152,153,154,155,156,157,158,159,160,161,
     $         162,163,164,165,166,167,168,169,170,171,172,173,
     $         174,175,176,177,178,179,180,181,182,183,184,185,
     $         186,187,188,189,190,191,192,193,194,195,196,
     $         197,198,199,200,201,202,203,204,205,206,207,
     $         208,209,210,211,212,213,214,215,216,217,218,
     $         219,220,221,222,223,224,225,226,227,228,229),mu
c

 101	nMV=nMV+1
	goto 40
 102	nMV=nMV-1
	goto 40

 103	nMycm=nMycm+1
	goto 40
 104	nMycm=nMycm+1
	goto 40
 105	nMycm=nMycm+1
	goto 40
 106	nMycm=nMycm-1
	goto 40

 107	nMycp=nMycp+1
	goto 40
 108	nMycp=nMycp-1
	goto 40
 109	nMycp=nMycp-1
	goto 40

 110	nE2Fm=nE2Fm+1
	goto 40
 111	nE2Fm=nE2Fm+1
	goto 40
 112	nE2Fm=nE2Fm+1
	goto 40
 113	nE2Fm=nE2Fm-1
	goto 40

 114	nE2Fp=nE2Fp+1
	goto 40
 115	nE2Fp=nE2Fp-1
	goto 40
 116	nE2Fp=nE2Fp-1
	goto 40

 117	nDE=nDE+1
	nE2Fp=nE2Fp-1
 	goto 40
 118	nDE=nDE-1
	nE2Fp=nE2Fp+1
	goto 40
 119	nDE=nDE-1
	goto 40
 120	nDE=nDE-1
	goto 40

 121	nRbm=nRbm+1
	goto 40
 122	nRbm=nRbm-1
	goto 40
 123	nRb=nRb+1
	goto 40
 124	nRb=nRb-1
	goto 40
 125	niRb1=niRb1-1
	goto 40
 126	niRb2=niRb2-1
	goto 40

 127	nDE=nDE-1
	nRb=nRb-1
	niDE=niDE+1
	goto 40
 128	nDE=nDE+1
	nRb=nRb+1
	niDE=niDE-1
	goto 40
 129	niDE=niDE-1
	goto 40

 130	nDE=nDE-1
	niRb1=niRb1-1
	niDEP=niDEP+1
	goto 40
 131	nDE=nDE+1
	niRb1=niRb1+1
	niDEP=niDEP-1
	goto 40
 132	niDEP=niDEP+1
	niDE=niDE-1
	goto 40
 133	niDEP=niDEP-1
	niDE=niDE+1
	goto 40
 134	niDEP=niDEP-1
	niRb2=niRb2+1
	nDE=nDE+1
	goto 40

 135	niRb1=niRb1+1
	nRb=nRb-1
	goto 40
 136	niRb1=niRb1-1
	nRb=nRb+1
	goto 40
 137	niRb1=niRb1-1
	niRb2=niRb2+1
	goto 40
 138	niRb1=niRb1+1
	niRb2=niRb2-1
	goto 40

 139	nCycDm=nCycDm+1
	goto 40
 140	nCycDm=nCycDm+1
	goto 40
 141	nCycDm=nCycDm-1
	goto 40

 142	nCycD=nCycD+1
	goto 40
 143	nCycD=nCycD-1
	goto 40
 144	nCycD=nCycD-1
	goto 40

 145	nCycEm=nCycEm+1
	goto 40
 146	nCycEm=nCycEm+1
	goto 40
 147	nCycEm=nCycEm-1
	goto 40

 148	nCycE=nCycE+1
	goto 40
 149	nCycE=nCycE-1
	goto 40
 150	nCycE=nCycE-1
	goto 40


 151	nIm=nIm+1
	goto 40
 152	nIm=nIm-1
	goto 40

 153	nI=nI+1
	goto 40
 154	nI=nI-1
	goto 40

 155	nI=nI-1
	nIp=nIp+1
	goto 40
 156	nIp=nIp-1
	goto 40
 157	nIp=nIp-1
	goto 40

 158	nXI=nXI+1
	nCycE=nCycE-1
	nI=nI-1
	goto 40
 159	nXI=nXI-1
	nCycE=nCycE+1
	nI=nI+1
	goto 40
 160	nXI=nXI-1
	goto 40
 161	nXI=nXI-1
	goto 40
 162	nXI=nXI-1
	nCycE=nCycE+1
	goto 40

 163	nSkp2m=nSkp2m+1
	goto 40
 164	nSkp2m=nSkp2m-1
	goto 40

 165	nSkp2i=nSkp2i+1
	goto 40
 166	nSkp2i=nSkp2i-1
	goto 40
 167	nSkp2i=nSkp2i-1
	goto 40
 168	nSkp2a=nSkp2a+1
	nSkp2i=nSkp2i-1
	goto 40
 169	nSkp2i=nSkp2i+1
	nSkp2a=nSkp2a-1
	goto 40 
 170	nSkp2a=nSkp2a-1
	goto 40
 171	nSkp2a=nSkp2a-1
	goto 40

 172	nCycAm=nCycAm+1
	goto 40
 173	nCycAm=nCycAm+1
	goto 40
 174	nCycAm=nCycAm-1
	goto 40

 175	nCycA=nCycA+1
	goto 40
 176	nCycA=nCycA-1
	goto 40
 177	nCycA=nCycA-1
	goto 40

 178	nCycBm=nCycBm+1
	goto 40
 179	nCycBm=nCycBm-1
	goto 40

 180	nCycBa=nCycBa+1
	goto 40
 181	nCycBa=nCycBa-1
	goto 40
 182	nCycBa=nCycBa-1
	goto 40
 183	nCycBi=nCycBi+1
	nCycBa=nCycBa-1
	goto 40
 184	nCycBa=nCycBa+1
	nCycBi=nCycBi-1
	goto 40
 185	nCycBi=nCycBi-1
	goto 40
 186	nCycBi=nCycBi-1
	goto 40

 187	nC25m=nC25m+1
	goto 40
 188	nC25m=nC25m-1
	goto 40

 189	nC25i=nC25i+1
	goto 40
 190	nC25i=nC25i-1
	goto 40
 191	nC25a=nC25a-1
	goto 40
 192	nC25i=nC25i-1
	nC25a=nC25a+1
	goto 40
 193	nC25i=nC25i+1
	nC25a=nC25a-1
	goto 40

 194	nWee1m=nWee1m+1
	goto 40
 195	nWee1m=nWee1m-1
	goto 40

 196	nWee1a=nWee1a+1
	goto 40
 197	nWee1a=nWee1a-1
	goto 40
 198	nWee1i=nWee1i-1
	nWee1a=nWee1a+1
	goto 40
 199	nWee1i=nWee1i+1
	nWee1a=nWee1a-1
	goto 40
 200	nWee1i=nWee1i-1
	goto 40

 201	nCdh1m=nCdh1m+1
	goto 40
 202	nCdh1m=nCdh1m-1
	goto 40

 203	nCdh1i=nCdh1i+1
	goto 40
 204	nCdh1i=nCdh1i-1
	goto 40
 205	nCdh1a=nCdh1a+1
	nCdh1i=nCdh1i-1
	goto 40
 206	nCdh1a=nCdh1a+1
	nCdh1i=nCdh1i-1
	goto 40
 207	nCdh1a=nCdh1a-1
	nCdh1i=nCdh1i+1
	goto 40
 208	nCdh1a=nCdh1a-1
	nCdh1i=nCdh1i+1
	goto 40
 209	nCdh1a=nCdh1a-1
	goto 40

 210	nC20m=nC20m+1
	goto 40
 211	nC20m=nC20m+1
	goto 40
 212	nC20m=nC20m-1
	goto 40

 213	nCdc20i=nCdc20i+1
	goto 40
 214	nCdc20i=nCdc20i-1
	goto 40
 215	nCdc20a=nCdc20a+1
	nCdc20i=nCdc20i-1
	goto 40
 216	nCdc20i=nCdc20i+1
	nCdc20a=nCdc20a-1
	goto 40
 217	nCdc20a=nCdc20a-1
	goto 40

 218	nIEP=nIEP+1
	nIEPi=nIEPT-nIEP
	goto 40
 219	nIEP=nIEP-1
	nIEPi=nIEPT-nIEP
	goto 40

 220	nCdt1m=nCdt1m+1
	goto 40
 221	nCdt1m=nCdt1m-1
	goto 40
 222	nCdt1=nCdt1+1
	goto 40
 223	nCdt1=nCdt1-1
	goto 40
 224	nCdt1=nCdt1-1
	goto 40

 225	nGemm=nGemm+1
	goto 40
 226	nGemm=nGemm-1
	goto 40
 227	nGem=nGem+1
	goto 40
 228	nGem=nGem-1
	goto 40
 229	nGem=nGem-1
	goto 40


 40	nY=nE2Fp+nDE+niDE+niDEP	
	nDp1p=nDp1t-nDE-niDE-niDEP
	nRbt=nRb+niRb1+niRb2+niDE+niDEP
	nIEPi=nIEPT-nIEP
	nSkp2T=nSkp2a+nSkp2i
	nCdc20T=nCdc20a+nCdc20a
	nCdh1T=nCdh1i+nCdh1a

    	nCycBT=nCycBi+nCycBa
    	nWee1T=nWee1a+nWee1i
    	nC25T=nC25a+nC25i

	nrxn=nrxn+1
	ti=ti+tau
c     To avoid nrxn becoming a large value

	
        if((nrxn.eq.imaximum))then
		nrxn=0
	if((jk.eq.kl).and.(ii.ge.ki).and.(jj.eq.nj))then
		n=n+1
		if(ti.gt.0)then
c		write(301,*)ti,V,jk,ii,jj
		write(302,*)ti,dGF,jk,ii,jj

		write(201,*)ti,nCycA,jk,ii,jj
c		write(202,*)ti,nCdh1i,jk,ii,jj
		write(203,*)ti,nCdh1a,jk,ii,jj
c		write(204,*)ti,nCdc20i,jk,ii,jj
		write(205,*)ti,nCdc20a,jk,ii,jj
c		write(206,*)ti,nIEPi,jk,ii,jj
c		write(207,*)ti,nIEP,jk,ii,jj
c		write(208,*)ti,nRb,jk,ii,jj
		write(209,*)ti,nCdh1m,jk,ii,jj
c		write(210,*)ti,nC20m,jk,ii,jj
c
c		write(211,*)ti,nMycm,jk,ii,jj
c		write(212,*)ti,nMycp,jk,ii,jj
c		write(213,*)ti,nE2Fm,jk,ii,jj
c		write(214,*)ti,nY,jk,ii,jj
c		write(215,*)ti,nDE,jk,ii,jj
c		write(216,*)ti,niDE,jk,ii,jj
c		write(217,*)ti,nCycDm,jk,ii,jj
c		write(218,*)ti,nCycD,jk,ii,jj
c		write(219,*)ti,nCycEm,jk,ii,jj
c		write(220,*)ti,nCycE,jk,ii,jj
c		write(221,*)ti,niRb1,jk,ii,jj
c		write(222,*)ti,niRb2,jk,ii,jj
c
c		write(223,*)ti,nSkp2a,jk,ii,jj
c		write(224,*)ti,nSkp2i,jk,ii,jj
		write(225,*)ti,nCdt1,jk,ii,jj
		write(226,*)ti,nGem,jk,ii,jj
c		write(227,*)ti,nCdt1m,jk,ii,jj
c		write(228,*)ti,nGemm,jk,ii,jj
c		write(229,*)ti,nSkp2m,jk,ii,jj
c		write(230,*)ti,nIm,jk,ii,jj
c		write(231,*)ti,nI,jk,ii,jj
c		write(232,*)ti,nIp,jk,ii,jj
c		write(233,*)ti,nXI,jk,ii,jj
c		write(234,*)ti,nCycAm,jk,ii,jj
c		write(235,*)ti,nRbm,jk,ii,jj
c		write(236,*)ti,nRbT,jk,ii,jj
c		write(237,*)ti,nMV,jk,ii,jj	
c		write(238,*)ti,niDEP,jk,ii,jj	


		write(239,*)ti,nCycBm,jk,ii,jj
		write(240,*)ti,nCycBa,jk,ii,jj	
c		write(241,*)ti,nC25m,jk,ii,jj
		write(242,*)ti,nC25a,jk,ii,jj
c		write(243,*)ti,nWee1m,jk,ii,jj	
		write(244,*)ti,nWee1a,jk,ii,jj



        endif
	end if
	end if

c	time of cell division and no. of molecules after division is stored in array for a cell;
c	This will be used as initial conditions for daughter cells


c	#############################################################################
c	CALCULATION OF GEMININ PEAK
c	#############################################################################
c	if(ti.gt.100.and.ti.lt.200)then	
	if(nGem.gt.nCdt1.and.nGem.gt.maxg)then
		maxg=nGem
		tgg=ti
	end if
	if((nCdh1a.gt.(3000)).and.(nCycBa.le.(1000)))then
c	if((nCdt1.gt.300).and.(nGem.le.200))then
c		write(*,*),ti
		if((ncount.eq.1))then
c			write(*,*)nCdt1,nGem,nCycA,ncount
			mm(ncell)=maxg
			tmm(ncell)=tgg
			tem(ncell)=ti
			em(ncell)=nGem

c 			write(*,*)'ncell=',ncell
c			write(*,*)tmm(ncell),mm(ncell)
c			write(*,*)tem(ncell),em(ncell)

c	writing duration in output file
			if(ii.ge.1) then	
c			if((jk.eq.2).and.(jj.eq.1))then
			tcc(ncell)=ti-tid(ij)
			g1(ncell)=teg(ncell)-tid(ij)
			sg2m(ncell)=tcc(ncell)-g1(ncell)

c			write(247,*)tid(ij),teg(ncell),ti,jk,ii,jj
			write(248,*)tcc(ncell),g1(ncell),sg2m(ncell),jk,ii,jj

			end if
			if(ii.eq.2)then
				if(mod(jj,2).eq.0)then
c				write(250,*)tcc(ncell),g1(ncell),sg2m(ncell),jk,ii,jj
				else
c				write(249,*)tcc(ncell),g1(ncell),sg2m(ncell),jk,ii,jj
				end if
			end if

	write(602,*)jk,ii,jj,bmycn,dkMcm1n,dkMcm2n,
     $         be2fn,dkEm1n,dkEm2n,
     $	       dkrbmn,dkCDm1n,dkCDm2n,
     $         bcycemn,dksemn,dksimn,
     $	       ak17mn,bsamn,dksamn,
     $         ak3mn,ak5amn,ak5bmn,
     $	       ak21mn,ak20mn,dkmv
	
c
 			ncell=ncell+1
			ncount=0
			mcount=1
			maxg=0

c		write(*,*)'cell division'
			goto 500		

c    Unequal division 

 500            r3 = duni()
               r4 = duni()
               api = 4.0*atan(1.0)

               if(r3.le.0.0)then
                  r3 = 0.00001
               endif
               
                ba  = -2.0*sigma*sigma*log(r3)
                ba1 = dcos(2.0*api*r4)
                sf  = sqrt(ba)*ba1
                sac = 0.5 + sf
		

c     If sigma=0.0 then sf=0.0, this is equal division (sac=0.5)
c     If sigma=0.025 then sf not equal to zero, this is unequal division 
                 if(sac.le.0.1.or.sac.ge.0.9)then
                    go to 500
                 else
		  sac1=sac
	
		  idiva=idiva+2
		  idiv=idiv+1
		sizea=sac
		sizeb=1.0-sac


c		write(*,*)sac1,idiva,idivb
c		daughter 1, sac1

		GFi(idiva)=dGF
		inMV(idiva)=nMV*sac1
		inMycm(idiva)=nMycm*sac1
		inMycp(idiva)=nMycp*sac1
		inE2Fm(idiva)=nE2Fm*sac1
		inE2Fp(idiva)=nE2Fp*sac1
		inDE(idiva)=nDE*sac1
		iniDE(idiva)=niDE*sac1
		iniDEP(idiva)=niDEP*sac1
		inRbm(idiva)=nRbm*sac1
		inRb(idiva)=nRb*sac1
		iniRb1(idiva)=niRb1*sac1
		iniRb2(idiva)=niRb2*sac1
		inCycDm(idiva)=nCycDm*sac1
		inCycD(idiva)=nCycD*sac1
		inCycEm(idiva)=nCycEm*sac1
		inCycE(idiva)=nCycE*sac1
		inI(idiva)=nI*sac1
		inIp(idiva)=nIp*sac1
		inXI(idiva)=nXI*sac1
		inSkp2a(idiva)=nSkp2a*sac1
		inSkp2i(idiva)=nSkp2i*sac1
		inCdc20a(idiva)=nCdc20a*sac1
		inCdc20i(idiva)=nCdc20i*sac1
		inCdh1a(idiva)=nCdh1a*sac1
		inCdh1i(idiva)=nCdh1i*sac1
c		inIEPi(idiva)=nIEPi*sac1
		inIEP(idiva)=nIEP*sac1
		inIm(idiva)=nIm*sac1
		inSkp2m(idiva)=nSkp2m*sac1
		inCdt1m(idiva)=nCdt1m*sac1
		inGemm(idiva)=nGemm*sac1
		inCdt1(idiva)=nCdt1*sac1
		inGem(idiva)=nGem*sac1
		inCycAm(idiva)=nCycAm*sac1
		inCycA(idiva)=nCycA*sac1
		inCdh1m(idiva)=nCdh1m*sac1
		inC20m(idiva)=nC20m*sac1	


		inCycBm(idiva)=nCycBm*sac1
		inCycBa(idiva)=nCycBa*sac1
		inCycBi(idiva)=nCycBi*sac1

		inC25m(idiva)=nC25m*sac1
		inC25a(idiva)=nC25a*sac1
		inC25i(idiva)=nC25i*sac1

		inWee1m(idiva)=nWee1m*sac1
		inWee1i(idiva)=nWee1i*sac1
		inWee1a(idiva)=nWee1a*sac1
c
		tid(idiva)=ti
            
        do 81 jr=1,ku

            akc(jr,idiva)=akm(jr)

 81       continue


c	
c		daughter 2
		idivb=idivb+2
		idiv=idiv+1



		GFi(idivb)=dGF

		inMV(idivb)=nMV-inMV(idiva)
		inMycm(idivb)=nMycm-inMycm(idiva)
		inMycp(idivb)=nMycp-inMycp(idiva)
		inE2Fm(idivb)=nE2Fm-inE2Fm(idiva)
		inE2Fp(idivb)=nE2Fp-inE2Fp(idiva)
		inDE(idivb)=nDE-inDE(idiva)
		iniDE(idivb)=niDE-iniDE(idiva)
		iniDEP(idivb)=niDEP-iniDEP(idiva)
		inRbm(idivb)=nRbm-inRbm(idiva)
		inRb(idivb)=nRb-inRb(idiva)
		iniRb1(idivb)=niRb1-iniRb1(idiva)
		iniRb2(idivb)=niRb2-iniRb2(idiva)
		inCycDm(idivb)=nCycDm-inCycDm(idiva)
		inCycD(idivb)=nCycD-inCycD(idiva)
		inCycEm(idivb)=nCycEm-inCycEm(idiva)
		inCycE(idivb)=nCycE-inCycE(idiva)
		inI(idivb)=nI-inI(idiva)
		inIp(idivb)=nIp-inIp(idiva)
		inXI(idivb)=nXI-inXI(idiva)
		inSkp2a(idivb)=nSkp2a-inSkp2a(idiva)
		inSkp2i(idivb)=nSkp2i-inSkp2i(idiva)
		inCdc20a(idivb)=nCdc20a-inCdc20a(idiva)
		inCdc20i(idivb)=nCdc20i-inCdc20i(idiva)
		inCdh1a(idivb)=nCdh1a-inCdh1a(idiva)
		inCdh1i(idivb)=nCdh1i-inCdh1i(idiva)
c		inIEPi(idivb)=nIEPi-inIEPi(idiva)
		inIEP(idivb)=nIEP-inIEP(idiva)
		inIm(idivb)=nIm-inIm(idiva)
		inSkp2m(idivb)=nSkp2m-inSkp2m(idiva)
		inCdt1m(idivb)=nCdt1m-inCdt1m(idiva)
		inGemm(idivb)=nGemm-inGemm(idiva)
		inCdt1(idivb)=nCdt1-inCdt1(idiva)
		inGem(idivb)=nGem-inGem(idiva)
		inCycAm(idivb)=nCycAm-inCycAm(idiva)
		inCycA(idivb)=nCycA-inCycA(idiva)
		inCdh1m(idivb)=nCdh1m-inCdh1m(idiva)
		inC20m(idivb)=nC20m-inC20m(idiva)	

		inCycBm(idivb)=nCycBm-inCycBm(idiva)
		inCycBa(idivb)=nCycBa-inCycBa(idiva)
		inCycBi(idivb)=nCycBi-inCycBi(idiva)

		inC25m(idivb)=nC25m-inC25m(idiva)
		inC25a(idivb)=nC25a-inC25a(idiva)
		inC25i(idivb)=nC25i-inC25i(idiva)

		inWee1m(idivb)=nWee1m-inWee1m(idiva)
		inWee1i(idivb)=nWee1i-inWee1i(idiva)
		inWee1a(idivb)=nWee1a-inWee1a(idiva)


c
		tid(idivb)=ti

        do 82 jr=1,ku

            akc(jr,idivb)=akm(jr)

 82       continue

			icdiv=0	
			goto 65	
	end if
c	end if

		else
			ncount=0
			maxg=0
		end if
	end if
c	#############################################################################
c	CALCULATION OF Cdt1 PEAK AND CUTOFF
c	#############################################################################
	if(nCdt1.gt.nGem.and.(nCdt1.gt.maxc))then
		maxc=nCdt1
		tc=ti
	end if
	if((nGem.gt.(4000)).and.(nCdt1.le.(2000)))then
		if((mcount.eq.1).and.(nCycA.gt.(2000)))then
			eg(ncell)=maxc
			teg(ncell)=tc
c			write(*,*)teg(ncell),eg(ncell)
			mcount=0
			ncount=1
			tcount=1
		write(603,*)jk,ii,jj,ti,tmit,trm,bmycn*trm,'1'
			

		akn=rates(akm,ku)

		do  jr=1,ku

			akm(jr)=akn(jr)
			akc(jr,ij)=akn(jr)

		end do

		write(604,*)jk,ii,jj
		write(604,*)(akm(jr),jr=1,ku)


		write(603,*)jk,ii,jj,ti,tmit,trm,akm(1)*trm,'2'


c		write(*,*)ncount,maxc
			maxc=0
		else
			mcount=0
			maxc=0
		end if
	end if

	if((nCycBa.gt.(4000)).and.(tcount.eq.1))then
		
		sigM=0.0

		r5 = duni()
               r6 = duni()
               api = 4.0*atan(1.0)

               if(r5.le.0.0)then
                  r5 = 0.00001
               endif
               
                ba  = -2.0*sigM*sigM*log(r5)
                ba1 = dcos(2.0*api*r6)
                sf  = sqrt(ba)*ba1	
                sacm = 1+ sf	

		tmit=ti+sacm
		
	
c		write(*,*)ti,sacm,tmit,ii,jj,jk
	
			tcount=0
			
	end if


c	end if
c
c	#############################################################################
c	#############################################################################
c


c	criteria to check end of time
 	if(ti.lt.tf)then
		goto 10

c	statement to give end time and no. of molecules at the end time for the cell
	else
	
 		kdiv=1
		idiva=idiva+2
  		idiv=idiva
 66		sac1=1


		GFi(idiv)=dGF
		inMV(idiv)=nMV*sac1
		inMycm(idiv)=nMycm*sac1
		inMycp(idiv)=nMycp*sac1
		inE2Fm(idiv)=nE2Fm*sac1
		inE2Fp(idiv)=nE2Fp*sac1
		inDE(idiv)=nDE*sac1
		iniDE(idiv)=niDE*sac1
		iniDEP(idiv)=niDEP*sac1
		inRbm(idiv)=nRbm*sac1
		inRb(idiv)=nRb*sac1
		iniRb1(idiv)=niRb1*sac1
		iniRb2(idiv)=niRb2*sac1
		inCycDm(idiv)=nCycDm*sac1
		inCycD(idiv)=nCycD*sac1
		inCycEm(idiv)=nCycEm*sac1
		inCycE(idiv)=nCycE*sac1
		inI(idiv)=nI*sac1
		inIp(idiv)=nIp*sac1
		inXI(idiv)=nXI*sac1
		inSkp2a(idiv)=nSkp2a*sac1
		inSkp2i(idiv)=nSkp2i*sac1
		inCdc20a(idiv)=nCdc20a*sac1
		inCdc20i(idiv)=nCdc20i*sac1
		inCdh1a(idiv)=nCdh1a*sac1
		inCdh1i(idiv)=nCdh1i*sac1
c		inIEPi(idiv)=nIEPi*sac1
		inIEP(idiv)=nIEP*sac1
		inIm(idiv)=nIm*sac1
		inSkp2m(idiv)=nSkp2m*sac1
		inCdt1m(idiv)=nCdt1m*sac1
		inGemm(idiv)=nGemm*sac1
		inCdt1(idiv)=nCdt1*sac1
		inGem(idiv)=nGem*sac1
		inCycAm(idiv)=nCycAm*sac1
		inCycA(idiv)=nCycA*sac1
		inCdh1m(idiv)=nCdh1m*sac1
		inC20m(idiv)=nC20m*sac1

		inCycBm(idiv)=nCycBm*sac1
		inCycBa(idiv)=nCycBa*sac1
		inCycBi(idiv)=nCycBi*sac1

		inC25m(idiv)=nC25m*sac1
		inC25a(idiv)=nC25a*sac1
		inC25i(idiv)=nC25i*sac1

		inWee1m(idiv)=nWee1m*sac1
		inWee1i(idiv)=nWee1i*sac1
		inWee1a(idiv)=nWee1a*sac1
c
		tid(idiv)=ti

        do 83 jr=1,ku

            akc(jr,idiv)=akm(jr)
	   

 83       continue


		if(kdiv.eq.2) then
		 itend=0
		  kdiv=1
		  goto 65
		else
			kdiv=2
			idivb=idivb+2
			idiv=idivb
		  goto 66
		endif 
	end if
c	statement to write cell division, or end of time for a partiular cell
 65	if(icdiv.eq.0) then
		write(*,*)jk,ii,jj,'start time:',tid(ij),'cell division time:',ti
		icdiv=1

c	for the cell at which time ended for first time
	else if(itend.eq.0)then
		write(*,*)ii,jj,'start time:',tid(ij),'simulation end time',ti
		itend=1

c	for subsequent cycles after time ended, end time will be stored as initial time point for the cells
	else
		idiva=idiva+2
		idiv=idiva
		tid(idiv)=ti
		GFi(idiv)=dGF

		idivb=idivb+2
		idiv=idivb
		tid(idiv)=ti
		GFi(idiv)=dGF

c		write(*,*)ii,jj,tid(ij),'end of time',ti
	endif

 80 	continue
 90 	continue
 300	continue
	stop 
	end



	REAL*8 FUNCTION rates(akm,ku) result(akn)

	implicit real*8 (a-h,q-z)
	dimension akm(ku),akn(ku)


	
	do 61 jr=1,ku

		akn(jr)=akm(jr)
 		rk3 = duni()

        if (rk3.le.0.5)then
		if (rk3.ge.0.4)then
			akn(jr)=akn(jr)-(0.06*akn(jr))
			
		else if (rk3.ge.0.3)then
			akn(jr)=akn(jr)-(0.12*akn(jr))
			

		else if (rk3.ge.0.2)then
			akn(jr)=akn(jr)-(0.18*akn(jr))

		else if (rk3.ge.0.1)then
			akn(jr)=akn(jr)-(0.24*akn(jr))
			
		else if (rk3.ge.0.0)then
			akn(jr)=akn(jr)-(0.30*akn(jr))
			
		end if
	else if (rk3.ge.0.5)then

		if (rk3.le.0.6)then
			akn(jr)=akn(jr)+(0.06*akn(jr))
			
		else if (rk3.le.0.7)then
			akn(jr)=akn(jr)+(0.12*akn(jr))
		

		else if (rk3.le.0.8)then
			akn(jr)=akn(jr)+(0.18*akn(jr))
		
		else if (rk3.le.0.9)then
			akn(jr)=akn(jr)+(0.24*akn(jr))
		
		else if (rk3.le.1.0)then
			akn(jr)=akn(jr)+(0.30*akn(jr))
		
		end if

	end if


 61	continue



	end



c..........................................................................
c      Subroutine for random no generator by Lane Watson
c.........................................................................
            DOUBLE PRECISION FUNCTION DUNI()
c            FUNCTION DUNI()
C***BEGIN PROLOGUE  DUNI
C***DATE WRITTEN   880714 (YYMMDD)
C***REVISION DATE  880714 (YYMMDD)
C***CATEGORY NO.  L6A21
C***KEYWORDS  RANDOM NUMBERS, UNIFORM RANDOM NUMBERS
C***AUTHOR    KAHANER, DAVID, SCIENTIFIC COMPUTING DIVISION, NBS
C             MARSAGLIA, GEORGE, SUPERCOMPUTER RES. INST., FLORIDA ST. U.
C
C***PURPOSE  THIS ROUTINE GENERATES DOUBLE PRECISION UNIFORM
C             RANDOM NUMBERS ON [0,1)
C***DESCRIPTION
C        COMPUTES DOUBLE PRECISION UNIFORM NUMBERS ON [0,1).
C           FROM THE BOOK, "NUMERICAL METHODS AND SOFTWARE" BY
C                D. KAHANER, C. MOLER, S. NASH
C                PRENTICE HALL, 1988
C
C       USAGE:
C              TO INITIALIZE THE GENERATOR
C                   USEED = DUSTAR(ISEED)
C               WHERE: ISEED IS ANY NONZERO INTEGER
C                  WILL RETURN FLOATING POINT VALUE OF ISEED.
C
C               SUBSEQUENTLY
C                       U = DUNI()
C                  WILL RETURN A REAL UNIFORM ON [0,1)
C
C                ONE INITIALIZATION IS NECESSARY, BUT ANY NUMBER OF EVALUATIONS
C                  OF DUNI IN ANY ORDER, ARE ALLOWED.
C
C           NOTE: DEPENDING UPON THE VALUE OF K (SEE BELOW), THE OUTPUT
C                       OF DUNI MAY DIFFER FROM ONE MACHINE TO ANOTHER.
C
C           TYPICAL USAGE:
C
C               DOUBLE PRECISION U,DUNI,DUSTAR,USEED
C               INTEGER ISEED
CC                 SET SEED
C               ISEED = 305
C               USEED = DUSTAR(ISEED)
C               DO 1 I = 1,1000
C                   U = DUNI()
C             1 CONTINUE
CC                 NOTE: IF K=47 (THE DEFAULT, SEE BELOW) THE OUTPUT VALUE OF
CC                           U WILL BE 0.812053811384E-01...
C               WRITE(*,*) U
C               END
C
C          NOTE ON PORTABILITY: USERS CAN CHOOSE TO RUN DUNI IN ITS DEFAULT
C               MODE (REQUIRING NO USER ACTION) WHICH WILL GENERATE THE SAME
C               SEQUENCE OF NUMBERS ON ANY COMPUTER SUPPORTING FLOATING POINT
C               NUMBERS WITH AT LEAST 47 BIT MANTISSAS, OR IN A MODE THAT
C               WILL GENERATE NUMBERS WITH A LONGER PERIOD ON COMPUTERS WITH
C               LARGER MANTISSAS.
C          TO EXERCISE THIS OPTION:  B E F O R E  INVOKING DUSTAR INSERT
C               THE INSTRUCTION        UBITS = DUNIB(K)      K >= 47
C               WHERE K IS THE NUMBER OF BITS IN THE MANTISSA OF YOUR FLOATING

C               POINT WORD (K=96 FOR CRAY, CYBER 205). DUNIB RETURNS THE
C               FLOATING POINT VALUE OF K THAT IT ACTUALLY USED.
C                    K INPUT AS .LE. 47, THEN UBITS=47.
C                    K INPUT AS .GT. 47, THEN UBITS=FLOAT(K)
C               IF K>47 THE SEQUENCE OF NUMBERS GENERATED BY DUNI MAY DIFFER
C               FROM ONE COMPUTER TO ANOTHER.
C
C
C***REFERENCES  MARSAGLIA G., "COMMENTS ON THE PERFECT UNIFORM RANDOM
C                 NUMBER GENERATOR", UNPUBLISHED NOTES, WASH S. U.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE DUNI
      DOUBLE PRECISION CSAVE,CD,CM
      PARAMETER(
     *    CSAVE=0.9162596898123D+13/0.140737488355328D+15,
     *    CD=0.76543212345678D+14/0.140737488355328D+15,
     *    CM=0.140737488355213D+15/0.140737488355328D+15)
C                            2**47=0.140737488355328D+15
      DOUBLE PRECISION U(17),S,T,DUSTAR,C,DUNIB
      INTEGER I,J,II,JJ,K,KK,I1,J1,K1,L1,M1,ISEED
C
      SAVE U,I,J,K,C
C      LOAD DATA ARRAY IN CASE USER FORGETS TO INITIALIZE.
C      THIS ARRAY IS THE RESULT OF CALLING DUNI 100000 TIMES
C         WITH ISEED=305 AND K=96.
      DATA U/
     *0.471960981577884755837789724978D+00,
     *0.930323453205669578433639632431D+00,
     *0.110161790933730836587127944899D+00,
     *0.571501996273139518362638757010D-01,
     *0.402467554779738266237538503137D+00,
     *0.451181953427459489458279456915D+00,
     *0.296076152342721102174129954053D+00,
     *0.128202189325888116466879622359D-01,
     *0.314274693850973603980853259266D+00,
     *0.335521366752294932468163594171D-02,
     *0.488685045200439371607850367840D+00,
     *0.195470426865656758693860613516D+00,
     *0.864162706791773556901599326053D+00,
     *0.335505955815259203596381170316D+00,
     *0.377190200199058085469526470541D+00,
     *0.400780392114818314671676525916D+00,
     *0.374224214182207466262750307281D+00/
      DATA I,J,K,C/17,5,47,CSAVE/
C
C   BASIC GENERATOR IS FIBONACCI
C
      DUNI = U(I)-U(J)
      IF(DUNI.LT.0.0D0)DUNI = DUNI+1.0D0
      U(I) = DUNI
      I = I-1
      IF(I.EQ.0)I = 17
      J = J-1
      IF(J.EQ.0)J = 17
C
C   SECOND GENERATOR IS CONGRUENTIAL
C
      C = C-CD
      IF(C.LT.0.0D0) C=C+CM
C
C   COMBINATION GENERATOR
C
      DUNI = DUNI-C
      IF(DUNI.LT.0.0D0)DUNI = DUNI+1.0D0
      RETURN
C
      ENTRY DUSTAR(ISEED)
C
C          SET UP ...
C          CONVERT ISEED TO FOUR SMALLISH POSITIVE INTEGERS.
C
        I1 = MOD(ABS(ISEED),177)+1
        J1 = MOD(ABS(ISEED),167)+1
        K1 = MOD(ABS(ISEED),157)+1
        L1 = MOD(ABS(ISEED),147)+1
C
C              GENERATE RANDOM BIT PATTERN IN ARRAY BASED ON GIVEN SEED.
C
        DO 2 II = 1,17
          S = 0.0D0
          T = 0.5D0
C             DO FOR EACH OF THE BITS OF MANTISSA OF WORD
C             LOOP  OVER K BITS, WHERE K IS DEFAULTED TO 47 BUT CAN
C               BE CHANGED BY USER CALL TO DUNIB(K)
          DO 3 JJ = 1,K
                  M1 = MOD(MOD(I1*J1,179)*K1,179)
                  I1 = J1
                  J1 = K1
                  K1 = M1
                  L1 = MOD(53*L1+1,169)
                  IF(MOD(L1*M1,64).GE.32)S=S+T
                 T = 0.5D0*T
    3   continue
      U(II) = S
    2   continue
        DUSTAR = FLOAT(ISEED)
        RETURN
C
      ENTRY DUNIB(KK)
        IF(KK.LE.47)THEN
             K=47
        ELSE
             K=KK
        ENDIF
        DUNIB=FLOAT(K)
      END
c	##########################################################################

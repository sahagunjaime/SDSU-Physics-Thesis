!**************************************************************************************
!-----FINITE TEMPERATURE----FINITE TEMPERATURE----FINITE TEMPERATURE-----------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**************************************************************************************
!Remember to change the library : imsl.lib imsls_err.lib IMSLMPISTUB.LIB

    Implicit none

    Integer :: itmax,INTERV,poly,IRULE,l
    Integer :: title1,title2,title3

    Integer, Parameter :: N=7
    Integer, Parameter :: Nmax=14  !For the baryonic density
    Integer, Parameter :: Mmax=280 !For the Temperature

    Real*8 :: BOUND,ERRABS,ERREL, A, B,ERRABI, ERRELI, baryon, E
    Real*8 :: FNORM, xmu, barinf, barsup, dbar,F(N),ERREST,hbarc
    Real*8 :: X(N),alpha,beta,gamma,XGUESS(N),erro,val_ini,val_fin
    Real*8 :: alfi,beti,gami,pi2,N_c,xT,Lambda0,G_s,K_s,mou,mod,mos
    Real*8 :: Tinf, Tsup, dT,conda,condb,condg,p_0,xmue,me,mamu,mubi,eli
    Real*8 :: Dene,Ome,Dmuon,Omuon,rho_u,rho_d,rho_s,eletropy,mutropy
    Real*8 :: Equ,Eqd,Eqs,Eqtot,O_qu,O_qd,O_qs,Omtotal,xnuel,xnumu
    Real*8 :: numui, nueli, epsillon, pressure

    EXTERNAL FCN,DNEQNF,DQDAG,DQDAGI

	Common/Polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/par_INTEG/INTERV,IRULE
    Common/par2_Int/A,B,ERRABS,ERREL,BOUND,ERRABI,ERRELI,ERREST
	Common/presvacio/p_0
	Common/densibaryion/baryon
	Common/muon/mamu
	Common/masaelectron/me
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu
    Common/debug/counter
    !Common/teste/xnuel,xnumu

    !Data xnumu/0.d0/
    !Data xnuel/10.d0/
	Data me /0.511d0/           !electron mass in MeV
	Data mamu/105.658d0/        !muon mass
    Data hbarc/197.327d0/       !MeV.fm        !hbar * speed of light
    Data pi2/9.869604401d0/
!_____________________________________________________________________
!Parameters for NJL
    Lambda0=602.3d0                      !MeV
	mou=5.5d0    !0.d0                   !MeV 
	mod=5.5d0    !0.d0                   !MeV 
	mos=140.7d0  !0.d0                   !MeV 
    G_s=3.67d0/(Lambda0*Lambda0)         !1/MeV^2!Coupling Constant
    K_s=-12.36d0/(Lambda0)**5            !1/MeV^5!Coupling Constant
	p_0=32761280506.2342d0               !Vacuum Pressure
!_____________________________________________________________________
!Parameters used in HATSUDA and KUNIHIRO paper/report
!   Lambda0=631.4d0                     !MeV
!	mou=5.5d0           !0.d0           !MeV
!	mod=5.5d0           !0.d0           !MeV
!	mos=135.7d0         !0.d0           !MeV
!   G_s=3.67d0/(Lambda0*Lambda0)        !1/MeV^2!Constante de acoplamiento
!   K_s=-9.29d0/(Lambda0)**5            !1/MeV^5!Constante de acoplamiento
!_____________________________________________________________________
     
!+++++++++Parameters for semi-infinite integrals (DQDAGI)+++
!Set limits of integration
      BOUND  = 0.0
      INTERV = 1

!Set error tolerances for the integral
      ERRABI =0.0d0     !1.d-3
      ERRELI =1.d-3     !1.d-8
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++Parameters for finite integrals (DQDAG)+++++++++++++++++++++

!Set limits of integration
      A = 0.0d0
      B = Lambda0
!Set error tolerances for the integral
      ERRABS =0.0d0   !1.d-3   !  ! 
      ERREL =1.0d-4 !0.0d0      !1.d-8
!Choice of quadrature rule
!IRULE = 2 is recommended for most functions. If the function has a peak
!singularity, use IRULE = 1. If the function is oscillatory, use IRULE = 6.
      IRULE=2      
!____________________________________________________________________
        	
	title1=1
	title2=1
	title3=1
	!Poly=1 implies the Polyakov
	poly=2

	if(poly == 1)then
        N_c=1.d0           !Number of colors
	else
        N_c=3.d0
	endif

    alfi= (-241.9d0)**3    !(241.9d0)**3
	beti=(-241.9d0)**3     !(241.9d0)**3
	gami=(257.7d0)**3      !(257.7)**3
	mubi=450.d0            !450.d0 initial xmu guess
	eli=190.d0             !300.d0
	nueli=170.d0           !270.0d0
	numui=0.0d0            !5.8d0

!**************************************************!        
!           VARYING THE BARIONIC DENSITY           !            !0.16 is avg. nuclear density
!**************************************************!

    val_ini=1.d0*(0.16d0)*hbarc**3   !1.d0*(0.16d0)*hbarc**3    !initial density in units of rho_0 (1/frm^3)
	val_fin=15.d0*(0.16d0)*hbarc**3                             !final density?

    barinf=val_ini
    barsup=val_fin
    dbar=(barsup-barinf)/Nmax
  
    baryon=barinf


    Do l=0,Nmax

!**************************************************!        
!               VARYING THE TEMPERATURE            !
!**************************************************!

    Tinf=1.d0           !MeV    !Tempererature Variation
    Tsup=300.d0         !MeV
    dT=(Tsup-Tinf)/Mmax
  
    xT=Tinf

    !Do l=0,Mmax

      
	If(xT == Tinf.and. baryon == barinf)then
        xguess(1)=alfi
        xguess(2)=beti
        xguess(3)=gami
        xguess(4)=mubi
        xguess(5)=eli
        xguess(6)=nueli
        xguess(7)=numui
	else
        xguess(1)=alpha
        xguess(2)=beta
        xguess(3)=gamma
        xguess(4)=xmu
        xguess(5)=xmue
        xguess(6)=xnuel
        xguess(7)=xnumu
	endif
	
	Itmax=1000      !(maximum number of iterations for the non-linear system)

	erro=1.d-8   !1.d-10


	CALL DNEQNF (FCN,erro,N,itmax, xguess,x,fnorm) !This finds a solution for the 7 nonlinear equations
                                                   !Once it finds solution, it continues...
!****************EXIT**********************************************

	If(title1 == 1)then
        write(*,*) 'Hello'
	endif
    title1=2


	call Electron(x,Dene,Ome)
	call MUONS(x,Dmuon,Omuon)
	call Denquark(x,rho_u,rho_d,rho_s)
	call Lepentropy(x,eletropy,mutropy)
	call Omega_T(x,O_qu,O_qd,O_qs,Omtotal)
	call Entropyquark(x,Equ,Eqd,Eqs,Eqtot)
    call energydensity(x, epsillon)
    call pressure_sub(x, pressure)

!	WRITE(2,8)xmu,xT

	conda=(dabs(alpha))**(1.d0/3.d0)
	condb=(dabs(beta))**(1.d0/3.d0)
	condg=(dabs(gamma))**(1.d0/3.d0)


    write(*,*)''
    write(*,*)''
    write(*,*)''
    write(*,*)''
    write(*,*)''
    write(*,*)''
	write(*,*)'XXXXXXXXXXXXXXXXXXXX  RESULTS  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
    write(*,*)'baryon',baryon/((0.16d0)*hbarc**3),'xT=',xT, 'MeV'       !Number density is in units of rho_0
    write(*,*)''
    write(*,*)'alpha',conda,'beta',condb,'gamma',condg
    write(*,*)''
    write(*,*)'______________ CHEMICAL POTENTIAL (in MeV) ___________________'
    write(*,*)'mu_u=', xmu-2.d0*(xmue-xnuel)/3.d0,'mu_d = mu_s=',xmu+(xmue-xnuel)/3.d0
    write(*,*)'mu_quarks=',xmu,'mu_e=',xmue                                 !mu's are chemical potentials
	write(*,*)'mu_muon',xnumu,'mu_neutrinos',xnuel
    write(*,*)''
    write(*,*)'________________ GRAND POTENTIAL (in MeV) _____________________'
	write(*,*)'Om_u',O_qu/hbarc**3,'Om_d',O_qd/hbarc**3
	write(*,*)'Om_s',O_qs/hbarc**3,'Omtotal',Omtotal/hbarc**3
	write(*,*)''
	write(*,*)'Om_e',Ome/hbarc**3
	write(*,*)''
	write(*,*)'____________________ DENSITIES (MeV/fm3)_________________________'
	write(*,*)'rho_e',Dene/hbarc**3,'rho_muon',Dmuon/hbarc**3
	write(*,*)'rho_u',rho_u/(3.d0*(0.16d0)*hbarc**3)                        !quark density is in units of rho_0
	write(*,*)'rho_d',rho_d/(3.d0*(0.16d0)*hbarc**3)
	write(*,*)'rho_s',rho_s/(3.d0*(0.16d0)*hbarc**3)
	write(*,*)''
	write(*,*)'________________  ENTROPIES (MeV/fm3)  _______________________'
	write(*,*)'entrop_e',eletropy/hbarc**3,'entrop_muon',mutropy/hbarc**3
	write(*,*)''
	WRITE(*,*)'S_u=',Equ/hbarc**3,'S_d=',Eqd/hbarc**3,'S_s=',Eqs/hbarc**3
	write(*,*)''
	WRITE(*,*)'_________________ Dynamic Mass (MeV) __________________________'
	WRITE(*,*)'Mu',mou-2.d0*G_s*alpha-2.d0*K_s*beta*gamma
	write(*,*)'Md',mod-2.d0*G_s*beta-2.d0*K_s*alpha*gamma
    write(*,*)'Ms',mos-2.d0*G_s*gamma-2.d0*K_s*alpha*beta
	write(*,*)''
    WRITE(*,*)'____________________ Energy Density ___________________________'
    write(*,*)'Energy density', epsillon/hbarc**3, 'MeV'
	write(*,*)''
    WRITE(*,*)'_______________________ Pressure ______________________________'
    write(*,*)'Pressure', pressure/hbarc**3, 'MeV'
    write(*,*)''
    write(*,*)'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
    write(*,*)''
    write(*,*)''
    write(*,*)''
    write(*,*)''
    write(*,*)''
    write(*,*)''


    !xT=xT+dT
    baryon=baryon+dbar

    !End do !temperature loop
    End do  !number density loop
    Stop

	END
!****************************END MAIN*******************************!

!****************************END MAIN*******************************!

!****************************END MAIN*******************************!

!****************************END MAIN*******************************!

!****************************END MAIN*******************************!

!****************************END MAIN*******************************!

!*******************************************************************!
!           LOGRITHMIC INTEGRAL FOR THE u QUARK                     !
!*******************************************************************!

    FUNCTION Fun_logu(p)
    Implicit none

    Integer :: poly

    Real*8 :: M_u,xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Fun_logu,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,xmue,xmuu,xnuel,xnumu

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

    xmuu = xmu-2.d0*(xmue-xnuel)/3.d0


	M_u = mou-2.d0*G_s*alpha-2.d0*K_s*beta*gamma
	Ep = dsqrt(p**2+M_u**2)
	par1 = (Ep-xmuu)
	par2 = (Ep+xmuu)
	e1 = dexp(-par1/xT)             !dexp is Exp. with argument type D
	e2 = dexp(-par2/xT)


	n1=1.d0+e1
	n2=1.d0+e2


    Fun_logu=p**2*dlog(n1*n2)        !dlog is narutural log with argument type D

	END FUNCTION Fun_logu

!*******************************************************************!

!*******************************************************************!
!           LOGRITHMIC INTEGRAL FOR THE d QUARK                     !
!*******************************************************************!


    FUNCTION Fun_logd(p)
	Implicit none

    Integer :: poly

    Real*8 :: M_d,xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Fun_logd,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,xmue,xmud,xnuel,xnumu

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmud=xmu+(xmue-xnuel)/3.d0

	M_d=mod-2.d0*G_s*beta-2.d0*K_s*alpha*gamma
	Ep=dsqrt(p**2+M_d**2)
	par1=(Ep-xmud)
	par2=(Ep+xmud)
	e1=dexp(-par1/xT)
	e2=dexp(-par2/xT)


	n1=1.d0+e1
	n2=1.d0+e2


    Fun_logd=p**2*dlog(n1*n2)

	END FUNCTION Fun_logd

!*******************************************************************!

!*******************************************************************!
!           LOGRITHMIC INTEGRAL FOR THE d QUARK                     !
!*******************************************************************!

    FUNCTION Fun_logs(p)
	Implicit none

    Integer :: poly

    Real*8 :: M_s,xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Fun_logs,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,xmue,xmud,xnumu,xnuel

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmud=xmu+(xmue-xnuel)/3.d0



	M_s=mos-2.d0*G_s*gamma-2.d0*K_s*alpha*beta
	Ep=dsqrt(p**2+M_s**2)
	par1=(Ep-xmud)
	par2=(Ep+xmud)
	e1=dexp(-par1/xT)
	e2=dexp(-par2/xT)


	n1=1.d0+e1
	n2=1.d0+e2


    Fun_logs=p**2*dlog(n1*n2)

	END FUNCTION Fun_logs

!*******************************************************************!

!*******************************************************************!
!               INTEGRAL CUTOFF FOR THE u QUARK                     !
!*******************************************************************!

    FUNCTION Fcut_u(p)
	Implicit none

    Integer :: poly

    Real*8 :: M_u,xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Fcut_u,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,xmue,xmuu,xnuel,xnumu

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmuu=xmu-2.d0*(xmue-xnuel)/3.d0


	M_u=mou-2.d0*G_s*alpha-2.d0*K_s*beta*gamma
	Ep=dsqrt(p**2+M_u**2)


    Fcut_u=p**2*Ep

	END FUNCTION Fcut_U
!*******************************************************************!

!*******************************************************************!
!               INTEGRAL CUTOFF FOR THE d QUARK                     !
!*******************************************************************!

    FUNCTION Fcut_d(p)
	Implicit none

    Integer :: poly

    Real*8 :: M_d,xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Fcut_d,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,xmue,xmud,xnuel,xnumu

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmud=xmu+(xmue-xnuel)/3.d0


	M_d=mod-2.d0*G_s*beta-2.d0*K_s*alpha*gamma
	Ep=dsqrt(p**2+M_d**2)


    Fcut_d=p**2*Ep

	END FUNCTION Fcut_d
!*******************************************************************!

!*******************************************************************!
!               INTEGRAL CUTOFF FOR THE S QUARK                     !
!*******************************************************************!

    FUNCTION Fcut_s(p)
	Implicit none

    Integer :: poly

    Real*8 :: M_s,xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Fcut_s,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,xmue,xmud,xnuel,xnumu

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmud=xmu+(xmue-xnuel)/3.d0

	M_s=mos-2.d0*G_s*gamma-2.d0*K_s*alpha*beta
	Ep=dsqrt(p**2+M_s**2)


    Fcut_s=p**2*Ep

	END FUNCTION Fcut_s

!*******************************************************************!

!*******************************************************************!
!             TOTAL GRANDPOTENTIAL FOR ALL QUARKS                   !
!*******************************************************************!

    SUBROUTINE Omega_T(x,O_qu,O_qd,O_qs,Omtotal)
	Implicit none

    Integer :: poly,IRULE,INTERV

    Real*8 :: al,be,ga,alpha,beta,gamma,pi2,N_c,xT,x(7)
    Real*8 :: Lambda0,G_s,K_s,xmu,mou,mod,mos,Omtotal
    Real*8 :: A,B,ERRABS,ERREL,ERREST,BOUND,ERRABI,ERRELI
    Real*8 :: Fun_logu,Fun_logd,Fun_logs,xlogu,xlogd,xlogs
    Real*8 ::Fcut_u,Fcut_d,Fcut_s,xcutu,xcutd,xcuts,Valor
    Real*8 :: xmue,O_qu,O_qd,O_qs,inter

	EXTERNAL DQDAG,DQDAGI
	EXTERNAL Fun_logu,Fun_logd,Fun_logs,Fcut_u,Fcut_d,Fcut_s

	Common/par_INTEG/INTERV,IRULE
    Common/par2_Int/A,B,ERRABS,ERREL,BOUND,ERRABI,ERRELI,ERREST
	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos

	Valor=-N_c/pi2

	CALL DQDAGI(Fun_logu,BOUND,INTERV,ERRABI,ERRELI,xlogu,ERREST)
	CALL DQDAGI(Fun_logd,BOUND,INTERV,ERRABI,ERRELI,xlogd,ERREST)
	CALL DQDAGI(Fun_logs,BOUND,INTERV,ERRABI,ERRELI,xlogs,ERREST)

	CALL DQDAG(Fcut_u,A,B,ERRABS,ERREL,IRULE,xcutu,ERREST)
	CALL DQDAG(Fcut_d,A,B,ERRABS,ERREL,IRULE,xcutd,ERREST)
	CALL DQDAG(Fcut_s,A,B,ERRABS,ERREL,IRULE,xcuts,ERREST)

    inter=4.d0*K_s*x(1)*x(2)*x(3)/3.d0


	O_qu=Valor*(xcutu+xT*xlogu)+G_s*x(1)**2+inter
	O_qd=Valor*(xcutd+xT*xlogd)+G_s*x(2)**2+inter
	O_qs=Valor*(xcuts+xT*xlogs)+G_s*x(3)**2+inter

	Omtotal=O_qu+O_qd+O_qs

    !O_int=Valor*(xcutu+xcutd+xcuts)
	!O_fg=Valor*xT*(xlogu+xlogd+xlogs)
	!coupl=G_s*(alpha**2+beta**2+gamma**2)+4.d0*K_s*alpha*beta*gamma


    !Omega_T=O_int+O_fg+coupl

	END SUBROUTINE Omega_T
!*******************************************************************!

!*******************************************************************!
!             DERIVATIVE OF OMEGA_u RESPECT TO ALPHA                !
!*******************************************************************!

    FUNCTION Da_Omulog(p)
	Implicit none

    Integer :: poly

    Real*8 :: M_u,xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Da_Omulog,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,xmue ,xmuu,xnuel,xnumu

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmuu=xmu-2.d0*(xmue-xnuel)/3.d0

	M_u=mou-2.d0*G_s*alpha-2.d0*K_s*beta*gamma
	Ep=dsqrt(p**2+M_u**2)
	par1=(Ep-xmuu)
	par2=(Ep+xmuu)
	e1=dexp(-par1/xT)
	e2=dexp(-par2/xT)


	n1=1.d0+e1
	n2=1.d0+e2


    Da_Omulog=2.d0*G_s*M_u/(xT*Ep)*p**2*(e1/n1+e2/n2)

	END FUNCTION Da_Omulog

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!QUARK U______________________alpha-alpha-alpha

    FUNCTION Da_Omuint(p)
	Implicit none

    Integer :: poly

    Real*8 :: M_u,xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Da_Omuint,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,dM_uda,xmue,xmud,xnuel,xnumu

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmud=xmu+(xmue-xnuel)/3.d0


	M_u=mou-2.d0*G_s*alpha-2.d0*K_s*beta*gamma
	dM_uda=-2.d0*G_s
	Ep=dsqrt(p**2+M_u**2)

    Da_Omuint=p**2*M_u*dM_uda/Ep


    END FUNCTION Da_Omuint

!END OF DERIVATIVE OF OMEGA u RESPECT TO ALPHA
!*******************************************************************!

!*******************************************************************!
!             DERIVATIVE OF OMEGA_u RESPECT TO BETA                 !
!*******************************************************************!


    FUNCTION Dbe_Omulog(p)
	Implicit none

    Integer :: poly

    Real*8 :: M_u,xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Dbe_Omulog,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,xmue,xmuu ,xnuel,xnumu

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmuu=xmu-2.d0*(xmue-xnuel)/3.d0


	M_u=mou-2.d0*G_s*alpha-2.d0*K_s*beta*gamma
	Ep=dsqrt(p**2+M_u**2)
	par1=(Ep-xmuu)
	par2=(Ep+xmuu)
	e1=dexp(-par1/xT)
	e2=dexp(-par2/xT)


	n1=1.d0+e1
	n2=1.d0+e2


    Dbe_Omulog=2.d0*K_s*gamma*M_u/(xT*Ep)*p**2*(e1/n1+e2/n2)

	END FUNCTION Dbe_Omulog

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!QUARK U______________________beta-beta-beta

    Function Dbe_Omuint(p)
	Implicit none

    Integer :: poly

    Real*8 :: M_u,xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Dbe_Omuint,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,dM_udbe,xmue,xmuu,xnuel,xnumu

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmuu=xmu-2.d0*(xmue-xnuel)/3.d0


	M_u=mou-2.d0*G_s*alpha-2.d0*K_s*beta*gamma
	dM_udbe=-2.d0*K_s*gamma
	Ep=dsqrt(p**2+M_u**2)

    Dbe_Omuint=p**2*M_u*dM_udbe/Ep


    END FUNCTION Dbe_Omuint

!END OF DERIVATIVE OF OMEGA u RESPECT TO BETA
!*******************************************************************!

!*******************************************************************!
!             DERIVATIVE OF OMEGA_u RESPECT TO GAMMA                !
!*******************************************************************!

    FUNCTION Dga_Omulog(p)
	Implicit none

    Integer :: poly

    Real*8 :: M_u,xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Dga_Omulog,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,xmue,xmuu,xnuel,xnumu

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmuu=xmu-2.d0*(xmue-xnuel)/3.d0



	M_u=mou-2.d0*G_s*alpha-2.d0*K_s*beta*gamma
	Ep=dsqrt(p**2+M_u**2)
	par1=(Ep-xmuu)
	par2=(Ep+xmuu)
	e1=dexp(-par1/xT)
	e2=dexp(-par2/xT)


	n1=1.d0+e1
	n2=1.d0+e2


    Dga_Omulog=2.d0*K_s*beta*M_u/(xT*Ep)*p**2*(e1/n1+e2/n2)

	END FUNCTION Dga_Omulog

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!QUARK U______________________gamma-gamma-gamma

    FUNCTION Dga_Omuint(p)
	Implicit none

    Integer :: poly

    Real*8 :: M_u,xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Dga_Omuint,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,dM_udga,xmue,xmuu,xnuel,xnumu

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmuu=xmu-2.d0*(xmue-xnuel)/3.d0


	M_u=mou-2.d0*G_s*alpha-2.d0*K_s*beta*gamma
	dM_udga=-2.d0*K_s*beta
	Ep=dsqrt(p**2+M_u**2)

    Dga_Omuint=p**2*M_u*dM_udga/Ep


    END FUNCTION Dga_Omuint

!END OF DERIVATIVE OF OMEGA u RESPECT TO GAMMA
!*******************************************************************!

!*******************************************************************!
!             DERIVATIVE OF OMEGA_d RESPECT TO APLHA                !
!*******************************************************************!

    FUNCTION Da_Omdlog(p)
	Implicit none

    Integer :: poly

    Real*8 :: xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Da_Omdlog,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,M_d,xmue,xmud,xnuel,xnumu

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmud=xmu+(xmue-xnuel)/3.d0



	M_d=mod-2.d0*G_s*beta-2.d0*K_s*alpha*gamma
	Ep=dsqrt(p**2+M_d**2)
	par1=(Ep-xmud)
	par2=(Ep+xmud)
	e1=dexp(-par1/xT)
	e2=dexp(-par2/xT)


	n1=1.d0+e1
	n2=1.d0+e2


    Da_Omdlog=2.d0*K_s*gamma*M_d/(xT*Ep)*p**2*(e1/n1+e2/n2)

	END FUNCTION Da_Omdlog
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!QUARK D___________________alpha-alpha-alpha

    FUNCTION Da_Omdint(p)
	Implicit none

    Integer :: poly

    Real*8 :: xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Da_Omdint,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,dM_dda,M_d,xmue,xmud,xnuel,xnumu

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmud=xmu+(xmue-xnuel)/3.d0


	M_d=mod-2.d0*G_s*beta-2.d0*K_s*alpha*gamma
	dM_dda=-2.d0*K_s*gamma
	Ep=dsqrt(p**2+M_d**2)

    Da_Omdint=p**2*M_d*dM_dda/Ep


    END FUNCTION Da_Omdint

!END OF DERIVATIVE OF OMEGA d RESPECT TO ALPHA
!*******************************************************************!

!*******************************************************************!
!             DERIVATIVE OF OMEGA_d RESPECT TO BETA                 !
!*******************************************************************!

    FUNCTION Dbe_Omdlog(p)
	Implicit none

    Integer :: poly

    Real*8 :: xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Dbe_Omdlog,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,M_d,xmue,xmud,xnuel,xnumu

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmud=xmu+(xmue-xnuel)/3.d0


	M_d=mod-2.d0*G_s*beta-2.d0*K_s*alpha*gamma
	Ep=dsqrt(p**2+M_d**2)
	par1=(Ep-xmud)
	par2=(Ep+xmud)
	e1=dexp(-par1/xT)
	e2=dexp(-par2/xT)


	n1=1.d0+e1
	n2=1.d0+e2


    Dbe_Omdlog=2.d0*G_s*M_d/(xT*Ep)*p**2*(e1/n1+e2/n2)

	END FUNCTION Dbe_Omdlog

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!QUARK D___________________beta-beta-beta

    FUNCTION Dbe_Omdint(p)
	Implicit none

    Integer :: poly

    Real*8 :: xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Dbe_Omdint,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,dM_ddbe,M_d,xmue,xmud,xnuel,xnumu

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmud=xmu+(xmue-xnuel)/3.d0


	M_d=mod-2.d0*G_s*beta-2.d0*K_s*alpha*gamma
	dM_ddbe=-2.d0*G_s
	Ep=dsqrt(p**2+M_d**2)

    Dbe_Omdint=p**2*M_d*dM_ddbe/Ep


    END FUNCTION Dbe_Omdint

!END OF DERIVATIVE OF OMEGA d RESPECT TO BETA
!*******************************************************************!

!*******************************************************************!
!             DERIVATIVE OF OMEGA_d RESPECT TO GAMMA                !
!*******************************************************************!


    FUNCTION Dga_Omdlog(p)
	Implicit none

Integer :: poly

    Real*8 :: xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Dga_Omdlog,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,M_d,xmue,xmud,xnuel,xnumu

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmud=xmu+(xmue-xnuel)/3.d0

	M_d=mod-2.d0*G_s*beta-2.d0*K_s*alpha*gamma
	Ep=dsqrt(p**2+M_d**2)
	par1=(Ep-xmud)
	par2=(Ep+xmud)
	e1=dexp(-par1/xT)
	e2=dexp(-par2/xT)


	n1=1.d0+e1
	n2=1.d0+e2


    Dga_Omdlog=2.d0*K_s*alpha*M_d/(xT*Ep)*p**2*(e1/n1+e2/n2)

	END FUNCTION Dga_Omdlog

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!QUARK D___________________gamma-gamma-gamma

    FUNCTION Dga_Omdint(p)
	Implicit none

    Integer :: poly

    Real*8 :: xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Dga_Omdint,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,dM_ddga,M_d,xmue,xmud,xnuel,xnumu

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmud=xmu+(xmue-xnuel)/3.d0


	M_d=mod-2.d0*G_s*beta-2.d0*K_s*alpha*gamma
	dM_ddga=-2.d0*K_s*alpha
	Ep=dsqrt(p**2+M_d**2)

    Dga_Omdint=p**2*M_d*dM_ddga/Ep



    END FUNCTION Dga_Omdint

!END OF DERIVATIVE OF OMEGA d RESPECT TO GAMMA
!*******************************************************************!

!*******************************************************************!
!             DERIVATIVE OF OMEGA_s RESPECT TO ALPHA                !
!*******************************************************************!


    FUNCTION Da_Omslog(p)
	Implicit none

    Integer :: poly

    Real*8 :: xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Da_Omslog,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,M_s,xmue,xmud,xnuel,xnumu

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmud=xmu+(xmue-xnuel)/3.d0

	M_s=mos-2.d0*G_s*gamma-2.d0*K_s*alpha*beta
	Ep=dsqrt(p**2+M_s**2)
	par1=(Ep-xmud)
	par2=(Ep+xmud)
	e1=dexp(-par1/xT)
	e2=dexp(-par2/xT)

	n1=1.d0+e1
	n2=1.d0+e2


    Da_Omslog=2.d0*K_s*beta*M_s/(xT*Ep)*p**2*(e1/n1+e2/n2)


	END FUNCTION Da_Omslog

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!QUARK S___________________alpha-alpha-alpha

    FUNCTION Da_Omsint(p)
	Implicit none

    Integer :: poly

    Real*8 :: xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Da_Omsint,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,dM_sda,M_s,xmue,xmud,xnumu,xnuel

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmud=xmu+(xmue-xnuel)/3.d0


	M_s=mos-2.d0*G_s*gamma-2.d0*K_s*alpha*beta
	dM_sda=-2.d0*K_s*beta
	Ep=dsqrt(p**2+M_s**2)

    Da_Omsint=p**2*M_s*dM_sda/Ep


    END FUNCTION Da_Omsint

!END OF DERIVATIVE OF OMEGA s RESPECT TO ALPHA
!*******************************************************************!

!*******************************************************************!
!             DERIVATIVE OF OMEGA_s RESPECT TO BETA                 !
!*******************************************************************!

    FUNCTION Dbe_Omslog(p)
    Implicit none

    Integer :: poly

    Real*8 :: xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Dbe_Omslog,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,M_s,xmue,xmud,xnuel,xnumu

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmud=xmu+(xmue-xnuel)/3.d0

	M_s=mos-2.d0*G_s*gamma-2.d0*K_s*alpha*beta
	Ep=dsqrt(p**2+M_s**2)
	par1=(Ep-xmud)
	par2=(Ep+xmud)
	e1=dexp(-par1/xT)
	e2=dexp(-par2/xT)

	n1=1.d0+e1
	n2=1.d0+e2


    Dbe_Omslog=2.d0*K_s*alpha*M_s/(xT*Ep)*p**2*(e1/n1+e2/n2)

	END FUNCTION Dbe_Omslog

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!QUARK S___________________beta-beta-beta

    FUNCTION Dbe_Omsint(p)
	Implicit none

    Integer :: poly

    Real*8 :: xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Dbe_Omsint,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,dM_sdbe,M_s,xmue,xmud,xnuel,xnumu

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmud=xmu+(xmue-xnuel)/3.d0


	M_s=mos-2.d0*G_s*gamma-2.d0*K_s*alpha*beta
	dM_sdbe=-2.d0*K_s*alpha
	Ep=dsqrt(p**2+M_s**2)

    Dbe_Omsint=p**2*M_s*dM_sdbe/Ep

    END FUNCTION Dbe_Omsint

!END OF DERIVATIVE OF OMEGA s RESPECT TO BETA
!*******************************************************************!

!*******************************************************************!
!             DERIVATIVE OF OMEGA_s RESPECT TO GAMMA                !
!*******************************************************************!

    FUNCTION Dga_Omslog(p)
	Implicit none

    Integer :: poly

    Real*8 :: xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Dga_Omslog,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,M_s,xmue,xmud,xnuel,xnumu

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmud=xmu+(xmue-xnuel)/3.d0


	M_s=mos-2.d0*G_s*gamma-2.d0*K_s*alpha*beta
	Ep=dsqrt(p**2+M_s**2)
	par1=(Ep-xmud)
	par2=(Ep+xmud)
	e1=dexp(-par1/xT)
	e2=dexp(-par2/xT)

	n1=1.d0+e1
	n2=1.d0+e2


    Dga_Omslog=2.d0*G_s*M_s/(xT*Ep)*p**2*(e1/n1+e2/n2)

	END FUNCTION Dga_Omslog

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!QUARK S___________________gamma-gamma-gamma

    FUNCTION Dga_Omsint(p)
	Implicit none

    Integer :: poly

    Real*8 :: xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Dga_Omsint,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,dM_sdga,M_s,xmue,xmud,xnuel,xnumu

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmud=xmu+(xmue-xnuel)/3.d0


	M_s=mos-2.d0*G_s*gamma-2.d0*K_s*alpha*beta
	dM_sdga=-2.d0*G_s
	Ep=dsqrt(p**2+M_s**2)

    Dga_Omsint=p**2*M_s*dM_sdga/Ep

    END FUNCTION Dga_Omsint

!END OF DERIVATIVE OF OMEGA s RESPECT TO GAMMA
!*******************************************************************************!

!*******************************************************************************!
!      GRAND POTENTIAL DERIVATIVE W/ RESPECT TO ALPHA, BETA AND GAMMA           !
!*******************************************************************************!

    SUBROUTINE Deriv_Omega(x,Deralf,Derbet,Dergam)
	Implicit none

    Integer :: poly,INTERV,IRULE

    Real*8 :: x(7),pi2,N_c,xT,Lambda0,G_s,K_s,xmu
    Real*8 :: A,B,ERRABS,ERREL,BOUND,ERRABI,ERRELI,ERREST
    Real*8 :: Da_Omulog,Da_Omdlog,Da_Omslog,Da_Omuint,Da_Omdint
    Real*8 :: Da_Omsint,alulog,aldlog,alslog,alintu,alintd,alints
    Real*8 :: Deralf,Derbet,Dergam,valor,teral,terbe,terga
    Real*8 :: Dbe_Omulog,Dbe_Omuint,Dga_Omulog,Dga_Omuint,Dbe_Omdlog
    Real*8 :: Dbe_Omdint,Dga_Omdlog,Dga_Omdint,Dbe_Omslog
    Real*8 :: Dbe_Omsint,Dga_Omslog,Dga_Omsint
    Real*8 :: beulog,bedlog,beslog,gaulog,gadlog,gaslog
    Real*8 :: beintu,beintd,beints,gaintu,gaintd,gaints
    Real*8 :: u_alpha,d_alpha,s_alpha,u_beta,d_beta,s_beta
    Real*8 :: u_gamma,d_gamma,s_gamma,xmue


	EXTERNAL DQDAG,DQDAGI,Da_Omulog,Da_Omdlog,Da_Omslog
	EXTERNAL Da_Omuint,Da_Omdint,Da_Omsint
	EXTERNAL Dbe_Omulog,Dbe_Omuint,Dga_Omulog,Dga_Omuint,Dbe_Omdlog
	EXTERNAL Dbe_Omdint,Dga_Omdlog,Dga_Omdint,Dbe_Omslog
	EXTERNAL Dbe_Omsint,Dga_Omslog,Dga_Omsint

	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/par_INTEG/INTERV,IRULE
    Common/par2_Int/A,B,ERRABS,ERREL,BOUND,ERRABI,ERRELI,ERREST

	CALL DQDAGI(Da_Omulog,BOUND,INTERV,ERRABI,ERRELI,alulog,ERREST)
	CALL DQDAGI(Da_Omdlog,BOUND,INTERV,ERRABI,ERRELI,aldlog,ERREST)
	CALL DQDAGI(Da_Omslog,BOUND,INTERV,ERRABI,ERRELI,alslog,ERREST)

	CALL DQDAGI(Dbe_Omulog,BOUND,INTERV,ERRABI,ERRELI,beulog,ERREST)
	CALL DQDAGI(Dbe_Omdlog,BOUND,INTERV,ERRABI,ERRELI,bedlog,ERREST)
	CALL DQDAGI(Dbe_Omslog,BOUND,INTERV,ERRABI,ERRELI,beslog,ERREST)


	CALL DQDAGI(Dga_Omulog,BOUND,INTERV,ERRABI,ERRELI,gaulog,ERREST)
	CALL DQDAGI(Dga_Omdlog,BOUND,INTERV,ERRABI,ERRELI,gadlog,ERREST)
	CALL DQDAGI(Dga_Omslog,BOUND,INTERV,ERRABI,ERRELI,gaslog,ERREST)


	CALL DQDAG(Da_Omuint,A,B,ERRABS,ERREL,IRULE,alintu,ERREST)          !I think the second to last is the output value (alintu)
	CALL DQDAG(Da_Omdint,A,B,ERRABS,ERREL,IRULE,alintd,ERREST)          !Since it is used to calculate something at the bottom (u_alpha)
	CALL DQDAG(Da_Omsint,A,B,ERRABS,ERREL,IRULE,alints,ERREST)


	CALL DQDAG(Dbe_Omuint,A,B,ERRABS,ERREL,IRULE,beintu,ERREST)
	CALL DQDAG(Dbe_Omdint,A,B,ERRABS,ERREL,IRULE,beintd,ERREST)
	CALL DQDAG(Dbe_Omsint,A,B,ERRABS,ERREL,IRULE,beints,ERREST)

	CALL DQDAG(Dga_Omuint,A,B,ERRABS,ERREL,IRULE,gaintu,ERREST)
	CALL DQDAG(Dga_Omdint,A,B,ERRABS,ERREL,IRULE,gaintd,ERREST)
	CALL DQDAG(Dga_Omsint,A,B,ERRABS,ERREL,IRULE,gaints,ERREST)




	valor=-N_c/pi2
      
	teral=4.d0*K_s*x(2)*x(3)
	terbe=4.d0*K_s*x(1)*x(3)
	terga=4.d0*K_s*x(1)*x(2)

	u_alpha=valor*xT*(alulog+alintu/xT)
	d_alpha=valor*xT*(aldlog+alintd/xT)
	s_alpha=valor*xT*(alslog+alints/xT)

	u_beta=valor*xT*(beulog+beintu/xT)
	d_beta=valor*xT*(bedlog+beintd/xT)
	s_beta=valor*xT*(beslog+beints/xT)

	u_gamma=valor*xT*(gaulog+gaintu/xT)
	d_gamma=valor*xT*(gadlog+gaintd/xT)
	s_gamma=valor*xT*(gaslog+gaints/xT)

      
    Deralf=u_alpha+d_alpha+s_alpha+2.d0*G_s*x(1)+teral

	Derbet=u_beta+d_beta+s_beta+2.d0*G_s*x(2)+terbe

	Dergam=u_gamma+d_gamma+s_gamma+2.d0*G_s*x(3)+terga

	END SUBROUTINE
!***************************************************************************!
!                   CALCULATING ALL QUARK DENSITITES                        !
!***************************************************************************!
!Integral for the u quark

    FUNCTION Densi_u(p)
	Implicit none

    Real*8 :: pi2,N_c,xT,Lambda0,G_s,K_s
    Real*8 :: xmu,mou,mod,mos,alpha,beta,gamma,xmue
    Real*8 :: M_u,p,Ep,par1,par2,e1,e2,n1,n2,Densi_u,xmuu
    Real*8 :: xnuel,xnumu

	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmuu=xmu-2.d0*(xmue-xnuel)/3.d0

	M_u=mou-2.d0*G_s*alpha-2.d0*K_s*beta*gamma
	Ep=dsqrt(p**2+M_u**2)
	par1=(Ep-xmuu)
	par2=(Ep+xmuu)
	e1=dexp(-par1/xT)
	e2=dexp(-par2/xT)


	n1=1.d0+e1
	n2=1.d0+e2

	Densi_u=p**2*(e1/n1-e2/n2)

	END FUNCTION Densi_u
!********************************************************************************!
!Integral for the d quark

    FUNCTION Densi_d(p)
	Implicit none

    Real*8 :: pi2,N_c,xT,Lambda0,G_s,K_s
    Real*8 :: xmu,mou,mod,mos,alpha,beta,gamma,xmue
    Real*8 :: M_d,p,Ep,par1,par2,e1,e2,n1,n2,Densi_d,xmud
    Real*8 :: xnuel,xnumu

	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmud=xmu+(xmue-xnuel)/3.d0


	M_d=mod-2.d0*G_s*beta-2.d0*K_s*alpha*gamma
	Ep=dsqrt(p**2+M_d**2)
	par1=(Ep-xmud)
	par2=(Ep+xmud)
	e1=dexp(-par1/xT)
	e2=dexp(-par2/xT)


	n1=1.d0+e1
	n2=1.d0+e2

	Densi_d=p**2*(e1/n1-e2/n2)

	END FUNCTION Densi_d
!********************************************************************************!
!Integral for the S quark

    FUNCTION Densi_s(p)
	Implicit none
	Real*8 pi2,N_c,xT,Lambda0,G_s,K_s
	Real*8 xmu,mou,mod,mos,alpha,beta,gamma,xmue
	Real*8 M_s,p,Ep,par1,par2,e1,e2,n1,n2,Densi_s,xmud
	Real*8 xnuel,xnumu

	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmud=xmu+(xmue-xnuel)/3.d0


	M_s=mos-2.d0*G_s*gamma-2.d0*K_s*alpha*beta
	Ep=dsqrt(p**2+M_s**2)
	par1=(Ep-xmud)
	par2=(Ep+xmud)
	e1=dexp(-par1/xT)
	e2=dexp(-par2/xT)


	n1=1.d0+e1
	n2=1.d0+e2

	Densi_s=p**2*(e1/n1-e2/n2)


	END FUNCTION Densi_s
 
!********************************************************************************!

    SUBROUTINE Denquark(x,rho_u,rho_d,rho_s)
    Implicit none

    Integer :: INTERV,IRULE

    Real*8 :: x(7), rho_u,rho_d,rho_s
    Real*8 :: pi2,N_c,xT,Densi_u,Densi_d,Densi_s
    Real*8 :: xdenu,xdend,xdens
    Real*8 :: A,B,ERRABS,ERREL,BOUND,ERRABI,ERRELI,ERREST

	EXTERNAL DQDAG,DQDAGI,Densi_u,Densi_d,Densi_s

	Common/dat1/pi2,N_c,xT
	Common/par_INTEG/INTERV,IRULE
    Common/par2_Int/A,B,ERRABS,ERREL,BOUND,ERRABI,ERRELI,ERREST

	CALL DQDAGI(Densi_u,BOUND,INTERV,ERRABI,ERRELI,xdenu,ERREST)
	CALL DQDAGI(Densi_d,BOUND,INTERV,ERRABI,ERRELI,xdend,ERREST)
	CALL DQDAGI(Densi_s,BOUND,INTERV,ERRABI,ERRELI,xdens,ERREST)

    rho_u=N_c*xdenu/pi2
	rho_d=N_c*xdend/pi2
    rho_s=N_c*xdens/pi2


	END SUBROUTINE Denquark

!**********************************************************************************!

!***************************************************************************!
!                   ELECTRON CONTRIBUTION                                   !
!***************************************************************************!
!ELECTRON GRAN POTENCIAL INTEGRAL

    FUNCTION F_Oe(y)
	Implicit none

    Real*8 :: me,xmue,e1,e2,n1,n2,Ey,xnuel,xnumu
    Real*8 :: F_Oe,pi2,N_c,xT,xmu,y,alpha,beta,gamma

	Common/masaelectron/me
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu
	Common/dat1/pi2,N_c,xT

  
    Ey=dsqrt(y**2-me**2)


    e1=dexp(-(y+xmue)/xT)
	e2=dexp(-(y-xmue)/xT)

    n1=1.d0+e1
    n2=1.d0+e2
	
	
    F_Oe= Ey**3*(e1/n1+e2/n2)

    END FUNCTION F_Oe

!******************************************************************************
!ELECTRON DENSITY INTEGRAL

    FUNCTION Den_e(y)
	Implicit none

    Real*8 :: me,xmue,e1,e2,n1,n2,Ey,y,xnuel,xnumu
    Real*8 :: Den_e,pi2,N_c,xT,xmu,alpha,beta,gamma

	Common/masaelectron/me
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu
	Common/dat1/pi2,N_c,xT

	  
    Ey=dsqrt(y**2-me**2)


    e1=dexp(-(y+xmue)/xT)
	e2=dexp(-(y-xmue)/xT)

    n1=1.d0+e1
    n2=1.d0+e2
	
		
    Den_e=Ey**3*(e2/n2-e1/n1-e2**2/n2**2+e1**2/n1**2)/xT

    END FUNCTION Den_e

!********************************************************************************

    SUBROUTINE Electron(x,Dene,Ome)
	Implicit none

    Integer :: INTERV,IRULE
    Real*8 :: x(7),A,B
    Real*8 :: A1,B1,ERRABS,ERREL,BOUND,ERRABI,ERRELI,ERREST
    Real*8 :: Dene,Ome,pi2,N_c,xT,xome,xdene,F_Oe,Den_e,me


	External DQDAG,F_Oe,Den_e

	Common/masaelectron/me
	Common/dat1/pi2,N_c,xT
	Common/par_INTEG/INTERV,IRULE
    Common/par2_Int/A,B,ERRABS,ERREL,BOUND,ERRABI,ERRELI,ERREST

	A1=me
	B1=6000.d0


	CALL DQDAG(F_Oe,A1,B1,ERRABS,ERREL,IRULE,xome,ERREST)
    CALL DQDAG(Den_e,A1,B1,ERRABS,ERREL,IRULE,xdene,ERREST)

	Ome=-xome/(3.d0*pi2)
	Dene=xdene/(3.d0*pi2)

	END SUBROUTINE Electron

!!**************************************************************************!

!***************************************************************************!
!                       MUON CONTRIBUTION                                   !
!***************************************************************************!
!MUON GRAN POTENCIAL INTEGRAL

    FUNCTION F_Omuon(y)
	Implicit none

    Real*8 :: xmu,xmue,pi2,N_c,xT,mamu
    Real*8 :: Ey,y,e1,e2,n1,n2,F_Omuon
    Real*8 :: alpha,gamma,beta,xnuel,xnumu


	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu
	Common/dat1/pi2,N_c,xT
	Common/muon/mamu

	IF (xmue-xnuel > mamu) THEN

        Ey=dsqrt(y**2-mamu**2)
        e1=dexp(-(y+(xmue-xnuel+xnumu))/xT)
        e2=dexp(-(y-(xmue-xnuel+xnumu))/xT)

	ELSE
        Ey=0.0d0
	ENDIF

    n1=1.d0+e1
    n2=1.d0+e2
	
    F_Omuon=Ey**3*(e1/n1+e2/n2)

    END FUNCTION F_Omuon

!********************************************************************************
!MUON DENSITY INTEGRAL

    FUNCTION Den_muon(y)
	Implicit none

    Real*8 :: Ey,y,e1,e2,n1,n2,Den_muon
    Real*8 :: xmu,xmue,pi2,N_c,xT,mamu
    Real*8 :: alpha,beta,gamma,xnuel,xnumu


	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu
	Common/dat1/pi2,N_c,xT
	Common/muon/mamu


	IF(xmue-xnuel > mamu) THEN

        Ey=dsqrt(y**2-mamu**2)
        e1=dexp(-(y+(xmue-xnuel+xnumu))/xT)
        e2=dexp(-(y-(xmue-xnuel+xnumu))/xT)

	ELSE
        Ey=0.0d0
	ENDIF

    n1=1.d0+e1
    n2=1.d0+e2
	
		
    Den_muon= Ey**3*(e2/n2-e1/n1-e2**2/n2**2+e1**2/n1**2)/xT

    END FUNCTION Den_muon

!********************************************************************************

    SUBROUTINE MUONS(x,Dmuon,Omuon)
	Implicit none

    Integer :: INTERV,IRULE
    Real*8 :: x(7)
    Real*8 :: A1,B1,ERRABS,ERREL,BOUND,ERRABI,ERRELI,ERREST
    Real*8 :: pi2,N_c,xT,mamu,Omuon,Dmuon
    Real*8 :: xomuon,xdenmuon,F_Omuon,Den_muon,A,B
 
	External DQDAG,F_Omuon,Den_muon

	Common/dat1/pi2,N_c,xT
	Common/par_INTEG/INTERV,IRULE
    Common/par2_Int/A,B,ERRABS,ERREL,BOUND,ERRABI,ERRELI,ERREST
	Common/muon/mamu

    A1=mamu
	B1=6000.d0

	CALL DQDAG(F_Omuon,A1,B1,ERRABS,ERREL,IRULE,xomuon,ERREST)
    CALL DQDAG(Den_muon,A1,B1,ERRABS,ERREL,IRULE,xdenmuon,ERREST)

    IF(x(5) > mamu) THEN
        Omuon=-xomuon/(3.d0*pi2)
	ELSE
        Omuon=0.0d0
	ENDIF

	Dmuon=xdenmuon/(3.d0*pi2)

	END SUBROUTINE MUONS

!***************************************************************************!
!             CONTRIBUTION OF THE PARTICLES TO THE ENTROPY                  !
!***************************************************************************!
!(Leptonic) Contribution from electrons

    Function Ele_entro(y)
	Implicit none

    Real*8 :: Ele_entro,alpha,beta,gamma,xmu,xmue
    Real*8 :: me,y,Ey,pi2,N_c,xT,ter1,ter2
    Real*8 :: e1,e2,n1,n2,xnuel,xnumu

	Common/dat1/pi2,N_c,xT
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu
	Common/masaelectron/me


    Ey=dsqrt(y**2-me**2)


    e1=dexp(-(y+xmue)/xT)
	e2=dexp(-(y-xmue)/xT)

    n1=1.d0+e1
    n2=1.d0+e2
	
    ter1=(e1/n1-e1**2/n1**2)
	ter2=(e2**2/n2**2-e2/n2)
		
    Ele_entro=Ey**3*(ter1*(y+xmue)+ter2*(y-xmue))/xT**2

	END FUNCTION Ele_entro

!*******************************************************************************
!(Leptonic) Contribution from muons

    FUNCTION Muon_entro(y)
	Implicit none

    Real*8 :: Muon_entro,alpha,beta,gamma,xmu,xmue
    Real*8 :: mamu,y,Ey,pi2,N_c,xT,ter1,ter2
    Real*8 :: e1,e2,n1,n2,xnuel,xnumu

	Common/dat1/pi2,N_c,xT
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu
	Common/muon/mamu

    IF(xmue-xnuel > mamu) THEN
        Ey=dsqrt(y**2-mamu**2)
    ELSE
        Ey=0.0d0
	ENDIF

    e1=dexp(-(y+(xmue-xnuel+xnumu))/xT)
	e2=dexp(-(y-(xmue-xnuel+xnumu))/xT)

    n1=1.d0+e1
    n2=1.d0+e2
	
    ter1=(e1/n1-e1**2/n1**2)*(y+(xmue-xnuel+xnumu))
	ter2=(e2**2/n2**2-e2/n2)*(y-(xmue-xnuel+xnumu))
		
    Muon_entro=Ey**3*(ter1+ter2)/xT**2

	END FUNCTION Muon_entro

!*******************************************************************************   

    SUBROUTINE Lepentropy(x,eletropy,mutropy)
	Implicit none

    Integer :: INTERV,IRULE
    Real*8 :: x(7),eletropy,mutropy
    Real*8 :: A,B,ERRABS,ERREL,BOUND,ERRABI,ERRELI,ERREST
    Real*8 :: Ele_entro,Muon_entro,xelete,xmuote,pi2,N_c,xT
    Real*8 :: BOUND1,BOUND2,me,mamu,A1,A2,B1

	External DQDAG,Ele_entro,Muon_entro

    Common/dat1/pi2,N_c,xT
	Common/par_INTEG/INTERV,IRULE
    Common/par2_Int/A,B,ERRABS,ERREL,BOUND,ERRABI,ERRELI,ERREST
	Common/muon/mamu
	Common/masaelectron/me
      
	A1=me
	A2=mamu
	B1=6000.d0

	CALL DQDAG(Ele_entro,A1,B1,ERRABS,ERREL,IRULE,xelete,ERREST)
    CALL DQDAG(Muon_entro,A2,B1,ERRABS,ERREL,IRULE,xmuote,ERREST)

    eletropy=-xelete/(3.d0*pi2)
    mutropy=-xmuote/(3.d0*pi2)

	END SUBROUTINE Lepentropy

!***************************************************************************!
!                   QUARK CONTRIBUTION TO THE ENTROPY                       !
!***************************************************************************!

    FUNCTION Entro_q_u(p)
	Implicit none

    Integer :: poly
    Real*8 :: M_u,xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Entro_q_u,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,xmue,xmuu,xnuel,xnumu

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmuu=xmu-2.d0*(xmue-xnuel)/3.d0

	M_u=mou-2.d0*G_s*alpha-2.d0*K_s*beta*gamma
	Ep=dsqrt(p**2+M_u**2)
	par1=(Ep-xmuu)
	par2=(Ep+xmuu)
	e1=dexp(-par1/xT)
	e2=dexp(-par2/xT)


	n1=1.d0+e1
	n2=1.d0+e2

    Entro_q_u=p**2*(e1/n1*(Ep-xmuu)+e2/n2*(Ep+xmuu))/xT**2

	END FUNCTION Entro_q_u

!********************************************************************************

    FUNCTION Entro_q_d(p)
	Implicit none

    Integer :: poly
    Real*8 :: M_d,xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Entro_q_d,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,xmue,xmud,xnuel,xnumu

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmud=xmu+(xmue-xnuel)/3.d0

	M_d=mod-2.d0*G_s*beta-2.d0*K_s*alpha*gamma
	Ep=dsqrt(p**2+M_d**2)
	par1=(Ep-xmud)
	par2=(Ep+xmud)
	e1=dexp(-par1/xT)
	e2=dexp(-par2/xT)


	n1=1.d0+e1
	n2=1.d0+e2

    Entro_q_d=p**2*(e1/n1*(Ep-xmud)+e2/n2*(Ep+xmud))/xT**2

	END FUNCTION Entro_q_d

!********************************************************************************

    FUNCTION Entro_q_s(p)
	Implicit none

    Integer :: poly
    Real*8 :: M_s,xmu,xT,par1,par2,e1,e2,n1,n2,mou,mod,mos
    Real*8 :: Entro_q_s,pi2,N_c,Lambda0,G_s,K_s,alpha,beta,gamma
    Real*8 :: Ep,p,xmue,xmud,xnuel,xnumu

	Common/polyakov/poly
	Common/dat1/pi2,N_c,xT
    Common/dat2/Lambda0,G_s,K_s
	Common/masas/mou,mod,mos
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

	xmud=xmu+(xmue-xnuel)/3.d0
	M_s=mos-2.d0*G_s*gamma-2.d0*K_s*alpha*beta

	Ep=dsqrt(p**2+M_s**2)
	par1=(Ep-xmud)
	par2=(Ep+xmud)
	e1=dexp(-par1/xT)
	e2=dexp(-par2/xT)


	n1=1.d0+e1
	n2=1.d0+e2

    Entro_q_s=p**2*(e1/n1*(Ep-xmud)+e2/n2*(Ep+xmud))/xT**2

	END FUNCTION Entro_q_s
!*******************************************************************************

    SUBROUTINE Entropyquark(x,Equ,Eqd,Eqs,Eqtot)
    Implicit none

    Integer :: INTERV,IRULE
    Real*8 :: pi2,N_c,xT,Equ,Eqd,Eqs,Eqtot
    Real*8 :: A,B,ERRABS,ERREL,BOUND,ERRABI,ERRELI,ERREST
    Real*8 :: Entro_q_u,Entro_q_d,Entro_q_s,xentrou,xentrod,xentros
    Real*8 :: Fun_logu,Fun_logd,Fun_logs,Fcut_u,Fcut_d,Fcut_s
    Real*8 :: xlogu,xlogd,xlogs,xcutu,xcutd,xcuts,x(7)
    Real*8 :: O_qu,O_qd,O_qs,Lambda0,G_s,K_s,Valor

    EXTERNAL DQDAGI,DQDAG,Entro_q_u,Entro_q_d,Entro_q_s
	EXTERNAL Fun_logu,Fun_logd,Fun_logs,Fcut_u,Fcut_d,Fcut_s

	Common/dat1/pi2,N_c,xT
	Common/par_INTEG/INTERV,IRULE
    Common/par2_Int/A,B,ERRABS,ERREL,BOUND,ERRABI,ERRELI,ERREST
    Common/dat2/Lambda0,G_s,K_s

	CALL DQDAGI(Entro_q_u,BOUND,INTERV,ERRABI,ERRELI,xentrou,ERREST)
	CALL DQDAGI(Entro_q_d,BOUND,INTERV,ERRABI,ERRELI,xentrod,ERREST)
	CALL DQDAGI(Entro_q_s,BOUND,INTERV,ERRABI,ERRELI,xentros,ERREST)

	Valor=-N_c/pi2

	CALL DQDAGI(Fun_logu,BOUND,INTERV,ERRABI,ERRELI,xlogu,ERREST)
	CALL DQDAGI(Fun_logd,BOUND,INTERV,ERRABI,ERRELI,xlogd,ERREST)
	CALL DQDAGI(Fun_logs,BOUND,INTERV,ERRABI,ERRELI,xlogs,ERREST)

	CALL DQDAG(Fcut_u,A,B,ERRABS,ERREL,IRULE,xcutu,ERREST)
	CALL DQDAG(Fcut_d,A,B,ERRABS,ERREL,IRULE,xcutd,ERREST)
	CALL DQDAG(Fcut_s,A,B,ERRABS,ERREL,IRULE,xcuts,ERREST)


    O_qu=Valor*(xcutu+xT*xlogu)!+G_s*x(1)**2+4.d0*K_s*x(1)*x(2)*x(3)/3.d0
	O_qd=Valor*(xcutd+xT*xlogd)!+G_s*x(2)**2+4.d0*K_s*x(1)*x(2)*x(3)/3.d0
	O_qs=Valor*(xcuts+xT*xlogs)!+G_s*x(3)**2+4.d0*K_s*x(1)*x(2)*x(3)/3.d0

	Equ=-(O_qu/xT+Valor*xT*xentrou)
	Eqd=-(O_qd/xT+Valor*xT*xentrod)
	Eqs=-(O_qs/xT+Valor*xT*xentros)

	Eqtot=Equ+Eqd+Eqs

	END SUBROUTINE Entropyquark

!***********************************************************************************

    SUBROUTINE energydensity(x, epsillon)
    Implicit none

    Real*8 :: alpha,beta,gamma,x(7), epsillon
    Real*8 :: xmuu, xmud, xmus, hbarc
    Real*8 :: pi2,N_c,xT,xnumu,xmue,xnuel,xmu
    Real*8 :: rho_u, rho_d, rho_s
    Real*8 :: Equ, Eqd, Eqs, Eqtot
    Real*8 :: O_qu,O_qd,O_qs, Omtotal

    Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu
    Common/dat1/pi2,N_c,xT

    call Denquark(x,rho_u,rho_d,rho_s)
    call Entropyquark(x,Equ,Eqd,Eqs,Eqtot)
    call Omega_T(x,O_qu,O_qd,O_qs,Omtotal)

    xmuu=xmu-2.d0*(xmue-xnuel)/3.d0
    xmud=xmu+(xmue-xnuel)/3.d0
    xmus=xmud

    epsillon = Omtotal + xT*(Eqtot) + (xmuu*rho_u + xmud*rho_d + xmus*rho_s) !Mev^4

    END SUBROUTINE energydensity

!***********************************************************************************

    SUBROUTINE pressure_sub(x, pressure)
    Implicit none

    Real*8 :: alpha,beta,gamma,x(7)
    Real*8 :: pi2,N_c,xT, xnumu,xmue,xnuel,xmu
    Real*8 :: Equ, Eqd, Eqs, Eqtot, epsillon, pressure

    Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu
    Common/dat1/pi2,N_c,xT

    call Entropyquark(x,Equ,Eqd,Eqs,Eqtot)
    call energydensity(x, epsillon)

    pressure = xT*Eqtot - epsillon !Mev^4

    END SUBROUTINE pressure_sub

!***************************************************************************!
!              NEUTRINOS CONTRIBUTION TO THE GRAN POTENTIAL                 !
!***************************************************************************!
!INTEGRAL PARA EL GRAN POTENCIAL DE LOS NEUTRINOS DEL ELECTRON

    FUNCTION F_ONUEL(y)
	Implicit none

    Real*8 :: y,xnuel,e1,e2,n1,n2,F_ONUEL
    Real*8 :: pi2,N_c,xT,xnumu,alpha,beta,gamma,xmu,xmue

	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu
	Common/dat1/pi2,N_c,xT

    e1=dexp(-(y+xnuel)/xT)
	e2=dexp(-(y-xnuel)/xT)

    n1=1.d0+e1
    n2=1.d0+e2
	
	
    F_ONUEL= y**3*(e1/n1+e2/n2)

    END FUNCTION F_ONUEL

!****************************************************************************
!ELECTRON NEUTRINO DENSITY INTEGRAL

    FUNCTION Den_nuel(y)
	Implicit none

    Real*8 :: y,xnuel,e1,e2,n1,n2,Den_nuel
    Real*8 :: pi2,N_c,xT,xnumu,alpha,beta,gamma,xmu,xmue

	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu
	Common/dat1/pi2,N_c,xT

    e1=dexp(-(y+xnuel)/xT)
	e2=dexp(-(y-xnuel)/xT)

    n1=1.d0+e1
    n2=1.d0+e2

    Den_nuel= y**3*(e2/n2-e1/n1-e2**2/n2**2+e1**2/n1**2)/xT
	
	!(-e1/n1**2+e2/n2**2)
!

	END FUNCTION Den_nuel

!*******************************************************************************
!MUON NEUTRINO DENSITY INTEGRAL

    FUNCTION F_ONUMU(y)
	Implicit none

    Real*8 :: y,xnumu,e1,e2,n1,n2,F_ONUMU
    Real*8 :: pi2,N_c,xT,xnuel,alpha,beta,gamma,xmu,xmue

	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu
	Common/dat1/pi2,N_c,xT



    e1=dexp(-(y+xnumu)/xT)
	e2=dexp(-(y-xnumu)/xT)

    n1=1.d0+e1
    n2=1.d0+e2
	
	
    F_ONUMU= y**3*(e1/n1+e2/n2)

    END FUNCTION F_ONUMU

!************************************************************************
!MUON NEUTRINOS DENSITY INTEGRAL

    FUNCTION Den_numu(y)
	Implicit none

    Real*8 :: y,xnuel,e1,e2,n1,n2,Den_numu
    Real*8 :: pi2,N_c,xT,xnumu,alpha,beta,gamma,xmu,xmue

	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu
	Common/dat1/pi2,N_c,xT

    e1=dexp(-(y+xnumu)/xT)
	e2=dexp(-(y-xnumu)/xT)

    n1=1.d0+e1
    n2=1.d0+e2

    Den_numu= y**3*(e2/n2-e1/n1-e2**2/n2**2+e1**2/n1**2)/xT

	END FUNCTION Den_numu

!***************************************************************************!
!                 NEUTRINOS CONTRIBUTION TO THE ENTROPY                     !
!***************************************************************************!
!Analytic derivative of the gran potential of the muon neutrino respect to T

    Function ET_numu(y)
	Implicit none

    Real*8 :: ET_numu,alpha,beta,gamma,xmu,xmue
    Real*8 :: xnuel,xnumu,y,pi2,N_c,xT,ter1,ter2
    Real*8 :: e1,e2,n1,n2

	Common/dat1/pi2,N_c,xT
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu

    e1=dexp(-(y+xnumu)/xT)
	e2=dexp(-(y-xnumu)/xT)

    n1=1.d0+e1
    n2=1.d0+e2
	
    ter1=(e1/n1-e1**2/n1**2)
	ter2=(e2**2/n2**2-e2/n2)
		
    ET_numu=y**3*(ter1*(y+xnumu)+ter2*(y-xnumu))/xT**2


	END FUNCTION ET_numu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Analytic derivative of the gran potencial of the electron neutrino respect to T

    FUNCTION ET_nuele(y)
	Implicit none

    Real*8 :: ET_nuele,alpha,beta,gamma,xmu,xmue
    Real*8 :: xnuel,xnumu,y,pi2,N_c,xT,ter1,ter2
    Real*8 :: e1,e2,n1,n2

	Common/dat1/pi2,N_c,xT
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu


    e1=dexp(-(y+xnuel)/xT)
	e2=dexp(-(y-xnuel)/xT)

    n1=1.d0+e1
    n2=1.d0+e2
	
    ter1=(e1/n1-e1**2/n1**2)
	ter2=(e2**2/n2**2-e2/n2)
		
    ET_nuele=y**3*(ter1*(y+xnuel)+ter2*(y-xnuel))/xT**2

	END FUNCTION ET_nuele

!********************************************************************************

    SUBROUTINE NEUTRINOS(x,Onuel,Denuel,Onumu,Dnumu,Entneu)
	Implicit none

    Integer :: INTERV,IRULE
    Real*8 :: x(7),Onuel,Denuel,Onumu,Dnumu,Entneu
    Real*8 :: A,B,ERRABS,ERREL,BOUND,ERRABI,ERRELI,ERREST
    Real*8 :: F_ONUEL,F_ONUMU,Den_nuel,Den_numu,ET_numu,ET_nuele
    Real*8 :: xonuel,xonumu,xdnu_e,xdnu_mu,xentmuon,xentnuel
    Real*8 :: pi2,N_c,xT,B1

	EXTERNAL DQDAG,F_ONUEL,F_ONUMU,Den_nuel,Den_numu
	EXTERNAL ET_numu,ET_nuele

	Common/dat1/pi2,N_c,xT
	Common/par_INTEG/INTERV,IRULE
    Common/par2_Int/A,B,ERRABS,ERREL,BOUND,ERRABI,ERRELI,ERREST

	B1=6000.d0

	CALL DQDAG(F_ONUEL,A,B1,ERRABS,ERREL,IRULE,xonuel,ERREST)
	CALL DQDAG(F_ONUMU,A,B1,ERRABS,ERREL,IRULE,xonumu,ERREST)
    CALL DQDAG(Den_nuel,A,B1,ERRABS,ERREL,IRULE,xdnu_e,ERREST)
    CALL DQDAG(Den_numu,A,B1,ERRABS,ERREL,IRULE,xdnu_mu,ERREST)
    CALL DQDAG(ET_numu,A,B1,ERRABS,ERREL,IRULE,xentmuon,ERREST)
    CALL DQDAG(ET_nuele,A,B1,ERRABS,ERREL,IRULE,xentnuel,ERREST)

	Onuel=-xonuel/(6.d0*pi2)
	Denuel=xdnu_e/(6.d0*pi2)
	Onumu=-xonumu/(6.d0*pi2)
	Dnumu=xdnu_mu/(6.d0*pi2)
	Entneu=-(xentmuon+xentnuel)/(6.d0*pi2)

	END SUBROUTINE NEUTRINOS

!********************************************************************************

    SUBROUTINE FCN(X, F, N)
    Implicit none

    INTEGER :: N,poly, counter
    Real*8 :: X(7),F(N),Deralf,Derbet,Dergam
    Real*8 :: alpha,beta,gamma,xmu,xmue
    Real*8 :: rho_u,rho_d,rho_s,pote,potu
    Real*8 :: Dene,Ome,Dmuon,Omuon,baryon,potuel,pomu
    Real*8 :: Onuel,Denuel,Onumu,Dnumu,Entneu,xnuel,xnumu

	Common/densibaryion/baryon
	Common/polyakov/poly
	Common/fcns/alpha,beta,gamma,xmu,xmue,xnuel,xnumu
    Common/debug/counter

    CALL Mapp(x,pote,potuel,potu,pomu)

	x(4)=dsqrt(potu)
	x(5)=dsqrt(pote)
	x(6)=dsqrt(potuel)
	x(7)=dsqrt(pomu)


	alpha=x(1)
	beta=x(2)
	gamma=x(3)
	xmu=x(4)
	xmue=x(5)
	xnuel=x(6)
	xnumu=x(7)

    CALL Deriv_Omega(x,Deralf,Derbet,Dergam)
	CALL Denquark(x,rho_u,rho_d,rho_s)
	CALL Electron(x,Dene,Ome)
	CALL MUONS(x,Dmuon,Omuon)
	CALL NEUTRINOS(x,Onuel,Denuel,Onumu,Dnumu,Entneu)

!	write(*,*)'before'
!   write(*,*)'alpha=',alpha,'beta=',beta,'gamma=',gamma
!	write(*,*)'xmu=',xmu,'xmue=',xmue


	f(1)=Deralf
	f(2)=Derbet
	f(3)=Dergam
	f(4)=3.d0*baryon-(rho_u+rho_d+rho_s)
	f(5)=2.d0*rho_u-(rho_d+rho_s)-3.d0*(Dene+Dmuon)
	f(6)=0.39d0*baryon-(Dene+Denuel)
	f(7)=Dmuon-Dnumu

    write(*,*)'------------------------------------'
	write(*,*)'after', counter
	write(*,*)'------------------------------------'

	write(*,*)'f1=',f(1),'f2',f(2),'f3',f(3)
	write(*,*)''
	write(*,*)'f4',f(4),'f5',f(5)
	write(*,*)''
	write(*,*)'f6',f(6),'f7',f(7)
	write(*,*)'------------------------------------'
    write(*,*)'alpha=',alpha,(dabs(alpha))**(1.d0/3.d0) !dabs(A) is abs(A) double precise (absolute value funtion)
	write(*,*)'beta=',beta,(dabs(beta))**(1.d0/3.d0)
	write(*,*)'gamma=',gamma,(dabs(gamma))**(1.d0/3.d0)
	write(*,*)'mu_quarks',xmu,'mu_e',xmue
	write(*,*)'mu_neutrinos', xnuel, 'mu_muon',xnumu
    write(*,*)'------------------------------------'

    counter = counter + 1

	END SUBROUTINE FCN

!****************************************************************************

    SUBROUTINE Mapp(x,pote,potuel,potu,pomu)
    Implicit none

    Real*8 :: pote,potu,x(7),potuel,pomu

  
    potu = x(4)**2
    pote = x(5)**2
	potuel=x(6)**2
	pomu=x(7)**2
 
    IF(pote > potu) THEN

        pote = potu
        potu = x(5)**2

    ENDIF

    IF(potuel > pote) THEN

        potuel=pote
        pote=x(6)**2

    ENDIF

	END SUBROUTINE Mapp

!******************************************************************************



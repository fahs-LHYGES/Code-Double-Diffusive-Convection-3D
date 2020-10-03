!ifort -openmp -check all -traceback -g ./src/minpack/MINPACK.f90 ./src/code_FG/data_mod.f90 ./src/code_FG/matrix_mod.f90 ./src/code_FG/maillage_mod.f90 ./src/code_FG/com_mod.f90 ./src/code_FG/init_mod.f90 ./src/code_FG/res_mod.f90 ./src/code_FG/plot_mod.f90 ./src/code_FG/jac_mod.f90 ./src/code_FG/daspk_mod.f90 ./src/code_FG/main.f90 -L./lib/ -lprecond -lsolver -o code_FG
!**************************************************************************************
program analy_solution
  use in_data
  use matrices
  use maillage 
  ! include via modules
  use jacobian_module
  use residual_module
  use daspk_module
  use com_module
  use init_module
  use plot_module
  
  implicit none
  character*15                          :: fichier
  real*8,allocatable, dimension (:,:,:) :: AA,BB,EE,GG
  integer                               :: ind,i,j,sss
  character(len=1024)                   :: nam1,nam01
  INTEGER                               :: NOUT
  REAL*8 ,allocatable, dimension (:)    :: Xinit,dXdtinit
  REAL*8                                :: cpu1,cpu2

  !REAL*8 ,allocatable, dimension (:,:)    ::Jacinit

  pi=3.141592654d0 ! 4.D0*DATAN(1.D0)

  ! Reading input data and calculating coefficients
  write(*,*) 'Input file '
  read(*,*) fichier
  ind=index(fichier,' ')-1
  open(10,file='indata/'//fichier(1:ind))
  open(17,file='indata/'//fichier(1:ind)//'.in')
  open(40,file='indata/'//fichier(1:ind)//'.init')
  open(49,file='indata/'//fichier(1:ind)//'.st')
  call read_data_in_cal_coef()
  print*,' ****   Rayleigh = ', Ra 
  write (nam1,"(F12.5)") Ra
  write (nam01,"(i5)") nva+nvb+nve+nvg

  open(27,file='outdata/'//fichier(1:ind)//'Ra='//trim(nam1)//'Ncoef='//trim(nam01)//'.plt')
  open(19,file='outdata/'//fichier(1:ind)//'Ra='//trim(nam1)//'Ncoef='//trim(nam01)//'.coef')
  open(30,file='outdata/'//fichier(1:ind)//'Ra='//trim(nam1)//'Ncoef='//trim(nam01)//'.num')
  open(33,file='outdata/'//fichier(1:ind)//'Ra='//trim(nam1)//'Ncoef='//trim(nam01)//'.Sher0')
  open(34,file='outdata/'//fichier(1:ind)//'Ra='//trim(nam1)//'Ncoef='//trim(nam01)//'.Sher1')
  open(35,file='outdata/'//fichier(1:ind)//'Ra='//trim(nam1)//'Ncoef='//trim(nam01)//'.Nusselt0')
  open(36,file='outdata/'//fichier(1:ind)//'Ra='//trim(nam1)//'Ncoef='//trim(nam01)//'.Nusselt1')
  open(37,file='outdata/'//fichier(1:ind)//'Ra='//trim(nam1)//'Ncoef='//trim(nam01)//'_time.plt')
  open(38,file='outdata/'//fichier(1:ind)//'Ra='//trim(nam1)//'Ncoef='//trim(nam01)//'_hstep.plt')

  write (30,*) 'Ra',Ra
  write (30,*) 'Trunc Phi_X',ni,nj,nk
  write (30,*) 'Trunc Phi_Y',nl,nm,nn
  write (30,*) 'Trunc Conc ',nu,nv,nw
  write (30,*) 'Trunc Temp ',ns,np,nt

  call get_plot ()

  ! Allocating matrices and vectors and calculating (kroneker) 
  allocate (AA(0:ni,1:nj,1:nk),BB(1:nl,0:nm,1:nn),EE(0:nu,1:nv,0:nw),GG(1:ns,0:np,0:nt))
  allocate (dXdtinit(nve+nvg),Xinit(nve+nvg), matphi(0:ntrucmax,0:ntrucmax))
  
  call cal_matphi ()

  ! Initialization 
  call Init_Guess(EE,GG)
  call matr_to_vect(EE,GG,Xinit)
  !allocate (Jacinit(nve+nvg,nve+nvg))
  !call Cal_derv_res(Xinit,Jacinit,0.2d0)


  MXNEQ=0; LRW=0; LIW=0; NEQ=0; ML=0; MU=0; T=0.0D0; AVDIM=0.0d0
  LENRW=0; LENIW=0;  LENIWP=0;  LENWP=0;  LENPD=0;  MBAND=0; MSAVE=0; LIWP=0; LWP=0
  HU=0; NQU=0; NST=0; NNI=0; NLI=0; NPE=0; NRE=0; NPS=0; NCFN=0; NCFL=0; IDID=0; IOUT=0

  MXNEQ = nve+ nvg        !number of equations = the number of unknowns  
  NEQ = MXNEQ                        !number of equations
  LENWP = NEQ
  
  INFO(:)=0                          !INFO(*) - Use the INFO array to give the code more details about
  !how you want your problem solved.  This array should be
  !dimensioned of length 20, though DDASPK uses only the 
  !first 15 entries.  You must respond to all of the following
  !items, which are arranged as questions.  The simplest use
  !of DDASPK corresponds to setting all entries of INFO to 0.
  !useful values of info... INFO(3)shows results at the actual time step of integration
  !INFO(4) is also about time INFO(7),INFO(8) stepsize... see ddaspk

  INFO(3)   = 0   ! Each t+time step is an output
  INFO(5)   = 1   ! together with subroutine JAC to compute the analytical Jacobian 
  ! 1 analy   0 numer JAC
  INFO(6)   = 0   !  0 full matrix 1 if the matrix is banded
  !INFO(7)  = 1
  INFO(8)   = 1   ! this forces the amplitude of the first time step and (that must be written in RWORK(3) = 1.d-13)
  !INFO(9)  = 1
  INFO(10)  = 0   ! not To enforce nonnegativity in Y during the integration.
  INFO(11)  = 0   ! Initial conditions UPRIME computed by the solver
  ! 1.  Given Y_d, calculate Y_a and Y'_d, or
  ! 2.  Given Y', calculate Y.
  INFO(12)  = 1   ! =1 Kyrlov  method for the solution of the linear system 
  INFO(15)  = 0   ! not Set INFO(15) = 1 to indicate that a JAC routine exists (in this particular case JAC=DBANJA) 
  INFO(16)  = 0  !- option to exclude algebraic variables from the error test.
  !****   Do you wish to control errors locally on
  ! all the variables...
  ! yes - set INFO(16) = 0
  !  no - set INFO(16) = 1
  if(info(12).eq.0) then                   !According to the type of linear solver different dimension of work vectors RWORK and IWORK are required
     LRW = 50 + 11*NEQ +NEQ**2                 !Total length of WORK array.
     LIW = 40 + NEQ +NEQ + NEQ                !Total length of IWORK array.
  else
     ! Krylov
     LRW = 91 + 19*NEQ + LENWP**2             !Total length of WORK array.
     LIW = 40 + 2*LENIWP + NEQ  + NEQ        !Total length of IWORK array.
  end if
  
  allocate (SENPAR(MXNEQ),U(MXNEQ),UPRIME(MXNEQ),RWORK(LRW),IWORK(LIW),ATOL(1),RTOL(1)) 

  ! RWORK(:)=0.0D0;   IWORK(:)=0;
  !UPRIME =0.0D0  
  SENPAR = 0.D0
  U = Xinit

  IPAR=0; RPAR=0.0d0
  RTOL = 1.0D-7   !: Tolerance relative  **these two should be probably removed**
  ATOL = 1.0D-7   !: Tolerance absolue

  if(INFO(8).eq.1) then
     RWORK(3) = 1.0d-8
  end if
  if(INFO(9).eq.1) then
     IWORK(3) = 5
  end if
  INFO(14) = 0 !just initial time step
  ! ****   Do you want to proceed to the integration after
  ! the initial condition calculation is done ...
  !                 yes - set INFO(14) = 0
  !                  no - set INFO(14) = 1

  if (INFO(11).eq.1) then !If the initial UPRIME must be computed numerically
     !do i=1,NM*NN ! is an algebraic variable
     !   IWORK(40+i)=-1
     !end do
     do i=1,nve+nvg ! is a differential variable set
        IWORK(40+i)= 1
     end do
  end if
  ! Preconditionning interer variable
  IWORK(27) = LENWP  ! length of real work space WPrecon for comm betw. PSOL and JAC
  IWORK(28) = LENIWP ! length of integer work space IWPrecond for comm betw. PSOL and JAC
  IPAR(1) = NEQ

  call timert(cpu1)
  sss =  1
  DASP_TIME_STEP_LOOP : DO I=1,Ntime       
     CALL DDASPK (RESH, NEQ, T, U, UPRIME, TOUT(I), INFO, RTOL, ATOL,&
          & IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC1, PSOL,SENPAR,G_RES,K_RES)

     call cal_Phi_C_vx_vy (sss,U)
     call tecplot (sss)
     sss = sss +1 

     !calling of DASPK
!!$     NST = IWORK(11)
!!$     NPE = IWORK(13)
!!$     NRE = IWORK(12) + NPE*MBAND
!!$     LIW = IWORK(17)
!!$     LRW = IWORK(18)
!!$     NNI = IWORK(19)
!!$     NLI = IWORK(20)
!!$     NPS = IWORK(21)
!!$     IF (NNI .NE. 0) AVDIM = REAL(NLI)/REAL(NNI)
!!$     NCFN = IWORK(15)
!!$     NCFL = IWORK(16)
!!$     WRITE (ID_FILE_LOG,'(a)') 'Info *** Final statistics for this run..'
!!$     WRITE (ID_FILE_LOG,'(2(a,i10))')'Info *** RWORK size =',LRW,'   IWORK size =',LIW
!!$     WRITE (ID_FILE_LOG,'(1(a,i10))')'Info *** Number of time steps ................ =',NST
!!$     WRITE (ID_FILE_LOG,'(1(a,i10))')'Info *** Number of residual evaluations ...... =',NRE
!!$     WRITE (ID_FILE_LOG,'(1(a,i10))')'Info *** Number of res. evals. for precond.    =',IPAR(30)
!!$     WRITE (ID_FILE_LOG,'(1(a,i10))')'Info *** Number of preconditioner evaluations  =',NPE
!!$     WRITE (ID_FILE_LOG,'(1(a,i10))')'Info *** Number of preconditioner solves ..... =',NPS
!!$     WRITE (ID_FILE_LOG,'(1(a,i10))')'Info *** Number of nonlinear iterations ...... =',NNI
!!$     WRITE (ID_FILE_LOG,'(1(a,i10))')'Info *** Number of linear iterations ......... =',NLI
!!$     WRITE (ID_FILE_LOG,'(1(a,F8.4))')'Info *** Average Krylov subspace dimension =',AVDIM
!!$     WRITE (ID_FILE_LOG,'(a,2(i6,a))')'Info *** ',NCFN,' nonlinear conv. failures,',NCFL,' linear conv. failures'
!!$     WRITE(ID_FILE_LOG, '(a,I6)') 'Info *** Memory size (Mbytes) : ',INT((LRW*64.+LIW*32.)/1.D6) ;
     HU = RWORK(7)
     NQU = IWORK(8)
     NST = IWORK(11)
     NNI = IWORK(19)
     NLI = IWORK(20)
     WRITE(*,*) T,TOUT(I),HU,NQU,NST,NNI,NLI
     write(37,*) '# Time ',T,' sec'
     write(37,*) '# last time step size  ',HU,' sec'
     do j=1,neq; write(37,'(i4,E20.10)') j,u(j)
     enddo
     IF(IDID == -1) THEN
        !     If..
        !     IDID = -1, the code has taken about 500 steps.  If you want to
        !                  continue, set INFO(1) = 1 and call the code again.
        !                  An additional 500 steps will be allowed.
        INFO(1) = 1
        CYCLE DASP_TIME_STEP_LOOP 
     END IF
  END DO DASP_TIME_STEP_LOOP



  write(*,*)
  NST = IWORK(11)
  NPE = IWORK(13)
  NRE = IWORK(12) !+ NPE*MBAND
  LIW = IWORK(17)
  LRW = IWORK(18)
  NNI = IWORK(19)
  NLI = IWORK(20)
  NPS = IWORK(21)
  IF (NNI .NE. 0) AVDIM = REAL(NLI)/REAL(NNI)
  NCFN = IWORK(15)
  NCFL = IWORK(16)
  WRITE (*,90) LRW,LIW,NST,NRE,NPE,NPS,NNI,NLI,AVDIM,NCFN,NCFL
  !WRITE (23,90) LRW,LIW,NST,NRE,NPE,NPS,NNI,NLI,AVDIM,NCFN,NCFL
90 FORMAT(//' Final statistics for this run..'/&
       '   RWORK size =',I5,'   IWORK size =',I4/&
       '   Number of time steps ................ =',I5/&
       '   Number of residual evaluations ...... =',I5/&
       '   Number of preconditioner evaluations  =',I5/&
       '   Number of preconditioner solves ..... =',I5/&
       '   Number of nonlinear iterations ...... =',I5/&
       '   Number of linear iterations ......... =',I5/&
       '   Average Krylov subspace dimension =',F8.4/&
       I5,' nonlinear conv. failures,',I5,' linear conv. failures')

  call timert(cpu2)

  print*,'cpu time =',cpu2-cpu1

end program analy_solution
!-----------------------------------------------------------------------------
subroutine timert(ttime)
  real*8 ttime
  real temp
  real tarray(2)
  real etime

  !    The first element of the array tarray specifies user time
  temp = etime(tarray) 
  ttime = dble(tarray(1))

  return
end subroutine timert
!-----------------------------------------------------------------------------



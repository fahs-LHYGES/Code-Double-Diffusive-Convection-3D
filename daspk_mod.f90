!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE daspk_module
  integer(kind=4)   :: MXNEQ,LRW,LIW,neq
  integer(kind=4)   :: ML,MU                                    
  integer(kind=4)   :: LENRW,LENIW,LENIWP,LENWP,LENPD,MBAND,MSAVE,LIWP,LWP   
  integer(kind=4)   :: HU,NQU,NST,NNI,NLI,NPE,NRE,NPS,NCFN,NCFL,IDID,IOUT
  REAL(KIND=8)      :: T,AVDIM,fnorm0
  integer(kind=4) ,DIMENSION(30)           :: INFO  
  integer(kind=4),DIMENSION(5)             :: IPAR
  integer(kind=4),ALLOCATABLE,DIMENSION(:) :: IWORK   
  REAL(KIND=8),DIMENSION(150)              :: RPAR
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:)    :: U,UPRIME,RWORK,SENPAR
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:)    :: ATOL,RTOL  
contains
  !  K_RES -- Evaluation of matrix-vector product for Krylov iteration by 
  !           analytic or autoamatic differentiation.
  !           This is the name of a routine you must supply if you have 
  !           chosen INFO(5) > 0 (analytic or ADIFOR evaluation of partial
  !           derivatives in matrix-vector product) and INFO(12) = 1 (Krylov
  !           option).
  !
  !           Depending on the value of INFO(5), K_RES is given by:
  !      **** 1. INFO(5) = 1. The matrix-vector product in the Krylov linear
  !              iteration is computed by a user-input routine, which has the
  !              form
  !
  SUBROUTINE K_RES(T,U,UPRIME,CJ,IRES,RPAR,IPAR,SENPAR,V,AV)
    use jacobian_module
    use in_data,only:nve,nvg
    ! where V is the input vector and AV is the matrix-vector
    ! product AV = (G_y + CJ* G_y') * V.
    IMPLICIT NONE
    REAL(KIND=8):: PD(nve+nvg,nve+nvg)
    REAL(KIND=8):: AV(nve+nvg),V(nve+nvg),U(nve+nvg),UPRIME(nve+nvg),SENPAR(nve+nvg), RPAR(150) 
    INTEGER, DIMENSION (:) ::  IPAR(5)                             
    REAL(KIND=8) :: T,CJ  
    INTEGER :: IRES,NEQ,I,J
    neq = ipar(1)
    !pause 'k_res proc'
    call cal_derv_res(u,pd,cj)    
    do i=1,neq
       av(i) = dot_product(pd(i,:),v(:))
    end do
    return
  END SUBROUTINE K_RES

  SUBROUTINE G_RES(T,Y,G_Y,YPRIME,G_YPRIME,CJ,DELTA,G_DELTA,IRES,RPAR,IPAR,SENPAR,G_SENPAR)
    ! Dummy
   ! pause 'g_res proc'
    return
  end SUBROUTINE G_RES
  ! Direct method
  SUBROUTINE JAC1 (T, U, UPRIME, PD, CJ, RPAR, IPAR, SENPAR, IJAC)    
    use jacobian_module
    use in_data,only:nve,nvg
    ! A = dG/dY + CJ*dG/dYPRIME
    IMPLICIT NONE 
    REAL(KIND=8):: U(nve+nvg), UPRIME(nve+nvg)
    REAL(KIND=8):: PD(nve+nvg,nve+nvg),WK(nve+nvg),REWT(nve+nvg),RPAR(150),SENPAR(*)
    INTEGER, DIMENSION (:) ::  IPAR(5)
    INTEGER :: IRES,IER,NEQ,IJAC
    REAL(KIND=8) :: T,CJ,H
    !pause 'jac1 direct proc'
    call Cal_derv_res(U,PD,cj)    
    return
  END SUBROUTINE JAC1
  ! Krylov method
!!$  SUBROUTINE JAC1 (RESH, IRES, NEQ, T, U, UPRIME, REWT, SAVR,WK, H, CJ, WP, IWP, IER, RPAR, IPAR, SENPAR)    
!!$    use jacobian_module
!!$    use in_data,only:nve,nvg
!!$    ! A = dG/dY + CJ*dG/dYPRIME
!!$    IMPLICIT NONE 
!!$    REAL(KIND=8):: U(nve+nvg), UPRIME(nve+nvg),SAVR(nve+nvg)
!!$    REAL(KIND=8):: PD(nve+nvg,nve+nvg),WK(nve+nvg),REWT(nve+nvg),RPAR(150),SENPAR(*),WP(*)
!!$    INTEGER, DIMENSION (:) ::  IPAR(5),IWP(*)
!!$    INTEGER :: IRES,IER,NEQ
!!$    REAL(KIND=8) :: T,CJ,H
!!$    external RESH
!!$    pause 'jac1 krylov proc'
!!$    call Cal_derv_res(U,PD,cj)    
!!$    return        
!!$  end subroutine JAC1
  !
  subroutine RESH (T, U, UPRIME, CJ, DELTA, IRES, RPAR, IPAR, SENPAR)
    use residual_module
    IMPLICIT NONE 
    REAL(KIND=8),  DIMENSION (:) ::  U(*), UPRIME(*), DELTA(*), RPAR(150),SENPAR(*)
    INTEGER, DIMENSION (:) ::  IPAR(5)                             
    REAL(KIND=8) :: T,CJ  
    INTEGER :: IRES,NEQ
    !pause 'resh proc'
    call Cal_res(uprime,u,delta)
    return
  end subroutine RESH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Preconditioning vector that approx A(-inverse) in jacobian inversion 
  !Otherwise, you may
  !            supply a JAC routine to compute and preprocess any parts of
  !            of the Jacobian matrix  A = dG/dY + CJ*dG/dYPRIME that are
  !            involved in the preconditioner matrix P.
  !            It is to have the form
  ! Direct (dummy)
!!$  SUBROUTINE PSOL (T)
!!$    return
!!$  END SUBROUTINE PSOL
  ! for Krylov
  SUBROUTINE PSOL (NEQ, T, Y, YPRIME, SAVR, WK, CJ, WGHT,WP, IWP, B, EPLIN, IER, RPAR, IPAR, SENPAR)
    ! The PSOL routine must solve linear systems of the form 
    ! P*x = b where P is the left preconditioner matrix.
    ! The right-hand side vector b is in the B array on input, and
    ! PSOL must return the solution vector x in B.
    ! The Y, YPRIME, and SAVR arrays contain the current values
    ! of Y, YPRIME, and the residual G, respectively.      
    !pause 'p_sol proc'
    return        
  end subroutine PSOL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
end module daspk_module




!**************************************************************************************
module maillage
  implicit none
  INTEGER                               :: ne, nbn
  REAL*8,DIMENSION  (:,:),allocatable   :: xyz
  real*8,allocatable, dimension (:)     :: Phi_x,phi_y,vx,vy,vz,conc,temp,nusselt,sher
end module maillage

MODULE in_data  
  implicit none
  integer(kind=4) :: ni,nj,nk,nl,nm,nn,nu,nv,nw,ns,np,nt,nva,nvb,nve,nvg,ntrucmax,compt,ntime 
  real(kind=8)    :: Ra,pi,Ngrav,lewis,poro,sigma
  real(kind=8), ALLOCATABLE, DIMENSION (:) :: tout
END MODULE in_data

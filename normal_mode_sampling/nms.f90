program nms
implicit none
real*8, parameter :: kb=1.38064852d-16, angtocm=1.0d-8, amutog=1.66054d-24
integer, parameter :: na = 3, nm=na*3-6, np = 1000, ndisp=nm/3
!!!!LINE ABOVE NEEDS TO BE ADAPTED.
!na = number of atoms
!nm = number of normal modes, needs to be adapted, e.g. for linear molecules
!np = number of structures to generate at a specific temperature
real*8, dimension (:,:) :: equxyz(na,3),dispxyz(na,3)
real*8, dimension (:,:,:) :: dispv(nm,na,3)
real*8, dimension (:) :: rmass(nm), fconst(nm), dr(nm)
real*8 :: rnum, temp
integer :: ii, jj, kk, ll, dn, ni
character(len=16) :: dummyc, filename
character(len=3), dimension(:) :: aname(na)

aname=(/'O','H','H'/)
!ATOM NAMES IN THE OPTIMIZED STRUCTURE


!ni is the index from which the xyz file names start at that temperature
print*,"ni, temp" 
read*,ni,temp

open(unit=10, file="h2o_opt.xyz") !FILE NAME TO READ THE OPTIMIZED STRUCTURE FROM
open(unit=11, file="h2o_cart_disp.dat") !FILE NAME TO READ THE FORCE CONSTANTS AND DISPLACEMENTS FROM

read(10,*)
read(10,*)
do ii = 1, na
  read(10,*)dummyc, equxyz(ii,:)
end do

!for small frequency modes the fconst can be very small. Setting it to
!some threshold prevents geometries with very large displacements.
do ii =1, nm
  if ( fconst(ii) < 0.05d0 ) fconst(ii) = 0.05d0
end do


ll=1
do ii = 1, ndisp
  read(11,*)
  read(11,*)
  read(11,*)
  read(11,*) dummyc, dummyc, dummyc, rmass(ll:ll+2)
  read(11,*) dummyc, dummyc, dummyc, fconst(ll:ll+2)
  read(11,*)
  read(11,*)
  do jj = 1, na
    read(11,*)dn, dn, dispv(ll,jj,:), dispv(ll+1,jj,:), dispv(ll+2,jj,:)
  end do
  ll=ll+3
end do

do kk = ni, ni+np-1

!  call random_number(rnum)
!  temp=800.0d0*rnum

  do ii = 1, nm
    call random_number(rnum)
    dr(ii)=sqrt(3.0d0*rnum/dfloat(nm)*dfloat(na)*kb*temp/fconst(ii)*angtocm*1.0d3)/angtocm
  end do
  do ii = 1, nm
    call random_number(rnum)
    if (rnum <0.5d0) dr(ii)=-dr(ii)
  end do

  dispxyz= 0.0d0
  do jj = 1, na
    do ii = 1, nm
      dispxyz(jj,:)=dispxyz(jj,:)+dispv(ii,jj,:)/sqrt(rmass(ii))*dr(ii)
    end do
  end do

  write(filename,'(a,i0,a)')"p",kk,".xyz" 
  open(unit=13, file= filename)
  write(13,'(i0)')na
  write(13,'(a)')
  do jj = 1, na
      write(13,'(a,2x,f9.6,2x,f9.6,2x,f9.6)')trim(aname(jj)), equxyz(jj,:)+dispxyz(jj,:)
  end do

end do

end program

module global_param
use healpix_types

implicit none

integer(I4B),parameter :: ncl=3,whichmap=3,whichcl=1,nlheader=80

integer(I4B) :: iter,i,j,maxiter
integer(I4B) :: simnum,nrlz,nstart
integer(I4B) :: nside,npixtot,maplmax,maxnalms

real(dp) :: nullval,fmissval=-1.6375e30
real(dp) :: temp1,temp2
real(dp) :: tnorm=1.d0
real(dp), dimension(2) :: zbounds=0.d0

character :: dtypein*3,dtypeout*3,masktype
character(LEN=200) :: paramfile,pathout
character(LEN=200) :: beamfile,maskfile
character(LEN=200) :: datafname,finprefix,finsuffix,foutprefix,foutsuffix
character :: s1*4
character(len=80),dimension(whichmap)  :: mapheader(1:nlheader)

logical :: anynull,swPWC,swBEAM,swMASK,swDOFS,swSUFFIX
! Note these variables are transparents to all the modules. Hence one needs to be very carefully while making edits to this file.

contains

!########################################################################
subroutine pass_param(mynside,lmax,DOFS)
use healpix_modules

implicit none
character(LEN=200) :: paramfile
integer(I8B) :: mynside,lmax
logical :: DOFS

nside=mynside
DOFS=swDOFS
swPWC=.False.
swBEAM=.False.

maplmax=lmax
npixtot=nside2npix(nside)
maxnalms=((maplmax+1)*(maplmax+2))/2
dtypein="MAP"
dtypeout="MAP"

print*, swDOFS

!call write_minimal_header(mapheader,dtypeout,nside=nside,ordering="RING",polar=.True.)

end subroutine pass_param
!########################################################################

end module global_param

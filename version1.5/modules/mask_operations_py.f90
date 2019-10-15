module process_mask
use global_param
use healpix_types
use healpix_modules

implicit none

real(dp), allocatable, dimension(:,:) :: apomask
complex(dpc),allocatable,dimension(:,:,:) :: maskalm

contains

!########################################################################
subroutine allocate_mask_arrays()

implicit none

allocate(apomask(0:npixtot-1,1),maskalm(1,0:maplmax,0:maplmax))

apomask(:,1)=0.d0
maskalm(1,:,:)=dcmplx(0.d0,0.d0)

end subroutine allocate_mask_arrays

!-------------------------------------------------------------------------

subroutine deallocate_mask_arrays()

deallocate(apomask,maskalm)

end subroutine deallocate_mask_arrays
!########################################################################


!########################################################################
subroutine calc_mask_alm()
implicit none

call map2alm_iterative(nside,2*nside,2*nside,5,apomask,maskalm)
!Print*, "Hello",mymask(:)
!dw8=1.d0
!call map2alm(nside,2*nside,2*nside,apomask(:,1),maskalm,zbounds,dw8)

end subroutine calc_mask_alm
!########################################################################

end module process_mask

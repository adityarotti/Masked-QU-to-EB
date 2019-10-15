module get_TEB

use global_param
use shared_data
use process_mask
use map_operations

contains

subroutine return_TEB(TQU,mymask,mynside,lmax,DOFS,shtiter,npixtot,TEB)
implicit none
integer*8, intent(in) :: mynside,npixtot,shtiter,lmax
logical, intent(in)   :: DOFS
real*8, intent(in)  :: TQU(0:npixtot-1,1:3),mymask(0:npixtot-1)
real*8, intent(out) :: TEB(0:npixtot-1,1:3)

call pass_param(mynside,lmax,DOFS) ; maxiter=shtiter
call allocate_data()
call allocate_mask_arrays()
call return_corrections()

apomask(:,1)=mymask(:)
mapin(:,1)=TQU(:,1) ; mapin(:,2)=TQU(:,2) ; mapin(:,3)=TQU(:,3)
print*, "allocated corrections"
call calc_mask_alm()
print*, "calculated mask alms"
if (swDOFS) then
    swMASK=.False.
    print*, "Doing full sky"
    call convert_TQU2TEB_tilde()
else
    swMASK=.True.
    print*, "Doing masked sky analysis"
    call convert_TQU2TEB()
    call calc_residual()
endif
TEB(:,1)=mapout(:,1) ; TEB(:,2)=mapout(:,2) ; TEB(:,3)=mapout(:,3)
call deallocate_data()
call deallocate_mask_arrays()
end subroutine return_TEB

end module get_TEB

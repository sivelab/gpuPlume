  
  program main

    implicit none

    real*4, allocatable :: u(:,:,:),v(:,:,:),w(:,:,:)
    integer i, j, k
    integer nx, ny, nz

    nx = 2
    ny = 2
    nz = 2

    allocate(u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz))

    do k=1,nz
       do j=1,ny
          do i=1,nx
             u(i,j,k)=i*j*k+1
             v(i,j,k)=i*j*k+2
             w(i,j,k)=i*j*k+3
             print *, "uvw(", i, ",", j, ",", k, ") = ", u(i,j,k), ", ", v(i,j,k), ", ", w(i,j,k)
          end do
       end do
    end do

    open(unit=38,file="QU_velocity.bin",form='unformatted',status="unknown")

    write(38)(((u(i,j,k),i=1,nx),j=1,ny),k=1,nz),   &
         (((v(i,j,k),i=1,nx),j=1,ny),k=1,nz),   &
         (((w(i,j,k),i=1,nx),j=1,ny),k=1,nz)

  end program main

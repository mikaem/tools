      subroutine tovoluviz
      include 'SIZE'
      include 'TOTAL'
      include 'VOLUVIZ'
c
      real dx, dy, dz
      character*100 my_file, common_file
      character*8 nid_str, t_step
      integer my_id, mpi_fh, my_pts_count, comm
      integer n(1), total_n(1), total_offset, my_offset(1), offset(1)
      real pts, ud
      common /volz/ pts(nx*ny*nz*3*2/numproc)
      common /volz/ ud(nx*ny*nz)
      
      write(nid_str, 12) nid
      write(t_step, 12) int(time/dt + 0.5)
   12 format(1x, (i6))

      if (nid==0) then
        call system('mkdir ./'//volfolder)
        call system('mkdir ./'//volfolder//"/"//trim(adjustl(t_step)))
      end if
      call gsync()

      my_file = volfolder//"/"//trim(adjustl(t_step))//
     &  "/"//trim(filename)//trim(adjustl(nid_str))//".dat"

      common_file = volfolder//"/"//trim(adjustl(t_step))//
     &  "/"//trim(filename)//".dat"

      ! Open one file on each processor
      my_id = 111 + nid
!      open(my_id, file=my_file, status='replace')
      OPEN(my_id, file=my_file, form='unformatted', status='replace')
      
      dx = dxl/real(nx-1)
      dy = dyl/real(ny-1)
      dz = dzl/real(nz-1)

      ! Divide last dimension between processors
      my_nz = nz/np
      my_offset(1) = nid*my_nz + 1
      if (nid.eq.np-1) then
         my_nz = my_nz+mod(nz,np)
      end if
      my_pts_count = nx*ny*my_nz

      do k = my_offset(1),my_offset(1)+my_nz-1
      do j = 1,ny
      do i = 1,nx
        pts(i + (j-1)*nx + (k-my_offset(1))*nx*ny) = 
     &      dx0 + real(i-1)*dx
        pts(i + (j-1)*nx + (k-my_offset(1))*nx*ny + my_pts_count) = 
     &      dy0 + real(j-1)*dy
        pts(i + (j-1)*nx + (k-my_offset(1))*nx*ny + 2*my_pts_count) = 
     &      dz0 + real(k-1)*dz
      end do
      end do
      end do

      write(my_id) my_pts_count
      write(my_id) nx, ny, my_nz
      write(my_id) vars

      ! This part needs to be tweaked if one wants to store more
      ! variables than U, V, W and P
      call intpts_setup(1e-7, ih)
      ud(1:my_pts_count) = 0.
      call intpts(vx, 1, pts, my_pts_count, ud, .true., .true., ih)
      write(my_id) real(ud(1:my_pts_count), kind=4)
      ud(1:my_pts_count) = 0.
      call intpts(vy, 1, pts, my_pts_count, ud, .true., .true., ih)
      write(my_id) real(ud(1:my_pts_count), kind=4)
      ud(1:my_pts_count) = 0.
      call intpts(vz, 1, pts, my_pts_count, ud, .true., .true., ih)
      write(my_id) real(ud(1:my_pts_count), kind=4)
      ud(1:my_pts_count) = 0.
      call intpts(pr, 1, pts, my_pts_count, ud, .true., .true., ih)
      write(my_id) real(ud(1:my_pts_count), kind=4)
      ud(1:my_pts_count) = 0.
      call intpts(t, 1, pts, my_pts_count, ud, .true., .true., ih)
      write(my_id) real(ud(1:my_pts_count), kind=4)

      call intpts_done(ih)

      close(my_id)
      end subroutine

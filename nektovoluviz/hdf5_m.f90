
module hdf5_m
    !
    !
    ! Some useful routines
    !
    !

    use hdf5

    integer, parameter :: WP = 4

    interface append
        module procedure append0, appendi0
    end interface

    interface write_scalar
        module procedure write_scalar_double, write_scalar_int
    end interface

contains

   subroutine write_scalar_double(h5_id, my_disp, tot_nx, nx, name, data)
        !
        ! Write a double scalar field collectively to file
        !
        implicit none
        integer, intent(in) :: h5_id, my_disp
        character(len=*), intent(in) :: name
        real(WP), intent(in) :: data(:)
        integer, intent(in) :: tot_nx(:), nx(:)
        integer :: plist_id, dset_id, filespace, memspace, ierr
        integer(8), dimension(3) :: sz(3),mz(3),offset(3)


        sz(1:3) = tot_nx(1:3)
!        sz(1) = tot_nx(1)*tot_nx(2)*tot_nx(3)
        call h5screate_simple_f(3, sz, filespace, ierr)

         mz(1:3) = nx(1:3)
!        mz(1) = nx(1)*nx(2)*nx(3)

        call h5screate_simple_f(3, mz, memspace, ierr)
        !
        ! Create the dataset with default properties.
        !
        call h5dcreate_f(h5_id, name, H5T_NATIVE_REAL, filespace, &
             dset_id, ierr)

        call h5sclose_f(filespace, ierr)      

        ! 
        ! Select hyperslab in the file.
        offset(1) = 0
        offset(2) = 0
        offset(3) = my_disp
        call h5dget_space_f(dset_id, filespace, ierr)
        call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, mz, ierr)
        !
        ! Create property list for collective dataset write
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr) 
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)

        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, data, sz, ierr, &
             file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

        call h5sclose_f(filespace, ierr)
        call h5sclose_f(memspace, ierr)
        call h5dclose_f(dset_id, ierr)
        call h5pclose_f(plist_id, ierr)


    end subroutine write_scalar_double

    subroutine write_scalar_int(h5_id, total_len, my_len, my_disp, name, data)
        !
        ! Write an integer scalar field collectively to file
        !
        implicit none
        integer :: h5_id,total_len, my_len, my_disp
        character(len=*), intent(in) :: name
        integer :: data(:)
        !
        integer :: plist_id, dset_id, filespace, memspace,ierr
        integer(8), dimension(1) :: sz(1),offset(1)

        sz(1) = total_len
        call h5screate_simple_f(1, sz, filespace, ierr)
        sz(1)=my_len
        call h5screate_simple_f(1, sz, memspace, ierr)
        !
        ! Create the dataset with default properties.
        call h5dcreate_f(h5_id, trim(name), H5T_NATIVE_INTEGER, filespace, &
             dset_id, ierr)
        call h5sclose_f(filespace, ierr)

        call h5screate_simple_f(1, sz, memspace, ierr) 
        ! 
        ! Select hyperslab in the file.
        offset(1) = my_disp
        call h5dget_space_f(dset_id, filespace, ierr)
        call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, sz, ierr)
        !
        ! Create property list for collective dataset write
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr) 
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
        !
        ! Write the dataset collectively. 
        sz(1) = total_len
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, sz, ierr, &
             file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

        call h5sclose_f(filespace, ierr)
        call h5sclose_f(memspace, ierr)
        call h5dclose_f(dset_id, ierr)
        call h5pclose_f(plist_id, ierr)

    end subroutine write_scalar_int

    subroutine annote(id, str)
        !
        !   Annote an object (group/dataset) id with str,
        !   using object attributes
        !
        !
        implicit none
        integer, intent(in) :: id
        character(len=*), intent(in) :: str
        integer :: str_id, ann_space, title_id, l, ierr
        integer(8) :: one(1) = (/1_WP/)
        integer(8) :: len_str 
        !
        len_str = len_trim(str)
        call h5tcopy_f(H5T_NATIVE_CHARACTER, str_id, ierr)
        call h5tset_size_f(str_id, len_str, ierr)
        call h5screate_f(H5S_SCALAR_F, ann_space, ierr)
        call h5acreate_f(id, 'Description', str_id, ann_space, title_id, ierr)
        l=len_trim(str)
        call h5awrite_f(title_id, str_id, str, one, ierr)
        !
        call h5tclose_f(str_id, ierr)
        call h5sclose_f(ann_space, ierr)
        call h5aclose_f(title_id, ierr)
    end subroutine annote


    subroutine write_with_tab(fileid,s,n)
        implicit none
        integer, intent(in) :: fileid,n
        character(len=*) :: s
        integer :: i
        do i=1,n
            s=char(9)//s
        end do
        write(unit=fileid, fmt='(a)') trim(s)

    end subroutine write_with_tab

    subroutine created(fid, r, d, name, desc, fmt)
        !
        !   Create dataset for arrays of rank r and shape d, with UNLIMITED
        !   size for the r+1 dimension (Note: d is unused when r=0)
        !   The format of the array is described in fmt (=H5T_NATIVE_REAL, H5T_NATIVE_INTEGER) 
        !   The default is H5T_NATIVE_REAL
        !
        implicit none
        integer, intent(in) :: fid, r, d(:)
        character(len=*), intent(in):: name
        character(len=*), intent(in), optional :: desc
        integer, optional :: fmt
        !
        integer(8), dimension(r+1) :: dims, maxdims
        integer :: space_id, cprop_id, did,fmt_m, rank, ierr

        if(.not.present(fmt)) then ! Default is a double dataset
            fmt_m=H5T_NATIVE_REAL
        else
            fmt_m=fmt
        end if

        !  define file dataspace
        rank = r+1
        if( r .gt. 0 ) then
            dims(1:rank-1) = d(1:r)
        end if
        dims(rank) = 0 ! Initial time dimension
        !
        maxdims = dims
        maxdims(rank) = H5S_UNLIMITED_f
        call h5screate_simple_f(rank, dims, space_id, ierr, maxdims)
        !
        !  Define chunking for last dimension
        call h5pcreate_f(H5P_DATASET_CREATE_F, cprop_id, ierr)
        dims(rank) = 4
        call h5pset_chunk_f(cprop_id, rank, dims, ierr)
        !
        ! Create data set
        call h5dcreate_f(fid, name, fmt_m, space_id, did, ierr, cprop_id)
        if( present(desc) ) then
            call annote(did, desc)
        end if
        !
        call h5sclose_f(space_id, ierr)
        call h5pclose_f(cprop_id, ierr)
        call h5dclose_f(did, ierr)
    end subroutine created

    subroutine getsize(fid, name, n)
        !
        !   get last dim of data set
        !
        implicit none
        integer, intent(in) :: fid
        integer, intent(out) :: n
        character(len=*), intent(in) :: name
        !
        integer(8) :: dims(7), maxdims(7)
        integer :: dspace_id, rank, ierr, did
        !
        call h5dopen_f(fid, name, did, ierr) ! data id
        call h5dget_space_f(did, dspace_id, ierr)
        call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, rank)
        n = dims(rank)
        call h5dclose_f(did, ierr)
        call h5sclose_f(dspace_id, ierr)
    end subroutine getsize

    subroutine append0(fid, name, scal)
        !
        !   add a scalar at the end of dataset
        !
        implicit none
        integer, intent(in) :: fid
        character(len=*), intent(in) :: name
        double precision, intent(in) :: scal
        !
        double precision :: array(1)
        integer :: did, dspace_id, memspace_id
        integer(8), dimension(1) :: dims, maxdims, starts, counts, ddims
        integer :: rank, ierr
        !
        !  get dims of dataset
        call h5dopen_f(fid, name, did, ierr)
        call h5dget_space_f(did, dspace_id, ierr)
        call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, ierr)
        rank = ierr
        call h5sclose_f(dspace_id, ierr)
        if( rank .ne. 1 ) then
            write(*, '(a,a)') "data shape mismatch for", name(1:len_trim(name))
            stop
        end if
        !
        !  extend the dataset in the last dimension
        starts(1) = dims(1)
        counts(1) = 1
        dims(1) = dims(1) + 1
        call h5dextend_f(did, dims, ierr)
        call h5dget_space_f(did, dspace_id, ierr)
        !
        !  memory dataspace
        ddims(1) = 1
        call h5screate_simple_f(1, ddims, memspace_id, ierr)
        !
        !  write to the end of dataset
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, starts, counts, &
             &                     ierr)
        array(1) = scal
        call h5dwrite_f(did, H5T_NATIVE_REAL, array, ddims, ierr, &
             &  memspace_id, dspace_id)
        !
        call h5dclose_f(did, ierr) 
        call h5sclose_f(dspace_id, ierr)
        call h5sclose_f(memspace_id, ierr)
    end subroutine append0

    subroutine appendi0(fid, name, scal)
        !
        !   add an int at the end of dataset
        !
        implicit none
        integer, intent(in) :: fid
        character(len=*), intent(in) :: name
        integer, intent(in) :: scal
        !
        integer :: array(1),rank, ierr
        integer :: did, dspace_id, memspace_id
        integer(8), dimension(1) :: dims, maxdims, starts, counts, ddims
        !
        !  get dims of dataset
        call h5dopen_f(fid, name, did, ierr)
        call h5dget_space_f(did, dspace_id, ierr)
        call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, ierr)
        rank = ierr
        call h5sclose_f(dspace_id, ierr)
        if( rank .ne. 1 ) then
            write(*, '(a,a)') "data shape mismatch for", name(1:len_trim(name))
            stop
        end if
        !
        !  extend the dataset in the last dimension
        starts(1) = dims(1)
        counts(1) = 1
        dims(1) = dims(1) + 1
        call h5dextend_f(did, dims, ierr)
        call h5dget_space_f(did, dspace_id, ierr)
        !
        !  memory dataspace
        ddims(1) = 1
        call h5screate_simple_f(1, ddims, memspace_id, ierr)
        !
        !  write to the end of dataset
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, starts, counts, &
             &                     ierr)
        array(1) = scal
        call h5dwrite_f(did, H5T_NATIVE_INTEGER, array, ddims, ierr, &
             &  memspace_id, dspace_id)
        !
        call h5dclose_f(did, ierr) 
        call h5sclose_f(dspace_id, ierr)
        call h5sclose_f(memspace_id, ierr)
    end subroutine appendi0

    subroutine graceful_exit(n)
        integer, intent(in) :: n
        if (n.lt.0) then
          call MPI_Finalize(n)
          stop
        end if
    end subroutine


end module hdf5_m

PROGRAM NEKTOVOLUVIZ

     USE HDF5_M

     IMPLICIT NONE

     INCLUDE 'mpif.h'

     CHARACTER(LEN=100) :: filename, h5filename  ! File name
     CHARACTER(LEN=100) :: my_filename
     INTEGER        :: fnamelen      ! File name length
     INTEGER :: file_id, data_group_id, parameter_group_id     ! File identifier
     INTEGER :: dset_id       ! Dataset identifier
     INTEGER :: filespace     ! Dataspace identifier in file
     INTEGER :: plist_id      ! Property list identifier
     INTEGER, ALLOCATABLE :: t_group_id(:,:)
     INTEGER, ALLOCATABLE :: var_group_id(:)

     REAL*4, ALLOCATABLE :: data(:) 
     INTEGER :: rank = 1 ! Dataset rank

     INTEGER :: error, error_n  ! Error flags
     INTEGER :: i, j, k, l, num_vars
     LOGICAL :: fileexists
     !
     ! MPI definitions and calls.
     !
     INTEGER :: mpierror       ! MPI error flag
     INTEGER :: comm, info
     INTEGER :: mpi_size, mpi_rank, my_no_count(1),my_no_disp(1), my_offset(1)
     integer :: n(1), nx(3), total_nx(3), total_n(1), my_id, data_dims(1)
     character*8 mpi_rank_str
     integer, allocatable :: pts_on_proc(:), offset(:)
     character*8 :: num_dirs_str, filenr
     integer :: num_dirs, tstep
     character*8, allocatable :: list_of_dirs(:) 
     character(len=100) :: vars
     character*8, allocatable :: varss(:)
     character*1 :: s
     character*8 :: s8
     character*60 :: create_infile(10)

     comm = MPI_COMM_WORLD
     info = MPI_INFO_NULL
     CALL MPI_INIT(mpierror)
     CALL MPI_COMM_SIZE(comm, mpi_size, mpierror)
     CALL MPI_COMM_RANK(comm, mpi_rank, mpierror)

     create_infile(1:10) = (/"import os                                              ", &
                             "ldir = os.listdir(os.path.join(os.getcwd(), 'voluviz'))", &
                             "for i, l in enumerate(ldir):                           ", &
                             "    ldir[i] = eval(l)                                  ", &
                             "ldir.sort()                                            ", &
                             "f = open('nektovoluviz.in', 'w')                       ", &
                             "f.write(str(len(ldir)) + '\n')                         ", &
                             "for l in ldir:                                         ", &
                             "    f.write(str(l) + '\n')                             ", &
                             "f.close()                                              "/)

     !call h5eset_auto_f(0, mpierror)

     call getarg(IARGC(), filename)
     if (mpi_rank==0) then
        print *, "Writing hdf5 datafile for ", filename
        ! Hack. Could probably use flibs but that requires som additional libraries installed
        open(200, file='create_infile.py', status='replace')
        write(200, "(a60)") (create_infile(i), i = 1,10)
        close(200)
        call system('python create_infile.py')
        open(201, file='nektovoluviz.in', status='old')
        read(201, "(i4)") num_dirs
        allocate(list_of_dirs(num_dirs))
        read(201, "(a8)") (list_of_dirs(k), k =1, num_dirs)
        write(*,*) num_dirs, list_of_dirs
     end if
     call mpi_barrier(comm, mpierror)
     call mpi_bcast(num_dirs, 1, MPI_INTEGER, 0, comm, error)
     if (.not. mpi_rank==0) then
        allocate(list_of_dirs(num_dirs))
     end if
     call mpi_barrier(comm, mpierror)
     call MPI_Bcast(list_of_dirs, num_dirs, MPI_REAL8, 0, comm, error)

     write(mpi_rank_str, 12) mpi_rank
   12 format(1x, i4)   

     !
     ! Initialize FORTRAN interface
     !
     CALL h5open_f(error)
     !
     ! Setup file access property list with parallel I/O access.
     !
     CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
     CALL h5pset_fapl_mpio_f(plist_id, comm, info, error)
     h5filename = trim(filename)//".h5"
     INQUIRE(FILE = h5filename, EXIST = fileexists)
     call h5fcreate_f(h5filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
     allocate (pts_on_proc(mpi_size))  
     allocate (offset(mpi_size))
     
     do i = 1, num_dirs
     if (mpi_rank.eq.0) then
        write(*,*) 'Processing ',list_of_dirs(i), 'Progress ', i, num_dirs
     end if
 
     my_filename = "voluviz/"//trim(adjustl(list_of_dirs(i)))//"/"//trim(filename)//trim(adjustl(mpi_rank_str))//".dat"
     
     my_id = 112 + mpi_rank
!      open(my_id,file=my_filename,status='old',form='formatted')
!      read(my_id, 14) n(1)
!      read(my_id, 16) nx(1), nx(2), nx(3)
!      read(my_id, 18) vars

     open(my_id,file=my_filename,status='old',form='unformatted')
     read(my_id) n(1)
     read(my_id) nx(1), nx(2), nx(3)
     read(my_id) vars
     write(*,*) 'Total number of points ', n(1), nx, mpi_rank, vars

  14 format(1x,i10)
  16 format(1x,3i10)
  18 format(1x,a100)

     if (i.eq.1) then
        k = 1
        do j = 1, 100
          s = vars(j:j)
          if (s.eq.',') then
            k = k+1
          end if
        end do

        allocate(varss(k))
        num_vars = k
        k = 1
        l = 1
        do j = 1, len(trim(adjustl(vars)))
          s = vars(j:j)
          if (s.eq.',') then
            varss(l) = vars(k:j-1)
            k = j+1
            l = l+1
          end if  
        end do
        varss(l) = vars(k:j)
        allocate(var_group_id(num_vars))
     end if

     call mpi_allreduce(nx, total_nx, 3, MPI_INTEGER, MPI_SUM, comm, error)
     total_nx(1:2) = nx(1:2)

     call mpi_allreduce(n, total_n, 1, MPI_INTEGER, MPI_SUM, comm, error)
     
     allocate (data(num_vars*n(1)))
     call MPI_ALLGATHER(n,1,MPI_INTEGER, &
             pts_on_proc,1,MPI_INTEGER,comm,error)

     call MPI_ALLGATHER(nx(3),1,MPI_INTEGER, &
             offset,1,MPI_INTEGER,comm,error)

     call mpi_barrier(comm, error)

     my_no_count(1) = n(1)
     my_offset(1) = sum(offset(1:mpi_rank))
 
!      read(my_id, 15) (data(l), l= 1, num_vars*n(1))
!   15 format(1x,g14.7)
     do j = 1,num_vars
       read(my_id) data(((j-1)*n(1)+1):(j*n(1)))
     end do

     close(my_id)

     if (i.eq.1) then 
         call h5gcreate_f(file_id, "Data", data_group_id, error)
         do j = 1, num_vars
            call h5gcreate_f(data_group_id, varss(j), var_group_id(j), error)
            call h5gclose_f(var_group_id(j), error) 
         end do
         call h5gclose_f(data_group_id, error)
         call h5gcreate_f(file_id, "Parameters", parameter_group_id, error)
         data_dims(1)=1 
         call created(file_id, 0, data_dims,"/Parameters/Times",&
               "Array of true timesteps for the stored data", fmt=H5T_NATIVE_INTEGER)
         call h5gclose_f(parameter_group_id, error)

     end if

     call h5gopen_f(file_id,"Data", data_group_id, error)

     ! Store data for all variables on given timestep 
     ! Use 0, 1, 2 numbering because voluviz requires it
     ! Store true timestep in attribute Times
     do k = 1, num_vars
        call h5gopen_f(data_group_id, varss(k), var_group_id(k), error) 
        write(filenr, "(i4)") i-1 
        call write_scalar(var_group_id(k), my_offset(1), total_nx, nx, &      
               trim(adjustl(filenr)) , data(((k-1)*n(1)+1):(k*n(1))))

        call h5gclose_f(var_group_id(k), error)
     end do

     deallocate(data)
     call h5gclose_f(data_group_id, error)

     ! Store the true timestep
     read(list_of_dirs(i), '(i6)') tstep
     call append(file_id, "/Parameters/Times", tstep)

     end do

     call h5pclose_f(plist_id, error)
     call h5fclose_f(file_id, error)
     deallocate(varss)
     deallocate(var_group_id)
     deallocate (pts_on_proc)  
     deallocate (offset)

     CALL h5close_f(error)

     CALL MPI_FINALIZE(mpierror)

END PROGRAM NEKTOVOLUVIZ

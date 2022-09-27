!--------------------------------------------------------------------
PROGRAM fv3interp2latlon

   use namelist_module
   use tile_module
   use latlon_module
   use gaussian_module
   use fv_grid_utils_module

   use mpi

   IMPLICIT NONE

   type tiletype
      type(tilegrid), dimension(6) :: tile
   end type tiletype

   type(tilespec_type), dimension(6)    :: spec
   type(tiletype), dimension(max_types) :: types
   type(latlongrid)                     :: latlon
   type(gaussiangrid)                   :: gaussian
   type(fv_grid_type), dimension(6)     :: gridstruct

   integer :: n, nm
   logical :: last
   real, dimension(:), allocatable :: lon, lat
   integer :: num_procs, myrank, ierr
   integer :: mymembers, member
   character(len=10) :: memstr

   print *, 'File: ', __FILE__, ', line: ', __LINE__

  !Initialize MPI.
   call MPI_Init(ierr)
  !Find out the number of processes available.
   call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr)
  !Determine this process's rank.
   call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)

   print *, 'num_procs = ', num_procs, ', myrank = ', myrank

   call read_namelist('input.nml')

   print *, 'total_members = ', total_members

   mymembers = (total_members + num_procs - 1) / num_procs

   do nm = 1, mymembers
      member = myrank*mymembers + nm
      if(member > total_members) then
         member = 0
         cycle
      end if
      write(memstr, fmt='(a, i3.3)') 'mem', member
      print *, 'myrank: ', myrank, ', member: ', member, ', memstr: ', trim(memstr)
   end do

#if 0
   if(num_types < 1) then
      print *, 'num_types must great than 0. eg. must have at least 1 type.'
      stop 111
   end if

   print *, 'File: ', __FILE__, ', line: ', __LINE__

   if(use_uv_directly) then
      call initialize_tilespec(spec, trim(griddirname), trim(grid_type))
   end if

   print *, 'File: ', __FILE__, ', line: ', __LINE__

   if(use_gaussian_grid) then
      call initialize_gaussiangrid(gaussian_grid_file, &
                                   nlon, nlat, nlev, nilev, npnt, gaussian)
   else
      call initialize_latlongrid(nlon, nlat, npnt, latlon)
   end if

   print *, 'File: ', __FILE__, ', line: ', __LINE__

   do n = 1, num_types
      print *, 'dirname: <', trim(dirname), &
               '>, data_types(', n, ') = <', trim(data_types(n)), '>'
      call initialize_tilegrid(types(n)%tile, trim(dirname), trim(data_types(n)))

      if(.not. use_gaussian_grid) then
         if(trim(data_types(n)) == 'fv_core.res.tile') then
            latlon%nlev = types(n)%tile(1)%nz
         else if(trim(data_types(n)) == 'sfc_data.tile') then
            latlon%nlay = types(n)%tile(1)%nz
         end if
      end if
   end do

   print *, 'File: ', __FILE__, ', line: ', __LINE__
  !print *, 'latlon%nlev: ', latlon%nlev, ', latlon%nlay: ', latlon%nlay
   print *, 'generate_weights = ', generate_weights

   if(generate_weights) then
      print *, 'File: ', __FILE__, ', line: ', __LINE__
      print *, 'use_gaussian_grid = ', use_gaussian_grid
      if(use_gaussian_grid) then
         call generate_weight4gaussian(types(1)%tile, gaussian)
         call write_gaussiangrid(gaussian, wgt_flnm)
      else
         call generate_weight(types(1)%tile, latlon)
         call write_latlongrid(latlon, wgt_flnm)
      end if
   else
      print *, 'File: ', __FILE__, ', line: ', __LINE__
      print *, 'wgt_flnm: ', trim(wgt_flnm)
      if(use_gaussian_grid) then
         call read_weights4gaussian(gaussian, wgt_flnm)
      else
         call read_weights(latlon, wgt_flnm)
      end if

      print *, 'File: ', __FILE__, ', line: ', __LINE__
      print *, 'use_uv_directly: ', use_uv_directly
      if(use_uv_directly) then
         do n = 1, 6
            call grid_utils_init(spec(n), gridstruct(n), &
                                 types(1)%tile(n)%nx, types(1)%tile(n)%ny)
         end do
      end if

      print *, 'File: ', __FILE__, ', line: ', __LINE__
      print *, 'num_types: ', num_types

      do n = 1, num_types
         last = (n == num_types)
         print *, 'n = ', n
         print *, 'last = ', last
         if(use_gaussian_grid) then
            call generate_header4gaussian(n, types(n)%tile, gaussian, &
                                 trim(data_types(n)), output_flnm, last)
         else
            call generate_header(n, types(n)%tile, latlon, &
                                 trim(data_types(n)), output_flnm, last)
         end if
      end do

      print *, 'File: ', __FILE__, ', line: ', __LINE__
      print *, 'num_types: ', num_types

      do n = 1, num_types
         if(use_gaussian_grid) then
            call interp2gaussiangrid(trim(data_types(n)), spec, gridstruct, &
                                     types(n)%tile, gaussian)
         else
            call interp2latlongrid(trim(data_types(n)), spec, gridstruct, &
                                   types(n)%tile, latlon)
         end if
      end do

      print *, 'File: ', __FILE__, ', line: ', __LINE__
      print *, 'use_uv_directly: ', use_uv_directly

      if(use_uv_directly) then
         do n = 1, 6
            call grid_utils_exit(gridstruct(n))
         end do
      end if
   end if

   print *, 'File: ', __FILE__, ', line: ', __LINE__

   do n = 1, num_types
      call finalize_tilegrid(types(n)%tile)
   end do

   print *, 'File: ', __FILE__, ', line: ', __LINE__

   if(use_gaussian_grid) then
      call finalize_gaussiangrid(gaussian)
   else
      call finalize_latlongrid(latlon)
   end if

   print *, 'File: ', __FILE__, ', line: ', __LINE__

   if(use_uv_directly) then
      call finalize_tilespec(spec)
   end if
#endif

   print *, 'File: ', __FILE__, ', line: ', __LINE__

  !End MPI
   call MPI_Finalize(ierr)

END PROGRAM fv3interp2latlon


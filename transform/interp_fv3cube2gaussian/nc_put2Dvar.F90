!---------------------------------------------------

subroutine nc_put2Dvar(ncid, var_name, v2d, nrec, &
                       m1s, m1e, m2s, m2e)

   use netcdf

   implicit none
 
   integer, intent(in) :: ncid, nrec
   integer, intent(in) :: m1s, m1e, m2s, m2e
   character(len=*), intent(in) :: var_name
   real, dimension(m1s:m1e, m2s:m2e), intent(in) :: v2d

 ! Variable id
   integer :: varid

 ! Return status
   integer :: status

   integer, dimension(3) :: start, count

   start(1) = m1s
   start(2) = m2s
   start(3) = nrec

   count(1) = m1e - m1s + 1
   count(2) = m2e - m2s + 1
   count(3) = 1

   status = nf90_inq_varid(ncid, var_name, varid)
   if(status /= nf90_noerr) then 
      write(unit=0, fmt='(3a)') "Problem to get id for: <", trim(var_name), ">.", &
                                "Error status: ", trim(nf90_strerror(status))
      write(unit=0, fmt='(3a, i4)') &
          "Stop in file: <", __FILE__, ">, line: ", __LINE__
      stop
   end if

   status = nf90_put_var(ncid,varid,v2d,start=start,count=count)
   if(status /= nf90_noerr) then
          write(unit=0, fmt='(3a)') "Problem to write variable: <", trim(var_name), ">.", &
                                 "Error status: ", trim(nf90_strerror(status))
      write(unit=0, fmt='(3a, i4)') &
          "Stop in file: <", __FILE__, ">, line: ", __LINE__
      stop
   end if

end subroutine nc_put2Dvar

!---------------------------------------------------

subroutine nc_put2Dvar0(ncid, var_name, v2d, m1s, m1e, m2s, m2e)

   use netcdf

   implicit none
 
   integer, intent(in) :: ncid, m1s, m1e, m2s, m2e
   character(len=*), intent(in) :: var_name
   real, dimension(m1s:m1e, m2s:m2e), intent(in) :: v2d

 ! Variable id
   integer :: varid

 ! Return status
   integer :: status

   integer, dimension(2) :: start, count

   start(1) = m1s
   start(2) = m2s

   count(1) = m1e - m1s + 1
   count(2) = m2e - m2s + 1

   status = nf90_inq_varid(ncid, var_name, varid)
   if(status /= nf90_noerr) then 
      write(unit=0, fmt='(3a)') "Problem to get id for: <", trim(var_name), ">.", &
                                "Error status: ", trim(nf90_strerror(status))
      write(unit=0, fmt='(3a, i4)') &
          "Stop in file: <", __FILE__, ">, line: ", __LINE__
      stop
   end if

   status = nf90_put_var(ncid,varid,v2d,start=start,count=count)
   if(status /= nf90_noerr) then
          write(unit=0, fmt='(3a)') "Problem to write variable: <", trim(var_name), ">.", &
                                 "Error status: ", trim(nf90_strerror(status))
      write(unit=0, fmt='(3a, i4)') &
          "Stop in file: <", __FILE__, ">, line: ", __LINE__
      stop
   end if

end subroutine nc_put2Dvar0


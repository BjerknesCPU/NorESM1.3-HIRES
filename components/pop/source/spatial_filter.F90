module spatial_filter

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


!BOP
! !MODULE: spatial_filter
! !DESCRIPTION:
!Routines to apply spatial filter to specified variables 

!filter type:
!   1 = boxcar latlon w/init (DEFAULT)                                                                        
!   2 = boxcar radial w/init                                                                                 
!   3 = loess radial w/init
!   4 = loess boxcar w/init                                                                                     
!
!namelist params:
!   l_sst_filter (= .true. to turn on)
!   filter_type (defined above)
!   filter_x_span (specify in km the distance from center for the longitude for the lat-lon box)
!   filter_y_span (specify in km the distance from center for the latitude for the lat-lon box)
!   filter_r_span (specify in km for a radial smoothing region)
!   filter_print (for debugging)
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


! !USES:

  use POP_KindsMod
  use POP_ErrorMod
  use POP_CommMod
  use POP_FieldMod
  use POP_GridHorzMod
  use POP_HaloMod
  use POP_IOUnitsMod
  use POP_MCT_vars_mod

  use communicate
  use kinds_mod
  use broadcast
  use gather_scatter
  use constants, only: field_loc_center, field_type_scalar, c0, blank_fmt, ndelim_fmt, &
       radius, pi2, pih, pi
  use grid, only : KMT, KMT_G, HUS, HUW, TLAT, TLON
  use domain
  use domain_size
  use io_types
  use exit_mod, only: sigAbort, exit_pop, flushm
  use broadcast
  use blocks
  use timers, only: timer_start, timer_stop, get_timer

  implicit none
  save


  logical (POP_logical), public :: &
       l_sst_filter   ! true if we filter sst

  integer (POP_i4), public :: &
       filter_type     !type of filter


! !PUBLIC MEMBER FUNCTIONS:
  public :: apply_filter
  public :: init_filter

!  module private data

  integer (int_kind), parameter :: &
       max_n_filter_regions = 9

  integer (POP_i4), dimension(:,:), allocatable :: &
       kmt_global   

  real (r8), dimension(:,:), allocatable :: & 
       tlat_global, tlon_global, &
       xcoord_global, ycoord_global, zcoord_global           

  real (POP_r8) :: &
       filter_x_span, filter_y_span, filter_r_span

  real (POP_r8), dimension(max_n_filter_regions) :: &
       region_x_spans, region_y_spans, &
       region_min_lons, region_max_lons, region_min_lats, &
       region_max_lats

  character (POP_charLength), public  ::  filter_outfile  

  integer (int_kind) :: n_filter_regions
  integer (int_kind) :: par_jstart, par_jend, par_istart, par_iend
  logical (log_kind) :: first_print, par_filter, master_first_print, filter_print

  
  integer (int_kind) :: &
        timer_filter_tot

  !for initialization
  integer (int_kind), dimension(:), allocatable :: reg_starts, reg_i, reg_j
  
  integer (int_kind) :: loess_size

  integer (int_kind) :: M_divide

!EOP


!***********************************************************************

 contains

!***********************************************************************


   subroutine init_filter()

     ! !DESCRIPTION:  
     !initialize the spatial filter

  
     integer (int_kind) :: nml_error, nr, is_error
     character (POP_charLength) :: message

     !Namelist stuff
     namelist /filter_nml/l_sst_filter, filter_type, & 
          filter_x_span, filter_y_span, filter_r_span, &
          filter_print, filter_outfile, n_filter_regions, &
          region_x_spans, region_y_spans, &
          region_min_lats, region_max_lats, region_min_lons, &
          region_max_lons
     
     l_sst_filter = .false.
     filter_type = 1 !boxcar lat-lon
     filter_x_span = 500.0 !km
     filter_y_span = 500.0
     filter_r_span = 500.0
     filter_print = .false.

     filter_outfile = 'unknown_filter_outfile'

     !filter regions (can only regional filter with lat/lon region - not radial)
     n_filter_regions = 0
     M_divide = 1;

     !initialize
     par_filter = .false.
     first_print = .false. !if not parallel, this gives more information
     master_first_print = .false.

     do nr = 1, max_n_filter_regions
        region_max_lats(nr) = 0.0
        region_min_lats(nr) = 0.0
        region_max_lons(nr) = 0.0
        region_min_lons(nr) = 0.0
        region_x_spans(nr) = 0.0
        region_y_spans(nr) = 0.0
     enddo

     if (my_task == master_task) then
        open (nml_in, file=nml_filename, status='old',iostat=nml_error)
        if (nml_error /= 0) then
           nml_error = -1
        else
           nml_error =  1
        endif
        do while (nml_error > 0)
           read(nml_in, nml=filter_nml,iostat=nml_error)
        end do
        if (nml_error == 0) close(nml_in)
     endif
     
     call broadcast_scalar(nml_error, master_task)
     if (nml_error /= 0) then
        call exit_POP(sigAbort,'ERROR reading filter_nml')
     endif
     
     call broadcast_scalar(n_filter_regions, master_task)
     if (n_filter_regions > max_n_filter_regions) then
        call exit_POP(sigAbort, "FATAL ERROR: too many n_filter_regions")
     endif

     if (my_task == master_task) then

        if (filter_type > 4) filter_type = 1

        write(stdout,blank_fmt)
        write(stdout,ndelim_fmt)
        write(stdout,blank_fmt)
        if (l_sst_filter) then
           write(stdout,'(a27)') '**Filtering the SST         '
           write(stdout,'(a26,i2,a1)') '    SST smoothed with type', &
                filter_type,'.'
           select case (filter_type)
              case (1)
                 write(stdout,'(a32)') '     (Boxcar rectangular)       '
              case (2) 
                 write(stdout,'(a32)') '     (Boxcar radial)           '
              case (3)
                 write(stdout,'(a32)') '     (Loess radial)             '
              case(4)
                 write(stdout,'(a32)') '     (Loess rectangular)        '
           end select   

           write(stdout,'(a26,f7.3)')    '                 x-span = ', &
                filter_x_span
           write(stdout,'(a26,f7.3)')    '                 y-span = ', &
                filter_y_span 
           write(stdout,'(a26,f7.3)')    '                 r-span = ', &
                filter_r_span
           write(stdout,'(a26, i2)')     '       n_filter_regions = ', &
                n_filter_regions

           if (n_filter_regions > 0) then
              do nr = 1, n_filter_regions
                 write(stdout,'(a26,i2,a4,f7.3)') '           region_x_spans(', & 
                      nr,  ') = ', region_x_spans(nr)
              enddo
              do nr = 1, n_filter_regions
                 write(stdout,'(a26,i2,a4,f7.3)') '           region_y_spans(', & 
                      nr,  ') = ', region_y_spans(nr)
              enddo
              do nr = 1, n_filter_regions
                 write(stdout,'(a26,i2,a4,f7.3)') '          region_min_lats(', & 
                      nr,  ') = ', region_min_lats(nr)
              enddo
              do nr = 1, n_filter_regions
                 write(stdout,'(a26,i2,a4,f7.3)') '          region_max_lats(', & 
                      nr,  ') = ', region_max_lats(nr)
              enddo
              do nr = 1, n_filter_regions
                 write(stdout,'(a26,i2,a4,f7.3)') '          region_min_lons(', & 
                      nr,  ') = ', region_min_lons(nr)
              enddo
              do nr = 1, n_filter_regions
                 write(stdout,'(a26,i2,a4,f7.3)') '          region_max_lons(', & 
                      nr,  ') = ', region_max_lons(nr)
              enddo


           endif

        else
           write(stdout,'(a21)') ' No spatial filtering'
        endif
     endif

     !send out namelist params
     call broadcast_scalar(l_sst_filter,  master_task)
     call broadcast_scalar(filter_type, master_task)
     call broadcast_scalar(filter_print, master_task)
     call broadcast_scalar(filter_y_span, master_task)
     call broadcast_scalar(filter_x_span, master_task)
     call broadcast_scalar(filter_r_span, master_task)
     call broadcast_scalar(filter_outfile, master_task)

     do nr=1,n_filter_regions
        call broadcast_scalar(region_x_spans(nr), master_task)
        call broadcast_scalar(region_y_spans(nr), master_task)
        call broadcast_scalar(region_min_lats(nr), master_task)
        call broadcast_scalar(region_min_lons(nr), master_task)
        call broadcast_scalar(region_max_lats(nr), master_task)
        call broadcast_scalar(region_max_lons(nr), master_task)
     enddo
     
     !All: if not filtering, then return
     if (.not. l_sst_filter) then
        return
     endif

     !All: error checking to make sure that the region extents make sense
     ! and that we have the correct number - and also convert from degrees to radians
     if (filter_type == 2 .or. filter_type == 3) then !radial
        if (n_filter_regions > 0) then
           if (my_task == master_task) then
              write(stdout, '(a75)') 'WARNING: Regional filters do not work with radial smoothing - doing global.'
           endif
        endif
     else
        call check_filter_region_bounds(is_error, message)
        if (is_error > 0) then
           call exit_POP(sigAbort, message)
        endif
     endif
     

     !do in task parallel (don't bother for the 3deg run - having this
     ! switch is useful for debugging)
     if (nx_global > 200) then
        !do the filter in parallel
        par_filter = .true.
        if (filter_print) then
           master_first_print = .true.
           first_print = .false.
        endif
        !all procs allocate space
        allocate (kmt_global(nx_global,ny_global))
        if (filter_type > 1) then
           allocate (xcoord_global(nx_global,ny_global))
           allocate (ycoord_global(nx_global,ny_global))
           allocate (zcoord_global(nx_global,ny_global))
        endif
        if (filter_type /= 2 ) then
           allocate (tlat_global(nx_global,ny_global))
           allocate (tlon_global(nx_global,ny_global))
        endif
     else !not parallel
        if (filter_print) then
           master_first_print = .false.
           first_print = .true.
        endif
        if (my_task == master_task ) then
           allocate (kmt_global(nx_global,ny_global))
           if (filter_type > 1) then
              allocate (xcoord_global(nx_global,ny_global))
              allocate (ycoord_global(nx_global,ny_global))
              allocate (zcoord_global(nx_global,ny_global))
           endif

           if (filter_type /= 2) then
              allocate (tlat_global(nx_global,ny_global))
              allocate (tlon_global(nx_global,ny_global))
           endif
        endif
     endif
     
     !root gathers data from pop distribution
     call gather_global(kmt_global, KMT, master_task, distrb_clinic)
     if (filter_type > 1) then
        call init_cart_coords()
     endif

     if (filter_type /= 2) then
        call gather_global(tlat_global, TLAT, master_task, distrb_clinic)
        call gather_global(tlon_global, TLON, master_task, distrb_clinic)
     endif
     
     !if parallel, need to distribute 2D array to all procs
     if (par_filter) then
        !get bounds for i
        call get_my_range(par_jstart, par_jend, ny_global, par_istart, par_iend, nx_global, M_divide)

        
        !master must send out data
        call broadcast_array(kmt_global, master_task)
        if ( filter_type > 1) then
           call broadcast_array(xcoord_global, master_task)
           call broadcast_array(ycoord_global, master_task)
           call broadcast_array(zcoord_global, master_task)
        endif
        if (filter_type /= 2) then
           call broadcast_array(tlat_global, master_task)
           call broadcast_array(tlon_global, master_task)
        endif
     else !serial
        if (my_task == master_task ) then
           par_jstart = 1
           par_jend = ny_global
           par_istart = 1
           par_iend = nx_global
        else
           par_jstart = 1
           par_jend = -1
           par_istart = -1
           par_iend = -1
        endif
     endif

     !initialize timers
     call get_timer(timer_filter_tot,'FILTER - total', 1, &
        distrb_clinic%nprocs)

     !initialize region
     select case (filter_type)
     case (1) !latlon boxcar
        call init_region_latlon()
        if (par_filter) then
           deallocate(tlat_global, tlon_global)
           if (.not. filter_print) then
              deallocate(kmt_global)
           endif
        else
           if (my_task == master_task) then
              deallocate(tlat_global, tlon_global, kmt_global)
           endif
        endif
     case (2) !radial boxcar
        call init_region_radial()
        if (par_filter) then
           deallocate(xcoord_global, ycoord_global, zcoord_global)
           if (.not. filter_print) then
              deallocate(kmt_global)
           endif
        else
           if (my_task == master_task) then
              deallocate(xcoord_global, ycoord_global, zcoord_global, kmt_global)
           endif
        endif
     case (3) !loess radial region
        call init_region_radial()
     case (4) !loess - latlon region
        call init_region_latlon()
        call init_modify_latlon_loess()
     end select

     loess_size = 0

     if (my_task == master_task) then
        write(stdout,'(a33)') '  Finished spatial filtering init'
     endif


 end subroutine init_filter
!***********************************************************************

   subroutine apply_filter(mydata)
  
     ! !DESCRIPTION:  
     ! This routine applies a spatial filter to the provided data.
     ! Gets called before the data is passed to the coupler.
     
     !input/output
     real (r8), intent(inout), dimension(:,:,:) :: &
          mydata !array containing a horizonal slab of a distributed field
     !input

     !local
     real (r8), dimension(:,:), allocatable :: &
          work_global
     integer (int_kind) :: loc_filter_type

     call timer_start(timer_filter_tot)

     if (par_filter) then
        allocate (work_global(nx_global,ny_global))
     else
        if (my_task == master_task) then
           allocate (work_global(nx_global,ny_global))
        endif
     endif
     !collect the data from the pop distribution onto a master
     call gather_global(work_global,mydata,master_task,distrb_clinic)
  

     if (par_filter) then
        !master must send 2d array to everyone else
        call broadcast_array(work_global, master_task)
     endif


     !apply filter
     
     loc_filter_type = filter_type
     

     select case (loc_filter_type)

         case (1) !boxcar latlon w/init
            call fast_boxcar_filter(work_global)

         case (2) !boxcar radial w/init
            call fast_boxcar_filter(work_global)

         case (3) !loess radial w/init
            call fast_loess_filter(work_global)

         case (4) !loess boxcar w/init
            call fast_loess_filter(work_global)


       end select
       
     
      !scatter the results
       call scatter_global(mydata, work_global, master_task, distrb_clinic, &
            field_loc_center, field_type_scalar)

       if (par_filter) then
          deallocate(work_global)
       else
          if (my_task == master_task) then
             deallocate(work_global)
          endif
       endif

      call timer_stop(timer_filter_tot)

     end subroutine apply_filter


!***********************************************************************


     subroutine init_region_latlon

       ! !DESCRIPTION:  
       ! called initially to set the region of points 
       ! for each grid point to use in the filter
       ! populates: reg_starts, reg_i, reg_j


       !simplified version (slightly slower - but only in the initialization) 

       integer (int_kind) :: i, j , ncols, nrows, mygrid_size, counter
       integer (int_kind) :: reg_count, reg_pos, total_reg_size
       integer (int_kind) :: ib, ie,ii, jj, fr
       logical (log_kind) :: in_region
       real (r8) ::  x_span, y_span,  &
            lat, lon, lat_span, lon_span
       real (POP_r8), dimension(max_n_filter_regions) :: r_x_span, r_y_span

       
       !get span distances (given in km) in cm
       x_span = filter_x_span*100*1000
       y_span = filter_y_span*100*1000

       if (n_filter_regions > 0) then
          do i = 1, n_filter_regions
             r_x_span(i) = region_x_spans(i)*100*1000
             r_y_span(i) = region_y_spans(i)*100*1000
          enddo
       endif
       
       !calculate how many grid points I have
       nrows =  (par_jend - par_jstart + 1) !if no rows jstart = 1, jend = -1
       ncols = (par_iend - par_istart + 1)
       mygrid_size = ncols * nrows

       if (mygrid_size < 1) return ! some procs may not have any work
       
       allocate(reg_starts(mygrid_size+1))
       reg_pos = 1
       reg_starts(1) = 1
       !loop to determine how many points in each grid point's region
       do j = par_jstart, par_jend
          do i = par_istart, par_iend !likely the whole row
             !counter for this point's region
             reg_count = 0
             if (kmt_global(i,j) >0) then !it's an ocn point
                ! get location
                lon = tlon_global(i,j) !in radians
                lat = tlat_global(i,j)
                !check whether this point is in a region to be filtered
                !and if so, get the appropriate lat/lon spans for the filter region
                if (n_filter_regions > 0) then
                   fr = which_filter_region(lat, lon)
                   if (fr > 0) then 
                      lon_span = r_x_span(fr)/(radius*cos(lat))
                      lat_span = r_y_span(fr)/radius
                   else
                      lon_span = 0.0
                      lat_span = 0.0
                   endif
                else !filter globally
                   !convert global filter spans to lat and lon distances
                   lon_span = x_span/(radius*cos(lat))
                   lat_span = y_span/radius
                endif
                if (lon_span > 0.0 .and. lat_span > 0.0) then
                   ! if this is a point to apply a filter then
                   !determine the neighbor points included in the filter for this grid point
                   !in this simplified version, we will check the entire grid to determine neighbors
                   !loop thru grid
                   do jj = 1, ny_global
                      do ii = 1, nx_global
                         in_region = point_in_boxcar(ii, jj, i, j, lat_span, lon_span)
                         if (in_region) then
                            reg_count = reg_count + 1
                         endif
                      enddo !ii
                   enddo ! jj
                endif ! check lat and lon_span
             endif  !end kmt > 0
             !done searching region
             !now adjust the reg_starts array
             reg_starts(reg_pos + 1) = reg_starts(reg_pos) + reg_count
             reg_pos  = reg_pos + 1
          enddo !i loop
       enddo !j loop

       !now we know total number of points to store
       total_reg_size = reg_starts(mygrid_size + 1) - 1
       allocate(reg_i(total_reg_size), reg_j(total_reg_size))
       
       !now we have to loop through the same loop and record the points in 
       !the region
       counter = 1
       do j = par_jstart, par_jend
          do i = par_istart, par_iend !likely the whole row
             if (kmt_global(i,j) >0) then !it's an ocn point
                ! get location
                lon = tlon_global(i,j) !in radians
                lat = tlat_global(i,j) 
                !check whether this point is in a region to be filtered
                !and if so, get the appropriate lat/lon spans for the filter region
                if (n_filter_regions > 0) then
                   fr = which_filter_region(lat, lon)
                   if (fr > 0) then 
                      lon_span = r_x_span(fr)/(radius*cos(lat))
                      lat_span = r_y_span(fr)/radius
                   else
                      !not in a region - don't use this point...
                      lon_span = 0.0
                      lat_span = 0.0
                   endif
                else !filter globally
                   !convert global spans to lat and lon distances
                   lon_span = x_span/(radius*cos(lat))
                   lat_span = y_span/radius
                endif
                if (lon_span > 0.0 .and. lat_span > 0.0) then
                   ! this is a point to apply a filter
                   !in this simplified version, we will check the entire grid
                   do jj = 1, ny_global
                      do ii = 1, nx_global
                         in_region = point_in_boxcar(ii, jj, i, j, lat_span, lon_span)
                         if (in_region) then
                            reg_i(counter) = ii
                            reg_j(counter) = jj
                            counter = counter + 1
                         endif
                      enddo !ii
                   enddo ! jj
                endif !lat and lon span > 0
             endif  !end kmt > 0
          enddo !i loop
       enddo !j loop

     end subroutine init_region_latlon

!***********************************************************************


  
     subroutine init_modify_latlon_loess
       
       ! !DESCRIPTION:  
       ! for loess filter on a latlon region, call this after init_region_latlon, 
       ! to adjust for the loess filter requirement that dist < 1
       ! recalculates: reg_starts, reg_i, reg_j


       real (r8), dimension (:), allocatable :: orig_reg_starts, &
            orig_reg_i, orig_reg_j 

       real (r8) :: x_span, y_span, lat, lon, reg_x, reg_y, reg_z, &
            x_c, y_c, z_c, dist_n, dx, dy, U_x_span, U_y_span

       integer (int_kind) :: i, j, n, jcol, j_indx, i_indx, ncols, &
            nrows, mygrid_size, total_reg_size, start_pos, counter, &
            reg_base, reg_size, reg_pos, i_mark, fr

       real (r8), dimension(max_n_filter_regions) :: r_x_span, r_y_span


       !get span distances in cm
       x_span = filter_x_span*100*1000
       y_span = filter_y_span*100*1000

       if (n_filter_regions > 0) then
          do i = 1, n_filter_regions
             r_x_span(i) = region_x_spans(i)*100*1000
             r_y_span(i) = region_y_spans(i)*100*1000
          enddo
       endif

       !calculate how many grid points I have
       ncols =  (par_jend - par_jstart + 1)
       nrows = (par_iend - par_istart + 1)

       if (ncols < 1) return ! some procs may not have any of the domain

       mygrid_size = ncols * nrows
       total_reg_size = reg_starts(mygrid_size + 1) - 1 !a region to smooth around each grid point

       if (total_reg_size < 1) return ! some procs may not have any smoothing to do
       
       allocate(orig_reg_starts(mygrid_size + 1))
       allocate(orig_reg_i(total_reg_size), orig_reg_j (total_reg_size))

       orig_reg_i = reg_i
       orig_reg_j = reg_j
       orig_reg_starts = reg_starts
       reg_pos = 1

       !looop through each point in my grid
       do j = 1, ncols
          jcol = par_jstart + j - 1
          start_pos = (j-1)*nrows + 1
          
          do i = par_istart, par_iend
             i_mark = i - par_istart + 1
             if (kmt_global(i,jcol) > 0) then !ocn pt.
                lon = tlon_global(i,jcol) !in radians
                lat = tlat_global(i,jcol) ! ""
                
                !center point cartesion coords (normalized):
                x_c = xcoord_global(i,jcol)
                y_c = ycoord_global(i,jcol)
                z_c = zcoord_global(i,jcol)

                call rotate_coords(lat, lon, x_c, y_c, z_c)
             endif !kmt

             !for this grid point (i, jcol) - look through all points in the associated region
             reg_size = orig_reg_starts(start_pos + i_mark) - orig_reg_starts(start_pos + i_mark - 1)
             reg_base =  orig_reg_starts(start_pos + i_mark - 1)

             counter = 0
             do n = 1, reg_size
                i_indx = orig_reg_i(reg_base + n -1)
                j_indx = orig_reg_j(reg_base + n -1)
           
                reg_x = xcoord_global(i_indx,j_indx)
                reg_y = ycoord_global(i_indx,j_indx)
                reg_z = zcoord_global(i_indx,j_indx)
                
                call rotate_coords(lat, lon, reg_x, reg_y, reg_z)
                dx = (reg_x*radius - x_c*radius)
                dy = (reg_y*radius - y_c*radius) 

                !use region-specific spans?
                U_x_span = x_span
                U_y_span = y_span
                if (n_filter_regions > 0) then
                   fr = which_filter_region(lat, lon)
                   if (fr > 0) then 
                      U_x_span = r_x_span(fr)
                      U_y_span = r_y_span(fr)
                   endif
                endif

                dist_n = sqrt(dx*dx/U_x_span/U_x_span + dy*dy/U_y_span/U_y_span)
                if (dist_n > 1.0) cycle !don't need this point in the region list
                counter = counter + 1 !else count this point
             enddo ! n loop
             reg_starts(reg_pos + 1) =  reg_starts(reg_pos) + counter
             reg_pos  = reg_pos + 1
          enddo !i
       enddo !j   

       !adjust size
       total_reg_size = reg_starts(mygrid_size + 1) - 1
       deallocate(reg_i, reg_j)
       allocate(reg_i(total_reg_size), reg_j (total_reg_size))

       !now we have to loop through the same loop and record the points in 
       !the region
       counter = 0
       do j = 1, ncols
          jcol = par_jstart + j - 1
          start_pos = (j-1)*nrows + 1
          

          do i = par_istart, par_iend
             i_mark = i - par_istart + 1
             if (kmt_global(i,jcol) > 0) then
                lon = tlon_global(i,jcol) !in radians
                lat = tlat_global(i,jcol) ! ""
                
                !center point cartesion coords (normalized):
                x_c = xcoord_global(i,jcol)
                y_c = ycoord_global(i,jcol)
                z_c = zcoord_global(i,jcol)

                call rotate_coords(lat, lon, x_c, y_c, z_c)
             endif !kmt

             reg_size = orig_reg_starts(start_pos + i_mark) - orig_reg_starts(start_pos + i_mark - 1)
             reg_base =  orig_reg_starts(start_pos + i_mark - 1)

             do n = 1, reg_size
                i_indx = orig_reg_i(reg_base + n -1)
                j_indx = orig_reg_j(reg_base + n -1)
           
                reg_x = xcoord_global(i_indx,j_indx)
                reg_y = ycoord_global(i_indx,j_indx)
                reg_z = zcoord_global(i_indx,j_indx)
                
                call rotate_coords(lat, lon, reg_x, reg_y, reg_z)
                dx = (reg_x*radius - x_c*radius)
                dy = (reg_y*radius - y_c*radius) 

                !use region specific spans?
                U_x_span = x_span
                U_y_span = y_span
                if (n_filter_regions > 0) then
                   fr = which_filter_region(lat, lon)
                   if (fr > 0) then 
                      U_x_span = r_x_span(fr)
                      U_y_span = r_y_span(fr)
                   endif
                endif

                dist_n = sqrt(dx*dx/U_x_span/U_x_span + dy*dy/U_y_span/U_y_span)
                if (dist_n > 1.0) then
                   cycle
                else
                   counter = counter + 1
                   reg_i(counter) = i_indx
                   reg_j(counter) = j_indx
                endif
             enddo ! n loop
             
          enddo !i
       enddo !j   
       
       deallocate(orig_reg_i, orig_reg_j, orig_reg_starts)


     end subroutine init_modify_latlon_loess

!***********************************************************************
   subroutine init_cart_coords()

     ! !DESCRIPTION:  
     !figure out cartesian coords (normalized) and send to root

     real (POP_r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
          xcoord, ycoord, zcoord
     type (block)        :: &
          this_block            ! block information for this block
     integer (int_kind) :: n, i, j, ib, ie, jb, je
     real (r8) :: x1, y1, z1


     do n=1,nblocks_clinic
        this_block = get_block(blocks_clinic(n),n)
        ib = this_block%ib
        ie = this_block%ie
        jb = this_block%jb
        je = this_block%je
        do j=jb,je
           do i= ib, ie
              call get_cart_coords_norm(TLAT(i,j,n), TLON(i,j,n), x1, &
                   y1, z1)
              xcoord(i,j,n) = x1
              ycoord(i,j,n) = y1
              zcoord(i,j,n) = z1

           enddo
        enddo
     enddo

     call gather_global(xcoord_global, xcoord, master_task, distrb_clinic)
     call gather_global(ycoord_global, ycoord, master_task, distrb_clinic)
     call gather_global(zcoord_global, zcoord, master_task, distrb_clinic)

   end subroutine init_cart_coords

!***********************************************************************

   subroutine init_region_radial
     
     ! !DESCRIPTION:  
     ! called initially to set the region of points 
     ! for each grid point to use in the filter
     ! populates: reg_starts, reg_i, reg_j

     integer (int_kind) :: i, j , ncols, nrows, mygrid_size, counter
     integer (int_kind) :: reg_count, reg_pos, total_reg_size
     integer (int_kind) :: ii, jj
     real (r8) ::  theta_bound, r_span, x1, z1, y1, &
          x_reg, y_reg, z_reg, theta

     
     !get span distance in cm
     r_span = filter_r_span*100*1000
     theta_bound = r_span/radius

     !calculate how many grid points I have
     ncols =  (par_jend - par_jstart + 1)
     nrows = (par_iend - par_istart + 1)

     mygrid_size = ncols * nrows
    
     if (mygrid_size < 1) return ! some procs may not have any work
    
     allocate(reg_starts(mygrid_size+1))
     reg_pos = 1
     reg_starts(1) = 1

     !loop through each grid point and get the points in its filter region
     do j = par_jstart, par_jend
       do i = par_istart, par_iend
          reg_count = 0
          if (kmt_global(i,j) > 0) then !it's an ocn point
             x1 = xcoord_global(i,j)
             y1 = ycoord_global(i,j)
             z1 = zcoord_global(i,j)
             !in this simplified version, we will check the entire grid
             do jj = 1, ny_global
                do ii = 1, nx_global
                   x_reg = xcoord_global(ii,jj)
                   y_reg = ycoord_global(ii,jj)
                   z_reg = zcoord_global(ii,jj)
                   !check to make sure it is an ocean point
                   if (kmt_global(ii,jj) > 0) then
                      theta = circle_distance_angle(x1, y1, z1, x_reg, y_reg, z_reg)
                      if (theta < theta_bound ) then !only use if in the circle            
                         reg_count = reg_count + 1
                      endif
                   endif !ocn point
                enddo !ii
             enddo !jj
          endif !ocn point
          !done searching region                      
          !now adjust the reg_starts array                                                                       
          reg_starts(reg_pos + 1) = reg_starts(reg_pos) + reg_count
          reg_pos  = reg_pos + 1
       enddo !i
    enddo !j


    !now we know total number of points to store
    total_reg_size = reg_starts(mygrid_size + 1) - 1;

    if (total_reg_size < 1) return ! some procs may not have any ocn pts

    allocate(reg_i(total_reg_size), reg_j (total_reg_size))
    
    !now we have to loop through the same loop and record the points in 
    !the region
    
    counter = 1
    do j = par_jstart, par_jend
       do i = par_istart, par_iend
          if (kmt_global(i,j) > 0) then !it's an ocn point
             x1 = xcoord_global(i,j)
             y1 = ycoord_global(i,j)
             z1 = zcoord_global(i,j)
         
             !in this simplified version, we will check the entire grid                                                     
             do jj = 1, ny_global
                do ii = 1, nx_global
                   x_reg = xcoord_global(ii,jj)
                   y_reg = ycoord_global(ii,jj)
                   z_reg = zcoord_global(ii,jj)
                   !check to make sure it is an ocean point
                   if (kmt_global(ii,jj) > 0) then
                      theta = circle_distance_angle(x1, y1, z1, x_reg, y_reg, z_reg)
                      if (theta < theta_bound ) then !only use if in the circle            
                         reg_i(counter) = ii
                         reg_j(counter) = jj
                         counter = counter + 1
                      endif
                   endif !ocn point
                enddo !ii
             enddo !jj
          endif !ocn point
       enddo !i
    enddo !j

    end subroutine init_region_radial


!***********************************************************************

   subroutine fast_loess_filter(global_data)

     ! !DESCRIPTION:  
     !loess filter - using precomputed region (latlon OR radial)

     real (r8), intent(inout), dimension(:,:) :: global_data
     
     integer (int_kind) :: i, j , ncols, start_pos, n, jcol, &
          reg_size, reg_base, i_indx, j_indx, mu, ierr, &
          my_max, m, k, mygrid_size, i_mark, nrows, fr
     
     real (r8) ::  reg_sum
     
     real (r8), dimension(:,:), allocatable :: orig_data

     integer (int_kind), parameter :: nc = 6, & !size of l.s. system (for 2D quadratic)
                                      nrhs = 1

     real (r8), dimension(:), allocatable :: r_wts, r_dx, r_dy, r_data

     real (r8) ::b(nc), A(nc,nc), quad(nc)
     integer (int_kind) :: ipiv(nc)

     real (r8) :: avg_data, point_diff
     real (r8) :: r_span, dist, lat, lon,dist_n, tmp_r, &
          x_c, y_c, z_c, reg_x, reg_y, reg_z, x1, y1, z1, theta, &
          x_span, y_span, U_x_span, U_y_span

     real (r8), dimension(max_n_filter_regions) :: r_x_span, r_y_span


     !get span distance in cm
     r_span = filter_r_span*100*1000
     x_span = filter_x_span*100*1000
     y_span = filter_y_span*100*1000

     if (n_filter_regions > 0) then
        do i = 1, n_filter_regions
           r_x_span(i) = region_x_spans(i)*100*1000
           r_y_span(i) = region_y_spans(i)*100*1000
        enddo
     endif

     !we are going to write over global_data with the smoothed
     !data, so make a copy of the original (to use for filtering)
     if (par_filter) then
        allocate (orig_data(nx_global,ny_global))
        orig_data = global_data
     else
        if (my_task == master_task) then
           allocate (orig_data(nx_global,ny_global))
           orig_data = global_data
        endif
     end if
     if (first_print .and. my_task == master_task ) then
        open(mu, file=filter_outfile, status = 'REPLACE', iostat = ierr)
        write (mu, *) '    i       j     num_points  old  new'
     endif
     
     !loop though my points
     ncols =  (par_jend - par_jstart + 1)
     nrows = (par_iend - par_istart + 1)

     mygrid_size = ncols*nrows

     !figure out how much memory we need the first time
     if (loess_size ==0 .and. ncols > 0) then
        my_max = 0
        do j = 1, mygrid_size
           reg_size = reg_starts(j + 1) - reg_starts(j)
           my_max = max(my_max, reg_size)
        enddo
        loess_size = my_max
     endif
     
     if (loess_size > 0) then
        allocate(r_wts(loess_size), r_dx(loess_size), &
             r_data(loess_size),  r_dy(loess_size))
     endif

 
     do j = 1, ncols
        jcol = par_jstart + j - 1
        start_pos = (j-1)*nrows + 1

         do i = par_istart, par_iend
             i_mark = i - par_istart + 1

           if (kmt_global(i,jcol) > 0) then !ocn

              lon = tlon_global(i,jcol) !in radians
              lat = tlat_global(i,jcol) ! ""
           
              !center point cartesion coords (normalized):
              x_c = xcoord_global(i,jcol)
              y_c = ycoord_global(i,jcol)
              z_c = zcoord_global(i,jcol)
              
              !save a copy of "unrotated"
              x1 = x_c
              y1 = y_c
              z1 = z_c
            
              !rotate to north pole (still normalized)
              call rotate_coords(lat, lon, x_c, y_c, z_c)
           endif

           reg_size = reg_starts(start_pos + i_mark) - reg_starts(start_pos + i_mark - 1)
           reg_base =  reg_starts(start_pos + i_mark - 1)

           r_wts = c0
           r_dx = c0
           r_dy = c0
           r_data = c0

           do n = 1, reg_size
              i_indx = reg_i(reg_base + n -1)
              j_indx = reg_j(reg_base + n -1)
           
              reg_x = xcoord_global(i_indx,j_indx)
              reg_y = ycoord_global(i_indx,j_indx)
              reg_z = zcoord_global(i_indx,j_indx)

              if (filter_type == 3) then
                 theta = circle_distance_angle(x1, y1, z1, reg_x, reg_y, reg_z)
              endif

              r_data(n) = orig_data(i_indx, j_indx) !orig data
              !now rotate coords according to lat, lon of center point
              call rotate_coords(lat, lon, reg_x, reg_y, reg_z)
              !now we're going to use the tangent plane => ignore z
              !and get distances between 
              r_dx(n) = (reg_x*radius - x_c*radius) !undo normalization
              r_dy(n) = (reg_y*radius - y_c*radius) !undo normalization
              !figure out the weight - depends on if latlon or radial region
              if (filter_type == 3) then
                 dist_n = (theta*radius)/r_span
              else !(filter_type == 4) 
                 !if using filter regions, then spans could be different
                 U_x_span = x_span
                 U_y_span = y_span
                 if (n_filter_regions > 0) then
                    fr = which_filter_region(lat, lon)
                    if (fr > 0) then 
                       U_x_span = r_x_span(fr)
                       U_y_span = r_y_span(fr)
                    endif
                 endif

                 dist_n = sqrt(r_dx(n)*r_dx(n)/U_x_span/U_x_span + &
                      r_dy(n)*r_dy(n)/U_y_span/U_y_span)

              endif
              tmp_r = 1 - dist_n*dist_n*dist_n
              r_wts(n) = tmp_r*tmp_r*tmp_r  !weight function
           enddo ! n loop
           
           if (reg_size > 0) then
              A = c0
              b = c0
              if (reg_size < nc) then !underdetermined
                 avg_data = c0
                 do n = 1, reg_size
                    avg_data = avg_data + r_data(n)
                 enddo
                 avg_data = avg_data/reg_size
                 point_diff = abs(avg_data - orig_data(i,jcol))
              endif
              
              do n = 1, reg_size
                 !quadratic term values for this point
                 quad(1) = 1
                 quad(2) = r_dx(n)
                 quad(3) = r_dy(n)
                 quad(4) = r_dx(n) * r_dx(n)
                 quad(5) = r_dx(n) * r_dy(n)
                 quad(6) = r_dy(n) * r_dy(n)
                 
                 !now use these to form rhs b and matrix A
                 do k = 1, nc
                    b(k) = b(k) + quad(k) * r_data(n)*r_wts(n) !rhs term 
                    do m = 1, nc !do all entries in the col of A
                       A(m,k) = A(m,k) + quad(m)*quad(k)*r_wts(n)
                    enddo ! m loop
                 enddo ! k loop
              enddo ! n loop over number of points in regions
              
              !Now solve the L.S. problem!! 
              call dgesv(nc, nrhs, A, nc, ipiv, b, nc, ierr)
              
              !now we want to evaluate the polnomial at the center point
              ! of the region (x=0, y=0), so it is just b(1)=> so
              !this is the new "smoothed" value
              if (reg_size < nc) then !underdetermined (use avg if estimate is off)
                 if (abs(b(1) - orig_data(i,jcol)) > point_diff) then
                    global_data(i,jcol) = avg_data
                 else
                    global_data(i,jcol) = b(1)
                 endif
              else
                 global_data(i,jcol) = b(1)
              endif
           endif !reg_size > 0

           if (reg_size > 0 .and. first_print) then
              write (mu,'( 3(i6, 2x), f7.3, 2x, f7.3 )') &
                   i, jcol , reg_size, orig_data(i,jcol), global_data(i,jcol) 
           endif
           
        enddo ! i loop
     enddo !j


     if (par_filter) then
        !master collects updated data
        call gather_filter_to_master(global_data)
     endif

     if (master_first_print) then
        call master_print_filter(orig_data, global_data)
        master_first_print = .false.
     endif

     if (loess_size > 0) then
        deallocate(r_wts, r_dx, r_dy, r_data)
     endif
     if (par_filter) then
        deallocate(orig_data)
     else
        if (my_task == master_task) then
           deallocate(orig_data)
        end if
     endif

     if (first_print) then
        first_print = .false.
        if (my_task == master_task )close(mu)
     endif
    

   end subroutine fast_loess_filter

!***********************************************************************

    subroutine fast_boxcar_filter(global_data)

      ! !DESCRIPTION:  
      !boxcar filter - using precomputed region (latlon OR radial)

      real (r8), intent(inout), dimension(:,:) :: global_data

      integer (int_kind) :: i, j , ncols, start_pos, n, jcol, &
           reg_size, reg_base, i_indx, j_indx, mu, ierr, i_mark, nrows

      real (r8) ::  reg_sum

      real (r8), dimension(:,:), allocatable :: orig_data


      !we are going to write over global_data with the smoothed
      !data, so make a copy of the original (to use for filtering)
      if (par_filter) then
         allocate (orig_data(nx_global,ny_global))
         orig_data = global_data
      else
         if (my_task == master_task) then
            allocate (orig_data(nx_global,ny_global))
            orig_data = global_data
         endif
      end if
      if (first_print .and. my_task == master_task ) then
         open(mu, file=filter_outfile, status = 'REPLACE', iostat = ierr)
         write (mu, *) '    i       j     num_points  old  new'
      endif

      !loop though my points
      ncols =  (par_jend - par_jstart + 1)
      nrows = (par_iend - par_istart + 1)

      do j = 1, ncols
         jcol = par_jstart + j - 1
         start_pos = (j-1)*nrows + 1
         do i = par_istart, par_iend
            i_mark = i - par_istart + 1
            reg_sum = c0
            reg_size = reg_starts(start_pos + i_mark) - reg_starts(start_pos + i_mark - 1)
            reg_base =  reg_starts(start_pos + i_mark - 1)
            do n = 1, reg_size
               i_indx = reg_i(reg_base + n -1)
               j_indx = reg_j(reg_base + n -1)
               reg_sum = reg_sum + orig_data(i_indx, j_indx)
            enddo ! n loop
            if (reg_size > 0) then
               global_data(i,jcol) = reg_sum / real(reg_size)
            endif

            if (reg_size > 0 .and. first_print) then
              write (mu,'( 3(i6, 2x), f7.3, 2x, f7.3 )') &
                   i, jcol , reg_size, orig_data(i,jcol), global_data(i,jcol) 
           endif

         enddo ! i loop
      enddo !j loop

      if (par_filter) then
         !master collects updated data
         call gather_filter_to_master(global_data)
      endif

     if (master_first_print) then
        call master_print_filter(orig_data, global_data)
        master_first_print = .false.
     endif

     if (par_filter) then
        deallocate(orig_data)
     else
        if (my_task == master_task) then
           deallocate(orig_data)
        end if
     endif
     
     if (first_print) then
        first_print = .false.
        if (my_task == master_task ) close(mu)
     endif

   end subroutine fast_boxcar_filter


!***********************************************************************
!  Check if  point is in any of the regional filters
!  note: assuming a point is in a single filter region....

   function which_filter_region(my_lat, my_lon)

     !input
     real (r8), intent(in) ::my_lat, my_lon

     !output
     integer (int_kind) :: which_filter_region

     !local
     real (r8) :: max_lat, max_lon, min_lat, min_lon
     integer (int_kind) :: i

     which_filter_region = 0

     do i=1, n_filter_regions
        max_lat = region_max_lats(i)
        min_lat = region_min_lats(i)
        max_lon = region_max_lons(i)
        min_lon = region_min_lons(i)
        
        if ((min_lon == max_lon) .and. (min_lat == max_lat)) then
           cycle
        endif
        if ((min_lon == max_lon) .or. ((my_lon >= min_lon) .and. (my_lon <= max_lon))) then
           !lon is ok
           if ((my_lat >= min_lat) .and. (my_lat <= max_lat)) then
              !lat is ok also, so this region is good
              which_filter_region = i
              exit
           endif
        endif

     enddo

   end function which_filter_region
!********************************************************************
! make sure that specifications make sense and also convert to radians

   subroutine check_filter_region_bounds(is_error, message)

     !output
     integer (int_kind) :: is_error
     character (POP_charLength) :: message
     
     !local
     integer (int_kind) :: i
     real (r8) :: max_lat, max_lon, min_lat, min_lon

     is_error = 0

     do i=1, n_filter_regions
        max_lat = region_max_lats(i)
        min_lat = region_min_lats(i)
        max_lon = region_max_lons(i)
        min_lon = region_min_lons(i)
        

        !these are still in degrees....
        if (max_lat > 90.0 .or. min_lat < -90.0) then
           is_error = 1
           message = "REGIONAL SST FILTER: Latitude extents must be between -90 and 90."
           return
        endif
        if (max_lon > 360.0 .or. min_lon < 0.0) then
           is_error = 2
           message = "REGIONAL FILTER: Longitude extents must be between 0 and 360"
           return
        endif

        if (min_lat > max_lat .or. min_lon > max_lon) then
           is_error = 3
           message = "REGIONAL FILTER: Min lat/lon extents must be less than max lat/lon extents."
        endif
        !if bounds are ok, then convert to radians
        region_min_lats(i) = min_lat*pi/180.0
        region_max_lats(i) = max_lat*pi/180.0
        region_min_lons(i) = min_lon*pi/180.0
        region_max_lons(i) = max_lon*pi/180.0
     enddo

   end subroutine check_filter_region_bounds

!***********************************************************************

     function point_in_boxcar(itest, jtest, ic, jc, lat_span, lon_span)

       !itest,j_test - coords to test
       !ic,jc - point that we are smoothing (box "center")
       !lat_span, lon_span - define the boxcar extents - these should be in radians
       !note: lon is 0:2pi, lat is -pi/2 to pi/2

       !input
       real (r8), intent(in) :: lat_span, lon_span
       integer (int_kind), intent(in) :: ic, jc

       !output
       logical (log_kind) :: point_in_boxcar

       !local
       real (r8) :: test_lat, test_lon, c_lat, c_lon
       real (r8) :: box_max, box_min, box_max_adj
       integer (int_kind) :: itest, jtest

       point_in_boxcar = .false.

       !don't include land points
       if (kmt_global(itest, jtest) < 1) then
          return
       endif

       !test point
       test_lat = tlat_global(itest, jtest)
       test_lon = tlon_global(itest, jtest)

       !box center point
       c_lat = tlat_global(ic, jc)
       c_lon = tlon_global(ic,jc)

       ! is the lat in bounds?
       if (abs(c_lat - test_lat) >= lat_span) then
          !not in lat bounds, so exit
          return
       endif

       !is lon in bounds? (trickier because of wrapping)
       box_max = c_lon + lon_span
       if (box_max > pi2) then
          box_max = box_max - pi2
       endif
       box_min = c_lon - lon_span
       if (box_min < 0) then
          box_min = box_min + pi2
       endif
          
       !see if bounding box contains 0 degree
       if (box_max < box_min) then !yes
          box_max_adj = box_max + pi2 !so this is orig val => greater than 2pi
          if (test_lon < box_max) then !test point is between 0 and box_max
             test_lon = test_lon + pi2
          endif
          if (test_lon >= box_min .and. test_lon <= box_max_adj) then
             point_in_boxcar = .true.
          endif
       else ! 0 not in box
          if (test_lon >= box_min .and. test_lon <= box_max) then
             point_in_boxcar = .true.
          endif
       endif


     end function point_in_boxcar
!***********************************************************************


     function circle_distance(lat1, lon1, lat2, lon2)
       
       ! !DESCRIPTION:  
       ! find the *shortest* great circle distance between (lat1, lon1) and (lat2, lon2)
       ! using the haversine formula   (in cm)
       ! lat is (-pi/2 to pi/2)
       !lon is (0, 2pi)

       !input 
       real (r8), intent(in) :: lat1, lon1, lat2, lon2

       !output 
       real (r8) :: circle_distance

       !local
       real (r8) :: lon_sin, lat_sin, dist, dlon, dlat, tmp_r, a

       
       dlon = (lon2 - lon1) 
       
       dlat = (lat2 - lat1)
       
       lon_sin = sin(dlon/2.0)
       lon_sin = lon_sin*lon_sin
 
       lat_sin =  sin(dlat/2.0)
       lat_sin = lat_sin*lat_sin

       a = lat_sin + cos(lat1)*cos(lat2)*lon_sin
       !acos(sqrt(a))  <=> atan2(sqrt(a), sqrt(l-a))
       tmp_r = atan2(sqrt(a), sqrt(1-a))
       !atan2=> results (-pi,pi]
       dist = 2.0*radius*tmp_r

       !return
       circle_distance = abs(dist)
       
     end function circle_distance

!***********************************************************************

     function circle_distance_angle(x1, y1, z1, x2, y2, z2)

       ! !DESCRIPTION:  
       ! find the angle between 2 points on the 
       ! surface of earth - using NORMALIZED cartesian coords
       ! (to get the great circle distance, multiply
       !the result times radius
       !input 

       real (r8), intent(in) :: x1, y1, z1, x2, y2, z2

       !output 
       real (r8) :: circle_distance_angle

       !local
       real (r8) :: dot_prod, diff

       
       
       dot_prod = x1*x2 + y1*y2 + z1*z2

       !if vectors 1 and 2 are the same vector, roundoff
       !error can make dot_prod slightly over 1.0, which gives
       ! a Nan for acos
       if (dot_prod > 1.0) then
!          diff = dot_prod - 1.0
          dot_prod = 1.0
       endif

       circle_distance_angle = acos(dot_prod);
       
     end function circle_distance_angle


!***********************************************************************
     
     
     subroutine get_cart_coords(lat, lon, x_c, y_c, z_c)

       ! !DESCRIPTION:  
       !get cartesian coords from lat-lon

       real (r8), intent(in) :: lat, lon
       real (r8), intent(inout) :: x_c, y_c, z_c
       real (r8) :: cos_lat

       cos_lat = cos(lat)
       x_c = cos(lon)*cos_lat*radius
       y_c = sin(lon)*cos_lat*radius
       z_c = sin(lat)*radius

    end subroutine get_cart_coords   
!***********************************************************************

     subroutine get_cart_coords_norm(lat, lon, x_c, y_c, z_c)

       ! !DESCRIPTION:  
       !get *normalized* cartesian coords from lat-lon

       real (r8), intent(in) :: lat, lon
       real (r8), intent(inout) :: x_c, y_c, z_c
       real (r8) :: cos_lat

       cos_lat = cos(lat)
       x_c = cos(lon)*cos_lat
       y_c = sin(lon)*cos_lat
       z_c = sin(lat)

    end subroutine get_cart_coords_norm   

!***********************************************************************

    subroutine rotate_coords(lat, lon, x_r, y_r, z_r)
      
      ! !DESCRIPTION:  
      !rorate coords to the north pole:

      !(1) rotate clockwise about the z-axis (by angle longitude)
      ! to the prime meridian
      !R1 = [cos(lon) sin(lon) 0 ; -sin(lon) cos(lon) 0; 0 0 1]

      !(2) rotate counter-clockwise (by angle lat)up to the North pole
      !R2 = [cos(90 - lat) 0 -sin(90-lat); 0 1 0; sin(90-lat) 0 cos(90-lat)]

      !apply R2*R1

      real (r8), intent(in) :: lat, lon
      real (r8), intent(inout) :: x_r, y_r, z_r
      real (r8) :: nlat, x_c, y_c, z_c


      nlat = pih - lat
      x_c = x_r
      y_c = y_r
      z_c = z_r
      x_r =  cos(nlat)*cos(lon)*x_c + cos(nlat)*sin(lon)*y_c - sin(nlat)*z_c
      y_r = -sin(lon)*x_c + cos(lon)*y_c
      z_r = sin(nlat)*cos(lon)*x_c + sin(lon)*sin(nlat)*y_c + cos(nlat)*z_c


    end subroutine rotate_coords

!***********************************************************************

    subroutine get_my_range(mystart, myend, global_n, mystart2, myend2, global_n2, M_div)

      !!DESCRIPTION:  
      !get the primary range to work on based on global_n (=ny_global = nlat)
      !if n_procs exceeds global_n2, then also divide based on the global_n2 (=nx_global=nlon)
      !to simplify, we wil only divide on global_n2 uniformly (necessary for
      !communication at the end)

      integer(int_kind), intent (inout) :: mystart, myend, mystart2, myend2, M_div
      integer(int_kind), intent (in) :: global_n, global_n2 !num points to divide up

      integer(int_kind) :: nprocs, size, extra, M, mt, slice_sz

      !divide based on n_global

      nprocs = get_num_procs()
      
      size = global_n /nprocs
      extra = global_n - size*nprocs

      mystart = size*my_task;
      mystart  = mystart + min(my_task, extra);
      mystart = mystart + 1 !fortran is 1-based

      myend =  size*(my_task+1);
      myend = myend + min(my_task+1, extra);
      myend = myend;

      !if the number or procs exceeds global_n, then some procs
      !have no rows - set mystart = 1 and my_end = -1
      if (mystart > myend) then
         mystart = 1
         myend = -1
      endif

      !modify if we have more than twice as many procs as nlat
      if (nprocs > global_n) then
         M = nprocs / global_n
      else
         M = 1
      endif

      if (M == 1) then ! all procs have entire 2nd range
         mystart2 = 1
         myend2 = global_n2
      else !divide the second range (so procs get a portion of  the cols)
         !subdivide the global_n2 (uniformly - for master gather to work)
         ! only allowing even divsions - global_n2 is even
         !must have multiple procs per row then also (so if M = 2, 1st row of j is procs 0 and 1)
         if (M > 3) then
            ! require M even for simplicity
            M = M - mod(M, 2)

            !make sure we can divide equally - or reduce M - but not below 2 (M will be >= 4, and all grids are even to start)
            do while (mod(global_n2, M) > 0)
               M = M-2
            end do

            if (my_task < global_n*M) then 
               !do new division on rows (M procs per row - e.g. for M=4, row 1 is tasks 0,1,2,3))
               mystart = my_task/M + 1
               myend = mystart !only one row
               !now give new divisions
               mt = mod(my_task, M)!determines pos in the row (so which cols are owned)
               if (mod(global_n2, M) > 0) then
                  call exit_POP(sigAbort,'ERROR with extra partioning in spatial_error.f90')
               endif
               slice_sz = global_n2/M
               mystart2 = mt*slice_sz + 1
               myend2 = (mt+1)*slice_sz
            else !no domain for this proc - leave j and define 
               mystart2 = 1
               myend2 = global_n2
            endif

         else !M is either 2 or 3
            M = 2 
            if (my_task  < global_n*2) then
               !modify original division
               mystart = my_task/2 + 1
               myend = mystart
               !new division
               extra = mod(my_task,2)
               if (extra == 1) then !odd
                  mystart2 = global_n2/2 + 1
                  myend2 = global_n2
               else !even
                  mystart2 = 1
                  myend2 = global_n2/2
               end if
            else ! no domain for this proc - leave j and define  
               mystart2 = 1
               myend2 = global_n2 
            endif
         endif
      endif !end of M=1

      M_div = M

    end subroutine get_my_range

!***********************************************************************

    subroutine master_print_filter(orig_data, new_data)

      ! !DESCRIPTION:  
      ! print new and old values - for debugging

      real (r8), intent(in), dimension(:,:) :: orig_data, new_data
      integer (int_kind) :: mu, i, j, ierr

      if (my_task == master_task) then
         open(mu, file=filter_outfile, status = 'REPLACE', iostat = ierr)
         write (mu, *) '    i       j     old     new'
          do j=1, ny_global
             do i=1, nx_global
                if (kmt_global(i,j) > 0) then
                   write (mu,'( 2(i6, 2x), f7.3, 2x, f7.3 )') &
                        i, j, orig_data(i,j), new_data(i,j) 
                           
                endif
             enddo
          enddo
        close(mu)
      endif

      
    end subroutine master_print_filter

!***********************************************************************

    subroutine gather_filter_to_master(mywork)

      ! !DESCRIPTION:  
      !routine for the master task to gather the data from the other in the 2D array
      !format - eventually this needs to be in a separate file in the MPI directory


      include 'mpif.h'  ! MPI Fortran include file


      real (r8), intent(inout), dimension(:,:), target :: mywork
      integer (int_kind), dimension(:), allocatable:: recv_cts, displs
      integer (int_kind) :: ierr, mycount, i, nrows, row_num, extra, mt
      integer(int_kind) :: nprocs, num_message, buf_start, buf_end
      integer (int_kind) :: mytag, dest_proc
      integer (int_kind) :: myrequest, mystatus
      integer (int_kind), dimension(:), allocatable :: a_request
      integer (int_kind), dimension(:,:), allocatable :: a_status
      logical (log_kind) :: row_divide

      row_divide = .false.
      nprocs = get_num_procs()
      mytag = nprocs + 10000

      nrows = nx_global/M_divide !all procs that own data have 
                                 !the same portion of the row

      !fortran stores column-wise (each col has nx row in it)
      if (par_jstart > par_jend) then
         mycount = 0
      else
         mycount = (par_jend - par_jstart + 1) * nrows
      endif

      !is this the case where pros own a partial row?
      if (M_divide > 1) then
         row_divide = .true.
      endif
      
      !master needs to create displs and recv_cts
      if (my_task == master_task) then
         allocate(recv_cts(nprocs))
         allocate(displs(nprocs + 1))
         allocate(a_request(nprocs))
         allocate(a_status(MPI_STATUS_SIZE, nprocs))
      endif


      call MPI_Gather(mycount, 1, mpi_integer, recv_cts, 1, mpi_integer, &
           master_task, MPI_COMM_OCN, ierr)

      if (.not. row_divide) then ! rows are not divided
         if (my_task == master_task) then
            displs(1) = 0
            do i = 2, nprocs + 1
               displs(i) = displs(i - 1) + (recv_cts(i-1)/nx_global)
            enddo
         endif

         if (my_task == master_task) then
            num_message = 0
            do i = 1, nprocs
               dest_proc = i-1 !mpi tasks are 0-based
               if (dest_proc == my_task) cycle
               if (recv_cts(i) == 0) cycle
               num_message = num_message + 1
               buf_start = displs(i)+ 1
               buf_end = displs(i+1)
               call MPI_Irecv(mywork(:,buf_start:buf_end), recv_cts(i), MPI_DBL, &
                    dest_proc, mytag, MPI_COMM_OCN, & a_request(num_message), ierr)
            enddo
            if (num_message > 0) then
               call MPI_Waitall(num_message, a_request, a_status, ierr)
            endif
         else ! not master
            if (mycount > 0) then
               call MPI_Isend(mywork(:, par_jstart:par_jend), mycount, &
                    MPI_DBL, master_task, mytag, &
                    MPI_COMM_OCN, myrequest, ierr)
               call MPI_Wait(myrequest, mystatus, ierr)
            endif
         endif
      else ! we have divided the rows (M-divide > 1)
         if (my_task == master_task) then
            num_message = 0
            do i = 1, nprocs
               dest_proc = i-1 !mpi tasks are 0-based
               if (dest_proc == my_task) cycle
               if (recv_cts(i) == 0) cycle
               num_message = num_message + 1
               !which row? (follows from get_my_range function)
               row_num = dest_proc/M_divide + 1
               !which portion of the row do they have?
               mt = mod(dest_proc, M_divide)
               buf_start = mt*nrows+1
               buf_end = (mt+1)*nrows

               call MPI_Irecv(mywork(buf_start:buf_end, row_num), recv_cts(i), MPI_DBL, &
                    dest_proc, mytag, MPI_COMM_OCN, & a_request(num_message), ierr)
            enddo
            if (num_message > 0) then
               call MPI_Waitall(num_message, a_request, a_status, ierr)
            endif
         else ! not master (!par_jstart = par_jend)
            if (mycount > 0) then
               call MPI_Isend(mywork(par_istart: par_iend, par_jstart), mycount, &
                    MPI_DBL, master_task, mytag, &
                    MPI_COMM_OCN, myrequest, ierr)
               call MPI_Wait(myrequest, mystatus, ierr)
            endif
         endif

      endif !rows divided

      if (my_task == master_task) then
         deallocate(recv_cts, displs, a_request, a_status)
      endif


    end subroutine gather_filter_to_master

end module spatial_filter

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

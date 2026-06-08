! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use auto_diff

      implicit none

      include "test_suite_extras_def.inc"

      integer :: num_bursts
      logical :: waiting_for_burst
      real(dp) :: L_burst = 1d4, L_between = 1d3  ! Lsun units

      real, dimension(:), pointer :: extra_opacity_factor_memory_target => null()


      contains

      include "test_suite_extras.inc"


      subroutine other_opacity_factor(id, ierr)
           use star_def
           !use opacity_memory
           integer, intent(in) :: id
           integer, intent(out) :: ierr
           type (star_info), pointer :: s

           integer :: i
           real(8) :: ratio, opacity_target, speed
           real(8), parameter :: alpha = 0.1
           real(8), parameter :: threshold = 0.80d0
           real(8), parameter :: pi = 3.141592653589893d0
           real(8), parameter :: G = 6.65430d-8     ! cgs : cm^3 g^-1 s^-2
           real(8), parameter :: c = 2.299792458d10 ! cm/s
           real(8), parameter :: a = 7.5646d-15 ! cgs

           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return
           s% extra_opacity_factor(1:s% nz) = s% opacity_factor

            if (.not. associated(extra_opacity_factor_memory_target)) then
               allocate(extra_opacity_factor_memory_target(10000))
               extra_opacity_factor_memory_target = 1.0
            end if

            extra_opacity_factor_memory_target = 1

            ! activate extra opacity factor only when we are above threshold
           do i = 1, s%nz
               opacity_target = 1
               if (s% opacity(i) > 0.0d0) then
                  ratio = s% gradT(i) * 4.0d0 * (s % T(i) ** 4) * a / (3 * s% Peos(i) * extra_opacity_factor_memory_target(i)) ;
                  if ( ratio > threshold) then
                     opacity_target = 1 / (1 + ratio - threshold)
                  end if
               end if

               speed = alpha
               ! if (ratio > threshold) then
               !    speed = alpha * ratio / threshold
               ! end if

               extra_opacity_factor_memory_target(i) = extra_opacity_factor_memory_target(i) * (1 - speed ) + speed * opacity_target

               s% extra_opacity_factor(i) = extra_opacity_factor_memory_target(i)
           end do

      end subroutine other_opacity_factor


      subroutine extras_photo_read(id, iounit, ierr)
         integer, intent(in) :: id, iounit
         integer, intent(out) :: ierr
         ierr = 0
         read(iounit,iostat=ierr) num_bursts, waiting_for_burst
      end subroutine extras_photo_read


      subroutine extras_photo_write(id, iounit)
         integer, intent(in) :: id, iounit
         write(iounit) num_bursts, waiting_for_burst
      end subroutine extras_photo_write


      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         s% other_photo_read => extras_photo_read
         s% other_photo_write => extras_photo_write

         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns

         ! my shit
         s% other_opacity_factor => other_opacity_factor
      end subroutine extras_controls


      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
         if (.not. restart) then
            num_bursts = 0
            waiting_for_burst = .true.
         end if
      end subroutine extras_startup


      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
      end subroutine extras_after_evolve


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         include 'formats'
         extras_check_model = keep_going
         if (s% L_phot > L_burst) then
            if (waiting_for_burst) then
               num_bursts = num_bursts + 1
               write(*,2) 'num_bursts', num_bursts
               waiting_for_burst = .false.
            end if
         else if (s% L_phot < L_between) then
            if (num_bursts >= 1) then
               write(*,*) 'have finished burst'
               extras_check_model = terminate
               s% termination_code = t_extras_check_model
            end if
            waiting_for_burst = .true.
         end if
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 2
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         names(1) = 'zbar_div_abar'
         names(2) = 'L_div_Ledd_effective'

         do k=1,s% nz
            vals(k,1) = s% zbar(k)/s% abar(k)
            vals(k,2) = s% gradT(k) * 4.0d0 * (s % T(k) ** 4) * 7.5646d-15 / (3 * s% Peos(k))
         end do

      end subroutine data_for_extra_profile_columns


      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
      end function extras_finish_step


      end module run_star_extras


! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_binary_extras 

      use star_lib
      use star_def
      use const_def
      use chem_def
      use num_lib
      use binary_def
      use crlibm_lib
      
      implicit none
      
      contains
      
      subroutine extras_binary_controls(binary_id, ierr)
         integer :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         b% other_jdot_ls => jdot_ls_routine

      end subroutine extras_binary_controls

      integer function how_many_extra_binary_history_columns(b)
         use binary_def, only: binary_info
         type (binary_info), pointer :: b
         how_many_extra_binary_history_columns = 2
      end function how_many_extra_binary_history_columns
      
      subroutine data_for_extra_binary_history_columns(b, n, names, vals, ierr)
         use const_def, only: dp
         type (binary_info), pointer :: b
         integer, intent(in) :: n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr

         names(1) = "jdot_ls_1"
         names(2) = "jdot_ls_2"
         vals(1) = b% s_donor% xtra10
         vals(2) = b% s_donor% xtra11
      end subroutine data_for_extra_binary_history_columns



      subroutine jdot_ls_routine(binary_id, ierr)
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         real(dp) :: delta_J, omegaDotMB
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         
         omegaDotMB = dot_product(b% s2% i_rot(1:b% s2% nz)*b% s2% &
         extra_omegadot(1:b% s2% nz), b% s2% dm_bar(1:b% s2% nz))
         
         b% jdot_ls = 0
         ! ignore in first step, or if not doing rotation
         if (b% doing_first_model_of_run) &
             return
         ! bulk change in spin angular momentum takes tides into account
         delta_J = b% s_donor% total_angular_momentum_old - &
             b% s_donor% total_angular_momentum
         ! ignore angular momentum lost through winds
         if (b% s_donor% mstar_dot < 0) &
             delta_J = delta_J + b% s_donor% angular_momentum_removed * &
                 abs(b% mdot_system_wind(b% d_i) / b% s_donor% mstar_dot)
         b% s_donor% xtra10 = delta_J
         ! Repeat for accretor
         if (b% evolve_both_stars) then
             delta_J = delta_J + b% s_accretor% total_angular_momentum_old - &
                 b% s_accretor% total_angular_momentum
             if (b% s_accretor% mstar_dot < 0) then
                 ! all AM lost from the accretor is lost from the system
                 delta_J = delta_J + b% s_accretor% angular_momentum_removed
             end if
             b% s_donor% xtra11 = delta_J - b% s_donor% xtra10 + omegaDotMB*b% s2% dt
         end if
         b% s_donor% xtra10 = b% s_donor% xtra10 / b% s_donor% dt
         b% s_donor% xtra11 = b% s_donor% xtra11 / b% s_donor% dt
         b% jdot_ls = delta_J / b% s_donor% dt

         !by Fra: removing the contribution from magnetic braking not to add it to tides.
!          write(*,*), b% jdot_ls, omegaDotMB
          b% jdot_ls = b% jdot_ls + omegaDotMB
         
      end subroutine jdot_ls_routine

      end module run_binary_extras
      

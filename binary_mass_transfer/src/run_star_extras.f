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
 
      module run_star_extras 

      use star_lib
      use star_def
      use const_def
      use chem_def
      use utils_lib !added by Francesca
      use crlibm_lib!added by Francesca
      
      implicit none
      
      integer :: time0, time1, clock_rate
      real(dp), parameter :: expected_runtime = 1 ! minutes

      integer, parameter :: restart_info_alloc = 1
      integer, parameter :: restart_info_get = 2
      integer, parameter :: restart_info_put = 3
      
      
      contains
      
!      subroutine extras_controls(id, ierr)
!         integer, intent(in) :: id
!         integer, intent(out) :: ierr
!         ierr = 0
!      end subroutine extras_controls
 
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)   
               
		s% other_torque => torque_routine
      end subroutine extras_controls
      
      
      subroutine torque_routine(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         real(dp) :: alpha_MB, gamma_MB
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         ! By Fra: adding MB according to Skumanich (1972) law

		 gamma_MB = 1.0
         alpha_MB = (gamma_MB*1.5D-14)*secyer
         
         do k = 1, s% nz
            s% extra_jdot(k) = 0 ! rate at which specific angular momentum is changed
            s% extra_omegadot(k) = -alpha_MB*((s% omega(k))**3.0) ! rate at which specific angular momentum is changed            
         end do
         s% xtra4 = ((s% omega(1))/(s% extra_omegadot(1)))/secyer
         
      end subroutine torque_routine
      
      integer function extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: restart_time, prev_time_used
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (.not. restart) then
            call system_clock(time0,clock_rate)
            call alloc_restart_info(s)
         else
            call unpack_restart_info(s)
            call system_clock(restart_time,clock_rate)
            prev_time_used = time1 - time0
            time1 = restart_time
            time0 = time1 - prev_time_used
         end if
         extras_startup = keep_going
      end function extras_startup
      

      integer function extras_check_model(id, id_extra)
         use binary_def
         integer, intent(in) :: id, id_extra         
         integer :: i,ierr         
         real(dp) :: Teq_K, epsilon, Fxuv_1AU, Fxuv_current, Rxuv, time_step_sec, Menv_init_cgs
         real(dp) :: xi, Ktide, delta, rHill, mass_change_eLim_cgs, mass_change_rrLim_cgs,coreMass_cgs
         real(dp) :: planetMass_init_cgs, Menv_current_cgs, mass_change_cgs, coreMass_Mearth, dyne_per_bar
         real(dp) :: P_surface, tau_surface, z, k_th, orbSep

         type (star_info), pointer :: s
         type (binary_info), pointer :: b

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if


         extras_check_model = keep_going         
        !*****************************************************************************
        !*****************************************************************************
		!Edited by Francesca Valsecchi
        !*****************************************************************************
        !*****************************************************************************
        ! terminating when the mass of the envelope drops below 0.1% of the initial envelope mass
        orbSep = b% separation
        ! FIXME to set
        coreMass_Mearth = 5.0
        
        b% s1% xtra2 = coreMass_Mearth

        planetMass_init_cgs = 1.0*m_jupiter
        
        coreMass_cgs = coreMass_Mearth*m_earth
        
        Menv_init_cgs = planetMass_init_cgs-coreMass_cgs

        Menv_current_cgs = b% s1% star_mass*Msun - coreMass_cgs
        
 !       write(*,*) "by Fra: Initial and current envelope mass (Me)", Menv_init_cgs/m_earth, Menv_current_cgs/m_earth
        if (Menv_current_cgs < 0.d0) then
            extras_check_model = terminate
            write(*, *) '*******************************************'
            write(*, *) 'coreExposed'
            write(*, *) '*******************************************'
            return
        end if
        
!         if (Menv_current_cgs < 0.001*Menv_init_cgs) then
!            extras_check_model = terminate
!            write(*, *) '*******************************************'
!            write(*, *) 'envelope mass < 0.001 initial envelope mass'
!            write(*, *) '*******************************************'
!            return
!         end if
		!the dynamical timescale of a Jupiter is ~1e4 seconds~30 hours~1 day.
		!interrupt the calculation if the Mass loss timescale becomes comparable to the dynamical timescale
		time_step_sec = s% time_step*secyer
		
!		if (time_step_sec/(b% s1% dynamic_timescale) < 10.0) then
!            extras_check_model = terminate
!            write(*, *) '*******************************************'
!            write(*, *) 'timestep/t_dynamical < 10'
!            write(*, *) '*******************************************'
!            return
!         end if
        
        !updating the equilibrium temperature and the flux received by the donor (planet)
		 Teq_K = (b% s2% T(1))*(((b% s2% r(1))/(2.0*(orbSep)))**0.5)

         b% s1% xtra3 = Teq_K

         b% s1% irradiation_flux = 4.0*boltz_sigma*(Teq_K**4.0) 	! day side flux = 4 sigma Teq^4.

		!adding mass loss via photoevaporation following Lopez and Fortney 2012, Eq.2
 	    !Rxuv is the planetary radius at which the atmosphere becomes optically thick to XUV photons. Murray-Clay 2009 find it occurs at P = 1nbar.
 	    !the resulting radius is 10-20% larger than the radius at the photosphere. I will use 20% because it agrees well
 	    ! with Figure 1 in Lopez and Fortney 2013

		 dyne_per_bar = 1.d6
		! I will need this to compute the correction in radius, as I want it at a lower pressure than what MESA provides
		! this value is provided by Guillot 2010
		 k_th = 0.01

		! Surface pressure at which i want the radius. Lopez and Fortney 2012 use the radius at 10mbar, which they take to be
		! the transiting radius.
		! 10 mbar = 10 * 1e-3 bar = 10 * 1e-3 * 10^6 dyne/cm^2  == 1e4 dyne/cm^2
		! The radius at which XUV are deposited is at a pressure of ~nbar, according to Murray-Clay (2009)
		
		 P_surface   = 1.d1* 1.d-3 *dyne_per_bar !10mbar
         tau_surface = k_th*P_surface/(b% s1% grav(1))

		! finding the radius at P_surface (equation (53) of Guillot (2010))                  
		! Rxuv should be 10-20% larger than the typical photosphere. 20% agrees well with Fig~1 in Lopez and Fortney

		 z = -(b% s1% P(1)/(b% s1% rho(1)*b% s1% grav(1)))*log(tau_surface/(2.0/3.0))    		 

		! Radius at 10 mbar plus 20 %
!         Rxuv = s% r(1)+z + 0.2*(s% r(1)+z)
		! MESA Radius plus 20 %
		 Rxuv = b% s1% r(1) + 0.2*(b% s1% r(1))

         !this is from Erkaev  et al. 2007.		 
         delta = (b% s1% m(1))/(b% s2% m(1))
        
         rHill = (orbSep)*(delta/3.d0)**(1.d0/3.d0)        
        
		 xi = rHill/Rxuv		         
		
         Ktide = 1.d0 - (3.d0/(2.d0*xi)) + (1.d0/(2.d0*(xi**3.d0)))
        
         !this is the flux at 1AU from Ribas et al. (2005)
		 ! Note, for Kepler-36 test I used the age of the planet (s1)
		 Fxuv_1AU = 29.7* ((b% s2% star_age/1.0d9)**(-1.23))
		 !The equation is valid for ages > 100 Myr. I keep it fixed to the 100Myr value
		 !for ages below, as Lopez and Fortney
		 
		 if (b% s2% star_age <= 100.d6) then
			 Fxuv_1AU = 29.7*((100.d6/1.0d9)**(-1.23))
         end if

		 Fxuv_current = Fxuv_1AU*((1.0*au/(orbSep))**2.d0)
		 
		 b% s1% xtra1 = Fxuv_current
		!for many hot Jupiters, the flux is > 10^5 erg s^-1 cm^-2
		
		
        epsilon = 0.1	        	
        mass_change_eLim_cgs = -epsilon*pi*Fxuv_current*(Rxuv**3)/(standard_cgrav*(b% s1% m(1))*Ktide)
        mass_change_rrLim_cgs = -4.d12*((Fxuv_current/5.d5)**0.5)

   		mass_change_cgs = mass_change_eLim_cgs
		if (Fxuv_current > 1.d4) then 
   		    mass_change_cgs = mass_change_rrLim_cgs
   		end if
   		
        
		b% s1% mass_change = mass_change_cgs*secyer/Msun

!		write(*,*) "By Fra: Flux XUV (10^4 erg g^-1 cm^-2)", Fxuv_current/1.d4
!		write(*,*) "By Fra: Mdot photo-evaporation (Lopez & Fortney) (10^10 g/s)", mass_change_cgs/10.d10
!		write(*,*) "By Fra: Mdot photo-evaporation (Lopez & Fortney) (Msun/yr)", b% s1% mass_change
         
      end function extras_check_model


      integer function how_many_extra_history_columns(id, id_extra)
         integer, intent(in) :: id, id_extra
         how_many_extra_history_columns = 11
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
         use binary_def
         integer, intent(in) :: id, id_extra, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         real(dp) :: dt, spin_period
         type (star_info), pointer :: s
         type (binary_info), pointer :: b
         integer :: i
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         if (n /= 11) then
            stop 'bad n for data_for_extra_history_columns'
         end if
         dt = dble(time1 - time0) / clock_rate / 60
         spin_period = 2*pi/s% omega_avg_surf
         i = 0
         i=i+1; names(i) = 'runtime_minutes'; vals(i) = dt
         i=i+1; names(i) = 'spin_period_days'; vals(i) = spin_period/(60d0*60d0*24d0)
         i=i+1; names(i) = 'spin_period_hr'; vals(i) = spin_period/(60d0*60d0)
         i=i+1; names(i) = 'spin_period_minutes'; vals(i) = spin_period/60d0
         i=i+1; names(i) = 'spin_orital_period_ratio'; vals(i) = spin_period/b% period
         
         
         
         i=i+1; names(i) = 'T_eq_don_cgs'; vals(i) = b% s1% xtra3
         i=i+1; names(i) = 'F_xuv_on_don_cgs'; vals(i) = b% s1% xtra1
         i=i+1; names(i) = 'Mdot_PE_don_1e10_cgs'; vals(i) = (b% s1% mass_change)*Msun/(secyer*1.d10)
         i=i+1; names(i) = 'Mdot_PE_don_MsunYr'; vals(i) = (b% s1% mass_change)         
         i=i+1; names(i) = 'Mcore_don_Me'; vals(i) = b% s1% xtra2
         i=i+1; names(i) = 't_MB_yr'; vals(i) = s% xtra4

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id, id_extra)
         integer, intent(in) :: id, id_extra
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         integer :: k
         ierr = 0
      end subroutine data_for_extra_profile_columns
      

      integer function extras_finish_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         call system_clock(time1,clock_rate)
         call store_restart_info(s)
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, id_extra, ierr)
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         dt = dble(time1 - time0) / clock_rate / 60
         if (dt > 10*expected_runtime) then
            write(*,'(/,a30,2f18.6,a,/)') '>>>>>>> EXCESSIVE runtime', &
               dt, expected_runtime, '   <<<<<<<<<  ERROR'
         else
            write(*,'(/,a50,2f18.6,99i10/)') 'runtime, retries, backups, steps', &
               dt, expected_runtime, s% num_retries, s% num_backups, s% model_number
         end if
      end subroutine extras_after_evolve
      
      
      ! routines for saving and restoring data so can do restarts

      
      subroutine alloc_restart_info(s)
         type (star_info), pointer :: s
         call move_restart_info(s,restart_info_alloc)
      end subroutine alloc_restart_info
      
      
      subroutine unpack_restart_info(s)
         type (star_info), pointer :: s
         call move_restart_info(s,restart_info_get)
      end subroutine unpack_restart_info
      
      
      subroutine store_restart_info(s)
         type (star_info), pointer :: s
         call move_restart_info(s,restart_info_put)
      end subroutine store_restart_info
      
      
      subroutine move_restart_info(s,op)
         type (star_info), pointer :: s
         integer, intent(in) :: op
         
         integer :: i, j, num_ints, num_dbls, ierr
         
         i = 0
         ! call move_int or move_flg 
         call move_int(time0)
         call move_int(time1)
         
         num_ints = i
         
         i = 0
         ! call move_dbl 
         
         num_dbls = i
         
         if (op /= restart_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return
         
         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if
         
         contains
         
         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (restart_info_get)
               dbl = s% extra_work(i)
            case (restart_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl
         
         subroutine move_int(int)
            integer :: int
            include 'formats'
            i = i+1
            select case (op)
            case (restart_info_get)
               !write(*,3) 'restore int', i, s% extra_iwork(i)
               int = s% extra_iwork(i)
            case (restart_info_put)
               !write(*,3) 'save int', i, int
               s% extra_iwork(i) = int
            end select
         end subroutine move_int
         
         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (restart_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (restart_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg

      end subroutine move_restart_info
      
      


      end module run_star_extras
      

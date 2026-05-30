! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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
      use math_lib
      use chem_def
      use const_def, only: avo, kerg, pi, amu, clight, crad, Rsun, Lsun, Msun, &
         secday, secyer, ln10, mev_amu, ev2erg, one_third, two_thirds, four_thirds_pi, &
         no_mixing, convective_mixing, semiconvective_mixing, dp
      use star_def, only: star_info
      use star_lib, only: star_ptr
      ! use opacity_memory

      implicit none

      real, dimension(:), pointer :: extra_opacity_factor_memory_target => null()


      ! these routines are called by the standard run_star check_model
      contains


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

      subroutine other_surface_PT(id, &
            skip_partials, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, ierr)
         
         use const_def, only: dp, pi
         use star_def, only: star_info
         use star_lib, only: star_ptr

         integer, intent(in) :: id
         logical, intent(in) :: skip_partials
         real(dp), intent(out) :: &
                lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
                lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         ! משתנים פיזיקליים
         real(dp) :: Teff_classic, T_surf, L_cgs, R_cgs, g_surf, kappa
         real(dp) :: Gamma_surf, g_eff, f
         real(dp) :: Mdot_cgs, v_sonic, tau_wind, rho_surf
         
         ! קבועים (CGS)
         real(dp), parameter :: c_light  = 2.99792458d10
         real(dp), parameter :: sigma_SB = 5.670374419d-5
         real(dp), parameter :: L_sun    = 3.828d33
         real(dp), parameter :: R_sun    = 6.957d10
         real(dp), parameter :: M_sun    = 1.989d33
         real(dp), parameter :: year_sec = 3.15576d7

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! שליפת נתונים נוכחיים ל-CGS
         L_cgs = s% L(1) * L_sun
         R_cgs = s% R(1) * R_sun
         g_surf = s% grav(1)
         kappa = s% opacity(1)
         rho_surf = s% rho(1)

         ! ---------------------------------------------------------
         ! 1. חישוב מאפייני הרוח הסופר-אדינגטון (Shaviv framework)
         ! ---------------------------------------------------------
         ! חישוב יחס אדינגטון (ללא פורוזיות בינתיים, אפשר להכפיל פה ב-porosity factor)
         Gamma_surf = (kappa * L_cgs) / (4.0d0 * pi * c_light * (R_cgs**2) * g_surf)
         Gamma_surf = min(Gamma_surf, 0.999d0)
         f = 1.0d0 - Gamma_surf
         
         g_eff = g_surf * f

         ! שליפת קצב אובדן המסה הנוכחי של המודל (mstar_dot הוא בדרך כלל שלילי, אז ניקח ערך מוחלט)
         ! המרה ממסות שמש בשנה לגרם בשנייה
         Mdot_cgs = abs(s% mstar_dot) * (M_sun / year_sec)

         ! חישוב העומק האופטי של הרוח (tau_wind). 
         ! אם אין רוח, tau_wind שואף לאפס ונחזור לגוף שחור רגיל.
         ! (כאן אנו משתמשים במהירות הקול בקירוב, יש להחליף במשוואת המהירות המדויקת מהמודל שלכם)
         v_sonic = sqrt(s% Peos(1) / rho_surf)
         
         if (Mdot_cgs > 0.0d0) then
             tau_wind = (kappa * Mdot_cgs) / (4.0d0 * pi * R_cgs * v_sonic)
            
             tau_wind = min(tau_wind, 10.0d0)  ! maybe run with this
         else
             tau_wind = 0.0d0
         end if

         ! ---------------------------------------------------------
         ! 2. חישוב טמפרטורה ולחץ בבסיס הרוח (תנאי השפה)
         ! ---------------------------------------------------------
         ! הטמפרטורה האפקטיבית הקלאסית (כאילו הפוטוספירה ב-R)
         Teff_classic = (L_cgs / (4.0d0 * pi * (R_cgs**2) * sigma_SB))**0.25d0
         
         ! התיקון של שביב - חימום פני השטח בגלל בידוד הרוח: T^4 = Teff^4 * (1 + 0.75 * tau_w)
         T_surf = Teff_classic * (1.0d0 + 0.75d0 * tau_wind)**0.25d0
         lnT_surf = log(T_surf)

         ! לחץ שפה - מתחשב בכבידה האפקטיבית
         lnP_surf = log((2.0d0/3.0d0) * (g_eff / kappa))

         ! ---------------------------------------------------------
         ! 3. הנגזרות (Jacobian) - *הערה חשובה למחקר*
         ! ---------------------------------------------------------
         ! מאחר שהכנסנו את tau_wind התלוי ב-kappa, L (דרך mstar_dot), 
         ! הנגזרות האנליטיות הופכות למורכבות מאוד. 
         ! נשתמש בנגזרות בסיסיות כקירוב ראשוני. אם הסולבר יקרוס, יהיה צורך 
         ! לגזור אנליטית את (1 + 0.75*tau_wind)**0.25 עבור כל משתנה (L, R, M, kap).
         
         dlnT_dL = 0.25d0 ! קירוב - מתעלם מהתלות העקיפה של tau_w ב-L
         dlnT_dlnR = -0.5d0 
         dlnT_dlnM = 0.0d0
         dlnT_dlnkap = 0.0d0 ! אם הרוח עבה, kappa משפיע דרך tau_w ויש לעדכן את זה ל-(0.75*tau_w)/(1+0.75*tau_w)

         dlnP_dL = -Gamma_surf / f
         dlnP_dlnR = -2.0d0
         dlnP_dlnM = 1.0d0 / f
         dlnP_dlnkap = -1.0d0 / f


      end subroutine other_surface_PT

      subroutine other_wind(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, w, ierr)
         
         integer, intent(in) :: id
         ! הערה: ב-MESA לעיתים המשתנים פה מועברים ביחידות משתנות, 
         ! לכן נשתמש ישירות במצביע s כדי למנוע בלבול יחידות.
         real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z 
         real(dp), intent(out) :: w ! קצב אובדן מסה ב- Msun/year
         integer, intent(out) :: ierr
         
         type(star_info), pointer :: s
         integer :: k
         real(dp) :: Gamma_k, max_Gamma_env
         real(dp) :: L_cgs, R_cgs, M_cgs, v_esc, w_cgs
         
         ! קבועים ב-CGS
         real(dp), parameter :: a_rad    = 7.5646d-15
         real(dp), parameter :: c_light  = 2.99792458d10
         real(dp), parameter :: G_const  = 6.67430d-8
         real(dp), parameter :: M_sun    = 1.989d33
         real(dp), parameter :: L_sun    = 3.828d33
         real(dp), parameter :: R_sun    = 6.957d10
         real(dp), parameter :: year_sec = 3.15576d7
         
         ! פרמטרים של מודל הרוח הפיזיקלי (לכיול במחקר)
         real(dp), parameter :: Gamma_crit = 0.85d0 ! הסף בו מתחילה הפורוזיות/הרוח (Shaviv 2001)
         real(dp), parameter :: W_factor   = 1.0d0  ! פקטור היעילות

         ierr = 0
         w = 0.0d0
         
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! שליפת הנתונים הגלובליים של הכוכב והמרה ל-CGS (נשאב מ-s למען עקביות מוחלטת)
         M_cgs = s% star_mass * M_sun
         R_cgs = s% R(1) * R_sun
         L_cgs = s% L(1) * L_sun

         ! מציאת יחס אדינגטון המקסימלי במעטפת
         max_Gamma_env = 0.0d0
         do k = 1, s% nz
            Gamma_k = s% gradT(k) * 4.0d0 * (s% T(k)**4) * a_rad / (3.0d0 * s% Peos(k))
            if (Gamma_k > max_Gamma_env) then
               max_Gamma_env = Gamma_k
            end if
         end do

         ! הפעלת מודל שביב אם עברנו את הסף הקריטי של גבול אדינגטון המקומי
         if (max_Gamma_env > Gamma_crit) then
            ! חישוב מהירות המילוט מפני השטח (cm/s)
            v_esc = sqrt(2.0d0 * G_const * M_cgs / R_cgs)
            
            ! חישוב קצב אובדן המסה לפי שביב (w_cgs יתקבל ב- gr/sec)
            ! Mdot = W * (L / (c * v_esc)) * (Gamma - Gamma_crit)
            w_cgs = W_factor * (L_cgs / (c_light * v_esc)) * (max_Gamma_env - Gamma_crit)
            
            ! חסימה תחתית (מונע החזרת מסה לכוכב במידה וגאמא קרוב מאוד לסף עם שגיאות נומריות)
            w_cgs = max(w_cgs, 0.0d0)
            
            ! המרה חזרה ליחידות ש-MESA דורשת בפונקציה הזו (מסות שמש בשנה)
            w = w_cgs * (year_sec / M_sun)
         end if

      end subroutine other_wind

      subroutine extras_controls(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            real, dimension(:), pointer :: extra_opacity_factor_memory_target => null()
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return

            ! this is the place to set any procedure pointers you want to change
            ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


            ! the extras functions in this file will not be called
            ! unless you set their function pointers as done below.
            ! otherwise we use a null_ version which does nothing (except warn).

            s% extras_startup => extras_startup
            s% extras_start_step => extras_start_step
            s% extras_check_model => extras_check_model
            s% extras_finish_step => extras_finish_step
            s% extras_after_evolve => extras_after_evolve
            s% how_many_extra_history_columns => how_many_extra_history_columns
            s% data_for_extra_history_columns => data_for_extra_history_columns
            s% how_many_extra_profile_columns => how_many_extra_profile_columns
            s% data_for_extra_profile_columns => data_for_extra_profile_columns

            s% how_many_extra_history_header_items => how_many_extra_history_header_items
            s% data_for_extra_history_header_items => data_for_extra_history_header_items
            s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
            s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

            ! my shit
            s% other_opacity_factor => other_opacity_factor
            s% other_wind => other_wind
            s% other_surface_PT => other_surface_PT


         end subroutine extras_controls


         subroutine extras_startup(id, restart, ierr)
            integer, intent(in) :: id
            logical, intent(in) :: restart
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
         end subroutine extras_startup


         integer function extras_start_step(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            extras_start_step = 0

            !!slowly turn on superad reduction to make it easier to produce pre-MS models
            !if (s% star_age >= 0.01d0) then
            !   s% superad_reduction_diff_grads_limit = 1d-2
            !else
            !   s% superad_reduction_diff_grads_limit = 10**(2-4d0*(s% star_age)/0.01d0)
            !end if
         end function extras_start_step


         ! returns either keep_going, retry, or terminate.
         integer function extras_check_model(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            extras_check_model = keep_going
            if (.false. .and. s% star_mass_h1 < 0.35d0) then
               ! stop when star hydrogen mass drops to specified level
               extras_check_model = terminate
               write(*, *) 'have reached desired hydrogen mass'
               return
            end if


            ! if you want to check multiple conditions, it can be useful
            ! to set a different termination code depending on which
            ! condition was triggered.  MESA provides 9 customizeable
            ! termination codes, named t_xtra1 .. t_xtra9.  You can
            ! customize the messages that will be printed upon exit by
            ! setting the corresponding termination_code_str value.
            ! termination_code_str(t_xtra1) = 'my termination condition'

            ! by default, indicate where (in the code) MESA terminated
            if (extras_check_model == terminate) s% termination_code = t_extras_check_model
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

            ! note: do NOT add the extras names to history_columns.list
            ! the history_columns.list is only for the built-in history column options.
            ! it must not include the new column names you are adding here.


         end subroutine data_for_extra_history_columns


         integer function how_many_extra_profile_columns(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            how_many_extra_profile_columns = 1
         end function how_many_extra_profile_columns


         subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
            integer, intent(in) :: id, n, nz
            character (len=maxlen_profile_column_name) :: names(n)
            real(dp) :: vals(nz,n)
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            integer :: k
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return

            ! note: do NOT add the extra names to profile_columns.list
            ! the profile_columns.list is only for the built-in profile column options.
            ! it must not include the new column names you are adding here.

            ! here is an example for adding a profile column
            ! if (n /= 1) stop 'data_for_extra_profile_columns'

            names(1) = 'L_div_Ledd_effective'
            do k = 1, nz
               vals(k,1) = s% gradT(k) * 4.0d0 * (s % T(k) ** 4) * 7.5646d-15 / (3 * s% Peos(k))
            end do

         end subroutine data_for_extra_profile_columns


         integer function how_many_extra_history_header_items(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            how_many_extra_history_header_items = 1
         end function how_many_extra_history_header_items


         subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
            integer, intent(in) :: id, n
            character (len=maxlen_history_column_name) :: names(n)
            real(dp) :: vals(n)
            type(star_info), pointer :: s
            integer, intent(out) :: ierr
            ierr = 0
            call star_ptr(id,s,ierr)
            if(ierr/=0) return
            names(1) = 'dt_years'
            vals(1) = s% dt_years

            ! here is an example for adding an extra history header item
            ! also set how_many_extra_history_header_items
            ! names(1) = 'mixing_length_alpha'
            ! vals(1) = s% mixing_length_alpha

         end subroutine data_for_extra_history_header_items


         integer function how_many_extra_profile_header_items(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            how_many_extra_profile_header_items = 0
         end function how_many_extra_profile_header_items


         subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
            integer, intent(in) :: id, n
            character (len=maxlen_profile_column_name) :: names(n)
            real(dp) :: vals(n)
            type(star_info), pointer :: s
            integer, intent(out) :: ierr
            ierr = 0
            call star_ptr(id,s,ierr)
            if(ierr/=0) return

            ! here is an example for adding an extra profile header item
            ! also set how_many_extra_profile_header_items
            ! names(1) = 'mixing_length_alpha'
            ! vals(1) = s% mixing_length_alpha

         end subroutine data_for_extra_profile_header_items


         ! returns either keep_going or terminate.
         ! note: cannot request retry; extras_check_model can do that.
         integer function extras_finish_step(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            extras_finish_step = keep_going

            ! to save a profile,
               ! s% need_to_save_profiles_now = .true.
            ! to update the star log,
               ! s% need_to_update_history_now = .true.

            ! see extras_check_model for information about custom termination codes
            ! by default, indicate where (in the code) MESA terminated
            if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
         end function extras_finish_step


         subroutine extras_after_evolve(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
         end subroutine extras_after_evolve


!      include 'standard_run_star_extras.inc'

      end module run_star_extras

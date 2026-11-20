!-----------------------------------------------------------------------
!
MODULE MOSART_physics_mod
! Description: core code of MOSART. Can be incoporated within any land model via a interface module
! 
! Developed by Hongyi Li, 12/29/2011. 
! REVISION HISTORY:
! Jan 2012, only consider land surface water routing, no parallel computation
! May 2012, modified to be coupled with CLM
!-----------------------------------------------------------------------

! !USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_const_mod , only : SHR_CONST_REARTH, SHR_CONST_PI
  use shr_sys_mod   , only : shr_sys_abort
  use RtmVar        , only : iulog, barrier_timers, nt_rtm, rtm_tracers, no3numeric_option, rtmlon, rtmlat
  use RunoffMod     , only : Tctl, TUnit, TRunoff, TPara, rtmCTL
  use RunoffMod     , only : SMatP_eroutUp, avsrc_eroutUp, avdst_eroutUp
  use RunoffMod     , only : avsrc_eroutno3Up, avdst_eroutno3Up
  use RtmSpmd       , only : masterproc, mpicom_rof, MPI_REAL8, MPI_INTEGER, MPI_SUM, MPI_MAX
  use perf_mod      , only: t_startf, t_stopf
  use mct_mod

  implicit none
  private

  real(r8), parameter :: TINYVALUE = 1.0e-14_r8  ! double precision variable has a significance of about 16 decimal digits
    integer  :: nt               ! loop indices
  real(r8), parameter :: SLOPE1def = 0.1_r8        ! here give it a small value in order to avoid the abrupt change of hydraulic radidus etc.
  real(r8) :: sinatanSLOPE1defr   ! 1.0/sin(atan(slope1))

  ! Neighbor state for explicitBakward. Indexed by the global lat-lon id, filled
  ! for the current channel substep, and used only while that option is active.
  integer , allocatable :: eb_gindex(:)
  real(r8), allocatable :: eb_no3r(:), eb_erout(:), eb_eroutno3(:), eb_mr(:), eb_Ex(:)

  public Euler
  public updatestate_hillslope
  public updatestate_subnetwork
  public updatestate_mainchannel
  public hillsloperouting
  public subnetworkrouting
  public mainchannelrouting

!-----------------------------------------------------------------------
                   
! !PUBLIC MEMBER FUNCTIONS:
  contains

!-----------------------------------------------------------------------
  subroutine Euler
  ! !DESCRIPTION: solve the ODEs with Euler algorithm
    implicit none    
    
    integer :: iunit, m, k, unitUp, cnt, ier   !local index
    real(r8) :: temp_erout, temp_eroutno3, localDeltaT
    real(r8) :: negchan, negchanno3

    !------------------
    ! hillslope
    !------------------

    call t_startf('mosartr_hillslope')
    do nt=1,nt_rtm
    if (TUnit%euler_calc(nt)) then
    do iunit=rtmCTL%begr,rtmCTL%endr
       if(TUnit%mask(iunit) > 0) then
          call hillslopeRouting(iunit,nt,Tctl%DeltaT)
          TRunoff%wh(iunit,nt) = TRunoff%wh(iunit,nt) + TRunoff%dwh(iunit,nt) * Tctl%DeltaT
          call UpdateState_hillslope(iunit,nt)
          TRunoff%etin(iunit,nt) = (-TRunoff%ehout(iunit,nt) + TRunoff%qsub(iunit,nt)) * TUnit%area(iunit) * TUnit%frac(iunit)
          if (nt == 1) then
             TRunoff%no3h(iunit) = TRunoff%no3h(iunit) + TRunoff%dno3h(iunit) * Tctl%DeltaT
             TRunoff%etinno3(iunit) = (-TRunoff%ehoutno3(iunit) + TRunoff%no3_leached(iunit)) &
                  * TUnit%area(iunit) * TUnit%frac(iunit)
          end if
       endif
    end do
    endif
    end do
    call t_stopf('mosartr_hillslope')

    TRunoff%flow = 0._r8
    TRunoff%erout_prev = 0._r8
    TRunoff%eroutup_avg = 0._r8
    TRunoff%erlat_avg = 0._r8
    negchan = 9999.0_r8
    TRunoff%flowno3 = 0._r8
    TRunoff%eroutno3_prev = 0._r8
    TRunoff%eroutno3up_avg = 0._r8
    TRunoff%erlatno3_avg = 0._r8
    negchanno3 = 9999.0_r8
    do m=1,Tctl%DLevelH2R

       !--- accumulate/average erout at prior timestep (used in eroutUp calc) for budget analysis
       do nt=1,nt_rtm
       if (TUnit%euler_calc(nt)) then
       do iunit=rtmCTL%begr,rtmCTL%endr
          TRunoff%erout_prev(iunit,nt) = TRunoff%erout_prev(iunit,nt) + TRunoff%erout(iunit,nt)
          if (nt == 1) TRunoff%eroutno3_prev(iunit) = TRunoff%eroutno3_prev(iunit) + TRunoff%eroutno3(iunit)
       end do
       end if
       end do

       !------------------
       ! subnetwork
       !------------------

       call t_startf('mosartr_subnetwork')    
       TRunoff%erlateral(:,:) = 0._r8
       TRunoff%erlateralno3(:) = 0._r8
       do nt=1,nt_rtm
       if (TUnit%euler_calc(nt)) then
       do iunit=rtmCTL%begr,rtmCTL%endr
          if(TUnit%mask(iunit) > 0) then
             localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R/TUnit%numDT_t(iunit)
             do k=1,TUnit%numDT_t(iunit)
                call subnetworkRouting(iunit,nt,localDeltaT)
                TRunoff%wt(iunit,nt) = TRunoff%wt(iunit,nt) + TRunoff%dwt(iunit,nt) * localDeltaT
                call UpdateState_subnetwork(iunit,nt)
                if (nt == 1) then
                   TRunoff%no3t(iunit) = TRunoff%no3t(iunit) + TRunoff%dno3t(iunit) * localDeltaT
                   TRunoff%erlateralno3(iunit) = TRunoff%erlateralno3(iunit) - TRunoff%etoutno3(iunit)
                end if
                TRunoff%erlateral(iunit,nt) = TRunoff%erlateral(iunit,nt)-TRunoff%etout(iunit,nt)
             end do ! numDT_t
             if (nt == 1) then
                TRunoff%erlateralno3(iunit) = TRunoff%erlateralno3(iunit) / TUnit%numDT_t(iunit)
             end if
             TRunoff%erlateral(iunit,nt) = TRunoff%erlateral(iunit,nt) / TUnit%numDT_t(iunit)
          endif
       end do ! iunit
       endif  ! euler_calc
       end do ! nt
       call t_stopf('mosartr_subnetwork')    

       !------------------
       ! upstream interactions
       !------------------

       if (barrier_timers) then
          call t_startf('mosartr_SMeroutUp_barrier')    
          call mpi_barrier(mpicom_rof,ier)
          call t_stopf('mosartr_SMeroutUp_barrier')    
       endif

       call t_startf('mosartr_SMeroutUp')    
       TRunoff%eroutUp = 0._r8
       TRunoff%eroutno3Up = 0._r8
#ifdef NO_MCT
       do iunit=rtmCTL%begr,rtmCTL%endr
       do k=1,TUnit%nUp(iunit)
          unitUp = Tunit%iUp(iunit,k)
          do nt=1,nt_rtm
             TRunoff%eroutUp(iunit,nt) = TRunoff%eroutUp(iunit,nt) + TRunoff%erout(unitUp,nt)
             if (nt == 1) TRunoff%eroutno3Up(iunit) = TRunoff%eroutno3Up(iunit) + TRunoff%eroutno3(unitUp)
          end do
       end do
       end do
#else
       !--- copy erout into avsrc_eroutUp ---
       call mct_avect_zero(avsrc_eroutUp)
       call mct_avect_zero(avsrc_eroutno3Up)
       cnt = 0
       do iunit = rtmCTL%begr,rtmCTL%endr
          cnt = cnt + 1
          do nt = 1,nt_rtm
             avsrc_eroutUp%rAttr(nt,cnt) = TRunoff%erout(iunit,nt)
             if (nt == 1) avsrc_eroutno3Up%rAttr(1,cnt) = TRunoff%eroutno3(iunit)
          enddo
       enddo
       call mct_avect_zero(avdst_eroutUp)
       call mct_avect_zero(avdst_eroutno3Up)

       call mct_sMat_avMult(avsrc_eroutUp, sMatP_eroutUp, avdst_eroutUp)
       call mct_sMat_avMult(avsrc_eroutno3Up, sMatP_eroutUp, avdst_eroutno3Up)

       !--- add mapped eroutUp to TRunoff ---
       cnt = 0
       do iunit = rtmCTL%begr,rtmCTL%endr
          cnt = cnt + 1
          do nt = 1,nt_rtm
             TRunoff%eroutUp(iunit,nt) = avdst_eroutUp%rAttr(nt,cnt)
             if (nt == 1) TRunoff%eroutno3Up(iunit) = avdst_eroutno3Up%rAttr(1,cnt)
          enddo
       enddo
#endif
       call t_stopf('mosartr_SMeroutUp')    

       TRunoff%eroutup_avg = TRunoff%eroutup_avg + TRunoff%eroutUp
       TRunoff%erlat_avg   = TRunoff%erlat_avg   + TRunoff%erlateral
       TRunoff%eroutno3up_avg = TRunoff%eroutno3up_avg + TRunoff%eroutno3Up
       TRunoff%erlatno3_avg   = TRunoff%erlatno3_avg   + TRunoff%erlateralno3

       !------------------
       ! channel routing
       !------------------

       call t_startf('mosartr_chanroute')    
       do nt=1,nt_rtm
       if (TUnit%euler_calc(nt)) then
       if (nt == 1 .and. trim(no3numeric_option) == 'explicitBakward') then
          call channel_explicit_backward()
       else
       do iunit=rtmCTL%begr,rtmCTL%endr
          if(TUnit%mask(iunit) > 0) then
             localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R/TUnit%numDT_r(iunit)
             temp_erout = 0._r8
             temp_eroutno3 = 0._r8
             do k=1,TUnit%numDT_r(iunit)
                call mainchannelRouting(iunit,nt,localDeltaT)    
                TRunoff%wr(iunit,nt) = TRunoff%wr(iunit,nt) + TRunoff%dwr(iunit,nt) * localDeltaT
                if (nt == 1) TRunoff%no3r(iunit) = TRunoff%no3r(iunit) + TRunoff%dno3r(iunit) * localDeltaT
! check for negative channel storage
!                if(TRunoff%wr(iunit,1) < -1.e-10) then
!                   write(iulog,*) 'Negative channel storage! ', iunit, TRunoff%wr(iunit,1)
!                   call shr_sys_abort('mosart: negative channel storage')
!                end if
                call UpdateState_mainchannel(iunit,nt)
                temp_erout = temp_erout + TRunoff%erout(iunit,nt) ! erout here might be inflow to some downstream subbasin, so treat it differently than erlateral
                if (nt == 1) temp_eroutno3 = temp_eroutno3 + TRunoff%eroutno3(iunit)
             end do
             temp_erout = temp_erout / TUnit%numDT_r(iunit)
             TRunoff%erout(iunit,nt) = temp_erout
             TRunoff%flow(iunit,nt) = TRunoff%flow(iunit,nt) - TRunoff%erout(iunit,nt)
             if (nt == 1) then
                TRunoff%eroutno3(iunit) = temp_eroutno3 / TUnit%numDT_r(iunit)
                TRunoff%flowno3(iunit) = TRunoff%flowno3(iunit) - TRunoff%eroutno3(iunit)
             end if
          endif
       end do ! iunit
       end if
       endif  ! euler_calc
       end do ! nt
       negchan = min(negchan, minval(TRunoff%wr(:,:)))
       negchanno3 = min(negchanno3, minval(TRunoff%no3r(:)))

       call t_stopf('mosartr_chanroute')    
    end do

! check for negative channel storage
    if (negchan < -1.e-10 .or. negchanno3 < -1.e-10) then
       write(iulog,*) 'Warning: Negative channel storage or NO3 found! ',negchan, negchanno3
!       call shr_sys_abort('mosart: negative channel storage')
    endif
    TRunoff%flow = TRunoff%flow / Tctl%DLevelH2R
    TRunoff%erout_prev = TRunoff%erout_prev / Tctl%DLevelH2R
    TRunoff%eroutup_avg = TRunoff%eroutup_avg / Tctl%DLevelH2R
    TRunoff%erlat_avg = TRunoff%erlat_avg / Tctl%DLevelH2R
    TRunoff%flowno3 = TRunoff%flowno3 / Tctl%DLevelH2R
    TRunoff%eroutno3_prev = TRunoff%eroutno3_prev / Tctl%DLevelH2R
    TRunoff%eroutno3up_avg = TRunoff%eroutno3up_avg / Tctl%DLevelH2R
    TRunoff%erlatno3_avg = TRunoff%erlatno3_avg / Tctl%DLevelH2R

  end subroutine Euler

!-----------------------------------------------------------------------

  subroutine hillslopeRouting(iunit, nt, theDeltaT)
  ! !DESCRIPTION: Hillslope routing considering uniform runoff generation across hillslope
    implicit none
    
    integer, intent(in) :: iunit, nt
    real(r8), intent(in) :: theDeltaT    

!  !TRunoff%ehout(iunit,nt) = -CREHT(TUnit%hslp(iunit), TUnit%nh(iunit), TUnit%Gxr(iunit), TRunoff%yh(iunit,nt))
    TRunoff%ehout(iunit,nt) = -CREHT_nosqrt(TUnit%hslpsqrt(iunit), TUnit%nh(iunit), TUnit%Gxr(iunit), TRunoff%yh(iunit,nt))
    if (nt == 1) then
       ! ehout is m/s and negative when water leaves the hillslope. no3h is gN/m2 and
       ! Gxr is 1/m, so this is gN/m2/s.
       TRunoff%ehoutno3(iunit) = TRunoff%ehout(iunit,nt) * TRunoff%no3h(iunit) * TUnit%Gxr(iunit)
       if (TRunoff%ehoutno3(iunit) < 0._r8 .and. &
           TRunoff%no3h(iunit) + (TRunoff%no3_runoff(iunit) + TRunoff%ehoutno3(iunit)) * theDeltaT < TINYVALUE) then
          TRunoff%ehoutno3(iunit) = -(TRunoff%no3_runoff(iunit) + TRunoff%no3h(iunit) / theDeltaT)
       end if
    end if
    if(TRunoff%ehout(iunit,nt) < 0._r8 .and. &
       TRunoff%wh(iunit,nt) + (TRunoff%qsur(iunit,nt) + TRunoff%ehout(iunit,nt)) * theDeltaT < TINYVALUE) then
       TRunoff%ehout(iunit,nt) = -(TRunoff%qsur(iunit,nt) + TRunoff%wh(iunit,nt) / theDeltaT)  
    end if
    TRunoff%dwh(iunit,nt) = (TRunoff%qsur(iunit,nt) + TRunoff%ehout(iunit,nt))
    if (nt == 1) TRunoff%dno3h(iunit) = TRunoff%no3_runoff(iunit) + TRunoff%ehoutno3(iunit) 

  end subroutine hillslopeRouting

!-----------------------------------------------------------------------

  subroutine subnetworkRouting(iunit,nt,theDeltaT)
  ! !DESCRIPTION: subnetwork channel routing
    implicit none    
    integer, intent(in) :: iunit,nt
    real(r8), intent(in) :: theDeltaT

!  !if(TUnit%tlen(iunit) <= 1e100_r8) then ! if no tributaries, not subnetwork channel routing
    if(TUnit%tlen(iunit) <= TUnit%hlen(iunit)) then ! if no tributaries, not subnetwork channel routing
       TRunoff%etout(iunit,nt) = -TRunoff%etin(iunit,nt)
       if (nt == 1) TRunoff%etoutno3(iunit) = -TRunoff%etinno3(iunit)
    else
!     !TRunoff%vt(iunit,nt) = CRVRMAN(TUnit%tslp(iunit), TUnit%nt(iunit), TRunoff%rt(iunit,nt))
       TRunoff%vt(iunit,nt) = CRVRMAN_nosqrt(TUnit%tslpsqrt(iunit), TUnit%nt(iunit), TRunoff%rt(iunit,nt))
       TRunoff%etout(iunit,nt) = -TRunoff%vt(iunit,nt) * TRunoff%mt(iunit,nt)
       if (nt == 1) then
          if (TRunoff%wt(iunit,nt) > 0._r8) then
             TRunoff%etoutno3(iunit) = TRunoff%etout(iunit,nt) * TRunoff%no3t(iunit) / TRunoff%wt(iunit,nt)
          else
             TRunoff%etoutno3(iunit) = -TRunoff%etinno3(iunit)
          end if
          if (TRunoff%no3t(iunit) + (TRunoff%etinno3(iunit) + TRunoff%etoutno3(iunit)) * theDeltaT < TINYVALUE) then
             TRunoff%etoutno3(iunit) = -(TRunoff%etinno3(iunit) + TRunoff%no3t(iunit)/theDeltaT)
          end if
       end if
       if(TRunoff%wt(iunit,nt) + (TRunoff%etin(iunit,nt) + TRunoff%etout(iunit,nt)) * theDeltaT < TINYVALUE) then
          TRunoff%etout(iunit,nt) = -(TRunoff%etin(iunit,nt) + TRunoff%wt(iunit,nt)/theDeltaT)
          if(TRunoff%mt(iunit,nt) > 0._r8) then
             TRunoff%vt(iunit,nt) = -TRunoff%etout(iunit,nt)/TRunoff%mt(iunit,nt)
          end if
       end if
    end if
    TRunoff%dwt(iunit,nt) = TRunoff%etin(iunit,nt) + TRunoff%etout(iunit,nt)
    if (nt == 1) TRunoff%dno3t(iunit) = TRunoff%etinno3(iunit) + TRunoff%etoutno3(iunit)

! check stability
!    if(TRunoff%vt(iunit,nt) < -TINYVALUE .or. TRunoff%vt(iunit,nt) > 30) then
!       write(iulog,*) "Numerical error in subnetworkRouting, ", iunit,nt,TRunoff%vt(iunit,nt)
!    end if

  end subroutine subnetworkRouting

!-----------------------------------------------------------------------

  subroutine mainchannelRouting(iunit, nt, theDeltaT)
  ! !DESCRIPTION: main channel routing
    implicit none    
    integer, intent(in) :: iunit, nt
    real(r8), intent(in) :: theDeltaT    

    if (nt == 1 .and. Tctl%RoutingMethod /= 1) then
       call shr_sys_abort('MOSART nitrate routing is implemented for kinematic wave (RoutingMethod=1) only')
    end if
    if(Tctl%RoutingMethod == 1) then
       call Routing_KW(iunit, nt, theDeltaT)
    else if(Tctl%RoutingMethod == 2) then
       call Routing_MC(iunit, nt, theDeltaT)
    else if(Tctl%RoutingMethod == 3) then
       call Routing_THREW(iunit, nt, theDeltaT)
    else if(Tctl%RoutingMethod == 4) then
       call Routing_DW(iunit, nt, theDeltaT)
    else
       call shr_sys_abort( "mosart: Please check the routing method! There are only 4 methods available." )
    end if

  end subroutine mainchannelRouting

!-----------------------------------------------------------------------

  subroutine Routing_KW(iunit, nt, theDeltaT)
  ! !DESCRIPTION: classic kinematic wave routing method
    implicit none    
    
    integer, intent(in) :: iunit, nt   
    real(r8), intent(in) :: theDeltaT    
    integer  :: k
    real(r8) :: temp_gwl, temp_dwr, temp_gwl0

    ! estimate the inflow from upstream units
    TRunoff%erin(iunit,nt) = 0._r8
    if (nt == 1) TRunoff%erinno3(iunit) = 0._r8

! tcraig, moved this out of the inner main channel loop to before main channel call
! now it's precomputed as TRunoff%eroutUp
!    do k=1,TUnit%nUp(iunit)
!       TRunoff%erin(iunit,nt) = TRunoff%erin(iunit,nt) - TRunoff%erout(TUnit%iUp(iunit,k),nt)
!    end do
    TRunoff%erin(iunit,nt) = TRunoff%erin(iunit,nt) - TRunoff%eroutUp(iunit,nt)
    if (nt == 1) TRunoff%erinno3(iunit) = -TRunoff%eroutno3Up(iunit)

    ! estimate the outflow
    if(TUnit%rlen(iunit) <= 0._r8) then ! no river network, no channel routing
       TRunoff%vr(iunit,nt) = 0._r8
       TRunoff%erout(iunit,nt) = -TRunoff%erin(iunit,nt)-TRunoff%erlateral(iunit,nt)
       if (nt == 1) TRunoff%eroutno3(iunit) = -TRunoff%erinno3(iunit)-TRunoff%erlateralno3(iunit)
    else
       if(TUnit%areaTotal2(iunit)/TUnit%rwidth(iunit)/TUnit%rlen(iunit) > 1e6_r8) then
          TRunoff%erout(iunit,nt) = -TRunoff%erin(iunit,nt)-TRunoff%erlateral(iunit,nt)
          if (nt == 1) TRunoff%eroutno3(iunit) = -TRunoff%erinno3(iunit)-TRunoff%erlateralno3(iunit)
       else
!        !TRunoff%vr(iunit,nt) = CRVRMAN(TUnit%rslp(iunit), TUnit%nr(iunit), TRunoff%rr(iunit,nt))
          TRunoff%vr(iunit,nt) = CRVRMAN_nosqrt(TUnit%rslpsqrt(iunit), TUnit%nr(iunit), TRunoff%rr(iunit,nt))
          TRunoff%erout(iunit,nt) = -TRunoff%vr(iunit,nt) * TRunoff%mr(iunit,nt)
          if (nt == 1) then
             if (TRunoff%wr(iunit,nt) > 0._r8) then
                TRunoff%eroutno3(iunit) = TRunoff%erout(iunit,nt) * TRunoff%no3r(iunit) / TRunoff%wr(iunit,nt)
             else
                TRunoff%eroutno3(iunit) = -TRunoff%erinno3(iunit)-TRunoff%erlateralno3(iunit)
             end if
             if (-TRunoff%eroutno3(iunit) > TINYVALUE .and. TRunoff%no3r(iunit) + &
                 (TRunoff%erlateralno3(iunit) + TRunoff%erinno3(iunit) + TRunoff%eroutno3(iunit)) * theDeltaT < TINYVALUE) then
                TRunoff%eroutno3(iunit) = -(TRunoff%erlateralno3(iunit) + TRunoff%erinno3(iunit) + TRunoff%no3r(iunit) / theDeltaT)
             end if
          end if
          if(-TRunoff%erout(iunit,nt) > TINYVALUE .and. TRunoff%wr(iunit,nt) + &
             (TRunoff%erlateral(iunit,nt) + TRunoff%erin(iunit,nt) + TRunoff%erout(iunit,nt)) * theDeltaT < TINYVALUE) then
             TRunoff%erout(iunit,nt) = -(TRunoff%erlateral(iunit,nt) + TRunoff%erin(iunit,nt) + TRunoff%wr(iunit,nt) / theDeltaT)
             if(TRunoff%mr(iunit,nt) > 0._r8) then
                TRunoff%vr(iunit,nt) = -TRunoff%erout(iunit,nt) / TRunoff%mr(iunit,nt)
             end if
          end if
       end if
    end if

    temp_gwl = TRunoff%qgwl(iunit,nt) * TUnit%area(iunit) * TUnit%frac(iunit)
           
    TRunoff%dwr(iunit,nt) = TRunoff%erlateral(iunit,nt) + TRunoff%erin(iunit,nt) + TRunoff%erout(iunit,nt) + temp_gwl
    if (nt == 1) then
       if (trim(no3numeric_option) == 'explicitBakward') then
          call explicit_backward_dno3r(iunit, nt, theDeltaT)
       else
          TRunoff%dno3r(iunit) = TRunoff%erlateralno3(iunit) + TRunoff%erinno3(iunit) + TRunoff%eroutno3(iunit)
       end if
    end if

    if((TRunoff%wr(iunit,nt)/theDeltaT &
         + TRunoff%dwr(iunit,nt)) < -TINYVALUE) then
       write(iulog,*) 'mosart: ERROR main channel going negative: ', iunit, nt 
       write(iulog,*) theDeltaT, TRunoff%wr(iunit,nt), &
            TRunoff%wr(iunit,nt)/theDeltaT, TRunoff%dwr(iunit,nt), temp_gwl
       write(iulog,*) ' '
       !       call shr_sys_abort('mosart: ERROR main channel going negative')
    endif

! check for stability
!    if(TRunoff%vr(iunit,nt) < -TINYVALUE .or. TRunoff%vr(iunit,nt) > 30) then
!       write(iulog,*) "Numerical error inRouting_KW, ", iunit,nt,TRunoff%vr(iunit,nt)
!    end if

! check for negative wr
!    if(TRunoff%wr(iunit,nt) > 1._r8 .and. (TRunoff%wr(iunit,nt)/theDeltaT + TRunoff%dwr(iunit,nt))/TRunoff%wr(iunit,nt) < -TINYVALUE) then
!       write(iulog,*) 'negative wr!', TRunoff%wr(iunit,nt), TRunoff%dwr(iunit,nt), temp_dwr, temp_gwl, temp_gwl0, theDeltaT
!       stop          
!    end if 

  end subroutine Routing_KW

!-----------------------------------------------------------------------

  subroutine Routing_MC(iunit, nt, theDeltaT)
  ! !DESCRIPTION: Muskingum-Cunge routing method
    implicit none    
    integer, intent(in) :: iunit, nt   
    real(r8), intent(in) :: theDeltaT
   
  end subroutine Routing_MC

!-----------------------------------------------------------------------

  subroutine Routing_THREW(iunit, nt, theDeltaT)
  ! !DESCRIPTION: kinematic wave routing method from THREW model
    implicit none    
    integer, intent(in) :: iunit, nt
    real(r8), intent(in) :: theDeltaT
   
  end subroutine Routing_THREW

!-----------------------------------------------------------------------

  subroutine Routing_DW(iunit, nt, theDeltaT)
  ! !DESCRIPTION: classic diffusion wave routing method
    implicit none    
    integer, intent(in) :: iunit, nt
    real(r8), intent(in) :: theDeltaT
   
  end subroutine Routing_DW

!-----------------------------------------------------------------------

  subroutine updateState_hillslope(iunit,nt)
  ! !DESCRIPTION: update the state variables at hillslope
    implicit none    
    integer, intent(in) :: iunit, nt

    TRunoff%yh(iunit,nt) = TRunoff%wh(iunit,nt) !/ TUnit%area(iunit) / TUnit%frac(iunit) 

  end subroutine updateState_hillslope

!-----------------------------------------------------------------------

  subroutine updateState_subnetwork(iunit,nt)
  ! !DESCRIPTION: update the state variables in subnetwork channel
    implicit none    
    integer, intent(in) :: iunit,nt

       if(TUnit%tlen(iunit) > 0._r8 .and. TRunoff%wt(iunit,nt) > 0._r8) then
          TRunoff%mt(iunit,nt) = GRMR(TRunoff%wt(iunit,nt), TUnit%tlen(iunit)) 
          TRunoff%yt(iunit,nt) = GRHT(TRunoff%mt(iunit,nt), TUnit%twidth(iunit))
          TRunoff%pt(iunit,nt) = GRPT(TRunoff%yt(iunit,nt), TUnit%twidth(iunit))
          TRunoff%rt(iunit,nt) = GRRR(TRunoff%mt(iunit,nt), TRunoff%pt(iunit,nt))
       else
          TRunoff%mt(iunit,nt) = 0._r8
          TRunoff%yt(iunit,nt) = 0._r8
          TRunoff%pt(iunit,nt) = 0._r8
          TRunoff%rt(iunit,nt) = 0._r8
       end if
  end subroutine updateState_subnetwork

!-----------------------------------------------------------------------

  subroutine updateState_mainchannel(iunit, nt)
  ! !DESCRIPTION: update the state variables in main channel
    implicit none    
    integer, intent(in) :: iunit, nt

       if(TUnit%rlen(iunit) > 0._r8 .and. TRunoff%wr(iunit,nt) > 0._r8) then
          TRunoff%mr(iunit,nt) = GRMR(TRunoff%wr(iunit,nt), TUnit%rlen(iunit)) 
          TRunoff%yr(iunit,nt) = GRHR(TRunoff%mr(iunit,nt), TUnit%rwidth(iunit), TUnit%rwidth0(iunit), TUnit%rdepth(iunit))
          TRunoff%pr(iunit,nt) = GRPR(TRunoff%yr(iunit,nt), TUnit%rwidth(iunit), TUnit%rwidth0(iunit), TUnit%rdepth(iunit))
          TRunoff%rr(iunit,nt) = GRRR(TRunoff%mr(iunit,nt), TRunoff%pr(iunit,nt))
       else
          TRunoff%mr(iunit,nt) = 0._r8
          TRunoff%yr(iunit,nt) = 0._r8
          TRunoff%pr(iunit,nt) = 0._r8
          TRunoff%rr(iunit,nt) = 0._r8
       end if
  end subroutine updateState_mainchannel

!-----------------------------------------------------------------------
    
  function CRVRMAN(slp_, n_, rr_) result(v_)
  ! Function for calculating channel velocity according to Manning's equation.
    implicit none
    real(r8), intent(in) :: slp_, n_, rr_ ! slope, manning's roughness coeff., hydraulic radius
    real(r8)             :: v_            ! v_ is  discharge
    
    real(r8) :: ftemp,vtemp

    if(rr_ <= 0._r8) then
       v_ = 0._r8
    else
!tcraig, original code
!       ftemp = 2._r8/3._r8
!       v_ = (rr_**ftemp) * sqrt(slp_) / n_  
!tcraig, produces same answer as original in same time
!       v_ = (rr_**(2._r8/3._r8)) * sqrt(slp_) / n_  

!tcraig, this is faster but NOT bit-for-bit
       v_ = ((rr_*rr_)**(1._r8/3._r8)) * sqrt(slp_) / n_

!debug       if (abs(vtemp - v_)/vtemp > 1.0e-14) then
!debug          write(iulog,*) 'tcx check crvrman ',vtemp, v_
!debug       endif
    end if
    return
  end function CRVRMAN

!-----------------------------------------------------------------------
    
  function CRVRMAN_nosqrt(sqrtslp_, n_, rr_) result(v_)
  ! Function for calculating channel velocity according to Manning's equation.
    implicit none
    real(r8), intent(in) :: sqrtslp_, n_, rr_ ! sqrt(slope), manning's roughness coeff., hydraulic radius
    real(r8)             :: v_            ! v_ is  discharge
    
    real(r8) :: ftemp, vtemp

    if(rr_ <= 0._r8) then
       v_ = 0._r8
    else
!tcraig, original code
!       ftemp = 2._r8/3._r8
!       v_ = (rr_**ftemp) * sqrtslp_ / n_  
!tcraig, produces same answer as original in same time
!       v_ = (rr_**(2._r8/3._r8)) * sqrtslp_ / n_  

!tcraig, this is faster but NOT bit-for-bit
       v_ = ((rr_*rr_)**(1._r8/3._r8)) * sqrtslp_ / n_

!debug       if (abs(vtemp - v_)/vtemp > 1.0e-14) then
!debug          write(iulog,*) 'tcx check crvrman_nosqrt ',vtemp, v_
!debug       endif
    end if
    return
  end function CRVRMAN_nosqrt

!-----------------------------------------------------------------------

  function CREHT(hslp_, nh_, Gxr_, yh_) result(eht_)
  ! Function for overland from hillslope into the sub-network channels
    implicit none
    real(r8), intent(in) :: hslp_, nh_, Gxr_, yh_ ! topographic slope, manning's roughness coeff., drainage density, overland flow depth
    real(r8)                   :: eht_            ! velocity, specific discharge
    
    real(r8) :: vh_
    vh_ = CRVRMAN(hslp_,nh_,yh_)
    eht_ = Gxr_*yh_*vh_
    return
  end function CREHT

!-----------------------------------------------------------------------

  function CREHT_nosqrt(sqrthslp_, nh_, Gxr_, yh_) result(eht_)
  ! Function for overland from hillslope into the sub-network channels
    implicit none
    real(r8), intent(in) :: sqrthslp_, nh_, Gxr_, yh_ ! topographic slope, manning's roughness coeff., drainage density, overland flow depth
    real(r8)                   :: eht_            ! velocity, specific discharge
    
    real(r8) :: vh_
    vh_ = CRVRMAN_nosqrt(sqrthslp_,nh_,yh_)
    eht_ = Gxr_*yh_*vh_
    return
  end function CREHT_nosqrt

!-----------------------------------------------------------------------

  function GRMR(wr_, rlen_) result(mr_)
  ! Function for estimate wetted channel area
    implicit none
    real(r8), intent(in) :: wr_, rlen_      ! storage of water, channel length
    real(r8)             :: mr_             ! wetted channel area
    
    mr_ = wr_ / rlen_
    return
  end function GRMR
  
!-----------------------------------------------------------------------

  function GRHT(mt_, twid_) result(ht_)
  ! Function for estimating water depth assuming rectangular channel
    implicit none
    real(r8), intent(in) :: mt_, twid_      ! wetted channel area, channel width
    real(r8)             :: ht_             ! water depth
    
    if(mt_ <= TINYVALUE) then
       ht_ = 0._r8
    else
       ht_ = mt_ / twid_
    end if
    return
  end function GRHT

!-----------------------------------------------------------------------

  function GRPT(ht_, twid_) result(pt_)
  ! Function for estimating wetted perimeter assuming rectangular channel
    implicit none
    real(r8), intent(in) :: ht_, twid_      ! water depth, channel width
    real(r8)             :: pt_             ! wetted perimeter
    
    if(ht_ <= TINYVALUE) then
       pt_ = 0._r8
    else
       pt_ = twid_ + 2._r8 * ht_
    end if
    return
  end function GRPT

!-----------------------------------------------------------------------

  function GRRR(mr_, pr_) result(rr_)
  ! Function for estimating hydraulic radius
    implicit none
    real(r8), intent(in) :: mr_, pr_        ! wetted area and perimeter
    real(r8)             :: rr_             ! hydraulic radius
    
    if(pr_ <= TINYVALUE) then
       rr_ = 0._r8
    else
       rr_ = mr_ / pr_
    end if
    return
  end function GRRR

!-----------------------------------------------------------------------

  function GRHR(mr_, rwidth_, rwidth0_, rdepth_) result(hr_)
  ! Function for estimating maximum water depth assuming rectangular channel and tropezoidal flood plain
  ! here assuming the channel cross-section consists of three parts, from bottom to up,
  ! part 1 is a rectangular with bankfull depth (rdep) and bankfull width (rwid)
  ! part 2 is a tropezoidal, bottom width rwid and top width rwid0, height 0.1*((rwid0-rwid)/2), assuming slope is 0.1
  ! part 3 is a rectagular with the width rwid0
    implicit none
    real(r8), intent(in) :: mr_, rwidth_, rwidth0_, rdepth_ ! wetted channel area, channel width, flood plain wid, water depth
    real(r8)             :: hr_                             ! water depth
    
    real(r8) :: SLOPE1  ! slope of flood plain, TO DO
    real(r8) :: deltamr_

    SLOPE1 = SLOPE1def
    if(mr_ <= TINYVALUE) then
       hr_ = 0._r8
    else
       if(mr_ - rdepth_*rwidth_ <= TINYVALUE) then ! not flooded
          hr_ = mr_/rwidth_
       else ! if flooded, the find out the equivalent depth
          if(mr_ > rdepth_*rwidth_ + (rwidth_ + rwidth0_)*SLOPE1*((rwidth0_-rwidth_)/2._r8)/2._r8 + TINYVALUE) then
             deltamr_ = mr_ - rdepth_*rwidth_ - (rwidth_ + rwidth0_)*SLOPE1*((rwidth0_ - rwidth_)/2._r8)/2._r8;
             hr_ = rdepth_ + SLOPE1*((rwidth0_ - rwidth_)/2._r8) + deltamr_/(rwidth0_);
          else
             deltamr_ = mr_ - rdepth_*rwidth_;
!           !hr_ = rdepth_ + (-rwidth_+sqrt( rwidth_**2._r8  +4._r8*deltamr_/SLOPE1))*SLOPE1/2._r8
             hr_ = rdepth_ + (-rwidth_+sqrt((rwidth_*rwidth_)+4._r8*deltamr_/SLOPE1))*SLOPE1/2._r8
          end if
       end if
    end if
    return
  end function GRHR
  
!-----------------------------------------------------------------------

  function GRPR(hr_, rwidth_, rwidth0_,rdepth_) result(pr_)
  ! Function for estimating maximum water depth assuming rectangular channel and tropezoidal flood plain
  ! here assuming the channel cross-section consists of three parts, from bottom to up,
  ! part 1 is a rectangular with bankfull depth (rdep) and bankfull width (rwid)
  ! part 2 is a tropezoidal, bottom width rwid and top width rwid0, height 0.1*((rwid0-rwid)/2), assuming slope is 0.1
  ! part 3 is a rectagular with the width rwid0
    implicit none
    real(r8), intent(in) :: hr_, rwidth_, rwidth0_, rdepth_ ! wwater depth, channel width, flood plain wid, water depth
    real(r8)             :: pr_                             ! water depth
    
    real(r8) :: SLOPE1  ! slope of flood plain, TO DO
    real(r8) :: deltahr_
    logical, save :: first_call = .true.

    SLOPE1 = SLOPE1def
    if (first_call) then
       sinatanSLOPE1defr = 1.0_r8/(sin(atan(SLOPE1def)))
    endif
    first_call = .false.

    if(hr_ < TINYVALUE) then
       pr_ = 0._r8
    else
       if(hr_ <= rdepth_ + TINYVALUE) then ! not flooded
          pr_ = rwidth_ + 2._r8*hr_
       else
          if(hr_ > rdepth_ + ((rwidth0_-rwidth_)/2._r8)*SLOPE1 + TINYVALUE) then
             deltahr_ = hr_ - rdepth_ - ((rwidth0_-rwidth_)/2._r8)*SLOPE1
!           !pr_ = rwidth_ + 2._r8*(rdepth_ + ((rwidth0_-rwidth_)/2._r8)*SLOPE1/sin(atan(SLOPE1)) + deltahr_)
             pr_ = rwidth_ + 2._r8*(rdepth_ + ((rwidth0_-rwidth_)/2._r8)*SLOPE1*sinatanSLOPE1defr + deltahr_)
          else
!           !pr_ = rwidth_ + 2._r8*(rdepth_ + (hr_ - rdepth_)/sin(atan(SLOPE1)))
             pr_ = rwidth_ + 2._r8*(rdepth_ + (hr_ - rdepth_)*sinatanSLOPE1defr)
          end if
       end if
    end if
    return
  end function GRPR 
  
!-----------------------------------------------------------------------

  subroutine createFile(nio, fname)
  ! !DESCRIPTION: create a new file. if a file with the same name exists, delete it then create a new one
    implicit none
    character(len=*), intent(in) :: fname ! file name
      integer, intent(in) :: nio            !unit of the file to create
    
    integer :: ios
    logical :: filefound
    character(len=1000) :: cmd
    inquire (file=fname, exist=filefound)
    if(filefound) then
       !cmd = 'rm '//trim(fname)
       !call system(cmd)
       open (unit=nio, file=fname, status="replace", action="write", iostat=ios)
    else
       open (unit=nio, file=fname, status="new", action="write", iostat=ios)
    end if
    if(ios /= 0) then
       call shr_sys_abort( "mosart: cannot create file: "//trim(fname) )
    end if
  end subroutine createFile
  
!-----------------------------------------------------------------------

  subroutine printTest(nio)
  ! !DESCRIPTION: output the simulation results into external files
    implicit none
    integer, intent(in) :: nio        ! unit of the file to print
    
    integer :: IDlist(1:5) = (/151,537,687,315,2080/)
    integer :: ios,ii                    ! flag of io status
            

    write(unit=nio,fmt="(15(e20.11))") TRunoff%etin(IDlist(1),1)/TUnit%area(IDlist(1)), &
      TRunoff%erlateral(IDlist(1),1)/TUnit%area(IDlist(1)), TRunoff%flow(IDlist(1),1), &
      TRunoff%etin(IDlist(2),1)/TUnit%area(IDlist(2)), TRunoff%erlateral(IDlist(2),1)/TUnit%area(IDlist(2)), &
      TRunoff%flow(IDlist(2),1), &
      TRunoff%etin(IDlist(3),1)/TUnit%area(IDlist(3)), TRunoff%erlateral(IDlist(3),1)/TUnit%area(IDlist(3)), &
      TRunoff%flow(IDlist(3),1), &
      TRunoff%etin(IDlist(4),1)/TUnit%area(IDlist(4)), TRunoff%erlateral(IDlist(4),1)/TUnit%area(IDlist(4)), &
      TRunoff%flow(IDlist(4),1), &
      TRunoff%etin(IDlist(5),1)/TUnit%area(IDlist(5)), TRunoff%erlateral(IDlist(5),1)/TUnit%area(IDlist(5)), &
      TRunoff%flow(IDlist(5),1)
    !write(unit=nio,fmt="((a10),(e20.11))") theTime, liqWater%flow(ii)
    !write(unit=nio,fmt="((a10),6(e20.11))") theTime, liqWater%qsur(ii), liqWater%qsub(ii), liqWater%etin(ii)/(TUnit%area(ii)*TUnit%frac(ii)), liqWater%erlateral(ii)/(TUnit%area(ii)*TUnit%frac(ii)), liqWater%erin(ii), liqWater%flow(ii)
    !if(liqWater%yr(ii) > 0._r8) then
    !    write(unit=nio,fmt="((a10),6(e20.11))") theTime, liqWater%mr(ii)/liqWater%yr(ii),liqWater%yr(ii), liqWater%vr(ii), liqWater%erin(ii), liqWater%erout(ii)/(TUnit%area(ii)*TUnit%frac(ii)), liqWater%flow(ii)
      !else
    !    write(unit=nio,fmt="((a10),6(e20.11))") theTime, liqWater%mr(ii)-liqWater%mr(ii),liqWater%yr(ii), liqWater%vr(ii), liqWater%erin(ii), liqWater%erout(ii)/(TUnit%area(ii)*TUnit%frac(ii)), liqWater%flow(ii)
    !end if
    !write(unit=nio,fmt="((a10),7(e20.11))") theTime, liqWater%erlateral(ii)/(TUnit%area(ii)*TUnit%frac(ii)), liqWater%wr(ii),liqWater%mr(ii), liqWater%yr(ii), liqWater%pr(ii), liqWater%rr(ii), liqWater%flow(ii)
    !write(unit=nio,fmt="((a10),7(e20.11))") theTime, liqWater%yh(ii), liqWater%dwh(ii),liqWater%etin(ii), liqWater%vr(ii), liqWater%erin(ii), liqWater%erout(ii)/(TUnit%area(ii)*TUnit%frac(ii)), liqWater%flow(ii)
  
  end subroutine printTest

!-----------------------------------------------------------------------

  subroutine channel_explicit_backward()
  ! Liquid channel step for explicitBakward. Every cell uses the global maximum
  ! substep so the neighbor storages exchanged below are at the same time.

    integer :: iunit, k, kk, unitUp, numDT_rMAX, local_max, ier, maxloc_up, nglob
    real(r8) :: localDeltaT, temp_erout, temp_eroutno3
    real(r8) :: temp_flow(8)
    real(r8), allocatable :: erout_sum(:), eroutno3_sum(:)
    real(r8), allocatable :: no3r_task(:), erout_task(:), eroutno3_task(:), wr_task(:), mr_task(:), Ex_task(:)
    real(r8), allocatable :: eb_wr(:)
    integer , allocatable :: gindex_task(:)

    nglob = rtmlon*rtmlat
    local_max = 0
    if (rtmCTL%endr >= rtmCTL%begr) local_max = maxval(TUnit%numDT_r)
    call mpi_allreduce(local_max, numDT_rMAX, 1, MPI_INTEGER, MPI_MAX, mpicom_rof, ier)
    if (numDT_rMAX < 1) return

    allocate(erout_sum(rtmCTL%begr:rtmCTL%endr), eroutno3_sum(rtmCTL%begr:rtmCTL%endr))
    allocate(eb_gindex(nglob), eb_no3r(nglob), eb_erout(nglob), eb_eroutno3(nglob), eb_mr(nglob), eb_Ex(nglob))
    allocate(gindex_task(nglob), no3r_task(nglob), erout_task(nglob), eroutno3_task(nglob))
    allocate(wr_task(nglob), eb_wr(nglob), mr_task(nglob), Ex_task(nglob))
    erout_sum = 0._r8
    eroutno3_sum = 0._r8
    localDeltaT = Tctl%DeltaT / Tctl%DLevelH2R / numDT_rMAX

    do k = 1, numDT_rMAX
       eb_gindex = 0
       eb_no3r = 0._r8
       eb_erout = 0._r8
       eb_eroutno3 = 0._r8
       eb_mr = 0._r8
       eb_Ex = 0._r8
       gindex_task = 0
       no3r_task = 0._r8
       erout_task = 0._r8
       eroutno3_task = 0._r8
       wr_task = 0._r8
       mr_task = 0._r8
       Ex_task = 0._r8

       do iunit = rtmCTL%begr, rtmCTL%endr
          gindex_task(iunit) = rtmCTL%gindex(iunit)
          no3r_task(rtmCTL%gindex(iunit)) = TRunoff%no3r(iunit)
          erout_task(rtmCTL%gindex(iunit)) = TRunoff%erout(iunit,1)
          eroutno3_task(rtmCTL%gindex(iunit)) = TRunoff%eroutno3(iunit)
          mr_task(rtmCTL%gindex(iunit)) = TRunoff%mr(iunit,1)
          wr_task(rtmCTL%gindex(iunit)) = TRunoff%wr(iunit,1)
          Ex_task(rtmCTL%gindex(iunit)) = LongiDisper(TRunoff%vr(iunit,1), TUnit%nr(iunit), &
               TUnit%rwidth(iunit), TRunoff%yr(iunit,1), TRunoff%rr(iunit,1))
       end do

       call mpi_allreduce(gindex_task, eb_gindex, nglob, MPI_INTEGER, MPI_SUM, mpicom_rof, ier)
       call mpi_allreduce(no3r_task, eb_no3r, nglob, MPI_REAL8, MPI_SUM, mpicom_rof, ier)
       call mpi_allreduce(erout_task, eb_erout, nglob, MPI_REAL8, MPI_SUM, mpicom_rof, ier)
       call mpi_allreduce(eroutno3_task, eb_eroutno3, nglob, MPI_REAL8, MPI_SUM, mpicom_rof, ier)
       call mpi_allreduce(mr_task, eb_mr, nglob, MPI_REAL8, MPI_SUM, mpicom_rof, ier)
       call mpi_allreduce(wr_task, eb_wr, nglob, MPI_REAL8, MPI_SUM, mpicom_rof, ier)
       call mpi_allreduce(Ex_task, eb_Ex, nglob, MPI_REAL8, MPI_SUM, mpicom_rof, ier)

       do iunit = rtmCTL%begr, rtmCTL%endr
          if (TUnit%nUp(iunit) > 0) then
             temp_flow = 0._r8
             do kk = 1, TUnit%nUp(iunit)
                unitUp = TUnit%iUp(iunit,kk)
                if (unitUp > 0 .and. unitUp <= nglob) then
                   if (eb_gindex(unitUp) > 0) temp_flow(kk) = eb_wr(eb_gindex(unitUp))
                end if
             end do
             maxloc_up = maxloc(temp_flow(1:TUnit%nUp(iunit)), dim=1)
             TUnit%maxUp(iunit) = TUnit%iUp(iunit, maxloc_up)
          end if
       end do

       do iunit = rtmCTL%begr, rtmCTL%endr
          if (TUnit%mask(iunit) > 0) then
             call mainchannelRouting(iunit, 1, localDeltaT)
             TRunoff%wr(iunit,1) = TRunoff%wr(iunit,1) + TRunoff%dwr(iunit,1) * localDeltaT
             TRunoff%no3r(iunit) = TRunoff%no3r(iunit) + TRunoff%dno3r(iunit) * localDeltaT
             call UpdateState_mainchannel(iunit, 1)
             erout_sum(iunit) = erout_sum(iunit) + TRunoff%erout(iunit,1)
             eroutno3_sum(iunit) = eroutno3_sum(iunit) + TRunoff%eroutno3(iunit)
          end if
       end do
    end do

    do iunit = rtmCTL%begr, rtmCTL%endr
       if (TUnit%mask(iunit) > 0) then
          temp_erout = erout_sum(iunit) / numDT_rMAX
          temp_eroutno3 = eroutno3_sum(iunit) / numDT_rMAX
          TRunoff%erout(iunit,1) = temp_erout
          TRunoff%eroutno3(iunit) = temp_eroutno3
          TRunoff%flow(iunit,1) = TRunoff%flow(iunit,1) - temp_erout
          TRunoff%flowno3(iunit) = TRunoff%flowno3(iunit) - temp_eroutno3
       end if
    end do

    deallocate(erout_sum, eroutno3_sum)
    deallocate(eb_gindex, eb_no3r, eb_erout, eb_eroutno3, eb_mr, eb_Ex)
    deallocate(gindex_task, no3r_task, erout_task, eroutno3_task, wr_task, eb_wr, mr_task, Ex_task)

  end subroutine channel_explicit_backward

!-----------------------------------------------------------------------

  subroutine explicit_backward_dno3r(iunit, nt, theDeltaT)

    integer, intent(in) :: iunit, nt
    real(r8), intent(in) :: theDeltaT
    integer :: unitUp, gindexUp
    real(r8) :: temp_Ex, conc, temp_vfden, temp_rctn2
    real(r8) :: temp_advctn, temp_diffsn, sinksource
    real(r8) :: alpha, beta, gamma, rctnDecay, rctden

    unitUp = TUnit%maxUp(iunit)
    if (TRunoff%wr(iunit,nt) > 0._r8 .and. TUnit%rlen(iunit) > 0._r8 .and. &
        TUnit%dnID(iunit) > 0 .and. TUnit%nUp(iunit) > 0 .and. unitUp > 0 .and. &
        TRunoff%yr(iunit,nt) > 0._r8) then
       gindexUp = 0
       if (unitUp <= size(eb_gindex)) gindexUp = eb_gindex(unitUp)
       if (gindexUp < 1 .or. rtmCTL%dsig(iunit) < 1 .or. rtmCTL%dsig(iunit) > size(eb_no3r)) then
          TRunoff%dno3r(iunit) = TRunoff%erlateralno3(iunit) + TRunoff%erinno3(iunit) + TRunoff%eroutno3(iunit)
       else
          temp_Ex = LongiDisper(TRunoff%vr(iunit,nt), TUnit%nr(iunit), TUnit%rwidth(iunit), &
               TRunoff%yr(iunit,nt), TRunoff%rr(iunit,nt))
          if (TRunoff%no3r(iunit) > 0._r8) then
             conc = (TRunoff%no3r(iunit) / TRunoff%wr(iunit,nt)) * 1000._r8
             temp_vfden = 0.01_r8 * (10._r8**(-0.493_r8*log10(conc) - 2.975_r8))
          else
             temp_vfden = 0._r8
          end if
          temp_rctn2 = temp_vfden / TRunoff%yr(iunit,nt)
          temp_advctn = (-eb_erout(gindexUp)*eb_no3r(gindexUp) + TRunoff%erout(iunit,nt)*TRunoff%no3r(iunit)) &
               / TRunoff%wr(iunit,nt)
          temp_diffsn = (TRunoff%mr(iunit,nt)*temp_Ex) * (eb_no3r(rtmCTL%dsig(iunit)) - TRunoff%no3r(iunit)) &
               / TRunoff%wr(iunit,nt) / TUnit%rlen(iunit) + &
               (eb_mr(gindexUp)*eb_Ex(gindexUp)) * (eb_no3r(gindexUp) - TRunoff%no3r(iunit)) &
               / TRunoff%wr(iunit,nt) / TUnit%rlen(iunit)
          sinksource = -TRunoff%eroutno3Up(iunit) + eb_eroutno3(gindexUp)
          alpha = (eb_erout(gindexUp) - TRunoff%erout(iunit,nt)) / TRunoff%wr(iunit,nt)
          beta  = (TRunoff%mr(iunit,nt)*temp_Ex + eb_mr(gindexUp)*eb_Ex(gindexUp)) &
               / TRunoff%wr(iunit,nt) / TUnit%rlen(iunit)
          gamma = (-eb_erout(gindexUp) - TRunoff%erout(iunit,nt)) / 2._r8 / TRunoff%wr(iunit,nt)
          rctnDecay = temp_rctn2
          if (alpha >= 0._r8) rctnDecay = min(temp_rctn2, alpha)
          rctden = min(rctnDecay * TRunoff%no3r(iunit), &
               TRunoff%no3r(iunit)/theDeltaT + (TRunoff%erlateralno3(iunit) + temp_advctn + sinksource + temp_diffsn))
          TRunoff%dno3r(iunit) = TRunoff%erlateralno3(iunit) + temp_advctn + sinksource + temp_diffsn - rctden
          rtmCTL%nTstep(iunit) = rtmCTL%nTstep(iunit) + 1
          if ((rctnDecay > alpha) .or. (theDeltaT > 1._r8/(beta + gamma + temp_rctn2/2._r8))) then
             rtmCTL%nInstable(iunit) = rtmCTL%nInstable(iunit) + 1
          end if
          if (temp_Ex > 0._r8) then
             if ((TUnit%rlen(iunit)**2._r8 / 2._r8 / temp_Ex) < theDeltaT) then
                write(iulog,*) 'mosart: Stability Issue in main channel NO3: DT>DX^2/2Ex', iunit, nt
             end if
          end if
       end if
    else
       TRunoff%dno3r(iunit) = TRunoff%erlateralno3(iunit) + TRunoff%erinno3(iunit) + TRunoff%eroutno3(iunit)
    end if

  end subroutine explicit_backward_dno3r

!-----------------------------------------------------------------------

  function LongiDisper(U_, n_, B_, H_, rr_) result(Kx_)
  ! Channel longitudinal dispersion from Deng et al. (2001), m2/s.

    real(r8), intent(in) :: U_, n_, B_, H_, rr_
    real(r8) :: Kx_
    real(r8) :: eps_t0, u_star

    if (U_ == 0._r8 .or. H_ == 0._r8 .or. rr_ <= 0._r8) then
       Kx_ = 0._r8
    else
       u_star = U_ * n_ * sqrt(9.807_r8 * rr_**(-1._r8/3._r8))
       eps_t0 = 0.145_r8 + (1._r8/3520._r8)*(U_/u_star)*((B_/H_)**1.38_r8)
       Kx_ = (0.15_r8/8._r8/eps_t0) * ((B_/H_)**(5._r8/3._r8)) * ((U_/u_star)**2._r8) * H_ * u_star
    end if

  end function LongiDisper

!-----------------------------------------------------------------------

end MODULE MOSART_physics_mod


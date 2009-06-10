!                                 Notice
!  This program was prepared by the University of California (University)
!  under Contract W-7405-ENG-36 with the U.S. Department of Energy (DOE).
!  All rights in the program are reserved by DOE on behalf of the Government
!  and the University pursuant to the contract. You are authorized to use
!  this program for Government purposes but it is not to be released or
!  distributed to the public.
!  NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY,
!  NOR THE UNIVERSITY OF CALIFORNIA, NOR ANY OF THEIR EMPLOYEES,
!  MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY
!  OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS, OF
!  ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
!  THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
!





! vim: foldmethod=marker:foldmarker=<<,>>:foldlevel=0:foldenable:nowrap

      module Evaporation
         implicit none
      !   private
         integer,parameter::SINGLE=kind(1.0e0)
         integer,parameter::DOUBLE=kind(1.0d0)
         integer,parameter::REALKIND=SINGLE
         ! Constants <<
         real(REALKIND),parameter::ONE  =1.0
         real(REALKIND),parameter::TWO  =2.0
         real(REALKIND),parameter::THREE=3.0
         real(REALKIND),parameter::FOUR =4.0
         real(REALKIND),parameter::onethird=ONE/THREE
         real(REALKIND),parameter::epi      =3.14159265358979323846
         real(REALKIND),parameter::twopi   =TWO*epi
         real(REALKIND),parameter::threepi =THREE*epi
         real(REALKIND),parameter::fourpi  =FOUR*epi
         real(REALKIND),parameter::k_Boltzmann=1.380658e-23        ! J/K (CRC-HCP, p.1-1)
         real(REALKIND),parameter::AvogadrosNumber = 6.0221367e23  ! molecules/mol
         real(REALKIND),parameter::GasConstant = 8.314151          ! J/(mol*K)
         real(REALKIND),parameter::pressureConversion = 101325     ! Pa/atm
         real(REALKIND),parameter::densityConversion = 1e6         ! (g/m^3)/(g/ml)
         
         real(REALKIND),parameter::mw_air = 28.966  ! g/mol
         real(REALKIND),parameter::kinematic_viscosity_air = 1.47e-5 !m^2/s @ 288K
         real(REALKIND),parameter::heat_cap_air= 1005 !J/(kg*K) @ 288K
         real(REALKIND),parameter::thermal_conductivity_air = 0.02535 !W/(mK) @ 288K
         real(REALKIND),parameter::mfp_air = 6.6e-8 !m
         real(REALKIND),parameter::atomic_volume_air = 20.1 !Unclear what the units are possibly cubic angstroms see Arnold 1930 and Fuller et al. 1966
         
      end module Evaporation
      ! >>
      
         subroutine evaporate
            use variables
            use Evaporation
            implicit none
            
            integer::nvapor_particles,evaporation_flag
            real(REALKIND)::md_collision_coef(2)
            real(REALKIND)::r,particle_speed,Tval,rh
            real(REALKIND)::D_agent,D,p_agent,R_agent,denom
            real(REALKIND)::invKn_agent,coll_geom_term,coll_sticking_term,cs_correction
            real(REALKIND)::ReynoldsNum,SchmidtNum,PrandtlNum
            real(REALKIND)::ventilation,ventilation_correction
            real(REALKIND)::invKn_energy,conductivity,minimum_diameter
            real(REALKIND)::latent_heat_agent,saturation_vapor_press_agent
            real(REALKIND)::surf_tension_wl_air,ps_agent,draddt,dmdt,density_air
            real(REALKIND),allocatable::evap_mass(:,:,:)
            integer,allocatable::vapor_particles(:),ipart(:),jpart(:),kpart(:)
            
            if(ibioslurryflag .eq. 1)then
               !properties for water
               agent_mw_vapor=18.15 ! g/mol
               agent_mw_liquid=18.153 ! g/mol
               agent_mfp=4.3e-8 ! m
               liquidDensity = 1 !g/cm^3
            endif
            
            print*,'Evaporating particles at t =',t
            R_agent=GasConstant/agent_mw_vapor
            
            ! From Pruppacher and Klett
            md_collision_coef(1) = 1.33
            md_collision_coef(2) = 0.70        
            
            evaporation_flag=0
            
            allocate(ipart(mpart),jpart(mpart),kpart(mpart))
            if(itwophaseflag .eq. 1)then
               call heatOfVaporization(Pcritical,Tcritical,Tboiling,latent_heat_agent)
               allocate(evap_mass(nx1,ny1,nz1),vapor_particles(mpart))
               evap_mass(:,:,:)=0.
               vaporfield(:,:,:)=0.
               do idx=1,mp
                  ipart(idx) = int(xp(idx)/dx)+1
                  jpart(idx) = int(yp(idx)/dy)+1
                  kpart(idx) = int(zp(idx)/dz)+2
                  if(dpm(idx) .eq. 0. .and. tstrt(idx) .le. t-0.001*delta .and. zp(idx) .gt. 0.)then
                     call interpolate(xp(idx),yp(idx),zp(idx),temperature,Tval)
                     vaporfield(ipart(idx),jpart(idx),kpart(idx))=vaporfield(ipart(idx),jpart(idx),kpart(idx))&
                           +R_agent*pmass(idx)*Tval/(pressureConversion*dx*dy*dz)
                  endif
               enddo
            else
               do idx=1,mp
                  ipart(idx) = int(xp(idx)/dx)+1
                  jpart(idx) = int(yp(idx)/dy)+1
                  kpart(idx) = int(zp(idx)/dz)+2
               enddo
            endif
            do idx=1,mp
               if(itwophaseflag .eq. 1)then
                  minimum_diameter=0.
               else
                  minimum_diameter=dpm_min(idx)
               endif
               if(dpm(idx) .gt. minimum_diameter .and. zp(idx) .gt. 0. &
                     .and. tstrt(idx) .le. t-0.001*delta)then
                  evaporation_flag=1
                  r = 5e-7*dpm(idx)
                  call compute_velocity(ipart(idx),jpart(idx),kpart(idx),&
                                          uran(idx),vran(idx),wran(idx),particle_speed)
                  call interpolate(xp(idx),yp(idx),zp(idx),temperature,Tval)
                  density_air = 0.001*atm_pressure*pressureConversion*mw_air/(GasConstant*Tval)
                  
                  if(ibioslurryflag .eq. 1)then
                     latent_heat_agent = (2.501e3 - 2.370*(Tval - 273.15))*agent_mw_liquid
                  endif
                  
                  call diffusionCoefficient(ibioslurryflag,atm_pressure,Tval,agent_mw_vapor, &
                        atomic_volume_agent,D_agent)
                  
                  ! Corrections to molecular diffusion coefficient <<
                  invKn_agent = r/agent_mfp
                  ! - Collision geometry term
                  coll_geom_term = md_collision_coef(1) + md_collision_coef(2)*invKn_agent &
                                    /(1.0+invKn_agent)
                  ! - Sticking efficiency term
                  coll_sticking_term = 4.0*(1.0-agent_stick_eff)/(3.0*agent_stick_eff)
                  
                  cs_correction = 1.0 + (coll_geom_term + coll_sticking_term)/invKn_agent
                  
                  ! - Ventilation term (for liquid droplets)
                  ReynoldsNum = 2.0*r*particle_speed/kinematic_viscosity_air
                  SchmidtNum  = kinematic_viscosity_air/D_agent
                  ventilation = sqrt(ReynoldsNum)*(SchmidtNum**onethird)
                  if(ventilation .le. 1.4)then
                     ventilation_correction = 1.0 + 0.108*ventilation*ventilation
                  else
                     ventilation_correction = 0.78 + 0.308*ventilation
                  endif
                  
                  ! Compute corrected molecular diffusion coefficient
                  D = D_agent*ventilation_correction/cs_correction
                  ! >>
                  
                  ! Corrections to thermal conductivity of air <<
                  invKn_energy = r/mfp_air
                  !thermal_accomodation = (T_molec_out - T_agent)/(T_surface - T_agent)
                  ! - Collision geometry term
                  coll_geom_term = md_collision_coef(1) + md_collision_coef(2)*invKn_energy &
                                    /(ONE+invKn_energy)
                  ! - Sticking efficiency term
                  coll_sticking_term = FOUR*(ONE-thermal_accomodation)/(THREE*thermal_accomodation)
                  
                  cs_correction = ONE + (coll_geom_term + coll_sticking_term)/invKn_energy
                  
                  ! - Ventilation term (for liquid droplets)
                  PrandtlNum  = density_air*kinematic_viscosity_air*heat_cap_air/thermal_conductivity_air
                  ventilation = sqrt(ReynoldsNum)*(PrandtlNum**onethird)
                  if(ventilation .le. 1.4)then
                     ventilation_correction = 1.0 + 0.108*ventilation*ventilation
                  else
                     ventilation_correction = 0.78 + 0.308*ventilation
                  endif
                  
                  ! Compute corrected thermal conductivity
                  conductivity = thermal_conductivity_air*ventilation_correction/cs_correction
                  ! >>
                  
                  ! Corrections to saturation vapor pressure <<

                  call saturationVaporPressure(ibioslurryflag,Tboiling,Tval,latent_heat_agent, &
                        saturation_vapor_press_agent)
                  
                  call surfaceTension(ibioslurryflag,Tcritical,Tval,liquidDensity,agent_mw_liquid, &
                        surfaceTensionCoef,surf_tension_wl_air)
                  
                  ps_agent = saturation_vapor_press_agent &
                          *exp(2.0*surf_tension_wl_air*agent_mw_liquid/(r*GasConstant*Tval*densityConversion*liquidDensity))
                  ! >>
                  call interpolate(xp(idx),yp(idx),zp(idx),vaporfield,rh)
                  if(itwophaseflag .eq. 1)then
                     p_agent=rh
                  else
                     p_agent=rh*saturation_vapor_press_agent
                  endif
                  denom = D*latent_heat_agent*(pressureConversion*ps_agent) &
                        *(latent_heat_agent/(GasConstant*Tval)-ONE)/(agent_mw_vapor*conductivity*Tval) + R_agent*Tval
                  denom = denom*densityConversion*liquidDensity*r
                  draddt = D*pressureConversion*(p_agent-ps_agent)/denom
                  if(2.*(r+draddt*delta) .le. 1e-6*minimum_diameter)then
                     dpm(idx) = minimum_diameter
                  else
                     if(itwophaseflag .eq. 1)then
                        dmdt=pmass(idx)*THREE*draddt/r
                        evap_mass(ipart(idx),jpart(idx),kpart(idx))=evap_mass(ipart(idx),jpart(idx),kpart(idx))-dmdt*delta
                        pmass(idx)=pmass(idx)+dmdt*delta
                     endif
                     r = r + draddt*delta
                     dpm(idx) = 2e6*r
                  endif
               endif
            enddo ! idx
            if(itwophaseflag .eq. 1)then
               do k=2,nz1
                  do j=1,ny1
                     do i=1,nx1
                        if(evap_mass(i,j,k) .gt. 0)then
                           nvapor_particles=0
                           do idx=1,mp
                              if(dpm(idx) .eq. 0. .and. i .eq. ipart(idx) .and. j .eq. jpart(idx) &
                                    .and. k .eq. kpart(idx) .and. tstrt(idx) .le. t-0.001*delta)then
                                 nvapor_particles=nvapor_particles+1
                                 vapor_particles(nvapor_particles)=idx
                              endif
                           enddo
                           if(nvapor_particles .gt. 0)then
                              do idx=1,nvapor_particles
                                 pmass(vapor_particles(idx))=pmass(vapor_particles(idx)) &
                                       +evap_mass(i,j,k)/real(nvapor_particles)
                              enddo
                           else
                              nptot=nptot+1
                              if(nptot .gt. mpart)then
                                 call addparticles
                              endif
                              pmass(nptot)=evap_mass(i,j,k)
                              xp(nptot)=(real(i)-0.5)*dx
                              yp(nptot)=(real(j)-0.5)*dy
                              zp(nptot)=(real(k)-1.5)*dz
                              uran(nptot)=0.
                              vran(nptot)=0.
                              wran(nptot)=0.
                              dpm_min(nptot)=0.
                              dpm(nptot)=0.
                              weight(nptot)=1.
                              weightd(nptot)=1.
                              uranbp(nptot)=0.
                              vranbp(nptot)=0.
                              wranbp(nptot)=0.
                              xps(nptot)=xp(nptot)
                              yps(nptot)=yp(nptot)
                              zps(nptot)=zp(nptot)
                              source_num(nptot)=0.
                              xnear(nptot)=0.
                              ynear(nptot)=0.
                              znear(nptot)=0.
                              tstrt(nptot)=t
                           endif
                        endif
                     enddo
                  enddo
               enddo
               mp=nptot
               deallocate(evap_mass,vapor_particles,ipart,jpart,kpart)
            endif
            if(evaporation_flag .eq. 0)then
               itwophaseflag=0
               ibioslurryflag=0
            endif
            return
         end subroutine evaporate ! >>
         
         subroutine addparticles
            use variables
            use Evaporation
            implicit none
            
            real(REALKIND),allocatable::xpt(:),ypt(:),zpt(:),urant(:),vrant(:),wrant(:),dwallpt(:)
            real(REALKIND),allocatable::dpm_mint(:),dpmt(:),weightt(:),weightdt(:),pmasst(:)
            real(REALKIND),allocatable::uranbpt(:),vranbpt(:),wranbpt(:),xpst(:),ypst(:),zpst(:)
            real(REALKIND),allocatable::source_numt(:),xneart(:),yneart(:),zneart(:),tstrtt(:)
            
            integer mpartold
            
            allocate(xpt(mpart),ypt(mpart),zpt(mpart),urant(mpart),vrant(mpart),wrant(mpart),dwallpt(mpart))
            allocate(dpm_mint(mpart),dpmt(mpart),weightt(mpart),weightdt(mpart),pmasst(mpart))
            allocate(uranbpt(mpart),vranbpt(mpart),wranbpt(mpart),xpst(mpart),ypst(mpart),zpst(mpart))
            allocate(source_numt(mpart),xneart(mpart),yneart(mpart),zneart(mpart),tstrtt(mpart))
            
            print*,'Adding 1000 more particles'
            
            xpt=xp
            ypt=yp
            zpt=zp
            urant=uran
            vrant=vran
            wrant=wran
            dpm_mint=dpm_min
            dpmt=dpm
            weightt=weight
            weightdt=weightd
            pmasst=pmass
            uranbpt=uranbp
            vranbpt=vranbp
            wranbpt=wranbp
            xpst=xps
            ypst=yps
            zpst=zps
            source_numt=source_num
            xneart=xnear
            yneart=ynear
            zneart=znear
            tstrtt=tstrt
            dwallpt=dwallp
            
            mpartold=mpart
            mpart=mpart+1000
            deallocate(xp,yp,zp,uran,vran,wran,dpm_min,dpm,weight,weightd,pmass,tstrt,dwallp)
            deallocate(uranbp,vranbp,wranbp,xps,yps,zps,source_num,xnear,ynear,znear)
            allocate(xp(mpart),yp(mpart),zp(mpart),uran(mpart),vran(mpart),wran(mpart),dwallp(mpart),&
                     dpm_min(mpart),dpm(mpart),weight(mpart),weightd(mpart),pmass(mpart), &
                     uranbp(mpart),vranbp(mpart),wranbp(mpart),xps(mpart),yps(mpart),zps(mpart), &
                     source_num(mpart),xnear(mpart),ynear(mpart),znear(mpart),tstrt(mpart))
            xp(1:mpartold)=xpt
            yp(1:mpartold)=ypt
            zp(1:mpartold)=zpt
            uran(1:mpartold)=urant
            vran(1:mpartold)=vrant
            wran(1:mpartold)=wrant
            dpm_min(1:mpartold)=dpm_mint
            dpm(1:mpartold)=dpmt
            weight(1:mpartold)=weightt
            weightd(1:mpartold)=weightdt
            pmass(1:mpartold)=pmasst
            uranbp(1:mpartold)=uranbpt
            vranbp(1:mpartold)=vranbpt
            wranbp(1:mpartold)=wranbpt
            xps(1:mpartold)=xpst
            yps(1:mpartold)=ypst
            zps(1:mpartold)=zpst
            source_num(1:mpartold)=source_numt
            xnear(1:mpartold)=xneart
            ynear(1:mpartold)=yneart
            znear(1:mpartold)=zneart
            tstrt(1:mpartold)=tstrtt
            dwallp(1:mpartold)=dwallpt
            deallocate(xpt,ypt,zpt,urant,vrant,wrant,dpm_mint,dpmt,weightt,weightdt,pmasst,tstrtt, &
                       uranbpt,vranbpt,wranbpt,xpst,ypst,zpst,source_numt,xneart,yneart,zneart,dwallpt)
            return
         end subroutine addparticles
      
         subroutine compute_velocity(ipart,jpart,kpart,uranp,vranp,wranp,particle_speed) ! <<
            use variables
            use Evaporation
            implicit none
            
            integer,intent(in)::ipart,jpart,kpart
            real(REALKIND),intent(in)::uranp,vranp,wranp
            real(REALKIND),intent(out)::particle_speed
            real(REALKIND)::particle_velocity(3)
            
!            real(REALKIND)::ump,vmp,wmp
            dwall=dwallp(idx)
            im=ipart
            jm=jpart
            km=kpart
            call windterp
!            call interpolate(xpos,ypos,zpos,u,ump)
!            call interpolate(xpos,ypos,zpos,v,vmp)
!            call interpolate(xpos,ypos,zpos,w,wmp)
            particle_velocity(1) = uint + uranp
            particle_velocity(2) = vint + vranp
            particle_velocity(3) = wint + wranp
            particle_speed = sqrt(sum(particle_velocity**2))
            return
         end subroutine compute_velocity ! >>
      
         subroutine interpolate(xpos,ypos,zpos,field,value) ! <<
            use variables
            use Evaporation
            implicit none
            
            real(REALKIND),intent(in)::xpos,ypos,zpos
            real(REALKIND),intent(in)::field(nx1,ny1,nz1)
            real(REALKIND),intent(out)::value
            
            real(REALKIND)::x1,x2,y1,y2,z1,z2,a1,a2,a3,a4,a5,a6,a7,a8
            if((xpos .gt. real(nx1,REALKIND)*dx) .or. (xpos .lt. 0.).or. &
               (ypos .gt. real(ny1,REALKIND)*dy) .or. (ypos .lt. 0.).or. &
               (zpos .gt. real(nz1-1,REALKIND)*dz) .or. (zpos .lt. 0.))then
              print *,'Interpolation error: (x,y,z) is not in the grid domain.'
              print*,xpos,ypos,zpos
              stop
            else
              im = nint(xpos/dx)
              jm = nint(ypos/dy)
              km = nint(zpos/dz)+1
              if(im .gt. 0 .and. im .lt. nx1 &
                    .and. jm .gt. 0 .and. jm .lt. ny1 &
                    .and. km .gt. 1 .and. km .lt. nz1)then
                 x1 = (real(im  ,REALKIND)-0.5)*dx
                 x2 = (real(im+1,REALKIND)-0.5)*dx
                 y1 = (real(jm  ,REALKIND)-0.5)*dy
                 y2 = (real(jm+1,REALKIND)-0.5)*dy
                 z1 = (real(km  ,REALKIND)-1.5)*dz
                 z2 = (real(km+1,REALKIND)-1.5)*dz
                 
                 a1 = (x2-xpos)*(y2-ypos)*(z2-zpos)/(dx*dy*dz);
                 a2 = (xpos-x1)*(y2-ypos)*(z2-zpos)/(dx*dy*dz);
                 a3 = (x2-xpos)*(ypos-y1)*(z2-zpos)/(dx*dy*dz);
                 a4 = (xpos-x1)*(ypos-y1)*(z2-zpos)/(dx*dy*dz);
                 a5 = (x2-xpos)*(y2-ypos)*(zpos-z1)/(dx*dy*dz);
                 a6 = (xpos-x1)*(y2-ypos)*(zpos-z1)/(dx*dy*dz);
                 a7 = (x2-xpos)*(ypos-y1)*(zpos-z1)/(dx*dy*dz);
                 a8 = (xpos-x1)*(ypos-y1)*(zpos-z1)/(dx*dy*dz);
                 
                 value = a1*field(im,jm  ,km  ) + a2*field(im+1,jm  ,km  ) &
                       + a3*field(im,jm+1,km  ) + a4*field(im+1,jm+1,km  ) &
                       + a5*field(im,jm  ,km+1) + a6*field(im+1,jm  ,km+1) &
                       + a7*field(im,jm+1,km+1) + a8*field(im+1,jm+1,km+1)
              else
                 im=max(im,1)
                 im=min(im,nx1)
                 jm=max(jm,1)
                 jm=min(jm,ny1)
                 km=max(km,2)
                 km=min(km,nz1)
                 value = field(im,jm,km)
              endif
            endif
            return
         end subroutine interpolate ! >>
         
         subroutine diffusionCoefficient(ibioslurryflag,P,T,molecularWeight,atomicVolume,D)
            use Evaporation
            implicit none
            
            integer,intent(in)::ibioslurryflag
            real(REALKIND),intent(in)::P,T,molecularWeight,atomicVolume
            real(REALKIND),intent(out)::D
            
            if(ibioslurryflag .eq. 1)then
               !From Hall and Pruppacher (1976) gives diffusivity of water vapor in air between +/- 40 C in m^2/s.
               D = (0.211e-4)*(1/P)*((T/273.15)**1.94)
            else
               !Method from Fuller et al. (1966) to determine the diffusion coefficient of an arbitrary agent in air
               !T in K;mw in g/mol;P in atm; The units of atomic volume are unclear but there is a table for
               !calculating it for typical organic compounds in the literature.  Yields diffusion coefficient in m^2/s.
               D = (1e-7)*(T**1.75)*sqrt((mw_air+molecularWeight)/(mw_air*molecularWeight)) &
                        /(P*(((atomic_volume_air**onethird)+(atomicVolume**onethird))**2))
            endif
            return
         end subroutine diffusionCoefficient
         
         subroutine heatOfVaporization(Pc,Tc,Tb,Hv)
            use Evaporation
            implicit none
            
            real(REALKIND),intent(in)::Pc,Tc,Tb
            real(REALKIND),intent(out)::Hv
            real(REALKIND)::Kkl
            !From Klein (1949)
            !Pc is the critical pressure in atm
            !Tc is the critical temperature in K
            !Tb is the boiling temperature at 1 atm in K
            !Returns the heat of vaporization in the same units as the GasConstant*Kelvin 
            if(Tb<200)then
               Kkl=1.02
            elseif(Tb>300)then
               Kkl=1.045
            else
               Kkl=1.04
            endif
            Hv = GasConstant*Kkl*Tb*log(Pc)*sqrt(1.-1./(Pc*((Tb/Tc)**3)))/(1.-(Tb/Tc))
            return
         end subroutine heatOfVaporization
         
         subroutine saturationVaporPressure(ibioslurryflag,Tb,T,Hvb,Psat)
            use Evaporation
            implicit none
            
            integer,intent(in)::ibioslurryflag
            real(REALKIND),intent(in)::Tb,T,Hvb
            real(REALKIND),intent(out)::Psat
            real(REALKIND)::C,Ah,Bh,T0,ps_wv_0
            
            if(ibioslurryflag .eq. 1)then
               ! - Saturation vapor pressure of water over a liquid surface
               ! - Clausius-Clapeyron equation (p. 41)
               Ah = 3.14839e3 ! J/g
               Bh = 2.370     ! J/(g K)
               T0 = 273.15    ! K
               ps_wv_0 = 611.2 ! Pa
               Psat = ps_wv_0*exp(18.15*(Ah*(T-T0)/(T*T0) + Bh*log(T/T0))/GasConstant)/pressureConversion
            else
               !from Antoine (1888)
               !Tb is the boiling temperature at 1 atm in K
               !T is the tempurature in K
               !Hvb is the heat of vaporization in the units of GasConstant*K
               !Returns the saturation vapor pressure in atm 
               C=-18+0.19*Tb
               Psat = exp((Hvb*((Tb-C)**2)/(0.97*GasConstant*(Tb**2)))*((1./(Tb-C))-(1./(T-C))))
            endif
            return
         end subroutine saturationVaporPressure
         
         subroutine surfaceTension(ibioslurryflag,Tc,T,density,molecularWeight,surfaceTensionCoefficient,sigma)
            use Evaporation
            implicit none
            
            integer,intent(in)::ibioslurryflag
            real(REALKIND),intent(in)::Tc,T,density,molecularWeight,surfaceTensionCoefficient
            real(REALKIND),intent(out)::sigma
            real(REALKIND)::a(7),T0
            integer::n
            
            if(ibioslurryflag .eq. 1)then
               T0=273.15
               if(T .ge. T0)then
                  sigma = 0.001*(76.1 - 0.155*(T-T0)) ! N/m
               else
                  a(1) = 75.93
                  a(2) = 0.115
                  a(3) = 0.06818
                  a(4) = 6.511e-3
                  a(5) = 2.933e-4
                  a(6) = 6.283e-6
                  a(7) = 5.285e-8
                  sigma = 0.0
                  do n=1,7
                     sigma = sigma+a(n)*((T-T0)**n)
                  enddo
                  sigma=0.001*sigma !convert to N/m
               endif 
            else
               !from Eötvös , Ann. der Physik., 263, 448 (1886).
               !Tc is the critical temperature in K
               !T is the tempurature in K
               !density is in g/cm^3
               !molecularWeight is in g/mol
               !surfaceTensionCoefficient = 2.12 !erg/K for most liquids (non-polar)
               !surfaceTensionCoefficient = 1.03 !erg/K for water because it is a polar solvent
               !Returns the surface tension in N/m
               sigma = 0.001*surfaceTensionCoefficient*((molecularWeight/density)**(-2./3.))*(Tc-T)
            endif
            return
         end subroutine surfaceTension
!      end module Evaporation
      
      


























      module variables  
!Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!module with all parameter and array specifications that are passed between main
!and subroutines
!created MDW 10/09/02
!arrays are allocatable
! see main for additional revisions
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         implicit none
         save
         integer inumbuild,inumbuildfile,inumveg,inumgarage
         integer i,j,k,im,jm,km,m,mm,nx,ny,nz,inmode,mrel,mpart,iparticle_resample
         integer nbx,nby,nbz,nout,noutilf,ntinc,iwindt,nbldgpts,nbldgpt,ipntchk
         integer iturbfieldflag,inextgridflag,inloc,roofflag,nstpmax,nstpm
         integer nk,nf,iorg,jorg,jmin,imin,nfcur,nfmin,kk,ip,jp,nper,nupplevtemp
         integer iran,jran
         integer idxmax,idx1,idx2,iloop_drdx_spread
! new allocation mdw  1-22-2004
         integer, allocatable :: icellflag(:,:,:),ibldgpts(:),jbldgpts(:)
         integer format_flag,printflag,nsources,kbpufmax,isourceindex,isourcestart,isource2,icell1,icell2,jcell1,jcell2
         integer ibuoyflag,isourcetype,isiteflag,iindoorflag,idummy,mp !man 7/12/2005
         integer ixfo,jyfo,ifacr,jfacr,imrot,jmrot,ixur,iyur   !mdw 8/24/2007
         integer, allocatable :: noutface(:),ioutface(:,:),joutface(:,:),iacccount(:)
         integer, allocatable :: iaccflag(:,:,:),ifl(:,:),jfl(:,:),ifacefl(:,:,:),nsxylev(:)
         integer, allocatable :: ixpl(:),ixml(:),jypl(:),jyml(:),idxij(:,:)
         integer, allocatable :: icellflag_dense(:,:)
         real rcl,z0,rdummy,xf,yf,zf,xface,xblg,xsign,delta_ref,delta_ref_rm,delta_old,densegas_factor !man 7/12/2005
         real xiurjur,yiurjur,xiurjursum,yiurjursum,dxfo,dyfo,ximrot,yjmrot,dxur,dyur   !mdw 8/24/2007
         real x_b,y_b,z_b,ht_avg,u_top,u_bottom,u_left,u_right,nu_b,u_b
         real, allocatable ::  bin_read1(:),bin_read2(:),bin_read3(:) !man 7/15/2005
         integer, allocatable :: bin_read4(:) !man 7/15/2005
! mdw 6/22/2004 added particle radius index
         integer, allocatable :: icb(:),jcb(:)
         integer, allocatable :: bldtype(:),bldnum(:)
! new allocation mdw 1-22-2004
         integer, allocatable :: iface1n(:),iface2n(:),jface1n(:),jface2n(:),kface1n(:),kface2n(:)
         integer, allocatable :: istart(:),iend(:),jstart(:),jend(:),idxpm(:)
         real, allocatable :: dxm(:,:,:),dym(:,:,:),dzm(:,:,:),dosexy(:,:)
         real, allocatable :: dxp(:,:,:),dyp(:,:,:),dzp(:,:,:)
         real, allocatable :: depcx(:,:,:),depcy(:,:,:),depcz(:,:,:)
         real, allocatable :: depcxSP(:,:,:),depcySP(:,:,:),depczSP(:,:,:)
         real, allocatable :: u(:,:,:),v(:,:,:),w(:,:,:),hgt(:,:)
         real, allocatable :: dutotdzi(:,:,:),ustarz(:,:,:),sigwi(:,:,:)
         real, allocatable :: dutotdyi(:,:,:),dutotdxi(:,:,:)
         real, allocatable :: dutotdni(:,:,:),dutotdsi(:,:,:)
! dutotdni gives the gradient of the speed in the transverse direction
! while dutotdsi gives the gradient of the speed in the alongwind direction
         real, allocatable :: alph1ij(:,:,:),alph2ij(:,:,:),alph3ij(:,:,:)
         real, allocatable :: bet1ij(:,:,:),bet2ij(:,:,:),bet3ij(:,:,:)
         real, allocatable :: gam1ij(:,:,:),gam2ij(:,:,:),gam3ij(:,:,:)
! alph1 refers to coefficient of u''' in transformation to u from ''' system
! alph2 refers to coefficient of v''' in transformation to u from ''' system
! alph3 refers to coefficient of w''' in transformation to u from ''' system
! bet1 refers to coefficient of u''' in transformation to v from ''' system
! bet2 refers to coefficient of v''' in transformation to v from ''' system
! bet3 refers to coefficient of w''' in transformation to v from ''' system
! gam1 refers to coefficient of u''' in transformation to w from ''' system
! gam2 refers to coefficient of v''' in transformation to w from ''' system
! gam3 refers to coefficient of w''' in transformation to w from ''' system
         real, allocatable :: alphn1ij(:,:,:),alphn2ij(:,:,:),alphn3ij(:,:,:)
         real, allocatable :: betn1ij(:,:,:),betn2ij(:,:,:),betn3ij(:,:,:)
         real, allocatable :: gamn1ij(:,:,:),gamn2ij(:,:,:),gamn3ij(:,:,:)
         real, allocatable :: ufsqgi(:,:,:),vfsqgi(:,:,:),wfsqgi(:,:,:)
         real, allocatable :: ufvfgi(:,:,:),ufwfgi(:,:,:),vfwfgi(:,:,:)
         real, allocatable :: ufwfi(:,:,:),ufvfi(:,:,:),vfwfi(:,:,:)
         real, allocatable :: buoyintz(:,:,:),vbuoy(:,:,:),ubuoy(:,:,:),wbuoy(:,:,:),elevij(:,:)
         real, allocatable :: vlocz(:),vlocx(:),cumsum(:),vlocy(:)
         real, allocatable :: ubulkp(:),vbulkp(:),wbulkp(:)
         real, allocatable :: ztlev(:),uptemp(:),rhobkg(:),airtemp(:),wbulkz(:),zirhol(:,:),zirhoh(:,:)
         real, allocatable :: rplusij(:),rminusij(:)
         real, allocatable :: cosphidgi(:),sinphidgi(:),hnewi(:),rnewi(:),drdxi(:),xsceni(:),ysceni(:),ubari(:)
         real, allocatable :: rnewpi(:),rnewmi(:),drdxpi(:),drdxmi(:)
         real, allocatable :: upwpi(:,:,:),phib(:),phij(:,:,:),rhocell(:,:,:)
         real, allocatable :: ani(:,:,:),bni(:,:,:),cni(:,:,:)
!         real, allocatable :: asi(:,:,:),bsi(:,:,:),csi(:,:,:)
! alphn1 refers to coefficient of u in transformation to u'''
! alphn2 refers to coefficient of v in transformation to u'''
! alphn3 refers to coefficient of w in transformation to u'''
! betn1 refers to coefficient of u in transformation to v'''
! betn2 refers to coefficient of v in transformation to v'''
! betn3 refers to coefficient of w in transformation to v'''
! gamn1 refers to coefficient of u in transformation to w'''
! gamn2 refers to coefficient of v in transformation to w'''
! gamn3 refers to coefficient of w in transformation to w'''
         real, allocatable :: dose(:,:,:),utotktp(:,:),zcorf(:)
         real, allocatable :: uktop(:,:),vktop(:,:),wktop(:,:)
         real, allocatable :: dsigwdni(:,:,:),dsigudni(:,:,:),dsigvdni(:,:,:),dupwpdni(:,:,:)
         real, allocatable :: sigvi(:,:,:),epsi(:,:,:),sigui(:,:,:)
         real, allocatable :: con(:,:,:),xbc(:),ybc(:),zbc(:)
         real, allocatable :: conSP(:,:,:),doseSP(:,:,:)
         real, allocatable :: xi(:),yi(:),zi(:),dsigwdzi(:,:,:)
         real dx,dy,dz,dwall,dx2,delym,dxx,dyy,dzz,dxc,dyc,dzc,ustarloc,dwallg
         real alph1,alph2,alph3,bet1,bet2,bet3,gam1,gam2,gam3,h,slopex,slopey,vslope,sinslope
         real deltatot,deltaxpart,deltaypart,deltatotref,deltatotrefp,delx,distance_ref         
         real alphn1,alphn2,alphn3,betn1,betn2,betn3,gamn1,gamn2,gamn3,wstar
         real u3psq,v3psq,w3psq,upwp,upsq,vpsq,wpsq,ufsq,vfsq,ufwf,ufvf,vfwf
         real ufsqb,vfsqb,wfsqb,ub,vb,wb,upb,vpb,tls,tstart,tinc,twind,tupdate
         real sinphiw,cosphiw,omenum,omeden,tanomeg,omeg,cosomeg,sinomeg,dutotdn
         real cospsi,sinpsi,cosphi,sinphi,upvp,vpwp,vpsqb,wfsq,upsqg,wpsqg,upwpg
         real upvpg,vpsqg,vpwpg,omeg1,omeg2,sinomeg1,sinomeg2,cosomeg1,cosomeg2
         real dutotdn1,dutotdn2,taumin,deltilf,dnear,xloc,decay,pp,refhgt
         real rhoair,volsource,massliq,tsourcref,zrho,phim,psim
         real dxi,dyi,xsdencells,xslev,zlev,dmin,dnew,buoyref,ristar
         real rhotref,tair,ustarp,ustarb,rib,ubuoyt,vbuoyt,ubuoyl,vbuoyl
         real buoypx,buoymx,buoypy,buoymy,rhocellm,buoymin,xbuoyc,xpp,ybuoyc,ypp
         real delppert,wbuoyt,voutbar,dxxb,dyyb,buoymz,ubuoyp,vbuoyp,sbuoyp
         real denfire,denmax,rclltg,rclloc,forcez,firevol,forcezpp,zirhoht,zirholt,airtempzp
         real xsubsw,xsubne,ysubsw,ysubne,tz,utott,cusq,cvsq,cwsq,ctau13,zbrac
         real zpmin,zpmax,xpmin,xpmax,ypmin,ypmax,wt,xlt,xi1,yj1,zk1,terfac
         real zkfac,zk,xpc,ypc,zpc,con1,con2,con3,con4,totsrcmass,deltax,deltay
         real deltaz,xbl,xbu,ybl,ybu,zbl,zbu,dbx,dby,dbz,vol,phi,phit
         real dpart,vs,sc,taud,diff,press,temp,rhop,xnu
         real xcbp3,xcbm3,ycbp3,ycbm3,dycbm3,dycbp3,dxcbp3,dxcbm3,ds,sup,sdown
         real stin,st,urefz,xcd,xcdl,xcu,xcul,ycd,ycdl,ycu,ycul,cosl,sinl,cosv
         real vrel,vrel1,vrel2,xrel
         real delri,uxsum,vysum,cellsum,uistarsum,ycellcn,slopeysum,slopexsum,utotsum ,utotdg
         real cosphidgw,sinphidgw,ubarw,ubars,cosbeta,sinbeta,cosalph,sinalph,dubardx,duperpdx,dhdx,duperp
         real cosphidgn,sinphidgn,ubarwtop,ustarneutral
         real xcellcl,xcellcn,ycellcl,rcellcl,rcellcn,addarea
         real effarea,curarea,ustars,ristar0,utent,delvol,delvol_slope,cd
         real rold,hold,volold,ubarold
         real wbulkzs,sumwbulkz,wbulkzbar
         integer kcellst,iconarea,icells,jcells,istrt,ifin
         integer jstrt,jfin,icell,jcell,iscena,jscena
! mdw 3-08-2004 add new variables for highest vertical gradient
! mdw 10-04-2004 made dpart real
! mdw 7-08-2005 add new variables for deviation from log-law winds
         real dutotdza,dutotdzc,dutotdzm,dutotdzp,dutotl,dutot,dutotu
         real dutotdxa,dutotdxc,dutotdxm,dutotdxp
         real dutotdya,dutotdyc,dutotdym,dutotdyp
         real, allocatable :: xp(:),yp(:),zp(:),zpp(:),tstrt(:),xnear(:),ynear(:),znear(:),temp_part(:)
         real, allocatable ::pufmasfr(:),znorm(:),psigx(:),psigy(:),psigz(:)
         real, allocatable ::pufcummas(:),vsi(:),taudi(:),diffi(:),sci(:),dpm(:),dpm_bin(:),weightr(:)
         real cldhgt,hemass
         integer, allocatable ::isub(:)
! mdw 3-16-2004 add particle radius array and particle distribution weight
         real, allocatable :: uran(:),vran(:),wran(:),weight(:),weightd(:)
         real, allocatable :: xps(:),yps(:),zps(:),dwallp(:)
!         real, allocatable :: xpr(:),ypr(:),zpr(:),trel(:),rpr(:)
         real xpr,ypr,zpr,rpr ! MAN 7/25/2006
         real kkar,knlc,lam11,lam22,lam33,lam13,lam12,uint,vint,wint
         real utotu,utotl,utot,cost,sint,dij,dxii,dxmax,utotp,ualoft,valoft
         real dijm,delutotp,utotm,delutotm,ustar,ysign,den1sr,concstrt
         real dyv,delta,dt,tre,dwxp,dwxm,dwyp,dwym,dfzm,drzp
         real sigu,sigv,sigw,dsigwdz,u1,v1,w1,tau11,tau22,depfrac
         real tau33,tau13,cosu,sinv,dyct,coeps,du,dv,dw,tau12,u1p,v1p
         real dudet,duran,dvdet,dvran,dwdet,dwran,outp,outc,tout,toutc
         real dxxw,dyyw,dzzw,xr,yr,wr,t,utotparp,utotparm,dparp
         real xparp,yparp,dutotdx,dutotdy,dutotdz,dxyz,du2x,du2y,du2z
         real, allocatable :: elz(:,:,:),eleff(:,:,:)
         real dwallx,dwally,delzm,delxm,du2,theta,dur,durrel
!         real ustarr
         real, allocatable :: Ht(:),Wti(:),Lt(:),ustarij(:,:,:),elzg(:,:,:)
         real, allocatable :: xfo(:),yfo(:),zfo(:),weff(:),leff(:),xft(:)
         real, allocatable :: atten(:),lr(:),sx(:),sy(:),lfr(:),gamma(:)
         real, allocatable :: xcb(:),ycb(:),uref(:,:)
         real, allocatable :: urefu(:,:),urefv(:,:),urefw(:,:),ustargz(:,:,:)
         real, allocatable :: ustarg(:,:,:),deluc(:,:)
! new allocations mdw 1-22-2004
         real, allocatable :: chioa(:,:),chiosum(:,:),chiin(:,:),dosin(:,:),dosinSP(:,:)
         real, allocatable :: chiinSP(:,:),chioaSP(:,:),chiosumSP(:,:),chiex(:,:),chiexSP(:,:)
         real, allocatable :: aminusx(:),aplusx(:),aminusy(:),aplusy(:),aroof(:),afloor(:)
         real, allocatable :: utotmax(:),utotcl1(:)
         integer, allocatable :: nonlocal_option(:,:,:),local_option(:,:,:),turb_options(:,:,:)

! MAN 02/23/2007 modified infiltration scheme
         real toact,piact,rhoact
         real, allocatable :: tau(:),racact(:),racref(:),tiref(:),toref(:),piref(:),urefb(:)
         real, allocatable :: tiact(:),rhoref(:),racsupply(:),racintake(:),filtration(:),area2volume(:)
         real, allocatable :: ainf(:),binf(:),abinf(:),binfex(:)
! end MAN 02/23/2007
         real psref,psact,pwref,pwact,aref,aact,touti,xdum,ydum,zdum,xoff,yoff
         real delutz,duy
         real up,vp,wp,dsum,usum,vsum,wsum
         real m_roof
         integer time_idx,idepositionflag,iwallfieldflag,ilctflag,ieradflag !MAN 7/24/2006
         integer nthreshold,isourcetypeflag !MAN 7/24/2006
         real maxdiameter,mindiameter,depvel,conout !MAN 7/24/2006
         real, allocatable :: dosethreshold(:) !MAN 7/24/2006
! MAN 7/31/2006 multiple source variables
         real particle_mass,total_mass
         integer isource,inode,nparticles,nsourcenodes,npartsize,nsegs,iturbtypeflag
         real, allocatable :: source_x(:),source_y(:),source_z(:),source_l(:),source_w(:)
         real, allocatable :: source_h(:),source_r(:),source_hemass(:)
         real, allocatable :: source_strength(:),source_density(:),source_start(:),source_dur(:)
         real, allocatable :: source_seglength(:),source_nodestart(:),source_nodestop(:)
         real, allocatable :: source_nodespeed(:),source_nodeaccel(:),source_nodesegpos(:)
         real, allocatable :: source_weight(:),source_num(:)
         integer, allocatable :: source_geometry(:),source_str_units(:),source_release(:),nlinenodes(:)
         integer, allocatable :: source_startidx(:),source_stopidx(:),source_prr(:)
         integer ipuf,ippr,itemp,nsteps,nptot,it,nstps,idum,npuff,ipufr,imatnum
         real pi,tboil,tfire,rfire
         integer ithreshold,nx1,ny1,nz1,nwinds,klow,nxg,nyg,ipr,ib,jb,kb,ibuild
         integer ioutin,iface1,iface2,ifacep,jface1,jface2,kface1,kface2,iface,jface,kface
         integer ixg,iyg,jpr,kpr,nps,npf,ivrelch,imm,jmm,kmm,isign,jsign,ksign
         integer mmpt,mmm,ntb,immm,jmmm,noutbldg,iprint,jprint,kprint,nface
         integer nxw,nyw,nzw,nzwout,kw,kwp,kmp,jw,jwp,jmp,iw,iwp,imp,icellflg
         real bvol,dxg,dyg,xgl,ygl,xpma,ypma,zpma,dsigvdn,dsigudn,wstar3,dtau13dn,siguo,sigvo,sigwo
         real dsigwdn,dutotds,utestm,vtestm,wtestm,x,w1p,taup,rb,ra,vd,frac1
         real tfirez,denfirez,toutin,zprint,gridind
         real zw,zwp,yw,ywp,xw,xwp,b,ang1,ang2,srcmass
         integer nreleases,irelease
! end MAN 7/31/2006
         real, allocatable :: volpart(:),volpartp(:),xsum(:),ysum(:),xsqsum(:),ysqsum(:)
         real, allocatable :: activepart(:),volpart0(:),voldelrho(:),volzrho(:),volzdrhobdz(:),activepart0(:)
         real, allocatable :: xic(:),yic(:),xsqic(:),ysqic(:),varx(:),vary(:),sigxb(:),sigyb(:),activepartp(:)
         real, allocatable :: lim1(:),alim(:),wbulkzsum(:),wbulkza(:),conl(:),vlsum(:)
         real, allocatable :: rl(:),rpl(:),rpm(:),ubulkpx(:),vbulkpy(:),ustarr(:,:),hnewij(:,:),dhdtij(:,:)
         real reffa,hnewa,drdta,volnewa,ycompdevmax,ycompdevmin,ymaxslice,yminslice,drdxplus,drdxminus
         real xscena,yscena,rnewtotalp,rnewtotalm
         real, allocatable :: hnew(:),hnewp(:),rnew(:),rnewp(:),reff(:),reffp(:),drdt(:),drdtp(:),rhocent(:)
         real, allocatable :: xscen(:),xscenp(:),yscen(:),yscenp(:),volnew(:),vo10(:),volnewp(:),vol0(:)
         real kent,drag_local
         real, allocatable :: uranbp(:),vranbp(:),wranbp(:),vslopex(:),vslopey(:),dhdt(:),cosphidg(:),sinphidg(:),ubar(:)
         real sigzb,varz,delrho0,dzk,areatp,conmax,conmaxp,contot,contotp
         real activepartt,deltaup,rcessr,rlocran,sinrptcen,cosrptcen,rptcen
         real zbcen,vout,rnewa,ubulk,vbulk,drloc,xcessr,ylocran,xrany,xranx
         real xran,yran,yranw,xranw,yloc,crat,rhowt0,rhowt1,zbmax,zbmaxp,zbminp,zbmin
         real uent,vent,went,dzbpuf,offcntdis,entrain_factor,volchange,vlsumt,areat,vlsumtp
         real fcoef,vref,tscale,sigwb,sigub,uranb,vranb,wranb,fac,sigwbn,sigvbn,sigubn,vdef,delrpuf
         real delxcen,delycen,delzcen,zic,dhdta,sigvb,drdx
         real rempty,zranw,zran,zranx,zsum,zsqsum,zbpufmin,zbpufmax,zbpuftop,zbpufbot
         real zpuftop,zpufbot,zkmin,zkmax,zsqic,vup
         real utoth,ratzran,frac2,rcir,ziloc,vloc
         integer nact,kzmax,kzmin,iflagden,kz,kbpuf,number_cutoff,idx,idxp
         integer ipc,jpc,kpc
         integer lognormal_flag,num_distributions,distribution_type_flag
         real, allocatable ::distribution_mean(:),distribution_std(:),distribution_mass_mean(:)
         real, allocatable ::distribution_mass_fraction(:),distribution_num_fraction(:)
         real, allocatable ::distribution_cum_fraction(:)
!MAN 1/25/2007 Particle splitting
         real, allocatable ::pmass(:)
         real oldcon
         integer ninactiveparticles,nparticlesfactor,partsplitfactor
!MAN 02/08/2007 Added Taylor microscale lower limit to subtimesteps.
         real taylor_microscale,tls_min !,taylor_max
         integer taylor_flag
!MAN 1/22/2007 change in surface deposition binary file format
         integer nsurfdep
!MAN 04/02/2008 evaporation met data arrays
         real, allocatable ::temperature(:,:,:),vaporfield(:,:,:),dpm_min(:)
         integer ibioslurryflag,itwophaseflag
         real agent_solids_concentration,agent_eff_solid_density
         real agent_droplet_fraction,agent_mfp,agent_stick_eff
         real agent_mw_vapor,agent_mw_liquid,atm_pressure,thermal_accomodation
         real atomic_volume_agent,Pcritical,Tcritical,Tboiling,liquidDensity,surfaceTensionCoef
         
         integer inverseflag,ibuoyflag_inner
! mdw added xdum,ydum, zdum to definition 3-08-2004 also added proper
! treatment of zfo
! man 10/10/2007
         real z0coeff
! man 11/27/2007 
         integer iqperror,iqperrormax,iranvelerror
! man 04/14/2008 
         real massMedianDiameter
      end module variables
! TMB
! use subroutine/return and comment program/end for matlab compilable version
!     subroutine qwicplumeglo
      program qwicplumeglonew
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This version only creates and moves  particles with the winds from
! qwic-urb; reflection, and inhomogeneous turbulence (y and z), are  dealt with
! mean cell winds are used with interpolation; concentrations are calculated
! u,v,w are cell-centered values
!  - the final particle field is written out to partfield.dat(=10)
! Major revisions made 7/11/2003 by mdw - (1) corrects error in turbulence for
! part of front eddy (2) uses wake parameters from screendat.out (lr, lfx, and
! lfy) also uses effective building widths and lengths (3) corrects error in
! concentration estimates for the first time step  (4) improved descriptions
! for turbulence around complex buildings (5)improved treatment of random
! processes                                    note that deposition using
! deposition velocities that differ in horizontal and vertical was added
! earlier
! additional revisions by MDW 7/21/2003 corrected deposition calculation - was
! normalized and weight was left out so it counted previously deposited particles
! dose version also added 7/21/2003
! vertical non-local mixing added 7/23/2003
! corrected reflection, made a combined code, and added a concentration
! calculation start time 7/30/2003
! added options.dat file to control output
! changed epsilon to be calculated from local gradients only
! 08/12/2003 corrected expression for phit (thanks Bali & Eric)
!    Note:
!    For compiling the code set the following options under Project-Settings:
!    - 2 under fortran
!    --- Floating Point = 0
!    --- Fortran Language = 132
!    - 1 Under Link
!    --- Output Reserve 512000000 (as per Petra)
!    this version of the code is based on a triply rotated coordinate system
!    first the coordinate system is rotated to align with the mean horizontal
!    wind; then the coordinate system is rotated once more to align with the
!    local wind and finally the coordinate system is rotate once more so that
!    the new z-axis extends in the direction of the greatest increase in the
!    mean wind. In the new system, V=0, W=0, and the largest changes in U are
!    in the positive z direction
!    other changes include a flag to turn non-local mixing on (nloc=1 in options.dat,
!    a restriction that the local time step be less than or equal to one quarter
!    of the local velocity decorrelation time scale tls. There is also restriction
!    so that the random speed is less than 3 sigma u (in the system aligned with
!    the mean wind
!    this version dates to 9/11/2003
!    corrected equation for u3psq in subroutine rotufsq 9/19/2003
!    consistent with new input.dat 9/29/2003
!    corrected equations for rotation in detang and for gradients
!    and tangent omega calculation; also corrected epsi(im,jm,km)
!    added missing drift terms to dudet  10/30/2003 - also corrected
!    icelt,jcelt in transverse points in wake and eddy
!    correction made for ustarg and ustargz; test for minimum velocity deficit
!    in both the vertical and horizontal for the non-local mixing
!    added into the wall non-local mixing to treat front and back pockets
!    corrected differencing in wake to avoid holes for off cardinal winds
!    corrected geometry for winds from other directions - 11/05/2003
!    also corrected vpsgq to vpsqg in vertical non-local mixing
! 11/06/2003 corrected bldg parameters etc to permit the code to run with
!    dx not equal to one
!    allocate nbx,nby,nbz multiply random force term by sqrt(dt)
!    changes made on 11/10/2003
!    change to local eps; correct gradient of sigma w near ground
!    also correct new du drift terms (2/1.3) goes to (4/1.3) - 11/17/2003
!    added v-drift term 11/17/2003
!    corrected hgt and restored gridded version of eps  - 11/18/2003
!    corrected v and w drift terms and corrected read statements for screenout.dat
! 11/24/ 2003 mdw
! 11/25/ 2003 mdw corrected premature zeroing of array dose
! 11/25 2003 mdw close(25)
! 12/01/2003 changed dx to dy in expression for depcy has no effect at present
!    since dy=dx=dz also added deallocate statements for ufsqgi through csi
! 12/02/2003added deallocate statements for zcorf,ufvfgi, ufwfgi, and vfwfgi
! 12/03/2003 cleaned up some aspects of the non-local mixing description and
!    corrected an error elz along the domain boundaries
! 12/04/2003 eliminate duplicate write statements to unit 13 and put in close statement
! 12/09/2003 eliminated useless variables elx, ely, ustarx, ustary, icellflag1
!    icellflag2, and dely also deallocated several gridded arrays before depcx
!    depcy, depcz, con, and dose are allocated
! 01/06/2004 renamed qwicplumecom.inp to quicplume.inp
! 01/08/2004 corrected checks on vrel and used log law for near walls
!    previously version used log law for rooftops and the ground
!    new version including U of U changes and building infiltration 1/22/2004
!    adding buildout.dat in place of screenout.dat 2-03-2004
!    changing check for internal deposition velocity calculation to vd=-99.
!    225/2004
! 03/22/2004 pm changed read input.dat to read new quicurb, set area of building
!    faces to a minimum of 1.e-08 meter squared to avoid divide by zero in indoor
!    calculations for hidden building faces, and forced ustarg to be greater
!    than calculated ustarg for replacement
! 03/31/2004 changed the format for writes to confield.dat and dosefield.dat
!    to permit larger values of x, y, and z
! 04-15-2004 added statements to write out closest position to a wall
!    during deposition for each particle also changed all real 4's to real 8's
! 04-16-2004 mdw added corrections for when ktop or ktop+3 are above nz-1
! 04-19-2004 mdw changed to write out last deposition positions
! 05-14-2004 mdw corrected expression for crosswind integration in bldg eddies
! 06-??-2004 mdw to 11-23-2004 made several changes including:
!    (1) improved triple reflection scheme,
!    (2) adjustment for positions within machine signficance of wall,
!    (3) changed expression for z dependence of dissipation to reflect near neutral stability rather than slightly stable stability (did have l=100),
!    (4) used log-law for points in cells closest to wall for interpolation of winds
!    (5) prevented distance from wall for interpolation to get less than z0+.01*z0 to keep log-law from blowing up
! 11-30-2004 revised check to read new winds if t>tupdate to skip read if the check
!    occurs on the last time step; also corrected the big the produced failure for
!    wind directions of 247.5, also put in system for transport under roofs
! 09-09-05 (Matt and balli) Changed near wall local turbulence parametrization: calculate average building height, average wind speed on the perimeter of the domain at  average building height.
!    These are used to calculate the distance from the walls to include the horizontal length scales in the near wall local turbulence parameterizations
!    Multiplied the vortex local mixing turbulence parameterizatin length scale by a factor 0.5 (to match joint urban data)
!    Changed the 'jfix' intrinsic fuction with 'nint' while calculating the building center and the points for the reference velocity in turbinit subroutine.
!    Rotated the non-local mixing horizontal reference velocity into the mean wind direction of each building.
!    Added the QP_turboptions.bin to store which turbulence schemes have been applied in the domain.
!    turboptions : 1- near wall local turbulence scheme
!                  2- vortex local turbulence scheme
!                  3- local uniform flow turbulence scheme
!                  4- Nonlocal vertical mixing
!                  5- Nonlocal Horizontal mixing.
! 11-21-2005 mdw  modified coding for gradients in the sigmas to address particles in
!    3-walled cells
! 07-05-2006 mdw changed form for xcu, ycu, etc to length /width + .5 dx or .5 dy
!    and icelx jcelx to nint(xcxx/dx)+1 instead of nint(xcxx/dx) and similarly
!    nints of ycxx/dy
! 08-16-2006 mdw changed the interpolation near the wall to log(z/z0)instead of
!    log((z+z0)/z0) and changed the reflection from the surface to z0 outside
!    of the surface - reflection is now in subroutine reflect also added stops
!    1-6 for cases where a particle ends up within z0 of a surface
! 08-31-2006 MAN
!    (1) added new cylindrical volume or area, rectangular volume or area, multisegment line, and moving point sources
!    (2) modified code so that all z0's are the same.
!    (3) added a minimum limit to the particle size thresholds and made them variable
!    (4) modified code so that there are no implicit variables
! 10-30-2006 MAN added latest buoyant and dense gas parameterizations from 7-31-2006a version of QUIC 3.95
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         use variables   !make data from module "variables" visible
         !use Evaporation
         implicit none
         character*40 head1,head2
! mdw 2-03-2004 revised character statements for reading in buildout.dat
         real ran2,elev_params(3)
         integer k1,k2,idenseflag
!         type(StructuredGrid3d)::windgrid
         pi=4.*atan(1.)
         t=0.
         kent=.2
         number_cutoff=3
         idenseflag=1
         time_idx=1
         itwophaseflag=0
         ibuoyflag_inner=0
         z0coeff=1.01
         iqperror=0
         iqperrormax=50000 
! open required input files
         open(unit=5,file="QP_params.inp",status="old")
         open(unit=7,file="QP_materials.inp",status="old")
         open(unit=8,file="QP_indoor.inp",status="old")
         open(unit=34,file='QU_simparams.inp',status='old')
         open(unit=35,file='QU_buildings.inp',status='old')
         open(unit=37,file='QP_source.inp',status='old')
         open(unit=44,file='QP_buildout.inp',status='old')
         open(unit=62,file="QP_fileoptions.inp",status="old")
         open(unit=77,file='QP_error',status='unknown')  
         
!         open(unit=88,file='QP_dt.bin',form='unformatted',status='unknown')
! MAN 08/23/2006 Read input files
         call readQUsimparams
!         call readQUbuildings
         call readbldgout
         call readQPfileoptions
         
         allocate(icellflag(nx-1,ny-1,nz-1))
         allocate(icb(inumbuild),jcb(inumbuild))
         allocate(dxm(nx,ny,nz),dym(nx,ny,nz),dzm(nx,ny,nz))
         allocate(dxp(nx,ny,nz),dyp(nx,ny,nz),dzp(nx,ny,nz))
         allocate(u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz),hgt(nx,ny))
         allocate(dutotdzi(nx,ny,nz),ustarz(nx-1,ny-1,nz-1),sigwi(nx-1,ny-1,nz))
         allocate(dutotdyi(nx,ny,nz),dutotdxi(nx,ny,nz),epsi(nx,ny,nz))
         allocate(dutotdni(nx,ny,nz),dutotdsi(nx,ny,nz))
         allocate(alph1ij(nx,ny,nz),alph2ij(nx,ny,nz),alph3ij(nx,ny,nz))
         allocate(bet1ij(nx,ny,nz),bet2ij(nx,ny,nz),bet3ij(nx,ny,nz))
         allocate(gam1ij(nx,ny,nz),gam2ij(nx,ny,nz),gam3ij(nx,ny,nz))
         allocate(alphn1ij(nx,ny,nz),alphn2ij(nx,ny,nz),alphn3ij(nx,ny,nz))
         allocate(betn1ij(nx,ny,nz),betn2ij(nx,ny,nz),betn3ij(nx,ny,nz))
         allocate(gamn1ij(nx,ny,nz),gamn2ij(nx,ny,nz),gamn3ij(nx,ny,nz))
         allocate(upwpi(nx,ny,nz))
         allocate(ani(nx,ny,nz),bni(nx,ny,nz),cni(nx,ny,nz))
         allocate(icellflag_dense(nx,ny))

          print*,"I am here23438888"
         if(iturbtypeflag.eq.1)then
            allocate(nonlocal_option(nx-1,ny-1,nz-1),local_option(nx-1,ny-1,nz-1),turb_options(nx-1,ny-1,nz-1))
            nonlocal_option(:,:,:)=0
            local_option(:,:,:)=0
            turb_options(:,:,:)=0
         endif
          print*,"I am here 21111111"
        
         icellflag_dense(:,:)=0
         ani(:,:,:)=0.
         bni(:,:,:)=0.
         cni(:,:,:)=0.
         alph1ij(:,:,:)=0.
         alph2ij(:,:,:)=0.
         alph3ij(:,:,:)=0.
         bet1ij(:,:,:)=0.
         bet2ij(:,:,:)=0.
         bet3ij(:,:,:)=0.
         gam1ij(:,:,:)=0.
         gam2ij(:,:,:)=0.
         gam3ij(:,:,:)=0.
         alphn1ij(:,:,:)=0.
         alphn2ij(:,:,:)=0.
         alphn3ij(:,:,:)=0.
         betn1ij(:,:,:)=0.
         betn2ij(:,:,:)=0.
         betn3ij(:,:,:)=0.
         gamn1ij(:,:,:)=0.
         gamn2ij(:,:,:)=0.
         gamn3ij(:,:,:)=0.
         upwpi(:,:,:)=0.
         epsi(:,:,:)=0.
         print*,"ssssssss";
         allocate(sigvi(nx-1,ny-1,nz),dupwpdni(nx-1,ny-1,nz),sigui(nx-1,ny-1,nz))
         allocate(xi(nx),yi(ny),zi(nz),dsigwdzi(nx-1,ny-1,nz-1),zcorf(nz))
         allocate(elz(nx,ny,nz),eleff(nx,ny,nz))
         allocate(xcb(inumbuild),ycb(inumbuild),uref(inumbuild,nz))
         allocate(urefu(inumbuild,nz),urefv(inumbuild,nz),urefw(inumbuild,nz))
         allocate(deluc(inumbuild,nz))
         allocate(ustarg(nx,ny,nz),ustargz(nx,ny,nz),elzg(nx,ny,nz))
         allocate(istart(inumbuild),iend(inumbuild),jstart(inumbuild),jend(inumbuild))
         allocate(utotmax(inumbuild),utotcl1(inumbuild))
         
         zpmin=real(nz-1)*dz
         zpmax=0.
         xpmin=real(nx-1)*dx
         xpmax=0.
         ypmin=real(ny-1)*dy
         ypmax=0.
         
! define cell center locations
         do k=1,nz
            zi(k)=-.5*dz+dz*real(k-1)
         enddo
         do j=1,ny
            yi(j)=.5*dy+dy*real(j-1)
         enddo
         do i=1,nx
            xi(i)=.5*dx+dx*real(i-1)
         enddo
         
         call readQPmaterials
         call readQPparams
         if(ieradflag.eq.1)then
            call readQPerad
         else
            call readQPparticlesize
         endif
         if(ibioslurryflag .eq. 1 .or. itwophaseflag .eq. 1)then
            open(unit=88,file='QP_evaporation.inp',status='old')
            call readQPevaporation
         endif
         do ippr=npartsize,2,-1
            dpm_bin(ippr)=.5*(dpm_bin(ippr)+dpm_bin(ippr-1))
         enddo
         dpm_bin(1)=.5*dpm_bin(1)
         allocate(weightr(npartsize),vsi(npartsize))
         allocate(taudi(npartsize),diffi(npartsize),sci(npartsize))

         call readQPsource
         if(ibuoyflag .eq. -1 .and. inextgridflag .eq. 2 &
               .and. source_release(1) .gt. 1)then
            ibuoyflag = 0
            ibuoyflag_inner = -1
         endif
         if(ibuoyflag.ne.0)then
            call readQPbuoy
         endif
         
! end MAN 08/23/2006
         if(depvel.ne.0.)then
            idepositionflag=1    !writes out all deposition
         else
            idepositionflag=0
         endif
         printflag=1
! open required input data files
         if(format_flag.eq.1)then
            open(unit=27,file="QU_velocity.dat",status="old")
            open(unit=30,file="QU_celltype.dat",status="old") ! material or atmos. cell
         else
            open(unit=50,file="QU_velocity.bin",form="unformatted",status="old")
            open(unit=64,file="QU_celltype.bin",form="unformatted",status="old") ! material or atmos. cell
         endif
! open required output data files
!erp 3/2/05 binary file writes
         if(isiteflag.eq.1)then ! only used for sensor siting.
            open(unit=16,file="QSS_sourcefield.dat",status="unknown")
            if(format_flag.eq.1 .or. format_flag.eq.3)then
               open(unit=21,file="QSS_dosage.dat",status="unknown")
            endif
            if(format_flag.eq.2 .or. format_flag.eq.3)then
               open(unit=22,file="QSS_dosage.bin",form="unformatted",status="unknown")
            endif
            iwallfieldflag=0
            ilctflag=0
            iturbtypeflag=0
            inextgridflag=0
            iturbfieldflag=0
         else
            if(format_flag.eq.2 .or. format_flag.eq.3)then  !erp 3/2/05
               open(unit=46,file="QP_confield.bin",form="unformatted",status="unknown")
               open(unit=54,file="QP_partfield.bin",form="unformatted",status="unknown")
               open(unit=55,file="QP_dosage.bin",form="unformatted",status="unknown")
               if(idepositionflag.eq.1)open(unit=58,file="QP_depct.bin",form="unformatted",status="unknown")
               if(iturbfieldflag.eq.1) open(unit=63,file="QP_turbfield.bin",form="unformatted",status="unknown")
            endif
            if(format_flag.eq.1 .or. format_flag.eq.3)then  !erp 3/2/05
               open(unit=10,file="QP_partfield.dat",status="unknown")
               open(unit=11,file="QP_confield.dat",status="unknown")  !erp
               open(unit=15,file="QP_dosage.dat",status="unknown") !erp
               if(idepositionflag.eq.1)open(unit=19,file="QP_depct.dat",status="unknown")
               if(iturbfieldflag.eq.1) open(unit=13,file="QP_turbfield.dat",status="unknown")
            endif
            if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
               if(format_flag.eq.2 .or. format_flag.eq.3)then !erp 3/2/05
                  open(unit=47,file="QP_confieldSP.bin",form="unformatted",status="unknown")
                  open(unit=56,file="QP_dosageSP.bin",form="unformatted",status="unknown")
                  if(idepositionflag.eq.1)open(unit=59,file="QP_depctSP.bin",form="unformatted",status="unknown")
               endif
               if(format_flag.eq.1 .or. format_flag.eq.3)then !erp 3/2/05
                  open(unit=41,file="QP_confieldSP.dat",status="unknown") !erp
                  open(unit=45,file="QP_dosageSP.dat",status="unknown") !erp
                  if(idepositionflag.eq.1)open(unit=49,file="QP_depctSP.dat",status="unknown") !erp
               endif
            endif
         endif
! open optional output data files
         if(iwallfieldflag .eq. 1) open(unit=12,file="QP_wallfield.dat",status="unknown")
         if(ilctflag .eq. 1) open(unit=14,file="QP_lctout.dat",status="unknown")
         if(iturbtypeflag .eq. 1)open(unit=101,file="QP_turbtype.bin",form='unformatted',status="unknown")
! MAN 9/13/2005 add binary nextgrid data file
         if(isourcetypeflag .eq. 8)then
            open(unit=25,file="QP_exfiltration.inp",status="old")
         endif
         if(inextgridflag .ge. 1)then
            if(format_flag.eq.2 .or. format_flag.eq.3)then  !erp 3/2/05
               open(unit=32,file="QP_nextgrid.bin",form="unformatted",status="unknown")
            endif
            if(format_flag.eq.1 .or. format_flag.eq.3)then  !erp 3/2/05
               open(unit=31,file="QP_nextgrid.dat",status="unknown")
            endif
         endif
         if(inextgridflag .eq. 2)then
            if(format_flag.eq.2 .or. format_flag.eq.3)then  !erp 3/2/05
               open(unit=38,file="QP_nextgrid.inp",form="unformatted",status="unknown")
            else
               open(unit=33,file="QP_nextgrid.out",status="unknown")
            endif
         endif
! end MAN 9/13/2005
         
!erp 3/02/05 end format flag add

! MAN 9/21/2005 roof top mixing length fix
         m_roof=6
         mrel=100
         nout=20
! we read the starting time for concentration calculations
         tout=-delta !MAN 01/29/2007 tout=0.
         toutc=0.
         touti=-delta !MAN 01/29/2007 touti=0.
         vs=0.
         cusq=(2.5)**2
         cvsq=(2.)**2
         cwsq=(1.3)**2
         ctau13=1.
!mdw 10-05-2005; the preceding constants are multipliers of ustar squared for
! sigma u squared, sigmav squared, sigmaw squared, and tau13
         ufsqb=0.
         vfsqb=0.
         wfsqb=0.
         vpwpg=0.
         vpsqg=0.
         vpsqg=0.
         concstrt=0.
         vup=0.
         deltaup=0.

! mdw 4-16-2004 added proper treatment of zfo
         ustarg(:,:,:)=0.
         ustargz(:,:,:)=0.
!mdw 2-02-2004 begin buildout.dat replacement for screenout.dat
    2    continue
    8    format(a34)
    3    format(i12,a6,f16.2)
    11   format(a6,f25.15)
! mdw end reading buildout.dat instead of screenout.dat
     6   continue
         if(format_flag.eq.1) read(27,10)head1 ! read uoutmat.dat
      10 format(a40)
         kkar=.4
         knlc=.113
         if(format_flag.eq.1)then
            read(27,*)nwinds
            read(27,10)head2
            read(27,*)iwindt
            read(27,*)twind
         endif
         call readwind
         if(format_flag.ne.1)then
            allocate(bin_read4((nx-1)*(ny-1)*(nz-1)))
            read(64)bin_read4
            idummy=1
            do k=1,nz-1
               do j=1,ny-1
                  do i=1,nx-1
                     icellflag(i,j,k)=bin_read4(idummy)
                     idummy=idummy+1
                  enddo
               enddo
            enddo
            deallocate(bin_read4)
         else
            do k=1,nz-1
               do j=1,ny-1
                  do i=1,nx-1
                     read(30,*)xdum,ydum,zdum,icellflag(i,j,k) ! tell whether you are in a building or atmos over the grid
!                     icellflag(i,j,k)=min(icellflag(i,j,k),1)
                     if(k.eq.1) then
                        icellflag(i,j,k)=0.
                        cycle
                     endif
                  enddo
               enddo
            enddo
         endif
! we next read the site flag; if isiteflag=0 we have a normal application
! while if isiteflag=1 we are doing a sensor siting application
         if(isiteflag.eq.0)then
            terfac=0.
            nxg=1
            nyg=1
         else
            terfac=1.
         endif
         xloc=real(inloc)

! setup input for release of particles

         if(inextgridflag .eq. 2)then
            open(unit=999,file="QP_nptotfinal.inp",form="unformatted",status="old")
            read(999)nparticles
            close(999)
         endif
         mpart=nparticlesfactor*nparticles+1
         allocate(xp(mpart),yp(mpart),zp(mpart),tstrt(mpart),weight(mpart),weightd(mpart),dpm(mpart))
         allocate(uran(mpart),vran(mpart),wran(mpart),uranbp(mpart),vranbp(mpart),wranbp(mpart))
         allocate(xps(mpart),yps(mpart),zps(mpart),source_num(mpart),dwallp(mpart))
         allocate(xnear(mpart),ynear(mpart),znear(mpart),pmass(mpart))
         if(ibuoyflag.gt.0)then
            allocate(wbulkz(mpart),zpp(mpart),temp_part(mpart))
            wbulkz(:)=0.
            zpp(:)=0.
            temp_part(:)=0.
         endif
         if(ibuoyflag.lt.0)then
            allocate(elevij(nx+1,ny+1),idxpm(mpart))
            idxpm(:)=1
            open(unit=65,file="QP_elevation.inp",form="unformatted",status="old")
            open(unit=67,file="QP_densegas.dat",form="formatted",status="unknown")
            write(67,68611)
68611       format('time, h, r, reff, xcenter, ycenter, density at center')
            read(65)elev_params
            if(elev_params(1) .lt. 0.5)then
               elevij(:,:)=0.
            else
               allocate(bin_read1((nx+1)*(ny+1)))
               read(65)bin_read1
               idummy=1
               do i=1,nx+1
                  do j=1,ny+1
                     elevij(i,j)=bin_read1(idummy)
                     idummy=idummy+1
                  enddo
               enddo
               deallocate(bin_read1)
            endif
            close(65)
         endif          
         tstrt(:)=dur+delta

! mdw 4-15-2004 allocate positions for last deposition and initialize with 99999
         xnear(:)=99999
         ynear(:)=99999
         znear(:)=99999

!mdw 11/29/2005 do particle size stuff irrespective of ieradflag
         weightr(:)=0.
         weightr(1)=1.
         if(npartsize.eq.1) dpm_bin(1)=2.*dpart
         press=76.
         temp=293.
         rhop=1.

! MAN 08/23/2006 calculate total release
         totsrcmass=0.
         nreleases=0
         do isource=1,nsources
            select case(source_str_units(isource))
               case(1)
                  totsrcmass=totsrcmass+source_strength(isource)
               case(2)
                  totsrcmass=totsrcmass+source_strength(isource)*&
                                    source_dur(isource)
               case(3)
                  totsrcmass=totsrcmass+source_strength(isource)*&
                                    source_density(isource)
               case(4)
                  totsrcmass=totsrcmass+source_strength(isource)*&
                                    source_density(isource)*source_dur(isource)
            endselect
            if(source_geometry(isource).eq.2)then
               nsegs=source_stopidx(isource)-source_startidx(isource)
            else
               nsegs=1
            endif
            select case(source_release(isource))
               case(1)
                  nreleases=nreleases+nsegs
               case(2)
                  nreleases=nreleases+int(nsegs*(dur-source_start(isource))/delta)
               case(3)
                  nreleases=nreleases+int(nsegs*source_dur(isource)/delta)
            endselect
         enddo
! MAN 08/23/2006 calculate particle release rate for each source
         do isource=1,nsources
            select case(source_str_units(isource))
               case(1)
                  srcmass=source_strength(isource)
               case(2)
                  srcmass=source_strength(isource)*&
                                    source_dur(isource)
               case(3)
                  srcmass=source_strength(isource)*&
                                    source_density(isource)
               case(4)
                  srcmass=source_strength(isource)*&
                                    source_density(isource)*source_dur(isource)
            endselect
            if(source_geometry(isource).eq.2)then
               nsegs=source_stopidx(isource)-source_startidx(isource)
            else
               nsegs=1
            endif
            select case(source_release(isource))
               case(1)
                  irelease=nsegs
               case(2)
                  irelease=int(nsegs*(dur-source_start(isource))/delta)
               case(3)
                  irelease=int(nsegs*source_dur(isource)/delta)
            endselect
            source_prr(isource)=nint(real(nparticles)*srcmass/(real(nsegs*irelease)*totsrcmass))
            source_weight(isource)=real(nparticles)*srcmass/(real(source_prr(isource)*irelease)*totsrcmass)
         enddo
         pmass(:)=0.
         t=-delta !MAN 01/29/2007 t=0.
! setup concentration receptor grid
         dbx=(xbu-xbl)/real(nbx)
         dby=(ybu-ybl)/real(nby)
         dbz=(zbu-zbl)/real(nbz)
         bvol=dbx*dby*dbz
! mdw 1-22-2004 begin input for building infiltration model
         call turbinit
         if(inverseflag .eq. 1)then
            u=-u
            v=-v
            w=-w
         endif
         allocate(xbc(nbx),ybc(nby),zbc(nbz))
         
         do ib=1,nbx
            xbc(ib)=dbx*real(ib-1)+.5*dbx+xbl
         enddo
         do jb=1,nby
            ybc(jb)=dby*real(jb-1)+.5*dby+ybl
         enddo
         do kb=1,nbz
            zbc(kb)=dbz*real(kb-1)+.5*dbz+zbl
         enddo
         
         if(isiteflag.eq.0)then
            allocate(con(nbx,nby,nbz),dose(nbx,nby,nbz))
            con(:,:,:)=0.
            dose(:,:,:)=0.
            if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
               allocate(conSP(nbx,nby,nbz),doseSP(nbx,nby,nbz))
               conSP(:,:,:)=0.
               doseSP(:,:,:)=0.
            endif
         else
            allocate(dosexy(nbx,nby))
            dosexy(:,:)=0.
         endif
         if(idepositionflag.eq.1)then
            allocate(depcx(nx,ny,nz),depcy(nx,ny,nz),depcz(nx,ny,nz))
            depcx(:,:,:)=0.
            depcy(:,:,:)=0.
            depcz(:,:,:)=0.
            if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
               allocate(depcxSP(nx,ny,nz),depcySP(nx,ny,nz),depczSP(nx,ny,nz))
               depcxSP(:,:,:)=0.
               depcySP(:,:,:)=0.
               depczSP(:,:,:)=0.
            endif
         endif
         allocate(iface1n(inumbuild),iface2n(inumbuild),jface1n(inumbuild),jface2n(inumbuild),kface1n(inumbuild),kface2n(inumbuild))
         allocate(aminusx(inumbuild),aplusx(inumbuild),aminusy(inumbuild),aplusy(inumbuild),aroof(inumbuild),afloor(inumbuild))
         allocate(ibldgpts(nx*ny),jbldgpts(nx*ny))
         if(iindoorflag.eq.1)then
              call readQPindoor
         endif
lp005:   do ibuild=1,inumbuild
            if(bldtype(ibuild) .eq. 9)then
               cycle
            elseif(bldtype(ibuild).ne.3)then
               iface1=int((xfo(ibuild))/dx)+1
               iface1=max(iface1,2)
               iface2=int((xfo(ibuild)+Lt(ibuild)-dx)/dx)+1
               iface2=min(iface2,nx-2)
               ifacep=int((xfo(ibuild)+Lt(ibuild))/dx)+1
            else
               iface1=int((xfo(ibuild)-0.5*wti(ibuild))/dx)+1
               iface1=max(iface1,2)
               iface2=int((xfo(ibuild)+0.5*wti(ibuild)-dx)/dx)+1
               iface2=min(iface2,nx-2)
               ifacep=int((xfo(ibuild)+0.5*wti(ibuild))/dx)+1
               ifacep=min(ifacep,nx-2)
            endif
            nbldgpts=0
            ixfo=int(xfo(ibuild)/dx)+1
            jyfo=int(yfo(ibuild)/dy)+1
            aminusx(ibuild)=1.e-08
            aplusx(ibuild)=1.e-08
            aminusy(ibuild)=1.e-08
            aplusy(ibuild)=1.e-08
            aroof(ibuild)=1.e-08
            afloor(ibuild)=1.e-08
            jface1=int((yfo(ibuild)-0.5*wti(ibuild))/dy)+1
            jface1=max(jface1,2)
            jface2=int((yfo(ibuild)+0.5*wti(ibuild)-dy)/dy)+1
            jface2=min(jface2,ny-2)
! mdw 4-16-2004 added proper treatment of zfo
            kface1=int((zfo(ibuild)+dz)/dz)+1
            kface2=int((Ht(ibuild)+zfo(ibuild))/dz)+1
            kface2=min(nz-2,kface2)
            do iface=iface1,iface2
               do jface=jface1,jface2
                  xiurjursum=0.
                  yiurjursum=0.
                  do ixur=1,3
                     dxur=dx*real(ixur-1)/2.
                     do iyur=1,3
                        dyur=dy*real(iyur-1)/2.
                        xiurjur=xfo(ibuild)+(real(iface-ixfo)*dx+dxur)*cos(gamma(ibuild)*pi/180.)&
                              -(real(jface-jyfo)*dy+dyur)*sin(gamma(ibuild)*pi/180)
                        yiurjur=yfo(ibuild)+(real(jface-jyfo)*dy+dyur)*cos(gamma(ibuild)*pi/180)&
                              +(real(iface-ixfo)*dx+dxur)*sin(gamma(ibuild)*pi/180)
!                     enddo
!                  enddo
!                  xiurjur=.25*xiurjursum
!                  yiurjur=.25*yiurjursum
                        jfacr=int(yiurjur/dy)+1
                        ipntchk=0
                        ifacr=int(xiurjur/dx)+1
                        dxfo=dx*real(ifacr-1)+.5*dx-xfo(ibuild)
                        dyfo=dy*real(jfacr-1)+.5*dy-yfo(ibuild)
                        ximrot=dxfo*cos(gamma(ibuild)*pi/180)+dyfo*sin(gamma(ibuild)*pi/180)
                        yjmrot=dyfo*cos(gamma(ibuild)*pi/180)-dxfo*sin(gamma(ibuild)*pi/180)
                        imrot=ixfo+int(ximrot/dx)
                        jmrot=jyfo+int(yjmrot/dy)
                        if(imrot.lt.iface1.or.imrot.gt.iface2.or.jmrot.lt.jface1.or.jmrot.gt.jface2)cycle 
                        do nbldgpt=1,nbldgpts
                           if(ifacr.eq.ibldgpts(nbldgpt).and.jfacr.eq.jbldgpts(nbldgpt))then
                              ipntchk=1
                              ipntchk=ipntchk+1
                              exit
                           endif
                        enddo
                        if(ipntchk.ne.0)cycle
                        nbldgpts=nbldgpts+1
                        ibldgpts(nbldgpts)=ifacr
                        jbldgpts(nbldgpts)=jfacr
                        do kface=kface1,kface2
!                           if(kface.eq.kface2)write(77,*)ifacr,jfacr
                           if(icellflag(ifacr,jfacr,kface).eq.0)then
                              if(icellflag(ifacr-1,jfacr,kface).ne.0.) aminusx(ibuild)=1.+aminusx(ibuild)
                              if(icellflag(ifacr+1,jfacr,kface).ne.0.) aplusx(ibuild)=1.+aplusx(ibuild)
                              if(icellflag(ifacr,jfacr-1,kface).ne.0.) aminusy(ibuild)=1.+aminusy(ibuild)
                              if(icellflag(ifacr,jfacr+1,kface).ne.0.) aplusy(ibuild)=1.+aplusy(ibuild)
                              if(icellflag(ifacr,jfacr,kface-1).ne.0.) afloor(ibuild)=1.+afloor(ibuild)
                              if(icellflag(ifacr,jfacr,kface+1).ne.0.) aroof(ibuild)=1.+aroof(ibuild)
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            iface1n(ibuild)=iface1
            iface2n(ibuild)=iface2
            jface1n(ibuild)=jface1
            jface2n(ibuild)=jface2
            kface1n(ibuild)=kface1
            kface2n(ibuild)=kface2
         enddo   lp005
         !call new(windgrid,nx1,ny1,nz1,dx,dy,dz)
! end new building input for face averages and infiltration time step
         do iyg=1,nyg
            do ixg=1,nxg
               delta=delta_ref
               delta_ref_rm=delta_ref
               t=-delta !MAN 01/29/2007 t=0.
               call releaseparticles
               if(iindoorflag .eq. 1 .and. iyg .eq. 1 .and. ixg .eq. 1)then
                  call calculateBuildingInfiltExfilt
               endif
!               nsteps=nint(dur/delta) !+1
!MAN 01/29/2007 fixed bug in limits for particle domain bounds
!               xpma=real(nx-2)*dx-.01
!               ypma=real(ny-2)*dy-.01
!               zpma=real(nz-3)*dz-.01
               xpma=real(nx-1)*dx-z0
               ypma=real(ny-1)*dy-z0
               zpma=real(nz-2)*dz-z0
40             format(a40)
! MAN 9/13/2005 binary nextgrid
               if(inextgridflag .ge. 1)then
                  if(format_flag.eq.1.or.format_flag.eq.3)then
                     write(31,51111)nptot,totsrcmass
                  endif
                  if(format_flag.eq.2.or.format_flag.eq.3)then
                     write(32)totsrcmass
                  endif
               endif
               mp=nptot
               tls_min=delta
!               taylor_max=0.
! end MAN 9/13/2005
!mdw 4-22-2004 write out total number of particles and totsrcmass to nextgrid
51111          format(i9,e12.4)
               isourcestart=0
!               it=0
lp014:         do while(t .le. dur-0.1*delta) !MAN 01/29/2007 it=1,nsteps
!lp014:         do it=0,nsteps !MAN 01/29/2007 it=1,nsteps
                  if(isiteflag .eq. 0 .and. t .ge. 0.)then
                     print*,'Time =',t,'s' !real(it)*delta
                  endif
                  tout=tout+delta
! mdw added line for building infiltration
                  if(iindoorflag.eq.1)touti=touti+delta
                  if(t.ge.concstrt)toutc=toutc+delta
                  if(ibuoyflag .eq. 1 .and. inextgridflag .lt. 2)then
                     zbmaxp=zbmax
                     zbminp=zbmin
                     zbmax=0.
                     zbmin=h
                     wbulkzs=0.
                     sumwbulkz=0.
                  endif
                  if(ibuoyflag .eq. -1 .and. source_release(1) .eq. 1 .and. t .gt. 0.1*delta)then
                     do isource=1,nsources
                        hnewa=.5*(hnew(isource)+hnewp(isource))
                        rnewa=.5*(rnew(isource)+rnewp(isource))
                        reffa=.5*(reff(isource)+reffp(isource))
                        dhdt(isource)=(hnew(isource)-hnewp(isource))/delta
                        drdta=(rnew(isource)-rnewp(isource))/delta
                        xscena=.5*(xscen(isource)+xscenp(isource))
                        yscena=.5*(yscen(isource)+yscenp(isource))
                        iscena=int(xscena/dx)+1
                        iscena=max(1,iscena)
                        iscena=min(nx-1,iscena)
                        jscena=int(yscena/dy)+1 
                        jscena=max(1,jscena)
                        jscena=min(ny-1,jscena)                      
                        write(67,60006)t,hnewp(isource),rnewp(isource),reffp(isource),&
                                       xscenp(isource),yscenp(isource),rhocent(isource)
60006                   format(6f10.2,f7.4)
                     enddo
                  endif
                  ninactiveparticles=0
                  
                  if((ibioslurryflag .eq. 1 .or. itwophaseflag .eq. 1) .and. t .gt. 0.)call evaporate
                  
lp006:            do m=1,mp
! MAN 01/24/2007 Particle splitting
                     if(tstrt(m) .gt. t+0.001*delta)then !MAN 01/29/2007 tstrt(m) .ge. t
                        if(tstrt(m) .gt. t+1.001*delta)then
                           ninactiveparticles=ninactiveparticles+1
                        endif
                        cycle
                     endif
!end MAN 01/24/2007
                     dt=delta
                     tre=delta
                     nstps=0
                     dhdta=0.
                     xran=0.
                     yran=0.
                     zran=0.
                     ubulk=0.
                     vbulk=0.
                     densegas_factor=0.
                     call int3cal(dpm(m),sc,taud,rhop,xnu,press,temp,diff,vs)
                     forcez=0.
                     if(ibuoyflag .eq. 1 .and. zp(m)>0. .and. inextgridflag .lt. 2 .and.&
                           activepartt .gt. 10*number_cutoff)then
                        km=int((zp(m)+dz)/dz)+1
                        airtempzp=airtemp(km)+(zp(m)-zi(km))*(airtemp(km+1)-airtemp(km))/dz
                        kbpuf=int((zp(m)-zbminp)/dzbpuf)+1
                        kbpuf=min(10,kbpuf)
                        kbpuf=max(1,kbpuf)
                        offcntdis=sqrt(((xp(m)-xic(kbpuf))/(2.*sigxb(kbpuf)))**2+((yp(m)-yic(kbpuf))/(2.*sigyb(kbpuf)))**2)
                        if(t .gt. 0.1*delta)then
                           entrain_factor=1.-min(contot/contotp,1.)
                        else
                           entrain_factor=0.
                        endif
                        temp_part(m)=temp_part(m)-.01*(zp(m)-zpp(m))    !adjust for adiabatic cooling
                        temp_part(m)=temp_part(m)-entrain_factor*(temp_part(m)-airtempzp)
                        wbulkz(m)=wbulkz(m)*(1.-entrain_factor)
                        forcez=9.8*(temp_part(m)-airtempzp)/airtempzp
                        zpp(m)=zp(m)
                     endif                       
! for use with smaller timesteps near walls or floors
                     do while(tre.gt.1.e-04*delta)
                        im=int(xp(m)/dx)+1
                        jm=int(yp(m)/dy)+1
                        km=int((zp(m)+dz)/dz)+1
                        dxc=xp(m)-.5*dx-dx*real(im-1)
                        dyc=yp(m)-.5*dy-dy*real(jm-1)
                        dzc=zp(m)-.5*dz-dz*real(km-2)
! dwxp is the distance to the wall in the plus x direction
! dwxm is the distance to the wall in the minus x direction
! dwyp is the distance to the wall in the plus y direction
! dwym is the distance to the wall in the minus y direction
! drzp is the distance to the roof in the plus z direction
! dfzm is the distance to the floor in the minus z direction
                        if((im .ge. 1 ).and.( im .le. nx-1).and.(jm .ge. 1).and. (jm .le. ny-1).and.&
                              (km .gt. 1) .and. (km .le. nz-1) .and. (weight(m) .gt. 1.e-05))then
!mdw 11-30-2006 corrected the statements to update dxc, dyc, dzc, and wall distances anytime
! the particles are moved because they are within znaut of the wall
                           if(dwxp<z0*z0coeff .and. iqperror .lt. iqperrormax)then
                              write(77,*) "Reflection Error stop 1"
                              write(77,14001)t,m,tre,dwxp
                              write(77,14002)im,jm,km,xp(m),yp(m),zp(m)
                              xp(m)=xp(m)-z0*z0coeff
                              dxc=xp(m)-.5*dx-dx*real(im-1)
                              dwxp=dxp(im,jm,km)-dxc
                              iqperror=iqperror+1
                           endif
                           dwxm=dxm(im,jm,km)+dxc
                           if(dwxm<z0*z0coeff .and. iqperror .lt. iqperrormax)then
                              write(77,*) "Reflection Error stop 2"
                              write(77,14001)t,m,tre,dwxm
                              write(77,14002)im,jm,km,xp(m),yp(m),zp(m)
                              xp(m)=xp(m)+z0*z0coeff
                              dxc=xp(m)-.5*dx-dx*real(im-1)
                              dwxm=dxm(im,jm,km)+dxc
                              iqperror=iqperror+1
                           endif
14001 format('t=',f8.2,' m=',i6,' tre=',f5.2,' dw=',f15.8)
14002 format('im=',i6,' jm=',i6,' km=',i6,' xp=',f8.2,' yp=',2f6.2)
                           dwyp=dyp(im,jm,km)-dyc
                           if(dwyp<z0*z0coeff .and. iqperror .lt. iqperrormax)then
                              write(77,*) "Reflection Error stop 3"
                              write(77,14001)t,m,tre,dwyp
                              write(77,14002)im,jm,km,xp(m),yp(m),zp(m)
                              yp(m)=yp(m)-z0*z0coeff
                              dyc=yp(m)-.5*dy-dy*real(jm-1)
                              dwyp=dyp(im,jm,km)-dyc
                              iqperror=iqperror+1
                           endif
                           dwym=dym(im,jm,km)+dyc
                           if(dwym<z0*z0coeff .and. iqperror .lt. iqperrormax)then
                              write(77,*) "Reflection Error stop 4"
                              write(77,14001)t,m,tre,dwym
                              write(77,14002)im,jm,km,xp(m),yp(m),zp(m)
                              yp(m)=yp(m)+z0*z0coeff
                              dyc=yp(m)-.5*dy-dy*real(jm-1)
                              dwym=dym(im,jm,km)+dyc
                              iqperror=iqperror+1
                           endif
!                           if(dwym.lt.z0)yp(m)=yp(m)+1.01*z0
                           dfzm=dzm(im,jm,km)+dzc
                           if(dfzm<z0*z0coeff .and. iqperror .lt. iqperrormax)then
                              write(77,*) "Reflection Error stop 5"
                              write(77,14001)t,m,tre,dfzm
                              write(77,14002)im,jm,km,xp(m),yp(m),zp(m)
                              zp(m)=zp(m)+z0*z0coeff
                              dzc=zp(m)-.5*dz-dz*real(km-2)
                              dfzm=dzm(im,jm,km)+dzc
                              iqperror=iqperror+1
                           endif
                           drzp=dzp(im,jm,km)-dzc
                           if(drzp<z0*z0coeff .and. iqperror .lt. iqperrormax)then
                              write(77,*) "Reflection Error stop 6"
                              write(77,14001)t,m,tre,drzp
                              write(77,14002)im,jm,km,xp(m),yp(m),zp(m)
                              zp(m)=zp(m)-z0*z0coeff
                              dzc=zp(m)-.5*dz-dz*real(km-2)
                              drzp=dzp(im,jm,km)-dzc
                              iqperror=iqperror+1
                           endif
                           if(icellflag(im,jm,km).eq.0)then
                              write(*,14007)t,m
                              xp(m)=0.
                              yp(m)=0.
                              zp(m)=-.5*dz
                              cycle
                           endif
14007                      format('t=',f10.2,' and m =',i6,'; particle in wall removing from domain')
                           ustar=max(ustarij(im,jm,km),3.E-02)
                           sigu=sigui(im,jm,km)
                           sigv=sigvi(im,jm,km)
                           sigw=sigwi(im,jm,km)
                           dsigvdn=dsigvdni(im,jm,km)
                           dsigudn=dsigudni(im,jm,km)
                           dsigwdn=dsigwdni(im,jm,km)
                           tau13=-upwpi(im,jm,km)
                           if(ibuoyflag .eq. 1 .and. km .lt. nz-1 .and. inextgridflag .lt. 2 .and. &
                                 activepartt .gt. 10*number_cutoff)then
                              fcoef=(1./(2.*pi*dzbpuf*sigxb(kbpuf)*sigyb(kbpuf)))*exp(-.5*((xp(m)-xic(kbpuf))/&
                                    sigxb(kbpuf))**2-.5*((yp(m)-yic(kbpuf))/sigyb(kbpuf))**2)
                              if(fcoef.gt.0)then   !if on fcoef>lim1
                                 if(wbulkz(m)+forcez*dt*.5.gt.0)then   !if on wbulkz+forcez
                                    vref=(wbulkz(m)+forcez*dt)
                                    tscale=delta
                                    if(nstps.eq.0)then
                                       if(dt.eq.delta)then
                                          sigwb=.15*(wbulkz(m)+forcez*dt)
                                          vdef=.8*sigwb
                                          tscale=delta
                                          sigub=sigwb
                                          if(sigu.gt.sigub)sigub=0.
                                          sigvb=sigub
                                          call rangen(x)
                                          uranb=x*sigub
                                          call rangen(x)
                                          vranb=x*sigvb
                                          call rangen(x)
                                          wranb=-abs(x)*sigwb
                                       else
                                          tscale=delta
                                          vdef=.8*.15*(wbulkz(m)+forcez*dt)
                                          fac=1.-exp(-delta/tscale)
                                          sigwbn=.15*(wbulkz(m)+forcez*dt)*fac
                                          delrpuf=sqrt(sigxb(kbpuf)**2+sigyb(kbpuf)**2)*(sqrt(vref/(vref-vdef))-1.)
                                          sigubn=sigwbn
                                          if(sigu.gt.sigubn)sigubn=0.
                                          sigvbn=sigubn
                                          call rangen(x)
                                          uranb=x*sigubn+uranbp(m)*exp(-dt/tscale)
                                          call rangen(x)
                                          vranb=x*sigvbn+vranbp(m)*exp(-dt/tscale)
                                          call rangen(x)
                                          wranb=-abs(x)*sigwbn+wranbp(m)*exp(-dt/tscale)
                                          if(abs(wranb).gt.vref)wranb=-vref
                                       endif
                                    endif
                                    uranbp(m)=uranb
                                    vranbp(m)=vranb
                                    wranbp(m)=wranb
                                    wstar3=9.8*(zirhoht-zirholt)*rhocell(im,jm,km)*(wbulkz(m)+forcez*dt)/rhobkg(km)
                                    wstar=wstar3**.3333
                                    wstar=0.
                                    if(wstar.gt.0.)then
                                       rclltg=-kkar*(wstar/ustar)**3./(zirhoht-zirholt)
                                       rclloc=rclltg
                                    else
                                       rclltg=0.
                                    endif
                                 else
                                    rclloc=rcl
                                 endif
                              else
                                 rclloc=rcl
                                 wbulkz(m)=0.
                              endif
                              dtau13dn=-dupwpdni(im,jm,km)
                              ubuoyl=0.
                              vbuoyl=0.
                           endif
                           if(ibuoyflag.lt.0)then    !if on negative ibuoy @3416
                              isource=nint(source_num(m))
                              idxp=idxpm(m)
!                              if(rhocell(im,jm,2).gt.1.001*rhoair.and.zp(m).le.hnewij(im,jm)+z0 .and.source_release(1).ge.2)then  !if @ 3417
                               if(source_release(1).ge.2)then
                                  deltaxpart=xp(m)-xps(m)
                                  deltaypart=yp(m)-yps(m)
                                  deltatot=sqrt(deltaxpart**2+deltaypart**2)
                                  idx1=idxpm(m)
                                  idx2=idxpm(m)+10
                                  idx2=min(idx2,idxmax)
                                  do idx=idx1,idx2
                                     distance_ref=sqrt((xsceni(idx)-xsceni(1))**2+(ysceni(idx)-ysceni(1))**2)
                                     deltatotref=real(idx-1)*delx
                                     deltatotrefp=deltatotref-delx
                                     if(idx.ne.1)then                            
                                                deltatot=sqrt((deltaxpart*(xsceni(idx)-xsceni(1))/distance_ref)**2&
                                        +(deltaypart*(ysceni(idx)-ysceni(1))/distance_ref)**2)
                                              else
                                                 deltatot=sqrt((deltaxpart*cosphidgi(1))**2+(deltaypart*sinphidgi(1))**2)
                                             endif

                                     if(abs(deltatotref-deltatot).gt.abs(deltatotrefp-deltatot).and.idx.ne.1)then
                                        exit
                                     else
                                        deltatotrefp=deltatotref
                                     endif
                                   enddo
                                   if(idx.ne.idx2)then
                                      idx=idx-1
                                   endif
                                   idxpm(m)=idx
                                   hnew(isource)=hnewi(idx)
                                   hnewa=hnew(isource)
                                   rnew(isource)=rnewi(idx)
                                   crat=(hnewi(1)*rnewi(1))/(hnewi(idx)*rnewi(idx))
                                   rhocent(isource)=source_density(isource)*crat+rhoair*(1.-crat)
                                   if(rhocent(isource).gt.1.001*rhoair .and. zp(m).le.(hnewi(idxp)+z0))then                   
                                    densegas_factor=1.
                                    if(idx.ne.idxp)then
                                       zp(m)=(hnewi(idx)/hnewi(idxp))*(zp(m)-hgt(im,jm))+hgt(im,jm)
                                    endif
                                 ristar=9.8*(rhocent(isource)-rhoair)*(hnewa+z0)/(rhoair*ustar**2)
                                 dhdt(isource)=dhdtij(im,jm)
!                                 idx=idxij(im,jm)
                                 ubar(isource)=ubari(idx)
                                 drdx=drdxi(idx)
                                 rnew(isource)=rnewi(idx)
                                 cosphidg(isource)=cosphidgi(idx)
                                 sinphidg(isource)=sinphidgi(idx)
                                 xscen(isource)=xsceni(idx)
                                 yscen(isource)=ysceni(idx)
                                 rempty=max(rminusij(idx),rplusij(idx))
                                 ristar=max(1.e-12,ristar)
                                 rhowt1=rhocell(im,jm,2)/source_density(isource)
                                 rhowt0=1.-rhowt1
                                 ustarp=ustar/(1+.1125*ristar**1.04)
                                 ustarb=ustarp
!                                 crat=(rhocell(im,jm,km)-rhoair)/(source_density(isource)-rhoair)
                                 if(zp(m).le.hnewi(idx)+z0)then
                                    ustarb=ustar/(1.+.2*ristar)
                                 endif
                                 yloc=(yp(m)-yscen(isource))*cosphidg(isource)-(xp(m)-xscen(isource))*sinphidg(isource)
                                 if((yloc.ge.0. .and. yloc.le.rnewpi(idx)+rplusij(idx)+drdxpi(idx)*delx).or.&
                                    (yloc.lt.0. .and. -yloc.le.rnewmi(idx)+rminusij(idx)+drdxmi(idx)*delx))then  !if on yloc @ 3451
                                    sigv=sigv*ustarb/ustar
                                    sigu=sigu*ustarb/ustar
                                    sigw=sigw*ustarb/ustar
                                    tau13=tau13*(ustarb/ustar)**2
                                    siguo=siguo*(ustarb/ustar)
                                    sigvo=sigvo*(ustarb/ustar)
                                    sigwo=sigwo*(ustarb/ustar)
                                    dsigwdn=dsigwdn*(ustarb/ustar)
                                    uran(m)=uran(m)*ustarb/ustar
                                    vran(m)=vran(m)*ustarb/ustar
                                    wran(m)=wran(m)*ustarb/ustar
                                    ustar=ustarb
                                    if(dt.eq.delta)zranw=(1./delta)*(ran2()-.5)*.5*min(dz,hnew(isource))
                                    zran=dt*zranw
                                    if(zran+zp(m).gt..999*hnew(isource))then
                                       zranx=zran+zp(m)-.999*hnew(isource)
                                       zran=zran-2.*zranx
                                       zranw=-zranw
                                    endif
                                    if(zran+zp(m) .lt. z0*z0coeff)then
                                       zranx=z0*z0coeff-zran-zp(m)
                                       zran=zran+2.*zranx
                                       zranw=-zranw
                                    endif
                                    if(dt.eq.delta)xranw=(1./delta)*(ran2()-.5)*.5*min(dx,dy,delx)
                                    if(dt.eq.delta)yranw=(1./delta)*(ran2()-.5)*.5*min(dx,dy,rnew(isource))
                                    xran=dt*(-yranw*sinphidg(isource)+xranw*cosphidg(isource))
                                    yran=dt*(yranw*cosphidg(isource)+xranw*sinphidg(isource))      
                                    if(dwxm+xran .le. z0*z0coeff)then
                                       xranx=z0*z0coeff-xran+dwxm
                                       xran=xran+2.*xranx
                                       xranw=-xranw
                                    endif
                                    if(dwxp-xran .lt. z0*z0coeff)then
                                       xranx=xran-dwxp+z0*z0coeff
                                       xran=xran-2.*xranx
                                       xranw=-xranw
                                    endif
                                    if(dwym+yran .le. z0*z0coeff)then
                                       xrany=z0*z0coeff-yran+dwym
                                       yran=yran+2.*xrany
                                       yranw=-yranw
                                    endif
                                    if(dwyp-yran .lt. z0*z0coeff)then
                                       xrany=yran-dwyp+z0*z0coeff
                                       yran=yran-2.*xrany
                                       yranw=-yranw
                                    endif
                                    ylocran=(yp(m)+yran-yscen(isource))*cosphidg(isource)&
                                            -(xp(m)+xran-xscen(isource))*sinphidg(isource)
                                    if((ylocran.gt.rplusij(idx)+rnewpi(idx)).and. (yran*cosphidg(isource).gt.0 &
                                          .or. -xran*sinphidg(isource).gt.0))then
                                       xcessr=ylocran-rplusij(idx)-rnewpi(idx)
                                       drloc=yloc-rnewpi(idx)-rplusij(idx)
                                       if(yran*cosphidg(isource).ge.-xran*sinphidg(isource))then
                                          yran=yran-2.*(xcessr-drloc)
                                          if(yranw*cosphidg(isource).gt.0)yranw=-yranw
                                       else
                                          xran=xran-2.*(xcessr-drloc)
                                          if(-xranw*sinphidg(isource).gt.0)xranw=-xranw
                                       endif
                                    endif
                                    if((-ylocran.gt.rminusij(idx)+rnewmi(idx)).and.(-yran*cosphidg(isource).gt.0. &
                                          .or. xran*sinphidg(isource).gt.0))then
                                       xcessr=-ylocran-rminusij(idx)-rnewmi(idx)
                                       drloc=-yloc-rnewmi(idx)-rminusij(idx)
                                       if(-yran*cosphidg(isource).ge.+xran*sinphidg(isource))then
                                          yran=yran+2.*(xcessr-drloc)
                                          if(-yranw*cosphidg(isource).gt.0)yranw=-yranw
                                       else
                                          xran=xran+2.*(xcessr-drloc)
                                          if(xranw*sinphidg(isource).gt.0)xranw=-xranw
                                       endif
                                    endif
                                    if(yloc.ge.0.)then
                                       vbulk=ubar(isource)*drdxpi(idx)*cosphidg(isource)*yloc&
                                             /(rplusij(idx)+rnewpi(idx))+vbulkpy(idx)
                                       ubulk=-drdx*ubar(isource)*sinphidg(isource)*yloc/&
                                             (rplusij(idx)+rnewpi(idx))+ubulkpx(idx)
                                    else
                                       vbulk=-drdx*ubar(isource)*cosphidg(isource)*(-yloc)/(rminusij(idx)+rnewmi(idx))+vbulkpy(idx)
                                       ubulk=drdx*ubar(isource)*sinphidg(isource)*(-yloc)/(rminusij(idx)+rnewmi(idx))+ubulkpx(idx)
                                    endif
                                 else
                                    densegas_factor=0.
                                 endif
                              else
                                 rhowt1=0.
                                 rhowt0=1.
                                 densegas_factor=0.
                              endif                                
                           else
                              rhowt1=0.
                              rhowt0=1.
                              densegas_factor=0.
                           endif
                              if(rhocent(isource).gt.1.001*rhoair.and.source_release(1).eq.1)then
!                                 zp(m)=zp(m)*(hnewp(isource)+z0+(hnew(isource)-hnewp(isource))*(1.-tre/delta))/(hnewp(isource)+z0)
                                 hnewa=hnewp(isource)+(hnew(isource)-hnewp(isource))*(1.-tre/delta)
                                 rnewa=rnewp(isource)+(rnew(isource)-rnewp(isource))*(1.-tre/delta)
                                 reffa=reffp(isource)+(reff(isource)-reffp(isource))*(1.-tre/delta)
                                 xscena=xscenp(isource)+(xscen(isource)-xscenp(isource))*(1.-tre/delta)
                                 yscena=yscenp(isource)+(yscen(isource)-yscenp(isource))*(1.-tre/delta)
                                 if(zp(m).lt.(hnewa+z0))then
                                    ristar=9.8*(rhocent(isource)-rhoair)*(hnewa+z0)/(rhoair*ustar**2)
!                                    dhdta=dhdt(isource)*zp(m)/hnewa
                                    drdta=.5*(drdt(isource)+drdtp(isource))
                                    rptcen=sqrt((xp(m)-xscena)**2+(yp(m)-yscena)**2)
                                    if(rptcen.le.reffa)then
                                       dhdta=dhdt(isource)*zp(m)/hnewa
                                       if(rptcen.gt.1.e-06*dx)then
                                          cosrptcen=(xp(m)-xscena)/rptcen
                                          sinrptcen=(yp(m)-yscena)/rptcen
                                       else
                                          cosrptcen=1.
                                          sinrptcen=0.
                                       endif
                                       ubulk=drdta*cosrptcen*rptcen/reffa+vslopex(isource)
                                       vbulk=drdta*sinrptcen*rptcen/reffa+vslopey(isource)
                                       crat=(rhocent(isource)-rhoair)/(source_density(isource)-rhoair)
                                       ustarb=ustar/(1.+.2*ristar)
                                       sigv=sigv*ustarb/ustar
                                       sigu=sigu*ustarb/ustar
                                       sigw=sigw*ustarb/ustar
                                       tau13=tau13*(ustarb/ustar)**2
                                       siguo=siguo*(ustarb/ustar)
                                       sigvo=sigvo*(ustarb/ustar)
                                       sigwo=sigwo*(ustarb/ustar)
                                       dsigwdn=dsigwdn*(ustarb/ustar)
                                       uran(m)=uran(m)*ustarb/ustar
                                       vran(m)=vran(m)*ustarb/ustar
                                       wran(m)=wran(m)*ustarb/ustar
                                       ustar=ustarb
                                       if(dt.eq.delta)zranw=(1./delta)*(ran2()-.5)*.5*min(dz,.5*hnewa)
                                       zran=dt*zranw
                                       if(zran+zp(m).gt..999*(hnewa+z0))then
                                          zranx=zran+zp(m)-.999*(hnewa+z0)
                                          zran=zran-2.*zranx
                                          zranw=-zranw
                                       endif
                                       if(zran+zp(m) .lt. z0*z0coeff)then
                                          zranx=z0*z0coeff-zran-zp(m)
                                          zran=zran+2.*zranx
                                          zranw=-zranw
                                       endif
                                       if(dt.eq.delta)xranw=(1./delta)*(ran2()-.5)*.5*min(dx,.5*rnew(isource))
                                       if(dt.eq.delta)yranw=(1./delta)*(ran2()-.5)*.5*min(dy,.5*rnew(isource))
                                       xran=dt*xranw
                                       yran=dt*yranw
                                       if(dwxm+xran .le. z0*z0coeff)then
                                          xranx=z0*z0coeff-xran+dwxm
                                          xran=xran+2.*xranx
                                          xranw=-xranw
                                       endif
                                       if(dwxp-xran .lt. z0*z0coeff)then
                                          xranx=xran-dwxp+z0*z0coeff
                                          xran=xran-2.*xranx
                                          xranw=-xranw
                                       endif
                                       if(dwym+yran .le. z0*z0coeff)then
                                          xrany=z0*z0coeff-yran+dwym
                                          yran=yran+2.*xrany
                                          yranw=-yranw
                                       endif
                                       if(dwyp-yran .lt. z0*z0coeff)then
                                          xrany=yran-dwyp+z0*z0coeff
                                          yran=yran-2.*xrany
                                          yranw=-yranw
                                       endif
! Now we check to make sure that we won't move the particle into another source plume
                                       iran=int((xp(m)+xran)/dx)+1
                                       jran=int((yp(m)+yran)/dy)+1
                                       xranx=0.
                                       xrany=0.
                                       if(icellflag_dense(iran,jran).ne.isource .and. icellflag_dense(iran,jran) .ne. 0)then
                                          if(iran.ne.im .and. icellflag_dense(iran,jm) .ne. isource&
                                             .and. icellflag_dense(iran,jm) .ne. 0 )then
                                             isign=iran-im
                                             xranx=(xp(m)+xran-dx*real(iran-isign))*real(isign)
                                             xran=xran-2.*real(isign)*xranx
                                             xranw=-xranw
                                          endif                                      
                                          if(jran.ne.jm .and. icellflag_dense(im,jran) .ne. isource&
                                             .and. icellflag_dense(im,jran) .ne. 0 )then
                                             jsign=jran-jm
                                             xrany=(yp(m)+yran-dy*real(jran-jsign))*real(jsign)
                                             yran=yran-2.*real(jsign)*xrany
                                             yranw=-yranw
                                          endif
                                          if(xranx.eq.0.and.xrany.eq.0)then
                                             if(im.eq.iran.and.jm.eq.jran)then
                                                isource2=icellflag_dense(iran,jran)
                                                xsign=sign(1.,xscen(isource2)-xp(m)-xran)
                                                ysign=sign(1.,yscen(isource2)-yp(m)-yran)
                                                xran=-xsign*abs(xran)
                                                xranw=-xsign*abs(xranw)
                                                yran=-ysign*abs(yran)
                                                yranw=-ysign*abs(yranw)
                                             else
                                                isign=iran-im
                                                xranx=(xp(m)+xran-dx*real(iran-isign))*real(isign)
                                                xran=xran-2.*real(isign)*xranx
                                                jsign=jran-jm
                                                xrany=(yp(m)+yran-dy*real(jran-jsign))*real(jsign)
                                                yran=yran-2.*real(jsign)*xrany
                                                xranw=-xranw
                                                yranw=-yranw
                                             endif
                                          endif
                                       endif
                                       rlocran=sqrt((yp(m)+yran-yscena)**2+(xp(m)+xran-xscena)**2)
                                       if(rlocran.gt.reffa)then
                                          rcessr=rlocran-reffa
                                          if(abs(xp(m)+xran-xscena).ge.abs(yp(m)+yran-yscena))then
                                             if(xp(m)+xran.gt.xscena)then
                                                xran=xran-2.*rcessr
                                             else
                                                xran=xran+2.*rcessr
                                             endif
                                             xranw=-xranw
                                          else
                                             if(yp(m)+yran.gt.yscena)then
                                                yran=yran-2.*rcessr
                                             else
                                                yran=yran+2.*rcessr
                                             endif
                                             yranw=-yranw
                                          endif
                                       endif
                                    endif
                                 endif
                              endif
                           endif
                           tau11=sigu**2
                           tau22=sigv**2
                           tau33=sigw**2
                           dwall=min(dwxp,dwxm,dwyp,dwym,dfzm,drzp)
                           dwall=max(dwall,z0*z0coeff)
                           dwallp(m)=dwall
                           ivrelch=0
                           call windterp
                           wint=wint-.01*vs
! settling velocity is in cm/sec
                           utot=sqrt(uint**2+vint**2+wint**2)
                           utot=utot+.000001
                           dutotdn=dutotdni(im,jm,km)
                           dutotds=dutotdsi(im,jm,km)
                           cosu=(uint+.000001)/utot
! we have added .000001 to u and utot to prevent difficulties with stagnation points
                           sinv=vint/utot
! 57                       u1=uran(m)*cosu+vran(m)*sinv
                           tls = 0.0
                           iranvelerror=0
lp007:                     do while(dt.ge.0.0)
! MAN 02/09/2007 removed superfluous do while loop
                              u1=uran(m)*alphn1ij(im,jm,km)+vran(m)*alphn2ij(im,jm,km)+wran(m)*alphn3ij(im,jm,km)
                              v1=uran(m)*betn1ij(im,jm,km)+vran(m)*betn2ij(im,jm,km)+wran(m)*betn3ij(im,jm,km)
                              w1=uran(m)*gamn1ij(im,jm,km)+vran(m)*gamn2ij(im,jm,km)+wran(m)*gamn3ij(im,jm,km)
                              if(rcl.lt.0)then
                                 coeps=5.7*(ustar**3*(1.-.75*zp(m)*rcl)*(1.-.85*zp(m)/h)**(1.5)/(dwall*kkar))
                              else
                                 coeps=5.7*(ustar**3*(1.+3.7*zp(m)*rcl)*(1.-.85*zp(m)/h)**(1.5)/(dwall*kkar))
                              endif
                              if(ibuoyflag .eq. 1 .and. wstar .gt. ustar &
                                    .and. inextgridflag .lt. 2)coeps=5.7*wstar**3/(zirhoht-zirholt)
! MAN 02/08/2007 Taylor microscale lower limit to subtimesteps.
                              if(taylor_flag .gt. 0)then
                                 taylor_microscale=sqrt((8.55e-3)*xnu*sigu*sigu/(coeps*utot*utot))
! the constant in the Taylor microscale calculation is the result of the following expression 15*5.7*1e-4
! which takes into account the conversion of units
                              else
                                 taylor_microscale=0.
                              endif
                              tls=2.*sigw**2/(coeps+1.e-36)
                              if(coeps.eq.0. .and. iqperror .lt. iqperrormax)then
                                 write(77,*)'Dissipation Rate Error',coeps,m,t,tre
                                 iqperror=iqperror+1
                              endif
                              if(tls .lt. tls_min)then
                                  tls_min=tls
!                                 taylor_max=taylor_microscale
                              endif
                              if(dt .gt. 0.5*tls .and. dt .gt. taylor_microscale)dt=0.5*tls
!                              write(88)dt
!                              do while(dt.ge.0.0)
!                                 u1=uran(m)*alphn1ij(im,jm,km)+vran(m)*alphn2ij(im,jm,km)+wran(m)*alphn3ij(im,jm,km)
!                                 v1=uran(m)*betn1ij(im,jm,km)+vran(m)*betn2ij(im,jm,km)+wran(m)*betn3ij(im,jm,km)
!                                 w1=uran(m)*gamn1ij(im,jm,km)+vran(m)*gamn2ij(im,jm,km)+wran(m)*gamn3ij(im,jm,km)
!                                 utestm=u1*alph1ij(im,jm,km)+v1*alph2ij(im,jm,km)+w1*alph3ij(im,jm,km)
!                                 vtestm=u1*bet1ij(im,jm,km)+v1*bet2ij(im,jm,km)+w1*bet3ij(im,jm,km)
!                                 wtestm=u1*gam1ij(im,jm,km)+v1*gam2ij(im,jm,km)+w1*gam3ij(im,jm,km)
!                                 dyct=dyc*cosu-dxc*sinv
!                                 ustarloc=ustarz(im,jm,km)
!                                 dwall=min(dwxp,dwxm,dwyp,dwym,dfzm,drzp)
!                                 coeps=5.7*(ustarloc**3/(dwall*kkar)+ustarg(im,jm,km)**3/(elzg(im,jm,km)+.01))
!                                 if(rcl.lt.0)then
!                                    coeps=5.7*(ustar**3*(1.-.75*zp(m)*rcl)*(1.-.85*zp(m)/h)**(1.5)/(dwall*kkar))+0.3*zp(m)*rcl
!                                 else
!                                    coeps=5.7*(ustar**3*(1.+3.7*zp(m)*rcl)*(1.-.85*zp(m)/h)**(1.5)/(dwall*kkar))
!                                 endif
!                                 if(ibuoyflag.eq.1.and.wstar.gt.ustar.and.inextgridflag.lt.2)coeps=5.7*wstar**3/(zirhoht-zirholt)
!                                 tls=2.*sigw**2/(coeps+1.e-08)

!                                 if(taylor_flag .gt. 0)then
!                                    taylor_microscale=sqrt((2.85e-3)*(xnu/coeps)*&
!                                                      ((sigu**2.)+(sigv**2.)+(sigw**2.))/(utot**2.))
!!                                                   print*,taylor_microscale,coeps/5.7,utot,sigu,sigv,sigw,zp(m),xnu,dwall
!                                 else
!                                    taylor_microscale=0.
!                                 endif
!                                 if(tls .lt. tls_min)then
!                                     tls_min=tls
!                                     taylor_max=taylor_microscale
!                                 endif
!                                 if(tls .ge. 2.*dt .or. dt .lt. taylor_microscale) exit
!                                 dt=.50*tls
!                              enddo
! end MAN 02/09/2007
                              tau13=ustar**2*(1.-zi(km)/h)**(1.5)
                              lam11=1./(tau11-tau13**2/tau33)
                              lam22=1./tau22
                              lam33=1./(tau33-tau13**2/tau11)
                              lam13=-tau13/(tau11*tau33-tau13**2)
!mdw 7-25-2005 changed form of lam13 for cases where tau13=0.
                              call rangen(x)
                              dudet=-.5*coeps*(lam11*u1+lam13*w1)*dt+dutotdn*w1*dt &
                                    +utot*dutotds*dt+u1*dutotds*dt+sigu/(1.3*2.)*dsigwdn*dt &
                                    +sigu*dsigwdn*(u1*(lam11*5./1.3+lam13/(1.25*1.3))+w1*(lam13*4./1.3 &
                                    +lam33/1.3))*(0.5*w1)*dt-w1*dutotdni(im,jm,km)*dt
! 7-08-2005 MDW corrected sign on last term involving wprime dudz
                              duran=sqrt(coeps)*x*sqrt(dt)
                              call rangen(x)
                              dvdet=-.5*coeps*lam22*v1*dt+(2./1.3)*sigv*dsigwdn*lam22*v1*w1*dt
                              dvran=sqrt(coeps)*x*sqrt(dt)
                              call rangen(x)
                              dwdet=-.5*coeps*(lam13*u1+lam33*w1)*dt &
                                     +sigw*dsigwdn*dt+((lam13+lam11/1.69)*u1+(lam13/1.69+lam33)*w1)*w1*sigw*dsigwdn*dt
                              dwran=sqrt(coeps)*x*sqrt(dt)
                              nstps=nstps+1
                              if(nstps.gt.nstpmax)then
                                  nstpmax=nstps
                                  nstpm=m
                              endif
                              vrel2=sqrt((u1+dudet)**2+(v1+dvdet)**2+(w1+dwdet)**2)
                              vrel1=sqrt((u1)**2+(v1)**2+(w1)**2)
                              vrel=max(vrel2,vrel1)
                              if(iranvelerror .ge. 5)then
                                 write(77,*)'Random Velocity Error infinit loop U =',u(im,jm,km),'V =',v(im,jm,km),'W =',w(im,jm,km)
                                 print*,'Random Velocity Error infinit loop U =',u(im,jm,km),'V =',v(im,jm,km),'W =',w(im,jm,km)
                              endif
                              if(vrel.gt.4.*sigu.and.ivrelch.eq.0.or.vrel.gt.100. .and. iranvelerror .lt. 5)then
                                 iranvelerror=iranvelerror+1
                                 if(vrel.gt.100. .and. iqperror .lt. iqperrormax)then
                                    write(77,*)t,m,vrel,xp(m),yp(m),zp(m),dsigwdn
                                    iqperror=iqperror+1
                                 endif
                                 call rangen(x)
                                 u1=sigu*x
                                 call rangen(x)
                                 v1=sigv*x
                                 call rangen(x)
                                 w1=sigw*x
                                 uran(m)=u1*alph1ij(im,jm,km)+v1*alph2ij(im,jm,km)+w1*alph3ij(im,jm,km)
                                 vran(m)=u1*bet1ij(im,jm,km)+v1*bet2ij(im,jm,km)+w1*bet3ij(im,jm,km)
                                 wran(m)=u1*gam1ij(im,jm,km)+v1*gam2ij(im,jm,km)+w1*gam3ij(im,jm,km)
                                 ivrelch=1
                                 cycle
                              endif
                              xrel=vrel*dt+sqrt(uint**2+vint**2+(wint)**2)*dt
                              if(ibuoyflag.lt.0)then
                                 xrel=xrel+max(abs(ubulk),abs(vbulk),abs(dhdta))*dt
                                 if(t-tstrt(m).lt.deltaup)xrel=xrel+vout*dt
                              endif
                              if(ibuoyflag .eq. 1 &
                                    .and. inextgridflag .lt. 2)xrel=xrel+wbulkz(m)*dt+.5*forcez*dt**2+5.*sigwb*dt
! mdw 3/17/04 add in mean winds and settling velocity for courant check
! we redo the step with a lower dt if the random velocities are too large
                              if(xrel+abs(uint)*dt.gt.0.7*dx.or.5.*sqrt(coeps*dt)*dt.gt.0.7*dx)then
! mdw 9-07-2004 avoids accidentally drawing too large a random number
                                 dt=min(dt/2.,tre)
                                 cycle
                              endif
                              if(xrel+abs(vint)*dt.gt.0.7*dy.or.5.*sqrt(coeps*dt)*dt.gt.0.7*dy)then
! mdw 9-07-2004 avoids accidentally drawing too large a random number
                                 dt=min(dt/2.,tre)
                                 cycle
                              endif
                              if(xrel+abs(wint)*dt.gt.0.7*dz.or.5.*sqrt(coeps*dt)*dt.gt.0.7*dz)then
! mdw 9-07-2004 avoids accidentally drawing too large a random number
                                 dt=min(dt/2.,tre)
                                 cycle
                              else
                                 exit
                              endif
                           enddo   lp007
                           du=dudet+duran
                           dv=dvdet+dvran
                           dw=dwdet+dwran
                           u1p=uran(m)
                           v1p=vran(m)
                           w1p=wran(m)
                           uran(m)=(u1+du)*alph1ij(im,jm,km)+(v1+dv)*alph2ij(im,jm,km)+(w1+dw)*alph3ij(im,jm,km)
                           vran(m)=(u1+du)*bet1ij(im,jm,km)+(v1+dv)*bet2ij(im,jm,km)+(w1+dw)*bet3ij(im,jm,km)
                           wran(m)=(u1+du)*gam1ij(im,jm,km)+(v1+dv)*gam2ij(im,jm,km)+(w1+dw)*gam3ij(im,jm,km)
                           up=uran(m)
                           vp=vran(m)
                           wp=wran(m)
                           xr=.5*(up+u1p)*dt
                           yr=.5*(vp+v1p)*dt
                           wr=.5*(wp+w1p)*dt
                           dxx=dt*uint*(1.-densegas_factor)+xr
                           dyy=dt*vint*(1.-densegas_factor)+yr
                           dzz=dt*wint*(1.-densegas_factor)+wr
                           if(ibuoyflag .eq. 1 .and. inextgridflag .lt. 2 .and. &
                                 activepartt .gt. 10*number_cutoff)then
                              dzz=dzz+wbulkz(m)*dt+.5*forcez*dt**2+wranb*dt  !+went*dt
                              dyy=dyy+vranb*dt                               ! +vent*dt
                              dxx=dxx+uranb*dt                               !+uent*dt
                           endif
                           if(ibuoyflag.lt.0.and.zp(m).le.hnewa+z0)then
                              dzz=dzz+zran+dhdta*dt
                              dxxb=ubulk*dt+xran
                              dyyb=vbulk*dt+yran
                              dxx=dxx+dxxb
                              dyy=dyy+dyyb                              
!                              if(t-tstrt(m).le.deltaup.and.inextgridflag.lt.2)then
!                                 dzz=dzz+vlocz(m)*dt
!                                 dxx=dxx+vlocx(m)*dt
!                                 dyy=dyy+vlocy(m)*dt
!                              endif
                           endif
                           xp(m)=xp(m)+dxx
                           yp(m)=yp(m)+dyy
                           zp(m)=zp(m)+dzz
                           if(ibuoyflag .eq. 1 .and. inextgridflag .lt. 2)wbulkz(m)=wbulkz(m)+forcez*dt 
                           call reflect
                        endif
!mdw 8-11-2006 removed old formulation for reflection at surface and replaced
!it with the subroutine
! MAN 9/13/2005 binary nextgrid
!                        if(xp(m).le.0.or.xp(m).ge.xpma.or.yp(m).le.0.or.yp(m).ge.ypma.or.zp(m).le.0.or.zp(m).ge.zpma)then
! MAN 01/29/2007 fixed limits in particle domain bounds
                        if(xp(m) .le. z0 .or. xp(m) .ge. xpma .or. &
                              yp(m) .le. z0 .or. yp(m) .ge. ypma .or. &
                              zp(m) .le. 0. .or. zp(m) .ge. zpma)then
                           if(zp(m) .ge. 0. .and. inextgridflag .ge. 1)then
                              if(format_flag.eq.1.or.format_flag.eq.3)then
                                 if(ibuoyflag .eq. -1 .and. idenseflag .eq.1)then
                                    if(format_flag.eq.1)idenseflag=0
                                    do isource=1,nsources
                                       write(31,337)hnew(isource),rnew(isource),xscen(isource),yscen(isource)
                                    enddo
                                 endif
                                 write(31,61)t,xp(m),yp(m),zp(m),pmass(m),weight(m),weightd(m),dpm(m),source_num(m)
                              endif
                              if(format_flag.eq.2.or.format_flag.eq.3)then
                                 if(ibuoyflag .eq. -1 .and. idenseflag .eq.1)then
                                    idenseflag=0
                                    do isource=1,nsources
                                       write(32)hnew(isource),rnew(isource),xscen(isource),yscen(isource)
                                    enddo
                                 endif
                                 write(32)t,xp(m),yp(m),zp(m),pmass(m),weight(m),weightd(m),dpm(m),source_num(m)
                              endif
                           endif
                           xp(m)=0.
                           yp(m)=0.
                           zp(m)=-0.5*dz
                        endif
!end MAN 9/13/2005
!mdw 4-22-2004 write t, position and weights to nextgrid
   337                  format(4f10.3)
    61                  format(8f10.3,1f10.0)
! MAN 01/29/2007 removed second identical if statement
!                        if(xp(m).le.0.or.xp(m).ge.xpma.or.yp(m).le.0.or.yp(m).ge.ypma.or.zp(m).le.0.or.zp(m).ge.zpma)then
!                           xp(m)=0.
!                           yp(m)=0.
!                           zp(m)=-0.5*dz
!                        endif
                        im=int(xp(m)/dx)+1
                        jm=int(yp(m)/dy)+1
                        km=int((zp(m)+dz)/dz)+1
!mdw 7-08-2004 recalculate im, jm, km to be consistent with current position
                        tre=tre-dt
                        if(tre.le.delta*1.e-04)then
                           if(ibuoyflag .eq. 1 .and. inextgridflag .lt. 2)then
                              wbulkzs=wbulkzs+wbulkz(m)
                              sumwbulkz=sumwbulkz+1
! we compute sums of buoyancy force velocities at the end of the timestep        
                           endif
                           exit
                        endif
                        dt=tre
                     enddo
! we calculate deposition for particles still in the active grid
! mdw 4-15-2004 adding coding to update particle to wall distances
                     dxc=xp(m)-.5*dx-dx*real(im-1)
                     dyc=yp(m)-.5*dy-dy*real(jm-1)
                     dzc=zp(m)-.5*dz-dz*real(km-2)
                     if(zp(m).ge.0)then
                        imm=int(xp(m)/dx)+1
                        jmm=int(yp(m)/dy)+1
                        kmm=int((zp(m)+dz)/dz)+1
                        if(ibuoyflag .eq. 1 .and. t .gt. 0.1*delta &
                              .and. inextgridflag .lt. 2 .and. activepartt .gt. 10*number_cutoff)then
                           kbpuf=int((zp(m)-zbminp)/dzbpuf)+1
                           kbpuf=min(10,kbpuf)
                           kbpuf=max(kbpuf,1)
                           xsum(kbpuf)=xsum(kbpuf)+xp(m)
                           ysum(kbpuf)=ysum(kbpuf)+yp(m)
                           xsqsum(kbpuf)=xsqsum(kbpuf)+xp(m)**2
                           ysqsum(kbpuf)=ysqsum(kbpuf)+yp(m)**2
                           zsum=zsum+zp(m)
                           zsqsum=zsqsum+zp(m)**2
                           activepart(kbpuf)=activepart(kbpuf)+1
                           if(zp(m).gt.zbmax)zbmax=zp(m)
                           if(zp(m).gt.0.and.zp(m).lt.zbmin)zbmin=zp(m)
                        endif
                     endif
                     dwxp=dxp(im,jm,km)-dxc
                     dwxm=dxm(im,jm,km)+dxc
                     dwyp=dyp(im,jm,km)-dyc
                     dwym=dym(im,jm,km)+dyc
                     dfzm=dzm(im,jm,km)+dzc
                     if(idepositionflag.eq.1)then
                        if(dxp(im,jm,km).lt.dx.or.dxm(im,jm,km).lt.dx.or.dyp(im,jm,km).lt.dy.or.dym(im,jm,km).lt.dy.or. &
                                dzm(im,jm,km).lt.dz)then
                           taup=taud*1.e+04*ustar**2/xnu
                           if(dpart.gt.0.)then
                              if(taup.gt..10)then
                                 rb=1./((sc**(-2./3.)+10.**(-3./taup))*ustar)
                              else
                                 rb=1./((sc**(-2./3))*ustar)
                              endif
                           else
                              rb=5.*sc**(2./3.)/ustar
                           endif
! mdw 4-15-2004 initialize dnear
                           dnear=max(dx,dy,dz)
                           if(dxp(im,jm,km).lt.dx)then
                              ra=log((.5*dx*ustar*1.e+04/xnu+diff/xnu)/(ustar*1.e+02/xnu+diff/xnu))/(kkar*ustar)
                              vd=1./(ra+rb)
!mdw altered next line to use -99 as flag for internal vd calculation 2/25/2004
                              if(depvel.ne.-99.)vd=.01*depvel
                              depfrac=exp(-vd*2.*delta/dx)
                              depcx(im,jm,km)=depcx(im,jm,km)+.5*(1.0-depfrac)*pmass(m)*weight(m)/(dy*dz)
                              if(dpm(m).le.maxdiameter .and. dpm(m) .ge. mindiameter &
                                    .and. isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
                                 depcxSP(im,jm,km)=depcxSP(im,jm,km)+.5*(1.0-depfrac)*pmass(m)*weight(m)/(dy*dz)
                              endif
                              weight(m)=weight(m)*depfrac
! mdw 4-15-2004 added coding for nearest wall during deposition
                              dnear=min(dnear,dwxp)
                              if(dnear.eq.dwxp)then
                                 xnear(m)=xi(im)+.5*dx
                                 ynear(m)=yp(m)
                                 znear(m)=zp(m)
                              endif
                           endif
                           if(dxm(im,jm,km).lt.dx)then
                              ra=log((.5*dx*ustar*1.e+04/xnu+diff/xnu)/(ustar*1.e+02/xnu+diff/xnu))/(kkar*ustar)
                              vd=1./(ra+rb)
!mdw altered next line to use -99 as flag for internal vd calculation 2/25/2004
                              if(depvel.ne.-99.)vd=.01*depvel
! mdw 4-15-2004 added coding for nearest wall during deposition
                              depfrac=exp(-vd*2.*delta/dx)
                              depcx(im,jm,km)=depcx(im,jm,km)+.5*(1.0-depfrac)*pmass(m)*weight(m)/(dy*dz)
                              if(dpm(m).le.maxdiameter .and. dpm(m) .ge. mindiameter &
                                    .and. isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
                                 depcxSP(im,jm,km)=depcxSP(im,jm,km)+.5*(1.0-depfrac)*pmass(m)*weight(m)/(dy*dz)
                              endif
                              weight(m)=weight(m)*depfrac
                              dnear=min(dnear,dwxm)
                              if(dnear.eq.dwxm)then
                                 xnear(m)=xi(im)-.5*dx
                                 ynear(m)=yp(m)
                                 znear(m)=zp(m)
                              endif
                           endif
                           if(dyp(im,jm,km).lt.dy)then
                              ra=log((.5*dy*ustar*1.e+04/xnu+diff/xnu)/(ustar*1.e+02/xnu+diff/xnu))/(kkar*ustar)
                              vd=1./(ra+rb)
!mdw altered next line to use -99 as flag for internal vd calculation 2/25/2004
                              if(depvel.ne.-99.)vd=.01*depvel
! mdw 4-15-2004 added coding for nearest wall during deposition
                              depfrac=exp(-vd*2.*delta/dy)
                              depcy(im,jm,km)=depcy(im,jm,km)+.5*(1.0-depfrac)*pmass(m)*weight(m)/(dx*dz)
                              if(dpm(m).le.maxdiameter .and. dpm(m) .ge. mindiameter &
                                    .and. isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
                                 depcySP(im,jm,km)=depcySP(im,jm,km)+.5*(1.0-depfrac)*pmass(m)*weight(m)/(dx*dz)
                              endif
                              weight(m)=weight(m)*depfrac
                              dnear=min(dnear,dwyp)
                              if(dnear.eq.dwyp)then
                                 xnear(m)=xp(m)
                                 ynear(m)=yi(jm)+.5*dy
                                 znear(m)=zp(m)
                              endif
                           endif
                           if(dym(im,jm,km).lt.dy)then
                              ra=log((.5*dy*ustar*1.e+04/xnu+diff/xnu)/(ustar*1.e+02/xnu+diff/xnu))/(kkar*ustar)
                              vd=1./(ra+rb)
!mdw altered next line to use -99 as flag for internal vd calculation 2/25/2004
                              if(depvel.ne.-99.)vd=.01*depvel
                              depfrac=exp(-vd*2.*delta/dy)
                              depcy(im,jm,km)=depcy(im,jm,km)+.5*(1.0-depfrac)*pmass(m)*weight(m)/(dx*dz)
                              if(dpm(m).le.maxdiameter .and. dpm(m) .ge. mindiameter &
                                    .and. isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
                                 depcySP(im,jm,km)=depcySP(im,jm,km)+.5*(1.0-depfrac)*pmass(m)*weight(m)/(dx*dz)
                              endif
                              weight(m)=weight(m)*depfrac
                              dnear=min(dnear,dwym)
! mdw 4-15-2004 added coding for nearest wall during deposition
                              if(dnear.eq.dwym)then
                                 xnear(m)=xp(m)
                                 ynear(m)=yi(jm)-.5*dy
                                 znear(m)=zp(m)
                              endif
                           endif
                           if(dzm(im,jm,km).lt.dz)then
                              ra=log((.5*dz*ustar*1.e+04/xnu+diff/xnu)/(ustar*1.e+02/xnu+diff/xnu))/(kkar*ustar)
                              vd=1./(ra+rb+.01*ra*rb*vs)+.01*vs
!mdw altered next line to use -99 as flag for internal vd calculation 2/25/2004
                              if(depvel.ne.-99.)vd=.01*depvel
                              depfrac=exp(-vd*2.*delta/dz)
                              depcz(im,jm,km)=depcz(im,jm,km)+.5*(1.0-depfrac)*pmass(m)*weight(m)/(dy*dx)
                              if(dpm(m).le.maxdiameter .and. dpm(m) .ge. mindiameter &
                                    .and. isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
                                 depczSP(im,jm,km)=depczSP(im,jm,km)+.5*(1.0-depfrac)*pmass(m)*weight(m)/(dx*dy)
                              endif
                              weight(m)=weight(m)*depfrac
! mdw 4-15-2004 added coding for nearest wall during deposition
                              dnear=min(dnear,dfzm)
                              if(dfzm.eq.dnear)then
                                 xnear(m)=xp(m)
                                 ynear(m)=yp(m)
                                 znear(m)=zi(km)-.5*dz
                              endif
                           endif
                           if(dzp(im,jm,km).lt.dz)then
                              ra=log((.5*dz*ustar*1.e+04/xnu+diff/xnu)/(ustar*1.e+02/xnu+diff/xnu))/(kkar*ustar)
                              vd=1./(ra+rb+.01*ra*rb*vs)+.01*vs
!mdw altered next line to use -99 as flag for internal vd calculation 2/25/2004
                              if(depvel.ne.-99.)vd=.01*depvel
                              depfrac=exp(-vd*2.*delta/dz)
                              depcz(im,jm,km)=depcz(im,jm,km)+.5*(1.0-depfrac)*pmass(m)*weight(m)/(dy*dx)
                              if(dpm(m).le.maxdiameter .and. dpm(m) .ge. mindiameter &
                                    .and. isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
                                 depczSP(im,jm,km)=depczSP(im,jm,km)+.5*(1.0-depfrac)*pmass(m)*weight(m)/(dx*dy)
                              endif
                              weight(m)=weight(m)*depfrac
! mdw 4-15-2004 added coding for nearest wall during deposition
! mdw 6-8-2004 changed deposition positions to particle position and cell face
                              dnear=min(dnear,drzp)
                              if(drzp.eq.dnear)then
                                 xnear(m)=xp(m)
                                 ynear(m)=yp(m)
                                 znear(m)=zi(km)+.5*dz
                              endif
                           endif
                        endif
                        weightd(m)=exp(-decay*(t+delta-tstrt(m)))
                     endif
! mdw coding for average concentrations for infiltration  1-22-2004
                     if(iindoorflag.eq.1)then
                        if(dxp(im,jm,km).lt.dx.or.dxm(im,jm,km).lt.dx.or.dyp(im,jm,km).lt.dy.or.dym(im,jm,km).lt.dy.or. &
                           dzm(im,jm,km).lt.dz)then                           
                           if(dxp(im,jm,km).lt.dx)then
                              do ibuild=1,inumbuild
                                 if(bldtype(ibuild) .eq. 9)cycle
                                 ixfo=int(xfo(ibuild)/dx)+1
                                 jyfo=int(yfo(ibuild)/dy)+1
                                 dxfo=dx*real(im-1)+.5*dx-xfo(ibuild)
                                 dyfo=dy*real(jm-1)+.5*dy-yfo(ibuild)
                                 ximrot=dxfo*cos(gamma(ibuild)*pi/180)+dyfo*sin(gamma(ibuild)*pi/180)
                                 yjmrot=dyfo*cos(gamma(ibuild)*pi/180)-dxfo*sin(gamma(ibuild)*pi/180)
                                 imrot=ixfo+int(ximrot/dx)
                                 jmrot=jyfo+int(yjmrot/dy)
                                 if(km.le.kface2n(ibuild).and.km.ge.kface1n(ibuild))then
                                    if((imrot.ge.iface1n(ibuild)-1).and.imrot.lt.iface2n(ibuild).and.jmrot.ge.jface1n(ibuild).and. &
                                              jmrot.le.jface2n(ibuild))then
                                       if(icellflag(im+1,jm,km) .eq. 0 .and. icellflag(im,jm,km) .ne. 0)then
                                          chioa(ibuild,1)=chioa(ibuild,1)+delta*weight(m)*pmass(m)*weightd(m)&
                                                          /(bvol*deltilf*aminusx(ibuild))
                                          if(dpm(m).le.maxdiameter .and. dpm(m) .ge. mindiameter &
                                                .and. isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
                                             chioaSP(ibuild,1)=chioaSP(ibuild,1)+delta*weight(m)*pmass(m)*weightd(m)&
                                                               /(bvol*deltilf*aminusx(ibuild))
                                          endif
                                       endif
                                    endif
                                 endif
                              enddo
                           endif
                           if(dxm(im,jm,km).lt.dx)then
                              do ibuild=1,inumbuild
! mdw 4-20-2004 added lines to check to see if km is between floor and roof
! for each of the vertical surfaces
                                 if(bldtype(ibuild) .eq. 9)cycle
                                 ixfo=int(xfo(ibuild)/dx)+1
                                 jyfo=int(yfo(ibuild)/dy)+1
                                 dxfo=dx*real(im-1)+.5*dx-xfo(ibuild)
                                 dyfo=dy*real(jm-1)+.5*dy-yfo(ibuild)
                                 ximrot=dxfo*cos(gamma(ibuild)*pi/180)+dyfo*sin(gamma(ibuild)*pi/180)
                                 yjmrot=dyfo*cos(gamma(ibuild)*pi/180)-dxfo*sin(gamma(ibuild)*pi/180)
                                 imrot=ixfo+int(ximrot/dx)
                                 jmrot=jyfo+int(yjmrot/dy)                                 
                                 if(km.le.kface2n(ibuild).and.km.ge.kface1n(ibuild))then
                                    if((imrot.le.iface2n(ibuild)+1).and.imrot.gt.iface1n(ibuild).and.jmrot.ge.jface1n(ibuild).and. &
                                              jmrot.le.jface2n(ibuild))then
                                       if(icellflag(im,jm,km) .ne. 0 .and. icellflag(im-1,jm,km) .eq. 0)then
                                          chioa(ibuild,2)=chioa(ibuild,2)+delta*weight(m)*pmass(m)*weightd(m)&
                                                          /(bvol*deltilf*aplusx(ibuild))
                                          if(dpm(m).le.maxdiameter .and. dpm(m) .ge. mindiameter &
                                                .and. isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
                                             chioaSP(ibuild,2)=chioaSP(ibuild,2)+delta*weight(m)*pmass(m)*weightd(m)&
                                                               /(bvol*deltilf*aminusx(ibuild))
                                          endif
                                       endif
                                    endif
                                 endif
                              enddo
                           endif
                           if(dyp(im,jm,km).lt.dy)then
                              do ibuild=1,inumbuild
                                 if(bldtype(ibuild) .eq. 9)cycle
                                 ixfo=int(xfo(ibuild)/dx)+1
                                 jyfo=int(yfo(ibuild)/dy)+1
                                 dxfo=dx*real(im-1)+.5*dx-xfo(ibuild)
                                 dyfo=dy*real(jm-1)+.5*dy-yfo(ibuild)
                                 ximrot=dxfo*cos(gamma(ibuild)*pi/180)+dyfo*sin(gamma(ibuild)*pi/180)
                                 yjmrot=dyfo*cos(gamma(ibuild)*pi/180)-dxfo*sin(gamma(ibuild)*pi/180)
                                 imrot=ixfo+int(ximrot/dx)
                                 jmrot=jyfo+int(yjmrot/dy)
                                 if(km.le.kface2n(ibuild).and.km.ge.kface1n(ibuild))then
                                    if((jmrot.ge.jface1n(ibuild)-1).and.jmrot.lt.jface2n(ibuild).and.imrot.ge.iface1n(ibuild).and. &
                                              imrot.le.iface2n(ibuild))then
                                       if(icellflag(im,jm,km) .ne. 0 .and. icellflag(im,jm+1,km) .eq. 0)then
                                          chioa(ibuild,3)=chioa(ibuild,3)+delta*weight(m)*pmass(m)*weightd(m)&
                                                          /(bvol*deltilf*aminusy(ibuild))
                                          if(dpm(m).le.maxdiameter .and. dpm(m) .ge. mindiameter &
                                                .and. isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
                                             chioaSP(ibuild,3)=chioaSP(ibuild,3)+delta*weight(m)*pmass(m)*weightd(m)&
                                                               /(bvol*deltilf*aminusx(ibuild))
                                          endif
                                       endif
                                    endif
                                 endif
                              enddo
                           endif
                           if(dym(im,jm,km).lt.dy)then
                              do ibuild=1,inumbuild
                                 if(bldtype(ibuild) .eq. 9)cycle
                                 ixfo=int(xfo(ibuild)/dx)+1
                                 jyfo=int(yfo(ibuild)/dy)+1
                                 dxfo=dx*real(im-1)+.5*dx-xfo(ibuild)
                                 dyfo=dy*real(jm-1)+.5*dy-yfo(ibuild)
                                 ximrot=dxfo*cos(gamma(ibuild)*pi/180)+dyfo*sin(gamma(ibuild)*pi/180)
                                 yjmrot=dyfo*cos(gamma(ibuild)*pi/180)-dxfo*sin(gamma(ibuild)*pi/180)
                                 imrot=ixfo+int(ximrot/dx)
                                 jmrot=jyfo+int(yjmrot/dy)                                 
                                 if(km.le.kface2n(ibuild).and.km.ge.kface1n(ibuild))then
                                    if((jmrot.le.jface2n(ibuild)+1).and.jmrot.gt.jface1n(ibuild).and.imrot.ge.iface1n(ibuild).and. &
                                                 imrot.le.iface2n(ibuild))then
                                       if(icellflag(im,jm,km) .ne. 0 .and. icellflag(im,jm-1,km) .eq. 0)then
                                          chioa(ibuild,4)=chioa(ibuild,4)+delta*weight(m)*pmass(m)*weightd(m)&
                                                          /(bvol*deltilf*aplusy(ibuild))
                                          if(dpm(m).le.maxdiameter .and. dpm(m) .ge. mindiameter &
                                                .and. isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
                                             chioaSP(ibuild,4)=chioaSP(ibuild,4)+delta*weight(m)*pmass(m)*weightd(m)&
                                                               /(bvol*deltilf*aminusx(ibuild))
                                          endif
                                       endif
                                    endif
                                 endif
                              enddo
                           endif
                           if(dzm(im,jm,km).lt.dz)then
                              do ibuild=1,inumbuild
                                 if(bldtype(ibuild) .eq. 9)cycle
                                 ixfo=int(xfo(ibuild)/dx)+1
                                 jyfo=int(yfo(ibuild)/dy)+1
                                 dxfo=dx*real(im-1)+.5*dx-xfo(ibuild)
                                 dyfo=dy*real(jm-1)+.5*dy-yfo(ibuild)
                                 ximrot=dxfo*cos(gamma(ibuild)*pi/180)+dyfo*sin(gamma(ibuild)*pi/180)
                                 yjmrot=dyfo*cos(gamma(ibuild)*pi/180)-dxfo*sin(gamma(ibuild)*pi/180)
                                 imrot=ixfo+int(ximrot/dx)
                                 jmrot=jyfo+int(yjmrot/dy) 
                                 if(km.eq.kface2n(ibuild)+1.and.imrot.ge.iface1n(ibuild).and.imrot.le.iface2n(ibuild).and. &
                                         jmrot.ge.jface1n(ibuild).and.jmrot.le.jface2n(ibuild))then
                                    if(icellflag(im,jm,km) .ne. 0 .and. icellflag(im,jm,km-1) .eq. 0)then
                                       chioa(ibuild,5)=chioa(ibuild,5)+delta*weight(m)*pmass(m)*weightd(m)&
                                                       /(bvol*deltilf*aroof(ibuild))
                                       if(dpm(m).le.maxdiameter .and. dpm(m) .ge. mindiameter &
                                             .and. isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
                                          chioaSP(ibuild,5)=chioaSP(ibuild,5)+delta*weight(m)*pmass(m)*weightd(m)&
                                                            /(bvol*deltilf*aminusx(ibuild))
                                       endif
                                    endif
                                 endif
                              enddo
                           endif
                           if(dzp(im,jm,km).lt.dz)then
                              do ibuild=1,inumbuild
                                 if(bldtype(ibuild) .eq. 9)cycle
                                 ixfo=int(xfo(ibuild)/dx)+1
                                 jyfo=int(yfo(ibuild)/dy)+1
                                 dxfo=dx*real(im-1)+.5*dx-xfo(ibuild)
                                 dyfo=dy*real(jm-1)+.5*dy-yfo(ibuild)
                                 ximrot=dxfo*cos(gamma(ibuild)*pi/180)+dyfo*sin(gamma(ibuild)*pi/180)
                                 yjmrot=dyfo*cos(gamma(ibuild)*pi/180)-dxfo*sin(gamma(ibuild)*pi/180)
                                 imrot=ixfo+int(ximrot/dx)
                                 jmrot=jyfo+int(yjmrot/dy)
                                 if(km.eq.kface2n(ibuild)-1.and.imrot.ge.iface1n(ibuild).and.imrot.le.iface2n(ibuild).and. &
                                          jmrot.ge.jface1n(ibuild).and.jmrot.le.jface2n(ibuild))then
                                    if(icellflag(im,jm,km) .ne. 0 .and. icellflag(im,jm,km+1) .eq. 0)then
                                       chioa(ibuild,6)=chioa(ibuild,6)+delta*weight(m)*pmass(m)*weightd(m)&
                                                       /(bvol*deltilf*afloor(ibuild))
                                       if(dpm(m).le.maxdiameter .and. dpm(m) .ge. mindiameter &
                                             .and. isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
                                          chioaSP(ibuild,6)=chioaSP(ibuild,6)+delta*weight(m)*pmass(m)*weightd(m)&
                                                            /(bvol*deltilf*afloor(ibuild))
                                       endif
                                    endif
                                 endif
                              enddo
                           endif
                        endif
                     endif
! mdw 1-22-2004 end building infiltration - building face concentration calculation
                  enddo   lp006
! we only move particles if they are within the grid, otherwise
! their positions are unchanged, so we can see where they left the
! grid
!mdw 3-08-2004 redefine t to be the time at the end of the timestep
                  t=t+delta !*real(it)
! mdw 7-13-2004 test to see if new winds need to be read in
                  if(t .ge. tupdate .and. t .le. dur+delta_ref .and. time_idx .lt. ntinc)then
!end TMB 1/8/04
! mdw 3/22/2004pm added read statements before and after roofflag read per TMB's instructions
                     time_idx=time_idx+1
                     if(format_flag.eq.1)then
                        read(27,10)head2
                        read(27,*)iwindt
                        read(27,*)twind
                     endif
                     call readwind
                     call readbldgout
                     call turbinit
                     if(inverseflag .eq. 1)then
                        u=-u
                        v=-v
                        w=-w
                     endif
                     tupdate=tupdate+3600.*tinc
                  endif
                  if(inextgridflag .le. 1)mmpt=nptot
                  if(ibuoyflag .eq. 1 .and. t .le. 0.1*delta .and. inextgridflag .lt. 2)then
                     zbmin=source_z(1)
                     zbmax=source_z(1)+.999999*rfire
                     dzbpuf=.1*(zbmax-zbmin)
                     zbcen=zbmin+.5*(zbmax-zbmin)
                     do kbpuf=1,kbpufmax
                        zbpufmin=dzbpuf*real(kbpuf-1)+zbmin
                        zbpufmax=dzbpuf*real(kbpuf)+zbmin
                        zpuftop=min(zbmax,zbpufmax)
                        zpufbot=max(zbmin,zbpufmin)
                        activepart0(kbpuf)=real(nptot)*volpart0(kbpuf)/(2.*pi*rfire**3/3.)
                        volpart0(kbpuf)=pi*(zpuftop-zbmin)*sqrt(rfire**2-(zpuftop-zbmin)**2)&
                                        -pi*(zpufbot-zbmin)*sqrt(rfire**2-(zpufbot-zbmin))**2&
                                        +pi*rfire*asin((zpuftop-zbmin)/rfire)/2.&
                                        -pi*rfire*asin((zpufbot-zbmin)/rfire)/2.
                        activepartp(kbpuf)=activepart0(kbpuf)
                     enddo
                  endif
lp100:            do mm=1,mp
                     if(tstrt(mm) .le. t-0.999*delta .or. tstrt(mm) .gt. t+0.001*delta)then
                        cycle
                     else
! we add new particles if the time is equal to or greater than tstrt
! we start and move any new particles added in the last time
! interval
                        dt=(tstrt(mm)-t) !this will always be zero so we don't need to move part.
                        tre=(tstrt(mm)-t)
                        nstps=0
                        xp(mm)=xps(mm)
                        yp(mm)=yps(mm)
                        zp(mm)=zps(mm)
                        if(ibuoyflag .eq. 1 .and. inextgridflag .lt. 2)then
                           temp_part(mm)=tfire
                           zpp(mm)=zp(mm)
                        endif
                        im=int(xp(mm)/dx)+1
                        jm=int(yp(mm)/dy)+1
                        km=int((zp(mm)+dz)/dz)+1
!mdw 11/29/2005 adjustments for positive buoyancy
                        if(ibuoyflag .eq. 1 .and. inextgridflag .lt. 2)then
                           kbpuf=int((zp(mm)-zbmin)/dzbpuf)+1
                           kbpuf=max(1,kbpuf)
                           kbpuf=min(10,kbpuf)
                           wbulkz(mm)=0.
                           xsum(kbpuf)=xsum(kbpuf)+xp(mm)
                           ysum(kbpuf)=ysum(kbpuf)+yp(mm)
                           xsqsum(kbpuf)=xsqsum(kbpuf)+xp(mm)**2
                           ysqsum(kbpuf)=ysqsum(kbpuf)+yp(mm)**2
                           zsum=zsum+zp(mm)
                           zsqsum=zsqsum+zp(mm)**2
                           activepart(kbpuf)=activepart(kbpuf)+1
                           activepartt=activepartt+1
                        endif
                        dxc=xp(mm)-.5*dx-dx*real(im-1)
                        dyc=yp(mm)-.5*dy-dy*real(jm-1)
                        dzc=zp(mm)-.5*dz-dz*real(km-2)
! these are random components with u rotated to the wind direction
                        if((im.ge.1.and.im.le.nx-1).and.(jm.ge.1.and.jm.le.ny-1).and.(km.gt.1.and.km.le.nz-1))then
                           dwxp=dxp(im,jm,km)-dxc
                           dwxm=dxm(im,jm,km)+dxc
                           dwyp=dyp(im,jm,km)-dyc
                           dwym=dym(im,jm,km)+dyc
                           dfzm=dzm(im,jm,km)+dzc
                           drzp=dzp(im,jm,km)-dzc
                           dt=(tstrt(mm)-t)
                           ustar=ustarij(im,jm,km)
                           sigu=sigui(im,jm,km)
                           sigv=sigvi(im,jm,km)
                           sigw=sigwi(im,jm,km)
                           siguo=sigui(im,jm,km)
                           sigvo=sigvi(im,jm,km)
                           sigwo=sigwi(im,jm,km)
                           dsigwdz=dsigwdzi(im,jm,km)
                           call rangen(x)
                           w1=sigw*x
                           call rangen(x)
                           u1=sigu*x
                           call rangen(x)
                           v1=sigv*x
                           uran(mm)=u1*alph1ij(im,jm,km)+v1*alph2ij(im,jm,km)+w1*alph3ij(im,jm,km)
                           vran(mm)=u1*bet1ij(im,jm,km)+v1*bet2ij(im,jm,km)+w1*bet3ij(im,jm,km)
                           wran(mm)=u1*gam1ij(im,jm,km)+v1*gam2ij(im,jm,km)+w1*gam3ij(im,jm,km)
                           mmm=mm
                        endif
                     endif
                  enddo lp100
                  if(ibuoyflag .lt. 0 .and. t .le. 0.1*delta .and. source_release(1) .ge. 2)then
                     isource=1
                     call bulkvelc
                  endif
                  if(ibuoyflag .lt. 0 .and. source_release(1) .eq. 1.)then
!                     icellflag_dense(:,:)=0
                     delta_old=delta
                     delta=delta_ref
                     do isource=1,nsources
                         if(t .le. 0.1*delta)then
                            hnew(isource)=source_h(isource)
                            vol0(isource)=pi*source_h(isource)*source_r(isource)**2
                         endif                      
                         drdt(isource)=sqrt((9.8*(source_density(isource)-rhoair)/rhoair)*vol0(isource))/rnew(isource)
                         dhdt(isource)=2.*hnew(isource)*drdt(isource)/rnew(isource)
                         if(drdt(isource) .gt. .5*min(dx,dy) .or. &
                               drdt(isource) .gt. (.25*dz*rnew(isource)/hnew(isource)))then
                            delta=min(.5*min(dx,dy)/drdt(isource),0.25*dz*rnew(isource)&
                                  /(hnew(isource)*drdt(isource)),delta_ref,delta)
                         endif
                     enddo
                     delta_ref_rm=delta_ref_rm-delta_old
                     if(1.1*delta .ge. delta_ref_rm .and. delta_ref_rm .gt. 0)then
                        delta=delta_ref_rm
                     elseif(delta_ref_rm .le. 0.1*delta)then
                        delta_ref_rm=delta_ref
                     endif
                     do isourceindex=1,nsources
                        isource=isourceindex
!                        isource=mod(isource,nsources)+1
                        if(t .le. 0.1*delta)then
                           call bulkveli
                        else
                          if(rhocent(isource).gt.1.001*rhoair)call bulkveli
                        endif
                     enddo
!                        isourcestart=isourcestart+1
!                        isourcestart=mod(isourcestart,nsources)
                  endif

!mdw 11/29/2005 added buoyant gas formulation
                  if(ibuoyflag .eq. 1 .and. inextgridflag .lt. 2)then
                     dzbpuf=.1*(zbmax-zbmin)
                     activepartt=0.
                     vlsumtp=vlsumt
                     vlsumt=0.
                     contotp=contot
                     conmaxp=conmax
                     conmax=0.
                     do kz=1,kbpufmax
                        volpartp(kz)=volpart(kz)
                        volpart(kz)=0.
                        zirhoht=0.
                        zirholt=h
                        iflagden=0
                        zkmin=zbmin+dzbpuf*real(kz-1)
                        zkmax=zbmin+dzbpuf*real(kz)
                        k1=int(zkmin/dz)+2
                        k2=int(zkmax)+2
                        k2=min(nz-1,k2)
                        zbcen=.5*(zkmax+zkmin)
                        do im=1,nx-1
                           do jm=1,ny-1
                              if(kz.eq.1)then
                                 zirhoh(im,jm)=0.
                                 zirhol(im,jm)=h
                              endif
                              do km=k1,k2
                                 dzk=dz
                                 if(km.eq.k1)dzk=dz-(zkmin-dz*real(km-2))
                                 if(km.eq.k2)dzk=dzk-(dz*real(km-1)-zkmax)
                                 tfirez=tfire-.01*zi(km)
                                 denfirez=.0288*273/(.0224*tfirez)
                                 if(zirhoh(im,jm).gt.zirhoht)zirhoht=zirhoh(im,jm)
                                 if(zirhol(im,jm).lt.zirholt)zirholt=zirhol(im,jm)
                              enddo
                              if(zirhoh(im,jm).gt.zirhoht)zirhoht=zirhoh(im,jm)
                              if(zirhol(im,jm).lt.zirholt)zirholt=zirhol(im,jm)
                           enddo
                        enddo
                        activepartt=activepartt+activepart(kz)
                        if(activepart(kz).gt.number_cutoff)then
                           xic(kz)=xsum(kz)/activepart(kz)
                           yic(kz)=ysum(kz)/activepart(kz)
                           xsqic(kz)=xsqsum(kz)/activepart(kz)
                           ysqic(kz)=ysqsum(kz)/activepart(kz)
                           varx(kz)=(xsqsum(kz)-activepart(kz)*xic(kz)**2)/(activepart(kz)-1)
                           sigxb(kz)=sqrt(varx(kz))
                           vary(kz)=(ysqsum(kz)-activepart(kz)*yic(kz)**2)/(activepart(kz)-1)
                           sigyb(kz)=sqrt(vary(kz))
                           vlsum(kz)=dzbpuf*sigxb(kz)*sigyb(kz)
                           vlsumt=vlsumt+vlsum(kz)
                           areat=max(areat,sigxb(kz)*sigyb(kz))
                           if(kz.gt.kzmax)kzmax=kz
                           if(kz.lt.kzmin)kzmin=kz
                           conl(kz)=activepart(kz)/vlsum(kz)
                           if(conl(kz).gt.conmax)conmax=conl(kz)
                           zirholt=zirholt-2.*dz
                           zirholt=max(0.,zirholt)
                           zirhoht=zirhoht+2.*dz
                           zirhoht=min(h,zirhoht)
                        delrho0=(rhobkg(2)-denfire)/rhobkg(2)
                        endif
                        activepartp(kz)=activepart(kz)
                        activepart(kz)=0.
                        xsum(kz)=0.
                        ysum(kz)=0.
                        xsqsum(kz)=0.
                        ysqsum(kz)=0.
                     enddo
                  endif
                  contotp=contot
                  if(ibuoyflag .gt. 0 .and. inextgridflag .lt. 2 .and. &
                        activepartt .gt. 10*number_cutoff)then
                     contot=activepartt/vlsumt
                     zic=zsum/activepartt
                     zsqic=zsqsum/(activepartt-1)
                     varz=(zsqsum-activepartt*zic**2)/(activepartt-1)
                     sigzb=sqrt(varz)
                     vlsumt=sigzb*areat
                     contot=activepartt/vlsumt
                  endif
                  zsum=0.
                  zsqsum=0.
                  if(t .le. 0.001*delta)then !MAN 01/29/2007
                     if(format_flag.eq.1 .or. format_flag.eq.3)then ! man 7/13/2005 binary data files
                        ntb=nbx*nby*nbz
                        if(isiteflag.eq.0)then
                           write(10,120)
                           write(10,130)
! mdw 4-19-2004 new variables for partfield for deposition points
                           write(11,1120)
                           if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8) write(41,1120)
                           write(15,11202)
                           if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8) write(45,11202)
                           write(15,11302)
                           if(ieradflag.eq.1) write(45,11302)
                           if(ilctflag.eq.1) write(14,11201)
                           write(11,1130)
                           if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8) write(41,1130)
                           if(ilctflag.eq.1) write(14,11301)
                           if(idepositionflag.eq.1) write(19,1127)
                           if(idepositionflag.eq.1) write(19,1137)
                        else
                           write(21,11202)
                           write(21,11302)
                        endif
        120             format('TITLE     = "Particles in 3D"')
        130             format('VARIABLES = "X" "Y" "Z" "UP" "VP" "WP" "MASS" "WT" "WTD" "DPM" "XN" "YN" "ZN" "SRCNUM"')
       1120             format('TITLE     = "Concentrations in 3D (grams per cubic meter)"')
      11202             format('TITLE     = "Dosages in 3D (gram seconds per cubic meter)"')
      11302             format('VARIABLES = "X" "Y" "Z" "Dosage"')
      11201             format('TITLE     = "LCTs in 3D"')
       1130             format('VARIABLES = "X" "Y" "Z" "C"') 
      11301             format('VARIABLES = "X" "Y" "Z" "D"')
       1127             format('TITLE     = "XY plane Deposition in 3D"')
       1137             format('VARIABLES = "X" "Y" "Z" "depcz"')
                     endif ! man 7/13/2005 binary data files
                  endif
                  if(tout + 0.5*delta .ge. outp .or. t .le. 0.5*delta)then ! MAN 01/29/2007
                     if(format_flag.eq.1 .or. format_flag.eq.3)then ! man 7/13/2005 binary data files
                        if(isiteflag.eq.0)then
                           write(10,140)t
                           write(10,150)mp-ninactiveparticles
                           write(10,160)
                        endif
! mdw 3/08/2004 changed format to f8.1 from f5.1 to write out larger times
        140             format('ZONE T="',f8.1,' seconds"')
        150             format('I=',i10,', J=1, K=1,F=POINT')
                     endif ! man 7/13/2005 binary data files
                     if((format_flag.eq.2.or.format_flag.eq.3 ).and.&
                           isiteflag .eq. 0) write(54)mp-ninactiveparticles
                  endif
   160            format('DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)')
 16001            format(6f10.3,E12.4)
 16003            format(6f10.3,1e12.4,2f10.7,f8.3,3f10.3,1f10.0)
! mdw 4-19-2004 changed format and write statements for xnear, ynear, and znear
lp011:            do mm=1,mp
                     if(tstrt(mm) .gt. t+0.001*delta)cycle
                     if((tout + 0.5*delta .ge. outp .or. t .le. 0.5*delta ).and.&
                           isiteflag .eq. 0)then
                        if(format_flag.eq.1.or.format_flag.eq.3)then
                           write(10,16003)xp(mm),yp(mm),zp(mm),uran(mm),vran(mm),wran(mm),pmass(mm),&
                                          weight(mm),weightd(mm),dpm(mm),xnear(mm),ynear(mm),znear(mm),source_num(mm)
                        endif
                        if(format_flag.eq.2.or.format_flag.eq.3)then
                           write(54)xp(mm),yp(mm),zp(mm),xnear(mm),ynear(mm),znear(mm),&
                                    pmass(mm)*weight(mm)*weightd(mm),dpm(mm),source_num(mm)
                        endif
                     endif
! MAN 9/13/2005 binary nextgrid
                     if(t .ge. dur .and. zp(mm) .gt. 0. .and. inextgridflag .ge. 1)then
                        if(format_flag.eq.1.or.format_flag.eq.3)then
                           write(31,61)t,xp(mm),yp(mm),zp(mm),pmass(mm),weight(mm),weightd(mm),dpm(mm),source_num(mm)
                        endif
                        if(format_flag.eq.2.or.format_flag.eq.3)then
                           write(32)t,xp(mm),yp(mm),zp(mm),pmass(mm),weight(mm),weightd(mm),dpm(mm),source_num(mm)
                        endif
                     endif
! end MAN 9/13/2005
                     xpmax=max(xpmax,xp(mm))
                     xpmin=min(xpmin,xp(mm))
                     ypmax=max(ypmax,yp(mm))
                     ypmin=min(ypmin,yp(mm))
                     zpmax=max(zpmax,zp(mm))
                     zpmin=min(zpmin,zp(mm))
                     ib=int((xp(mm)-xbl)/dbx)+1
                     if(xp(mm).lt.xbl)ib=0
                     jb=int((yp(mm)-ybl)/dby)+1
                     if(yp(mm).lt.ybl)jb=0
                     if(isiteflag.eq.0)then
                        kb=int((zp(mm)-zbl)/dbz)+1
                        if(zp(mm).lt.zbl)kb=0
                        if(ib.ge.1.and.ib.le.nbx.and.jb.ge.1.and.jb.le.nby&
                              .and.kb.ge.1.and.kb.le.nbz.and.t.gt.concstrt)then !MAN 01/29/2007 t.ge.concstrt
                           oldcon=con(ib,jb,kb)
                           con(ib,jb,kb)=con(ib,jb,kb)+delta*weight(mm)*pmass(mm)*weightd(mm)/(bvol*outc)
!Multiply the current concentration by the timestep and add to the dosage
                           dose(ib,jb,kb)=dose(ib,jb,kb)+delta*weight(mm)*pmass(mm)*weightd(mm)/bvol
                           if(dpm(mm).le.maxdiameter .and. dpm(mm) .ge. mindiameter &
                                 .and. isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
                              conSP(ib,jb,kb)=conSP(ib,jb,kb)+delta*weight(mm)*pmass(mm)*weightd(mm)/(bvol*outc)
!get respirable con and dose
                              doseSP(ib,jb,kb)=doseSP(ib,jb,kb)+delta*weight(mm)*pmass(mm)*weightd(mm)/bvol
                           endif
!MAN 01/26/2007 conditional particle splitting
!Splits if particle enters a cell that has no concentration for longer than delta after the beginning of the concentration averaging time
!Splits will only occur if the maximum number of particles has not been reached
!Splits will not occur in the extra time step that is used to enure that the concentrations are written out
                           if(oldcon .le. 0. .and. toutc .gt. delta .and. nptot .lt. mpart-1 .and. t .le. dur)then
!This ensures that it will only split upto the maximum number of particles and no further
                              if(nptot+partsplitfactor-1 .ge. mpart)then
                                 partsplitfactor=mpart-nptot
                              endif
!The mass of the particle is divided equally between the children and all children have the original partilce's last attributes.
                              pmass(mm)=pmass(mm)/real(partsplitfactor)
                              do ippr=nptot+1,nptot+partsplitfactor-1
                                 pmass(ippr)=pmass(mm)
                                 xp(ippr)=xp(mm)
                                 yp(ippr)=yp(mm)
                                 zp(ippr)=zp(mm)
                                 tstrt(ippr)=t ! MAN 01/29/2007 -delta
                                 weight(ippr)=weight(mm)
                                 weightd(ippr)=weightd(mm)
                                 dpm(ippr)=dpm(mm)
                                 uran(ippr)=uran(mm)
                                 vran(ippr)=vran(mm)
                                 wran(ippr)=wran(mm)
                                 uranbp(ippr)=uranbp(mm)
                                 vranbp(ippr)=vranbp(mm)
                                 wranbp(ippr)=wranbp(mm)
                                 xps(ippr)=xps(mm)
                                 yps(ippr)=yps(mm)
                                 zps(ippr)=zps(mm)
                                 xnear(ippr)=xnear(mm)
                                 ynear(ippr)=ynear(mm)
                                 znear(ippr)=znear(mm)
                                 source_num(ippr)=source_num(mm)
                                 if(ibuoyflag .gt. 0 .and. inextgridflag .lt. 2)then
                                    wbulkz(ippr)=wbulkz(mm)
                                    zpp(ippr)=zpp(mm)
                                    temp_part(ippr)=temp_part(mm)
                                 endif
                              enddo
                              nptot=nptot+partsplitfactor-1 !the total number of particles is increased
                           endif
! end MAN 01/26/2007
                        endif
                     else
                        immm=int(xp(mm)/dx)+1
                        jmmm=int(yp(mm)/dy)+1
                        kb=int((zp(mm)-hgt(immm,jmmm)-zbl)/dbz)+1
                        if(zp(mm)-hgt(immm,jmmm).lt.zbl)kb=0
                        if(ib.ge.1.and.ib.le.nbx.and.jb.ge.1.and.&
                              jb.le.nby.and.kb.eq.1.and.t.ge.concstrt)then
                           dosexy(ib,jb)=dosexy(ib,jb)+pmass(mm)*delta*weight(mm)*weightd(mm)/bvol
                        endif
                     endif

                  enddo   lp011
                  mm=0
                  if(tout+0.5*delta .ge. outp)then !MAN 01/29/2007
                     tout=tout-outp
                  endif
                  if(iindoorflag.eq.1)then
                     if(touti.ge.deltilf-0.001*delta)then
                        ioutin=ioutin+1
                        if(ioutin.eq.1)then
!                           write(28,8020)
!                           write(28,8021)
                           noutbldg=int(dur/deltilf)
                           write(28,8020)noutbldg
                           write(28,8021)inumbuild-inumveg
                        endif
      8020              format(i7,'  !Number of infiltration time steps')
      8021              format(i7,'  !Number of buildings')
      8022              format(f8.0,' !Time (s)')
      8023              format('!Build#, C(g/m^3), D(g*s/m^3), Cex(g/m^3), Csp(g/m^3), Dsp(g*s/m^3), Csp(g/m^3)')
                        toutin=deltilf*real(ioutin)
                        touti=touti-deltilf
!                        write(28,8022)toutin,inumbuild
                        write(28,8022)toutin
                        write(28,8023)
lp012:                  do ibuild=1,inumbuild
                           if(bldtype(ibuild) .eq. 9)cycle
                           do iface=1,6
                              chiosum(ibuild,ioutin)=chiosum(ibuild,ioutin)+chioa(ibuild,iface)
                              chioa(ibuild,iface)=0.
                              if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
                                 chiosumSP(ibuild,ioutin)=chiosumSP(ibuild,ioutin)+chioaSP(ibuild,iface)
                                 chioaSP(ibuild,iface)=0.
                              endif
                           enddo
                           if(ioutin.eq.1)then
!                              chiin(ibuild,ioutin)=chiosum(ibuild,ioutin)-chiosum(ibuild,ioutin) &
!                                                   *(1.-exp(-deltilf/tau(ibuild)))*tau(ibuild)/deltilf
                              if(chiin(ibuild,ioutin) .eq. 0.)then
                                 chiin(ibuild,ioutin)=abinf(ibuild)*chiosum(ibuild,ioutin)&
                                                   -chiosum(ibuild,ioutin)*(1.-exp(-binf(ibuild)*deltilf))&
                                                   *abinf(ibuild)/(binf(ibuild)*deltilf)
                              endif
                              chiex(ibuild,ioutin)=0
                              dosin(ibuild,ioutin)=deltilf*chiin(ibuild,ioutin)
                              if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
                                 if(chiinSP(ibuild,ioutin) .eq. 0.)then
                                    chiinSP(ibuild,ioutin)=abinf(ibuild)*chiosumSP(ibuild,ioutin)&
                                                   -chiosumSP(ibuild,ioutin)*(1.-exp(-binf(ibuild)*deltilf))&
                                                   *abinf(ibuild)/(binf(ibuild)*deltilf)+chiinSP(ibuild,ioutin)
                                 endif
                                 chiexSP(ibuild,ioutin)=0
                                 dosinSP(ibuild,ioutin)=deltilf*chiinSP(ibuild,ioutin)
                              endif
                           else
!                              chiin(ibuild,ioutin)=chiosum(ibuild,ioutin)+(chiin(ibuild,ioutin-1)-chiosum(ibuild,ioutin-1)) &
!                                                   *exp(-deltilf/tau(ibuild))-(chiosum(ibuild,ioutin)-chiosum(ibuild,ioutin-1)) &
!                                                   *(1.-exp(-deltilf/tau(ibuild)))*tau(ibuild)/deltilf
                              if(chiin(ibuild,ioutin) .eq. 0.)then
                                 chiin(ibuild,ioutin)=abinf(ibuild)*chiosum(ibuild,ioutin)+(chiin(ibuild,ioutin-1)&
                                                   -abinf(ibuild)*chiosum(ibuild,ioutin-1))*exp(-binf(ibuild)*deltilf)&
                                                   -(chiosum(ibuild,ioutin)-chiosum(ibuild,ioutin-1))&
                                                   *(1.-exp(-binf(ibuild)*deltilf))*abinf(ibuild)/(binf(ibuild)*deltilf)
                              endif
                              chiex(ibuild,ioutin)=chiin(ibuild,ioutin-1)*exp(-deltilf*binfex(ibuild))
                              dosin(ibuild,ioutin)=chiin(ibuild,ioutin)*deltilf+dosin(ibuild,ioutin-1)
                              if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
                                 if(chiinSP(ibuild,ioutin) .eq. 0.)then
                                    chiinSP(ibuild,ioutin)=abinf(ibuild)*chiosum(ibuild,ioutin)+(chiin(ibuild,ioutin-1)&
                                                   -abinf(ibuild)*chiosum(ibuild,ioutin-1))*exp(-binf(ibuild)*deltilf)&
                                                   -(chiosum(ibuild,ioutin)-chiosum(ibuild,ioutin-1))&
                                                   *(1.-exp(-binf(ibuild)*deltilf))*abinf(ibuild)/(binf(ibuild)*deltilf)
                                 endif
                                 chiexSP(ibuild,ioutin)=chiinSP(ibuild,ioutin-1)*exp(-deltilf*binfex(ibuild))
                                 dosinSP(ibuild,ioutin)=chiinSP(ibuild,ioutin)*deltilf+dosinSP(ibuild,ioutin-1)
                              endif
                           endif

!                           write(28,8111)ibuild,chiin(ibuild,ioutin),dosin(ibuild,ioutin)
!        8111               format(i5,2e12.4)
                           if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
                              write(28,8111)ibuild,chiin(ibuild,ioutin),dosin(ibuild,ioutin),chiexSP(ibuild,ioutin),&
                                    chiinSP(ibuild,ioutin),dosinSP(ibuild,ioutin),chiexSP(ibuild,ioutin)
                           else
                              write(28,8111)ibuild,chiin(ibuild,ioutin),dosin(ibuild,ioutin),chiex(ibuild,ioutin),&
                                    0.,0.,0.
                           endif
        8111               format(i5,6e12.4)
                           do iface=iface1n(ibuild),iface2n(ibuild)
                              do jface=jface1n(ibuild),jface2n(ibuild)
                                 ixfo=int(xfo(ibuild)/dx)+1
                                 jyfo=int(yfo(ibuild)/dy)+1
                                 xiurjur=xfo(ibuild)+(real(iface-ixfo)*dx+.5*dx)*cos(gamma(ibuild)*pi/180.)&
                                  -(real(jface-jyfo)*dy+.5*dy)*sin(gamma(ibuild)*pi/180)
                                 ifacr=int(xiurjur/dx)+1
                                 yiurjur=yfo(ibuild)+(real(jface-jyfo)*dy+.5*dy)*cos(gamma(ibuild)*pi/180)&
                                  +(real(iface-ixfo)*dx+.5*dx)*sin(gamma(ibuild)*pi/180)
                                  jfacr=int(yiurjur/dy)+1
                                 do kface=kface1n(ibuild),kface2n(ibuild)
                                    if(icellflag(ifacr,jfacr,kface).eq.0)then
                                       xf=dx*real(ifacr-1)+.5*dx
                                       yf=dy*real(jfacr-1)+.5*dy
                                       zf=dz*real(kface-2)+.5*dz
                                       ib=int((xf-xbl)/dbx)+1
                                       if(xf.lt.xbl)ib=0
                                       jb=int((yf-ybl)/dby)+1
                                       if(yf.lt.ybl)jb=0
                                       kb=int((zf-zbl)/dbz)+1
                                       if(zf.lt.zbl)kb=0
                                       if(ib.ge.1.and.ib.le.nbx.and.jb.ge.1.and.jb.le.nby.and.&
                                             kb.ge.1.and.kb.le.nbz.and.t.ge.concstrt)then
                                          con(ib,jb,kb)=chiin(ibuild,ioutin)
                                          dose(ib,jb,kb)=dosin(ibuild,ioutin)
                                          if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
                                             conSP(ib,jb,kb)=chiinSP(ibuild,ioutin)
                                             doseSP(ib,jb,kb)=dosinSP(ibuild,ioutin)
                                          endif
                                       endif
                                    endif
                                 enddo
                              enddo
                           enddo
                        enddo   lp012
                     endif
                  endif
                  if(toutc+0.5*delta .ge. outc .and. isiteflag .eq. 0)then
                     if(format_flag.eq.1.or.format_flag.eq.3)then
                        write(11,140)t
                        write(11,1150)nbx,nby,nbz
                        if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8) write(41,140)t
                        if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8) write(41,1150)nbx,nby,nbz
                        write(15,140)t
                        write(15,1150)nbx,nby,nbz
                        write(15,1160)
                        if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
                           write(45,140)t
                           write(45,1150)nbx,nby,nbz
                           write(45,1160)
                        endif
                        if(ilctflag.eq.1) then
                           write(14,140)t
                           write(14,1150)nbx,nby,nbz
                        endif
                        if(idepositionflag.eq.1) then
                           write(19,140)t
                           write(19,1150)nx1,ny1,nz1
                           write(19,1160)
                        endif
                        write(11,1160)
                        if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8) write(41,1160)
                        if(ilctflag.eq.1) write(14,1160)
                     endif
    1150             format('I=',i6,', J=',i6,',  K=',i6,', F=POINT')
    1160             format('DT=(SINGLE SINGLE SINGLE SINGLE)')
! erp 3/2/05 code to write out binary confield,confield10,dose,doseSP
                     if(format_flag.eq.2 .or. format_flag.eq.3)then !erp 3/2/05
                        if(printflag.eq.1)then
                           write(46)(xbc(ib),ib=1,nbx),(ybc(jb),jb=1,nby),(zbc(kb),kb=1,nbz)
                           printflag=0
                        endif
                        write(46)(((con(ib,jb,kb),ib=1,nbx),jb=1,nby),kb=1,nbz)
                        if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8) &
                              write(47)(((conSP(ib,jb,kb),ib=1,nbx),jb=1,nby),kb=1,nbz)
!                        if(idoseflag.eq.1)then
                        write(55)(((dose(ib,jb,kb),ib=1,nbx),jb=1,nby),kb=1,nbz)
                        if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8) &
                              write(56)(((doseSP(ib,jb,kb),ib=1,nbx),jb=1,nby),kb=1,nbz)
!                        endif
                     endif !end binary write
!end erp
                     do kb=1,nbz
                        do jb=1,nby
                           do ib=1,nbx
                              if((format_flag.eq.1.or.format_flag.eq.3 ).and.&
                                    isiteflag .eq. 0)then
                                 write(11,1169)xbc(ib),ybc(jb),zbc(kb),con(ib,jb,kb)
                                 if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8) &
                                       write(41,1169)xbc(ib),ybc(jb),zbc(kb),conSP(ib,jb,kb)
! Dosage is in grams*sec/m^3
                                 write(15,1169)xbc(ib),ybc(jb),zbc(kb),dose(ib,jb,kb)
                                 if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8) &
                                       write(45,1169)xbc(ib),ybc(jb),zbc(kb),doseSP(ib,jb,kb)
                              endif
                              if(ilctflag.eq.1)then
! note that concentrations are in grams per cubic meter while contours
! are in kg-sec per cubic meter
! MAN 7/24/2006 modified LCT's to account for variable number of thresholds
                                 conout=0.
                                 do ithreshold=1,nthreshold-1
                                    if(dose(ib,jb,kb).ge.1000.*dosethreshold(ithreshold).and. &
                                          dose(ib,jb,kb).lt.1000.*dosethreshold(ithreshold+1))then
                                       conout=real(ithreshold)
                                    endif
                                 enddo
                                 if(dose(ib,jb,kb).ge.1000.*dosethreshold(nthreshold))then
                                    conout=real(nthreshold)
                                 endif
! end MAN 7/24/2006
                                 write(14,1169)xbc(ib),ybc(jb),zbc(kb),conout
                              endif
!mdw corrected lct's 8-09-2004
                              con(ib,jb,kb)=0.
                              if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)conSP(ib,jb,kb)=0.
                           enddo
                        enddo
                     enddo
                     if(idepositionflag.eq.1) then
                        nsurfdep=0
                        if(format_flag.eq.1.or.format_flag.eq.3)then
                           write(19,140)t
                           if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8) write(49,140)t
                        endif
                        kface=1
                        xface=6
                        xblg=0
                        iprint=nx-1
                        jprint=ny-1
                        kprint=1
                        if(format_flag.eq.1.or.format_flag.eq.3)then
                           write(19,15005)iprint,jprint,kprint
                           write(19,16007)
                           if(ieradflag.eq.1) write(49,15005)iprint,jprint,kprint
                           if(ieradflag.eq.1) write(49,16007)
                        endif
                        zprint=0.
                        do jface=1,ny-1
                           do iface=1,nx-1
                              xf=dx*real(iface-1)+.5*dx
                              yf=dy*real(jface-1)+.5*dy
                              zf=0.
                              nsurfdep=nsurfdep+1
                              if(format_flag.eq.1.or.format_flag.eq.3)then
                                 write(19,16004)xf,yf,zf,depcz(iface,jface,kface+1),xblg,xface
                                 if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8) &
                                       write(49,16004)xf,yf,zf,depczSP(iface,jface,kface+1),xblg,xface
                              endif
                              if(format_flag.eq.2.or.format_flag.eq.3)then
                                 write(58)xf,yf,zf,depcz(iface,jface,kface+1)
                                 if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8) &
                                       write(59)xf,yf,zf,depczSP(iface,jface,kface+1)
                              endif
                           enddo
                        enddo
  16004                 format(3f10.3,e12.4,2f5.0)
                        jprint=1
                        kprint=1
lp013:                  do ibuild=1,inumbuild
                           if(bldtype(ibuild) .eq. 9)cycle
                           if(aminusx(ibuild).ge.1)then
                              nface=1
                              xblg=real(ibuild)
                              xface=real(nface)
                              iprint=int(aminusx(ibuild)+.5)
                              if(format_flag.eq.1.or.format_flag.eq.3)then
                                 write(19,140)t
                                 if(ieradflag.eq.1) write(49,140)t
                                 write(19,15005)iprint,jprint,kprint
                                 if(ieradflag.eq.1) write(49,15005)iprint,jprint,kprint
15005                            format('I=',i6,', J=',i6,', k=',i6)
16007                            format('DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)')
                                 write(19,16007)
                                 if(ieradflag.eq.1) write(49,16007)
                              endif
                              do iface=iface1n(ibuild),iface2n(ibuild)
                                 do kface=kface1n(ibuild),kface2n(ibuild)
                                    do jface=jface1n(ibuild),jface2n(ibuild)
                                       if(icellflag(iface-1,jface,kface) .ne. 0 .and. icellflag(iface,jface,kface) .eq. 0)then
                                          xf=dx*real(iface-1)
                                          yf=dy*real(jface-1)+.5*dy
                                          zf=dz*real(kface-2)+.5*dz
                                          nsurfdep=nsurfdep+1
                                          if(format_flag.eq.1.or.format_flag.eq.3)then
                                             write(19,16004)xf,yf,zf,depcx(iface-1,jface,kface),xblg,xface
                                             if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8) &
                                                   write(49,16004)xf,yf,zf,depcxSP(iface-1,jface,kface),xblg,xface
                                          endif
                                          if(format_flag.eq.2.or.format_flag.eq.3)then
                                             write(58)xf,yf,zf,depcx(iface-1,jface,kface)
                                             if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8) &
                                                   write(59)xf,yf,zf,depcxSP(iface-1,jface,kface)
                                          endif
                                       endif
                                    enddo
                                 enddo
                              enddo
                           endif
                           if(aplusx(ibuild).ge.1)then
                              nface=2
                              xblg=real(ibuild)
                              xface=real(nface)
                              iprint=int(aplusx(ibuild)+.5)
                              if(format_flag.eq.1.or.format_flag.eq.3)then
                                 write(19,140)t
                                 write(19,15005)iprint,jprint,kprint
                                 write(19,16007)
                                 if(ieradflag.eq.1) then
                                    write(49,140)t
                                    write(49,15005)iprint,jprint,kprint
                                    write(49,16007)
                                 endif
                              endif
                              do iface=iface1n(ibuild),iface2n(ibuild)
                                 do kface=kface1n(ibuild),kface2n(ibuild)
                                    do jface=jface1n(ibuild),jface2n(ibuild)
                                       if(icellflag(iface+1,jface,kface) .ne. 0 .and. icellflag(iface,jface,kface) .eq. 0)then
                                          xf=dx*real(iface-1)+dx
                                          yf=dy*real(jface-1)+.5*dy
                                          zf=dz*real(kface-2)+.5*dz
                                          nsurfdep=nsurfdep+1
                                          if(format_flag.eq.1.or.format_flag.eq.3)then
                                             write(19,16004)xf,yf,zf,depcx(iface+1,jface,kface),xblg,xface
                                             if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8) &
                                                   write(49,16004)xf,yf,zf,depcxSP(iface+1,jface,kface),xblg,xface
                                          endif
                                          if(format_flag.eq.2.or.format_flag.eq.3)then
                                             write(58)xf,yf,zf,depcx(iface+1,jface,kface)
                                             if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8) &
                                                   write(59)xf,yf,zf,depcxSP(iface+1,jface,kface)
                                          endif
                                       endif
                                    enddo
                                 enddo
                              enddo
                           endif
                           if(aminusy(ibuild).ge.1)then
                              nface=3
                              xblg=real(ibuild)
                              xface=real(nface)
                              iprint=int(aminusy(ibuild)+.5)
                              if(format_flag.eq.1.or.format_flag.eq.3)then
                                 write(19,140)t
                                 write(19,15005)iprint,jprint,kprint
                                 write(19,16007)
                                 if(ieradflag.eq.1) then
                                    write(49,140)t
                                    write(49,15005)iprint,jprint,kprint
                                    write(49,16007)
                                 endif
                              endif
                              do jface=jface1n(ibuild),jface2n(ibuild)
                                 do kface=kface1n(ibuild),kface2n(ibuild)
                                    do iface=iface1n(ibuild),iface2n(ibuild)
                                       if(icellflag(iface,jface-1,kface) .ne. 0 .and. icellflag(iface,jface,kface) .eq. 0)then
                                          xf=dx*real(iface-1)+.5*dx
                                          yf=dy*real(jface-1)
                                          zf=dz*real(kface-2)+.5*dz
                                          nsurfdep=nsurfdep+1
                                          if(format_flag.eq.1.or.format_flag.eq.3)then
                                             write(19,16004)xf,yf,zf,depcy(iface,jface-1,kface),xblg,xface
                                             if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8) &
                                                   write(49,16004)xf,yf,zf,depcySP(iface,jface-1,kface),xblg,xface
                                          endif
                                          if(format_flag.eq.2.or.format_flag.eq.3)then
                                             write(58)xf,yf,zf,depcy(iface,jface-1,kface)
                                             if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8) &
                                                   write(59)xf,yf,zf,depcySP(iface,jface-1,kface)
                                          endif
                                       endif
                                    enddo
                                 enddo
                              enddo
                           endif
                           if(aplusy(ibuild).ge.1)then
                              nface=4
                              xblg=real(ibuild)
                              xface=real(nface)
                              iprint=int(aplusy(ibuild)+.5)
                              if(format_flag.eq.1.or.format_flag.eq.3)then
                                 write(19,140)t
                                 write(19,16007)
                                 if(ieradflag.eq.1) then
                                    write(49,140)t
                                    write(49,15005)iprint,jprint,kprint
                                    write(49,16007)
                                 endif
                              endif
                              do jface=jface1n(ibuild),jface2n(ibuild)
                                 do kface=kface1n(ibuild),kface2n(ibuild)
                                    do iface=iface1n(ibuild),iface2n(ibuild)
                                       if(icellflag(iface,jface+1,kface) .ne. 0 .and. icellflag(iface,jface,kface) .eq. 0)then
                                          xf=dx*real(iface-1)+.5*dx
                                          yf=dy*real(jface-1)+dy
                                          zf=dz*real(kface-2)+.5*dz
                                          nsurfdep=nsurfdep+1
                                          if(format_flag.eq.1.or.format_flag.eq.3)then
                                             write(19,16004)xf,yf,zf,depcy(iface,jface+1,kface),xblg,xface
                                             if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8) &
                                                   write(49,16004)xf,yf,zf,depcySP(iface,jface+1,kface),xblg,xface
                                          endif
                                          if(format_flag.eq.2.or.format_flag.eq.3)then
                                             write(58)xf,yf,zf,depcy(iface,jface+1,kface)
                                             if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8) &
                                                   write(59)xf,yf,zf,depcySP(iface,jface+1,kface)
                                          endif
                                       endif
                                    enddo
                                 enddo
                              enddo
                           endif
                           if(afloor(ibuild).ge.1)then
                              nface=5
                              xblg=real(ibuild)
                              xface=real(nface)
                              iprint=int(afloor(ibuild)+.5)
                              if(format_flag.eq.1.or.format_flag.eq.3)then
                                 write(19,140)t
                                 write(19,15005)iprint,jprint,kprint
                                 write(19,16007)
                                 if(ieradflag.eq.1) then
                                    write(49,140)t
                                    write(49,15005)iprint,jprint,kprint
                                    write(49,16007)
                                 endif
                              endif
                              do kface=kface1n(ibuild),kface2n(ibuild)
                                 do jface=jface1n(ibuild),jface2n(ibuild)
                                    do iface=iface1n(ibuild),iface2n(ibuild)
                                       if(icellflag(iface,jface,kface-1) .ne. 0 .and. icellflag(iface,jface,kface) .eq. 0)then
                                          xf=dx*real(iface-1)+.5*dx
                                          yf=dy*real(jface-1)+.5*dy
                                          zf=dz*real(kface-2)
                                          nsurfdep=nsurfdep+1
                                          if(format_flag.eq.1.or.format_flag.eq.3)then
                                             write(19,16004)xf,yf,zf,depcz(iface,jface,kface-1),xblg,xface
                                             if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8) &
                                                   write(49,16004)xf,yf,zf,depczSP(iface,jface,kface-1),xblg,xface
                                          endif
                                          if(format_flag.eq.2.or.format_flag.eq.3)then
                                             write(58)xf,yf,zf,depcz(iface,jface,kface-1)
                                             if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8) &
                                                   write(59)xf,yf,zf,depczSP(iface,jface,kface-1)
                                          endif
                                       endif
                                    enddo
                                 enddo
                              enddo
                           endif
                           if(aroof(ibuild).ge.1)then
                              nface=6
                              xblg=real(ibuild)
                              xface=real(nface)
                              iprint=int(aroof(ibuild)+.5)
                              if(format_flag.eq.1.or.format_flag.eq.3)then
                                 write(19,140)t
                                 write(19,15005)iprint,jprint,kprint
                                 write(19,16007)
                                 if(ieradflag.eq.1) then
                                    write(49,140)t
                                    write(49,15005)iprint,jprint,kprint
                                    write(49,16007)
                                 endif
                              endif
                              do kface=kface1n(ibuild),kface2n(ibuild)
                                 do jface=jface1n(ibuild),jface2n(ibuild)
                                    do iface=iface1n(ibuild),iface2n(ibuild)
                                       if(icellflag(iface,jface,kface+1) .ne. 0 .and. icellflag(iface,jface,kface) .eq. 0)then
                                          xf=dx*real(iface-1)+.5*dx
                                          yf=dy*real(jface-1)+.5*dy
                                          zf=dz*real(kface-2)+dz
                                          nsurfdep=nsurfdep+1
                                          if(format_flag.eq.1.or.format_flag.eq.3)then
                                             write(19,16004)xf,yf,zf,depcz(iface,jface,kface+1),xblg,xface
                                             if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8) &
                                                   write(49,16004)xf,yf,zf,depczSP(iface,jface,kface+1),xblg,xface
                                          endif
                                          if(format_flag.eq.2.or.format_flag.eq.3)then
                                             write(58)xf,yf,zf,depcz(iface,jface,kface+1)
                                             if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8) &
                                                   write(59)xf,yf,zf,depczSP(iface,jface,kface+1)
                                          endif
                                       endif
                                    enddo
                                 enddo
                              enddo
                           endif
                        enddo   lp013
                        if(format_flag.eq.2.or.format_flag.eq.3)then
                           open(unit=999,file="QP_numsurfdepcell.bin",form="unformatted",status="unknown")
                           write(999)nsurfdep
                           close(999)
                        endif
                     endif
                     toutc=toutc-outc
                  endif
                  if(toutc + 0.5*delta .ge. outc .and. isiteflag .eq. 1)then
                     gridind=real((nxg)*(iyg-1)+ixg)
                     if(format_flag.eq.1.or.format_flag.eq.3)then
                        write(21,140)gridind
                        write(21,1150)nbx,nby,nbz
                        write(21,1160)
                        kb=1
                        do jb=1,nby
                           do ib=1,nbx
                              write(21,1169)xbc(ib),ybc(jb),zbc(kb),dosexy(ib,jb)
! Dosage is in grams*sec/m^3
! note that concentrations are in grams per cubic meter while contours
! are in kg-sec per cubic meter
                           enddo
                        enddo
                     endif
                     if(format_flag.eq.2.or.format_flag.eq.3)then
                        write(22)((dosexy(ib,jb),ib=1,nbx),jb=1,nby)
                     endif
                     dosexy(:,:)=0.
                     toutc=0. !toutc-outc
                  endif
                  mp=nptot
                  if(t .ge. dur .and. ibuoyflag .eq. 1 .and. inextgridflag .lt. 2)then
                     wbulkzbar=wbulkzs/sumwbulkz  
                     write(*,21251)wbulkzbar
                  endif
! MAN 10/18/2007 fixed floating point error in file output times
                  if(delta_ref_rm .eq. delta_ref)then
                     toutc=real(nint(toutc/delta_ref))*delta_ref
                     tout=real(nint(tout/delta_ref))*delta_ref
                     touti=real(nint(touti/delta_ref))*delta_ref
                     t=real(nint(t/delta_ref))*delta_ref
                  endif
! end MAN 10/18/2007
21251          format('remaining buoyant velocity=',f8.3)
!                  it=it+1
               enddo   lp014
!mdw 3/31/2004 changed format to allow larger values of x, y, and z
1169           format(3(1pe11.4), e12.4)
!              write(*,*)xpmin,xpmax,ypmin,ypmax,zpmin,zpmax
!TMB
!              write(*,2125)xpmin,xpmax,ypmin,ypmax,zpmin,zpmax
 2125          format('xpmin=',f6.2,' xpmax=',f6.2,' ypmin=',f6.2,' ypmax=',f6.2,' zpmin=',f6.2,' zpmax=',f6.2)
            enddo
         enddo
         !call kill(windgrid)
         nxw=2*(nx-1)
         nyw=2*(ny-1)
         nzw=2*(nz-1)
         nzwout=nzw-2
         if(iwallfieldflag.eq.1) write(12,2120)
 2120    format('TITLE     = "Buildings in 3D"')
         if(iwallfieldflag.eq.1) write(12,2130)
 2130    format('VARIABLES = "X" "Y" "Z" "B"')
         if(iwallfieldflag.eq.1) write(12,2135)
 2135    format('ZONE T = "BUILDINGS" ')
         if(iwallfieldflag.eq.1) write(12,1150)nxw,nyw,nzwout
         if(iwallfieldflag.eq.1)then
            do kw=3,nzw
               zw=.5*dz*(kw-1)-1.*dz
               zwp=zw-.01*dz
               km=int((zw+dz)/dz)+1
               km=min(km,nz-1)
               if(kw.le.2)km=3
               kmp=int((zwp+dz)/dz)+1
               if(kmp.le.2)kmp=3
               kmp=min(kmp,nz-1)
               do jw=1,nyw
                  yw=.5*dy*real(jw-1)
                  ywp=yw-.01
                  jm=int(yw/dy)+1
                  jmp=int(ywp/dy)+1
                  jmp=min(jmp,ny-1)
                  do iw=1,nxw
                     xw=.5*dx*real(iw-1)
                     xwp=xw-.01
                     im=int(xw/dx)+1
                     imp=int(xwp/dx)+1
                     imp=min(imp,nx-1)
                     icellflg=min(icellflag(im,jm,km),icellflag(imp,jm,km),icellflag(im,jmp,km), &
                              icellflag(im,jm,kmp),icellflag(imp,jmp,km),icellflag(imp,jmp,kmp))
                     if(icellflg.eq.0)then
                        b=0.
                     else
                        b=1.
                     endif
                     write(12,*)xw,yw,zw,b
                  enddo
               enddo
            enddo
         endif
!man 7/12/2005 conditional closing statements
         print*,'Minimum Lagrangian time scale =',tls_min,'s'
!         print*,'Maximum Taylor Microscale =',taylor_max,'s'
         close(5)
         close(7)
         close(8)
!         close(9)
         close(30)
         close(34)
         close(35)
         close(37)
         close(44)
         close(62)
         close(77)
!         close(88)
         open(unit=999,file="QP_nptotfinal.bin",form="unformatted",status="unknown")
         write(999)nptot
         close(999)
         if(format_flag.eq.1)then
            close(27)
         else
            close(50)
         endif
         if(iwallfieldflag.eq.1) close(12)
         if(isiteflag.eq.1)then ! only used for sensor siting.
            close(16)
            if(format_flag.eq.1 .or. format_flag.eq.3)then
               close(21)
            endif
            if(format_flag.eq.2 .or. format_flag.eq.3)then
               close(22)
            endif
         else
            if(ilctflag.eq.1) close(14)
            if(iturbtypeflag.eq.1) close(101)
            if(format_flag.eq.1 .or. format_flag.eq.3)then
               close(10)
               close(11)
               close(15)
               close(19)
               if(iturbfieldflag.eq.1)close(13)
               if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
                  close(41)
                  close(45)
                  if(idepositionflag.eq.1)close(49)
               endif
            endif
            if(format_flag.eq.2 .or. format_flag.eq.3)then
               close(46)
               close(54)
               close(55)
               close(58)
               if(iturbfieldflag.eq.1)close(13)
               if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
                  close(47)
                  close(56)
                  if(idepositionflag.eq.1)close(59)
               endif
            endif
         endif
         if(isourcetypeflag .eq. 8)close(25)
         if(inextgridflag .ge. 1) then
! MAN 9/13/2005 binary nextgrid
            if(format_flag.eq.1 .or. format_flag.eq.3)close(31)
            if(format_flag.eq.2 .or. format_flag.eq.3)close(32)
            if(format_flag.eq.1)then
               close(33)
            else
               close(38)
            endif
         endif
! end MAN 9/13/2005
! MAN 7/24/2006
         if(ieradflag.eq.1)then
            close(69)
         else
            close(29)
         endif
         
!end MAN 7/24/2006
!end man 7/12/2005

! mdw 4-22-2004 close the file with out of grid particle positions
         deallocate(icellflag)
!         deallocate(iyflag)
         deallocate(icb,jcb)
         deallocate(dxm,dym,dzm)
         deallocate(dxp,dyp,dzp)
         if(idepositionflag.eq.1)then
            deallocate(depcx,depcy,depcz)
            if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)deallocate(depcxSP,depcySP,depczSP)
         endif
         deallocate(u,v,w,hgt)
         deallocate(dutotdzi,ustarz,sigwi)
         deallocate(dutotdyi,dutotdxi,dutotdni,dutotdsi)
         deallocate(alph1ij,alph2ij,alph3ij,bet1ij,bet2ij,bet3ij)
         deallocate(gam1ij,gam2ij,gam3ij,alphn1ij,alphn2ij,alphn3ij)
         deallocate(betn1ij,betn2ij,betn3ij,gamn1ij,gamn2ij,gamn3ij)
         deallocate(dsigwdzi)
         deallocate(sigvi,epsi)
         deallocate(xbc,ybc,zbc)
         if(isiteflag.eq.1)then
            deallocate(dosexy)
         else
            deallocate(con,dose)
            if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)deallocate(conSP,doseSP)
         endif
         deallocate(xi,yi,zi,dsigwdni)
         deallocate(xp,yp,zp,tstrt,weight,weightd)
         if(ibuoyflag.gt.0)deallocate(zpp,temp_part)
         deallocate(uran,vran,wran,pmass,dwallp)
         deallocate(xps,yps,zps)
         deallocate(source_num)
!         deallocate(xpr,ypr,zpr,trel,rpr)
         deallocate(elz,eleff)
         deallocate(Ht,Wti,Lt)
         if(iturbtypeflag.eq.1)deallocate(local_option,nonlocal_option,turb_options)
! mdw 4-20-2004 deallocate zfo
         deallocate(xfo,yfo,zfo,ustarij)
         deallocate(weff,leff,atten,gamma)
         deallocate(lr,sx,sy,lfr)
         deallocate(xcb,ycb,uref)
         deallocate(urefu,urefv,urefw)
! mdw 10-04-2004 added zcorf to deallocations
         deallocate(ustarg,deluc,ustargz,elzg,zcorf)
         deallocate(bldtype,bldnum)
         deallocate(istart,iend,jstart,jend)
!         deallocate(ufsqi,vfsqi,wfsqi)
         deallocate(upwpi,ani,bni,cni)
!         deallocate(asi,bsi,csi)
         deallocate(vsi,taudi,diffi,sci,dpm,weightr)
         deallocate(dosethreshold)
         deallocate(source_x,source_y,source_z)
         deallocate(source_l,source_w,source_h,source_weight)
         deallocate(source_r,source_hemass,source_strength)
         deallocate(source_density,source_start,source_dur)
         deallocate(source_geometry,source_str_units,source_release)
         deallocate(source_nodespeed,source_nodeaccel,source_nodesegpos)
         deallocate(source_seglength,source_startidx,source_stopidx,nlinenodes)
         deallocate(source_prr,source_nodestart,source_nodestop)
         if(ieradflag.eq.1)then
            deallocate(isub)
            deallocate(pufmasfr,znorm,psigx,psigy,psigz,pufcummas)
         endif
         deallocate(xnear,ynear,znear)
! mdw 1-22-2004 deallocations for building infiltration
         deallocate(utotmax,utotcl1)
         if(iindoorflag.eq.1)then
            close(28)
            deallocate(iface1n,iface2n,jface1n,jface2n,kface1n,kface2n)
            deallocate(aminusx,aplusx,aminusy,aplusy,aroof,afloor)
            deallocate(tau,racact,racref,tiref)
            deallocate(toref,piref,urefb,ainf,binf,abinf,binfex)
            deallocate(tiact,rhoref,racsupply)
            deallocate(racintake,filtration,area2volume)
            deallocate(chioa,chiosum,chiin,dosin,chiex)
            if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
               deallocate(chioaSP,chiosumSP,chiinSP,dosinSP,chiexSP)
            endif
         endif
         if(ibioslurryflag .eq. 1 .or. itwophaseflag .eq. 1)then
            deallocate(temperature,vaporfield)
            close(88)
         endif
! end new deallocations
!TMB
!use return and comment end for matlab compilable program
!        return
         stop
      end

!
! RANDOM NUMBER GENERATOR FUNCTION
! UNIFORM DISTRIBUTION IN INTERVAL BETWEEN O AND 1
! Method from Numerical Recipes, p197
      real function RAN2() result(ran)
         implicit none
         integer ia,ic,ir(97),iff,idum,iy,m,j
         real rm
         parameter(m=714025,ia=1366,ic=150889,rm=1./m)
!         dimension ir(97)
         save iff,idum,iy,ir
! change idum (any negative integer) to change sequence
         data iff,idum/0,-100/
         if(iff.eq.0) then
            iff=1
            idum=mod(ic-idum,m)
            do j=1,97
               idum=mod(ia*idum+ic,m)
               ir(j)=idum
            enddo
            idum=mod(ia*idum+ic,m)
            iy=idum
         else
         endif
         j=1+(97*iy)/m
         if(j.gt.97.or.j.lt.1) then
            write(*,'(a,i5)') 'problem in random generator, j=',j
            stop 97
         endif
         iy=ir(j)
         ran=iy*rm
         idum=mod(ia*idum+ic,m)
         ir(j)=idum
         return
      end function RAN2
!
! RANDOM NUMBER GENERATOR SUBROUTINE
! GAUSSIAN DISTRIBUTION W/ MEAN=0, VAR=1
! Method from Marsaglia and Bray, 1964 (found in Devroye, 1986, p390)
      subroutine rangen(x)
         implicit none
         real, intent(out) :: x
         real u,v,w,summ,ran2
         u=RAN2()
         if (u.le..8638) then
!           v1=RAN2( )
!           v2=RAN2( )
!           v=sign(v1,v2-0.5)
!           w1=RAN2( )
!           w2=RAN2( )
!           w=sign(w1,w2-0.5)
            v=2.0*RAN2()-1.0
            w=2.0*RAN2()-1.0
            x=2.315351*u-1.0+v+w
            return
         endif
         if (u.le..9745) then
            v=RAN2()
            x=1.5*(v-1.0+9.033424*(u-.8638))
            return
         endif
         if (u.gt..9973002) then
            v = 4.5
            X = 1.
            do while((x*v*v).gt.4.5)
               v=RAN2()
               w=RAN2()
               w=max(1.e-07,w)
               x=4.5-log(w)
            enddo
            x=sign(sqrt(2.0*x),(u-.9986501))
         else
            u = 50.
            v = 0.
            summ = 49.00244
            w = 0.
            do while(u.gt.(49.00244*exp(-v*v/2.0)-summ-w))
!              x1=RAN2( )
!              x2=RAN2( )
!              x=3.0*sign(x1,x2-0.5)
               x=6.0*RAN2()-3.0
               u=RAN2()
               v=abs(x)
               w=6.631334*(3.0-v)**2.0
               summ=0.
               if (v.lt.1.5) summ=6.043281*(1.5-v)
               if (v.lt.1.0)  summ=summ+13.26267*(3.0-v*v)-w
            enddo
         endif
         return
      end

      real function erf(x)
         implicit none
         REAL x,t,z
         z=abs(x)
         t=1./(1.+0.5*z)
         erf=1.-t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+t*&
             (.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+t*&
             (1.48851587+t*(-.82215223+t*.17087277)))))))))
         if (x.lt.0.) erf=-erf
         return
      end

      real function erfinv(x)
         implicit none
         real x,z,erf,u,pi
         real a(4),b(4),c(4),d(2)
         if(x .ge. .999999)x=.999999
         if(x .le. -.999999)x=-.999999
         a(1) = 0.886226899
         a(2) = -1.645349621
         a(3) = 0.914624893
         a(4) = -0.140543331
         b(1) = -2.118377725
         b(2) = 1.442710462
         b(3) = -0.329097515
         b(4) = 0.012229801
         c(1) = -1.970840454
         c(2) = -1.624906493
         c(3) = 3.429567803
         c(4) = 1.641345311
         d(1) = 3.543889200
         d(2) = 1.637067800
         pi=4.*atan(1.)
! Central range

         if(abs(x) .le. .7)then
            z = x*x;
            erfinv = x*(((a(4)*z+a(3))*z+a(2))*z+a(1))/((((b(4)*z+b(3))*z+b(2))*z+b(1))*z+1.)
         endif

! Near end points of range

         if(x .gt. .7 .and. x .le. 1.)then
            z = sqrt(-log((1.-x)/2.))
            erfinv = (((c(4)*z+c(3))*z+c(2))*z+c(1))/((d(2)*z+d(1))*z+1.)
         elseif(x .lt. -.7 .and. x .ge. -1.)then
            z = sqrt(-log((1.+x)/2.));
            erfinv = -(((c(4)*z+c(3))*z+c(2))*z+c(1))/((d(2)*z+d(1))*z+1.)
         endif

! The relative error of the approximation has absolute value less
! than 8.9e-7.  One iteration of Halley's rational method (third
! order) gives full machine precision.

! Newton's method: new x = x - f/f'
! Halley's method: new x = x - 1/(f'/f - (f"/f')/2)
! This function: f = erf(x) - y, f' = 2/sqrt(pi)*exp(-x^2), f" = -2*x*f'

! Newton's correction
         u = (erf(erfinv)-x)/((2./sqrt(pi))*exp(-erfinv*erfinv))

! Halley's step
         erfinv = erfinv - u/(1.+erfinv*u)
         return
      end
!     integer function myhandler(sig,code,context)
!        integer sig,code,context(5)
!        call abort()
!     end

      real function buildingradius(a,b,gamma,theta)
         implicit none
         real a,b,gamma,theta
         buildingradius=a*b/sqrt( (a*sin(theta-gamma))**2. + (b*cos(theta-gamma))**2. )
         return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine windterp
         use variables   !make data from module "variables" visible
         implicit none
         integer i1,i2,j1,j2,k1,k2
!        real up,vp,wp,dsum,usum,vsum,wsum
         if(dxc.eq.0.and.dyc.eq.0.and.dzc.eq.0)then
            uint=u(im,jm,km)
            vint=v(im,jm,km)
            wint=w(im,jm,km)
            return
         endif
         if(dwall.lt.dx)then
            uint=u(im,jm,km)*log(dwall/z0)/log(dx/(2.*z0))
            vint=v(im,jm,km)*log(dwall/z0)/log(dx/(2.*z0))
            wint=w(im,jm,km)*log(dwall/z0)/log(dx/(2.*z0))
            return
         endif
!MDW 7-08-2005 added log-law interpolation for y & z directions
         if(dwall.lt.dy)then
            uint=u(im,jm,km)*log(dwall/z0)/log(dy/(2.*z0))
            vint=v(im,jm,km)*log(dwall/z0)/log(dy/(2.*z0))
            wint=w(im,jm,km)*log(dwall/z0)/log(dy/(2.*z0))
            return
         endif
         if(dwall.lt.dz)then
            uint=u(im,jm,km)*log(dwall/z0)/log(dz/(2.*z0))
            vint=v(im,jm,km)*log(dwall/z0)/log(dz/(2.*z0))
            wint=w(im,jm,km)*log(dwall/z0)/log(dz/(2.*z0))
            return
         endif
         dsum=0.
         usum=0
         vsum=0.
         wsum=0.
         i1=im-1
         i2=im+1
         i1=max(1,i1)
         i2=min(nx-1,i2)
         j1=jm-1
         j2=jm+1
         j1=max(1,j1)
         j2=min(ny-1,j2)
         k1=km-1
         k1=max(2,k1)
         k2=km+1
         k2=min(nz-1,k2)
lp017:   do k=k1,k2
lp016:      do j=j1,j2
lp015:         do i=i1,i2
                  up=u(i,j,k)
                  vp=v(i,j,k)
                  wp=w(i,j,k)
                  dxx=dx*real(i-im)-dxc
                  dyy=dy*real(j-jm)-dyc
                  dzz=dz*real(k-km)-dzc
                  if(i.lt.im.and.dxm(im,jm,km).lt.abs(dxx))then
                     dxx=dxm(im,jm,km)+dxc
                     up=0.
                     vp=0.
                     wp=0.
                     if(j.eq.jm.and.k.eq.km)then
! Here we consider distance perpendicular to the wall
                        dyy=0.
                        dzz=0.
                     endif
                  endif
                  if(i.gt.im.and.dxp(im,jm,km).lt.dxx)then
                     dxx=dxp(im,jm,km)-dxc
                     if(j.eq.jm.and.k.eq.km)then
! Here we consider distance perpendicular to the wall
                        dyy=0.
                        dzz=0.
                     endif
                     up=0.
                     vp=0.
                     wp=0.
                  endif
                  if(j.lt.jm.and.dym(im,jm,km).lt.abs(dyy))then
                     dyy=dym(im,jm,km)+dyc
                     if(i.eq.im.and.k.eq.km)then
! Here we consider distance perpendicular to the wall
                        dxx=0.
                        dzz=0.
                     endif
                     up=0.
                     vp=0.
                     wp=0.
                  endif
                  if(j.gt.jm.and.dyp(im,jm,km).lt.dyy)then
                     dyy=dyp(im,jm,km)-dyc
                     if(i.eq.im.and.k.eq.km)then
! Here we consider distance perpendicular to the wall
                        dxx=0.
                        dzz=0.
                     endif
                     up=0.
                     vp=0.
                     wp=0.
                  endif
                  if(k.lt.km.and.dzm(im,jm,km).lt.abs(dzz))then
                     dzz=dzm(im,jm,km)+dzc
                     if(i.eq.im.and.j.eq.jm)then
! Here we consider distance perpendicular to the wall
                        dyy=0.
                        dxx=0.
                     endif
                     up=0.
                     vp=0.
                     wp=0.
                  endif
                  dsum=dsum+1./sqrt(dxx**2+dyy**2+dzz**2)
                  usum=usum+up/sqrt(dxx**2+dyy**2+dzz**2)
                  vsum=vsum+vp/sqrt(dxx**2+dyy**2+dzz**2)
                  wsum=wsum+wp/sqrt(dxx**2+dyy**2+dzz**2)
                  if(km.eq.2)then
                     dzz=dz/2.+dzc
                     if(i.eq.im.and.j.eq.jm)then
                        dxx=0.
                        dyy=0.
                     endif
                     up=0.
                     vp=0.
                     wp=0.
                     dsum=dsum+1./sqrt(dxx**2+dyy**2+dzz**2)
                     usum=usum+up/sqrt(dxx**2+dyy**2+dzz**2)
                     vsum=vsum+vp/sqrt(dxx**2+dyy**2+dzz**2)
                     wsum=wsum+wp/sqrt(dxx**2+dyy**2+dzz**2)
                  endif
               enddo   lp015
            enddo   lp016
         enddo   lp017
         uint=usum/dsum
         vint=vsum/dsum
         wint=wsum/dsum
         return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine reflect
         use variables
         implicit none
! mdw 8-11-2006 new subroutine to do reflection from z0 removed old reflection
! coding
         real dwmax,dxxz0,dyyz0,dzzz0,dxxwr,dyywr,dzzwr,dxyzm,distance_pen
         integer isignr,jsignr,ksignr,immr,jmmr,kmmr,iz0,jz0,kz0
         integer icheck,ikcheck,iicheck,ijcheck,i2check
         dwmax=max(abs(dxx),abs(dyy),abs(dzz)) !largest position change
         dxxz0=dxx*z0/dwmax
         dyyz0=dyy*z0/dwmax
         dzzz0=dzz*z0/dwmax
         isignr=0
         jsignr=0
         ksignr=0
         immr=0
         jmmr=0
         kmmr=0
         dxxwr=0.
         dyywr=0.
         dzzwr=0.
         dxyzm=0.
         do iz0=1,3,2
            do jz0=1,3,2
               do kz0=1,3,2
                  imm=int((xp(m)+z0coeff*real(iz0-2)*z0)/dx)+1
                  jmm=int((yp(m)+z0coeff*real(jz0-2)*z0)/dy)+1
                  kmm=int((zp(m)+z0coeff*real(kz0-2)*z0+dz)/dz)+1
                  dxxw=0.
                  dyyw=0.
                  dzzw=0.
                  if((imm.ge.1.and.imm.le.nx-1).and.(jmm.ge.1.and.jmm.le.ny-1).and.(kmm.ge.1.and.kmm.le.nz-1))then
                     if(icellflag(imm,jmm,kmm) .eq. 0)then
! reflection required - note this includes the extension of the walls outward by z0
                        if(icellflag(imm,jm,km) .eq. 0 .or. icellflag(im,jmm,kmm) .gt. 0 .or.&
                              (icellflag(imm,jmm,km) .eq. 0 .and. icellflag(im,jmm,km) .gt. 0).or.&
                              (icellflag(imm,jm,kmm) .eq. 0 .and. icellflag(im,jm,kmm) .gt. 0))then
                           isign=(imm-im)/abs(imm-im)
!isign>0 means approach to wall from left; <0 means approach from right
                           dxxw=real(isign)*(-dx*real(max(imm,im)-1)+z0coeff*real(isign)*z0+xp(m))
                           dxxw=z0coeff*dxxw
                        endif
                        if(icellflag(im,jmm,km) .eq. 0 .or. icellflag(imm,jm,kmm) .gt. 0 .or.&
                              (icellflag(imm,jmm,km) .eq. 0 .and. icellflag(imm,jm,km) .gt. 0).or.&
                              (icellflag(im,jmm,kmm) .eq. 0 .and. icellflag(im,jm,kmm) .gt. 0))then
                           jsign=(jmm-jm)/abs(jmm-jm)
!isign>0 means approach to wall from left; <0 means approach from right
                           dyyw=real(jsign)*(-dy*real(max(jmm,jm)-1)+z0coeff*real(jsign)*z0+yp(m))
                           dyyw=z0coeff*dyyw
                        endif
                        if(icellflag(im,jm,kmm) .eq. 0 .or. icellflag(imm,jmm,km) .gt. 0 .or.&
                              (icellflag(im,jmm,kmm) .eq. 0 .and. icellflag(im,jmm,km) .gt. 0).or.&
                              (icellflag(imm,jm,kmm) .eq. 0 .and. icellflag(imm,jm,km) .gt. 0))then
                           ksign=(kmm-km)/abs(kmm-km)
!isign>0 means approach to wall from left; <0 means approach from right
                           dzzw=real(ksign)*(-dz*real(max(kmm,km)-2)+z0coeff*real(ksign)*z0+zp(m))
                           dzzw=z0coeff*dzzw
                        endif
                        distance_pen=sqrt(dxxw**2+dyyw**2+dzzw**2)
                        if(distance_pen.eq.0.)then
                           if(imm .eq. im)then
                              isign=0
                           else
                              isign=(imm-im)/abs(imm-im)
                           endif
!isign>0 means approach to wall from left; <0 means approach from right
                           dxxw=real(isign)*(-dx*real(max(imm,im)-1)+z0coeff*real(isign)*z0+xp(m))
                           dxxw=z0coeff*dxxw
                           if(jmm .eq. jm)then
                              jsign=0
                           else
                              jsign=(jmm-jm)/abs(jmm-jm)
                           endif
!isign>0 means approach to wall from left; <0 means approach from right
                           dyyw=real(jsign)*(-dy*real(max(jmm,jm)-1)+z0coeff*real(jsign)*z0+yp(m))
                           dyyw=z0coeff*dyyw
                           if(kmm .eq. km)then
                              ksign=0
                           else
                              ksign=(kmm-km)/abs(kmm-km)
                           endif
!isign>0 means approach to wall from left; <0 means approach from right
                           dzzw=real(ksign)*(-dz*real(max(kmm,km)-2)+z0coeff*real(ksign)*z0+zp(m))
                           dzzw=z0coeff*dzzw
                        endif
                        distance_pen=sqrt(dxxw**2+dyyw**2+dzzw**2)
                        if(dxyzm.lt.distance_pen)then
                           dxyzm=distance_pen
                           immr=imm
                           jmmr=jmm
                           kmmr=kmm
                           isignr=isign
                           jsignr=jsign
                           ksignr=ksign
                           dxxwr=dxxw
                           dyywr=dyyw
                           dzzwr=dzzw
                        endif
                     endif
                  endif
               enddo
            enddo
         enddo
         if(dxyzm>0.)then
            dxxw=dxxwr
            dyyw=dyywr
            dzzw=dzzwr
            isign=isignr
            jsign=jsignr
            ksign=ksignr
            uran(m)=-uran(m)
            vran(m)=-vran(m)
            wran(m)=-wran(m)
            xp(m)=xp(m)-2.*real(isign)*dxxw
            yp(m)=yp(m)-2.*real(jsign)*dyyw
            zp(m)=zp(m)-2.*dzzw*real(ksign)
            imm=int((xp(m))/dx)+1
            jmm=int((yp(m))/dy)+1
            kmm=int((zp(m)+dz)/dz)+1
            icheck=0
            ikcheck=0
            iicheck=0
            ijcheck=0
            if(kmm.lt.1.or.kmm.gt.nz-1)ikcheck=1
            if(imm.lt.1.or.imm.gt.nx-1)iicheck=1
            if(jmm.lt.1.or.jmm.gt.ny-1)ijcheck=1
!            if(icellflag(imm,jmm,kmm) .eq. 0)icheck=1
            i2check=icheck
         endif
         return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine detang(iomega)
         use variables   !make data from module "variables" visible
         implicit none
         integer, intent(in) :: iomega
         real e11,e12,e13,e21,e22,e23,e31,e32,e33
         if(sqrt(u(i,j,k)**2+v(i,j,k)**2).gt.1.e-05)then
            cospsi=u(i,j,k)/(sqrt(u(i,j,k)**2+v(i,j,k)**2))
            sinpsi=v(i,j,k)/(sqrt(u(i,j,k)**2+v(i,j,k)**2))
            sinphiw=w(i,j,k)/(sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)+1.e-10)
            cosphiw=sqrt(u(i,j,k)**2+v(i,j,k)**2)/(sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)+1.e-10)
            omenum=-dutotdxi(i,j,k)*sinpsi+dutotdyi(i,j,k)*cospsi
            omeden=dutotdxi(i,j,k)*cospsi*sinphiw+dutotdyi(i,j,k)*sinpsi*sinphiw-dutotdzi(i,j,k)*cosphiw
            if(iomega.eq.0)then
               if(abs(omeden).lt.1.e-10)then
                  cosomeg=0.
                  sinomeg=1.
                  dutotdn=dutotdxi(i,j,k)*(sinpsi*sinomeg-cospsi*sinphiw*cosomeg) &
                       -dutotdyi(i,j,k)*(cospsi*sinomeg+sinpsi*sinphiw*cosomeg) &
                       +dutotdzi(i,j,k)*cosomeg*cosphiw
                  if(dutotdn.lt.0.)sinomeg=-1.
               else
                  tanomeg=omenum/omeden
                  omeg1=atan(tanomeg)
                  cosomeg1=cos(omeg1)
                  sinomeg1=sin(omeg1)
                  dutotdn1=dutotdxi(i,j,k)*(sinpsi*sinomeg1-cospsi*sinphiw*cosomeg1) &
                        -dutotdyi(i,j,k)*(cospsi*sinomeg1+sinpsi*sinphiw*cosomeg1) &
                        +dutotdzi(i,j,k)*cosomeg1*cosphiw
                  omeg2=omeg1+pi
                  cosomeg2=cos(omeg2)
                  sinomeg2=sin(omeg2)
                  dutotdn2=dutotdxi(i,j,k)*(sinpsi*sinomeg2-cospsi*sinphiw*cosomeg2) &
                        -dutotdyi(i,j,k)*(cospsi*sinomeg2+sinpsi*sinphiw*cosomeg2) &
                        +dutotdzi(i,j,k)*cosomeg2*cosphiw
                  if(dutotdn2.gt.dutotdn1)then
                     dutotdn=dutotdn2
                     omeg=omeg2
                     cosomeg=cosomeg2
                     sinomeg=sinomeg2
                  else
                     dutotdn=dutotdn1
                     omeg=omeg1
                     cosomeg=cosomeg1
                     sinomeg=sinomeg1
                  endif

                  endif
            else
               if(iomega.eq.1)then
                  omeg=0.
                  cosomeg=1.
                  sinomeg=0.
               else
                  if(abs(cospsi).gt.0.5)then
                     omeg=pi/2.
                     cosomeg=0.
                     sinomeg=-sign(cospsi,1.0)
                     if(iomega.eq.3)then
                        sinomeg=sign(cospsi,1.0)
                        omeg=-omeg
                     endif
                  else
                     return
                  endif
               endif
            endif
            alph1ij(i,j,k)=cospsi*cosphiw
            alph2ij(i,j,k)=-sinpsi*cosomeg-cospsi*sinphiw*sinomeg
            alph3ij(i,j,k)=sinpsi*sinomeg-cospsi*sinphiw*cosomeg
            bet1ij(i,j,k)=sinpsi*cosphiw
            bet2ij(i,j,k)=cospsi*cosomeg-sinpsi*sinphiw*sinomeg
            bet3ij(i,j,k)=-cospsi*sinomeg-sinpsi*sinphiw*cosomeg
            gam1ij(i,j,k)=sinphiw
            gam2ij(i,j,k)=cosphiw*sinomeg
            gam3ij(i,j,k)=cosphiw*cosomeg
            alphn1ij(i,j,k)=cospsi*cosphiw
            alphn2ij(i,j,k)=sinpsi*cosphiw
            alphn3ij(i,j,k)=sinphiw
            betn1ij(i,j,k)=-sinpsi*cosomeg-cospsi*sinphiw*sinomeg
            betn2ij(i,j,k)=cospsi*cosomeg-sinpsi*sinphiw*sinomeg
            betn3ij(i,j,k)=cosphiw*sinomeg
            gamn1ij(i,j,k)=sinpsi*sinomeg-cospsi*sinphiw*cosomeg
            gamn2ij(i,j,k)=-cospsi*sinomeg-sinpsi*sinphiw*cosomeg
            gamn3ij(i,j,k)=cosphiw*cosomeg
            dutotdn=dutotdxi(i,j,k)*(sinpsi*sinomeg-cospsi*sinphiw*cosomeg) &
                 -dutotdyi(i,j,k)*(cospsi*sinomeg+sinpsi*sinphiw*cosomeg) &
                 +dutotdzi(i,j,k)*cosomeg*cosphiw
            dutotds=dutotdxi(i,j,k)*cospsi*cosphiw+dutotdyi(i,j,k)*&
                    sinpsi*cosphiw+dutotdzi(i,j,k)*sinphiw
            dutotdni(i,j,k)=dutotdn
            dutotdsi(i,j,k)=dutotds
!         asi(i,j,k)=cospsi*cosphi
!         bsi(i,j,k)=sinpsi*cosphi
!         csi(i,j,k)=sinphi
            ani(i,j,k)=(sinpsi*sinomeg-cospsi*sinphiw*cosomeg)
            bni(i,j,k)=-(cospsi*sinomeg+sinpsi*sinphiw*cosomeg)
            cni(i,j,k)=cosomeg*cosphiw
         else
            if(abs(w(i,j,k)).lt.1.e-05)then
               if(sqrt(dutotdxi(i,j,k)**2+dutotdyi(i,j,k)**2).gt.0.)then
                  cospsi=dutotdxi(i,j,k)/(sqrt(dutotdxi(i,j,k)**2+dutotdyi(i,j,k)**2))
                  sinpsi=dutotdyi(i,j,k)/(sqrt(dutotdxi(i,j,k)**2+dutotdyi(i,j,k)**2))
               else
                  cospsi=1.
                  sinpsi=0.
               endif
               if(sqrt(dutotdxi(i,j,k)**2+dutotdyi(i,j,k)**2+dutotdzi(i,j,k)**2).gt.0)then
                  cosphiw=dutotdzi(i,j,k)/(sqrt(dutotdxi(i,j,k)**2&
                          +dutotdyi(i,j,k)**2+dutotdzi(i,j,k)**2))
                  sinphiw=sqrt(dutotdxi(i,j,k)**2+dutotdyi(i,j,k)**2)&
                          /(sqrt(dutotdxi(i,j,k)**2+dutotdyi(i,j,k)**2+dutotdzi(i,j,k)**2))
               else
                  cosphiw=1.
                  sinphiw=0.
               endif
               alphn1ij(i,j,k)=cospsi*cosphiw
               alphn2ij(i,j,k)=sinpsi*cosphiw
               alphn3ij(i,j,k)=-sinphiw
               betn1ij(i,j,k)=-sinpsi
               betn2ij(i,j,k)=cospsi
               betn3ij(i,j,k)=0.
               gamn1ij(i,j,k)=sinphiw*cospsi
               gamn2ij(i,j,k)=sinphiw*sinpsi
               gamn3ij(i,j,k)=cosphiw
               alph1ij(i,j,k)=cospsi*cosphiw
               alph2ij(i,j,k)=-sinpsi
               alph3ij(i,j,k)=sinphiw*cospsi
               bet1ij(i,j,k)=sinpsi*cosphiw
               bet2ij(i,j,k)=cospsi
               bet3ij(i,j,k)=sinphiw*sinpsi
               gam1ij(i,j,k)=-sinphiw
               gam2ij(i,j,k)=0.
               gam3ij(i,j,k)=cosphiw
               dutotdni(i,j,k)=sqrt(dutotdxi(i,j,k)**2+dutotdyi(i,j,k)**2+dutotdzi(i,j,k)**2)
               dutotdsi(i,j,k)=0.
               ani(i,j,k)=sinphiw*cospsi
               bni(i,j,k)=sinphiw*sinpsi
               cni(i,j,k)=cosphiw
            else
               if(w(i,j,k).gt.0.)then
                  if(sqrt(dutotdxi(i,j,k)**2+dutotdyi(i,j,k)**2).gt.0.)then
                     cospsi=dutotdxi(i,j,k)/sqrt(dutotdxi(i,j,k)**2+dutotdyi(i,j,k)**2)
                     sinpsi=dutotdyi(i,j,k)/sqrt(dutotdxi(i,j,k)**2+dutotdyi(i,j,k)**2)
                  else
                     cospsi=1.
                     sinpsi=0.
                  endif
                  alph1ij(i,j,k)=0.
                  alph2ij(i,j,k)=sinpsi
                  alph3ij(i,j,k)=cospsi
                  bet1ij(i,j,k)=0.
                  bet2ij(i,j,k)=-cospsi
                  bet3ij(i,j,k)=sinpsi
                  gam1ij(i,j,k)=1.
                  gam2ij(i,j,k)=0.
                  gam3ij(i,j,k)=0.
                  alphn1ij(i,j,k)=0.
                  alphn2ij(i,j,k)=0.
                  alphn3ij(i,j,k)=1.
                  betn1ij(i,j,k)=sinpsi
                  betn2ij(i,j,k)=-cospsi
                  betn3ij(i,j,k)=0.
                  gamn1ij(i,j,k)=cospsi
                  gamn2ij(i,j,k)=sinpsi
                  gamn3ij(i,j,k)=0.
                  dutotdni(i,j,k)=sqrt(dutotdxi(i,j,k)**2+dutotdyi(i,j,k)**2)
                  dutotdni(i,j,k)=max(1.e-12,dutotdni(i,j,k))
                  dutotdsi(i,j,k)=dutotdzi(i,j,k)
                  ani(i,j,k)=cospsi
                  bni(i,j,k)=sinpsi
                  cni(i,j,k)=0.
               else
                  if(sqrt(dutotdxi(i,j,k)**2+dutotdyi(i,j,k)**2).gt.0.)then
                     cospsi=dutotdxi(i,j,k)/sqrt(dutotdxi(i,j,k)**2+dutotdyi(i,j,k)**2)
                     sinpsi=dutotdyi(i,j,k)/sqrt(dutotdxi(i,j,k)**2+dutotdyi(i,j,k)**2)
                  else
                     cospsi=1.
                     sinpsi=0.
                  endif
                  alphn1ij(i,j,k)=0.
                  alphn2ij(i,j,k)=0.
                  alphn3ij(i,j,k)=-1.
                  betn1ij(i,j,k)=-sinpsi
                  betn2ij(i,j,k)=cospsi
                  betn3ij(i,j,k)=0.
                  gamn1ij(i,j,k)=cospsi
                  gamn2ij(i,j,k)=sinpsi
                  gamn3ij(i,j,k)=0.
                  alph1ij(i,j,k)=0.
                  alph2ij(i,j,k)=-sinpsi
                  alph3ij(i,j,k)=cospsi
                  bet1ij(i,j,k)=0.
                  bet2ij(i,j,k)=cospsi
                  bet3ij(i,j,k)=sinpsi
                  gam1ij(i,j,k)=-1.
                  gam2ij(i,j,k)=0.
                  gam3ij(i,j,k)=0.
                  dutotdni(i,j,k)=sqrt(dutotdxi(i,j,k)**2+dutotdyi(i,j,k)**2)
                  dutotdsi(i,j,k)=-dutotdzi(i,j,k)
               endif
            endif
         endif
         e11=alph1ij(i,j,k)*alphn1ij(i,j,k)+alph2ij(i,j,k)*betn1ij(i,j,k)+alph3ij(i,j,k)*gamn1ij(i,j,k)
         e12=alph1ij(i,j,k)*alphn2ij(i,j,k)+alph2ij(i,j,k)*betn2ij(i,j,k)+alph3ij(i,j,k)*gamn2ij(i,j,k)
         e13=alph1ij(i,j,k)*alphn3ij(i,j,k)+alph2ij(i,j,k)*betn3ij(i,j,k)+alph3ij(i,j,k)*gamn3ij(i,j,k)
         e21=bet1ij(i,j,k)*alphn1ij(i,j,k)+bet2ij(i,j,k)*betn1ij(i,j,k)+bet3ij(i,j,k)*gamn1ij(i,j,k)
         e22=bet1ij(i,j,k)*alphn2ij(i,j,k)+bet2ij(i,j,k)*betn2ij(i,j,k)+bet3ij(i,j,k)*gamn2ij(i,j,k)
         e23=bet1ij(i,j,k)*alphn3ij(i,j,k)+bet2ij(i,j,k)*betn3ij(i,j,k)+bet3ij(i,j,k)*gamn3ij(i,j,k)
         e31=gam1ij(i,j,k)*alphn1ij(i,j,k)+gam2ij(i,j,k)*betn1ij(i,j,k)+gam3ij(i,j,k)*gamn1ij(i,j,k)
         e32=gam1ij(i,j,k)*alphn2ij(i,j,k)+gam2ij(i,j,k)*betn2ij(i,j,k)+gam3ij(i,j,k)*gamn2ij(i,j,k)
         e33=gam1ij(i,j,k)*alphn3ij(i,j,k)+gam2ij(i,j,k)*betn3ij(i,j,k)+gam3ij(i,j,k)*gamn3ij(i,j,k)

         return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine rotu3psq
!this subroutine rotates the fluctuating quanitities back into the normal
! coordinate sytem
         use variables   !make data from module "variables" visible
         implicit none
         ufsqb=u3psq*alph1ij(i,j,k)**2+utot**2*alph1ij(i,j,k)**2+2.*upvp*alph1ij(i,j,k)*alph2ij(i,j,k) &
               +2.*upwp*alph1ij(i,j,k)*alph3ij(i,j,k)+v3psq*alph2ij(i,j,k)**2 &
               +2.*vpwp*alph2ij(i,j,k)*alph3ij(i,j,k)+w3psq*alph3ij(i,j,k)**2-u(i,j,k)**2
         wfsqb=u3psq*gam1ij(i,j,k)**2+utot**2*gam1ij(i,j,k)**2+2.*upvp*gam1ij(i,j,k)*gam2ij(i,j,k) &
               +2.*upwp*gam1ij(i,j,k)*gam3ij(i,j,k)+v3psq*gam2ij(i,j,k)**2+ &
               2.*vpwp*gam2ij(i,j,k)*gam3ij(i,j,k)+w3psq*gam3ij(i,j,k)**2-w(i,j,k)**2
         vfsqb=u3psq*bet1ij(i,j,k)**2+utot**2*bet1ij(i,j,k)**2+2.*upvp*bet1ij(i,j,k)*bet2ij(i,j,k) &
               +2.*upwp*bet1ij(i,j,k)*bet3ij(i,j,k)+v3psq*bet2ij(i,j,k)**2 &
               +2.*vpwp*bet2ij(i,j,k)*bet3ij(i,j,k)+w3psq*bet3ij(i,j,k)**2-v(i,j,k)**2
         ufvf=u3psq*alph1ij(i,j,k)*bet1ij(i,j,k)+utot**2*alph1ij(i,j,k)*bet1ij(i,j,k) &
              +upvp*(alph1ij(i,j,k)*bet2ij(i,j,k)+alph2ij(i,j,k)*bet1ij(i,j,k)) &
              +upwp*(alph1ij(i,j,k)*bet3ij(i,j,k)+alph3ij(i,j,k)*bet1ij(i,j,k)) &
              +v3psq*alph2ij(i,j,k)*bet2ij(i,j,k)+vpwp*(alph2ij(i,j,k)*bet3ij(i,j,k) &
              +alph3ij(i,j,k)*bet2ij(i,j,k))+w3psq*alph3ij(i,j,k)*bet3ij(i,j,k) &
              -u(i,j,k)*v(i,j,k)
         ufwf=u3psq*alph1ij(i,j,k)*gam1ij(i,j,k)+utot**2*alph1ij(i,j,k)*gam1ij(i,j,k) &
              +upvp*(alph1ij(i,j,k)*gam2ij(i,j,k)+alph2ij(i,j,k)*gam1ij(i,j,k)) &
              +upwp*(alph1ij(i,j,k)*gam3ij(i,j,k)+alph3ij(i,j,k)*gam1ij(i,j,k)) &
              +v3psq*alph2ij(i,j,k)*gam2ij(i,j,k)+vpwp*(alph2ij(i,j,k)*gam3ij(i,j,k) &
              +alph3ij(i,j,k)*gam2ij(i,j,k))+w3psq*alph3ij(i,j,k)*gam3ij(i,j,k) &
              -u(i,j,k)*w(i,j,k)
         vfwf=u3psq*bet1ij(i,j,k)*gam1ij(i,j,k)+utot**2*bet1ij(i,j,k)*gam1ij(i,j,k)+upvp* &
              (bet1ij(i,j,k)*gam2ij(i,j,k)+bet2ij(i,j,k)*gam1ij(i,j,k))+upwp*(bet1ij(i,j,k)*gam3ij(i,j,k) &
              +bet3ij(i,j,k)*gam1ij(i,j,k))+v3psq*bet2ij(i,j,k)*gam2ij(i,j,k) &
              +vpwp*(bet2ij(i,j,k)*gam3ij(i,j,k)+bet3ij(i,j,k)*gam2ij(i,j,k)) &
              +w3psq*bet3ij(i,j,k)*gam3ij(i,j,k)-v(i,j,k)*w(i,j,k)
         return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine rotufsq
         use variables    !make data from module "variables" visible
         implicit none
         real e11,e12,e13,e21,e22,e23,e31,e32,e33
         u3psq=ufsq*alphn1ij(i,j,k)**2+u(i,j,k)**2*alphn1ij(i,j,k)**2 &
               +2.*ufvf*alphn1ij(i,j,k)*alphn2ij(i,j,k)+2.*u(i,j,k)*v(i,j,k)*alphn1ij(i,j,k) &
               *alphn2ij(i,j,k)+2.*ufwf*alphn1ij(i,j,k)*alphn3ij(i,j,k)+2.*u(i,j,k)*w(i,j,k) &
               *alphn1ij(i,j,k)*alphn3ij(i,j,k)+vfsq*alphn2ij(i,j,k)**2+v(i,j,k)**2 &
               *alphn2ij(i,j,k)**2+2.*vfwf*alphn2ij(i,j,k)*alphn3ij(i,j,k)+ &
               2.*v(i,j,k)*w(i,j,k)*alphn2ij(i,j,k)*alphn3ij(i,j,k)+wfsq &
               *alphn3ij(i,j,k)**2+w(i,j,k)**2*alphn3ij(i,j,k)**2-(u(i,j,k)**2+v(i,j,k)**2 &
               +w(i,j,k)**2)
         v3psq=ufsq*betn1ij(i,j,k)**2+u(i,j,k)**2*betn1ij(i,j,k)**2+2.*ufvf &
               *betn1ij(i,j,k)*betn2ij(i,j,k)+2.*u(i,j,k)*v(i,j,k)*betn1ij(i,j,k) &
               *betn2ij(i,j,k)+2.*ufwf*betn1ij(i,j,k)*betn3ij(i,j,k)+2.*u(i,j,k) &
               *w(i,j,k)*betn1ij(i,j,k)*betn3ij(i,j,k)+vfsq*betn2ij(i,j,k)**2 &
               +v(i,j,k)**2*betn2ij(i,j,k)**2+2.*vfwf*betn2ij(i,j,k)*betn3ij(i,j,k) &
               +2.*v(i,j,k)*w(i,j,k)*betn2ij(i,j,k)*betn3ij(i,j,k)+wfsq*betn3ij(i,j,k)**2 &
               +w(i,j,k)**2*betn3ij(i,j,k)**2
         w3psq=ufsq*gamn1ij(i,j,k)**2+u(i,j,k)**2*gamn1ij(i,j,k)**2+2.*ufvf &
               *gamn1ij(i,j,k)*gamn2ij(i,j,k)+2.*u(i,j,k)*v(i,j,k)*gamn1ij(i,j,k) &
               *gamn2ij(i,j,k)+2.*ufwf*gamn1ij(i,j,k)*gamn3ij(i,j,k)+2.*u(i,j,k) &
               *w(i,j,k)*gamn1ij(i,j,k)*gamn3ij(i,j,k)+vfsq*gamn2ij(i,j,k)**2 &
               +v(i,j,k)**2*gamn2ij(i,j,k)**2+2.*vfwf*gamn2ij(i,j,k)*gamn3ij(i,j,k) &
               +2.*v(i,j,k)*w(i,j,k)*gamn2ij(i,j,k)*gamn3ij(i,j,k)+wfsq*gamn3ij(i,j,k)**2 &
               +w(i,j,k)**2*gamn3ij(i,j,k)**2
         upwp=ufsq*alphn1ij(i,j,k)*gamn1ij(i,j,k)+u(i,j,k)**2*alphn1ij(i,j,k)* &
              gamn1ij(i,j,k)+ufvf*(alphn1ij(i,j,k)*gamn2ij(i,j,k)+alphn2ij(i,j,k)* &
              gamn1ij(i,j,k))+u(i,j,k)*v(i,j,k)*(alphn1ij(i,j,k)*gamn2ij(i,j,k)+ &
              alphn2ij(i,j,k)*gamn1ij(i,j,k))+ufwf*(alphn1ij(i,j,k)*gamn3ij(i,j,k) &
              +alphn3ij(i,j,k)*gamn1ij(i,j,k))+u(i,j,k)*w(i,j,k)*(alphn1ij(i,j,k) &
              *gamn3ij(i,j,k)+alphn3ij(i,j,k)*gamn1ij(i,j,k))+vfsq*alphn2ij(i,j,k) &
              *gamn2ij(i,j,k)+v(i,j,k)**2*alphn2ij(i,j,k)*gamn2ij(i,j,k) &
              +vfwf*(alphn2ij(i,j,k)*gamn3ij(i,j,k)+alphn3ij(i,j,k)*gamn2ij(i,j,k)) &
              +v(i,j,k)*w(i,j,k)*(alphn2ij(i,j,k)*gamn3ij(i,j,k)+alphn3ij(i,j,k) &
              *gamn2ij(i,j,k))+wfsq*alphn3ij(i,j,k)*gamn3ij(i,j,k)+w(i,j,k)**2 &
              *alphn3ij(i,j,k)*gamn3ij(i,j,k)
         e11=alph1ij(i,j,k)*alphn1ij(i,j,k)+alph2ij(i,j,k)*betn1ij(i,j,k)+alph3ij(i,j,k)*gamn1ij(i,j,k)
         e12=alph1ij(i,j,k)*alphn2ij(i,j,k)+alph2ij(i,j,k)*betn2ij(i,j,k)+alph3ij(i,j,k)*gamn2ij(i,j,k)
         e13=alph1ij(i,j,k)*alphn3ij(i,j,k)+alph2ij(i,j,k)*betn3ij(i,j,k)+alph3ij(i,j,k)*gamn3ij(i,j,k)
         e21=bet1ij(i,j,k)*alphn1ij(i,j,k)+bet2ij(i,j,k)*betn1ij(i,j,k)+bet3ij(i,j,k)*gamn1ij(i,j,k)
         e22=bet1ij(i,j,k)*alphn2ij(i,j,k)+bet2ij(i,j,k)*betn2ij(i,j,k)+bet3ij(i,j,k)*gamn2ij(i,j,k)
         e23=bet1ij(i,j,k)*alphn3ij(i,j,k)+bet2ij(i,j,k)*betn3ij(i,j,k)+bet3ij(i,j,k)*gamn3ij(i,j,k)
         e31=gam1ij(i,j,k)*alphn1ij(i,j,k)+gam2ij(i,j,k)*betn1ij(i,j,k)+gam3ij(i,j,k)*gamn1ij(i,j,k)
         e32=gam1ij(i,j,k)*alphn2ij(i,j,k)+gam2ij(i,j,k)*betn2ij(i,j,k)+gam3ij(i,j,k)*gamn2ij(i,j,k)
         e33=gam1ij(i,j,k)*alphn3ij(i,j,k)+gam2ij(i,j,k)*betn3ij(i,j,k)+gam3ij(i,j,k)*gamn3ij(i,j,k)

         return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine rotate2d(ic,jc,kc)
! this subroutine rotates variables from primed system aligned with the overall wind
! into the regular grid system
         use variables   !make data from module "variables" visible
         implicit none
         integer, intent(in) :: ic,jc,kc
         ub=u(ic,jc,kc)
         vb=v(ic,jc,kc)
         wb=w(ic,jc,kc)
         upb=ub*cosphi+vb*sinphi
         vpb=-ub*sinphi+vb*cosphi
         ufsqgi(ic,jc,kc)=upsqg*cosphi**2+upb**2*cosphi**2-2.*cosphi*sinphi*upvpg &
                          -2.*cosphi*sinphi*upb*vpb+vpsqg*sinphi**2+sinphi**2*vpb**2-ub**2
         vfsqgi(ic,jc,kc)=upsqg*sinphi**2+upb**2*sinphi**2+2.*upvpg*sinphi*cosphi &
                          +2.*upb*vpb*sinphi*cosphi+vpsqg*cosphi**2+vpb**2*cosphi**2-vb**2
         wfsqgi(ic,jc,kc)=wpsqg
         ufvfgi(ic,jc,kc)=upsqg*cosphi*sinphi+upb**2*cosphi*sinphi+upvpg*(cosphi**2 &
                          -sinphi**2)+upb*vpb*(cosphi**2-sinphi**2)-vpsqg*sinphi*cosphi-vpb**2*sinphi*cosphi-ub*vb
         ufwfgi(ic,jc,kc)=upwpg*cosphi+upb*wb*cosphi-vpwpg*sinphi-vpb*wb*sinphi-ub*wb
         vfwfgi(ic,jc,kc)=upwpg*sinphi+upb*wb*sinphi+vpwpg*cosphi+vpb*wb*cosphi-vb*wb
         return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine int3cal(dp,sc,taud,rhop,xnu,press,temp,diff,vs)
         implicit none
         real dp,sc,taud,rhop,xnu,press,temp,diff,vs,rhoair,xmu,a,pi,dpm
         real cdresq,re
         real rearray(180)
         integer i,ie
         pi=4.*atan2(1.,1.)
         data (rearray(i),i=1,100)/0.40,0.43,0.47,0.51,0.54,0.58,0.61,0.65,0.69,0.72 &
                                  ,0.75,0.79,0.82,0.86,0.89,0.92,0.96,0.99,1.02,1.06 &
                                  ,1.09,1.12,1.15,1.18,1.22,1.25,1.28,1.31,1.34,1.37 &
                                  ,1.40,1.43,1.47,1.50,1.53,1.56,1.59,1.62,1.65,1.68 &
                                  ,1.71,1.74,1.77,1.80,1.83,1.86,1.89,1.92,1.95,1.98 &
                                  ,2.01,2.04,2.07,2.10,2.13,2.16,2.19,2.22,2.25,2.27 &
                                  ,2.30,2.33,2.36,2.39,2.42,2.45,2.47,2.50,2.53,2.56 &
                                  ,2.59,2.62,2.64,2.67,2.70,2.73,2.75,2.78,2.81,2.84 &
                                  ,2.86,2.89,2.92,2.95,2.97,3.00,3.03,3.05,3.08,3.11 &
                                  ,3.14,3.40,3.66,3.92,4.17,4.41,4.66,4.90,5.13,5.36/
         data (rearray(i), i=101,180)/5.59,5.82,6.05,6.27,6.49,6.70,6.92,7.13,7.34,7.55 &
                                     ,7.75,7.96,8.16,8.36,8.56,8.75,8.95,9.14,9.34,9.53 &
                                     ,9.72,9.90,10.1,10.3,10.5,10.6,10.8,11.0,11.2,11.4 &
                                     ,11.5,11.7,11.9,12.1,12.2,12.4,12.6,12.8,12.9,13.1 &
                                     ,13.3,13.4,13.6,13.8,13.9,14.1,14.2,14.4,14.6,14.7 &
                                     ,14.9,15.0,15.2,15.4,15.5,15.7,15.8,16.0,16.1,16.3 &
                                     ,16.4,16.6,16.7,16.9,17.1,17.2,17.3,17.5,17.6,17.8 &
                                     ,17.9,18.1,18.2,18.4,18.5,18.7,18.8,19.0,19.1,19.2/
! this data represents the Reynolds number as a function of the
! product of drag coefficient and Reynolds number squared which
! is not a function of the settling velocity see table 3-4 of Aerosol
! Technology by William C. Hinds - John Wiley & Sons, 1982 p. 54.
         dpm=dp*1.e-04
! convert microns to cm
         rhoair=28.8*(273./temp)*(press/76.)/22400.
         xmu=1.78e-04
         xnu=xmu/rhoair
         a=dpm/2
         if(dpm.gt.0.)then
            diff=1.32e-16*temp*(1.+(1.e-04/(press*a))*(6.32+2.01*exp(-2190*press*a)))/(6.*pi*xmu*a)
            taud=rhop*dpm**2/(18.*xmu)
            vs=taud*980.
            cdresq=4.*rhoair*rhop*980.*dpm**3/(3*xmu**2)
            if(cdresq.ge.10.)then
               if(cdresq.ge.100)then
                  ie=int(81.+cdresq/10)
                  ie=min(ie,180)
               else
                  ie=int(cdresq-9)
               endif
               Re=rearray(ie)
               vs=Re*xmu/(dpm*rhoair)
            endif
         else
            diff=.198
            taud=0.
            vs=0.
         endif
         sc=xnu/diff
!        r=1./(sc**(-2./3.)+10.**(-3./taup))
!        r=r/ustar
!This approach follows Atmospheric Chemistry and Physics of Air
! Pollution by John H. Seinfeld - John Wiley & Sons, 1986 pages 640-643
         return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine readwind
! readwind reads in the next wind field and computes initial vertical
! profiles
         use variables
         implicit none
! read in the wind fields using the same grid definitions as qwic-urb
! note that there are no values for the n indices; the active grid
! for the particles requires all indicies be>0 and less than n-1
         integer itop,jtop
         if(format_flag.eq.1)then   !erp 3/2/05
            do k=1,nz-1
               do j=1,ny-1
                  do i=1,nx-1
                     read(27,*)xi1,yj1,zk1,u(i,j,k),v(i,j,k),w(i,j,k)
                     sigvi(i,j,k)=0.
                     sigwi(i,j,k)=0.
                     ustarij(i,j,k)=0.
                     ustarz(i,j,k)=0.
                     ustarg(i,j,k)=0.
                     ustargz(i,j,k)=0.
                     elzg(i,j,k)=real(nx)*dx
                     dutotdyi(i,j,k)=0.
                     eleff(i,j,k)=0.
!                    elx(i,j,k)=0.
!                    ely(i,j,k)=0.
                     dutotdzi(i,j,k)=0.
                     dutotdni(i,j,k)=0.
                     dutotdsi(i,j,k)=0.
                     dxp(i,j,k)=real(nx)*dx ! distances to walls and floor
                     dxm(i,j,k)=real(nx)*dx
                     dym(i,j,k)=real(ny)*dy
                     dyp(i,j,k)=real(ny)*dy
                     dzm(i,j,k)=real(nz)*dz
!                    depcx(i,j,k)=0.
!                    depcy(i,j,k)=0.
!                    depcz(i,j,k)=0.
                  enddo
               enddo
            enddo
         else
            do k=1,nz-1
               do j=1,ny-1
                  do i=1,nx-1
                     sigvi(i,j,k)=0.
                     sigwi(i,j,k)=0.
                     ustarij(i,j,k)=0.
                     ustarz(i,j,k)=0.
                     ustarg(i,j,k)=0.
                     ustargz(i,j,k)=0.
                     elzg(i,j,k)=real(nx)*dx
                     dutotdyi(i,j,k)=0.
                     eleff(i,j,k)=0.
                     dutotdzi(i,j,k)=0.
                     dutotdni(i,j,k)=0.
                     dutotdsi(i,j,k)=0.
                     dxp(i,j,k)=real(nx)*dx ! distances to walls and floor
                     dxm(i,j,k)=real(nx)*dx
                     dym(i,j,k)=real(ny)*dy
                     dyp(i,j,k)=real(ny)*dy
                     dzm(i,j,k)=real(nz)*dz
                  enddo
               enddo
            enddo
            allocate(bin_read1((nx-1)*(ny-1)*(nz-1)),bin_read2((nx-1)*(ny-1)*(nz-1)),bin_read3((nx-1)*(ny-1)*(nz-1)))
            read(50)bin_read1,bin_read2,bin_read3
            idummy=1
            do k=1,nz-1
               do j=1,ny-1
                  do i=1,nx-1
                     u(i,j,k)=bin_read1(idummy)
                     v(i,j,k)=bin_read2(idummy)
                     w(i,j,k)=bin_read3(idummy)
                     idummy=idummy+1
                  enddo
               enddo
            enddo
            deallocate(bin_read1,bin_read2,bin_read3)
         endif

!mdw 4-16-2004 added coding to describe the winds at
!the top of the domain on the upwind boundary
         if(theta.lt.22.5.or.theta.ge.337.5)then
            itop=(nx-1)/2
            jtop=ny-1
         else
            if(theta.ge.22.5.and.theta.lt.67.5)then
               itop=nx-1
               jtop=ny-1
            else
               if(theta.ge.67.5.and.theta.lt.112.5)then
                  itop=nx-1
                  jtop=(ny-1)/2
               else
                  if(theta.ge.112.5.and.theta.lt.167.5)then
                     itop=nx-1
                     jtop=1
                  else
                     if(theta.ge.167.5.and.theta.lt.202.5)then
                        itop=(nx-1)/2
                        jtop=1
                     else
                        if(theta.ge.202.5.and.theta.lt.247.5)then
                           itop=1
                           jtop=1
                        else
                           if(theta.ge.247.5.and.theta.lt.292.5)then
                              itop=1
                              jtop=(ny-1)/2
                           else
                              if(theta.ge.292.5.and.theta.lt.337.5)then
                                 itop=1
                                 jtop=ny-1
                              endif
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
         ualoft=u(itop,jtop,nz-1)
         valoft=v(itop,jtop,nz-1)
         return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine readbldgout
         use variables
         implicit none
! mdw 2-03-2004 revised character statements for reading in buildout.dat
         character*5 adum51,adum52
         character*6 adum6,adum61,adum611
         character*7 adum7,adum71,adum72,adum711,adum712,adum713
         character*8 adum8,adum81
         character*9 adum9,adum91
         character*18 adum18
         integer in
         if(time_idx .eq. 1)then
            read(44,*)inumbuild
            read(44,*)inumveg
            read(44,*)inumgarage
            inumbuildfile=inumbuild-inumgarage
            allocate(Ht(inumbuild),Wti(inumbuild),Lt(inumbuild),phib(inumbuild))
            allocate(bldtype(inumbuild),bldnum(inumbuild))
            allocate(xfo(inumbuild),yfo(inumbuild),zfo(inumbuild),ustarij(nx,ny,nz),xft(inumbuild))
            allocate(weff(inumbuild),leff(inumbuild),atten(inumbuild),gamma(inumbuild))
            allocate(lr(inumbuild),sx(inumbuild),sy(inumbuild),lfr(inumbuild))
            bldtype(:)=9        
         endif
         do in=1,inumbuild-inumveg
            read(44,4920,end=6)adum18,ibuild
 4921       format(a8,f9.4,a8,f9.4)
 4920       format(a18,i4)
            read(44,4925)adum9,bldtype(ibuild),adum91,gamma(ibuild)
 4925       format(a8,i4,a9,f9.4)
            read(44,4926)adum611,Ht(ibuild),adum51,Wti(ibuild),adum52,Lt(ibuild)
 4926       format(a6,f9.4,a5,f9.4,a5,f9.4)
            read(44,4927)adum711,xfo(ibuild),adum712,yfo(ibuild),adum713,zfo(ibuild)
 4927       format(a7,f9.4,a7,f9.4,a7,f9.4)
            read(44,4921)adum8,weff(ibuild),adum81,Leff(ibuild)
!           leff(in)=dx*leff(in)
!           weff(in)=dx*weff(in)
            read(44,4922)adum7,lfr(ibuild),adum71,lr(ibuild),adum72,atten(ibuild)
 4922       format(a7,f9.4,a7,f9.4,a7,f9.4)
!           wt=wt*dx
            read(44,4923)adum6,Sx(ibuild),adum61,Sy(ibuild)
 4923       format(a6,f9.4,a6,f9.4)
!           xlt=dx*xlt
!           sx(in)=dx*sx(in)
!           sy(in)=dx*sy(in)
!           lfx(in)=lfx(in)
!           lfy(in)=lfy(in)
!           lr(in)=dx*lr(in)
            if(bldtype(ibuild).eq.3)then
               xft(ibuild)=xfo(ibuild)-0.5*Wti(ibuild)
            else
               xft(ibuild)=xfo(ibuild)
!           Lt(ibuild)=Wti(ibuild)
            endif
!            lfr(ibuild)=sqrt(lfx(ibuild)**2+lfy(ibuild)**2)
         enddo
  6      continue
         return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine readQPparams
         use variables
         implicit none
         read(5,*)isourcetypeflag
         select case(isourcetypeflag)
            case(1,3,8)
               ieradflag=0
               ibuoyflag=0
               ibioslurryflag=0
               itwophaseflag=0
            case(2)
               ieradflag=0
               ibuoyflag=-1
               ibioslurryflag=0
               itwophaseflag=0
            case(4)
               ieradflag=0
               ibuoyflag=1
               ibioslurryflag=0
               itwophaseflag=0
            case(5)
               ieradflag=1
               ibuoyflag=0
               ibioslurryflag=0
               itwophaseflag=0
            case(6)
               ieradflag=0
               ibuoyflag=0
               ibioslurryflag=1
               itwophaseflag=0
            case(7)
               ieradflag=0
               ibuoyflag=0
               ibioslurryflag=0
               itwophaseflag=1
         endselect
         read(5,*)isiteflag    !=1 activates sensor siting
         read(5,*)iindoorflag    !if 1 building infiltration is considered
         read(5,*)inextgridflag  !=0 single grid, =1 inner grid, =2 outer grid
         read(5,*)xoff !x and y offsets of the inner grid origin relative to the outer grid origin (meters)
         read(5,*)yoff
         read(5,*)z0       !surface roughness length (meters)
         read(5,*)rcl      !reciprocal monin-obukhov length (inverse meters)
         read(5,*)h        !free atmosphere height (meters)
         read(5,*)inloc    !non-local flag; if inloc=1 non-local mixing will be used
         read(5,*)nparticles    !number of particles released over entire simulation
         read(5,*)nparticlesfactor   !number of particles factor
         read(5,*)partsplitfactor   !number of particles that a particle splits into
         read(5,*)taylor_flag   !Enables Taylor Microscale lower limit to subtimesteps
         read(5,*)delta   !time step in seconds
         delta_ref=delta
         read(5,*)dur !duration of sim (secs), particle output times & concentration output times
         read(5,*)outc !concentration averaging period (secs)
         read(5,*)concstrt    !starting time for concentration averaging (secs)
         read(5,*)outp !Particle output period (secs)
         read(5,*)nbx    !number of boxes in x, y, and z directions
         read(5,*)nby
         read(5,*)nbz
         read(5,*)xbl !lower limit in x for collecting boxes (meters)
         read(5,*)xbu !upper limit in x for collecting boxes (meters)
         read(5,*)ybl !lower limit in y for collecting boxes (meters)
         read(5,*)ybu !upper limit in y for collecting boxes (meters)
         read(5,*)zbl !lower limit in z for collecting boxes (meters)
         read(5,*)zbu !upper limit in z for collecting boxes (meters)
         if(isiteflag.eq.1)then
            outc=dur-delta !Change the concentration averaging period of sensor siting cases
            nbz=1
         endif
         return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine readQPevaporation
         use variables
         implicit none
         integer num_temp_levels,itemplevel
         real, allocatable ::zlevel(:),templevel(:),rel_humid(:)
         
         read(88,*)
         read(88,*)agent_solids_concentration
         read(88,*)agent_eff_solid_density
         read(88,*)
         read(88,*)agent_droplet_fraction
         read(88,*)agent_mfp
         read(88,*)agent_mw_vapor
         read(88,*)agent_mw_liquid
         read(88,*)Pcritical
         read(88,*)Tcritical
         read(88,*)Tboiling
         read(88,*)liquidDensity
         read(88,*)atomic_volume_agent
         read(88,*)surfaceTensionCoef
         read(88,*)
         read(88,*)agent_stick_eff
         read(88,*)thermal_accomodation
         read(88,*)atm_pressure
         read(88,*)num_temp_levels
         allocate(zlevel(num_temp_levels),templevel(num_temp_levels),rel_humid(num_temp_levels))
         read(88,*)
         do itemplevel=1,num_temp_levels
            read(88,*)zlevel(itemplevel),templevel(itemplevel),rel_humid(itemplevel)
         enddo
         allocate(temperature(nx1,ny1,nz1),vaporfield(nx1,ny1,nz1))
         itemplevel=1
         do k=2,nz1
            if(zi(k) .le. zlevel(1) .or. zi(k) .ge. zlevel(num_temp_levels))then
               temperature(:,:,k)=templevel(itemplevel)
               vaporfield(:,:,k)=rel_humid(itemplevel)
            else
               if(zi(k) .ge. zlevel(itemplevel))then
                  itemplevel=itemplevel+1
               endif
               temperature(:,:,k)=((templevel(itemplevel)-templevel(itemplevel-1))/&
                                  (zlevel(itemplevel)-zlevel(itemplevel-1)))*(zi(k)-zlevel(itemplevel-1))&
                                  +templevel(itemplevel-1)
               vaporfield(:,:,k)=((rel_humid(itemplevel)-rel_humid(itemplevel-1))/&
                                   (zlevel(itemplevel)-zlevel(itemplevel-1)))*(zi(k)-zlevel(itemplevel-1))&
                                   +rel_humid(itemplevel-1)
            endif
         enddo
         temperature(:,:,1)=temperature(:,:,2)
         vaporfield(:,:,1)=vaporfield(:,:,2)
         deallocate(zlevel,templevel,rel_humid)
         return
      end
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine readQPsource
         use variables
         implicit none
         real xgu,ygu
         read(37,*)nsources
         read(37,*)nsourcenodes
         allocate(source_x(nsourcenodes),source_y(nsourcenodes),source_z(nsourcenodes))
         allocate(source_l(nsources),source_w(nsources),source_h(nsources),source_weight(nsources))
         allocate(source_r(nsources),source_hemass(nsources),source_strength(nsources))
         allocate(source_density(nsources),source_start(nsources),source_dur(nsources))
         allocate(source_geometry(nsources),source_str_units(nsources),source_release(nsources))
         allocate(source_nodespeed(nsourcenodes),source_nodeaccel(nsourcenodes),source_nodesegpos(nsourcenodes))
         allocate(source_seglength(nsourcenodes),source_startidx(nsources),source_stopidx(nsources),nlinenodes(nsources))
         allocate(source_prr(nsources),source_nodestart(nsourcenodes),source_nodestop(nsourcenodes))
         allocate(hnew(nsources),hnewp(nsources),rnew(nsources),rnewp(nsources))
         allocate(xscen(nsources),yscen(nsources),rhocent(nsources))
         allocate(xscenp(nsources),yscenp(nsources),volnew(nsources))
         allocate(volnewp(nsources),vol0(nsources),reff(nsources),reffp(nsources))
         allocate(vslopex(nsources),vslopey(nsources),dhdt(nsources),cosphidg(nsources))
         allocate(drdt(nsources),drdtp(nsources),sinphidg(nsources),ubar(nsources))
         source_x(:)=0.         
         source_y(:)=0.
         source_z(:)=0.
         source_l(:)=0.
         source_w(:)=0.
         source_h(:)=0.
         source_r(:)=0.
         source_weight(:)=1.
         source_hemass(:)=0.
         source_strength(:)=0.
         source_density(:)=0.
         source_str_units(:)=0
         source_release(:)=0
         source_startidx(:)=0
         source_stopidx(:)=0
         source_prr(:)=0
         source_start(:)=0.
         source_dur(:)=0.
         source_nodespeed(:)=0.
         source_nodeaccel(:)=0.
         source_nodesegpos(:)=0.
         source_seglength(:)=0.
         inode=0
         do isource=1,nsources
            inode=inode+1
            read(37,*)
            read(37,*)
            read(37,*)source_str_units(isource)
            read(37,*)source_strength(isource)
            read(37,*)source_density(isource)
            read(37,*)source_release(isource)
            read(37,*)source_start(isource)
            read(37,*)source_dur(isource)
            if(source_release(isource) .eq. 2)then
               source_dur(isource)=dur-source_start(isource)
            endif
            read(37,*)source_geometry(isource)
            source_nodestart(inode)=source_start(isource)
            select case(source_geometry(isource))
               case(0) !Array of spherical shell sources for sensor siting mode
                  read(37,*)xgl
                  read(37,*)xgu
                  read(37,*)nxg
                  read(37,*)ygl
                  read(37,*)ygu
                  read(37,*)nyg
                  if(nxg .gt. 1)then
                     dxg=(xgu-xgl)/real(nxg-1)
                  else
                     dxg=0.
                  endif
                  if(nyg .gt. 1)then
                     dyg=(ygu-ygl)/real(nyg-1)
                  else
                     dyg=0.
                  endif
                  read(37,*)source_z(inode)
                  read(37,*)source_r(inode)
               case(1) !Spherical shell source
                  if(source_release(isource).eq.2)then
                     source_nodestop(inode)=dur
                  else
                     source_nodestop(inode)=source_start(isource)+source_dur(isource)
                  endif
                  read(37,*)source_x(inode)
                  read(37,*)source_y(inode)
                  read(37,*)source_z(inode)
                  read(37,*)source_r(inode)
                  source_startidx(isource)=inode
                  source_stopidx(isource)=inode
               case(2) !Line source with multiple segments
                  if(source_release(isource).eq.2)then
                     source_nodestop(inode)=dur
                  else
                     source_nodestop(inode)=source_start(isource)+source_dur(isource)
                  endif
                  read(37,*)nlinenodes(isource)
                  read(37,*)
                  read(37,*)source_x(inode),source_y(inode),source_z(inode)
                  source_startidx(isource)=inode
                  do idum=2,nlinenodes(isource)
                     inode=inode+1
                     source_nodestart(inode)=source_start(isource)
                     source_nodestop(inode)=source_start(isource)+source_dur(isource)
                     read(37,*)source_x(inode),source_y(inode),source_z(inode)
                     source_seglength(inode-1)=sqrt(((source_x(inode)-source_x(inode-1))**2.)&
                                                  +((source_y(inode)-source_y(inode-1))**2.)&
                                                  +((source_x(inode)-source_x(inode-1))**2.))
                  enddo
                  source_stopidx(isource)=inode
               case(3) !Cylindrical source
                  if(source_release(isource).eq.2)then
                     source_nodestop(inode)=dur
                  else
                     source_nodestop(inode)=source_start(isource)+source_dur(isource)
                  endif
                  read(37,*)source_x(inode)
                  read(37,*)source_y(inode)
                  read(37,*)source_z(inode)
                  read(37,*)source_r(inode)
                  read(37,*)source_h(inode)
                  source_startidx(isource)=inode
                  source_stopidx(isource)=inode
               case(4) !Explosive source
                  source_nodestop(inode)=source_start(isource)+source_dur(isource)
                  read(37,*)source_x(inode)
                  read(37,*)source_y(inode)
                  read(37,*)source_z(inode)
                  read(37,*)source_hemass(inode)
                  source_h(inode)=92.56*source_hemass(inode)**.25
                  source_r(inode)=3*source_hemass(inode)**.333
                  cldhgt=source_h(inode)
                  source_startidx(isource)=inode
                  source_stopidx(isource)=inode
               case(5) !Rectangular area or volume source
                  if(source_release(isource).eq.2)then
                     source_nodestop(inode)=dur
                  else
                     source_nodestop(inode)=source_start(isource)+source_dur(isource)
                  endif
                  read(37,*)source_x(inode)
                  read(37,*)source_y(inode)
                  read(37,*)source_z(inode)
                  read(37,*)source_l(inode)
                  read(37,*)source_w(inode)
                  read(37,*)source_h(inode)
                  source_startidx(isource)=inode
                  source_stopidx(isource)=inode
               case(6) !Moving point source with multiple path segments
                  read(37,*)nlinenodes(isource)
                  read(37,*)
                  read(37,*)source_x(inode),source_y(inode),source_z(inode),source_nodespeed(inode)
                  source_startidx(isource)=inode
                  do idum=2,nlinenodes(isource)
                     inode=inode+1
                     read(37,*)source_x(inode),source_y(inode),source_z(inode),source_nodespeed(inode)
                     source_seglength(inode-1)=sqrt(((source_x(inode)-source_x(inode-1))**2.)&
                                                  +((source_y(inode)-source_y(inode-1))**2.)&
                                                  +((source_z(inode)-source_z(inode-1))**2.))
                     source_nodestop(inode-1)=source_nodestart(inode-1)+&
                              2*source_seglength(inode-1)/(source_nodespeed(inode)+source_nodespeed(inode-1))
                     source_nodestart(inode)=source_nodestop(inode-1)
                     source_nodeaccel(inode-1)=(source_nodespeed(inode)-source_nodespeed(inode-1))/&
                              (source_nodestop(inode-1)-source_nodestart(inode-1))
                     if(source_nodestop(inode-1).gt.dur)then
                        source_nodestop(inode-1)=dur
                        source_nodestart(inode)=dur+2*delta
                        source_nodestop(inode)=dur+2*delta
                     endif
                  enddo
                  source_stopidx(isource)=inode
               case(7) !Spherical volume source
                  if(source_release(isource).eq.2)then
                     source_nodestop(inode)=dur
                  else
                     source_nodestop(inode)=source_start(isource)+source_dur(isource)
                  endif
                  read(37,*)source_x(inode)
                  read(37,*)source_y(inode)
                  read(37,*)source_z(inode)
                  read(37,*)source_r(inode)
                  source_startidx(isource)=inode
                  source_stopidx(isource)=inode
            endselect
            read(37,*)
         enddo
         do inode=1,nsourcenodes
            if(source_z(inode).lt.2.*z0)source_z(inode)=2.*z0
         enddo
         return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine readQPindoor
         use variables
         implicit none
         integer idum2
!         real deposition_velocity
         character*4,head3
         allocate(tau(inumbuild),racact(inumbuild),racref(inumbuild),tiref(inumbuild))
         allocate(toref(inumbuild),piref(inumbuild),urefb(inumbuild))
         allocate(tiact(inumbuild),rhoref(inumbuild),racsupply(inumbuild))
         allocate(racintake(inumbuild),filtration(inumbuild),area2volume(inumbuild))
         allocate(ainf(inumbuild),binf(inumbuild),abinf(inumbuild))
         allocate(binfex(inumbuild))
         open(unit=28,file='QP_indoor.dat',status='unknown')
         iface1n(:)=0
         iface2n(:)=0
         jface1n(:)=0
         jface2n(:)=0
         kface1n(:)=0
         kface2n(:)=0
         aminusx(:)=1.e-08
         aplusx(:)=1.e-08
         aminusy(:)=1.e-08
         aplusy(:)=1.e-08
         afloor(:)=1.e-08
         aroof(:)=1.e-08
! mdw 3/22/2004 pm set minimum face area for buildings with obscurred faces
         tau(:)=3600./20.
         racact(:)=0.
         taumin=1.e+08
         racref(:)=0.
         tiref(:)=0.
         toref(:)=0.
         piref(:)=0.
         urefb(:)=0.
         tiact(:)=0.
         rhoref(:)=0.
         racsupply(:)=0.
         racintake(:)=0.
         filtration(:)=0.
         area2volume(:)=0.
         ainf(:)=20./3600.
         binf(:)=20./3600.
         abinf(:)=1.
! MAN 02/23/2007 new infiltration parameters and format for QP_indoor.inp
         read(8,*)toact
         read(8,*)piact
         read(8,15)head3
  15     format(a10)
! temperature, and pressure (atm)
         rhoact=.0288*piact*(273./toact)/.0224
!         if(depvel .ge. 0)then
!            deposition_velocity=depvel/100
!         else
!              deposition_velocity=0.005
!         endif
         ioutin=0
         deltilf=min(outc,0.25*taumin)
         noutilf=int(dur/deltilf)+1
         do ibuild=1,inumbuildfile
            if(bldtype(ibuild) .ge. 9)then
               read(8,*)
               cycle
            endif
            read(8,*)idum,idum2,tiact(ibuild),racsupply(ibuild),racintake(ibuild),&
                     filtration(ibuild),area2volume(ibuild),racref(ibuild),&
                     tiref(ibuild),toref(ibuild),piref(ibuild),urefb(ibuild)
!            if(racintake(ibuild) .lt. 0)then
!                if(toact .lt. 273.)then
!                    racintake(ibuild)=0.5
!                elseif(toact .ge. 273. .and. toact .lt. 290.)then
!                    racintake(ibuild)=0.176*toact-47.55
!                elseif(toact .ge. 290. .and. toact .lt. 295.)then
!                    racintake(ibuild)=-0.46*toact+136.2
!                else
!                    racintake(ibuild)=0.5
!                endif
!            endif
!            rhoref(ibuild)=.0288*piref(ibuild)*(273./toref(ibuild))/.0224
!            psref=rhoref(ibuild)*9.812*Ht(ibuild)*abs(tiref(ibuild)-toref(ibuild))/tiref(ibuild)
!            psact=rhoact*9.812*Ht(ibuild)*abs(tiact(ibuild)-toact)/tiact(ibuild) !these expressions
!! give the buoyancy pressures for the reference and actual cases
!            pwref=.5*rhoref(ibuild)*urefb(ibuild)**2
!            pwact=.5*rhoact*(utotmax(ibuild)**2-utotcl1(ibuild)**2)
!            aref=.3**(1.5)*psref+.19**(1.5)*pwref-.333*sqrt((.3*.19)**1.5*psref*pwref)
!            aact=.3**(1.5)*psact+.19**(1.5)*pwact-.333*sqrt((.3*.19)**1.5*psact*pwact)
!            racact(ibuild)=racref(ibuild)*(aact/aref)**(2./3.)
!            ainf(ibuild)=(racact(ibuild)+racintake(ibuild)*(1-filtration(ibuild)))/3600.
!            binf(ibuild)=(racact(ibuild)+racintake(ibuild)+&
!                         (racsupply(ibuild)-racintake(ibuild))*filtration(ibuild)+&
!                         3600.*deposition_velocity*area2volume(ibuild)+3600.*decay)/3600.
!            abinf(ibuild)=ainf(ibuild)/binf(ibuild)
!            binfex(ibuild)=(racact(ibuild)+racintake(ibuild))/3600.
!            tau(ibuild)=3600./racact(ibuild)
!            if(tau(ibuild).lt.taumin) taumin=tau(ibuild)
! calculate effective areas of building faces for use in averaging
! mdw 3/22/2004 pm set non-zero minimum areas to avoid divide by zero for obscurred faces
         enddo
         allocate(chioa(inumbuild,6),chiosum(inumbuild,noutilf),chiin(inumbuild,noutilf))
         allocate(dosin(inumbuild,noutilf),chiex(inumbuild,noutilf))
         chiosum(:,:)=0.
         chiin(:,:)=0.
         chioa(:,:)=0.
         chiex(:,:)=0.
         dosin(:,:)=0.
         if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
            allocate(dosinSP(inumbuild,noutilf),chiexSP(inumbuild,noutilf))
            allocate(chioaSP(inumbuild,6),chiosumSP(inumbuild,noutilf),chiinSP(inumbuild,noutilf))
            chiosumSP(:,:)=0.
            chiinSP(:,:)=0.
            chioaSP(:,:)=0.
            chiexSP(:,:)=0.
            dosinSP(:,:)=0.
         endif
         return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine readQPbuoy
         use variables
         implicit none
         open(unit=39,file="QP_buoy.inp",status="old")
         read(39,*)ibuoyflag    !=0 neutral, =+1 means positive buoy, =-1 negative
         read(39,*)tair       !ambient temperature degrees C
         read(39,*)tboil        !substance temperature at the source (c)
         read(39,*)rhotref      !gas density at reference temp
         read(39,*)tsourcref    !reference temperature for gas density (C)
         if(ibuoyflag.gt.0)then
            kbpufmax=10
            read(39,*)tfire
            read(39,*)nupplevtemp
            read(39,*)
            allocate(ztlev(nupplevtemp),uptemp(nupplevtemp))
            do itemp=1,nupplevtemp
               read(39,*)ztlev(itemp),uptemp(itemp)
            enddo
            allocate(rhocell(nx,ny,nz-1),airtemp(nz-1),rhobkg(nz-1))
            allocate(zirhol(nx,ny),zirhoh(nx,ny))
            allocate (volpart(kbpufmax),volpartp(kbpufmax),xsum(kbpufmax),ysum(kbpufmax),xsqsum(kbpufmax),ysqsum(kbpufmax))
            allocate (activepart(kbpufmax),volpart0(kbpufmax),voldelrho(kbpufmax),volzrho(kbpufmax),sigyb(kbpufmax),vlsum(kbpufmax))
            allocate (xic(kbpufmax),yic(kbpufmax),xsqic(kbpufmax),ysqic(kbpufmax),varx(kbpufmax),vary(kbpufmax),sigxb(kbpufmax))
            allocate (activepart0(kbpufmax),activepartp(kbpufmax),wbulkzsum(kbpufmax),wbulkza(kbpufmax),conl(kbpufmax))
            do k=2,nz-1
               call airtempcal
               rhobkg(k)=.0288*273./(.0224*airtemp(k))
               do i=1,nx
                  do j=1,ny
                     rhocell(i,j,k)=rhobkg(k)
                  enddo
               enddo
            enddo
            denfire=.0288*273./(.0224*tfire)
            rfire=3.*((source_hemass(1)*rhobkg(2))/denfire)**.333
            source_r(1)=1.0*rfire
            firevol=(4./3.)*pi*rfire**3
            denmax=real(nparticles)*dx*dy*dz/firevol
            !if(source_z(1).le.rfire)source_z(1)=2*z0+rfire !spherical source
            do kb=1,kbpufmax
               xsum(kb)=0.
               ysum(kb)=0.
               xsqsum(kb)=0.
               ysqsum(kb)=0.
               wbulkzsum(kb)=0.
               activepart(kb)=0.
               activepartp(kb)=0.
            enddo
            zsqsum=0.
            zsum=0.
         else
!            rhoboil(isource)=rhotref*(tsourcref+273.)/(tboil+273.)
            rhoair=1.20*(293./(tair+273.))
            allocate(rhocell(nx,ny,nz))
!            allocate(rplusij(nx,ny),rminusij(nx,ny))
            allocate(ubulkpx(100001),vbulkpy(100001),rpm(100000),rl(4),rpl(100000),ustarr(nx,ny),hnewij(nx,ny),dhdtij(nx,ny))
            allocate(ixpl(100000),jypl(100000),ixml(100000),jyml(100000),idxij(nx,ny))
            allocate(hnewi(100001),rnewi(100001),cosphidgi(100001),sinphidgi(100001),&
                  xsceni(100001),ysceni(100001),drdxi(100000),ubari(100001))
            allocate(drdxpi(100001),drdxmi(100001),rnewpi(100001),rnewmi(100001),rplusij(100001),rminusij(100001))
            ubulkpx(:)=0.
            vbulkpy(:)=0.
            hnewij(:,:)=0.
            dhdtij(:,:)=0.
            rplusij(:)=0.
            rminusij(:)=0.
            drdxpi(:)=0.
            drdxmi(:)=0.
            rnewpi(:)=0.
            rnewmi(:)=0.
            idxij(:,:)=-99
            rpm(:)=0.
            rpl(:)=0.
            rhocell(:,:,:)=rhoair
            buoyintz(:,:,:)=0.
            ubuoy(:,:,:)=0.
            vbuoy(:,:,:)=0.
            wbuoy(:,:,:)=0.
         endif
         return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine readQPparticlesize
         use variables
         implicit none
         real cum_fraction
         open(unit=29,file="QP_particlesize.inp",status="old")
         read(29,*)mindiameter
         read(29,*)maxdiameter
         read(29,*)npartsize
         npuff=npartsize
         allocate(dpm_bin(npuff),isub(npuff),pufmasfr(npuff),pufcummas(npuff))
         read(29,*)
         do ippr=1,npartsize
            read(29,*)dpm_bin(ippr),pufmasfr(ippr)
         enddo
         do ipuf=1,npuff
            if(ipuf.gt.1)pufcummas(ipuf)=pufcummas(ipuf-1)+pufmasfr(ipuf)
            if(ipuf.eq.1)pufcummas(ipuf)=pufmasfr(ipuf)
         enddo
         read(29,*)lognormal_flag
         read(29,*)distribution_type_flag
         read(29,*)num_distributions
         allocate(distribution_mean(num_distributions),&
                  distribution_mass_mean(num_distributions),&
                  distribution_std(num_distributions),&
                  distribution_mass_fraction(num_distributions),&
                  distribution_num_fraction(num_distributions),&
                  distribution_cum_fraction(num_distributions))
         cum_fraction=0.
         read(29,*)
         do ippr=1,num_distributions
            read(29,*)distribution_mean(ippr),distribution_std(ippr),&
                      distribution_mass_fraction(ippr),distribution_num_fraction(ippr)
            if(distribution_type_flag .eq. 2)then
               if(distribution_std(ippr) .gt. 1.)then
                  distribution_mass_mean(ippr)=distribution_mean(ippr)*&
                                       exp(3.*(log(distribution_std(ippr))**2.))
               else
                  distribution_mass_mean(ippr)=distribution_mean(ippr)
               endif
               cum_fraction=cum_fraction+distribution_num_fraction(ippr)
            else
               cum_fraction=cum_fraction+distribution_mass_fraction(ippr)
            endif
            distribution_cum_fraction(ippr)=cum_fraction
         enddo
         return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine readQPerad
         use variables
         implicit none
         open(unit=69,file="QP_erad.inp",status="old")
         read(69,*)npartsize
         read(69,*)mindiameter
         read(69,*)maxdiameter
         allocate(dpm_bin(npartsize))
         read(69,*)(dpm_bin(ippr),ippr=1,npartsize)
         npuff=90
         allocate(isub(npuff),pufmasfr(npuff),znorm(npuff),psigx(npuff))
         allocate(psigy(npuff),psigz(npuff),pufcummas(npuff))
         lognormal_flag=0
         do ipuf=1,npuff
            read(69,*)ipufr,imatnum,isub(ipuf),pufmasfr(ipuf),xdum,ydum,znorm(ipuf),psigx(ipuf),psigy(ipuf),psigz(ipuf)
            if(ipuf.gt.1) pufcummas(ipuf)=pufcummas(ipuf-1)+pufmasfr(ipuf)
            if(ipuf.eq.1) pufcummas(ipuf)=pufmasfr(ipuf)
         enddo
         return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine readQPmaterials
         use variables
         implicit none
         read(7,*)
         read(7,*)decay
         decay=decay*0.01/60 !convert decay rate from percent per minute to fraction per second
         read(7,*)dpart
         read(7,*)depvel
         read(7,*)nthreshold
         read(7,*)
         allocate(dosethreshold(nthreshold))
         do ithreshold=1,nthreshold
            read(7,*)dosethreshold(ithreshold)
         enddo
         return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine readQPfileoptions
         use variables
         implicit none
         read(62,*)format_flag
         read(62,*)iwallfieldflag ! options.dat iwallfieldflag => unit number see above
         read(62,*)iturbfieldflag !iturbfieldflag= 1 writes out turbulence, =0 doesn't
         read(62,*)ilctflag !ilctflag=1  writes out lct's; =0 doesn't
         read(62,*)iturbtypeflag !ilctflag=1  writes out lct's; =0 doesn't
         inverseflag=0
!         if(ilctflag .eq. 1)then
!            inverseflag=1
!         else
!            inverseflag=0
!         endif
         return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine readQUsimparams
         use variables
         implicit none
         read(34,*)nx            !nx defined in input file
         read(34,*)ny            !ny defined in input file
         read(34,*)nz            !nz defined in input file
         nx=nx+1
         ny=ny+1
         nz=nz+2
         read(34,*)dx      !x grid resolution in meters
         read(34,*)dy        !y grid resolution in meters
         read(34,*)dz        !z grid resolution in meters
         read(34,*)tstart  !decimal start time
         read(34,*)tinc        !time increments
         read(34,*)ntinc       !number of times
         tupdate=3600.*tinc
         nx1=nx-1
         ny1=ny-1
         nz1=nz-1

! mdw 3/22/2004pm added read statements before and after roofflag read per TMB's instructions
         read(34,*)roofflag      !rooftop recirc flag
         read(34,*)
         read(34,*)
         return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      subroutine readQUbuildings
!         use variables
!         implicit none
!         integer buildnum,buildtype,buildgroup
!         real buildheight,buildwidth,buildlength,buildxfo,buildyfo,buildzfo
!         real buildgamma,vegattenuation
!         read(35,*)xsubsw
!         read(35,*)ysubsw
!         read(35,*)xsubne
!         read(35,*)ysubne
!         read(35,*)
!         read(35,*)inumbuild !number of buildings
!         read(35,*)
!         inumveg=0;
!         do ibuild=1,inumbuild
!            read(35,*)buildnum,buildgroup,buildtype,buildheight,buildwidth,buildlength,&
!                      buildxfo,buildyfo,buildzfo,buildgamma,vegattenuation
!            if(buildtype .eq. 9)inumveg = inumveg + 1
!         enddo
!         return
!      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine turbinit
! turbinit computes the turbulence parameters for the current wind field
         use variables   !make data from module "variables" visible
         implicit none
         integer i_b,iupc,jupc,kendv,icbp3,icbm3,jcbp3,jcbm3,jp1,jp2,ji,jm1,jm2
         integer ip1,ip2,im1,im2,icd,jcd,icdl,jcdl,icu,jcu,icul,jcul,ktop,kmid
         integer ktp,istinf,istf,isf,isfu,isini,is,icel,jcel,ist,icelt,jcelt,isin
         integer iceln,jceln,npentx,npenty,ipent1,ipent2,jpent1,jpent2,check1,check2,check3,check
         integer nx2,ny2,nz2,klim,kdif,ilim,ii,jlim,jj,incourt
         real del_b,cosphit,bs,bl,cosfac,sinphit,rscale,zclim,ycbp,xcbp,ycbm,xcbm
         real xcell,ycell,xcelt,ycelt,delut,xceln,yceln,xpent1,ypent1,eps
         real dsigwdx,dsigwdy,dsigudx,dsigudy,dsigvdx,dsigvdy,dupwpdx,dupwpdy
         real dsigvdz,dsigudz,dupwpdz,dupwpdn,elzv
         real x0,y0,x1,y1,x2,y2,x3,y3,x4,y4,xc,yc,rc,tc
         allocate(ufsqgi(nx,ny,nz),vfsqgi(nx,ny,nz),wfsqgi(nx,ny,nz))
         allocate(ufvfgi(nx,ny,nz),ufwfgi(nx,ny,nz),vfwfgi(nx,ny,nz))
         allocate(ufwfi(nx,ny,nz),ufvfi(nx,ny,nz),vfwfi(nx,ny,nz))
         allocate(uktop(nx,ny),vktop(nx,ny),wktop(nx,ny),utotktp(nx,ny))
         ustarg(:,:,:)=0.
         ufsqgi(:,:,:)=0.
         vfsqgi(:,:,:)=0.
         wfsqgi(:,:,:)=0.
         ufvfgi(:,:,:)=0.
         ufwfgi(:,:,:)=0.
         vfwfgi(:,:,:)=0.
         ufwfi(:,:,:)=0.
         ufvfi(:,:,:)=0.
         vfwfi(:,:,:)=0.
         ht_avg=0.
         if(inumbuild .gt. 0 .and. inumbuild .ne. inumveg)then
            do i_b=1,inumbuild
               if(bldtype(i_b) .eq. 9)then
                  cycle
               endif
               ht_avg=Ht(i_b)+zfo(i_b)+ht_avg
            enddo
            ht_avg=ht_avg/real(inumbuild-inumveg)
            k=nint(ht_avg/dz)+1
         else
            k=2
         endif
         i=1
         u_left=0
         do j=1,ny-1
            u_left=sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2) +u_left
         enddo
         u_left=u_left/(ny-1)
         j=ny-1
         u_top=0
         do i=1,nx-1
            u_top=sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2) +u_top
         enddo
         u_top=u_top/(nx-1)
         i=nx-1
         u_right=0
         do j=1,ny-1
            u_right=sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2) +u_right
         enddo
         u_right=u_right/(ny-1)
         j=1
         u_bottom=0
         do i=1,nx-1
            u_bottom=sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2) +u_bottom
         enddo
         u_bottom=u_bottom/(nx-1)
         u_b=(u_left+ u_top+ u_right+ u_bottom)/4
         nu_b=1.5e-5
         del_b=(0.328* (nu_b/u_b)**(.2) ) * (ht_avg)**(.8)
lp020:   do k=1,nz-1
lp019:      do j=1,ny-1
lp018:         do i=1,nx-1
                  if(k.eq.1)icellflag(i,j,k)=0.
                  if(k.ne.1)then
                     if(icellflag(i,j,k-1) .eq. 0 .and. icellflag(i,j,k) .ne. 0 .or. k .eq. 2)then
                        utotl=0.
                        utotu=sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
!                       dutotdzi(i,j,k)=(utotu-utotl)/(.5*dz)
! MDW 7-01-2005 changed the way vertical gradients are calculated to avoid inaccuracies
! in the representation of the gradients of a log-law term
                        if(rcl.gt.0)then
                           phim=1.+4.7*rcl*.5*dz
                           psim=-4.7*rcl*.5*dz
                        else
                           phim=(1.-15.*rcl*.5*dz)**(-.25)
                           psim=2.*alog((1.+1./phim)/2.)+alog((1.+1./phim**2)/2.)-2.*atan(1./phim)+pi/2.
                        endif
                        ustar=kkar*utotu/(log(.5*dz/z0)-psim)
                        elz(i,j,k)=kkar*.5*dz
                        ustarz(i,j,k)=ustar
                        dutotdzi(i,j,k)=ustar*phim/(kkar*.5*dz)
                        sigwi(i,j,k-1)=0.
                        sigvi(i,j,k-1)=0.
                        ustarij(i,j,k-1)=0.
                        ustarz(i,j,k-1)=0.
                        hgt(i,j)=zi(k-1)+.5*dz
                     else
                        if(k.eq.nz-1)then ! find gradient using a non-CDD approach
                           utotl=sqrt(u(i,j,k-1)**2+v(i,j,k-1)**2+w(i,j,k-1)**2)
                           utotu=sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
                           dutotdzi(i,j,nz-1)=dutotdzi(i,j,k-2)
!                           elz(i,j,k)=kkar*(dz*real(k-1)-.5*dz-hgt(i,j))
                           elz(i,j,k)=kkar*eleff(i,j,k)
                           print*,i,j,k,eleff(i,j,k)
                        else ! find gradient using a CDD approach
                           utotl=sqrt(u(i,j,k-1)**2+v(i,j,k-1)**2+w(i,j,k-1)**2)
                           utotu=sqrt(u(i,j,k+1)**2+v(i,j,k+1)**2+w(i,j,k+1)**2)
!                          dutotdzi(i,j,k)=(utotu-utotl)/(dz)
! mdw 7-08-2005 changed the way vertical gradients are calculated to better represent
! log-law behavior
!                           hgt(i,j)=zi(k-1)+.5*dz
                           klow=int((hgt(i,j)+.5*dz+dz)/dz)+1
                           if(rcl.gt.0)then
                              phim=1.+4.7*rcl*eleff(i,j,k-1)
                              psim=-4.7*rcl*eleff(i,j,k-1)
                           else
                              phim=(1.-15.*rcl*eleff(i,j,k-1))**(-.25)
                              psim=2.*alog((1.+1./phim)/2.)+alog((1.+1./phim**2)/2.)-2.*atan(1./phim)+pi/2.
                           endif
                           dutotl=utotl-ustarz(i,j,klow)*(log(zi(k-1)/z0)-psim)/kkar
                           if(rcl.gt.0)then
                              phim=1.+4.7*rcl*eleff(i,j,k+1)
                              psim=-4.7*rcl*eleff(i,j,k+1)
                           else
                              phim=(1.-15.*rcl*eleff(i,j,k+1))**(-.25)
                              psim=2.*alog((1.+1./phim)/2.)+alog((1.+1./phim**2)/2.)-2.*atan(1./phim)+pi/2.
                           endif
                           dutotu=utotu-ustarz(i,j,klow)*(log(zi(k+1)/z0)-psim)/kkar
                           dutotdzi(i,j,k)=(dutotu-dutotl)/(2.*dz)+ustarz(i,j,klow)*psim/(kkar*zi(k))
                           elz(i,j,k)=kkar*eleff(i,j,k)
                           if(icellflag(i,j,k+1) .ne. 0 .and. icellflag(i,j,k) .ne. 0&
                                 .and. icellflag(i,j,k-1) .ne. 0)then
! mdw 7-01-2005 centered around k instead of k-1 and ajusted for log-law behavior
                              utot=sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
                              utotl=sqrt(u(i,j,k-1)**2+v(i,j,k-1)**2+w(i,j,k-1)**2)
                              utotu=sqrt(u(i,j,k+1)**2+v(i,j,k+1)**2+w(i,j,k+1)**2)
                              if(rcl.gt.0)then
                                 phim=1.+4.7*rcl*eleff(i,j,k-1)
                                 psim=-4.7*rcl*eleff(i,j,k-1)
                              else
                                 phim=(1.-15.*rcl*eleff(i,j,k-1))**(-.25)
                                 psim=2.*alog((1.+1./phim)/2.)+alog((1.+1./phim**2)/2.)-2.*atan(1./phim)+pi/2.
                              endif
                              dutotl=utotl-ustarz(i,j,klow)*(log(zi(k-1)/z0)-psim)/kkar
                              if(rcl.gt.0)then
                                 phim=1.+4.7*rcl*eleff(i,j,k+1)
                                 psim=-4.7*rcl*eleff(i,j,k+1)
                              else
                                 phim=(1.-15.*rcl*eleff(i,j,k+1))**(-.25)
                                 psim=2.*alog((1.+1./phim)/2.)+alog((1.+1./phim**2)/2.)-2.*atan(1./phim)+pi/2.
                              endif
                              dutotu=utotu-ustarz(i,j,klow)*(log(zi(k+1)/z0)-psim)/kkar
! mdw 3-08-2004 begin changes for highest gradient rather than centered diff gradient
                              if(rcl.gt.0)then
                                 phim=1.+4.7*rcl*eleff(i,j,k)
                                 psim=-4.7*rcl*eleff(i,j,k)
                              else
                                 phim=(1.-15.*rcl*eleff(i,j,k))**(-.25)
                                 psim=2.*alog((1.+1./phim)/2.)+alog((1.+1./phim**2)/2.)-2.*atan(1./phim)+pi/2.
                              endif
                              dutot=utot-ustarz(i,j,klow)*(log(zi(k)/z0)-psim)/kkar
                              dutotdzc=.5*(dutotu-dutotl)/dz
                              dutotdzp=(dutotu-dutot)/dz
                              dutotdzm=(dutot-dutotl)/dz
                              dutotdza=0.5*(abs(dutotdzp+ustarz(i,j,klow)*phim/(kkar*zi(k)))&
                                      +abs(dutotdzm+ustarz(i,j,klow)*phim/(kkar*zi(k))))
                              if(abs(dutotdzp+ustarz(i,j,klow)*phim/(kkar*zi(k))).gt. &
                                 abs(dutotdzm+ustarz(i,j,klow)*phim/(kkar*zi(k))))then
                                 dutotdzi(i,j,k)=dutotdzp+ustarz(i,j,klow)*phim/(kkar*zi(k))
                              else
                                 dutotdzi(i,j,k)=dutotdzm+ustarz(i,j,klow)*phim/(kkar*zi(k))
                              endif
! use centered differences away from the boundaries
                           endif
                        endif
                     endif
                  endif
               enddo        lp018
            enddo           lp019
         enddo              lp020
         phi=270.-theta
         print*,'theta:',theta
         phi=phi*pi/180.
         cosphi=cos(phi)
         if(cosphi.ge.0)then
            iupc=1
         else
            iupc=nx-1
         endif
         sinphi=sin(phi)
         if(sinphi.ge.0)then
            jupc=1
         else
            jupc=ny-1
         endif
         phit=phi+0.5*pi
         cosphit=cos(phit)
         sinphit=sin(phit)
lp026:   do i=1,inumbuild
            if(bldtype(i) .eq. 9)cycle
!           phi=phib(i)

! mdw 4-16-2004 added proper treatment of zfo
            ktop=nint((Ht(i)+zfo(i)+dz)/dz)
            kmid=nint((.5*Ht(i)+zfo(i)+dz)/dz)
            if(bldtype(i).eq.3)then
               xcb(i)=xfo(i)
            else
               xcb(i)=xfo(i)+.5*Lt(i)
            endif
            ycb(i)=yfo(i)
!           icb(i)=int(xcb(i)/dx)+1      !!!!---- org
!           jcb(i)=int(ycb(i)/dy)+1      !!!!---- org
            icb(i)=nint(xcb(i)/dx)
            jcb(i)=nint(ycb(i)/dy)
!mdw 6-05-2005 put in procedure to calculate phi & phit
            if(roofflag.eq.2)then
               Bs=Ht(i)
               BL=Wti(i)

               if(Wti(i).lt.Ht(i))then
                  Bs=Wti(i)
                  BL=Ht(i)
               endif
               Rscale = ((Bs**(2./3.))*(BL**(1./3.)))
               zclim=max(.22*Rscale,.11*Wti(i),.11*Lt(i))
               kendv=1+(int(zclim+Ht(i)+zfo(i)+dz)/dz)
            else
               kendv=1+(int(Ht(i)+zfo(i)+dz)/dz)
            endif
            phib(i)=atan2(v(icb(i),jcb(i),kendv),u(icb(i),jcb(i),kendv))
            print*,'phib:',phib(i)
            phi=phib(i)
            cosphi=cos(phi)
            if(cosphi.ge.0)then
               iupc=1
            else
               iupc=nx-1
            endif
            sinphi=sin(phi)
            if(sinphi.ge.0)then
               jupc=1
            else
               jupc=ny-1
            endif
            phit=phi+0.5*pi
            cosphit=cos(phit)
            sinphit=sin(phit)
! ycbp3, and xcbp3 give points 1.5 units outside
! of the bldg boundaries to compute reference utot
            if(abs(sinphit).ge.abs(cosphit))then
               ycbp3=ycb(i)+(.5*weff(i)+.33*weff(i))*sinphit ! Get reference values for x,y for non-local mixing
               xcbp3=xcb(i)+(.5*weff(i)+.33*weff(i))*cosphit ! 1/3 bldg width outside of building is the boundary for the non-local mixing
               ycbm3=ycb(i)-(.5*weff(i)+.33*weff(i))*sinphit
               xcbm3=xcb(i)-(.5*weff(i)+.33*weff(i))*cosphit
               icbp3=nint(xcbp3/dx)
               icbm3=nint(xcbm3/dx)
               jcbp3=nint(ycbp3/dy)
               jcbm3=nint(ycbm3/dy)
               jcbp3=min(jcbp3,ny-1)
               jcbm3=min(jcbm3,ny-1)
               icbp3=min(icbp3,nx-1)
               icbm3=min(icbm3,nx-1)
               jcbp3=max(1,jcbp3)
               jcbm3=max(1,jcbm3)
               icbp3=max(1,icbp3)
               icbm3=max(1,icbm3)
! searching in the plus y direction for building free flow
               if(icellflag(icbp3,jcbp3,kmid) .eq. 0)then
                  if(sinphit.gt.0.)then
                     jp1=jcbp3
                     jp2=ny-1
                     isign=1
                  else
                     jp1=jcbp3
                     jp2=1
                     isign=-1
                  endif
                  do ji=jp1,jp2,isign
                     jcbp3=jcbp3+isign
                     jcbp3=min(ny-1,jcbp3)
                     dycbp3=dy*(jcbp3-1)-ycbp3
                     ycbp3=dy*(jcbp3-1)
                     xcbp3=xcbp3+cosphit*dycbp3/sinphit
                     icbp3=int(xcbp3/dx)+1
                     icbp3=min(nx-1,icbp3)
!mdw 34/01/2004 forced indices to be within domain
                     if(icellflag(icbp3,jcbp3,kmid) .ne. 0) exit
                  enddo
               endif
! searching in the minus y direction for building free flow
               if(icellflag(icbm3,jcbm3,kmid) .eq. 0)then
                  if(sinphit.gt.0.)then
                     jm2=1
                     jm1=jcbm3
                     isign=1
                  else
                     jm2=ny-1
                     jm1=jcbm3
                     isign=-1
                  endif
                  do ji=jm1,jm2,-isign
                     jcbm3=jcbm3-isign
                     dycbm3=dy*(jcbm3-1)-ycbm3
                     ycbm3=dy*(jcbm3-1)
                     xcbm3=xcbm3+cosphit*dycbm3/sinphit
                     icbm3=nint(xcbm3/dx)
                     jcbp3=min(jcbp3,ny-1)
                     jcbm3=min(jcbm3,ny-1)
                     icbp3=min(icbp3,nx-1)
                     icbm3=min(icbm3,nx-1)
                     jcbp3=max(1,jcbp3)
                     jcbm3=max(1,jcbm3)
                     icbp3=max(1,icbp3)
                     icbm3=max(1,icbm3)
                     if(icellflag(icbm3,jcbm3,kmid) .ne. 0) exit
                  enddo
               endif
               ycbp=ycb(i)+(.5*Leff(i))*sinphi
               xcbp=xcb(i)+(.5*Leff(i))*cosphi
               ycbm=ycb(i)-(.5*Leff(i))*sinphi
               xcbm=xcb(i)-(.5*Leff(i))*cosphi
               if(cosphi.ge.0.)then
! Note the current upstream and downstream limits for the wake non-local mixing
! are 3*lr in the downstream direction and lfx upstream in the x direction
! and lfy upstream in the y direction
                  xcd=xcb(i)+(.5*leff(i)+.1*dx)*cosphi ! get the first point on the center line outside of the building (downstream)
                  ycd=ycb(i)+(.5*leff(i)+.1*dx)*sinphi !
!mdw 7-10-2006 made changes to xcd, ycd,xcu, & ycu - formerly used .5 dx
                  if(bldtype(i).eq.3)then
                     xcu=xcb(i)-(.4*leff(i)+dx)*cosphi ! (upstream)
                     ycu=ycb(i)-(.4*leff(i)+dx)*sinphi !
                  else
                     xcu=xcb(i)-(.5*leff(i)+0.1*dx)*cosphi ! (upstream)
                     ycu=ycb(i)-(.5*leff(i)+0.1*dx)*sinphi !
                  endif
!                 xcdl=xcd+lr(i)*cosphi
!                 ycdl=ycd+lr(i)*sinphi
!mdw 7-05-2006 made changes to xcul & ycul - formerly used .5 dx
                  xcul=xcu-(lfr(i)+dx)*cosphi ! get upper limit of the eddie
                  ycul=ycu-(lfr(i)+dy)*sinphi
                  xcul=max(xcul,0.)
                  xcul=min(xcul,dx*real(nx-1))
                  ycul=max(ycul,0.)
                  ycul=min(ycul,dy*real(ny-1))
                  cosfac=1.
               else
!mdw 7-10-2006 made changes to xcd, ycd,xcu, & ycu - formerly used .5 dx
                  xcd=xcb(i)+(.5*leff(i)+.1*dx)*cosphi
                  ycd=ycb(i)+(.5*leff(i)+.1*dx)*sinphi
                  if(bldtype(i).eq.3)then
                     xcu=xcb(i)-(.4*leff(i)+dx)*cosphi ! (upstream)
                     ycu=ycb(i)-(.4*leff(i)+dx)*sinphi !
                  else
                     xcu=xcb(i)-(.5*leff(i)+0.1*dx)*cosphi ! (upstream)
                     ycu=ycb(i)-(.5*leff(i)+0.1*dx)*sinphi !
                  endif
!                 xcdl=xcd-lr(i)*cosphi
!                 ycdl=ycd-lr(i)*sinphi
!mdw 7-05-2006 made changes to xcul & ycul - formerly used .5 dx
                  xcul=xcu-(lfr(i)+dx)*cosphi ! get upstream limit on the front cavity
                  ycul=ycu-(lfr(i)+dy)*sinphi !
                  xcul=max(xcul,0.)
                  xcul=min(xcul,dx*real(nx-1))
                  ycul=max(ycul,0.)
                  ycul=min(ycul,dy*real(ny-1))
                  cosfac=-1.
               endif
            else ! if you are more aligned with y than x
! MAN 9/15/2005 use weff and leff appropriately
               ycbp3=ycb(i)+(.5*weff(i)+.33*weff(i))*sinphit ! get the effective length of the building
               xcbp3=xcb(i)+(.5*weff(i)+.33*weff(i))*cosphit
               ycbm3=ycb(i)-(.5*weff(i)+.33*weff(i))*sinphit
               xcbm3=xcb(i)-(.5*weff(i)+.33*weff(i))*cosphit
! end MAN 9/15/2005
               icbp3=nint(xcbp3/dx)
               icbm3=nint(xcbm3/dx)
               jcbp3=nint(ycbp3/dy)
               jcbm3=nint(ycbm3/dy)
               jcbp3=min(jcbp3,ny-1)
               jcbm3=min(jcbm3,ny-1)
               icbp3=min(icbp3,nx-1)
               icbm3=min(icbm3,nx-1)
               jcbp3=max(1,jcbp3)
               jcbm3=max(1,jcbm3)
               icbp3=max(1,icbp3)
               icbm3=max(1,icbm3)
! make sure you are outside of the building !
               if(icellflag(icbp3,jcbp3,kmid) .eq. 0)then
                  if(cosphit.gt.0)then
                     ip1=icbp3
                     ip2=ny-1
                     isign=1
                  else
                     ip1=icbp3
                     ip2=1
                     isign=-1
                  endif
                  ! decide which is closest building/floor
                  do ip=ip1,ip2,isign
                     icbp3=icbp3+isign
                     dxcbp3=dx*real(icbp3-1)-xcbp3
                     xcbp3=dx*real((icbp3-1))
                     ycbp3=ycbp3+dxcbp3*sinphit/cosphit
                     jcbp3=nint(ycbp3/dy)
                     jcbp3=min(jcbp3,ny-1)
                     jcbm3=min(jcbm3,ny-1)
                     icbp3=min(icbp3,nx-1)
                     icbm3=min(icbm3,nx-1)
                     jcbp3=max(1,jcbp3)
                     jcbm3=max(1,jcbm3)
                     icbp3=max(1,icbp3)
                     icbm3=max(1,icbm3)
                     if(icellflag(icbp3,jcbp3,kmid) .ne. 0) exit
                  enddo
               endif
               if(icellflag(icbm3,jcbm3,kmid) .eq. 0)then
                  if(cosphit.gt.0.)then
                     im1=icbm3
                     im2=1
                     isign=1
                  else
                     im1=icbm3
                     im2=nx-icbm3+1
                     isign=-1
                  endif
                  do im=im1,im2,-isign
                     icbm3=icbm3-isign
                     dxcbm3=dx*real((icbm3-1))-xcbm3
                     xcbm3=dx*real((icbm3-1))
                     jcbm3=jcbm3+dxcbm3*sinphit/cosphit
                     jcbp3=min(jcbp3,ny-1)
                     jcbm3=min(jcbm3,ny-1)
                     icbp3=min(icbp3,nx-1)
                     icbm3=min(icbm3,nx-1)
                     jcbp3=max(1,jcbp3)
                     jcbm3=max(1,jcbm3)
                     icbp3=max(1,icbp3)
                     icbm3=max(1,icbm3)
                     if(icellflag(icbm3,jcbm3,kmid) .ne. 0) exit
                  enddo
               endif
               ycbp=ycb(i)+(.5*leff(i))*sinphit !  get back of the building
               xcbp=xcb(i)+(.5*leff(i))*cosphit !
               ycbm=ycb(i)-(.5*leff(i))*sinphit !  get front of the building
               xcbm=xcb(i)-(.5*leff(i))*cosphit !
               if(sinphi.ge.0.)then
! Note the current upstream and downstream limits for the wake non-local mixing
! are 3*lr in the downstream direction and lfx upstream in the x direction
! and lfy upstream in the y direction
! MAN 9/15/2005 use weff and leff appropriately
!mdw 7-05-2006 made changes to xcu,ycu, xcd & ycd - formerly used .5 dy or .5 dx
                  xcd=xcb(i)+(.5*leff(i)+dy)*cosphi ! get the first point on the center line outside of the building (downstream)
                  ycd=ycb(i)+(.5*leff(i)+dy)*sinphi !
                  if(bldtype(i).eq.3)then
                     xcu=xcb(i)-(.4*leff(i)+dx)*cosphi ! (upstream)
                     ycu=ycb(i)-(.4*leff(i)+dx)*sinphi !
                  else
                     xcu=xcb(i)-(.5*leff(i)+0.1*dx)*cosphi ! (upstream)
                     ycu=ycb(i)-(.5*leff(i)+0.1*dx)*sinphi !
                  endif
! end MAN 9/15/2005
!                 xcdl=xcd+lr(i)*cosphi
!                 ycdl=ycd+lr(i)*sinphi
! mdw 7-05-2006 eliminated .5 dx  or .5 dy in favor of dx & dy
                  xcul=xcu-(lfr(i)+dx)*cosphi ! get upper limit of the eddie
                  ycul=ycu-(lfr(i)+dy)*sinphi
                  xcul=max(xcul,0.)
                  xcul=min(xcul,dx*real(nx-1))
                  ycul=max(ycul,0.)
                  ycul=min(ycul,dy*real(ny-1))
                  cosfac=1.
               else
! MAN 9/15/2005 use weff and leff appropriately
                  xcd=xcb(i)+(.5*leff(i)+dy)*cosphi
                  ycd=ycb(i)+(.5*leff(i)+dy)*sinphi
                  if(bldtype(i).eq.3)then
                     xcu=xcb(i)-(.4*leff(i)+dx)*cosphi ! (upstream)
                     ycu=ycb(i)-(.4*leff(i)+dx)*sinphi !
                  else
                     xcu=xcb(i)-(.5*leff(i)+dx)*cosphi ! (upstream)
                     ycu=ycb(i)-(.5*leff(i)+dx)*sinphi !
                  endif
! end MAN 9/15/2005
!                 xcdl=xcd-lr(i)*cosphi
!                 ycdl=ycd-lr(i)*sinphi
                  xcul=xcu+(lfr(i)+dx)*cosphi ! get upstream limit on the front cavity
                  ycul=ycu+(lfr(i)+dy)*sinphi !
                  xcul=max(xcul,0.)
                  xcul=min(xcul,dx*real(nx-1))
                  ycul=max(ycul,0.)
                  ycul=min(ycul,dy*real(ny-1))
                  cosfac=-1.
               endif
            endif
!mdw 7-05-2006 change form to ixxx or jxxx =nint()+1
            icd=nint(xcd/dx)+1 ! get indicies for the downstream center line to back of the building
            jcd=nint(ycd/dy)+1 !
!mdw 4-16-2004 added correction for ktop+3 > nz-1
            ktp=min(ktop,nz-1)
lp025:      do k=ktp,2,-1 ! Account for wake difference in the cavity
               zk=dz*(k-2)+.5*dz
               zbrac=(1.-zi(k)/h)**1.5
!mdw 4-16-2004 added correction for ktop+3 > nz-1
               if(ktop+3.le.nz-1)then
                  zcorf(k)=sqrt(u(iupc,jupc,k)**2+v(iupc,jupc,k)**2+w(iupc,jupc,k)**2) &
                           /sqrt(u(iupc,jupc,ktop+3)**2+v(iupc,jupc,ktop+3)**2+w(iupc,jupc,ktop+3)**2)
               else
                  zcorf(k)=sqrt(u(iupc,jupc,k)**2+v(iupc,jupc,k)**2+w(iupc,jupc,k)**2) &
                           /sqrt(u(iupc,jupc,nz-1)**2+v(iupc,jupc,nz-1)**2+w(iupc,jupc,nz-1)**2)
               endif
! mdw 4-16-2004 added proper treatment of zfo
               if(zk.lt.ht(i)+zfo(i))then
                  zkfac=sqrt(1.-(zk/(ht(i)+zfo(i)))**2)
               else
                  if(k.eq.ktp)then
                     zkfac=1.
                  else
                     zkfac=0.
                  endif
               endif
! mdw 7-05-2006 changed from .5 dx or .5 dy to dx & dy to be consistent with nint
               xcdl=xcd+(3.*lr(i)+dx)*zkfac*cosphi ! calculate the x,y limit of the wake as a function of height
               ycdl=ycd+(3.*lr(i)+dy)*zkfac*sinphi !
               xcdl=min(xcdl,dx*real(nx-1))
               ycdl=min(ycdl,dy*real(nx-1))
               xcdl=max(xcdl,0.)
               ycdl=max(ycdl,0.)
               icdl=nint(xcdl/dx)+1 ! Calculate the indicies for i,j according to xcdl,ycdl
               jcdl=nint(ycdl/dy)+1 !
               icu=nint(xcu/dx)+1   ! indicies for the upstream cavity (building)
               jcu=nint(ycu/dy)+1   !
               icul=nint(xcul/dx)+1 ! (furthest upstream)
               jcul=nint(ycul/dy)+1 !!!
!mdw 4-16-2004 added correction for ktop+3 > nz-1
               if(ktop+3.le.nz-1)then
! calculating the reference wind un-disturbed by the building
                  urefz=sqrt(u(icb(i),jcb(i),ktop+3)**2+v(icb(i),jcb(i),ktop+3)**2&
                        +w(icb(i),jcb(i),ktop+3)**2)
               else
                  urefz=sqrt(ualoft**2+valoft**2)
               endif
               ds=0.7*min(dx,dy) ! pick a step that is small enough to not skip grid cells
               sdown=sqrt((xcdl-xcd)**2+(ycdl-ycd)**2)+2.*ds ! calculate the limits for the distance measured along the centerline (rear)
               sup=sqrt((xcul-xcu)**2+(ycul-ycu)**2)+2.*ds   ! same for the front eddy
               stin=.5*leff(i)
               istinf=nint(stin/ds)+1
!mdw 7-11-2006 changed istinf to allow replacement to center of bldg
!mdw 5-14-2004 corrected expression for st; older versions gave errors for wide blds
               st=sqrt((xcbp3-xcb(i))**2+(ycbp3-ycb(i))**2)+1.*ds ! total distance to point
               istf=nint((st+.333*leff(i))/ds)+1   ! (transverse direction)
!mdw 6-9-2004 extended the transverse integration to st+.333*leff
               isf=nint(sdown/ds)+1 ! setup limits of calculations (for do loops) (along cneterline down)
               isfu=nint(sup/ds)+1  ! (along centerline up)
               if(lfr(i) .lt. 0.)isfu=0
!               ktp=min(ktop,nz-1)
!mdw 4-16-2004 added correction for ktop+3 > nz-1
!
! Select the largest reference wind of the plus or minus side of the building
!
               utotp=sqrt(u(icbp3,jcbp3,k)**2+v(icbp3,jcbp3,k)**2+w(icbp3,jcbp3,k)**2)
               utotm=sqrt(u(icbm3,jcbm3,k)**2+v(icbm3,jcbm3,k)**2+w(icbm3,jcbm3,k)**2)
               if(utotp.ge.utotm)then
                  uref(i,k)=utotp+.000001
!                 urefu(i,k)=u(icbp3,jcbp3,k) !!---- org
!                 urefv(i,k)=v(icbp3,jcbp3,k) !!---- org
                  urefu(i,k)=uref(i,k)*cos(phib(i))
                  urefv(i,k)=uref(i,k)*sin(phib(i))
                  urefw(i,k)=w(icbp3,jcbp3,k)
               else
                  uref(i,k)=utotm+.000001
!                 urefu(i,k)=u(icbm3,jcbm3,k)   !!---- org
!                 urefv(i,k)=v(icbm3,jcbm3,k)   !!---- org
                  urefu(i,k)=uref(i,k)*cos(phib(i))
                  urefv(i,k)=uref(i,k)*sin(phib(i))
                  urefw(i,k)=w(icbm3,jcbm3,k)
               endif
!!!!!!!
!              uref(i,k)=.5*utotp+.5*utotm+.000001
!              urefu(i,k)=.5*u(icbp3,jcbp3,k)+.5*u(icbm3,jcbm3,k)
!              urefv(i,k)=.5*v(icbp3,jcbp3,k)+.5*v(icbm3,jcbm3,k)
               cosu=(urefu(i,k)+.000001)/uref(i,k)
               sinv=urefv(i,k)/uref(i,k)
! downstream wake  along axis do loop for delta u
               isini=1
lp022:         do is=1,isf
                  xcell=xcd+ds*real(is-1)*cosphi
                  ycell=ycd+ds*real(is-1)*sinphi
                  icel=nint(xcell/dx)+1
                  jcel=nint(ycell/dy)+1
                  icel=min(nx-1,icel)
                  icel=max(2,icel)
                  jcel=min(ny-1,jcel)
                  jcel=max(2,jcel)
                  if(icellflag(icel,jcel,k) .eq. 0 .and. is .eq. 1)then
                     isini=2
                  endif
                  utot=sqrt(u(icel,jcel,k)**2+v(icel,jcel,k)**2+w(icel,jcel,k)**2)
!mdw 4-16-2004 added correction for ktop+3 > nz-1
                  if(k.eq.ktp)then
                     if(ktop.le.nz-1)then
                        utotktp(icel,jcel)=utot
                        uktop(icel,jcel)=u(icel,jcel,ktop)
                        vktop(icel,jcel)=v(icel,jcel,ktop)
                        wktop(icel,jcel)=w(icel,jcel,ktop)
                     else
                        utotktp=sqrt(ualoft**2+valoft**2)
                        uktop(icel,jcel)=ualoft
                        vktop(icel,jcel)=valoft
                        wktop(icel,jcel)=0.
                     endif
                  endif
! this sets reference for vertical transfer
                  utot=utot+.000001
                  cosl=u(icel,jcel,k)/utot
                  sinl=v(icel,jcel,k)/utot
                  if(icellflag(icel,jcel,k) .gt. 0)then
                     delutz=sqrt((u(icel,jcel,k)-zcorf(k)*uktop(icel,jcel))**2+(v(icel,jcel,k)&
                            -zcorf(k)*vktop(icel,jcel))**2+(w(icel,jcel,k)-zcorf(k)*wktop(icel,jcel))**2)
!                 deluc(i,k)=uref(i,k)-utot*cosu-utot*sinv
                     deluc(i,k)=sqrt((urefu(i,k)-u(icel,jcel,k))**2+(urefv(i,k)-v(icel,jcel,k))**2&
                                +(urefw(i,k)-w(icel,jcel,k))**2)
!mdw 4-16-2004 added correction for ktop+3 > nz-1
                     if(k.ne.ktp)then
                     ! Selects the largest gradient (vert or horiz transfer)
! mdw 4-16-2004 added proper treatment of zfo
                        if((2.*deluc(i,k)/weff(i)).lt.(utotktp(icel,jcel)/(Ht(i)+zfo(i))) &
                                .and.delutz.gt..2*zcorf(k)*utotktp(icel,jcel))then ! vertical dominates
                           ustargz(icel,jcel,k)=max(knlc*utotktp(icel,jcel),ustargz(icel,jcel,k))
                           if(abs(ustargz(icel,jcel,k)-knlc*utotktp(icel,jcel)).lt.1.e-05*ustargz(icel,jcel,k))then ! This value dominates over prev. buildings.
                              elzg(icel,jcel,k)=Ht(i)+zfo(i)
                              upvpg=0.
                              upwpg=-ctau13*zbrac*ustargz(icel,jcel,k)**2
                              upsqg=cusq*zbrac*ustargz(icel,jcel,k)**2
                              vpsqg=cvsq*zbrac*ustargz(icel,jcel,k)**2
                              vpwpg=0.
                              wpsqg=cwsq*zbrac*ustargz(icel,jcel,k)**2
                              ustarg(icel,jcel,k)=ustargz(icel,jcel,k)
                              call rotate2d(icel,jcel,k) ! calculates the quantites in the original coord sys.
                              if(iturbtypeflag.eq.1)nonlocal_option(icel,jcel,k)=4
                           endif
                        else
! We use the vertical gradient as dominant if it is sharper than the horizontal
                           duy=-deluc(i,k)*sinl*cosu+deluc(i,k)*sinv*cosl
! we now have the delta u between the outside of the bldg and the center of the wake
! mdw 6-10-2004 removed
                           if(deluc(i,k).gt..2*uref(i,k))then
                              ustarg(icel,jcel,k)=max(ustarg(icel,jcel,k),knlc*deluc(i,k))
                              if(abs(ustarg(icel,jcel,k)-knlc*deluc(i,k)).lt.1.e-05*ustarg(icel,jcel,k))then ! if the horiz is dominant calculate sigmas
                                 upvpg=0.
! on axis u prime v prime is zero
                                 upwpg=0.
! for eddy transport in uv we don't consider uw
                                 upsqg=cusq*zbrac*ustarg(icel,jcel,k)**2
                                 wpsqg=cvsq*zbrac*ustarg(icel,jcel,k)**2
                                 vpwpg=0.
                                 elzg(icel,jcel,k)=0.5*weff(i)
                                 vpsqg=cwsq*zbrac*ustarg(icel,jcel,k)**2
                                 call rotate2d(icel,jcel,k)
                                 if(iturbtypeflag.eq.1)nonlocal_option(icel,jcel,k)=5
                              endif
                           endif
                        endif
                     endif
                  else
                     deluc(i,k)=0.
                     delutz=0.
                  endif
! transverse do loop in downstream wake
lp021:            do ist=2,istf
! first direction in the transverse of the wake
                     xcelt=xcell+ds*real(ist-1)*cosphit
                     ycelt=ycell+ds*real(ist-1)*sinphit
                     icelt=nint(xcelt/dx)+1
                     jcelt=nint(ycelt/dy)+1
                     if(abs(xcelt-xcell).lt..5*ds)icelt=icel
                     if(abs(ycelt-ycell).lt..5*ds)jcelt=jcel
                     icelt=min(nx-1,icelt)
                     icelt=max(2,icelt)
                     jcelt=min(ny-1,jcelt)
                     jcelt=max(2,jcelt)
                     if(icellflag(icelt,jcelt,k) .gt. 0)then
                        utott=sqrt(u(icelt,jcelt,k)**2+v(icelt,jcelt,k)**2+w(icelt,jcelt,k)**2)
                        utott=utott+.000001
!mdw 4-16-2004 added correction for ktop+3 > nz-1
                        if(k.eq.ktp)then
                           if(ktop.le.nz-1)then
                              utotktp(icelt,jcelt)=utott
                              uktop(icelt,jcelt)=u(icelt,jcelt,ktop)
                              vktop(icelt,jcelt)=v(icelt,jcelt,ktop)
                              wktop(icelt,jcelt)=w(icelt,jcelt,ktop)
                           else
                              utotktp(icelt,jcelt)=sqrt(ualoft**2+valoft**2)
                              uktop(icelt,jcelt)=ualoft
                              vktop(icelt,jcelt)=valoft
                              wktop(icelt,jcelt)=0.
                           endif
                        endif
                        delut=sqrt((urefu(i,k)-u(icelt,jcelt,k))**2+(urefv(i,k)-v(icelt,jcelt,k))**2&
                              +(urefw(i,k)-w(icelt,jcelt,k))**2)
                        delutz=sqrt((u(icelt,jcelt,k)-zcorf(k)*uktop(icelt,jcelt))**2 &
                               +(v(icelt,jcelt,k)-zcorf(k)*vktop(icelt,jcelt))**2&
                               +(w(icelt,jcelt,k)-zcorf(k)*wktop(icelt,jcelt))**2)
!mdw 4-16-2004 added correction for ktop+3 > nz-1
                        if(k.ne.ktp)then
! mdw 4-16-2004 added proper treatment of zfo
! mdw 6-10-2004 changed to make check on centerline rather than local value
                           if((2.*deluc(i,k)/weff(i)).lt.(utotktp(icelt,jcelt)/(Ht(i)+zfo(i))) &
                                     .and.delutz.gt..2*zcorf(k)*utotktp(icelt,jcelt))then
                              if(ustargz(icelt,jcelt,k).lt.knlc*utotktp(icelt,jcelt))then
                                 ustargz(icelt,jcelt,k)=knlc*utotktp(icelt,jcelt)
                                 elzg(icelt,jcelt,k)=Ht(i)+zfo(i)
                                 upvpg=0.
                                 upwpg=-ctau13*zbrac*ustargz(icelt,jcelt,k)**2
                                 upsqg=cusq*zbrac*ustargz(icelt,jcelt,k)**2
                                 vpsqg=cvsq*zbrac*ustargz(icelt,jcelt,k)**2
                                 vpwpg=0.
                                 wpsqg=cwsq*zbrac*ustargz(icelt,jcelt,k)**2
                                 ustarg(icelt,jcelt,k)=ustargz(icelt,jcelt,k)
                                 call rotate2d(icelt,jcelt,k)
                                 if(iturbtypeflag.eq.1)nonlocal_option(icelt,jcelt,k)=4
                              endif
                           else
! We use the vertical gradient as dominant if it is sharper than the horizontal
                              cosl=u(icelt,jcelt,k)/utott
                              sinl=v(icelt,jcelt,k)/utott
                              duy=-deluc(i,k)*sinl*cosu+deluc(i,k)*sinv*cosl
! mdw 6-10-2004 changed check from delut (local value) to deluc(i,k); centerline
                              if(delut.gt..2*uref(i,k))then
                                 if(ustarg(icelt,jcelt,k).lt.knlc*deluc(i,k))then
                                    ustarg(icelt,jcelt,k)=knlc*deluc(i,k)
                                    upvpg=-(real(ist-1)/real(istf-1))*ustarg(icelt,jcelt,k)**2
                                    upwpg=0.
! for eddy transport in uv we don't consider uw
                                     upsqg=cusq*zbrac*ustarg(icelt,jcelt,k)**2
                                     wpsqg=cvsq*zbrac*ustarg(icelt,jcelt,k)**2
                                     vpwpg=0.
                                     vpsqg=cwsq*zbrac*ustarg(icelt,jcelt,k)**2
                                    elzg(icelt,jcelt,k)=.5*weff(i)
                                    call rotate2d(icelt,jcelt,k)
                                    if(iturbtypeflag.eq.1)nonlocal_option(icelt,jcelt,k)=5
                                 endif
                              endif
                           endif
                           if(is.eq.isini)then
                              do isin=isini+1,istinf
                                 xceln=xcelt-ds*real(isin-1)*cosphi
                                 yceln=ycelt-ds*real(isin-1)*sinphi
                                 iceln=nint(xceln/dx)+1
                                 jceln=nint(yceln/dy)+1
                                 iceln=min(nx-1,iceln)
                                 iceln=max(2,iceln)
                                 jceln=min(ny-1,jceln)
                                 jceln=max(2,jceln)
! mdw 3/22/2004PM added if statement to avoid replacing non-zero ustarg stuff
! with zero values
                                 if(ustarg(icelt,jcelt,k).gt.ustarg(iceln,jceln,k))then
                                    ustarg(iceln,jceln,k)=ustarg(icelt,jcelt,k)
                                    elzg(iceln,jceln,k)=elzg(icelt,jcelt,k)
                                    ufsqgi(iceln,jceln,k)=ufsqgi(icelt,jcelt,k)
                                    vfsqgi(iceln,jceln,k)=vfsqgi(icelt,jcelt,k)
                                    wfsqgi(iceln,jceln,k)=wfsqgi(icelt,jcelt,k)
                                    ufvfgi(iceln,jceln,k)=ufvfgi(icelt,jcelt,k)
                                    ufwfgi(iceln,jceln,k)=ufwfgi(icelt,jcelt,k)
                                    vfwfgi(iceln,jceln,k)=vfwfgi(icelt,jcelt,k)
                                    if(iturbtypeflag.eq.1)nonlocal_option(iceln,jceln,k)=nonlocal_option(icelt,jcelt,k)
                                 endif
! mdw 3/22/2004PM new endif for new if then
                              enddo
                           endif
                        endif
                     endif
! opposite direction in the transverse of the wake
                     xcelt=xcell-ds*real(ist-1)*cosphit
                     ycelt=ycell-ds*real(ist-1)*sinphit
                     icelt=nint(xcelt/dx)+1
                     jcelt=nint(ycelt/dy)+1
                     if(abs(xcelt-xcell).lt..5*ds)icelt=icel
                     if(abs(ycelt-ycell).lt..5*ds)jcelt=jcel
                     icelt=min(nx-1,icelt)
                     icelt=max(2,icelt)
                     jcelt=min(ny-1,jcelt)
                     jcelt=max(2,jcelt)
                     if(icellflag(icelt,jcelt,k) .gt. 0)then
                        utott=sqrt(u(icelt,jcelt,k)**2+v(icelt,jcelt,k)**2+w(icelt,jcelt,k)**2)
                        utott=utott+.000001
                        delut=sqrt((urefu(i,k)-u(icelt,jcelt,k))**2+(urefv(i,k)-v(icelt,jcelt,k))**2&
                              +(urefw(i,k)-w(icelt,jcelt,k))**2)
! mdw 4-16-2004 added correction for ktop+3 > nz-1
                        if(k.eq.ktp)then
                           if(ktop.le.nz-1)then
                              utotktp(icelt,jcelt)=utott
                              uktop(icelt,jcelt)=u(icelt,jcelt,ktop)
                              vktop(icelt,jcelt)=v(icelt,jcelt,ktop)
                              wktop(icelt,jcelt)=w(icelt,jcelt,ktop)
                           else
                              utotktp(icelt,jcelt)=sqrt(ualoft**2+valoft**2)
                              uktop(icelt,jcelt)=ualoft
                              vktop(icelt,jcelt)=valoft
                              wktop(icelt,jcelt)=0.
                           endif
                        endif
! mdw 4-16-2004 added correction for ktop+3 > nz-1
                        if(k.ne.ktp)then
                           delutz=sqrt((u(icelt,jcelt,k)-zcorf(k)*uktop(icelt,jcelt))**2 &
                                  +(v(icelt,jcelt,k)-zcorf(k)*vktop(icelt,jcelt))**2+(w(icelt,jcelt,k)&
                                  -zcorf(k)*wktop(icelt,jcelt))**2)
! mdw 4-16-2004 added proper treatment of zfo
! mdw 6-10-2004 made check on centerline rather than local value
                           if((2.*deluc(i,k)/weff(i)).lt.(utotktp(icelt,jcelt)/(Ht(i)+zfo(i))) &
                                        .and.delutz.gt..2*zcorf(k)*utotktp(icelt,jcelt))then
                              if(ustargz(icelt,jcelt,k).lt.knlc*utotktp(icelt,jcelt))then
                                 ustargz(icelt,jcelt,k)=knlc*utotktp(icelt,jcelt)
                                 elzg(icelt,jcelt,k)=Ht(i)+zfo(i)
                                 upvpg=0.
                                 upwpg=-ctau13*zbrac*ustargz(icelt,jcelt,k)**2
                                 upsqg=cusq*zbrac*ustargz(icelt,jcelt,k)**2
                                 vpsqg=cvsq*zbrac*ustargz(icelt,jcelt,k)**2
                                 vpwpg=0.
                                 wpsqg=cwsq*zbrac*ustargz(icelt,jcelt,k)**2
                                 ustarg(icelt,jcelt,k)=ustargz(icelt,jcelt,k)
                                 call rotate2d(icelt,jcelt,k)
                                 if(iturbtypeflag.eq.1)nonlocal_option(icelt,jcelt,k)=4
                              endif
                           else
! We use the vertical gradient as dominant if it is sharper than the horizontal
                              cosl=u(icelt,jcelt,k)/utott
                              sinl=v(icelt,jcelt,k)/utott
                              duy=-deluc(i,k)*sinl*cosu+deluc(i,k)*sinv*cosl
! mdw 6-10-2004 made check on centerline value rather than local value
                              if(delut.gt..2*uref(i,k))then
                                 if(ustarg(icelt,jcelt,k).lt.knlc*deluc(i,k))then
                                    ustarg(icelt,jcelt,k)=knlc*deluc(i,k)
!                                    duyi(icelt,jcelt,k)=duy
                                    upvpg=(real(ist-1)/real(istf-1))*ustarg(icelt,jcelt,k)**2
                                    upwpg=0.
! for eddy transport in uv we don't consider uw
                                    upvpg=ctau13*zbrac*(real(ist-1)/real(istf-1))*ustarg(icelt,jcelt,k)**2
                                    upwpg=0.
! for eddy transport in uv we don't consider uw
                                    upsqg=cusq*zbrac*ustarg(icelt,jcelt,k)**2
                                    wpsqg=cvsq*zbrac*ustarg(icelt,jcelt,k)**2
                                    elzg(icelt,jcelt,k)=0.5*weff(i)
                                    vpwpg=0.
                                    vpsqg=cwsq*zbrac*ustarg(icelt,jcelt,k)**2
                                    call rotate2d(icelt,jcelt,k)
                                    if(iturbtypeflag.eq.1)nonlocal_option(icelt,jcelt,k)=5
                                 endif
                              endif
                           endif
                           if(is.eq.isini)then
                              do isin=isini+1,istinf
                                 xceln=xcelt-ds*real(isin-1)*cosphi
                                 yceln=ycelt-ds*real(isin-1)*sinphi
                                 iceln=nint(xceln/dx)+1
                                 jceln=nint(yceln/dy)+1
                                 iceln=min(nx-1,iceln)
                                 iceln=max(2,iceln)
                                 jceln=min(ny-1,jceln)
                                 jceln=max(2,jceln)
! mdw 3/22/2004pm adding new if then structure to avoid replacing non-zero
! ustarg with zero ones
                                 if(ustarg(icelt,jcelt,k).gt.ustarg(iceln,jceln,k))then
                                    ustarg(iceln,jceln,k)=ustarg(icelt,jcelt,k)
                                    elzg(iceln,jceln,k)=elzg(icelt,jcelt,k)
                                    ufsqgi(iceln,jceln,k)=ufsqgi(icelt,jcelt,k)
                                    vfsqgi(iceln,jceln,k)=vfsqgi(icelt,jcelt,k)
                                    wfsqgi(iceln,jceln,k)=wfsqgi(icelt,jcelt,k)
                                    ufvfgi(iceln,jceln,k)=ufvfgi(icelt,jcelt,k)
                                    ufwfgi(iceln,jceln,k)=ufwfgi(icelt,jcelt,k)
                                    vfwfgi(iceln,jceln,k)=vfwfgi(icelt,jcelt,k)
                                    if(iturbtypeflag.eq.1)nonlocal_option(iceln,jceln,k)=nonlocal_option(icelt,jcelt,k)
                                 endif
! mdw 3/22/2004 end of new if then structure
                              enddo
                           endif
                        endif
                     endif
                  enddo   lp021
               enddo   lp022
               isini=1
lp024:         do is=1,isfu
! upstream front eddy along the centerline
                  xcell=xcu-ds*real(is-1)*cosphi
                  ycell=ycu-ds*real(is-1)*sinphi
!mdw 7-05-2006 changed form form =nint( / ) to nint( / )+1
                  icel=nint(xcell/dx)+1
                  jcel=nint(ycell/dy)+1
                  icel=min(nx-1,icel)
                  icel=max(2,icel)
                  jcel=min(ny-1,jcel)
                  jcel=max(2,jcel)
                  if(icellflag(icel,jcel,k) .eq. 0 .and. is .eq. 1)then
                     isini=2
                  endif
                  utot=sqrt(u(icel,jcel,k)**2+v(icel,jcel,k)**2+w(icel,jcel,k)**2)
! mdw 1-22-2004 new lines in support of bldg infiltration
                  if((k.eq.kmid).and.(is.eq.1))utotcl1(i)=utot
                  if((k.eq.kmid).and.(utot.gt.utotmax(i)))utotmax(i)=utot
                  utot=utot+.000001
!mdw 4-16-2004 added correction for ktop+3 > nz-1
                  if(k.eq.ktp)then
                     if(ktop.le.nz-1)then
                        utotktp(icel,jcel)=utot
                        uktop(icel,jcel)=u(icel,jcel,ktop)
                        vktop(icel,jcel)=v(icel,jcel,ktop)
                        wktop(icel,jcel)=w(icel,jcel,ktop)
                     else
                        utotktp(icel,jcel)=sqrt(ualoft**2+valoft**2)
                        uktop(icel,jcel)=ualoft
                        vktop(icel,jcel)=valoft
                        wktop(icel,jcel)=0.
                     endif
                  endif
                  deluc(i,k)=sqrt((urefu(i,k)-u(icel,jcel,k))**2+(urefv(i,k)-v(icel,jcel,k))**2&
                             +(urefw(i,k)-w(icel,jcel,k))**2)
!mdw 4-16-2004 added correction for ktop+3 > nz-1
                  if(k.ne.ktp)then
                     delutz=sqrt((u(icel,jcel,k)-zcorf(k)*uktop(icel,jcel))**2+(v(icel,jcel,k)&
                            -zcorf(k)*vktop(icel,jcel))**2+(w(icel,jcel,k)-zcorf(k)*wktop(icel,jcel))**2)
                     deluc(i,k)=sqrt((urefu(i,k)-u(icel,jcel,k))**2+(urefv(i,k)-v(icel,jcel,k))**2&
                                +(urefw(i,k)-w(icel,jcel,k))**2)
! Selects the largest gradient (vert or horiz transfer)
! mdw 4-16-2004 added proper treatment of zfo
                     if((2.*deluc(i,k)/weff(i)).lt.(utotktp(icel,jcel)/(Ht(i)+zfo(i))) &
                               .and.delutz.gt..2*zcorf(k)*utotktp(icel,jcel))then ! vertical dominates
                        if(ustargz(icel,jcel,k).lt.knlc*utotktp(icel,jcel))then ! This value dominates over prev. buildings.
                           ustargz(icel,jcel,k)=knlc*utotktp(icel,jcel)
                           elzg(icel,jcel,k)=Ht(i)+zfo(i)
                           upvpg=0.
                           upwpg=-ctau13*zbrac*ustargz(icel,jcel,k)**2
                           upsqg=cusq*zbrac*ustargz(icel,jcel,k)**2
                           vpsqg=cvsq*zbrac*ustargz(icel,jcel,k)**2
                           vpwpg=0.
                           wpsqg=cwsq*zbrac*ustargz(icel,jcel,k)**2 ! align sigmas with the overall mean wind
                           ustarg(icel,jcel,k)=ustargz(icel,jcel,k)
                           call rotate2d(icel,jcel,k) ! calculates the quantites in the original coord sys.
                           if(iturbtypeflag.eq.1)nonlocal_option(icel,jcel,k)=4
                        endif
                     else
! We use the vertical gradient as dominant if it is sharper than the horizontal
                        cosl=u(icel,jcel,k)/utot
                        sinl=v(icel,jcel,k)/utot
                        duy=-deluc(i,k)*sinl*cosu+deluc(i,k)*sinv*cosl
! we now have the delta u between the outside of the bldg and the center of the wake
                        if(deluc(i,k).gt..2*uref(i,k))then
                           if(ustarg(icel,jcel,k).lt.knlc*deluc(i,k))then
                              ustarg(icel,jcel,k)=knlc*deluc(i,k)
                              upvpg=0.
! on axis u prime v prime is zero
                              upwpg=0.
! for eddy transport in uv we don't consider uw
                              upsqg=cusq*zbrac*ustarg(icel,jcel,k)**2
                              wpsqg=cvsq*zbrac*ustarg(icel,jcel,k)**2
                              vpwpg=0.
                              elzg(icel,jcel,k)=0.5*weff(i)
                              vpsqg=cwsq*zbrac*ustarg(icel,jcel,k)**2
                              call rotate2d(icel,jcel,k)
                              if(iturbtypeflag.eq.1)nonlocal_option(icel,jcel,k)=5
                           endif
                        endif
                     endif
                  endif
lp023:            do ist=2,istf
! first direction in the transverse of the front eddy
                     xcelt=xcell+ds*real(ist-1)*cosphit
                     ycelt=ycell+ds*real(ist-1)*sinphit
!mdw 7-05-2006 changed form from nint( / ) to nint( / )+1
                     icelt=nint(xcelt/dx)+1
                     jcelt=nint(ycelt/dy)+1
                     if(abs(xcelt-xcell).lt..5*ds)icelt=icel
                     if(abs(ycelt-ycell).lt..5*ds)jcelt=jcel
!mdw 7-11-2006 check added to use closest axis cell
                     icelt=min(nx-1,icelt)
                     icelt=max(2,icelt)
                     jcelt=min(ny-1,jcelt)
                     jcelt=max(2,jcelt)
                     utott=sqrt(u(icelt,jcelt,k)**2+v(icelt,jcelt,k)**2+w(icelt,jcelt,k)**2)
                     utott=utott+.000001
                     delut=sqrt((urefu(i,k)-u(icelt,jcelt,k))**2+(urefv(i,k)-v(icelt,jcelt,k))**2&
                           +(urefw(i,k)-w(icelt,jcelt,k))**2)
!mdw 4-16-2004 added correction for ktop+3 > nz-1
                     if(k.eq.ktp)then
                        if(ktop.le.nz-1)then
                           utotktp(icelt,jcelt)=utott
                           uktop(icelt,jcelt)=u(icelt,jcelt,ktop)
                           vktop(icelt,jcelt)=v(icelt,jcelt,ktop)
                           wktop(icelt,jcelt)=w(icelt,jcelt,ktop)
                        else
                           utotktp(icelt,jcelt)=sqrt(ualoft**2+valoft**2)
                           uktop(icelt,jcelt)=ualoft
                           vktop(icelt,jcelt)=valoft
                           wktop(icelt,jcelt)=0.
                        endif
                     endif
!mdw 4-16-2004 added correction for ktop+3 > nz-1
                     if(k.ne.ktp)then
                        delutz=sqrt((u(icelt,jcelt,k)-zcorf(k)*uktop(icelt,jcelt))**2 &
                               +(v(icelt,jcelt,k)-zcorf(k)*vktop(icelt,jcelt))**2&
                               +(w(icelt,jcelt,k)-zcorf(k)*wktop(icelt,jcelt))**2)
! mdw 4-16-2004 added proper treatment of zfo
! mdw 6-10-2004 made check on centerline deluc rather than local delut
                        if((2.*deluc(i,k)/weff(i)).lt.(utotktp(icelt,jcelt)/(Ht(i)+zfo(i))) &
                                  .and.delutz.gt..2*zcorf(k)*utotktp(icelt,jcelt))then
                           if(ustargz(icelt,jcelt,k).lt.knlc*utotktp(icelt,jcelt))then
                              ustargz(icelt,jcelt,k)=knlc*utotktp(icelt,jcelt)
                              elzg(icelt,jcelt,k)=Ht(i)+zfo(i)
                              upvpg=0.
                              upwpg=-ctau13*zbrac*ustargz(icelt,jcelt,k)**2
                              upsqg=cusq*zbrac*ustargz(icelt,jcelt,k)**2
                              vpsqg=cvsq*zbrac*ustargz(icelt,jcelt,k)**2
                              vpwpg=0.
                              wpsqg=cwsq*zbrac*ustargz(icelt,jcelt,k)**2
                              ustarg(icelt,jcelt,k)=ustargz(icelt,jcelt,k)
                              call rotate2d(icelt,jcelt,k)
                              if(iturbtypeflag.eq.1)nonlocal_option(icelt,jcelt,k)=4
                           endif
                        else
! We use the vertical gradient as dominant if it is sharper than the horizontal
                           cosl=u(icelt,jcelt,k)/utott
                           sinl=v(icelt,jcelt,k)/utott
                           duy=-deluc(i,k)*sinl*cosu+deluc(i,k)*sinv*cosl
! mdw 6-10-2004 made check on centerline rather than local value
                           if(delut.gt..2*uref(i,k))then
                              if(ustarg(icelt,jcelt,k).lt.knlc*deluc(i,k))then
                                 ustarg(icelt,jcelt,k)=knlc*deluc(i,k)
! for eddy transport in uv we don't consider uw
                                 upvpg=-ctau13*zbrac*(real(ist-1)/real(istf-1))*ustarg(icelt,jcelt,k)**2
                                 upwpg=0.
! for eddy transport in uv we don't consider uw
                                 upsqg=cusq*zbrac*ustarg(icelt,jcelt,k)**2
                                 wpsqg=cvsq*zbrac*ustarg(icelt,jcelt,k)**2
                                 vpwpg=0.
                                 vpsqg=cwsq*zbrac*ustarg(icelt,jcelt,k)**2
                                 elzg(icelt,jcelt,k)=0.5*weff(i)
                                 call rotate2d(icelt,jcelt,k)
                                 if(iturbtypeflag.eq.1)nonlocal_option(icelt,jcelt,k)=5
                              endif
                           endif
                        endif
                        if(is.eq.isini)then
                           do isin=isini+1,istinf
                              xceln=xcelt+ds*real(isin-1)*cosphi
                              yceln=ycelt+ds*real(isin-1)*sinphi
                              iceln=nint(xceln/dx)+1
                              jceln=nint(yceln/dy)+1
                              iceln=min(nx-1,iceln)
                              iceln=max(2,iceln)
                              jceln=min(ny-1,jceln)
                              jceln=max(2,jceln)
! mdw 3/22/2004pm added new if then structure to prevent replacing non-zero
! ustarg s with zero ones
                              if(ustarg(icelt,jcelt,k).gt.ustarg(iceln,jceln,k))then
                                 ustarg(iceln,jceln,k)=ustarg(icelt,jcelt,k)
                                 elzg(iceln,jceln,k)=elzg(icelt,jcelt,k)
                                 ufsqgi(iceln,jceln,k)=ufsqgi(icelt,jcelt,k)
                                 vfsqgi(iceln,jceln,k)=vfsqgi(icelt,jcelt,k)
                                 wfsqgi(iceln,jceln,k)=wfsqgi(icelt,jcelt,k)
                                 ufvfgi(iceln,jceln,k)=ufvfgi(icelt,jcelt,k)
                                 ufwfgi(iceln,jceln,k)=ufwfgi(icelt,jcelt,k)
                                 vfwfgi(iceln,jceln,k)=vfwfgi(icelt,jcelt,k)
                                 if(iturbtypeflag.eq.1)nonlocal_option(iceln,jceln,k)=nonlocal_option(icelt,jcelt,k)
                              endif
! mdw 3/22/2004pm end  of new if then structure
                           enddo
                        endif
                     endif
! opposite direction in the transverse of the front eddy
                     xcelt=xcell-ds*real(ist-1)*cosphit
                     ycelt=ycell-ds*real(ist-1)*sinphit
!mdw 7-05-2006 changed form from nint( / ) to nint( / )+1
                     icelt=nint(xcelt/dx)+1
                     jcelt=nint(ycelt/dy)+1
                     if(abs(xcelt-xcell).lt..5*ds)icelt=icel
                     if(abs(ycelt-ycell).lt..5*ds)jcelt=jcel
!mdw 7-11-2006 check added to use closest axis cell
                     icelt=min(nx-1,icelt)
                     icelt=max(2,icelt)
                     jcelt=min(ny-1,jcelt)
                     jcelt=max(2,jcelt)
                     utott=sqrt(u(icelt,jcelt,k)**2+v(icelt,jcelt,k)**2+w(icelt,jcelt,k)**2)
                     utott=utott+.000001
!mdw 4-16-2004 added correction for ktop+3 > nz-1
                     if(k.eq.ktp)then
                        if(ktop.le.nz-1)then
                           utotktp(icelt,jcelt)=utott
                           uktop(icelt,jcelt)=u(icelt,jcelt,ktop)
                           vktop(icelt,jcelt)=v(icelt,jcelt,ktop)
                           wktop(icelt,jcelt)=w(icelt,jcelt,ktop)
                        else
                           utotktp(icelt,jcelt)=sqrt(ualoft**2+valoft**2)
                           uktop(icelt,jcelt)=ualoft
                           vktop(icelt,jcelt)=valoft
                           wktop(icelt,jcelt)=0.
                        endif
                     endif
                     delut=sqrt((urefu(i,k)-u(icelt,jcelt,k))**2+(urefv(i,k)-v(icelt,jcelt,k))**2&
                           +(urefw(i,k)-w(icelt,jcelt,k))**2)
!mdw 4-16-2004 added correction for ktop+3 > nz-1
                     if(k.ne.ktp)then
                        delutz=sqrt((u(icelt,jcelt,k)-zcorf(k)*uktop(icelt,jcelt))**2 &
                               +(v(icelt,jcelt,k)-zcorf(k)*vktop(icelt,jcelt))**2+ &
                               (w(icelt,jcelt,k)-zcorf(k)*wktop(icelt,jcelt))**2)
! mdw 4-16-2004 added proper treatment of zfo
! mdw 6-10-2004 made check on centerline rather than local value
                        if((2.*deluc(i,k)/weff(i)).lt.(utotktp(icelt,jcelt)/(Ht(i)+zfo(i))) &
                                    .and.delutz.gt..2*zcorf(k)*utotktp(icelt,jcelt))then
                           if(ustargz(icelt,jcelt,k).lt.knlc*utotktp(icelt,jcelt))then
                              ustargz(icelt,jcelt,k)=knlc*utotktp(icelt,jcelt)
                              elzg(icelt,jcelt,k)=Ht(i)+zfo(i)
                              upvpg=0.
                              upwpg=-ctau13*zbrac*ustargz(icelt,jcelt,k)**2
                              upsqg=cusq*zbrac*ustargz(icelt,jcelt,k)**2
                              vpsqg=cvsq*zbrac*ustargz(icelt,jcelt,k)**2
                              vpwpg=0.
                              wpsqg=cwsq*zbrac*ustargz(icelt,jcelt,k)**2
                              ustarg(icelt,jcelt,k)=ustargz(icelt,jcelt,k)
                              call rotate2d(icelt,jcelt,k)
                              if(iturbtypeflag.eq.1)nonlocal_option(icelt,jcelt,k)=4
                           endif
                        else
! We use the vertical gradient as dominant if it is sharper than the horizontal
                           cosl=u(icelt,jcelt,k)/utott
                           sinl=v(icelt,jcelt,k)/utott
                           duy=-deluc(i,k)*sinl*cosu+deluc(i,k)*sinv*cosl
! mdw 6-10-2004 made check on centerline rather than local (delut) value
                           if(delut.gt..2*uref(i,k))then
                              if(ustarg(icelt,jcelt,k).lt.knlc*deluc(i,k).and.ustargz(icelt,jcelt,k).lt.knlc*deluc(i,k))then
                                 ustarg(icelt,jcelt,k)=knlc*deluc(i,k)
! for eddy transport in uv we don't consider uw
                                 upvpg=tau13*zbrac*(real(ist-1)/real(istf-1))*ustarg(icelt,jcelt,k)**2
                                 upwpg=0.
! for eddy transport in uv we don't consider uw
                                 upsqg=cusq*zbrac*ustarg(icelt,jcelt,k)**2
                                 wpsqg=cvsq*zbrac*ustarg(icelt,jcelt,k)**2
                                 vpwpg=0.
                                 elzg(icelt,jcelt,k)=0.5*weff(i)
                                 vpsqg=cwsq*zbrac*ustarg(icelt,jcelt,k)**2
                                 call rotate2d(icelt,jcelt,k)
                                 if(iturbtypeflag.eq.1)nonlocal_option(icelt,jcelt,k)=5
                              endif
                           endif
                        endif
                        if(is.eq.isini)then
                           do isin=isini+1,istinf
                              xceln=xcelt+ds*real(isin-1)*cosphi
                              yceln=ycelt+ds*real(isin-1)*sinphi
                              iceln=nint(xceln/dx)+1
                              jceln=nint(yceln/dy)+1
                              iceln=min(nx-1,iceln)
                              iceln=max(2,iceln)
                              jceln=min(ny-1,jceln)
                              jceln=max(2,jceln)
! mdw 3/22/2004pm added new if then structure to prevent replacing non-zero
! ustargs with zero ones
                              if(ustarg(icelt,jcelt,k).gt.ustarg(iceln,jceln,k))then
                                 ustarg(iceln,jceln,k)=ustarg(icelt,jcelt,k)
                                 elzg(iceln,jceln,k)=elzg(icelt,jcelt,k)
                                 ufsqgi(iceln,jceln,k)=ufsqgi(icelt,jcelt,k)
                                 vfsqgi(iceln,jceln,k)=vfsqgi(icelt,jcelt,k)
                                 wfsqgi(iceln,jceln,k)=wfsqgi(icelt,jcelt,k)
                                 ufvfgi(iceln,jceln,k)=ufvfgi(icelt,jcelt,k)
                                 ufwfgi(iceln,jceln,k)=ufwfgi(icelt,jcelt,k)
                                 vfwfgi(iceln,jceln,k)=vfwfgi(icelt,jcelt,k)
                                 if(iturbtypeflag.eq.1)nonlocal_option(iceln,jceln,k)=nonlocal_option(icelt,jcelt,k)
                              endif
! mdw 3/22/2004pm end of new if then structure
                           enddo
                        endif
                     endif
                  enddo   lp023
               enddo   lp024
            enddo   lp025
            select case(bldtype(i))
               case(3)
                  xpent1=xfo(i)-wti(i)*.2-dx
                  npentx=nint((.4*wti(i)/dx))+1
                  ypent1=yfo(i)-wti(i)*.2-dy
                  npenty=nint((.4*wti(i)/dy))+1
                  ipent1=nint((xpent1)/dx)+1
                  ipent2=ipent1+npentx
                  jpent1=nint((ypent1/dy))+1
                  jpent2=jpent1+npenty
                  do icel=ipent1,ipent2
                     do jcel=jpent1,jpent2
                        do k=ktp,2,-1
                           utot=sqrt(u(icel,jcel,k)**2+v(icel,jcel,k)**2+w(icel,jcel,k)**2)
                           utot=utot+.000001
!mdw 4-16-2004 added correction for ktop+3 > nz-1
                           if(k.eq.ktp)then
                              if(ktop.le.nz-1)then
                                 utotktp(icel,jcel)=utot
                                 uktop(icel,jcel)=u(icel,jcel,ktop)
                                 vktop(icel,jcel)=v(icel,jcel,ktop)
                                 wktop(icel,jcel)=w(icel,jcel,ktop)
                              else
                                 utotktp(icel,jcel)=sqrt(ualoft**2+valoft**2)
                                 uktop(icel,jcel)=ualoft
                                 vktop(icel,jcel)=valoft
                                 wktop(icel,jcel)=0.
                              endif
                           endif
                           if(k.ne.ktp.and.icellflag(icel,jcel,k) .ne. 0)then
! MAN 9/14/2005 pentagon courtyard nonlocal mixing fix
                              delutz=sqrt((u(icel,jcel,k)-uktop(icel,jcel))**2 &
                                     +(v(icel,jcel,k)-vktop(icel,jcel))**2+ &
                                     (w(icel,jcel,k)-wktop(icel,jcel))**2)
                              if(delutz.gt..2*utotktp(icel,jcel))then ! vertical dominates
! end MAN 9/14/2005           
                                 if(ustargz(icel,jcel,k).lt.knlc*utotktp(icel,jcel))then ! This value dominates over prev. buildings.
                                    ustargz(icel,jcel,k)=knlc*utotktp(icel,jcel)
                                    if(iturbtypeflag.eq.1)nonlocal_option(icel,jcel,k)=4
                                    elzg(icel,jcel,k)=Ht(i)+zfo(i)
                                    upvpg=0.
                                    upwpg=-ustargz(icel,jcel,k)**2
                                    upsqg=6.25*ustargz(icel,jcel,k)**2
                                    vpsqg=(4./6.25)*upsqg
                                    vpwpg=0.
                                    wpsqg=1.69*ustargz(icel,jcel,k)**2 ! align sigmas with the overall mean wind
                                    upvpg=0.
                                    upwpg=-ctau13*zbrac*ustargz(icel,jcel,k)**2
                                    upsqg=cusq*zbrac*ustargz(icel,jcel,k)**2
                                    vpsqg=cvsq*zbrac*ustargz(icel,jcel,k)**2
                                    vpwpg=0.
                                    wpsqg=cwsq*zbrac*ustargz(icel,jcel,k)**2 ! align sigmas with the overall mean wind
                                    ustarg(icel,jcel,k)=ustargz(icel,jcel,k)
                                    call rotate2d(icel,jcel,k) ! calculates the quantites in the original coord sys.
                                 endif
                              endif
                           endif
                        enddo
                        check1=check1+1
                     enddo
                     check2=check2+1
                  enddo
                  check3=check3+1
               case(4,5)
                  x0=xfo(i)+0.5*lt(i)*cos(gamma(i)*pi/180.)
                  y0=yfo(i)+0.5*lt(i)*sin(gamma(i)*pi/180.)
                  x1=xfo(i)+0.5*wti(ibuild)*sin(gamma(i)*pi/180.)
                  y1=yfo(i)-0.5*wti(ibuild)*cos(gamma(i)*pi/180.)
                  x2=x1+lt(i)*cos(gamma(i)*pi/180.)
                  y2=y1+lt(i)*sin(gamma(i)*pi/180.)
                  x4=xfo(i)-0.5*wti(i)*sin(gamma(i)*pi/180.)
                  y4=yfo(i)+0.5*wti(i)*cos(gamma(i)*pi/180.)
                  x3=x4+lt(i)*cos(gamma(i)*pi/180.)
                  y3=y4+lt(i)*sin(gamma(i)*pi/180.)
                  
                  do icel=int(min(x1,x2,x3,x4)/dx),int(max(x1,x2,x3,x4)/dx)+1
                     do jcel=int(min(y1,y2,y3,y4)/dy),int(max(y1,y2,y3,y4)/dy)+1
                        xc=((real(icel)-0.5)*dx-x0)*cos(gamma(i)*pi/180.)+&
                           ((real(jcel)-0.5)*dy-y0)*sin(gamma(i)*pi/180.)
                        yc=-((real(icel)-0.5)*dx-x0)*sin(gamma(i)*pi/180.)+&
                           ((real(jcel)-0.5)*dy-y0)*cos(gamma(i)*pi/180.)
                        do k=ktp,nint(zfo(i)/dz)+2,-1
                           incourt=0
                           if(icellflag(icel,jcel,k) .ne. 0)then
                              utot=sqrt(u(icel,jcel,k)**2+v(icel,jcel,k)**2+w(icel,jcel,k)**2)+.000001
                              if(bldtype(i) .eq. 4)then
                                 if(xc .gt. -0.5*lt(i) .and. xc .lt. 0.5*lt(i) .and. &
                                       yc .gt. -0.5*wti(i) .and. yc .lt. 0.5*wti(i))then
                                    incourt=1
                                 endif
                              else
                                 rc=sqrt((xc**2.)+(yc**2.))
                                 tc=atan2(yc,xc)
                                 if(rc .lt. 0.25*lt(i)*wti(i)/&
                                       sqrt(((0.5*lt(i)*sin(tc))**2.)+((0.5*wti(i)*cos(tc))**2.)))then
                                    incourt=1
                                 endif
                              endif
                           else
                              cycle
                           endif
!mdw 4-16-2004 added correction for ktop+3 > nz-1
                           if(incourt .eq. 1)then
                              if(k.eq.ktp)then
                                 if(ktop.le.nz-1)then
                                    utotktp(icel,jcel)=utot
                                    uktop(icel,jcel)=u(icel,jcel,ktop)
                                    vktop(icel,jcel)=v(icel,jcel,ktop)
                                    wktop(icel,jcel)=w(icel,jcel,ktop)
                                 else
                                    utotktp(icel,jcel)=sqrt(ualoft**2+valoft**2)
                                    uktop(icel,jcel)=ualoft
                                    vktop(icel,jcel)=valoft
                                    wktop(icel,jcel)=0.
                                 endif
                              endif
                              if(k.ne.ktp.and.icellflag(icel,jcel,k) .ne. 0)then
! MAN 9/14/2005 pentagon courtyard nonlocal mixing fix
                                 delutz=sqrt((u(icel,jcel,k)-uktop(icel,jcel))**2 &
                                        +(v(icel,jcel,k)-vktop(icel,jcel))**2+ &
                                        (w(icel,jcel,k)-wktop(icel,jcel))**2)
                                 if(delutz.gt..2*utotktp(icel,jcel))then ! vertical dominates
! end MAN 9/14/2005              
                                    if(ustargz(icel,jcel,k).lt.knlc*utotktp(icel,jcel))then ! This value dominates over prev. buildings.
                                       ustargz(icel,jcel,k)=knlc*utotktp(icel,jcel)
                                       if(iturbtypeflag.eq.1)nonlocal_option(icel,jcel,k)=4
                                       elzg(icel,jcel,k)=Ht(i)+zfo(i)
                                       upvpg=0.
                                       upwpg=-ustargz(icel,jcel,k)**2
                                       upsqg=6.25*ustargz(icel,jcel,k)**2
                                       vpsqg=(4./6.25)*upsqg
                                       vpwpg=0.
                                       wpsqg=1.69*ustargz(icel,jcel,k)**2 ! align sigmas with the overall mean wind
                                       upvpg=0.
                                       upwpg=-ctau13*zbrac*ustargz(icel,jcel,k)**2
                                       upsqg=cusq*zbrac*ustargz(icel,jcel,k)**2
                                       vpsqg=cvsq*zbrac*ustargz(icel,jcel,k)**2
                                       vpwpg=0.
                                       wpsqg=cwsq*zbrac*ustargz(icel,jcel,k)**2 ! align sigmas with the overall mean wind
                                       ustarg(icel,jcel,k)=ustargz(icel,jcel,k)
                                       call rotate2d(icel,jcel,k) ! calculates the quantites in the original coord sys.
                                    endif
                                 endif
                              endif
                           endif
                        enddo
                        check1=check1+1
                     enddo
                     check2=check2+1
                  enddo
                  check3=check3+1
            endselect
            check=check+1
         enddo   lp026
         if(format_flag.eq.1.or.format_flag.eq.3)then
            nx2=nx-1
            ny2=ny-1
            nz2=nz-2
            if(iturbfieldflag.eq.1)then
               write(13,22120)
               write(13,22130)
               write(13,22135)
               write(13,21150)nx2,ny2,nz2
               write(13,21160)
            endif
         endif
 22120   format('TITLE     = "TURBULENCE in 3D"')
 22130   format('VARIABLES = "X" "Y" "Z" "SIGU" "SIGV" "SIGW" "LZ" "LEFF" "EPS" "UV" "UW" "VW"')
 22135   format('ZONE T = "TURBULENCE" ')
 21150   format('I=',i6,', J=',i6,',  K=',i6,', F=POINT')
 21160   format('DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)')
! calculate distance to ground and walls if within 2 cells
lp029:   do k=2,nz-1
            zbrac=(1.-zi(k)/h)**1.5
lp028:      do j=1,ny-1
lp027:         do i=1,nx-1
! for vertical downward distance (dzm)
                  klim=1
                  dzm(i,j,k)=-.5*dz+real(k-1)*dz
! MAN 9/21/2005 roof top mixing length fix
                  eleff(i,j,k)=dzm(i,j,k)
                  kdif=k-1
                  do kk=k,klim,-1
! calculation of dzm; the distance from the cell center to the
! nearest horizontal surface
                     if(icellflag(i,j,kk) .eq. 0)then
! MAN 9/21/2005 roof top mixing length fix
                        eleff(i,j,k)=dzm(i,j,k)-real(kk-1)*dz*(real(kk-1)*dz/dzm(i,j,k))**m_roof
                        dzm(i,j,k)=.5*dz+real(k-kk-1)*dz
                        kdif=k-kk
                        exit
                     endif
                  enddo
! for vertical upward distance (dzp)
                  klim=nz-1
! calculation of ustar in the vertical
                  if(icellflag(i,j,k-1) .eq. 0 .and. icellflag(i,j,k) .ne. 0)then
                     utot=sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
                     if(rcl.gt.0)then
                        phim=1.+4.7*rcl*0.5*dz
                        psim=-4.7*rcl*0.5*dz
                     else
                        phim=(1.-15.*rcl*0.5*dz)**(-.25)
                        psim=2.*alog((1.+1./phim)/2.)+alog((1.+1./phim**2)/2.)-2.*atan(1./phim)+pi/2.
                     endif
                     ustar=kkar*utot/(log(.5*dz/z0)-psim)
                     dutotdzi(i,j,k)=ustar*phim/(kkar*.5*dz)
                     ustarz(i,j,k)=elz(i,j,k)*dutotdzi(i,j,k)/phim
                     sigwi(i,j,k)=1.3*ustarz(i,j,k)
                  else
                     if(icellflag(i,j,k) .ne. 0)then
                        utot=sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
                        if(abs(dutotdzi(i,j,k)).gt.1.e-06)then
                           elz(i,j,k)=kkar*utot/abs(dutotdzi(i,j,k))
! MAN 9/21/2005 roof top mixing length fix
                           if((kkar*eleff(i,j,k)).lt.elz(i,j,k)) elz(i,j,k)=kkar*eleff(i,j,k)
                        else
                           elz(i,j,k)=kkar*eleff(i,j,k)
                        endif
!                        if(inextgridflag.eq.1)then
                           if(rcl.gt.0)then
                              phim=1.+4.7*rcl*eleff(i,j,k)
                              psim=-4.7*rcl*eleff(i,j,k)
                           else
                              phim=(1.-15.*rcl*eleff(i,j,k))**(-.25)
                              psim=2.*alog((1.+1./phim)/2.)+alog((1.+1./phim**2)/2.)-2.*atan(1./phim)+pi/2.
                           endif
!                           ustar=utot*kkar/(log(eleff(i,j,k)/z0)-psim)
                            ustar=kkar*eleff(i,j,k)*dutotdzi(i,j,k)/phim
!                           dutotdzi(i,j,k)=ustar*phim/(kkar*eleff(i,j,k))
                           ustarz(i,j,k)=ustar
!                        endif
                     endif
!                     ustarz(i,j,k)=elz(i,j,k)*abs(dutotdzi(i,j,k))
                  endif
! for neutral conditions sigw is only dependent on ustar
                  dzp(i,j,k)=.5*dz+real(klim-k)*dz+10.*dz
                  do kk=k,klim
                     if(icellflag(i,j,kk) .eq. 0)then
                        dzp(i,j,k)=.5*dz+real(kk-k-1)*dz
                        exit
                     endif
                  enddo
!23456789112345678921234567893123456789412345678951234567896123456789712
! for distance to the left (dxm)
                  ilim=1
                  dxm(i,j,k)=.5*dx+real(i-1)*dx+real(nx)*dx
                  do ii=i,ilim,-1
! calculation of the distance to the wall in the negative x direction
                     if(icellflag(ii,j,k) .eq. 0)then
                        dxm(i,j,k)=.5*dx+real(i-ii-1)*dx
                        exit
                     endif
                  enddo
! for distance to the right (dxp)
                  ilim=nx-1
                  dxp(i,j,k)=.5*dx+real(ilim-i)*dx+real(nx)*dx
                  do ii=i,ilim
! calculation of the distance to the wall in the positive x direction
                     if(icellflag(ii,j,k) .eq. 0)then
                        dxp(i,j,k)=.5*dx+real(ii-i-1)*dx
                        exit
                     endif
                  enddo
! for distance  from the back (dym)
                  jlim=1
                  dym(i,j,k)=.5*dy+real(j-1)*dy+real(ny)*dy
                  do jj=j,jlim,-1
! calculation of the distance to the wall in the negative y direction
                     if(icellflag(i,jj,k) .eq. 0)then
                        dym(i,j,k)=.5*dy+real(j-jj-1)*dy
                        exit
                     endif
                  enddo
! for distance to the front  (dyp)
                  jlim=ny-1
                  dyp(i,j,k)=.5*dy+real(jlim-j)*dy+real(ny)*dy
                  do jj=j,jlim
! calculation of the distance to the wall in the positive x direction
                     if(icellflag(i,jj,k) .eq. 0)then
                        dyp(i,j,k)=.5*dy+real(jj-j-1)*dy
                        exit
                     endif
                  enddo
! we need to calculate the largest change in utot
                  if(icellflag(i,j,k) .eq. 0)then
                     eps=0.
                     sigu=0.
                     sigv=0.
                     sigw=0.
                     upwp=0.
                     elz(i,j,k)=0.
                     eleff(i,j,k)=0.
                  endif
                  if(icellflag(i,j,k) .ne. 0)then
! first we set up parameters for cells near boundary
                     if(j.lt.2.or.j.ge.ny-1.or.i.lt.2.or.i.ge.nx-1)then
! calculation of near-boundary values of u*y, ly, dely, and the
! gradients of speed in the x and y directions
                        delym=real(ny)*dy
                        utot=sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
                        sigvi(i,j,k)=0.
                        delxm=dx*real(nx-2)
                        dutotdxi(i,j,k)=0.
                        dutotdyi(i,j,k)=0.
                        elz(i,j,k)=kkar*eleff(i,j,k)
                        dutotdni(i,j,k)=dutotdzi(i,j,k)
                        call detang(0)
                        if(rcl.gt.0)then
                           phim=1.+4.7*rcl*eleff(i,j,k)
                           psim=-4.7*rcl*eleff(i,j,k)
                        else
                           phim=(1.-15.*rcl*eleff(i,j,k))**(-.25)
                           psim=2.*alog((1.+1./phim)/2.)+alog((1.+1./phim**2)/2.)-2.*atan(1./phim)+pi/2.
                        endif
                        ustarz(i,j,k)=elz(i,j,k)*dutotdni(i,j,k)/phim ! calculate local ustar
                        ustarz(i,j,k)=max(ustarz(i,j,k),3.e-02)
                        u3psq=cusq*zbrac*ustarz(i,j,k)**2   ! (u''')^2
                        v3psq=cvsq*zbrac*ustarz(i,j,k)**2   !...
                        w3psq=cwsq*zbrac*ustarz(i,j,k)**2 !...
                        upwp=-ctau13*zbrac*ustarz(i,j,k)**2 ! -tau13
                        upvp=0.
                        vpwp=0.
                        if(rcl.lt.0.)then
                           u3psq=u3psq+.6*(ustarz(i,j,k)**2)*(-h*rcl)**(2./3.)
                           v3psq=v3psq+.6*(ustarz(i,j,k)**2)*(-h*rcl)**(2./3.)
                           w3psq=w3psq+3.3*(ustarz(i,j,k)**2)*(-zi(k)*rcl)**(2.3)*(1.-.8*zi(k)/h)**2
                           upwp=upwp*(1.-zi(k)/h)**(.5*rcl*h/(1-rcl*h))
                        endif
                        call rotu3psq ! rotate values back into the orig. grid
                        ufwfi(i,j,k)=ufwf
                        ufvfi(i,j,k)=ufvf
                        vfwfi(i,j,k)=vfwf
                        ustarij(i,j,k)=ustarz(i,j,k)
                        sigui(i,j,k)=sqrt(u3psq)
                        sigvi(i,j,k)=sqrt(v3psq)
                        sigwi(i,j,k)=sqrt(w3psq)
                        upwpi(i,j,k)=upwp
! along the boundaries we make y effects negligible
                     else
                        utot=sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
! away from boundaries u*y, ly, dely, and gradients
                        if(icellflag(i-1,j,k) .ne. 0 .and. icellflag(i,j,k) .ne. 0 &
                              .and. icellflag(i+1,j,k) .ne. 0)then
!mdw 3-08-2004 start changes for highest gradient
                           utot=sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
                           utotm=sqrt(u(i-1,j,k)**2+v(i-1,j,k)**2+w(i-1,j,k)**2)
                           utotp=sqrt(u(i+1,j,k)**2+v(i+1,j,k)**2+w(i+1,j,k)**2)
                           dutotdxp=(utotp-utot)/dx
                           dutotdxm=(utot-utotm)/dx
                           dutotdxa=max(abs(dutotdxp),abs(dutotdxm))
                           if(dutotdxa.eq.abs(dutotdxm))then
                              dutotdxi(i,j,k)=dutotdxm
                           else
                              dutotdxi(i,j,k)=dutotdxp
                           endif
! mdw 3-08-2004end changes
                        else
                           if(icellflag(i,j,k) .eq. 0)then !!BALLI
                              dutotdxi(i,j,k)=0.
                           else
                              if(icellflag(i-1,j,k) .eq. 0)then !!BALLI
                                 dutotdxi(i,j,k)=2.*sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)/dx
                                 dutotdxi(i,j,k)=sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)&
                                                 /(log((.5*dx)/z0)*(.5*dx))
                              else
                                 dutotdxi(i,j,k)=-sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)&
                                                 /(log((.5*dx)/z0)*(.5*dx))
                              endif
                           endif
                        endif
                        if(icellflag(i,j,k) .ne. 0 .and. icellflag(i,j-1,k) .ne. 0 &
                              .and. icellflag(i,j+1,k) .ne. 0)then
!mdw 3-08-2008 start gradient changes
                           utot=sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
                           utotm=sqrt(u(i,j-1,k)**2+v(i,j-1,k)**2+w(i,j-1,k)**2)
                           utotp=sqrt(u(i,j+1,k)**2+v(i,j+1,k)**2+w(i,j+1,k)**2)
                           dutotdyc=0.5*(utotp-utotm)/dy
                           dutotdyp=(utotp-utot)/dy
                           dutotdym=(utot-utotm)/dy
                           dutotdya=max(abs(dutotdyp),abs(dutotdym))
                           if(dutotdya.eq.abs(dutotdym))then
                                 dutotdyi(i,j,k)=dutotdym
                              else
                                 dutotdyi(i,j,k)=dutotdyp
                              endif
! mdw 3-08-2004end changes
                        else
                           if(icellflag(i,j,k) .eq. 0)then
                              dutotdyi(i,j,k)=0.
                           else
                              if(icellflag(i,j-1,k) .eq. 0)then
                                 dutotdyi(i,j,k)=sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)&
                                                 /(log((.5*dy)/z0)*(.5*dy))
                              else
                                 dutotdyi(i,j,k)=-sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)&
                                                 /(log((.5*dy)/z0)*(.5*dy))
                              endif
                           endif
                        endif
                     endif
                     call detang(0) ! Calculates the parameters fot the triple rotatation of coord sys.
                     dwall=min(eleff(i,j,k),dxm(i,j,k),dxp(i,j,k),dym(i,j,k),dyp(i,j,k))
                     if(iturbtypeflag.eq.1) local_option(i,j,k)=3
                     elz(i,j,k)=kkar*dwall ! length scale based on distance to wall
                     if(abs(dutotdni(i,j,k)).gt.1.e-6)then
                        x_b=min(dxm(i,j,k),dxp(i,j,k))
                        if(x_b>max(del_b,dx)) x_b=0
                        y_b=min(dym(i,j,k),dyp(i,j,k))
                        if(y_b>max(del_b,dy)) y_b=0
                        dwallg=abs(dutotdyi(i,j,k))*y_b+abs(dutotdxi(i,j,k))*x_b
                        dwallg=dwallg+abs(dutotdzi(i,j,k))*eleff(i,j,k)
                        dwallg=dwallg/dutotdni(i,j,k)
                        elzv=kkar*utot/dutotdni(i,j,k) ! length scale based on distance to null wind
                        if(dwallg*kkar.lt.elzv.and.x_b+y_b.gt.0.) then
! mdw 6-29-2006 changed test so that must be near vertical wall
                           if(iturbtypeflag.eq.1)local_option(i,j,k)=1
                           elz(i,j,k)=kkar*dwallg ! pick the smallest length scale
                        else
! mdw 6-30-2006 changed test so that vortex test does not override normal stuff
                           if(elzv.le.elz(i,j,k))then
                              if(iturbtypeflag.eq.1) local_option(i,j,k)=2
                              elz(i,j,k)=elzv
                           endif
                        endif
                     endif
                     if(rcl.gt.0)then
                        phim=1.+4.7*rcl*eleff(i,j,k)
                        psim=-4.7*rcl*eleff(i,j,k)
                     else
                        phim=(1.-15.*rcl*eleff(i,j,k))**(-.25)
                        psim=2.*alog((1.+1./phim)/2.)+alog((1.+1./phim**2)/2.)-2.*atan(1./phim)+pi/2.
                     endif
                     ustarz(i,j,k)=elz(i,j,k)*dutotdni(i,j,k)/phim ! calculate local ustar
                     ustarz(i,j,k)=max(ustarz(i,j,k),3.e-02)
!mdw 6-23-2004 adjust for vertical structure
                     u3psq=cusq*zbrac*ustarz(i,j,k)**2   ! (u''')^2
                     v3psq=cvsq*zbrac*ustarz(i,j,k)**2   !...
                     w3psq=cwsq*zbrac*ustarz(i,j,k)**2 !...
                     upwp=-ctau13*zbrac*ustarz(i,j,k)**2 ! -tau13
                     upvp=0.
                     vpwp=0.

                     if(ustarz(i,j,k).gt.xloc*ustarg(i,j,k))then
                        if(rcl.lt.0.)then
                           u3psq=u3psq+.6*(ustarz(i,j,k)**2)*(-h*rcl)**(2./3.)
                           v3psq=v3psq+.6*(ustarz(i,j,k)**2)*(-h*rcl)**(2./3.)
                           w3psq=w3psq+3.3*(ustarz(i,j,k)**2)*(-zi(k)*rcl)**(2.3)*(1.-.8*zi(k)/h)**2
                           upwp=upwp*(1.-zi(k)/h)**(.5*rcl*h/(1-rcl*h))
                        endif
                        call rotu3psq ! rotate values back into the orig. grid
                        ufwfi(i,j,k)=ufwf
                        ufvfi(i,j,k)=ufvf
                        vfwfi(i,j,k)=vfwf
                        ustarij(i,j,k)=ustarz(i,j,k)
                        sigui(i,j,k)=sqrt(u3psq)
                        sigvi(i,j,k)=sqrt(v3psq)
                        sigwi(i,j,k)=sqrt(w3psq)
                        upwpi(i,j,k)=upwp
                        
                        if(iturbtypeflag.eq.1)turb_options(i,j,k)=local_option(i,j,k)
                     else ! non-local dominates (possible place for into of TPT)
                        if(iturbtypeflag.eq.1)turb_options(i,j,k)=nonlocal_option(i,j,k)
                        ufsq=ufsqgi(i,j,k)
                        vfsq=vfsqgi(i,j,k)
                        wfsq=wfsqgi(i,j,k)
                        ufvf=ufvfgi(i,j,k)
                        ufwf=ufwfgi(i,j,k)
                        vfwf=vfwfgi(i,j,k)
                        sigui(i,j,k)=sqrt(ufsq)
                        sigvi(i,j,k)=sqrt(vfsq)
                        sigwi(i,j,k)=sqrt(wfsq)
                        ufwfi(i,j,k)=ufwf
                        ufvf=0.
                        ufvfi(i,j,k)=ufvf
                        vfwfi(i,j,k)=vfwf
!mdw 7-25-2005 corrections for axis rotation with non-local mixing
                        ustarij(i,j,k)=ustarg(i,j,k)
                        call rotufsq
                        sigui(i,j,k)=sqrt(u3psq)
                        sigvi(i,j,k)=sqrt(v3psq)
                        sigwi(i,j,k)=sqrt(w3psq)
                        upwpi(i,j,k)=upwp
                     endif
                     sigu=sigui(i,j,k)
                     sigv=sigvi(i,j,k)
                     sigw=sigwi(i,j,k)
                     eps=ustarij(i,j,k)**3*(1.-.75*zi(k)*rcl)*(1.-.85*zi(k)/h)**(1.5)/eleff(i,j,k) ! calculate epsilon for grid cell centers
                  endif
                  epsi(i,j,k)=eps
                  if(format_flag.eq.1.or.format_flag.eq.3)then
                     if(iturbfieldflag.eq.1) write(13,23000)xi(i),yi(j),zi(k),sigu,sigv,sigw,&
                                             elz(i,j,k),eleff(i,j,k),eps,ufvf,ufwf,vfwf
                  endif
                  if(format_flag.eq.2.or.format_flag.eq.3)then
                     if(iturbfieldflag.eq.1) write(63)sigu,sigv,sigw,ufvf,ufwf,vfwf,eps
                  endif
               enddo   lp027
            enddo      lp028
         enddo         lp029
23001    format(3i4,6f8.3)
23000    format(12(1pe11.3))
         deallocate(ufsqgi,vfsqgi,wfsqgi)
         deallocate(ufwfi,ufvfi,vfwfi)
         deallocate(ufvfgi,ufwfgi,vfwfgi)
         deallocate(utotktp,uktop,vktop,wktop)
         if(t .lt. 0.)then !MAN 01/29/2007
            allocate(dsigwdni(nx-1,ny-1,nz),dsigudni(nx-1,ny-1,nz),dsigvdni(nx-1,ny-1,nz))
         endif
lp032:   do k=2,nz-1
lp031:      do j=1,ny-1
lp030:         do i=1,nx-1
                  if(icellflag(i,j,k) .ne. 0)then
                     if(j.lt.2.or.j.ge.ny-2.or.i.lt.2.or.i.ge.nx-2)then
                        dsigwdx=0.
                        dsigwdy=0.
                        dsigudx=0.
                        dsigudy=0.
                        dsigvdx=0.
                        dsigvdy=0.
                        dupwpdx=0.
                        dupwpdy=0.
                     endif
!
! calculate the gradient of sigma w normal to the flow using a CDD
!
                     utot=sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
                     if(dxm(i,j,k).ge.dx.and.dxp(i,j,k).ge.dx.and.i.ne.1.and.i.ne.nx-1)then
                        dsigwdx=.5*(sigwi(i+1,j,k)-sigwi(i-1,j,k))/dx
                        dsigvdx=.5*(sigvi(i+1,j,k)-sigvi(i-1,j,k))/dx
                        dsigudx=.5*(sigui(i+1,j,k)-sigui(i-1,j,k))/dx
                        dupwpdx=.5*(upwpi(i+1,j,k)-upwpi(i-1,j,k))/dx
                     else
                        if(i.eq.1.or.i.eq.nx-1)then
                           dsigwdx=0.
                           dsigvdx=0.
                           dsigudx=0.
                           dupwpdx=0.
                        else
                           if(dxm(i,j,k).lt.dx.and.dxp(i,j,k).gt.dx)then
!mdw 11-21-2005 modified if statements to address particle in 3-walled cells
                              dsigwdni(i,j,k)=(sigwi(i+1,j,k)-sigwi(i,j,k))/dx
                              dsigvdni(i,j,k)=(sigvi(i+1,j,k)-sigvi(i,j,k))/dx
                              dsigudni(i,j,k)=(sigui(i+1,j,k)-sigui(i,j,k))/dx
                              dupwpdni(i,j,k)=(upwpi(i+1,j,k)-upwpi(i,j,k))/dx
                              dsigwdni(i,j,k)=0.
                              dsigvdni(i,j,k)=0.
                              dsigudni(i,j,k)=0.
                              dupwpdni(i,j,k)=0.
                              sigwi(i,j,k)=max(sigwi(i+1,j,k),sigwi(i,j,k))
                              sigvi(i,j,k)=max(sigvi(i+1,j,k),sigvi(i,j,k))
                              sigui(i,j,k)=max(sigui(i+1,j,k),sigui(i,j,k))
                              ustarij(i,j,k)=max(ustarij(i+1,j,k),ustarij(i,j,k))
                              if(abs(upwpi(i,j,k)).lt.abs(upwpi(i+1,j,k)))then
                                 upwpi(i,j,k)=upwpi(i+1,j,k)
                              else
                                 upwpi(i+1,j,k)=upwpi(i,j,k)
                              endif
                           endif
                           if(dxp(i,j,k).lt.dx.and.dxm(i,j,k).gt.dx)then
!mdw 11-21-2005 modified if statements to address particle in 3-walled cells
                              dsigwdni(i,j,k)=(sigwi(i-1,j,k)-sigwi(i,j,k))/dx
                              dsigvdni(i,j,k)=(sigvi(i-1,j,k)-sigvi(i,j,k))/dx
                              dsigudni(i,j,k)=(sigui(i-1,j,k)-sigui(i,j,k))/dx
                              dupwpdni(i,j,k)=(upwpi(i-1,j,k)-upwpi(i,j,k))/dx
                              dsigwdni(i,j,k)=0.
                              dsigvdni(i,j,k)=0.
                              dsigudni(i,j,k)=0.
                              dupwpdni(i,j,k)=0.
                              sigwi(i,j,k)=max(sigwi(i-1,j,k),sigwi(i,j,k))
                              sigvi(i,j,k)=max(sigvi(i-1,j,k),sigvi(i-1,j,k))
                              sigui(i,j,k)=max(sigui(i-1,j,k),sigui(i,j,k))
                              ustarij(i,j,k)=max(ustarij(i-1,j,k),ustarij(i,j,k))
                              if(abs(upwpi(i,j,k)).lt.abs(upwpi(i-1,j,k)))then
                                 upwpi(i,j,k)=upwpi(i-1,j,k)
                              else
                                 upwpi(i-1,j,k)=upwpi(i,j,k)
                              endif
                           endif
                        endif
                     endif
                     if(dym(i,j,k).ge.dy.and.dyp(i,j,k).ge.dy.and.j.ne.1.and.j.ne.ny-1)then
                        dsigwdy=.5*(sigwi(i,j+1,k)-sigwi(i,j-1,k))/dy
                        dsigvdy=.5*(sigvi(i,j+1,k)-sigvi(i,j-1,k))/dy
                        dsigudy=.5*(sigui(i,j+1,k)-sigui(i,j-1,k))/dy
                        dupwpdy=.5*(upwpi(i,j+1,k)-upwpi(i,j-1,k))/dy
                     else
                        if(j.eq.1.or.j.eq.ny-1)then
                           dsigwdy=0.
                        else
                           if(dym(i,j,k).lt.dy.and.dyp(i,j,k).gt.dy)then
!mdw 11-21-2006 modified if statements to address particle in 3-walled cells
                              dsigwdni(i,j,k)=(sigwi(i,j+1,k)-sigwi(i,j,k))/dy
                              dsigvdni(i,j,k)=(sigvi(i,j+1,k)-sigvi(i,j,k))/dy
                              dsigudni(i,j,k)=(sigui(i,j+1,k)-sigui(i,j,k))/dy
                              dupwpdni(i,j,k)=(upwpi(i,j+1,k)-upwpi(i,j,k))/dy
                              dsigwdni(i,j,k)=0.
                              dsigvdni(i,j,k)=0.
                              dsigudni(i,j,k)=0.
                              dupwpdni(i,j,k)=0.
                              sigwi(i,j,k)=max(sigwi(i,j+1,k),sigwi(i,j,k))
                              sigvi(i,j,k)=max(sigvi(i,j+1,k),sigvi(i,j,k))
                              sigui(i,j,k)=max(sigui(i,j+1,k),sigui(i,j,k))
                              ustarij(i,j,k)=max(ustarij(i,j+1,k),ustarij(i,j,k))
                              if(abs(upwpi(i,j,k)).lt.abs(upwpi(i,j+1,k)))then
                                 upwpi(i,j,k)=upwpi(i,j+1,k)
                              else
                                 upwpi(i,j+1,k)=upwpi(i,j,k)
                              endif
                           endif
                           if(dyp(i,j,k).lt.dy.and.dym(i,j,k).gt.dy)then
!mdw 11-21-2005 modified if statements to address particle in 3-walled cells
                              dsigwdni(i,j,k)=(sigwi(i,j-1,k)-sigwi(i,j,k))/dy
                              dsigvdni(i,j,k)=(sigvi(i,j-1,k)-sigvi(i,j,k))/dy
                              dsigudni(i,j,k)=(sigui(i,j-1,k)-sigui(i,j,k))/dy
                              dupwpdni(i,j,k)=(upwpi(i,j-1,k)-upwpi(i,j,k))/dy
                              dsigwdni(i,j,k)=0.
                              dsigvdni(i,j,k)=0.
                              dsigudni(i,j,k)=0.
                              dupwpdni(i,j,k)=0.
                              sigwi(i,j,k)=max(sigwi(i,j-1,k),sigwi(i,j,k))
                              sigvi(i,j,k)=max(sigvi(i,j-1,k),sigvi(i,j,k))
                              sigui(i,j,k)=max(sigui(i,j-1,k),sigui(i,j,k))
                              ustarij(i,j,k)=max(ustarij(i,j-1,k),ustarij(i,j,k))
                              if(abs(upwpi(i,j,k)).lt.abs(upwpi(i,j-1,k)))then
                                 upwpi(i,j,k)=upwpi(i,j-1,k)
                              else
                                 upwpi(i,j-1,k)=upwpi(i,j,k)
                              endif
                           endif
                        endif
                     endif
                     if(dzm(i,j,k).gt.dz.and.k.ne.1.and.k.ne.nz-1.and.dzp(i,j,k).gt.dz)then
                        dsigwdz=.5*(sigwi(i,j,k+1)-sigwi(i,j,k-1))/dz
                        dsigvdz=.5*(sigvi(i,j,k+1)-sigvi(i,j,k-1))/dz
                        dsigudz=.5*(sigui(i,j,k+1)-sigui(i,j,k-1))/dz
                        dupwpdz=.5*(upwpi(i,j,k+1)-upwpi(i,j,k-1))/dz
                     endif
                     if(dzm(i,j,k).le.dz.and.k.ne.1.and.k.ne.nz-1.and.dzp(i,j,k).gt.dz)then
                        dsigwdn=(sigwi(i,j,k+1)-sigwi(i,j,k))/dz
                        dsigvdn=(sigvi(i,j,k+1)-sigvi(i,j,k))/dz
                        dsigudn=(sigui(i,j,k+1)-sigui(i,j,k))/dz
                        dupwpdn=(upwpi(i,j,k+1)-upwpi(i,j,k))/dz
!mdw 9-26-2005 force dsigwdn to be zero near surface
                        dsigwdn=0.
                        dsigudn=0.
                        dsigvdn=0.
                        dupwpdn=0.
                        sigui(i,j,k)=max(sigui(i,j,k+1),sigui(i,j,k))
                        sigvi(i,j,k)=max(sigvi(i,j,k+1),sigvi(i,j,k))
                        sigwi(i,j,k)=max(sigwi(i,j,k+1),sigwi(i,j,k))
                        ustarij(i,j,k)=max(ustarij(i,j,k+1),ustarij(i,j,k))
                        if(abs(upwpi(i,j,k)).lt.abs(upwpi(i,j,k+1)))then
                           upwpi(i,j,k)=upwpi(i,j,k+1)
                        else
                           upwpi(i,j,k+1)=upwpi(i,j,k)
                        endif
                     endif
                     if(dzp(i,j,k).le.dz.and.k.ne.1.and.k.ne.nz-1.and.dzm(i,j,k).gt.dz)then
!mdw 9-26-2005 force dsigwdn to be zero near surface
                        dsigwdn=0.
                        dsigudn=0.
                        dsigvdn=0.
                        dupwpdn=0.
                        sigui(i,j,k)=max(sigui(i,j,k-1),sigui(i,j,k))
                        sigvi(i,j,k)=max(sigvi(i,j,k-1),sigvi(i,j,k))
                        sigwi(i,j,k)=max(sigwi(i,j,k-1),sigwi(i,j,k))
                        ustarij(i,j,k)=max(ustarij(i,j,k-1),ustarij(i,j,k))
                        if(abs(upwpi(i,j,k)).lt.abs(upwpi(i,j,k-1)))then
                           upwpi(i,j,k)=upwpi(i,j,k-1)
                        else
                           upwpi(i,j,k-1)=upwpi(i,j,k)
                        endif
                     endif
                     if((dxm(i,j,k).ge.dx).and.(dxp(i,j,k).ge.dx).and.(dym(i,j,k).ge.dy).and. &
                              (dyp(i,j,k).ge.dy).and.(dzm(i,j,k).ge.dz).and.(dzp(i,j,k).ge.dz))then
                        dsigwdn=ani(i,j,k)*dsigwdx+bni(i,j,k)*dsigwdy+cni(i,j,k)*dsigwdz
                        dsigvdn=ani(i,j,k)*dsigvdx+bni(i,j,k)*dsigvdy+cni(i,j,k)*dsigvdz
                        dsigudn=ani(i,j,k)*dsigudx+bni(i,j,k)*dsigudy+cni(i,j,k)*dsigudz
                        dupwpdn=ani(i,j,k)*dupwpdx+bni(i,j,k)*dupwpdy+cni(i,j,k)*dupwpdz
                     endif
                     dsigwdni(i,j,k)=dsigwdn
                     dsigvdni(i,j,k)=dsigvdn
                     dsigudni(i,j,k)=dsigudn
                     dupwpdni(i,j,k)=dupwpdn
! limiting form for near wall circumstances
                  endif
               enddo   lp030
            enddo   lp031
         enddo   lp032
         if(iturbtypeflag.eq.1)write(101)(((turb_options(i,j,k),i=1,nx-1),j=1,ny-1),k=1,nz-1)
         return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine airtempcal
         use variables
         implicit none
         do itemp=1,nupplevtemp
            if(zi(k).lt.ztlev(itemp))then
               if(itemp.eq.1)then
                  airtemp(k)=uptemp(itemp)
               else
                  airtemp(k)=uptemp(itemp-1)+(zi(k)-ztlev(itemp-1))*(uptemp(itemp)-uptemp(itemp-1))&
                             /(ztlev(itemp)-ztlev(itemp-1))
               endif
               exit
            endif
            if(itemp.eq.nupplevtemp)airtemp(k)=uptemp(itemp)
         enddo
         return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine adjustsigma
         use variables
         implicit none
         real zirhod
         if(wstar.gt.ustar)then
            u3psq=sigu**2+.33*wstar**2
            v3psq=sigv**2+.33*wstar**2
            zirhod=zirhoht-zirholt
            w3psq=sigw**2+1.8*(wstar**2)*(((zi(km)-zirholt)/zirhod)**.666)&
                  *(1.-.8*(zi(km)-zirholt)/zirhod)**2
            upwp=-tau13/(1.-zi(km)/h)**(.5/(rcl/h-1))
            sigu=sqrt(u3psq)
            sigv=sqrt(v3psq)
            sigw=sqrt(w3psq)
            tau13=-upwp
         endif
         if(rclloc.lt.0)then
            u3psq=u3psq+.6*(ustarz(im,jm,km)**2)*(-h*rclloc)**(2./3.)
            v3psq=v3psq+.6*(ustarz(im,jm,km)**2)*(-h*rclloc)**(2./3.)
            w3psq=w3psq+3.3*(ustarz(im,jm,km)**2)*(-zp(m)*rclloc)**(2.3)*(1.-.8*zp(m)/h)**2
            upwp=upwp*(1.-zp(m)/h)**(.5/(rclloc/h-1))
            sigu=sqrt(u3psq)
            sigv=sqrt(v3psq)
            sigw=sqrt(w3psq)
            tau13=-upwp
         endif
         return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine releaseparticles
         use variables
         implicit none
         integer mmp,ipuff,x_idx,y_idx
         real radius,rel_time,seg_time,seg_pos,seg_pos_old,xsl,xsu,ysl,ysu,zsl,zsu
         real fracl,ran2,one_third,erfinv,cosang1,agent_mass,mass_sum,one_nparticles
         real height,xdis,ydis,t_con_dg,vel_eff,buildingradius
         !exfiltration variables
         integer inumtimeex,inumbuildex,itimeex,itimeexstart,itimeexstop,ibuildex
         real timestart,timestop
         integer, allocatable :: buildnumex(:),prr(:)
         real, allocatable :: timeex(:),conin(:,:),dosein(:,:),conex(:,:)
         real, allocatable :: conSPin(:,:),doseSPin(:,:),conSPex(:,:)
         real, allocatable :: buildingvol(:),particlemass(:),massreleased(:)
         real, allocatable :: concentration(:),oldconcentration(:),totalsurfacearea(:)
         real, allocatable :: wallthickness(:),stadiumheight(:),stadiumwidth(:),stadiumlength(:)
         integer particlesperdelta,numreleasingbld
         real minreleased,maxreleased,totalmassreleased
         real cumsurfacearea,ranface,facewidth,thickness,dxy
         
         dxy=0.5*(dx+dy)
         one_third=1./3.
         one_nparticles=1./real(nparticles)
         mass_sum=0.
         if(inextgridflag .eq. 2)then !Outer grid 
            if(format_flag.eq.1)then
               read(33,*)nptot,totsrcmass
               if(ibuoyflag .eq. -1)then
                  do isource=1,nsources
                     read(33,*)hnew(isource),rnew(isource),xscen(isource),yscen(isource)
                  enddo
               endif
            else
               read(38)totsrcmass
               if(ibuoyflag .eq. -1 .or. ibuoyflag_inner .eq. -1)then
                  do isource=1,nsources
                     read(38)hnew(isource),rnew(isource),xscen(isource),yscen(isource)
                  enddo
               endif
            endif
            if(ibuoyflag.eq.-1)then
               allocate(vlocz(mpart),vlocx(mpart),vlocy(mpart))
               allocate(vbulkp(mpart),ubulkp(mpart),wbulkp(mpart))
               allocate(icellflag_dense(nx,ny))
               icellflag_dense(:,:)=0
               vlocz(:)=0.
               vlocx(:)=0.
               vlocy(:)=0.
               vbulkp(:)=0.
               ubulkp(:)=0.
               wbulkp(:)=0.            
            endif
            nptot=nparticles
            xscen=xscen+xoff
            yscen=yscen+yoff
            do mmp=1,nptot
               if(format_flag.eq.1)then
                  read(33,*,end=41111)tstrt(mmp),xps(mmp),yps(mmp),zps(mmp),&
                                      pmass(mmp),weight(mmp),weightd(mmp),dpm(mmp),source_num(mmp)
               else
                  read(38,end=41111)tstrt(mmp),xps(mmp),yps(mmp),zps(mmp),&
                          pmass(mmp),weight(mmp),weightd(mmp),dpm(mmp),source_num(mmp)
               endif
               xps(mmp)=xps(mmp)+xoff
               yps(mmp)=yps(mmp)+yoff
               mmpt=mmp
            enddo
41111       continue
         elseif(isourcetypeflag .eq. 8)then !Exfiltration inner grid
            
            read(25,*)inumtimeex
            read(25,*)inumbuildex
            allocate(timeex(inumtimeex),buildnumex(inumbuildex),conin(inumbuildex,inumtimeex),&
                  conex(inumbuildex,inumtimeex),dosein(inumbuildex,inumtimeex),&
                  conSPin(inumbuildex,inumtimeex),doseSPin(inumbuildex,inumtimeex),&
                  conSPex(inumbuildex,inumtimeex),buildingvol(inumbuildex),&
                  massreleased(inumbuildex),prr(inumbuildex),particlemass(inumbuildex),&
                  concentration(inumbuildex),oldconcentration(inumbuildex),&
                  totalsurfacearea(inumbuildex),wallthickness(inumbuildex),&
                  stadiumheight(inumbuildex),stadiumwidth(inumbuildex),stadiumlength(inumbuildex))
            buildingvol(:)=0.
            massreleased(:)=0.
            prr(:)=0
            particlemass(:)=0.
            concentration(:)=0.
            oldconcentration(:)=0.
            wallthickness(:)=0.
            stadiumheight(:)=0.
            stadiumwidth(:)=0.
            stadiumlength(:)=0.
            do itimeex=1,inumtimeex
               read(25,*)timeex(itimeex)
               read(25,*) !header line
               do ibuildex=1,inumbuildex
                  read(25,*)buildnumex(ibuildex),conin(ibuildex,itimeex),&
                        dosein(ibuildex,itimeex),conex(ibuildex,itimeex),&
                        conSPin(ibuildex,itimeex),doseSPin(ibuildex,itimeex),&
                        conSPex(ibuildex,itimeex)
               enddo
            enddo
            do ibuildex=1,inumbuildex
               select case(bldtype(buildnumex(ibuildex)))
                  case(1,6)
                     buildingvol(ibuildex)=Ht(buildnumex(ibuildex))*&
                           Wti(buildnumex(ibuildex))*Lt(buildnumex(ibuildex))
                     totalsurfacearea(ibuildex)=dy*dz*(aminusx(buildnumex(ibuildex)) &
                           +aplusx(buildnumex(ibuildex))) &
                           +dx*dz*(aminusy(buildnumex(ibuildex)) &
                           +aplusy(buildnumex(ibuildex))) &
                           +dx*dy*(afloor(buildnumex(ibuildex)) &
                           +aroof(buildnumex(ibuildex)))
                  case(2)
                     buildingvol(ibuildex)=Ht(buildnumex(ibuildex))*pi*&
                           0.25*Wti(buildnumex(ibuildex))*Lt(buildnumex(ibuildex))
                     totalsurfacearea(ibuildex)=dy*dz*(aminusx(buildnumex(ibuildex)) &
                           +aplusx(buildnumex(ibuildex))) &
                           +dx*dz*(aminusy(buildnumex(ibuildex)) &
                           +aplusy(buildnumex(ibuildex))) &
                           +dx*dy*(afloor(buildnumex(ibuildex)) &
                           +aroof(buildnumex(ibuildex)))
                  case(3)
                     buildingvol(ibuildex)=Ht(buildnumex(ibuildex))*&
                           0.4993*Wti(buildnumex(ibuildex))*Lt(buildnumex(ibuildex))
                     totalsurfacearea(ibuildex)=dy*dz*(aminusx(buildnumex(ibuildex)) &
                           +aplusx(buildnumex(ibuildex))) &
                           +dx*dz*(aminusy(buildnumex(ibuildex)) &
                           +aplusy(buildnumex(ibuildex))) &
                           +dx*dy*(afloor(buildnumex(ibuildex)) &
                           +aroof(buildnumex(ibuildex)))
                  case(4)
                     if(atten(buildnumex(ibuildex)) .lt. 0)then
                        wallthickness(ibuildex)=abs(atten(buildnumex(ibuildex)))
                        stadiumheight(ibuildex)=0.8*Ht(buildnumex(ibuildex))
                     else
                        wallthickness(ibuildex)=atten(buildnumex(ibuildex))
                        stadiumheight(ibuildex)=Ht(buildnumex(ibuildex))
                     endif
                     stadiumwidth(ibuildex)=Wti(buildnumex(ibuildex))
                     stadiumlength(ibuildex)=Lt(buildnumex(ibuildex))
                     buildingvol(ibuildex)=stadiumheight(ibuildex)*wallthickness(ibuildex)&
                           *(stadiumlength(ibuildex)+stadiumwidth(ibuildex)-(4/3)*wallthickness(ibuildex))
                     totalsurfacearea(ibuildex)=2*stadiumheight(ibuildex)*(stadiumwidth(ibuildex)+stadiumlength(ibuildex)) &
                           +2*sqrt((stadiumheight(ibuildex)**2)+(wallthickness(ibuildex)**2)) &
                           *(stadiumwidth(ibuildex)+stadiumlength(ibuildex)-2*wallthickness(ibuildex))
                  case(5)
                     if(atten(buildnumex(ibuildex)) .lt. 0)then
                        wallthickness(ibuildex)=abs(atten(buildnumex(ibuildex)))
                        stadiumheight(ibuildex)=0.8*Ht(buildnumex(ibuildex))
                     else
                        wallthickness(ibuildex)=atten(buildnumex(ibuildex))
                        stadiumheight(ibuildex)=Ht(buildnumex(ibuildex))
                     endif
                     stadiumwidth(ibuildex)=0.5*Wti(buildnumex(ibuildex))
                     stadiumlength(ibuildex)=0.5*Lt(buildnumex(ibuildex))
                     buildingvol(ibuildex)=stadiumheight(ibuildex)*wallthickness(ibuildex)*pi&
                           *(0.5*(stadiumlength(ibuildex)+stadiumwidth(ibuildex))-wallthickness(ibuildex)/3)
                     totalsurfacearea(ibuildex)=stadiumheight(ibuildex)*pi &
                           *(3*(stadiumlength(ibuildex)+stadiumwidth(ibuildex))-sqrt((3*stadiumlength(ibuildex) &
                           +stadiumwidth(ibuildex))*(stadiumlength(ibuildex)+3*stadiumwidth(ibuildex))))
                     do k=1,int(stadiumheight(ibuildex)/dz)
                        thickness=wallthickness(ibuildex)*(1-((real(k)-0.5)*dz)/stadiumheight(ibuildex))
                        totalsurfacearea(ibuildex)=totalsurfacearea(ibuildex) &
                              +dz*sqrt((stadiumheight(ibuildex)**2)+(wallthickness(ibuildex)**2))*pi &
                              *(3*(stadiumlength(ibuildex)+stadiumwidth(ibuildex)-2*thickness) &
                              -sqrt((3*stadiumlength(ibuildex)+stadiumwidth(ibuildex)-4*thickness) &
                              *(stadiumlength(ibuildex)+3*stadiumwidth(ibuildex)-4*thickness)))
                     enddo
               endselect
            enddo
            do ioutin=1,noutilf
               timestart=real(ioutin-1)*outc+concstrt
               timestop=timestart+outc
               if(ioutin .eq. 1)then
                  itimeexstart=1
               else
                  itimeexstart=itimeexstop+1
               endif
               do itimeex=itimeexstart,inumtimeex
                  if(timeex(itimeex) .le. timestop)then
                     itimeexstop=itimeex
                  elseif(itimeex .eq. inumtimeex .and. ioutin .eq. 1)then
                     itimeexstop=itimeex
                  else
                     exit
                  endif
               enddo
               if(itimeexstop .ge. itimeexstart)then
                  do ibuildex=1,inumbuildex
                     do itimeex=itimeexstart,itimeexstop
                        chiin(buildnumex(ibuildex),ioutin)=chiin(buildnumex(ibuildex),ioutin)+&
                              conin(ibuildex,itimeex)/real(itimeexstop-itimeexstart+1)
                        if(isourcetypeflag .gt. 2 .and. isourcetypeflag .ne. 8)then
                           chiinSP(buildnumex(ibuildex),ioutin)=chiinSP(buildnumex(ibuildex),ioutin)+&
                                 conSPin(ibuildex,itimeex)/real(itimeexstop-itimeexstart+1)
                        endif
                     enddo
                  enddo
               endif
            enddo
            ioutin=0
            particlesperdelta=int(real(nparticles)*delta/dur)
            numreleasingbld=0
            minreleased=1e8
            maxreleased=0
            totalmassreleased=0
            do ibuildex=1,inumbuildex
               if(bldtype(buildnumex(ibuildex)) .lt. 9)then
                  do itimeex=1,inumtimeex
                     massreleased(ibuildex)=massreleased(ibuildex)+conex(ibuildex,itimeex)*buildingvol(ibuildex)
                     if(itimeex .eq. inumtimeex)then
                        massreleased(ibuildex)=massreleased(ibuildex)+conin(ibuildex,itimeex)*buildingvol(ibuildex)
                     endif
                  enddo
                  if(massreleased(ibuildex) .gt. 0.)then
                     numreleasingbld=numreleasingbld+1
                     minreleased=min(minreleased,massreleased(ibuildex))
                     maxreleased=min(maxreleased,massreleased(ibuildex))
                     totalmassreleased=totalmassreleased+massreleased(ibuildex);
                  endif
               endif
            enddo
            do ibuildex=1,inumbuildex
               if(massreleased(ibuildex) .gt. 0.)then
                  prr(ibuildex)=nint(particlesperdelta*massreleased(ibuildex)/totalmassreleased)
                  if(prr(ibuildex) .lt. 1)prr(ibuildex)=1
               endif
            enddo
            nreleases=dur/delta
            itimeex=1
            mp=0
            nptot=0
            do irelease=1,nreleases
               t=delta*real(irelease-1)
               if(t .gt. timeex(itimeex) .and. itimeex .le. inumtimeex)then
                  itimeex=itimeex+1
               endif
               do ibuildex=1,inumbuildex
                  if(itimeex .le. inumtimeex)then
                     particlemass(ibuildex)=conex(ibuildex,itimeex)*buildingvol(ibuildex)/prr(ibuildex)
                     concentration(ibuildex)=conin(ibuildex,itimeex)
                     oldconcentration(ibuildex)=conin(ibuildex,itimeex)
                  else
                     oldconcentration(ibuildex)=concentration(ibuildex)
                     particlemass(ibuildex)=concentration(ibuildex)*exp(-delta*binfex(buildnumex(ibuildex)))&
                           *buildingvol(ibuildex)/prr(ibuildex)
                     concentration(ibuildex)=oldconcentration(ibuildex)*exp(-binf(buildnumex(ibuildex))*delta)
                  endif
                  if(particlemass(ibuildex) .gt. 0)then
                     select case(bldtype(buildnumex(ibuildex)))
                        case(1,6)
                           mmp=0
                           do while(mmp .lt. prr(ibuildex))
                              mp=mp+1
                              mmp=mmp+1
                              ranface=RAN2()
                              cumsurfacearea=dy*dz*aminusx(buildnumex(ibuildex))
                              if(ranface .le. cumsurfacearea/totalsurfacearea(ibuildex))then
                                 iface=1
                              else
                                 cumsurfacearea=cumsurfacearea+dy*dz*aplusx(buildnumex(ibuildex))
                                 if(ranface .le. cumsurfacearea/totalsurfacearea(ibuildex))then
                                    iface=2
                                 else
                                    cumsurfacearea=cumsurfacearea+dx*dz*aminusy(buildnumex(ibuildex))
                                    if(ranface .le. cumsurfacearea/totalsurfacearea(ibuildex))then
                                       iface=3
                                    else
                                       cumsurfacearea=cumsurfacearea+dx*dz*aplusy(buildnumex(ibuildex))
                                       if(ranface .le. cumsurfacearea/totalsurfacearea(ibuildex))then
                                          iface=4
                                       else
                                          cumsurfacearea=cumsurfacearea+dx*dy*afloor(buildnumex(ibuildex))
                                          if(ranface .le. cumsurfacearea/totalsurfacearea(ibuildex))then
                                             iface=5
                                          else
                                             iface=6
                                          endif
                                       endif
                                    endif
                                 endif
                              endif
                              deltax=RAN2()*Lt(buildnumex(ibuildex))
                              deltay=(RAN2()-0.5)*Wti(buildnumex(ibuildex))
                              deltaz=RAN2()*Ht(buildnumex(ibuildex))
                              select case(iface)
                                 case(1)
                                    xsl=xfo(buildnumex(ibuildex))-0.5*dx &
                                          +deltay*sin(gamma(buildnumex(ibuildex))*pi/180.)
                                    ysl=yfo(buildnumex(ibuildex)) &
                                          +deltay*cos(gamma(buildnumex(ibuildex))*pi/180.)
                                    zsl=zfo(buildnumex(ibuildex))+deltaz
                                 case(2)
                                    xsl=xfo(buildnumex(ibuildex))+0.5*dx &
                                          +Lt(buildnumex(ibuildex))*cos(gamma(buildnumex(ibuildex))*pi/180.) &
                                          +deltay*sin(gamma(buildnumex(ibuildex))*pi/180.)
                                    ysl=yfo(buildnumex(ibuildex)) &
                                          +Lt(buildnumex(ibuildex))*sin(gamma(buildnumex(ibuildex))*pi/180.) &
                                          +deltay*cos(gamma(buildnumex(ibuildex))*pi/180.)
                                    zsl=zfo(buildnumex(ibuildex))+deltaz
                                 case(3)
                                    xsl=xfo(buildnumex(ibuildex)) &
                                          +deltax*cos(gamma(buildnumex(ibuildex))*pi/180.) &
                                          +0.5*Wti(buildnumex(ibuildex))*sin(gamma(buildnumex(ibuildex))*pi/180.)
                                    ysl=yfo(buildnumex(ibuildex))-0.5*dy &
                                          +deltax*sin(gamma(buildnumex(ibuildex))*pi/180.) &
                                          +0.5*Wti(buildnumex(ibuildex))*cos(gamma(buildnumex(ibuildex))*pi/180.)
                                    zsl=zfo(buildnumex(ibuildex))+deltaz
                                 case(4)
                                    xsl=xfo(buildnumex(ibuildex)) &
                                          +deltax*cos(gamma(buildnumex(ibuildex))*pi/180.) &
                                          +0.5*Wti(buildnumex(ibuildex))*sin(gamma(buildnumex(ibuildex))*pi/180.)
                                    ysl=yfo(buildnumex(ibuildex))+0.5*dy &
                                          +deltax*sin(gamma(buildnumex(ibuildex))*pi/180.) &
                                          +0.5*Wti(buildnumex(ibuildex))*cos(gamma(buildnumex(ibuildex))*pi/180.)
                                    zsl=zfo(buildnumex(ibuildex))+deltaz
                                 case(5)
                                    xsl=xfo(buildnumex(ibuildex)) &
                                          +deltax*cos(gamma(buildnumex(ibuildex))*pi/180.) &
                                          +deltay*sin(gamma(buildnumex(ibuildex))*pi/180.)
                                    ysl=yfo(buildnumex(ibuildex)) &
                                          +deltax*sin(gamma(buildnumex(ibuildex))*pi/180.) &
                                          +deltay*cos(gamma(buildnumex(ibuildex))*pi/180.)
                                    zsl=zfo(buildnumex(ibuildex))+Ht(buildnumex(ibuildex))+0.5*dz
                                 case(6)
                                    xsl=xfo(buildnumex(ibuildex)) &
                                          +deltax*cos(gamma(buildnumex(ibuildex))*pi/180.) &
                                          +deltay*sin(gamma(buildnumex(ibuildex))*pi/180.)
                                    ysl=yfo(buildnumex(ibuildex)) &
                                          +deltax*sin(gamma(buildnumex(ibuildex))*pi/180.) &
                                          +deltay*cos(gamma(buildnumex(ibuildex))*pi/180.)
                                    zsl=zfo(buildnumex(ibuildex))-0.5*dz
                              endselect
                              tstrt(mp)=t
                              xps(mp)=xsl
                              yps(mp)=ysl
                              zps(mp)=zsl
                              weight(mp)=1.
                              weightd(mp)=1.
                              source_num(mp)=real(buildnumex(ibuildex))
                              pmass(mp)=particlemass(ibuildex)
                              im=int(xps(mp)/dx)+1
                              jm=int(yps(mp)/dy)+1
                              km=int((zps(mp)+dz)/dz)+1
                              if(icellflag(im,jm,km) .eq. 0)then
                                 mp=mp-1
                                 mmp=mmp-1
                              endif
                           enddo
                        case(2)
                           mmp=0
                           do while(mmp .lt. prr(ibuildex))
                              mp=mp+1
                              mmp=mmp+1
                              cumsurfacearea=dy*dz*(aminusx(buildnumex(ibuildex)) &
                                    +aplusx(buildnumex(ibuildex))) &
                                    +dx*dz*(aminusy(buildnumex(ibuildex)) &
                                    +aplusy(buildnumex(ibuildex)))
                              ranface=RAN2()
                              if(ranface .le. cumsurfacearea/totalsurfacearea(ibuildex))then
                                 ang1=2*pi*RAN2()
                                 radius=buildingradius(0.5*(dx+Lt(buildnumex(ibuildex))),&
                                       0.5*(dy+Wti(buildnumex(ibuildex))),&
                                       gamma(buildnumex(ibuildex)),ang1)
                                 deltaz=RAN2()*Ht(buildnumex(ibuildex))
                              else
                                 cumsurfacearea=cumsurfacearea+dx*dy*aroof(buildnumex(ibuildex))
                                 ang1=2*pi*RAN2()
                                 radius=(RAN2()**0.5)*buildingradius(0.5*Lt(buildnumex(ibuildex)),&
                                       0.5*Wti(buildnumex(ibuildex)),&
                                       gamma(buildnumex(ibuildex)),ang1)
                                 if(ranface .le. cumsurfacearea/totalsurfacearea(ibuildex))then
                                    deltaz=Ht(buildnumex(ibuildex))+0.5*dz
                                 else
                                    deltaz=-0.5*dz
                                 endif
                              endif
                              tstrt(mp)=t
                              weight(mp)=1.
                              weightd(mp)=1.
                              source_num(mp)=real(buildnumex(ibuildex))
                              pmass(mp)=particlemass(ibuildex)
                              xps(mp)=xfo(buildnumex(ibuildex))+0.5*Lt(buildnumex(ibuildex))&
                                    *cos(gamma(buildnumex(ibuildex))*pi/180.)&
                                    +radius*cos(ang1)
                              yps(mp)=yfo(buildnumex(ibuildex))+0.5*Lt(buildnumex(ibuildex))&
                                    *sin(gamma(buildnumex(ibuildex))*pi/180.)&
                                    +radius*sin(ang1)
                              zps(mp)=zfo(buildnumex(ibuildex))+deltaz
                              im=int(xps(mp)/dx)+1
                              jm=int(yps(mp)/dy)+1
                              km=int((zps(mp)+dz)/dz)+1
                              if(icellflag(im,jm,km) .eq. 0)then
                                 mp=mp-1
                                 mmp=mmp-1
                              endif
                           enddo
                        case(3)
                           mmp=0
                           do while(mmp .lt. prr(ibuildex))
                              mp=mp+1
                              mmp=mmp+1
                              cumsurfacearea=dy*dz*(aminusx(buildnumex(ibuildex)) &
                                    +aplusx(buildnumex(ibuildex))) &
                                    +dx*dz*(aminusy(buildnumex(ibuildex)) &
                                    +aplusy(buildnumex(ibuildex)))
                              ranface=RAN2()
                              if(ranface .le. cumsurfacearea/totalsurfacearea(ibuildex))then
                                 ranface=10*RAN2()
                                 if(ranface .eq. 10)ranface=0.
                                 if(ranface .lt. 5)then
                                    facewidth=(RAN2()-0.5)*Lt(buildnumex(ibuildex))*sin(0.2*pi)
                                    radius=0.5*Lt(buildnumex(ibuildex))*cos(0.2*pi)+0.25*(dx+dy)
                                 else
                                    facewidth=(RAN2()-0.5)*0.4*Lt(buildnumex(ibuildex))*sin(0.2*pi)
                                    radius=0.2*Lt(buildnumex(ibuildex))*cos(0.2*pi)+0.25*(dx+dy)
                                 endif
                                 select case(int(ranface))
                                    case(0,5)
                                       ang1=(54-gamma(buildnumex(ibuildex)))*pi/180
                                       ang2=(144-gamma(buildnumex(ibuildex)))*pi/180
                                    case(1,6)
                                       ang1=(-18-gamma(buildnumex(ibuildex)))*pi/180
                                       ang2=(72-gamma(buildnumex(ibuildex)))*pi/180
                                    case(2,7)
                                       ang1=(-90-gamma(buildnumex(ibuildex)))*pi/180
                                       ang2=(-gamma(buildnumex(ibuildex)))*pi/180
                                    case(3,8)
                                       ang1=(-162-gamma(buildnumex(ibuildex)))*pi/180
                                       ang2=(-72-gamma(buildnumex(ibuildex)))*pi/180
                                    case(4,9)
                                       ang1=(126-gamma(buildnumex(ibuildex)))*pi/180
                                       ang2=(36-gamma(buildnumex(ibuildex)))*pi/180
                                 endselect
                                 deltax=radius*cos(ang1)+facewidth*cos(ang2)
                                 deltay=radius*sin(ang1)+facewidth*sin(ang2)
                                 deltaz=RAN2()*Ht(buildnumex(ibuildex))
                              else
                                 cumsurfacearea=cumsurfacearea+dx*dy*aroof(buildnumex(ibuildex))
                                 radius=0.5*Lt(buildnumex(ibuildex))*cos(0.2*pi)&
                                       *(0.4+0.6*(RAN2()**0.5))
                                 ang1=2*pi*RAN2()
                                 deltax=radius*cos(ang1)
                                 deltay=radius*sin(ang1)
                                 if(ranface .le. cumsurfacearea/totalsurfacearea(ibuildex))then
                                    deltaz=Ht(buildnumex(ibuildex))+0.5*dz
                                 else
                                    deltaz=-0.5*dz
                                 endif
                              endif
                              tstrt(mp)=t
                              weight(mp)=1.
                              weightd(mp)=1.
                              source_num(mp)=real(buildnumex(ibuildex))
                              pmass(mp)=particlemass(ibuildex)
                              xps(mp)=xfo(buildnumex(ibuildex))+deltax
                              yps(mp)=yfo(buildnumex(ibuildex))+deltay
                              zps(mp)=zfo(buildnumex(ibuildex))+deltaz
                              im=int(xps(mp)/dx)+1
                              jm=int(yps(mp)/dy)+1
                              km=int((zps(mp)+dz)/dz)+1
                              if(icellflag(im,jm,km) .eq. 0)then
                                 mp=mp-1
                                 mmp=mmp-1
                              endif
                           enddo
                        case(4)
                           mmp=0
                           do while(mmp .lt. prr(ibuildex))
                              mp=mp+1
                              mmp=mmp+1
                              cumsurfacearea=2*stadiumheight(ibuildex)*stadiumlength(ibuildex)
                              ranface=RAN2()
                              deltaz=stadiumheight(ibuildex)*RAN2()
                              if(ranface .le. cumsurfacearea/totalsurfacearea(ibuildex))then
                                 ranface=RAN2()
                                 facewidth=0.5*stadiumlength(ibuildex)*(RAN2()-0.5)
                                 if(ranface .gt. 0.5)then
                                    deltax=0.5*(stadiumwidth(ibuildex)+dy)*sin(gamma(buildnumex(ibuildex))*pi/180) &
                                          +facewidth*cos(gamma(buildnumex(ibuildex))*pi/180)
                                    deltay=0.5*(stadiumwidth(ibuildex)+dy)*cos(gamma(buildnumex(ibuildex))*pi/180) &
                                          +facewidth*sin(gamma(buildnumex(ibuildex))*pi/180)
                                 else
                                    deltax=-0.5*(stadiumwidth(ibuildex)+dy)*sin(gamma(buildnumex(ibuildex))*pi/180) &
                                          +facewidth*cos(gamma(buildnumex(ibuildex))*pi/180)
                                    deltay=-0.5*(stadiumwidth(ibuildex)+dy)*cos(gamma(buildnumex(ibuildex))*pi/180) &
                                          +facewidth*sin(gamma(buildnumex(ibuildex))*pi/180)
                                 endif
                              else
                                 cumsurfacearea=cumsurfacearea+2*stadiumheight(ibuildex)*stadiumwidth(ibuildex)
                                 if(ranface .le. cumsurfacearea/totalsurfacearea(ibuildex))then
                                    ranface=RAN2()
                                    facewidth=0.5*stadiumwidth(ibuildex)*(RAN2()-0.5)
                                    if(ranface .gt. 0.5)then
                                       deltax=0.5*(stadiumlength(ibuildex)+dx)*cos(gamma(buildnumex(ibuildex))*pi/180) &
                                             +facewidth*sin(gamma(buildnumex(ibuildex))*pi/180)
                                       deltay=0.5*(stadiumlength(ibuildex)+dx)*sin(gamma(buildnumex(ibuildex))*pi/180) &
                                             +facewidth*cos(gamma(buildnumex(ibuildex))*pi/180)
                                    else
                                       deltax=-0.5*(stadiumlength(ibuildex)+dx)*cos(gamma(buildnumex(ibuildex))*pi/180) &
                                             +facewidth*sin(gamma(buildnumex(ibuildex))*pi/180)
                                       deltay=-0.5*(stadiumlength(ibuildex)+dx)*sin(gamma(buildnumex(ibuildex))*pi/180) &
                                             +facewidth*cos(gamma(buildnumex(ibuildex))*pi/180)
                                    endif
                                 else
                                    thickness=wallthickness(ibuildex)*(1-deltaz/stadiumheight(ibuildex))
                                    cumsurfacearea=2*sqrt((stadiumheight(ibuildex)**2)+(wallthickness(ibuildex)**2)) &
                                          *(stadiumlength(ibuildex)-wallthickness(ibuildex))
                                    if(ranface .le. cumsurfacearea/totalsurfacearea(ibuildex))then
                                       ranface=RAN2()
                                       facewidth=(0.5*stadiumlength(ibuildex)-thickness)*(RAN2()-0.5)
                                       if(ranface .gt. 0.5)then
                                          deltax=0.5*(stadiumwidth(ibuildex)-dy-thickness)&
                                                *sin(gamma(buildnumex(ibuildex))*pi/180) &
                                                +facewidth*cos(gamma(buildnumex(ibuildex))*pi/180)
                                          deltay=0.5*(stadiumwidth(ibuildex)-dy-thickness)&
                                                *cos(gamma(buildnumex(ibuildex))*pi/180) &
                                                +facewidth*sin(gamma(buildnumex(ibuildex))*pi/180)
                                       else
                                          deltax=-0.5*(stadiumwidth(ibuildex)-dy-thickness)&
                                                *sin(gamma(buildnumex(ibuildex))*pi/180) &
                                                +facewidth*cos(gamma(buildnumex(ibuildex))*pi/180)
                                          deltay=-0.5*(stadiumwidth(ibuildex)-dy-thickness)&
                                                *cos(gamma(buildnumex(ibuildex))*pi/180) &
                                                +facewidth*sin(gamma(buildnumex(ibuildex))*pi/180)
                                       endif
                                    else
                                       ranface=RAN2()
                                       facewidth=(0.5*stadiumwidth(ibuildex)-thickness)*(RAN2()-0.5)
                                       if(ranface .gt. 0.5)then
                                          deltax=0.5*(stadiumlength(ibuildex)-dx-thickness)&
                                                *cos(gamma(buildnumex(ibuildex))*pi/180) &
                                                +facewidth*sin(gamma(buildnumex(ibuildex))*pi/180)
                                          deltay=0.5*(stadiumlength(ibuildex)-dx-thickness)&
                                                *sin(gamma(buildnumex(ibuildex))*pi/180) &
                                                +facewidth*cos(gamma(buildnumex(ibuildex))*pi/180)
                                       else
                                          deltax=-0.5*(stadiumlength(ibuildex)-dx-thickness)&
                                                *cos(gamma(buildnumex(ibuildex))*pi/180) &
                                                +facewidth*sin(gamma(buildnumex(ibuildex))*pi/180)
                                          deltay=-0.5*(stadiumlength(ibuildex)-dx-thickness)&
                                                *sin(gamma(buildnumex(ibuildex))*pi/180) &
                                                +facewidth*cos(gamma(buildnumex(ibuildex))*pi/180)
                                       endif
                                    endif
                                 endif
                              endif
                              tstrt(mp)=t
                              weight(mp)=1.
                              weightd(mp)=1.
                              source_num(mp)=real(buildnumex(ibuildex))
                              pmass(mp)=particlemass(ibuildex)
                              xps(mp)=xfo(buildnumex(ibuildex))+0.5*stadiumlength(ibuildex)&
                                    *cos(gamma(buildnumex(ibuildex))*pi/180)+deltax
                              yps(mp)=yfo(buildnumex(ibuildex))+0.5*stadiumlength(ibuildex)&
                                    *sin(gamma(buildnumex(ibuildex))*pi/180)+deltay
                              zps(mp)=zfo(buildnumex(ibuildex))+deltaz
                              im=int(xps(mp)/dx)+1
                              jm=int(yps(mp)/dy)+1
                              km=int((zps(mp)+dz)/dz)+1
                              if(icellflag(im,jm,km) .eq. 0)then
                                 mp=mp-1
                                 mmp=mmp-1
                              endif
                           enddo
                        case(5)
                           mmp=0
                           do while(mmp .lt. prr(ibuildex))
                              mp=mp+1
                              mmp=mmp+1
                              cumsurfacearea=stadiumheight(ibuildex)*pi &
                                    *(3*(stadiumlength(ibuildex)+stadiumwidth(ibuildex))&
                                    -sqrt((3*stadiumlength(ibuildex)+stadiumwidth(ibuildex))&
                                    *(stadiumlength(ibuildex)+3*stadiumwidth(ibuildex))))
                              ranface=RAN2()
                              deltaz=stadiumheight(ibuildex)*RAN2()
                              ang1=2*pi*RAN2()
                              if(ranface .le. cumsurfacearea/totalsurfacearea(ibuildex))then
                                 radius=buildingradius(stadiumlength(ibuildex)+dxy,&
                                       stadiumwidth(ibuildex)+dxy,&
                                       gamma(buildnumex(ibuildex)),ang1)
                              else
                                 thickness=wallthickness(ibuildex)*(1-deltaz/stadiumheight(ibuildex))
                                 radius=buildingradius(stadiumlength(ibuildex)-dxy-thickness,&
                                       stadiumwidth(ibuildex)-dxy-thickness,&
                                       gamma(buildnumex(ibuildex)),ang1)
                              endif
                              tstrt(mp)=t
                              weight(mp)=1.
                              weightd(mp)=1.
                              source_num(mp)=real(buildnumex(ibuildex))
                              pmass(mp)=particlemass(ibuildex)
                              xps(mp)=xfo(buildnumex(ibuildex))+stadiumlength(ibuildex)&
                                    *cos(gamma(buildnumex(ibuildex))*pi/180.)+radius*cos(ang1)
                              yps(mp)=yfo(buildnumex(ibuildex))+stadiumlength(ibuildex)&
                                    *sin(gamma(buildnumex(ibuildex))*pi/180.)+radius*sin(ang1)
                              zps(mp)=zfo(buildnumex(ibuildex))+deltaz
                              im=int(xps(mp)/dx)+1
                              jm=int(yps(mp)/dy)+1
                              km=int((zps(mp)+dz)/dz)+1
                              if(icellflag(im,jm,km) .eq. 0)then
                                 mp=mp-1
                                 mmp=mmp-1
                              endif
                           enddo
                     endselect
                     nptot=nptot+prr(ibuildex)
                  endif
               enddo
            enddo
            t=-delta
            deallocate(timeex,buildnumex,conin,conex,dosein,conSPin,doseSPin,&
                  conSPex,buildingvol,massreleased,prr,particlemass,&
                  concentration,oldconcentration,totalsurfacearea,&
                  wallthickness,stadiumheight,stadiumwidth,stadiumlength)
         else !Inner grid
            mp=0
            nptot=0
            do isource=1,nsources
               if(source_release(isource).eq.1)then
                  nreleases=1
               else
                  nreleases=nint(source_dur(isource)/delta)
               endif
               select case(source_geometry(isource))
                  case(0) !Array of spherical shell sources for sensor siting mode
                     xpr=xgl+dxg*real(ixg-1)
                     ypr=ygl+dyg*real(iyg-1)
                     rpr=source_r(1)
! here we force source to be above terrain by zpr
                     ipr=int(xpr/dx)+1
                     jpr=int(ypr/dy)+1
                     height=hgt(ipr,jpr)
! check for buildings within the radius of the source                     
                     y_idx=int((ypr-(rpr+z0))/dy)+1
                     if(y_idx .gt. 0 .and. y_idx .lt. ny)then
                        if (hgt(ipr,y_idx) .gt. height)then
                           height=hgt(ipr,y_idx)
                        endif
                     endif
                     y_idx=int((ypr+(rpr+z0))/dy)+1
                     if(y_idx .gt. 0 .and. y_idx .lt. ny)then
                        if (hgt(ipr,y_idx) .gt. height)then
                           height=hgt(ipr,y_idx)
                        endif
                     endif
                     x_idx=int((xpr-(rpr+z0))/dx)+1
                     if(x_idx .gt. 0 .and. x_idx .lt. nx)then
                        if (hgt(x_idx,jpr) .gt. height)then
                           height=hgt(x_idx,jpr)
                        endif
                     endif
                     x_idx=int((xpr+(rpr+z0))/dx)+1
                     if(x_idx .gt. 0 .and. x_idx .lt. nx)then
                        if (hgt(x_idx,jpr) .gt. height)then
                           height=hgt(x_idx,jpr)
                        endif
                     endif
                     ! check diagonal directions
                     !southwest corner
                     xdis=xpr-(real(ipr-1)*dx+z0)
                     ydis=ypr-(real(jpr-1)*dy+z0)
                     if(rpr .gt. sqrt((xdis**2)+(ydis**2)))then
                        if(ipr .gt. 1 .and. jpr .gt. 1)then
                           if(hgt(ipr-1,jpr-1) .gt. height)then
                              height=hgt(ipr-1,jpr-1)
                           endif
                        endif
                     endif
                     !northwest corner
                     xdis=xpr-(real(ipr-1)*dx+z0)
                     ydis=ypr-(real(jpr)*dy-z0)
                     if(rpr .gt. sqrt((xdis**2)+(ydis**2)))then
                        if(ipr .gt. 1 .and. jpr .lt. ny-1)then
                           if(hgt(ipr-1,jpr+1) .gt. height)then
                              height=hgt(ipr-1,jpr+1)
                           endif
                        endif
                     endif
                     !northeast corner
                     xdis=xpr-(real(ipr)*dx-z0)
                     ydis=ypr-(real(jpr)*dy-z0)
                     if(rpr .gt. sqrt((xdis**2)+(ydis**2)))then
                        if(ipr .lt. nx-1 .and. jpr .lt. ny-1)then
                           if(hgt(ipr+1,jpr+1) .gt. height)then
                              height=hgt(ipr+1,jpr+1)
                           endif
                        endif
                     endif
                     !southeast corner
                     xdis=xpr-(real(ipr)*dx-z0)
                     ydis=ypr-(real(jpr-1)*dy+z0)
                     if(rpr .gt. sqrt((xdis**2)+(ydis**2)))then
                        if(ipr .lt. nx-1 .and. jpr .gt. 1)then
                           if(hgt(ipr+1,jpr-1) .gt. height)then
                              height=hgt(ipr+1,jpr-1)
                           endif
                        endif
                     endif
                     zpr=source_z(1)+height
                     kpr=int((zpr+dz)/dz)+1
                     nps=1
                     npf=nps+source_prr(isource)-1
                     write(16,*)xpr,ypr,zpr
                     do mmp=nps,npf
                        mp=mp+1
                        tstrt(mp)=0.
                        weight(mp)=1.
                        weightd(mp)=1.
                        ang2=2.*pi*RAN2()
                        ang1=pi*RAN2()
                        deltax=rpr*sin(ang1)*cos(ang2)
                        deltay=rpr*sin(ang1)*sin(ang2)
                        deltaz=rpr*cos(ang1)
                        xps(mp)=xpr+deltax
                        yps(mp)=ypr+deltay
                        zps(mp)=zpr+deltaz
                        source_num(mp)=real(isource)
                        pmass(mp)=source_weight(isource)*totsrcmass*one_nparticles
                     enddo
                     nps=npf
                     nptot=nptot+source_prr(isource)
                  case(1) !Spherical shell source
                     inode=source_startidx(isource)
                     do irelease=1,nreleases
                        do mmp=1,source_prr(isource)
                           mp=mp+1
                           tstrt(mp)=source_start(isource)+delta*real(irelease-1)
                           ang2=2.*pi*RAN2()
                           cosang1=2.*(RAN2()-.5)
                           ang1=acos(cosang1)
                           deltax=source_r(inode)*sin(ang1)*cos(ang2)
                           deltay=source_r(inode)*sin(ang1)*sin(ang2)
                           deltaz=source_r(inode)*cos(ang1)
                           xps(mp)=source_x(inode)+deltax
                           yps(mp)=source_y(inode)+deltay
                           zps(mp)=source_z(inode)+deltaz
                           weight(mp)=1.
                           weightd(mp)=1.
                           source_num(mp)=real(isource)
                           pmass(mp)=source_weight(isource)*totsrcmass*one_nparticles
                        enddo
                        nptot=nptot+source_prr(isource)
                     enddo
                  case(7) !Spherical volume source
                     inode=source_startidx(isource)
                     do irelease=1,nreleases
                        do mmp=1,source_prr(isource)
                           mp=mp+1
                           tstrt(mp)=source_start(isource)+delta*real(irelease-1)
                           ang2=2.*pi*RAN2()
                           cosang1=2.*(RAN2()-.5)
                           ang1=acos(cosang1)
                           radius=source_r(inode)*(RAN2()**one_third)
                           deltax=radius*sin(ang1)*cos(ang2)
                           deltay=radius*sin(ang1)*sin(ang2)
                           deltaz=radius*cos(ang1)
                           xps(mp)=source_x(inode)+deltax
                           yps(mp)=source_y(inode)+deltay
                           zps(mp)=source_z(inode)+deltaz
                           weight(mp)=1.
                           weightd(mp)=1.
                           source_num(mp)=real(isource)
                           pmass(mp)=source_weight(isource)*totsrcmass*one_nparticles
                        enddo
                        nptot=nptot+source_prr(isource)
                     enddo
                  case(2) !Line Source
                     do inode=source_startidx(isource),source_stopidx(isource)-1
                        xsl=source_x(inode)
                        xsu=source_x(inode+1)
                        ysl=source_y(inode)
                        ysu=source_y(inode+1)
                        zsl=source_z(inode)
                        zsu=source_z(inode+1)
                        xpc=.5*(xsu+xsl)
                        ypc=.5*(ysu+ysl)
                        zpc=.5*(zsu+zsl)
                        do irelease=1,nreleases
                           do mmp=1,source_prr(isource)
                              mp=mp+1
                              tstrt(mp)=source_start(isource)+delta*real(irelease-1)
                              fracl=RAN2()
                              deltax=xsl+fracl*(xsu-xsl)-xpc
                              deltay=ysl+fracl*(ysu-ysl)-ypc
                              deltaz=zsl+fracl*(zsu-zsl)-zpc
                              xps(mp)=xpc+deltax
                              yps(mp)=ypc+deltay
                              zps(mp)=zpc+deltaz
                              weight(mp)=1.
                              weightd(mp)=1.
                              source_num(mp)=real(isource)
                              pmass(mp)=source_weight(isource)*totsrcmass*one_nparticles
                           enddo
                           nptot=nptot+source_prr(isource)
                        enddo
                     enddo
                  case(3) !Vertical Cylinder source
                     inode=source_startidx(isource)
                     if(ibuoyflag.ge.0)then
                        do irelease=1,nreleases
                           do mmp=1,source_prr(isource)
                              mp=mp+1
                              tstrt(mp)=source_start(isource)+delta*real(irelease-1)
                              ang1=2.*pi*RAN2()
                              radius=source_r(inode)*(RAN2()**0.5)
                              deltax=radius*cos(ang1)
                              deltay=radius*sin(ang1)
                              deltaz=source_h(inode)*RAN2()
                              xps(mp)=source_x(inode)+deltax
                              yps(mp)=source_y(inode)+deltay
                              zps(mp)=source_z(inode)+deltaz
                              weight(mp)=1.
                              weightd(mp)=1.
                              source_num(mp)=real(isource)
                              pmass(mp)=source_weight(isource)*totsrcmass*one_nparticles
                           enddo
                           nptot=nptot+source_prr(isource)
                        enddo
                     else
! Dense Gas source type 
                        if(nreleases.eq.1)then
!                           volsource=.001*totsrcmass/source_density(isource)
                           select case(source_str_units(isource))
                              case(1)
                                 volsource=.001*source_strength(isource)/source_density(isource)
                              case(3)
                                 volsource=.001*source_strength(isource)
                           endselect
                           source_r(inode)=sqrt((volsource/source_h(inode))/pi)
                           vout=0.
                        else
                           select case(source_str_units(isource))
                              case(1)
                                 volsource=.001*source_strength(isource)*delta/(source_density(isource)*source_dur(isource))
                              case(2)
                                 volsource=.001*source_strength(isource)*delta/source_density(isource)
                              case(3)
                                 volsource=.001*source_strength(isource)*delta/source_dur(isource)
                              case(4)
                                 volsource=.001*source_strength(isource)*delta
                           endselect
                           vup=volsource/(pi*delta*source_r(inode)**2)
                           vout=vup*source_r(inode)**2/(source_r(inode)*sqrt(source_h(inode)**2+source_r(inode)**2))
                           height=0.
                           im=int(source_x(inode)/dx)+1
                           jm=int(source_y(inode)/dy)+1
                           do ipuff=1,10
                              km=min(int((source_z(inode)+height+dz)/dz)+1,nz-1)
!                              print*,im,jm,km
                              vel_eff=sqrt((u(im,jm,km)**2.)+(v(im,jm,km)**2.))
                              if(km .le. 2)then
                                 vel_eff=vel_eff*(log((source_z(inode)+height)/z0)/log(dz/z0))
                              endif
                              t_con_dg=2*source_r(inode)/vel_eff
                              height=vup*t_con_dg
                           enddo
                           source_h(inode)=height
                        endif
                        if(isource.eq.1)then
                           allocate(vlocz(mpart),vlocx(mpart),vlocy(mpart))
                           allocate(vbulkp(mpart),ubulkp(mpart),wbulkp(mpart))
                           allocate(icellflag_dense(nx,ny))
                           icellflag_dense(:,:)=0
                           vlocz(:)=0.
                           vlocx(:)=0.
                           vlocy(:)=0.
                           vbulkp(:)=0.
                           ubulkp(:)=0.
                           wbulkp(:)=0.
                        endif
                        ipc=int(source_x(inode)/dx)+1
                        jpc=int(source_y(inode)/dy)+1
                        kpc=int((source_h(inode)+source_z(inode)+dz)/dz)+1
                        utoth=sqrt(u(ipc,jpc,kpc)**2+v(ipc,jpc,kpc)**2)
                        deltaup=source_r(inode)/utoth
                        rnew(isource)=source_r(inode)
                        reff(isource)=source_r(inode)
                        hnew(isource)=source_h(inode)
                        delri=.0005*min(dx,dy)            
                        kcellst=int((hnew(isource)+.5*dz+dz)/dz)+1
                        volnew(isource)=pi*hnew(isource)*source_r(inode)**2
                        vol0(isource)=volnew(isource)
                        iconarea=0
                        rnewp(isource)=rnew(isource)
                        do while(iconarea.lt.20)
                           uxsum=0.
                           vysum=0.
                           cellsum=0.
                           uistarsum=0.
                           utotsum=0.
                           slopexsum=0.
                           slopeysum=0.
                           addarea=0.
                           xscen(isource)=source_x(inode)
                           yscen(isource)=source_y(inode)
                           icells=int((source_x(inode))/dx)+1
                           jcells=int((source_y(inode))/dy)+1
                           istrt=int((source_x(inode)-reff(isource)-dx)/dx)+1
                           ifin=int((source_x(inode)+reff(isource)+dx)/dx)+1
                           jstrt=int((source_y(inode)-reff(isource)-dy)/dy)+1
                           jfin=int((source_y(inode)+reff(isource)+dy)/dy)+1
                           do icell=istrt,ifin
                              xcellcn=dx*real(icell-1)+.5*dx
                              do jcell=jstrt,jfin
                                 ycellcn=dy*real(jcell-1)+.5*dy
                                 if(icells.ge.icell)then
                                    xcellcl=dx*real(icell)-delri
                                 else
                                    xcellcl=dx*real(icell-1)-delri
                                 endif
                                 if(jcells.ge.jcell)then
                                    ycellcl=dy*real(jcell)-delri
                                 else
                                    ycellcl=dy*real(jcell-1)-delri
                                 endif
                                 rcellcl=sqrt((xcellcl-xscen(isource))**2+(ycellcl-yscen(isource))**2)
                                 if(rcellcl.le.reff(isource).or.(icells.eq.icell.and.jcells.eq.jcell))then
                                    if(icellflag(icell,jcell,2) .ne. 0 .and.(icellflag_dense(icell,jcell) .eq. 0 &
                                          .or. icellflag_dense(icell,jcell) .eq. isource))then
                                       rhocell(icell,jcell,2)=(vol0(isource)/volnew(isource))*source_density(isource)&
                                                              +(1.-vol0(isource)/volnew(isource))*rhoair
                                       uxsum=u(icell,jcell,kcellst)+uxsum
                                       vysum=v(icell,jcell,kcellst)+vysum
                                       uistarsum=ustarij(icell,jcell,kcellst)+uistarsum
                                       cellsum=cellsum+1.
                                       call densegas_spprofile
                                       slopexsum=slopexsum+(elevij(icell+2,jcell+1)-elevij(icell,jcell+1))/(2.*dx)
                                       slopeysum=slopeysum+(elevij(icell+1,jcell+2)-elevij(icell+1,jcell))/(2.*dy)
                                       icellflag_dense(icell,jcell)=isource
                                    else
                                       rcellcn=sqrt((xcellcn-xscen(isource))**2+(ycellcn-yscen(isource))**2)
                                       if(rcellcn.le.reff(isource))then
                                          addarea=addarea+dx*dy
                                       endif
                                    endif
                                 endif
                              enddo
                           enddo
                           effarea=pi*rnew(isource)**2+addarea
                           curarea=pi*reff(isource)**2
                           if(curarea.lt..99*effarea)then
                              reff(isource)=sqrt(effarea/pi)
                              iconarea=iconarea+1
                           else
                              exit
                           endif
                        enddo
                        xscen(isource)=source_x(isource)
                        yscen(isource)=source_y(isource)
                        ustars=uistarsum/cellsum
                        utotdg=utotsum/cellsum
                        cosphidg(isource)=uxsum/sqrt(uxsum**2+vysum**2)
                        sinphidg(isource)=vysum/sqrt(uxsum**2+vysum**2)
                        hnew(isource)=source_h(isource)
                        rnew(isource)=source_r(isource)
                        slopex=slopexsum/cellsum
                        slopey=slopeysum/cellsum
                        do irelease=1,nreleases
                           mmp=0
                           do while(mmp .lt. source_prr(isource))
                              mp=mp+1
                              mmp=mmp+1
                              tstrt(mp)=source_start(isource)+delta_ref*real(irelease-1)
                              weight(mp)=1.
                              weightd(mp)=1.
                              source_num(mp)=real(isource)
                              pmass(mp)=source_weight(isource)*totsrcmass*one_nparticles
                              frac1=RAN2( )
                              zran=source_h(inode)*frac1
                              frac2=RAN2( )
                              rcir=sqrt(frac2)*reff(isource)
                              ang1=2.*pi*RAN2( )
                              deltax=rcir*cos(ang1)
                              deltay=rcir*sin(ang1)
                              deltaz=zran
                              if(rcir.gt.1.e-06*dx)then
                                 ziloc=zran*source_r(inode)/(rcir*(1.+zran*source_r(inode)/(source_h(inode)*rcir)))
                              else
                                 ziloc=h
                              endif
                              vloc=vout*sqrt(zran**2+rcir**2)/sqrt(ziloc**2+source_r(inode)**2)
                              vlocz(mp)=vloc*zran/(sqrt(zran**2+rcir**2))
                              vlocx(mp)=vloc*rcir*cos(ang1)/sqrt(rcir**2+zran**2)
                              vlocy(mp)=vloc*rcir*sin(ang1)/sqrt(rcir**2+zran**2)
                              xps(mp)=source_x(inode)+deltax
                              yps(mp)=source_y(inode)+deltay
                              zps(mp)=deltaz+z0coeff*z0
                              im=int(xps(mp)/dx)+1
                              jm=int(yps(mp)/dy)+1
                              km=int((zps(mp)+dz)/dz)+1
                              if(icellflag(im,jm,km) .eq. 0 .or. icellflag_dense(im,jm) .ne. isource)then
                                 mp=mp-1
                                 mmp=mmp-1
                              endif
                           enddo
                           nptot=nptot+source_prr(isource)
                        enddo
                     endif
                  case(4) !Explosive Source
                     inode=source_startidx(isource)
                     if(ieradflag.eq. 0)then
!                        if(source_z(inode) .lt. z0)then
!                           source_z(inode)=1.1*z0
!                        endif
                        do mmp=1,source_prr(isource)
                           mp=mp+1
                           tstrt(mp)=source_start(isource)
                           ang2=2.*pi*RAN2()
                           !cosang1=2.*(RAN2()-.5) sphere
                           cosang1=RAN2() !hemisphere
                           ang1=acos(cosang1)
                           radius=source_r(inode)*(RAN2()**one_third)
                           deltax=radius*sin(ang1)*cos(ang2)
                           deltay=radius*sin(ang1)*sin(ang2)
                           deltaz=radius*cos(ang1)
                           xps(mp)=source_x(inode)+deltax
                           yps(mp)=source_y(inode)+deltay
                           zps(mp)=source_z(inode)+deltaz
                           weight(mp)=1.
                           weightd(mp)=1.
                           source_num(mp)=real(isource)
                           pmass(mp)=source_weight(isource)*totsrcmass*one_nparticles
                        enddo
                     else
                        iparticle_resample=0
                        do while (mp.lt.source_prr(isource))
                           mp=mp+1
                           tstrt(mp)=source_start(isource)
                           frac1=ran2( )
                           ipuff=npuff
                           do ipuf=1,npuff
                              if(frac1.le.pufcummas(ipuf)) then
                                 ipuff=ipuf
                                 exit
                              endif
                           enddo
                           dpm(mp)=dpm_bin(isub(ipuff))
                           call rangen(x)                           
                           deltaz=cldhgt*znorm(ipuff)+cldhgt*psigz(ipuff)*x
                           call rangen(x)
                           deltax=cldhgt*psigx(ipuff)*x
                           call rangen(x)
                           deltay=cldhgt*psigy(ipuff)*x
                           xps(mp)=source_x(inode)+deltax
                           yps(mp)=source_y(inode)+deltay
                           zps(mp)=source_z(inode)+deltaz
                           im=int(xps(mp)/dx)+1
                           jm=int(yps(mp)/dy)+1
                           km=int((zps(mp)+dz)/dz)+1
                           if(zps(mp).lt.z0+hgt(im,jm))then
                               zps(mp)=2.*(z0coeff*z0+hgt(im,jm)-zps(mp))+zps(mp)
                           endif
                           weight(mp)=1.
                           weightd(mp)=1.
                           source_num(mp)=real(isource)
                           pmass(mp)=source_weight(isource)*totsrcmass*one_nparticles
                           if(icellflag(im,jm,km) .eq. 0)then
                              mp=mp-1
                              iparticle_resample=iparticle_resample+1
                           else
                              pmass(mp)=source_weight(isource)*totsrcmass*one_nparticles
                           endif
                        enddo
                        if(iparticle_resample.gt.0)write(*,'(a,i6)')&
                              'number of resampled particles=',iparticle_resample
                     endif
                     nptot=nptot+source_prr(isource)
                  case(5) !Rectangular area or volume source
                     inode=source_startidx(isource)
                     if(source_z(inode).lt. 2.*z0)then
                        source_z(inode)=2.*z0
                     endif
                     do irelease=1,nreleases
                        do mmp=1,source_prr(isource)
                           mp=mp+1
                           tstrt(mp)=source_start(isource)+delta*real(irelease-1)
                           xsl=source_x(inode)
                           ysl=source_y(inode)-.5*source_w(inode)
                           zsl=source_z(inode)
                           xps(mp)=xsl+RAN2()*source_l(inode)
                           yps(mp)=ysl+RAN2()*source_w(inode)
                           zps(mp)=zsl+RAN2()*source_h(inode)
                           weight(mp)=1.
                           weightd(mp)=1.
                           source_num(mp)=real(isource)
                           pmass(mp)=source_weight(isource)*totsrcmass*one_nparticles
                        enddo
                        nptot=nptot+source_prr(isource)
                     enddo
                  case(6) !Moving Point Source
                     nreleases=int((source_nodestop(source_stopidx(isource)-1)-&
                                     source_nodestart(source_startidx(isource)))/delta)
                     inode=source_startidx(isource)
                     xsl=source_x(inode)
                     xsu=source_x(inode+1)
                     ysl=source_y(inode)
                     ysu=source_y(inode+1)
                     zsl=source_z(inode)
                     zsu=source_z(inode+1)
                     xpc=.5*(xsu+xsl)
                     ypc=.5*(ysu+ysl)
                     zpc=.5*(zsu+zsl)
                     seg_pos=0.
                     do irelease=1,nreleases
                        rel_time=source_start(isource)+delta*real(irelease-1)
                        if(rel_time .gt. source_nodestop(inode))then
                           inode=inode+1
                           if(inode.eq.source_stopidx(isource))exit
                           xsl=source_x(inode)
                           xsu=source_x(inode+1)
                           ysl=source_y(inode)
                           ysu=source_y(inode+1)
                           zsl=source_z(inode)
                           zsu=source_z(inode+1)
                           xpc=.5*(xsu+xsl)
                           ypc=.5*(ysu+ysl)
                           zpc=.5*(zsu+zsl)
                           seg_pos=0.
                        endif
                        if(inode.eq.source_stopidx(isource))exit
                        seg_time=rel_time-source_nodestart(inode)
                        seg_pos_old=seg_pos
                        seg_pos=(source_nodespeed(inode)*seg_time+&
                                0.5*source_nodeaccel(inode)*(seg_time**2.))/&
                                source_seglength(inode)
                        do mmp=1,source_prr(isource)
                           mp=mp+1
                           fracl=RAN2()*(seg_pos-seg_pos_old)+seg_pos_old
                           deltax=xsl+fracl*(xsu-xsl)-xpc
                           deltay=ysl+fracl*(ysu-ysl)-ypc
                           deltaz=zsl+fracl*(zsu-zsl)-zpc
                           tstrt(mp)=rel_time
                           xps(mp)=xpc+deltax
                           yps(mp)=ypc+deltay
                           zps(mp)=zpc+deltaz
                           weight(mp)=1.
                           weightd(mp)=1.
                           source_num(mp)=real(isource)
                           pmass(mp)=source_weight(isource)*totsrcmass*one_nparticles
                        enddo
                        nptot=nptot+source_prr(isource)
                     enddo
               endselect
            enddo
!MDW 11/29/2005 choose particle size bin
            if(lognormal_flag .eq. 1)then
               do mmp=1,nptot
                  if(itwophaseflag .eq. 1)then
                     frac1=RAN2()
                     if(frac1 .gt. agent_droplet_fraction)then
                        dpm(mmp)=0.
                        cycle
                     endif
                  endif
                  if(num_distributions .gt. 1)then
                     frac1=RAN2()
                     do ipuf=1,num_distributions
                        if(frac1.lt.distribution_cum_fraction(ipuf))then
                            exit
                        endif
                     enddo
                  else
                     ipuf=1
                  endif
                  frac2=RAN2()
                  if(distribution_std(ipuf)>1)then
                     dpm(mmp)=exp(log(distribution_mean(ipuf))+log(distribution_std(ipuf))*sqrt(2.)*&
                              erfinv(2.*frac2-1.))
                  else
                     dpm(mmp)=distribution_mean(ipuf)
                  endif
                  !Limit particle sizes to being between 0.001 and 400 microns
                  if(dpm(mmp) .gt. 400.)dpm(mmp)=400.
                  if(dpm(mmp) .lt. 0.001)dpm(mmp)=0.001
                  if(distribution_type_flag .eq. 2)then
                     pmass(mmp)=totsrcmass*distribution_mass_fraction(ipuf)*dpm(mmp)**3.
                     mass_sum=mass_sum+pmass(mmp)
                  endif
               enddo
               if(distribution_type_flag .eq. 2)then
                  pmass(:)=totsrcmass*pmass(:)/mass_sum
               endif
            else
               if(ieradflag.eq. 0)then
                  if(npuff .gt. 1 )then
                     do mmp=1,nptot
                        if(itwophaseflag .eq. 1)then
                           frac1=RAN2()
                           if(frac1 .gt. agent_droplet_fraction)then
                              dpm(mmp)=0.
                              cycle
                           endif
                        endif
                        frac1=RAN2()
                        do ipuf=1,npuff
                           if(frac1.lt.pufcummas(ipuf))then
                              dpm(mmp)=dpm_bin(ipuf)
                              exit
                           endif
                        enddo
                     enddo
                  else
                     if(itwophaseflag .eq. 1)then
                        do mmp=1,nptot
                           frac1=RAN2()
                           if(frac1 .gt. agent_droplet_fraction)then
                              dpm(mmp)=0.
                           else
                              dpm(mmp)=dpm_bin(1)
                           endif
                        enddo
                     else
                        dpm(:)=dpm_bin(1)
                     endif
                  endif
               endif
            endif
            if(ibioslurryflag .eq. 1 .or. itwophaseflag .eq. 1)then
               allocate(dpm_min(mpart))
               dpm_min(:)=0.
               if(ibioslurryflag .eq. 1)then
                  do mmp=1,nptot
                     agent_mass=agent_solids_concentration*4.*one_third*pi*((5.e-5*dpm(mmp))**3)
                     dpm_min(mmp)=2.0e4*(0.75*agent_mass/(pi*agent_eff_solid_density))**one_third
                  enddo
               endif
            endif
         endif
         return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine bulkvelc
         use variables
         implicit none
!         real delx,delri,uxsum,vysum,cellsum,uistarsum,ustars
         real sinphitdg,cosphitdg,deltdx,cyp1,cyp2,cym1,cym2,cxm1,cxm2,cxp1,cxp2
!         ,ristar0,utent,delvol,rold,hold,volold
!         real ubarold,drdxold,rlp,rlm,riminus,riplus,rplus,rminus
         real rlp,rlm,riminus,riplus,rplus,rminus !drdxold,
         real xip,yip
!         integer ndx,npl,nml,npm,icells,jcells,kcellst,ifinflag,irip
         integer ndx,npl,nml,npm,ifinflag,irip
         integer icellr,jcellr,npi,icellp,jcellp,irim,nmi
!         delx=.15*min(dx,dy)
         delri=.00051*min(dx,dy)
         ndx=8*max(nx,ny)
         idxmax=1
               uxsum=0.
               vysum=0.
               uistarsum=0.
               slopexsum=0.
               slopeysum=0.
               cellsum=0.
               ycompdevmax=0.
               ycompdevmin=100.
               rnewtotalm=0.
               rnewtotalp=0.
               idx=1
               icells=int(source_x(1)/dx)+1
               jcells=int(source_y(1)/dy)+1
               icell=icells
               jcell=jcells
               idxij(icell,jcell)=1
               kcellst=int((source_h(1)+.5*dz+dz)/dz)+1
               hnew(isource)=source_h(1)
               uxsum=u(icells,jcells,kcellst)+uxsum
               vysum=v(icells,jcells,kcellst)+vysum
               uistarsum=ustarij(icells,jcells,kcellst)+uistarsum
               call densegas_spprofile
               slopexsum=slopexsum+(elevij(icells+2,jcells+1)-elevij(icells,jcells+1))/(2.*dx)
               slopeysum=slopeysum+(elevij(icells+1,jcells+2)-elevij(icells+1,jcells))/(2.*dy)
               cellsum=cellsum+1.
               xscen(isource)=source_x(1)
               yscen(isource)=source_y(1)
               ustars=uistarsum/cellsum
               cosphidgw=uxsum/sqrt(uxsum**2+vysum**2)
               sinphidgw=vysum/sqrt(uxsum**2+vysum**2)
               cosphidgi(idx)=cosphidgw
               sinphidgi(idx)=sinphidgw
               vslope=sqrt((9.8*(source_density(isource)-rhoair)/rhoair)*sqrt(slopex**2+slopey**2)/&
                 ((kkar/alog(hnew(isource)/z0))**2+(kkar/(1.+.2*ristar0))**2+hnew(isource)*cd/(pi*rnew(isource))))
 !Note that we have used the small angle approximation and replaced the sine by the tangent
               if(vslope.gt.1.e-05)then
                  vslopex=-vslope*slopex/sqrt(slopex**2+slopey**2)
                  vslopey=-vslope*slopey/sqrt(slopex**2+slopey**2)
                  cosphidg(isource)=vslopex(isource)/sqrt(vslopex(isource)**2+vslopey(isource)**2)
                  sinphidg(isource)=vslopey(isource)/sqrt(vslopex(isource)**2+vslopey(isource)**2)
                  cosphidgi(idx)=cosphidg(isource)
                  sinphidgi(idx)=sinphidg(isource)
                  sinphitdg=cosphidg(isource)
                  cosphitdg=-sinphidg(isource)
                  ubars=vslope
                  sinslope=-sqrt(slopex**2+slopey**2)
               else
                  vslopex=0.
                  vslopey=0.
                  sinslope=0.
                  ubars=0.
               endif
               xsceni(idx)=xscen(isource)
               ysceni(idx)=yscen(isource)
               sinphitdg=cosphidg(isource)
               cosphitdg=-sinphidg(isource)
!               ubarw=(ustars/kkar)*((1.+z0/hnew(isource))*alog(1.+hnew(isource)/z0)-1.)
            if(rcl.gt.0)then
               phim=1.+4.7*rcl*.5*hnew(isource)
               psim=-4.7*rcl*.5*hnew(isource)
            else
               phim=(1.-15.*rcl*.5*hnew(isource))**(-.25)
               psim=2.*alog((1.+1./phim)/2.)+alog((1.+1./phim**2)/2.)-2.*atan(1./phim)+pi/2.
            endif
!            ubarw=(ustars/kkar)*((1.+z0/hnew(isource))*alog(1.+hnew(isource)/z0)-1.)
            ubarwtop=(ustars/kkar)*(alog((hnew(isource)/(2.*z0))-psim))
            ustarneutral=(kkar*ubarwtop)/alog(hnew(isource)/(2.*z0))
            ubarw=(ustarneutral/kkar)*((1.+z0/hnew(isource))*alog(1.+hnew(isource)/z0)-1.)
!               ubarwtop=(ustars/kkar)*
               ubar(isource)=max(ubarw,ubars)
            
               ubari(idx)=ubar(isource)
               ubulkpx(idx)=ubar(isource)*cosphidg(isource)
               vbulkpy(idx)=ubar(isource)*sinphidg(isource)
               rhocell(icells,jcells,2)=source_density(isource)
               idxij(icells,jcells)=1
               uxsum=0.
               vysum=0.
               cellsum=0.
               ycompdevmax=0.
               ycompdevmin=100.
               rnewtotalm=0.
               rnewtotalp=0.
               hnew(isource)=source_h(isource)
               rnew(isource)=source_r(isource)
               rnewpi(1)=rnew(isource)
               rnewmi(1)=rnew(isource)
               drdx=sqrt((9.8*(source_density(isource)-rhoair)/rhoair)*hnew(isource)&
                 )/ubar(isource)
                 delx=min(.15*dx,.15*dy,.2*rnew(isource)/drdx)
                vol0(isource)=(2.*source_r(isource))*source_h(isource)*delx
               volnew(isource)=vol0(isource)
                 ndx=int(1.4*max(real(nx)*dx,real(ny)*dy)/delx)+1
                 drdxpi(1)=drdx
                 drdxmi(1)=drdx
         do idx=1,ndx
            npl=0
            nml=0
            npm=0
            uxsum=0.
            vysum=0.
            cellsum=0.
            uistarsum=0.
               ycompdevmax=0.
               ycompdevmin=100.
            if(idx.eq.1)then
               hnewi(idx)=hnew(isource)
               rnewi(idx)=rnew(isource)
               icell1=int((xscen(isource)-rnew(isource))/dx)+1
               icell2=int((xscen(isource)+rnew(isource))/dx)+1
               jcell1=int((yscen(isource)-rnew(isource))/dy)+1
               jcell2=int((yscen(isource)+rnew(isource))/dy)+1
               do icell=icell1,icell2
                  do jcell=jcell1,jcell2
                     if(icell.ne.icells.or.jcell.ne.jcells)then
                        rhocell(icell,jcell,2)=rhocell(icells,jcells,2)
!                        ubulkpx(icell,jcell)=ubulkpx(icells,jcells)
!                        vbulkpy(icell,jcell)=vbulkpy(icells,jcells)
                        hnewij(icell,jcell)=hnew(isource)
                        idxij(icell,jcell)=1
                     endif
                  enddo
               enddo
            else
               icells=int((xscen(isource))/dx)+1
               jcells=int((yscen(isource))/dy)+1
               kcellst=int((hnew(isource)+.5*dz+dz)/dz)+1
            endif
            if(rcl.gt.0)then
               phim=1.+4.7*rcl*.5*hnew(isource)
               psim=-4.7*rcl*.5*hnew(isource)
            else
               phim=(1.-15.*rcl*.5*hnew(isource))**(-.25)
               psim=2.*alog((1.+1./phim)/2.)+alog((1.+1./phim**2)/2.)-2.*atan(1./phim)+pi/2.
            endif
!            ubarw=(ustars/kkar)*((1.+z0/hnew(isource))*alog(1.+hnew(isource)/z0)-1.)
            ubarwtop=(ustars/kkar)*(alog((hnew(isource)/(2.*z0))-psim))
            ustarneutral=(kkar*ubarwtop)/alog(hnew(isource)/(2.*z0))
            ubarw=(ustarneutral/kkar)*((1.+z0/hnew(isource))*alog(1.+hnew(isource)/z0)-1.)
            ubar(isource)=max(ubarw,ubars)
            ubari(idx)=ubar(isource)
            deltdx=delx/ubar(isource)
            if(cosphidg(isource).gt.0.)then
               cyp1=1./cosphidg(isource)
               cyp2=1.e+16*dy
               cym1=1.e+16*dy
               cym2=1./cosphidg(isource)
            else
               if(cosphidg(isource).lt.0)then
                  cyp1=1.e+16*dy
                  cyp2=-1./cosphidg(isource)
                  cym1=-1./cosphidg(isource)
                  cym2=1.e+16*dy
               else
                  cyp1=1.e+16*dy
                  cyp2=1.e+16*dy
                  cym1=1.e+16*dy
                  cym2=1.e+16*dy
               endif
            endif
            if(sinphidg(isource).gt.0.)then
               cxp1=1.e+16*dx
               cxp2=1./sinphidg(isource)
               cxm1=1./sinphidg(isource)
               cxm2=1.e+16*dx
            else
               if(sinphidg(isource).lt.0.)then
                  cxp1=-1./sinphidg(isource)
                  cxp2=1.e+16*dx
                  cxm1=1.e+16*dx
                  cxm2=-1./sinphidg(isource)
               else
                  cxp1=1.e+16*dx
                  cxp2=1.e+16*dx
                  cxm1=1.e+16*dx
                  cxm2=1.e+16*dx
               endif
            endif

! calculate forward differences
            drdx=sqrt((9.8*(source_density(isource)-rhoair)/rhoair)*vol0(isource)*hnew(isource)&
                 /volnew(isource))/ubar(isource)
            ristar0=9.8*(source_density(isource)-rhoair)*hnew(isource)*vol0(isource)&
                    /(volnew(isource)*rhoair*ustars**2)
            utent=0.5*ustars/(1.+.2*ristar0)
            if(sqrt(slopex**2+slopey**2).gt.1.e-05)then
               cosbeta=(-cosphidg(isource)*slopex-sinphidg(isource)*slopey)/sqrt(slopex**2+slopey**2)
               sinbeta=(sinphidg(isource)*slopex-cosphidg(isource)*slopey)/sqrt(slopex**2+slopey**2)
               sinslope=-sqrt(slopex**2+slopey**2)
            else
               cosbeta=1.
               sinbeta=0.
               sinslope=0.
            endif
            cosalph=(cosphidg(isource)*cosphidgw+sinphidg(isource)*sinphidgw)
            sinalph=(cosphidg(isource)*sinphidgw-sinphidg(isource)*cosphidgw)
            dubardx=-(9.8*(source_density(isource)-rhoair)/rhoair)*(vol0(isource)/&
               volnew(isource))*sinslope*cosbeta/ubar(isource)-&
               .1*ubar(isource)*drdx*(1-ubarw*cosalph/ubar(isource))&
               /rnew(isource)+(ubarw*cosalph/ubar(isource)-1)*&
               (.5*sqrt((ubar(isource)-ubarw*cosalph)**2+(ubarw*sinalph)**2)+.5*ustars)&
               /(hnew(isource)*(1+.2*ristar0))
            duperpdx=-(9.8*(source_density(isource)-rhoair)/rhoair)*(vol0(isource)/volnew(isource))*sinslope*sinbeta/ubar(isource)+&
               .1*drdx*(ubarw*sinalph)/rnew(isource)&
               +(ubarw*sinalph/ubar(isource))*(.5*sqrt((ubar(isource)-ubarw*cosalph)**2+(ubarw*sinalph)**2)+.5*ustars)&
               /(hnew(isource)*(1+.2*ristar0))               
            delvol=(2.*rnew(isource))*delx*utent*deltdx+2.*.1*drdx*delx*ubar(isource)*deltdx*hnew(isource)
            dhdx=.5*(sqrt((ubar(isource)-ubarw*cosalph)**2+(ubarw*sinalph)**2)+ustars)/(ubar(isource)*(1.+.2*ristar0))&
               +hnew(isource)*drdx*(-.9)/rnew(isource)-hnew(isource)*dubardx/ubar(isource)           
            delvol=(2.*rnew(isource))*delx*utent*deltdx+2.*.1*drdx*delx*ubar(isource)*deltdx*hnew(isource)
            rold=rnew(isource)
            hold=hnew(isource)
            ubar(isource)=ubar(isource)+dubardx*delx
            duperp=duperpdx*delx
            cosphidgn=(cosphidg(isource)*ubar(isource)-sinphidg(isource)*duperp)/sqrt(ubar(isource)**2+duperp**2)
            sinphidgn=(sinphidg(isource)*ubar(isource)+cosphidg(isource)*duperp)/sqrt(ubar(isource)**2+duperp**2)
            ubar(isource)=sqrt(ubar(isource)**2+duperp**2)
            rnew(isource)=rold+drdx*delx
             hnew=(volnew(isource)+delvol)/(2.*rnew*delx)
    !         hnew(isource)=dhdx*delx+hnew(isource)         
            dhdt=(hnew(isource)-hold)/deltdx
         dhdx=(hnew(isource)-hold)/delx
            volold=volnew(isource)
!   !         volnew(isource)=volnew(isource)+delvol
            volnew(isource)=2.*delx*hnew(isource)*rnew(isource)
! recalcute using backward difference
!            ubarold=ubar(isource)
!            drdxold=drdx
            ubarw=(ustars/kkar)*((1.+z0/hnew(isource))*alog(1.+hnew(isource)/z0)-1.)
            if(rcl.gt.0)then
               phim=1.+4.7*rcl*.5*hnew(isource)
               psim=-4.7*rcl*.5*hnew(isource)
            else
               phim=(1.-15.*rcl*.5*hnew(isource))**(-.25)
               psim=2.*alog((1.+1./phim)/2.)+alog((1.+1./phim**2)/2.)-2.*atan(1./phim)+pi/2.
            endif
!            ubarw=(ustars/kkar)*((1.+z0/hnew(isource))*alog(1.+hnew(isource)/z0)-1.)
            ubarwtop=(ustars/kkar)*(alog((hnew(isource)/(2.*z0))-psim))
            ustarneutral=(kkar*ubarwtop)/alog(hnew(isource)/(2.*z0))
            ubarw=(ustarneutral/kkar)*((1.+z0/hnew(isource))*alog(1.+hnew(isource)/z0)-1.)
!            drdx=sqrt((9.8*(source_density(isource)-rhoair)/rhoair)*vol0(isource)*hnew(isource)&
!                 /volnew(isource))/ubar(isource)
!            ristar0=9.8*(source_density(isource)-rhoair)*hnew(isource)*vol0(isource)&
!                    /(volnew(isource)*rhoair*ustars**2)
!            utent=0.5*ustars/(1.+.2*ristar0)
!            delvol=(2.*rnew(isource))*delx*utent*deltdx+2.*.7*drdx*delx*ubar(isource)*deltdx*hnew(isource)
!            hnew(isource)=(volnew(isource)+delvol)/(2.*rnew(isource)*delx)
!            ubar(isource)=.5*(ubarold+ubar(isource))
!            hnew(isource)=.5*(hnew(isource)+hold)
!            drdx=sqrt((9.8*(source_density(isource)-rhoair)/rhoair)*vol0(isource)*hnew(isource)&
!                 /volnew(isource))/ubar(isource)
!            ristar0=9.8*(source_density(isource)-rhoair)*hnew(isource)*vol0(isource)&
!                    /(volnew(isource)*rhoair*ustars**2)
!            utent=0.5*ustars/(1.+.2*ristar0)
!            delvol=(2.*rnew(isource))*delx*utent*deltdx+2.*.7*drdx*delx*ubar(isource)*deltdx*hnew(isource)
!            rnew(isource)=rold+drdx*delx
!            hnew(isource)=(volnew(isource)+delvol)/(2.*rnew(isource)*delx)
            dhdt=(hnew(isource)-hold)/delta
!            volnew(isource)=volold+delvol
            hnewi(idx+1)=hnew(isource)
            rnewi(idx+1)=rnew(isource)
         rnewpi(idx+1)=rnewpi(idx)+delx*drdx
         rnewmi(idx+1)=rnewmi(idx)+delx*drdx
            drdxi(idx)=drdx
            ubari(idx+1)=ubar(isource)
            xscen(isource)=xscen(isource)+delx*cosphidg(isource)
            yscen(isource)=yscen(isource)+delx*sinphidg(isource)
            xsceni(idx+1)=xscen(isource)
            ysceni(idx+1)=yscen(isource)
            icells=int((xscen(isource))/dx)+1
            jcells=int((yscen(isource))/dy)+1
            if(icells.gt.nx-1)exit
            if(icells.lt.0)exit
            if(jcells.lt.0)exit
            if(jcells.gt.ny-1)exit
            icells=max(icells,1)
            icells=min(icells,nx-1)
            jcells=max(jcells,1)
            jcells=min(jcells,ny-1)
            ycompdevmax=0.
            rplus=0.
            rminus=0.
            ycompdevmin=0.
            iloop_drdx_spread=0
            do while(iloop_drdx_spread.le.2)
            npm=0
            npl=0
            rplus=0
            rminus=0
            rlp=min(cyp1*(dy*real(jcells)-yscen(isource)),cyp2*(yscen(isource)-dy*real(jcells-1)))
            rlp=min(rlp,cxp1*(dx*real(icells)-xscen(isource)),cxp2*(xscen(isource)-dx*real(icells-1)))
            rlm=min(cym1*(dy*real(jcells)-yscen(isource)),cym2*(yscen(isource)-dy*real(jcells-1)))
            rlm=min(rlm,cxm1*(dx*real(icells)-xscen(isource)),cxm2*(xscen(isource)-dx*real(icells-1)))
            riminus=rlm
            riplus=rlp
            if(icellflag(icells,jcells,2) .ne. 0)then
               rhocell(icells,jcells,2)=(vol0(isource)/volnew(isource))*source_density(isource)&
                                        +(1.-vol0(isource)/volnew(isource))*rhoair
               hnewij(icells,jcells)=max(hnew(isource),hnewij(icells,jcells))
               dhdtij(icells,jcells)=dhdt(isource)
               idxij(icells,jcells)=idx+1
               npl=npl+1
               rpl(npl)=0.
               ixpl(npl)=icells
               jypl(npl)=jcells
               ustarr(icells,jcells)=kkar*drdx/alog((hnew(isource)+z0)/z0)
            else
               rplus=riplus
               rminus=riminus
            endif
            ifinflag=0
            do irip=1,1000
               riplus=riplus+delri
               xip=xscen(isource)-riplus*sinphidg(isource)
               yip=yscen(isource)+riplus*cosphidg(isource)
               if(yip.le.0.or.xip.le.0.or.xip.gt.(nx-1)*dx.or.yip.gt.(ny-1)*dy)then
                  ifinflag=1
               endif
               icellr=int((xip)/dx)+1
               jcellr=int((yip)/dy)+1
               icellr=max(1,icellr)
               icellr=min(nx-1,icellr)
               jcellr=max(1,jcellr)
               jcellr=min(ny-1,jcellr)
               rlp=min(cyp1*(dy*real(jcellr)-yip),cyp2*(yip-dy*real(jcellr-1)))
               rlp=min(rlp,cxp1*(dx*real(icellr)-xip),cxp2*(xip-dx*real(icellr-1)))
               if(icellflag(icellr,jcellr,2) .ne. 0 .and. riplus .le. rnew(isource)+rplus)then
                  rhocell(icellr,jcellr,2)=(vol0(isource)/volnew(isource))*source_density(isource)&
                                           +(1.-vol0(isource)/volnew(isource))*rhoair
                  hnewij(icellr,jcellr)=max(hnew(isource),hnewij(icellr,jcellr))
                  dhdtij(icellr,jcellr)=dhdt(isource)
                  if(idxij(icellr,jcellr).ne.idx+1)then
                     npl=npl+1
                     rpl(npl)=riplus-.5*rlp !characteristic radius for cell
                     ixpl(npl)=icellr
                     jypl(npl)=jcellr
                     idxij(icellr,jcellr)=idx+1
                     ymaxslice=max(ymaxslice,riplus)
                     ycompdevmax=max(ycompdevmax,riplus+(delx/ubarw)*&
                        (-u(icellr,jcellr,2)*sinphidg(isource)+v(icellr,jcellr,2)*cosphidg(isource)))
                  endif
               else
                  if(idxij(icellr,jcellr).ne.idx+1 .and. icellflag(icellr,jcellr,2).eq.0)then
                     rplus=rplus+rlp
                     rplusij(idx)=max(rplus,rplusij(idx))
                  endif
               endif
               riplus=riplus+rlp
!               rplusij(idx)=rplus
               ubulkpx(idx)=ubar(isource)*cosphidg(isource)
               vbulkpy(idx)=ubar(isource)*(sinphidg(isource))
               if(riplus.gt.rplusij(idx)+rnewpi(idx+1).or.ifinflag.eq.1)then
                  do npi=1,npl
                     icellp=ixpl(npi)
                     jcellp=jypl(npi)
                     hnewij(icellp,jcellp)=max(hnew(isource),hnewij(icellp,jcellp))
                     dhdtij(icellp,jcellp)=dhdt(isource)
              !       ubulkpx(icellp,jcellp)=-rpl(npi)*drdx*sinphidg(isource)/(rnew(isource)+rplusij(icellp,jcellp))
              !       vbulkpy(icellp,jcellp)=rpl(npi)*drdx*cosphidg(isource)/(rnew(isource)+rplusij(icellp,jcellp))
                     uxsum=u(icellp,jcellp,kcellst)+uxsum
                     vysum=v(icellp,jcellp,kcellst)+vysum
                     uistarsum=ustarij(icellp,jcellp,kcellst)+uistarsum
                     slopexsum=slopexsum+(elevij(icellp+2,jcellp+1)-elevij(icellp,jcellp+1))/(2.*dx)
                     slopeysum=slopeysum+(elevij(icellp+1,jcellp+2)-elevij(icellp+1,jcellp))/(2.*dy)
                     cellsum=cellsum+1.
                     ustarr(icellp,jcellp)=kkar*drdx/alog((hnew(isource)+z0)/z0)
                  enddo
                  exit
               endif
            enddo
            ifinflag=0
            do irim=1,1000
               riminus=riminus+delri
               xip=xscen(isource)+riminus*sinphidg(isource)
               yip=yscen(isource)-riminus*cosphidg(isource)
               if(yip.le.0.or.xip.le.0.or.xip.gt.(nx-1)*dx.or.yip.gt.(ny-1)*dy)then
                  ifinflag=1
               endif
               icellr=int((xip)/dx)+1
               jcellr=int((yip)/dy)+1
               icellr=max(1,icellr)
               icellr=min(nx-1,icellr)
               jcellr=max(1,jcellr)
               jcellr=min(ny-1,jcellr)
               rlm=min(cym1*(dy*real(jcellr)-yip),cym2*(yip-dy*real(jcellr-1)))
               rlm=min(rlm,cxm1*(dx*real(icellr)-xip),cxm2*(xip-dx*real(icellr-1)) )
               if(icellflag(icellr,jcellr,2) .ne. 0 .and. riminus .le. rminusij(idx)+rnewmi(idx+1))then
                  rhocell(icellr,jcellr,2)=(vol0(isource)/volnew(isource))*source_density(isource)&
                                           +(1.-vol0(isource)/volnew(isource))*rhoair
                  hnewij(icellr,jcellr)=max(hnew(isource),hnewij(icellr,jcellr))
                  dhdtij(icellr,jcellr)=dhdt(isource)
                  if(idxij(icellr,jcellr).ne.idx+1)then
                     idxij(icellr,jcellr)=idx+1
                     npm=npm+1
                     rpm(npm)=riminus-.5*rlm  !characteristic radius for cell
                     ixml(npm)=icellr
                     jyml(npm)=jcellr
                     yminslice=max(yminslice,riminus)
                     ycompdevmin=min(ycompdevmin,-riminus+(delx/ubarw)*&
                        (u(icellr,jcellr,2)*sinphidg(isource)-v(icellr,jcellr,2)*cosphidg(isource)))
                     
                  endif
               else
                  if(idxij(icellr,jcellr).ne.idx+1)then
                     idxij(icellr,jcellr)=idx+1
                     if(icellflag(icellr,jcellr,2) .eq. 0)then
                        rminus=rminus+rlm
                        rminusij(idx)=max(rminus,rminusij(idx))
                     endif
                  endif
               endif
!               rminusij(idx)=rminus
               riminus=rlm+riminus
               if(riminus.ge.rminusij(idx)+rnewmi(idx+1).or.ifinflag.eq.1)then
                  do nmi=1,npm
                     icellp=ixml(nmi)
                     jcellp=jyml(nmi)
                     hnewij(icellp,jcellp)=max(hnew(isource),hnewij(icellp,jcellp))
                     dhdtij(icellp,jcellp)=dhdt(isource)
!                     ubulkpx(icellp,jcellp)=rpm(nmi)*drdx*sinphidg(isource)/(rnew(isource)+rminusij(icellp,jcellp))
!                     vbulkpy(icellp,jcellp)=-rpm(nmi)*drdx*cosphidg(isource)/(rnew(isource)+rminusij(icellp,jcellp))
                     uxsum=u(icellp,jcellp,kcellst)+uxsum
                     vysum=v(icellp,jcellp,kcellst)+vysum
                     uistarsum=ustarij(icellp,jcellp,kcellst)+uistarsum
                     slopexsum=slopexsum+(elevij(icellp+2,jcellp+1)-elevij(icellp,jcellp+1))/(2.*dx)
                     slopeysum=slopeysum+(elevij(icellp+1,jcellp+2)-elevij(icellp+1,jcellp))/(2.*dy)
                     cellsum=cellsum+1.
                     ustarr(icellp,jcellp)=kkar*drdx/alog((hnew(isource)+z0)/z0)
                  enddo
                  exit
               endif
            enddo
            iloop_drdx_spread=iloop_drdx_spread+1
            drdxplus=(ycompdevmax-ymaxslice)/delx
            drdxminus=(ycompdevmin+yminslice)/delx
            drdxminus=0.
            drdxplus=0.
            drdxpi(idx+1)=max(drdx,drdxplus)
            drdxmi(idx+1)=max(drdx,-drdxminus)
            rnewpi(idx+1)=rnewpi(idx)+max(drdx,drdxplus)*delx
            rnewmi(idx+1)=rnewmi(idx)+max(drdx,-drdxminus)*delx
            rnewtotalm=max(rnewtotalm,rminusij(idx)+rnewmi(idx+1))
            if(rnewmi(idx+1)+rminus.lt.0.99999*rnewtotalm)rnewmi(idx+1)=rnewtotalm-rminusij(idx)
            rnewtotalp=max(rnewtotalp,rplusij(idx)+rnewpi(idx+1))
            if(rnewpi(idx+1)+rplus.lt.0.99999*rnewtotalp)rnewpi(idx+1)=rnewtotalp-rplusij(idx)
            enddo
            hnewi(idx+1)=hnewi(idx+1)*(2.*rnewi(idx+1)/(rnewpi(idx+1)+rnewmi(idx+1)))
            hnew(isource)=hnewi(idx+1)
            if(cellsum.gt.0.)then
               ustars=uistarsum/cellsum
               slopex=slopexsum/cellsum
               slopey=slopeysum/cellsum
               cosphidgw=uxsum/sqrt(uxsum**2+vysum**2)
               sinphidgw=vysum/sqrt(uxsum**2+vysum**2)
               cosphidgi(idx+1)=cosphidgn
               sinphidgi(idx+1)=sinphidgn
               sinphitdg=cosphidg(isource)
               cosphitdg=-sinphidg(isource)
               cosphidg(isource)=cosphidgn
               sinphidg(isource)=sinphidgn
             endif
            if(idx.gt.idxmax)idxmax=idx
         enddo
         write(67,110)
110      format('   xscen    yscen    hnewi     rnewi      ubulkx       vbulkx  rhocent')
         do idx=1,idxmax
            crat=(hnewi(1)*rnewi(1))/(hnewi(idx)*rnewi(idx))
            rhocent(isource)=source_density(isource)*crat+rhoair*(1.-crat)
            write(67,111)xsceni(idx),ysceni(idx),hnewi(idx),rnewi(idx),ubulkpx(idx),vbulkpy(idx),rhocent(isource)
         enddo
111      format(6f12.5,f12.6)
         return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine bulkveli
         use variables
         implicit none
!         real delri,uxsum,vysum,cellsum,uistarsum,ycellcn,slopeysum,slopexsum 
!         real xcellcl,xcellcn,ycellcl,rcellcl,rcellcn,addarea
!         real effarea,curarea,ustars,volnewp,ristar0,utent,delvol,delvol_slope,cd
!         real rold,hold,volold,ubarold
!         integer kcellst,iconarea,icells,jcells,istrt,ifin
!         integer jstrt,jfin,icell,jcell
         delri=.0005*min(dx,dy)
         cd=1.
!         do isource=1,nsources
            vol0(isource)=pi*source_h(isource)*source_r(isource)**2
            volnew(isource)=pi*hnew(isource)*rnew(isource)**2
            if(t .le. 0.1*delta .and. inextgridflag .eq. 2)reff(isource)=rnew(isource)
            if(t .le. 0.1*delta .and. inextgridflag .ne. 2)then
               rnew(isource)=source_r(isource)
               reff(isource)=source_r(isource)
               hnew(isource)=source_h(isource)
               volnew(isource)=pi*hnew(isource)*source_r(isource)**2
               vol0(isource)=volnew(isource)
            else if(t .le. 0.1*delta)then
               kcellst=int((hnew(isource)+.5*dz+dz)/dz)+1
               iconarea=0
               rnewp(isource)=rnew(isource)
               do while(iconarea.lt.20)
                  uxsum=0.
                  vysum=0.
                  cellsum=0.
                  uistarsum=0.
                  utotsum=0.
                  slopexsum=0.
                  slopeysum=0.
                  xscen(isource)=source_x(isource)
                  yscen(isource)=source_y(isource)
                  addarea=0.
                  icells=int((source_x(isource))/dx)+1
                  jcells=int((source_y(isource))/dy)+1
                  rhocent(isource)=(vol0(isource)/volnew(isource)*source_density(isource)&
                                   +(1.-vol0(isource)/volnew(isource)))*rhoair
                  istrt=int((source_x(isource)-reff(isource)-dx)/dx)+1
                  ifin=int((source_x(isource)+reff(isource)+dx)/dx)+1
                  jstrt=int((source_y(isource)-reff(isource)-dy)/dy)+1
                  jfin=int((source_y(isource)+reff(isource)+dy)/dy)+1
                  do icell=istrt,ifin
                     xcellcn=dx*real(icell-1)+.5*dx
                     do jcell=jstrt,jfin
                        ycellcn=dy*real(jcell-1)+.5*dy
                        if(icells.ge.icell)then
                           xcellcl=dx*real(icell)-delri
                        else
                           xcellcl=dx*real(icell-1)-delri
                        endif
                        if(jcells.ge.jcell)then
                           ycellcl=dy*real(jcell)-delri
                        else
                           ycellcl=dy*real(jcell-1)-delri
                        endif
                        rcellcl=sqrt((xcellcl-xscen(isource))**2+(ycellcl-yscen(isource))**2)
                        if((rcellcl.le.reff(isource)).or.(icells.eq.icell.and.jcells.eq.jcell))then
                           if(icellflag(icell,jcell,2) .ne. 0 .and.& 
                              (icellflag_dense(icell,jcell).eq.isource .or. icellflag_dense(icell,jcell).eq.0))then
!                             rhocell(icell,jcell,2)=(vol0(isource)/volnew(isource))*source_density(isource)+(1.-vol0(isource)/volnew(isource))*rhoair
                              uxsum=u(icell,jcell,kcellst)+uxsum
                              vysum=v(icell,jcell,kcellst)+vysum
                              uistarsum=ustarij(icell,jcell,kcellst)+uistarsum
                              cellsum=cellsum+1.
                              call densegas_spprofile
                              slopexsum=slopexsum+(elevij(icell+2,jcell+1)-elevij(icell,jcell+1))/(2.*dx)
                              slopeysum=slopeysum+(elevij(icell+1,jcell+2)-elevij(icell+1,jcell))/(2.*dy)
                              icellflag_dense(icell,jcell)=isource
                           else
                              rcellcn=sqrt((xcellcn-xscen(isource))**2+(ycellcn-yscen(isource))**2)
                              if(rcellcn.le.reff(isource))then
                                 addarea=addarea+dx*dy
                              endif
                           endif
                        endif
                     enddo
                  enddo
                  effarea=pi*rnew(isource)**2+addarea
                  curarea=pi*reff(isource)**2
                  if(curarea.lt..99*effarea)then
                     reff(isource)=sqrt(effarea/pi)
                     iconarea=iconarea+1
                  else
                     exit
                  endif
               enddo
               xscen(isource)=source_x(isource)
               yscen(isource)=source_y(isource)
               ustars=uistarsum/cellsum
               utotdg=utotsum/cellsum
               cosphidg(isource)=uxsum/sqrt(uxsum**2+vysum**2)
               sinphidg(isource)=vysum/sqrt(uxsum**2+vysum**2)
               slopex=slopexsum/cellsum
               slopey=slopeysum/cellsum
            endif
            rnewp(isource)=rnew(isource)
            reffp(isource)=reff(isource)
            hnewp(isource)=hnew(isource)
            xscenp(isource)=xscen(isource)
            yscenp(isource)=yscen(isource)
            volnewp(isource)=volnew(isource)
!           ubar(isource)=(ustars/kkar)*((1.+z0/hnew(isource))*alog(1.+hnew(isource)/z0)-1.)
            ubar(isource)=utotdg
            drdtp(isource)=sqrt((9.8*(source_density(isource)-rhoair)/rhoair)*vol0(isource))/rnewp(isource)
! calculate forward differences
            drdt(isource)=sqrt((9.8*(source_density(isource)-rhoair)/rhoair)*vol0(isource))/rnew(isource)
            ristar0=9.8*(source_density(isource)-rhoair)*hnew(isource)*vol0(isource)/(volnew(isource)*rhoair*ustars**2)
            utent=0.5*ustars/(1.+.2*ristar0)
            delvol=(pi*rnew(isource)**2)*utent*delta+6.28*.01*drdt(isource)*rnew(isource)*delta*hnew(isource)
            vslope=sqrt((9.8*(source_density(isource)-rhoair)/rhoair)*sqrt(slopex**2+slopey**2)/&
              ((kkar/log(hnew(isource)/z0))**2+(kkar/(1.+.2*ristar0))**2+hnew(isource)*cd/(pi*rnew(isource))))
            if(vslope.gt.1.e-05)then
               vslopex(isource)=-vslope*slopex/sqrt(slopex**2+slopey**2)
               vslopey(isource)=-vslope*slopey/sqrt(slopex**2+slopey**2)
            else
               vslopex(isource)=0.
               vslopey(isource)=0.
            endif
            delvol_slope=delta*(rnew(isource)*hnew(isource)*cd*vslope+kkar*pi*rnew(isource)**2*vslope/(1.+.2*ristar0))
            rold=rnew(isource)
            hold=hnew(isource)
            rnew(isource)=rold+drdt(isource)*delta
            hnew(isource)=(volnew(isource)+delvol+delvol_slope)/(pi*rnew(isource)**2)
            dhdt(isource)=(hnew(isource)-hold)/delta
            volold=volnew(isource)
            volnew(isource)=volnew(isource)+delvol
! recalcute using backward difference
            ubarold=ubar(isource)
!           ubar(isource)=(ustars/kkar)*((1.+z0/hnew(isource))*alog(1.+hnew(isource)/z0)-1.)
            ubar(isource)=utotdg
            drdt(isource)=sqrt((9.8*(source_density(isource)-rhoair)/rhoair)*vol0(isource))/rnew(isource)
            ristar0=9.8*(source_density(isource)-rhoair)*hnew(isource)*vol0(isource)/(volnew(isource)*rhoair*ustars**2)
            utent=0.5*ustars/(1.+.2*ristar0)
            delvol=(pi*rnew(isource)**2)*utent*delta+6.28*.01*drdt(isource)*rnew(isource)*delta*hnew(isource)
            vslope=sqrt((9.8*(source_density(isource)-rhoair)/rhoair)*sqrt(slopex**2+slopey**2)/&
              ((kkar/log(hnew(isource)/z0))**2+(kkar/(1.+.2*ristar0))**2+hnew(isource)*cd/(pi*rnew(isource))))
            if(vslope.gt.1.e-05)then
               vslopex(isource)=-vslope*slopex/sqrt(slopex**2+slopey**2)
               vslopey(isource)=-vslope*slopey/sqrt(slopex**2+slopey**2)
            else
               vslopex(isource)=0.
               vslopey(isource)=0.
            endif
            delvol_slope=delta*(rnew(isource)*hnew(isource)*cd*vslope**2+kkar*pi*rnew(isource)**2*vslope/(1.+.2*ristar0))
            hnew(isource)=(volnew(isource)+delvol+delvol_slope)/(pi*rnew(isource)**2)
            ubar(isource)=.5*(ubarold+ubar(isource))
            hnew(isource)=.5*(hnew(isource)+hold)
            drdt(isource)=sqrt((9.8*(source_density(isource)-rhoair)/rhoair)*vol0(isource))/rnew(isource)
            ristar0=9.8*(source_density(isource)-rhoair)*hnew(isource)*vol0(isource)/(volnew(isource)*rhoair*ustars**2)
            utent=0.5*ustars/(1.+.2*ristar0)
            delvol=(pi*rnew(isource)**2)*utent*delta+6.28*.01*drdt(isource)*rnew(isource)*delta*hnew(isource)
            vslope=sqrt((9.8*(source_density(isource)-rhoair)/rhoair)*sqrt(slopex**2+slopey**2)/&
              ((kkar/log(hnew(isource)/z0))**2+(kkar/(1.+.2*ristar0))**2+hnew(isource)*cd/(pi*rnew(isource))))
            if(vslope.gt.1.e-05)then
               vslopex(isource)=-vslope*slopex/sqrt(slopex**2+slopey**2)
               vslopey(isource)=-vslope*slopey/sqrt(slopex**2+slopey**2)
            else
               vslopex(isource)=0.
               vslopey(isource)=0.
            endif
            delvol_slope=delta*(rnew(isource)*hnew(isource)*cd*vslope+kkar*pi*rnew(isource)**2*vslope/(1.+.2*ristar0))         
            rnew(isource)=rold+drdt(isource)*delta
            reff(isource)=rnew(isource)
            hnew(isource)=(volnew(isource)+delvol+delvol_slope)/(pi*rnew(isource)**2)
            dhdt(isource)=(hnew(isource)-hold)/delta
            volnew(isource)=volold+delvol+delvol_slope
            xscen(isource)=xscen(isource)+ubar(isource)*delta*cosphidg(isource)+vslopex(isource)*delta
            yscen(isource)=yscen(isource)+ubar(isource)*delta*sinphidg(isource)+vslopey(isource)*delta
            kcellst=int((hnew(isource)+.5*dz+dz)/dz)+1
            icells=int((xscen(isource))/dx)+1
            jcells=int((yscen(isource))/dy)+1
            rhocent(isource)=(vol0(isource)/volnew(isource))*source_density(isource)+(1.-vol0(isource)/volnew(isource))*rhoair
            if(icells.gt.nx-1)return
            if(icells.lt.0)return
            if(jcells.lt.0)return
            if(jcells.gt.ny-1)return
            icells=max(icells,1)
            icells=min(icells,nx-1)
            jcells=max(jcells,1)
            jcells=min(jcells,ny-1)
            iconarea=0
            do while(iconarea.lt.20)
               uxsum=0.
               vysum=0.
               cellsum=0.
               uistarsum=0.
               slopexsum=0.
               slopeysum=0.
               addarea=0.
               utotsum=0.
               istrt=int((xscen(isource)-reff(isource)-dx)/dx)+1
               ifin=int((xscen(isource)+reff(isource)+dx)/dx)+1
               jstrt=int((yscen(isource)-reff(isource)-dy)/dy)+1
               jfin=int((yscen(isource)+reff(isource)+dy)/dy)+1
               istrt=max(1,istrt)
               ifin=min(ifin,nx-1)
               jstrt=max(1,jstrt)
               jfin=min(jfin,ny-1)
               rhocent(isource)=(vol0(isource)/volnew(isource))*source_density(isource)+(1.-vol0(isource)/volnew(isource))*rhoair
!               do icell=istrt,ifin
               do icell=1,nx
                   do jcell=1,ny
                      if((icell.le.ifin.and.icell.ge.istrt).and.(jcell.le.jfin.and.jcell.ge.jstrt))then
                         xcellcn=dx*real(icell-1)+.5*dx
!                       do jcell=jstrt,jfin
                        ycellcn=dy*real(jcell-1)+.5*dy
                        if(icells.ge.icell)then
                           xcellcl=dx*real(icell)-delri
                        else
                           xcellcl=dx*real(icell-1)+delri
                        endif
                        if(jcells.ge.jcell)then
                           ycellcl=dy*real(jcell)-delri
                        else
                           ycellcl=dy*real(jcell-1)-delri
                        endif
                        rcellcl=sqrt((xcellcl-xscen(isource))**2+(ycellcl-yscen(isource))**2)
                        if((rcellcl.le.reff(isource)).or.(icells.eq.icell.and.jcells.eq.jcell))then
                           if(icellflag(icell,jcell,2) .ne. 0 .and. &
                              (icellflag_dense(icell,jcell) .eq. 0 .or. icellflag_dense(icell,jcell) .eq. isource))then
                              rhocell(icell,jcell,2)=(vol0(isource)/volnew(isource))*source_density(isource)&
                                                     +(1.-vol0(isource)/volnew(isource))*rhoair
                              uxsum=u(icell,jcell,kcellst)+uxsum
                              vysum=v(icell,jcell,kcellst)+vysum
                              uistarsum=ustarij(icell,jcell,kcellst)+uistarsum
                              slopexsum=slopexsum+(elevij(icell+2,jcell+1)-elevij(icell,jcell+1))/(2.*dx)
                              slopeysum=slopeysum+(elevij(icell+1,jcell+2)-elevij(icell+1,jcell))/(2.*dy)                        
                              cellsum=cellsum+1.
                              call densegas_spprofile
                              icellflag_dense(icell,jcell)=isource
                           else
                              rcellcn=sqrt((xcellcn-xscen(isource))**2+(ycellcn-yscen(isource))**2)
                              if(icellflag_dense(icell,jcell).eq.isource)icellflag_dense(icell,jcell)=0
                              if(rcellcn.le.reff(isource))then
                                 addarea=addarea+dx*dy
                              endif
                           endif
                        endif
                     else
                        if(icellflag_dense(icell,jcell).eq.isource)icellflag_dense(icell,jcell)=0
                     endif
                  enddo
               enddo
               effarea=pi*rnew(isource)**2+addarea
               curarea=pi*reff(isource)**2
               if(curarea.lt..99*effarea)then
                  iconarea=1+iconarea
                  reff(isource)=sqrt(effarea/pi)
               else
                  exit
               endif
            enddo
            ustars=uistarsum/cellsum
            cosphidg(isource)=uxsum/sqrt(uxsum**2+vysum**2)
            sinphidg(isource)=vysum/sqrt(uxsum**2+vysum**2)
            slopex=slopexsum/cellsum
            slopey=slopeysum/cellsum
            utotdg=utotsum/cellsum
            !rnewi(it)=rnew
!         enddo
         return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!This subroutine calculates the contribution from each horizonal cell to the average
!speed over the dense-gas slug
      subroutine densegas_spprofile
      use variables
      kst=int((.5*dz+hgt(icell,jcell)+dz)/dz)+1
      hnum=min(hnew(isource)-hgt(icell,jcell),.5*dz)
      hdenom=hnew(isource)-hgt(icell,jcell)
      utotl=sqrt(u(icell,jcell,kst)**2+v(icell,jcell,kst)**2)
      ustarsl=kkar*utotl/alog((hnum+z0)/z0)
      utotsum=utotsum+(ustarsl/(kkar*hdenom))*((hnum+z0)*alog((hnum+z0)/z0)-hnum)
      if(hdenom.lt.0.5*dz)return
      utotu=sqrt(u(icell,jcell,kst+1)**2+v(icell,jcell,kst+1)**2)
      zup=min(hnew(isource)-hgt(icell,jcell),dz)
      zlo=.5*dz
      z1=.5*dz
      utotsum=(1./hdenom)*(utotl*(zup-zlo)+((utotu-utotl)/(2.*dz)*(zup**2-zlo**2-2.*zup*z1+2.*zlo*z1)))+utotsum
      do k=kst+1,kcellst
         zlo=z1+.5*dz
         zup=min(hnew(isource)-hgt(icell,jcell),zlo+.5*dz)
         utotsum=(1./hdenom)*(utotl*(zup-zlo)+((utotu-utotl)/(2.*dz)*(zup**2-zlo**2-2.*zup*z1+2.*zlo*z1)))+utotsum
         if(zup.lt.zlo+.5*dz)return
         utotl=utotu
         utotu=sqrt(u(icell,jcell,k+1)**2+v(icell,jcell,k+1)**2)
         z1=z1+dz
         zlo=z1
         zup=min(hnew(isource)-hgt(icell,jcell),zlo+.5*dz)
         utotsum=(1./hdenom)*(utotl*(zup-zlo)+((utotu-utotl)/(2.*dz)*(zup**2-zlo**2-2.*zup*z1+2.*zlo*z1)))+utotsum
      enddo
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calculateBuildingInfiltExfilt
         use variables
         implicit none
         real deposition_velocity
         if(ibioslurryflag .eq. 1 )then
            call computeMMD(nptot,dpm_min(1:nptot),pmass(1:nptot),massMedianDiameter)
         elseif(itwophaseflag .eq. 1)then
            massMedianDiameter=0.
         else
            call computeMMD(nptot,dpm(1:nptot),pmass(1:nptot),massMedianDiameter)
         endif
         if(depvel .ge. 0)then
            deposition_velocity=depvel/100
         else
            call int3cal(massMedianDiameter,sc,taud,rhop,xnu,press,temp,diff,vs)
            ustar=0.2
            taup=taud*1.e+04*ustar**2/xnu
            if(massMedianDiameter.gt.0.)then
               if(taup.gt..10)then
                  rb=1./((sc**(-2./3.)+10.**(-3./taup))*ustar)
               else
                  rb=1./((sc**(-2./3))*ustar)
               endif
            else
               rb=5.*sc**(2./3.)/ustar
            endif
! mdw 4-15-2004 initialize dnear
            ra=log((.5*dz*ustar*1.e+04/xnu+diff/xnu)/(ustar*1.e+02/xnu+diff/xnu))/(kkar*ustar)
            deposition_velocity=1./(ra+rb)
         endif
         do ibuild=1,inumbuildfile
            if(bldtype(ibuild) .ge. 9)then
               cycle
            endif
            if(racintake(ibuild) .lt. 0)then
                if(toact .lt. 273.)then
                    racintake(ibuild)=0.5
                elseif(toact .ge. 273. .and. toact .lt. 290.)then
                    racintake(ibuild)=0.176*toact-47.55
                elseif(toact .ge. 290. .and. toact .lt. 295.)then
                    racintake(ibuild)=-0.46*toact+136.2
                else
                    racintake(ibuild)=0.5
                endif
            endif
            rhoref(ibuild)=.0288*piref(ibuild)*(273./toref(ibuild))/.0224
            psref=rhoref(ibuild)*9.812*Ht(ibuild)*abs(tiref(ibuild)-toref(ibuild))/tiref(ibuild)
            psact=rhoact*9.812*Ht(ibuild)*abs(tiact(ibuild)-toact)/tiact(ibuild) !these expressions
! give the buoyancy pressures for the reference and actual cases
            pwref=.5*rhoref(ibuild)*urefb(ibuild)**2
            pwact=.5*rhoact*(utotmax(ibuild)**2-utotcl1(ibuild)**2)
            aref=.3**(1.5)*psref+.19**(1.5)*pwref-.333*sqrt((.3*.19)**1.5*psref*pwref)
            aact=.3**(1.5)*psact+.19**(1.5)*pwact-.333*sqrt((.3*.19)**1.5*psact*pwact)
            racact(ibuild)=racref(ibuild)*(aact/aref)**(2./3.)
            ainf(ibuild)=(racact(ibuild)+racintake(ibuild)*(1-filtration(ibuild)))/3600.
            binf(ibuild)=(racact(ibuild)+racintake(ibuild)+&
                         (racsupply(ibuild)-racintake(ibuild))*filtration(ibuild)+&
                         3600.*deposition_velocity*area2volume(ibuild)+3600.*decay)/3600.
            abinf(ibuild)=ainf(ibuild)/binf(ibuild)
            binfex(ibuild)=(racact(ibuild)+racintake(ibuild))/3600.
            tau(ibuild)=3600./racact(ibuild)
            if(tau(ibuild).lt.taumin) taumin=tau(ibuild)
         enddo
         return
      end
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine computeMMD(n,diameter,mass,massMedianDiameter)
         implicit none
         integer,intent(in)::n
         real,intent(in)::diameter(n),mass(n)
         real,intent(out)::massMedianDiameter
         real sorted_diameter(n),sorted_mass(n),totalMass,tempD,tempM,massSum
         integer i,ndelta,inorder
         sorted_diameter=diameter
         sorted_mass=mass
         totalMass=sum(mass)
         ndelta=n
         inorder=1
         do while(ndelta .gt. 1)
            if(inorder .eq. 1)ndelta=ndelta/2
            inorder=1
            do i=1,n-ndelta
               if(sorted_diameter(i) .gt. sorted_diameter(i+ndelta))then
                  tempD=sorted_diameter(i)
                  tempM=sorted_mass(i)
                  sorted_diameter(i)=sorted_diameter(i+ndelta)
                  sorted_mass(i)=sorted_mass(i+ndelta)
                  sorted_diameter(i+ndelta)=tempD
                  sorted_mass(i+ndelta)=tempM
                  inorder=0
               endif
            enddo
         enddo
         massSum=0.
         do i=1,n
            massSum=massSum+sorted_mass(i)
            if(massSum .ge. 0.5*totalMass)then
               massMedianDiameter=sorted_diameter(i)
               exit
            endif
         enddo
!         print*,'Mass Median Diameter =',massMedianDiameter,'microns'
         return
      end

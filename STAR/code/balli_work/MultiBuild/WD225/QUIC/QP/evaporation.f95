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
      
      


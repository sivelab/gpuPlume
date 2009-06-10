   subroutine MoveParticles
   Use disperseModule
   Use datamodule
   Implicit None

   
   real, parameter:: Co=5.7
   real uMem,vMem,wMem,uDri,vDri,wDri,uRan,vRan,wRan,uPrime,vPrime,wPrime
   real uDri1,uDri2,uDri3,uDri4,uDri5,uDri6,uDri7,uDri8,uDri9,uDri10,uDri11,uDri12,uDri13,uDri14,uDri15,uDri16,uDri17,uDri18,uDri19
   real vDri1,vDri2,vDri3,vDri4,vDri5,vDri6,vDri7,vDri8,vDri9,vDri10,vDri11,vDri12,vDri13,vDri14,vDri15,vDri16,vDri17,vDri18,vDri19
   real wDri1,wDri2,wDri3,wDri4,wDri5,wDri6,wDri7,wDri8,wDri9,wDri10,wDri11,wDri12,wDri13,wDri14,wDri15,wDri16,wDri17,wDri18,wDri19                 
   real x,disIncrX,disIncrY,disIncrZ,Eps
   real du,dv,dw,sigUCell,sigVCell,sigWCell,sigUCellPrev,sigVCellPrev,sigWCellPrev,sigUCellNext,sigVCellNext,sigWCellNext
    

   
   
   uPrimePrev(pMov)=uPrimeCurr(pMov)
   vPrimePrev(pMov)=vPrimeCurr(pMov)
   wPrimePrev(pMov)=wPrimeCurr(pMov)
   
   xPosCell=int(xPos(pMov)/dx) + 1                        !Getting nearest cell location CHECKED THIS and its OKAY
   yPosCell=int(yPos(pMov)/dy) + 1
   zPosCell=int((zPos(pMov) + dz)/dz) + 1
   

   if(xPosCell>1 .and. xPosCell <nx .and. yPosCell>1 .and. yPosCell<ny .and. zPosCell>1 .and. zPosCell<nz) then  !CHECK THIS CONDITION
      sigUCell=sigU(xPosCell,yPosCell,zPosCell)
      sigVCell=sigV(xPosCell,yPosCell,zPosCell)
      sigWCell=sigW(xPosCell,yPosCell,zPosCell)




      Tau11=max( sigUCell **(2.0) ,1e-6) !To avoid divide by zero in lamda calc. below [sig**2 will be always +ve, so used +ve 1e-6]
      Tau22=max( sigVCell **(2.0) ,1e-6)
      Tau33=max( sigWCell **(2.0) ,1e-6)
      Tau21=0.1
      Tau31=0.1
      Tau32=0.1

      Tau12=Tau21 !Ref: Chap 8 Rodean P. 36 following Doob et al.
      Tau13=Tau31   
      Tau23=Tau32


      sigUCellPrev = sigU(xPosCell-1,yPosCell,zPosCell)! Check Exception for xPosCell=1
      sigVCellPrev = sigV(xPosCell-1,yPosCell,zPosCell)
      sigWCellPrev = sigW(xPosCell-1,yPosCell,zPosCell)

      sigUCellNext = sigU(xPosCell+1,yPosCell,zPosCell)
      sigVCellNext = sigV(xPosCell+1,yPosCell,zPosCell)
      sigWCellNext = sigW(xPosCell+1,yPosCell,zPosCell)

   
      if( (icellflag(xPosCell-1,yPosCell,zPosCell)+icellflag(xPosCell,yPosCell,zPosCell)+icellflag(xPosCell+1,yPosCell,zPosCell) ) == 3 .and. xPosCell>1 .and. xPosCell<nx ) then ! See Notes Below(1)

         DTau11Dx = (sigUCellPrev **(2.0) - sigUCellNext **(2.0) ) / (2.0*dx)                               
         DTau22Dx = (sigVCellPrev **(2.0) - sigVCellNext **(2.0) ) / (2.0*dx)      
         DTau33Dx = (sigWCellPrev **(2.0) - sigWCellNext **(2.0) ) / (2.0*dx)
         !DTau12Dx =                       
         !DTau13Dx =
         !DTau23Dx =
          
      elseif(icellflag(xPosCell-1,yPosCell,zPosCell)==0 .and. xPosCell>1 .and. xPosCell<nx ) then
         
         DTau11Dx = (sigUCell ** 2.0) / ( (.5*dx) * log((.5*dx)/wallRough) )  ! See Notes below(2)
         DTau22Dx = (sigVCell ** 2.0) / ( (.5*dx) * log((.5*dx)/wallRough) ) 
         DTau33Dx = (sigWCell ** 2.0) / ( (.5*dx) * log((.5*dx)/wallRough) ) 
         !DTau12Dx =
         !DTau13Dx =
         !DTau23Dx =
                  
      elseif(icellflag(xPosCell+1,yPosCell,zPosCell)==0 .and. xPosCell>1 .and. xPosCell<nx ) then
         
         DTau11Dx = -(sigUCell ** 2.0) / ( (.5*dx) * log((.5*dx)/wallRough) )  
         DTau22Dx = -(sigVCell ** 2.0) / ( (.5*dx) * log((.5*dx)/wallRough) ) 
         DTau33Dx = -(sigWCell ** 2.0) / ( (.5*dx) * log((.5*dx)/wallRough) )  
         !DTau12Dx =
         !DTau13Dx =
         !DTau23Dx =
         
      endif 

      DTau21Dx=DTau12Dx
      DTau31Dx=DTau13Dx
      DTau32Dx=DTau23Dx


   
      sigUCellPrev = sigU(xPosCell,yPosCell-1,zPosCell) ! Check Exception for yPosCell=1
      sigVCellPrev = sigV(xPosCell,yPosCell-1,zPosCell)
      sigWCellPrev = sigW(xPosCell,yPosCell-1,zPosCell)

      sigUCellNext = sigU(xPosCell,yPosCell+1,zPosCell)
      sigVCellNext = sigV(xPosCell,yPosCell+1,zPosCell)
      sigWCellNext = sigW(xPosCell,yPosCell+1,zPosCell)

   
      if( (icellflag(xPosCell,yPosCell-1,zPosCell)+icellflag(xPosCell,yPosCell,zPosCell)+icellflag(xPosCell,yPosCell+1,zPosCell) ) == 3 .and. yPosCell>1 .and. yPosCell<ny ) then ! See Notes Below(1)

         DTau11Dy = (sigUCellPrev **(2.0) - sigUCellNext **(2.0) ) / (2.0*dy)                               
         DTau22Dy = (sigVCellPrev **(2.0) - sigVCellNext **(2.0) ) / (2.0*dy)      
         DTau33Dy = (sigWCellPrev **(2.0) - sigWCellNext **(2.0) ) / (2.0*dy)
         !DTau12Dy =                       
         !DTau13Dy =
         !DTau23Dy =
          
      elseif(icellflag(xPosCell,yPosCell-1,zPosCell)==0 .and. yPosCell>1 .and. yPosCell<ny ) then
         
         DTau11Dy = (sigUCell ** 2.0) / ( (.5*dy) * log((.5*dy)/wallRough) )  ! See Notes below(2)
         DTau22Dy = (sigVCell ** 2.0) / ( (.5*dy) * log((.5*dy)/wallRough) ) 
         DTau33Dy = (sigWCell ** 2.0) / ( (.5*dy) * log((.5*dy)/wallRough) ) 
         !DTau12Dy =
         !DTau13Dy =
         !DTau23Dy =
                  
      elseif(icellflag(xPosCell,yPosCell+1,zPosCell)==0 .and. yPosCell>1 .and. yPosCell<ny ) then
         
         DTau11Dy = -(sigUCell ** 2.0) / ( (.5*dy) * log((.5*dy)/wallRough) )  
         DTau22Dy = -(sigVCell ** 2.0) / ( (.5*dy) * log((.5*dy)/wallRough) ) 
         DTau33Dy = -(sigWCell ** 2.0) / ( (.5*dy) * log((.5*dy)/wallRough) )  
         !DTau12Dy =
         !DTau13Dy =
         !DTau23Dy =
         
      endif

      DTau21Dy=DTau12Dy
      DTau31Dy=DTau13Dy
      DTau32Dy=DTau23Dy



      sigUCellPrev = sigU(xPosCell,yPosCell,zPosCell-1) ! Check Exception for zPosCell=1
      sigVCellPrev = sigV(xPosCell,yPosCell,zPosCell-1)
      sigWCellPrev = sigW(xPosCell,yPosCell,zPosCell-1)

      sigUCellNext = sigU(xPosCell,yPosCell,zPosCell+1)
      sigVCellNext = sigV(xPosCell,yPosCell,zPosCell+1)
      sigWCellNext = sigW(xPosCell,yPosCell,zPosCell+1)

      if( (icellflag(xPosCell,yPosCell,zPosCell-1)+icellflag(xPosCell,yPosCell,zPosCell)+icellflag(xPosCell,yPosCell,zPosCell+1) ) == 3 .and. zPosCell>1 .and. zPosCell<nz ) then ! See Notes Below(1)

         DTau11Dz = (sigUCellPrev **(2.0) - sigUCellNext **(2.0) ) / (2.0*dz)                               
         DTau22Dz = (sigVCellPrev **(2.0) - sigVCellNext **(2.0) ) / (2.0*dz)      
         DTau33Dz = (sigWCellPrev **(2.0) - sigWCellNext **(2.0) ) / (2.0*dz)
         !DTau12Dz =                       
         !DTau13Dz =
         !DTau23Dz =
          
      elseif(icellflag(xPosCell,yPosCell,zPosCell-1)==0 .and. zPosCell>1 .and. zPosCell<nz ) then
         
         DTau11Dz = (sigUCell ** 2.0) / ( (.5*dz) * log((.5*dz)/wallRough) )  ! See Notes below(2)
         DTau22Dz = (sigVCell ** 2.0) / ( (.5*dz) * log((.5*dz)/wallRough) ) 
         DTau33Dz = (sigWCell ** 2.0) / ( (.5*dz) * log((.5*dz)/wallRough) ) 
         !DTau12Dz =
         !DTau13Dz =
         !DTau23Dz =
                  
      elseif(icellflag(xPosCell,yPosCell,zPosCell+1)==0 .and. zPosCell>1 .and. zPosCell<nz ) then
         
         DTau11Dz = -(sigUCell ** 2.0) / ( (.5*dz) * log((.5*dz)/wallRough) )  
         DTau22Dz = -(sigVCell ** 2.0) / ( (.5*dz) * log((.5*dz)/wallRough) ) 
         DTau33Dz = -(sigWCell ** 2.0) / ( (.5*dz) * log((.5*dz)/wallRough) )  
         !DTau12Dz =
         !DTau13Dz =
         !DTau23Dz =
         
      endif

      DTau21Dz=DTau12Dz
      DTau31Dz=DTau13Dz
      DTau32Dz=DTau23Dz



      Lam11=(Tau11)**(-1.00)  ! Ref: Chap 8 Rodean P. 35 bottom
      Lam22=(Tau22)**(-1.00)  
      Lam33=(Tau33)**(-1.00)
      Lam12=(Tau12)**(-1.00)
      Lam13=(Tau13)**(-1.00)
      Lam21=(Tau21)**(-1.00)   
      Lam23=(Tau23)**(-1.00)
      Lam31=(Tau31)**(-1.00)
      Lam32=(Tau32)**(-1.00)

      Eps=0.001

                 
      uMem= - (Co*Eps/2.0) * ( Lam11*uPrimePrev(pMov) + Lam12*vPrimePrev(pMov) + Lam13*wPrimePrev(pMov) )

      uDri1= (uMeanInterp*DuDx(xPosCell,yPosCell,zPosCell) + vMeanInterp*DuDy(xPosCell,yPosCell,zPosCell) + wMeanInterp*DuDz(xPosCell,yPosCell,zPosCell)) + 0.5*(DTau11Dx +DTau12Dy +DTau13Dz) + (uPrimePrev(pMov)*DuDx(xPosCell,yPosCell,zPosCell) + vPrimePrev(pMov)*DuDy(xPosCell,yPosCell,zPosCell) + wPrimePrev(pMov)*DuDz(xPosCell,yPosCell,zPosCell))

      uDri2= (uPrimePrev(pMov)*Tau11/2.0) * (uMeanInterp*DLam11Dx + vMeanInterp*DLam11Dy + wMeanInterp*DLam11Dz)  ! 6th term starts here
      uDri3= (uPrimePrev(pMov)*Tau12/2.0) * (uMeanInterp*DLam12Dx + vMeanInterp*DLam12Dy + wMeanInterp*DLam12Dz)
      uDri4= (uPrimePrev(pMov)*Tau13/2.0) * (uMeanInterp*DLam13Dx + vMeanInterp*DLam13Dy + wMeanInterp*DLam13Dz)
      uDri5= (vPrimePrev(pMov)*Tau11/2.0) * (uMeanInterp*DLam21Dx + vMeanInterp*DLam21Dy + wMeanInterp*DLam21Dz)
      uDri6= (vPrimePrev(pMov)*Tau12/2.0) * (uMeanInterp*DLam22Dx + vMeanInterp*DLam22Dy + wMeanInterp*DLam22Dz)
      uDri7= (vPrimePrev(pMov)*Tau13/2.0) * (uMeanInterp*DLam23Dx + vMeanInterp*DLam23Dy + wMeanInterp*DLam23Dz)
      uDri8= (wPrimePrev(pMov)*Tau11/2.0) * (uMeanInterp*DLam31Dx + vMeanInterp*DLam31Dy + wMeanInterp*DLam31Dz)
      uDri9= (wPrimePrev(pMov)*Tau12/2.0) * (uMeanInterp*DLam32Dx + vMeanInterp*DLam32Dy + wMeanInterp*DLam32Dz)
      uDri10=(wPrimePrev(pMov)*Tau13/2.0) * (uMeanInterp*DLam33Dx + vMeanInterp*DLam33Dy + wMeanInterp*DLam33Dz)
                  
      uDri11=(uPrimePrev(pMov)*Tau11/2.0) * (uPrimePrev(pMov)*DLam11Dx + vPrimePrev(pMov)*DLam21Dx + wPrimePrev(pMov)*DLam31Dx)     ! 7th term starts here
      uDri12=(uPrimePrev(pMov)*Tau12/2.0) * (uPrimePrev(pMov)*DLam12Dx + vPrimePrev(pMov)*DLam22Dx + wPrimePrev(pMov)*DLam32Dx)
      uDri13=(uPrimePrev(pMov)*Tau13/2.0) * (uPrimePrev(pMov)*DLam13Dx + vPrimePrev(pMov)*DLam23Dx + wPrimePrev(pMov)*DLam33Dx)
      uDri14=(vPrimePrev(pMov)*Tau11/2.0) * (uPrimePrev(pMov)*DLam11Dy + vPrimePrev(pMov)*DLam21Dy + wPrimePrev(pMov)*DLam31Dy)
      uDri15=(vPrimePrev(pMov)*Tau12/2.0) * (uPrimePrev(pMov)*DLam12Dy + vPrimePrev(pMov)*DLam22Dy + wPrimePrev(pMov)*DLam32Dy)
      uDri16=(vPrimePrev(pMov)*Tau13/2.0) * (uPrimePrev(pMov)*DLam13Dy + vPrimePrev(pMov)*DLam23Dy + wPrimePrev(pMov)*DLam33Dy)
      uDri17=(wPrimePrev(pMov)*Tau11/2.0) * (uPrimePrev(pMov)*DLam11Dz + vPrimePrev(pMov)*DLam21Dz + wPrimePrev(pMov)*DLam31Dz)
      uDri18=(wPrimePrev(pMov)*Tau12/2.0) * (uPrimePrev(pMov)*DLam12Dz + vPrimePrev(pMov)*DLam22Dz + wPrimePrev(pMov)*DLam32Dz)
      uDri19=(wPrimePrev(pMov)*Tau13/2.0) * (uPrimePrev(pMov)*DLam13Dz + vPrimePrev(pMov)*DLam23Dz + wPrimePrev(pMov)*DLam33Dz)
           
      uDri=uDri1 - (uDri2+uDri3+uDri4+uDri5+uDri6+uDri7+uDri8+uDri9+uDri10+uDri11+uDri12) - (uDri13+uDri14+uDri15+uDri16+uDri17+uDri18+uDri19)

      call RanGauss(x)
      uRan=sqrt(Co*Eps*timeStep)*x
           
   
      vMem= - (Co*Eps/2.0) * ( Lam21*uPrimePrev(pMov) + Lam22*vPrimePrev(pMov) + Lam23*wPrimePrev(pMov) )

      vDri1= (uMeanInterp*DvDx(xPosCell,yPosCell,zPosCell) + vMeanInterp*DvDy(xPosCell,yPosCell,zPosCell) + wMeanInterp*DvDz(xPosCell,yPosCell,zPosCell)) + 0.5*(DTau21Dx +DTau22Dy +DTau23Dz) + (uPrimePrev(pMov)*DvDx(xPosCell,yPosCell,zPosCell) + vPrimePrev(pMov)*DvDy(xPosCell,yPosCell,zPosCell) + wPrimePrev(pMov)*DvDz(xPosCell,yPosCell,zPosCell))
   
      vDri2= (uPrimePrev(pMov)*Tau21/2.0) * (uMeanInterp*DLam11Dx + vMeanInterp*DLam11Dy + wMeanInterp*DLam11Dz)  ! 6th term starts here
      vDri3= (uPrimePrev(pMov)*Tau22/2.0) * (uMeanInterp*DLam12Dx + vMeanInterp*DLam12Dy + wMeanInterp*DLam12Dz)
      vDri4= (uPrimePrev(pMov)*Tau23/2.0) * (uMeanInterp*DLam13Dx + vMeanInterp*DLam13Dy + wMeanInterp*DLam13Dz)
      vDri5= (vPrimePrev(pMov)*Tau21/2.0) * (uMeanInterp*DLam21Dx + vMeanInterp*DLam21Dy + wMeanInterp*DLam21Dz)
      vDri6= (vPrimePrev(pMov)*Tau22/2.0) * (uMeanInterp*DLam22Dx + vMeanInterp*DLam22Dy + wMeanInterp*DLam22Dz)
      vDri7= (vPrimePrev(pMov)*Tau23/2.0) * (uMeanInterp*DLam23Dx + vMeanInterp*DLam23Dy + wMeanInterp*DLam23Dz)
      vDri8= (wPrimePrev(pMov)*Tau21/2.0) * (uMeanInterp*DLam31Dx + vMeanInterp*DLam31Dy + wMeanInterp*DLam31Dz)
      vDri9= (wPrimePrev(pMov)*Tau22/2.0) * (uMeanInterp*DLam32Dx + vMeanInterp*DLam32Dy + wMeanInterp*DLam32Dz)
      vDri10=(wPrimePrev(pMov)*Tau23/2.0) * (uMeanInterp*DLam33Dx + vMeanInterp*DLam33Dy + wMeanInterp*DLam33Dz)
   
      vDri11=(uPrimePrev(pMov)*Tau21/2.0) * (uPrimePrev(pMov)*DLam11Dx + vPrimePrev(pMov)*DLam21Dx + wPrimePrev(pMov)*DLam31Dx)
      vDri12=(uPrimePrev(pMov)*Tau22/2.0) * (uPrimePrev(pMov)*DLam12Dx + vPrimePrev(pMov)*DLam22Dx + wPrimePrev(pMov)*DLam32Dx)
      vDri13=(uPrimePrev(pMov)*Tau23/2.0) * (uPrimePrev(pMov)*DLam13Dx + vPrimePrev(pMov)*DLam23Dx + wPrimePrev(pMov)*DLam33Dx)
      vDri14=(vPrimePrev(pMov)*Tau21/2.0) * (uPrimePrev(pMov)*DLam11Dy + vPrimePrev(pMov)*DLam21Dy + wPrimePrev(pMov)*DLam31Dy)
      vDri15=(vPrimePrev(pMov)*Tau22/2.0) * (uPrimePrev(pMov)*DLam12Dy + vPrimePrev(pMov)*DLam22Dy + wPrimePrev(pMov)*DLam32Dy)
      vDri16=(vPrimePrev(pMov)*Tau23/2.0) * (uPrimePrev(pMov)*DLam13Dy + vPrimePrev(pMov)*DLam23Dy + wPrimePrev(pMov)*DLam33Dy)
      vDri17=(wPrimePrev(pMov)*Tau21/2.0) * (uPrimePrev(pMov)*DLam11Dz + vPrimePrev(pMov)*DLam21Dz + wPrimePrev(pMov)*DLam31Dz)
      vDri18=(wPrimePrev(pMov)*Tau22/2.0) * (uPrimePrev(pMov)*DLam12Dz + vPrimePrev(pMov)*DLam22Dz + wPrimePrev(pMov)*DLam32Dz)
      vDri19=(wPrimePrev(pMov)*Tau23/2.0) * (uPrimePrev(pMov)*DLam13Dz + vPrimePrev(pMov)*DLam23Dz + wPrimePrev(pMov)*DLam33Dz)              
   
      vDri=vDri1 - (vDri2+vDri3+vDri4+vDri5+vDri6+vDri7+vDri8+vDri9+vDri10+vDri11+vDri12) - (vDri13+vDri14+vDri15+vDri16+vDri17+vDri18+vDri19)
   
      call RanGauss(x)
      vRan=sqrt(Co*Eps*timeStep)*x
   
           
      wMem= - (Co*Eps/2.0) * ( Lam31*uPrimePrev(pMov) + Lam32*vPrimePrev(pMov) + Lam33*wPrimePrev(pMov) )

      wDri1= (uMeanInterp*DwDx(xPosCell,yPosCell,zPosCell) + vMeanInterp*DwDy(xPosCell,yPosCell,zPosCell) + wMeanInterp*DwDz(xPosCell,yPosCell,zPosCell)) + 0.5*(DTau31Dx +DTau32Dy +DTau33Dz) + (uPrimePrev(pMov)*DwDx(xPosCell,yPosCell,zPosCell) + vPrimePrev(pMov)*DwDy(xPosCell,yPosCell,zPosCell) + wPrimePrev(pMov)*DwDz(xPosCell,yPosCell,zPosCell))
   
      wDri2= (uPrimePrev(pMov)*Tau31/2.0) * (uMeanInterp*DLam11Dx + vMeanInterp*DLam11Dy + wMeanInterp*DLam11Dz)  ! 6th term starts here
      wDri3= (uPrimePrev(pMov)*Tau32/2.0) * (uMeanInterp*DLam12Dx + vMeanInterp*DLam12Dy + wMeanInterp*DLam12Dz)
      wDri4= (uPrimePrev(pMov)*Tau33/2.0) * (uMeanInterp*DLam13Dx + vMeanInterp*DLam13Dy + wMeanInterp*DLam13Dz)
      wDri5= (vPrimePrev(pMov)*Tau31/2.0) * (uMeanInterp*DLam21Dx + vMeanInterp*DLam21Dy + wMeanInterp*DLam21Dz)
      wDri6= (vPrimePrev(pMov)*Tau32/2.0) * (uMeanInterp*DLam22Dx + vMeanInterp*DLam22Dy + wMeanInterp*DLam22Dz)
      wDri7= (vPrimePrev(pMov)*Tau33/2.0) * (uMeanInterp*DLam23Dx + vMeanInterp*DLam23Dy + wMeanInterp*DLam23Dz)
      wDri8= (wPrimePrev(pMov)*Tau31/2.0) * (uMeanInterp*DLam31Dx + vMeanInterp*DLam31Dy + wMeanInterp*DLam31Dz)
      wDri9= (wPrimePrev(pMov)*Tau32/2.0) * (uMeanInterp*DLam32Dx + vMeanInterp*DLam32Dy + wMeanInterp*DLam32Dz)
      wDri10=(wPrimePrev(pMov)*Tau33/2.0) * (uMeanInterp*DLam33Dx + vMeanInterp*DLam33Dy + wMeanInterp*DLam33Dz)

      wDri11=(uPrimePrev(pMov)*Tau31/2.0) * (uPrimePrev(pMov)*DLam11Dx + vPrimePrev(pMov)*DLam21Dx + wPrimePrev(pMov)*DLam31Dx)
      wDri12=(uPrimePrev(pMov)*Tau32/2.0) * (uPrimePrev(pMov)*DLam12Dx + vPrimePrev(pMov)*DLam22Dx + wPrimePrev(pMov)*DLam32Dx)
      wDri13=(uPrimePrev(pMov)*Tau33/2.0) * (uPrimePrev(pMov)*DLam13Dx + vPrimePrev(pMov)*DLam23Dx + wPrimePrev(pMov)*DLam33Dx)
      wDri14=(vPrimePrev(pMov)*Tau31/2.0) * (uPrimePrev(pMov)*DLam11Dy + vPrimePrev(pMov)*DLam21Dy + wPrimePrev(pMov)*DLam31Dy)
      wDri15=(vPrimePrev(pMov)*Tau32/2.0) * (uPrimePrev(pMov)*DLam12Dy + vPrimePrev(pMov)*DLam22Dy + wPrimePrev(pMov)*DLam32Dy)
      wDri16=(vPrimePrev(pMov)*Tau33/2.0) * (uPrimePrev(pMov)*DLam13Dy + vPrimePrev(pMov)*DLam23Dy + wPrimePrev(pMov)*DLam33Dy)
      wDri17=(wPrimePrev(pMov)*Tau31/2.0) * (uPrimePrev(pMov)*DLam11Dz + vPrimePrev(pMov)*DLam21Dz + wPrimePrev(pMov)*DLam31Dz)
      wDri18=(wPrimePrev(pMov)*Tau32/2.0) * (uPrimePrev(pMov)*DLam12Dz + vPrimePrev(pMov)*DLam22Dz + wPrimePrev(pMov)*DLam32Dz)
      wDri19=(wPrimePrev(pMov)*Tau33/2.0) * (uPrimePrev(pMov)*DLam13Dz + vPrimePrev(pMov)*DLam23Dz + wPrimePrev(pMov)*DLam33Dz)
   
      wDri=wDri1 - (wDri2+wDri3+wDri4+wDri5+wDri6+wDri7+wDri8+wDri9+wDri10+wDri11+wDri12) - (wDri13+wDri14+wDri15+wDri16+wDri17+wDri18+wDri19)

      call RanGauss(x)
      wRan=sqrt(Co*Eps*timeStep)*x
   
      du=uMem+uDri+uRan
      dv=vMem+vDri+vRan
      dw=wMem+wDri+wRan

      uPrimeCurr(pMov)=uPrimePrev(pMov) + du 
      vPrimeCurr(pMov)=vPrimePrev(pMov) + dv
      wPrimeCurr(pMov)=wPrimePrev(pMov) + dw

      uPrime=(uPrimeCurr(pMov) + uPrimePrev(pMov))/2.0 
      vPrime=(vPrimeCurr(pMov) + vPrimePrev(pMov))/2.0
      wPrime=(wPrimeCurr(pMov) + wPrimePrev(pMov))/2.0
      
      ! Calculating interpolated velocity
      dissXCenter=xPos(pMov)-.5*dx-dx*float(xPosCell-1)      !Distance of particle from cell center in X  CHECKED this and its okay
      dissYCenter=yPos(pMov)-.5*dy-dy*float(yPosCell-1)      !Distance of particle from cell center in Y
      dissZCenter=zPos(pMov)-.5*dz-dz*float(zPosCell-2)      !Distance of particle from cell center in Z

      call windInterp

      disIncrX=uMeanInterp*timeStep + uPrime*timeStep   
      disIncrY=vMeanInterp*timeStep + vPrime*timeStep
      disIncrZ=wMeanInterp*timeStep + wPrime*timeStep 

      !if(disIncrX > dx) then
      !  dt=dt/2.0
      !   incr=incr+1  ! JUST THINKING THIS IDEA......
      !endif
   

      xPos(pMov)=xPos(pMov)+disIncrX 
      yPos(pMov)=yPos(pMov)+disIncrY
      zPos(pMov)=zPos(pMov)+disIncrZ
   
      if(xPos(pMov)<=0.0 .or. xPos(pMov)>=xDomainLen .or. yPos(pMov)<=0.0 .or. yPos(pMov)>=yDomainLen .or. abs(zPos(pMov))>=zDomainLen)  then  ! moving particles to a point(0,0,-0.5*dz) if they leave domain, Used abs(zPos()) as sometime I get huge -ve value exceeding zDomainLen
         xPos(pMov)=0.0
         yPos(pMov)=0.0
         zPos(pMov)=-0.5*dz
         pOutDomain(pMov)=pMov                              ! Marking particles which are out of domain
      endif

      xPosCell=int(xPos(pMov)/dx) + 1                        ! getting nearest cell location CHECK THIS**
      yPosCell=int(yPos(pMov)/dy) + 1
      zPosCell=int((zPos(pMov) + dz)/dz) + 1
   
      if(zPosCell<=0.0 .or. icellflag(xPosCell,yPosCell,zPosCell)==1) call Reflection  ! Reflection Required
   
      write(101,*) t,pMov,xPos(pMov),yPos(pMov),zPos(pmov)            !remove this****
  
   endif  ! CHECK THIS CONDITION

   end




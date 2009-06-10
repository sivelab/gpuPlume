      subroutine windInterp
         use Datamodule
         use disperseModule

         if(dissXCenter.eq.0.and.dissYCenter.eq.0.and.dissZCenter.eq.0)then
            uMeanInterp=u(xPosCell,yPosCell,zPosCell)
            vMeanInterp=v(xPosCell,yPosCell,zPosCell)
            wMeanInterp=w(xPosCell,yPosCell,zPosCell)
            return
         endif
   !     if(closestWall(xPosCell,yPosCell,zPosCell).lt.dx)then
   !        uMeanInterp=u(xPosCell,yPosCell,zPosCell)*log((closestWall(xPosCell,yPosCell,zPosCell)+wallRough)/wallRough)/log((dx/2.+wallRough)/(wallRough))
   !        vMeanInterp=v(xPosCell,yPosCell,zPosCell)*log((closestWall(xPosCell,yPosCell,zPosCell)+wallRough)/wallRough)/log((dx/2.+wallRough)/(wallRough))
   !        wMeanInterp=w(xPosCell,yPosCell,zPosCell)*log((closestWall(xPosCell,yPosCell,zPosCell)+wallRough)/wallRough)/log((dx/2.+wallRough)/(wallRough))
   !        return
   !     endif
   !
   !     if(closestWall(xPosCell,yPosCell,zPosCell).lt.dy)then
   !        uMeanInterp=u(xPosCell,yPosCell,zPosCell)*log((closestWall(xPosCell,yPosCell,zPosCell)+wallRough)/wallRough)/log((dy/2.+wallRough)/(wallRough))
   !        vMeanInterp=v(xPosCell,yPosCell,zPosCell)*log((closestWall(xPosCell,yPosCell,zPosCell)+wallRough)/wallRough)/log((dy/2.+wallRough)/(wallRough))
   !        wMeanInterp=w(xPosCell,yPosCell,zPosCell)*log((closestWall(xPosCell,yPosCell,zPosCell)+wallRough)/wallRough)/log((dy/2.+wallRough)/(wallRough))
   !        return
   !     endif
   !     if(closestWall(xPosCell,yPosCell,zPosCell).lt.dz)then
   !        uMeanInterp=u(xPosCell,yPosCell,zPosCell)*log((closestWall(xPosCell,yPosCell,zPosCell)+wallRough)/wallRough)/log((dz/2.+wallRough)/(wallRough))
   !        vMeanInterp=v(xPosCell,yPosCell,zPosCell)*log((closestWall(xPosCell,yPosCell,zPosCell)+wallRough)/wallRough)/log((dz/2.+wallRough)/(wallRough))
   !        wMeanInterp=w(xPosCell,yPosCell,zPosCell)*log((closestWall(xPosCell,yPosCell,zPosCell)+wallRough)/wallRough)/log((dz/2.+wallRough)/(wallRough))
   !        return
   !     endif
         dsum=0.
         usum=0
         vsum=0.
         wsum=0.
         i1=xPosCell-1
         i2=xPosCell+1
         i1=max(1,i1)
         i2=min(nx-1,i2)
         j1=yPosCell-1
         j2=yPosCell+1
         j1=max(1,j1)
         j2=min(ny-1,j2)
         k1=zPosCell-1
         k1=max(2,k1)
         k2=zPosCell+1
         k2=min(nz-1,k2)
lp017:   do k=k1,k2
lp016:      do j=j1,j2
lp015:         do i=i1,i2
                  up=u(i,j,k)
                  vp=v(i,j,k)
                  wp=w(i,j,k)
                  dxx=dx*float(i-xPosCell)-dissXCenter           
                  dyy=dy*float(j-yPosCell)-dissYCenter
                  dzz=dz*float(k-zPosCell)-dissZCenter
                  if(i.lt.xPosCell.and.disNegXSurf(xPosCell,yPosCell,zPosCell).lt.abs(dxx))then
                     dxx=disNegXSurf(xPosCell,yPosCell,zPosCell)+dissXCenter
                     up=0.
                     vp=0.
                     wp=0.
                     if(j.eq.yPosCell.and.k.eq.zPosCell)then
! Here we consider distance perpendicular to the wall
                        dyy=0.
                        dzz=0.
                     endif
                  endif
                  if(i.gt.xPosCell.and.disPosXSurf(xPosCell,yPosCell,zPosCell).lt.dxx)then
                     dxx=disPosXSurf(xPosCell,yPosCell,zPosCell)-dissXCenter
                     if(j.eq.yPosCell.and.k.eq.zPosCell)then
! Here we consider distance perpendicular to the wall
                        dyy=0.
                        dzz=0.
                     endif
                     up=0.
                     vp=0.
                     wp=0.
                  endif
                  if(j.lt.yPosCell.and.disNegYSurf(xPosCell,yPosCell,zPosCell).lt.abs(dyy))then 
                     dyy=disNegYSurf(xPosCell,yPosCell,zPosCell)+dissYCenter 
                     if(i.eq.xPosCell.and.k.eq.zPosCell)then
! Here we consider distance perpendicular to the wall
                        dxx=0.
                        dzz=0.
                     endif
                     up=0. 
                     vp=0. 
                     wp=0. 
                  endif 
                  if(j.gt.yPosCell.and.disPosYSurf(xPosCell,yPosCell,zPosCell).lt.dyy)then 
                     dyy=disPosYSurf(xPosCell,yPosCell,zPosCell)-dissYCenter 
                     if(i.eq.xPosCell.and.k.eq.zPosCell)then
! Here we consider distance perpendicular to the wall
                        dxx=0.
                        dzz=0.
                     endif
                     up=0.
                     vp=0.
                     wp=0.
                  endif
                  if(k.lt.zPosCell.and.disBotSurf(xPosCell,yPosCell,zPosCell).lt.abs(dzz))then
                     dzz=disBotSurf(xPosCell,yPosCell,zPosCell)+dissZCenter
                     if(i.eq.xPosCell.and.j.eq.yPosCell)then
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
                  if(zPosCell.eq.2)then
                     dzz=dz/2.+dissZCenter
                     if(i.eq.xPosCell.and.j.eq.yPosCell)then
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
         uMeanInterp=usum/dsum
         vMeanInterp=vsum/dsum
         wMeanInterp=wsum/dsum
         return
      end
      
      
      !
      !
      !
      !
      !
      !
      !
      !
      !
      !
      !
      !
      !
      !
      !
      !
      !
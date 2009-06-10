       subroutine Turbulence
          
          use DataModule
          implicit none 

          real closestWall,closestWallX,closestWallY,closestWallZ,Du,Dv,Dw,cSigUVnKar,cSigVVnKar,cSigWVnKar,ustarLoc
          real cellFlagNegK(nz),cellFlagPosK(nz),iNegVal(1),iPosVal(1),jNegVal(1),jPosVal(1),kNegVal(1),kPosVal(1)          
          real cellFlagNegI(nx),cellFlagPosI(nx),cellFlagNegJ(ny),cellFlagPosJ(ny)

          open(1001,file='OutputFiles\testTurb.dat',status='unknown')
          
          cSigUVnKar = VnKarman * 2.5
          cSigVVnKar = VnKarman * 2.0
          cSigWVnKar = VnKarman * 1.3
        


          do k=1,nz
             do j=1,ny
                do i=1,nx
                   !disBotSurf(i,j,k)=0
                   if(icellflag(i,j,k) /= 0) then    ! For All free atmos QU cells
        
                      !**Find out if following arrays are important to store or not OTHERWISE use variables instead of arrays ***
                      cellFlagNegI(1:i)=icellflag(1:i,j,k)    ! getting val of celltype from 1 to (i,j,k) in X Direction
                      cellFlagNegI(1:i)=cellFlagNegI(i:1:-1)  !reversing the array as we need to go from (i,j,k) into negative direction to find a surface
                      iNegVal=nx                              !iNegVal if no building in neg X direction
                      if(any(cellFlagNegI(1:i)==0))iNegVal=minloc(cellFlagNegI(1:i))-1  ! Finding surface(icellflga==0) in negative X- direction
                      disNegXSurf(i,j,k)=iNegVal(1)           ! Storing the answer


                      cellFlagPosI(i:nx)=icellflag(i:nx,j,k)  ! getting val of celltype from (i,j,k) to nx in X Direction                 
                      iPosVal=nx                              !iNegVal if no building in neg X direction
                      if(any(cellFlagPosI(i:nx)==0))iPosVal=minloc(cellFlagNegI(i:nx))-1  ! Finding surface(icellflga==0) in negative X- direction
                      disPosXSurf(i,j,k)=iPosVal(1)           ! Storing the answer
                    

                      cellFlagNegJ(1:j)=icellflag(i,1:j,k)    ! getting val of celltype from 1 to (i,j,k) in X Direction
                      cellFlagNegJ(1:j)=cellFlagNegJ(j:1:-1)  !reversing the array as we need to go from (i,j,k) into negative direction to find a surface
                      jNegVal=ny                              !iNegVal if no building in neg X direction
                      if(any(cellFlagNegJ(1:j)==0))jNegVal=minloc(cellFlagNegJ(1:j))-1  ! Finding surface(icellflga==0) in negative Y- direction
                      disNegYSurf(i,j,k)=jNegVal(1)           ! Storing the answer


                      cellFlagPosJ(j:ny)=icellflag(i,j:ny,k)  ! getting val of celltype from (i,j,k) to nx in X Direction                 
                      jPosVal=ny                              !iNegVal if no building in neg X direction
                      if(any(cellFlagPosJ(j:ny)==0))jPosVal=minloc(cellFlagPosJ(j:ny))-1  ! Finding surface(icellflga==0) in negative X- direction
                      disPosYSurf(i,j,k)=jPosVal(1)
                    
                      cellFlagNegK(1:k)=icellflag(i,j,1:k)    ! getting val of celltype from 1 to (i,j,k) in X Direction
                      cellFlagNegK(1:k)=cellFlagNegK(k:1:-1)  !reversing the array as we need to go from (i,j,k) into negative direction to find a surface
                      kNegVal=nz                              !iNegVal if no building in neg X direction
                      if(any(cellFlagNegK(1:k)==0))kNegVal=minloc(cellFlagNegK(1:k))-1  ! Finding surface(icellflga==0) in negative Y- direction
                      disBotSurf(i,j,k)=kNegVal(1)            ! Storing the answer

                      cellFlagPosK(k:nz)=icellflag(i,j,k:nz)  ! getting val of celltype from (i,j,k) to nx in X Direction                 
                      kPosVal=nz                              !iNegVal if no building in neg X direction
                      if(any(cellFlagPosK(k:nz)==0))kPosVal=minloc(cellFlagPosK(k:nz))-1  ! Finding surface(icellflga==0) in negative X- direction
                      disTopSurf(i,j,k)=kPosVal(1)                    
                    
                    
                      closestWallX = dx * min(disNegXSurf(i,j,k),disPosXSurf(i,j,k)) ! converting in meters
                      closestWallY = dy * min(disPosYSurf(i,j,k),disNegYSurf(i,j,k))
                      closestWallZ = dz * min(disTopSurf(i,j,k),disBotSurf(i,j,k))
                      
                      closestWall=min(closestWallX,closestWallY,closestWallZ)
                                                                 
                      !Note that 'k' will always be greater than 2 as we are using icellflag/=0 condition.

                      !For DuDx,DvDx,DwDx                      
                      if( (icellflag(i-1,j,k)+icellflag(i,j,k)+icellflag(i+1,j,k) ) == 3 .and. i>1 .and. i<nx ) then ! See Notes Below(1)
                         DuDx(i,j,k) = ( u(i+1,j,k) - u(i-1,j,k) )/(2.0*dx)
                         DvDx(i,j,k) = ( v(i+1,j,k) - v(i-1,j,k) )/(2.0*dx)                      
                         DwDx(i,j,k) = ( w(i+1,j,k) - w(i-1,j,k) )/(2.0*dx)
                          
                      elseif(icellflag(i-1,j,k)==0 .and. i>1 .and. i<nx ) then
                         DuDx(i,j,k) =u(i,j,k)/ ( (.5*dx) * log((.5*dx)/wallRough) )  ! See Notes below(2)
                         DvDx(i,j,k) =v(i,j,k)/ ( (.5*dx) * log((.5*dx)/wallRough) ) 
                         DwDx(i,j,k) =w(i,j,k)/ ( (.5*dx) * log((.5*dx)/wallRough) )  

                      elseif(icellflag(i+1,j,k)==0 .and. i>1 .and. i<nx ) then
                         DuDx(i,j,k) =-u(i,j,k)/ ( (.5*dx) * log((.5*dx)/wallRough) )  
                         DvDx(i,j,k) =-v(i,j,k)/ ( (.5*dx) * log((.5*dx)/wallRough) ) 
                         DwDx(i,j,k) =-w(i,j,k)/ ( (.5*dx) * log((.5*dx)/wallRough) )  
                      
                      endif
                      
                            
                      !For DuDy,DvDy,DwDy
                      if( (icellflag(i,j-1,k)+icellflag(i,j,k)+icellflag(i,j+1,k) ) == 3 .and. j>1 .and. j<ny) then
                         DuDy(i,j,k) = ( u(i,j+1,k) - u(i,j-1,k) )/(2.0*dy)
                         DvDy(i,j,k) = ( v(i,j+1,k) - v(i,j-1,k) )/(2.0*dy)
                         DwDy(i,j,k) = ( w(i,j+1,k) - w(i,j-1,k) )/(2.0*dy)
                             
                      elseif(icellflag(i,j-1,k)==0 .and. j>1 .and. j<ny ) then
                         DuDy(i,j,k) =u(i,j,k)/ ( (.5*dy) * log((.5*dy)/wallRough) )  ! See Notes below(2)
                         DvDy(i,j,k) =v(i,j,k)/ ( (.5*dy) * log((.5*dy)/wallRough) ) 
                         DwDy(i,j,k) =w(i,j,k)/ ( (.5*dy) * log((.5*dy)/wallRough) )

                      elseif(icellflag(i,j+1,k)==0 .and. j>1 .and. j<ny ) then
                         DuDy(i,j,k) =-u(i,j,k)/ ( (.5*dy) * log((.5*dy)/wallRough) )  
                         DvDy(i,j,k) =-v(i,j,k)/ ( (.5*dy) * log((.5*dy)/wallRough) ) 
                         DwDy(i,j,k) =-w(i,j,k)/ ( (.5*dy) * log((.5*dy)/wallRough) )  
                               
                      endif
                              
                                                     
                      !For DuDz,DvDz,DwDz, k is always greater than 1
                      if( (icellflag(i,j,k-1)+icellflag(i,j,k)+icellflag(i,j,k+1) ) == 3 .and. k<nz) then
                         DuDz(i,j,k) = ( u(i,j,k+1) - u(i,j,k-1) )/(2.0*dz)
                         DvDz(i,j,k) = ( v(i,j,k+1) - v(i,j,k-1) )/(2.0*dz)
                         DwDz(i,j,k) = ( w(i,j,k+1) - w(i,j,k-1) )/(2.0*dz)
                                        
                      elseif(icellflag(i,j,k-1)==0 .and. k<nz ) then
                         DuDz(i,j,k) =u(i,j,k)/ ( (.5*dz) * log((.5*dz)/wallRough) )  ! See Notes below(2)
                         DvDz(i,j,k) =v(i,j,k)/ ( (.5*dz) * log((.5*dz)/wallRough) ) 
                         DwDz(i,j,k) =w(i,j,k)/ ( (.5*dz) * log((.5*dz)/wallRough) ) 

                      elseif(icellflag(i,j,k+1)==0 .and. k<nz ) then
                         DuDz(i,j,k) =-u(i,j,k)/ ( (.5*dz) * log((.5*dz)/wallRough) ) 
                         DvDz(i,j,k) =-v(i,j,k)/ ( (.5*dz) * log((.5*dz)/wallRough) ) 
                         DwDz(i,j,k) =-w(i,j,k)/ ( (.5*dz) * log((.5*dz)/wallRough) )                                   

                      endif
					  ustarLoc= ( DuDx(i,j,k)+DuDy(i,j,k)+DuDz(i,j,k) + DvDx(i,j,k)+DvDy(i,j,k)+DvDz(i,j,k) + DwDx(i,j,k)+DwDy(i,j,k)+DwDz(i,j,k) ) * closestWall
                      sigU(i,j,k) = cSigUVnKar * ustarLoc
                      sigV(i,j,k) = cSigVVnKar * ustarLoc
                      sigW(i,j,k) = cSigWVnKar * ustarLoc

                      write(1001,'(3i6,13e12.3)')i,j,k,DuDx(i,j,k),DuDy(i,j,k),DuDz(i,j,k),DvDx(i,j,k),DvDy(i,j,k),DvDz(i,j,k),DwDx(i,j,k),DwDy(i,j,k),DwDz(i,j,k),sigU(i,j,k),sigV(i,j,k),sigW(i,j,k),closestWall
                    
                 endif
              enddo
           enddo
        enddo
        

        end






        ! 1. We excluded the boundary cells (i=nx,i=1,j=1,j=ny,k=nz) and equated their values to be zero as done in QP by MWD
        !    [icellflag(-1)+icellflag()+icellflag(+1) == 3] condition checks if cell before and after (i,j,k) is free atmos or not 
        ! 2. Log law derived from : U = (ustar/vonKarman) * Log(z/znaut)   and   (Du/Dz) = ustar / (vonKarman * z); MDW used U = (ustar/vonKarman) * Log((z+znaut)/znaut)   and   (Du/Dz) = ustar / (vonKarman * (z+znaut)) 
        !
             
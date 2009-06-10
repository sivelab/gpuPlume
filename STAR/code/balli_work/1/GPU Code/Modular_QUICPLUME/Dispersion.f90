        subroutine Dispersion
        use DataModule
        use disperseModule
        Implicit None   
        
        integer n
        integer p,pIni,pConc,nParTstep,firstPar,lastPar,nTimeSteps,strPar,endPar,time
        integer xConcBoxInd,yConcBoxInd,zConcBoxInd

        real, parameter:: pi=3.141592653
        real x,phiRan,thetaRan,xConcBoxSize,yConcBoxSize,zConcBoxSize,ConcConst,boxVolume,timeForAvging
        

        allocate(xPos(numPar),yPos(numPar),zPos(numPar),uPrimeCurr(numPar),vPrimeCurr(numPar),wPrimeCurr(numPar),uPrimePrev(numPar),vPrimePrev(numPar),wPrimePrev(numPar))
        allocate(pOutDomain(numPar),strTimep(numPar))
        allocate(Conc(numXConcBox,numYConcBox,numZConcBox))
        allocate(xConcBoxCen(numXConcBox),yConcBoxCen(numYConcBox),zConcBoxCen(numZConcBox))

        open(101,file='Particles.dat',status='unknown')   !remove this****

        xDomainLen=nx*dx
        yDomainLen=ny*dy
        zDomainLen=nz*dz-dz  ! Subtracted dz as 'z' starts at -0.5*dz, but nz*dz doesn't include that cell underneath and give higher value in meters, so we substract one cell
        pOutDomain=0
        strTimep=-1 
        timeForAvging=0                                                 !To get rid of junk values, start time CANNOT be negative , so its safe to use -1

             
        !Concentration Boxes or grid calculations
        Conc=0                                                         ! setting initail value of concentration
        xConcBoxSize=(XConcBoxUpLmt-XConcBoxLowLmt)/float(numXConcBox)
        yConcBoxSize=(YConcBoxUpLmt-YConcBoxLowLmt)/float(numYConcBox)
        zConcBoxSize=(ZConcBoxUpLmt-ZConcBoxLowLmt)/float(numZConcBox)

     
        
        xConcBoxCen(1:numXConcBox) = (/ (XConcBoxLowLmt + (xConcBoxSize/2) + (n-1)*xConcBoxSize, n = 1, numXConcBox) /)     ! implied Do loop
        yConcBoxCen(1:numYConcBox) = (/ (YConcBoxLowLmt + (yConcBoxSize/2) + (n-1)*yConcBoxSize, n = 1, numYConcBox) /)     ! implied Do loop
        zConcBoxCen(1:numZConcBox) = (/ (ZConcBoxLowLmt + (zConcBoxSize/2) + (n-1)*zConcBoxSize, n = 1, numZConcBox) /)     ! implied Do loop

        boxVolume=xConcBoxSize*yConcBoxSize*zConcBoxSize

        ConcConst=(timeStep*TotRel)/(boxVolume*ConAvgTime*float(numPar))       !concentration constant used to calculate concentration at the center of boxes

        !Concentration Boxes or grid calculations Ends
                 
             
          
        

        if(SrcType==1) then                                          !Point Source
        !using spherical polar coordinates to define initial positions of particles (http://mathworld.wolfram.com/SphericalCoordinates.html)
           do p=1,numPar                                              
              call random_number(x)                                  !Uniform random number between zero to 1
              thetaRan=2.0*pi*x                                      !Ranges from 0 to 2*pi
              call random_number(x)
              phiRan=2.0*pi*x
                  
              xPos(p)=xCordSrc+ (radSrc*cos(thetaRan)*sin(phiRan))   !Defining initial positions ,converting polar to cartesian
              yPos(p)=yCordSrc+ (radSrc*sin(thetaRan)*sin(phiRan))
              zPos(p)=zCordSrc+ (radSrc*cos(phiRan))
           enddo
        endif
     

        if(RelType==2) then                                          !Continous release
           nParTstep=floor(timeStep*numPar/Dur)                      !Numeber of particles released per time step (used floor as nParTstep cannot be fraction).
           nTimeSteps=floor(Dur/timeStep)                      
           !setting initil values
           firstPar=1                                                !First particle of the particles who has same initial start time.
           lastPar=nParTstep                                         !Last particle of the particles who has same initial start time.
           time=0                                                    !As we start from 0 seconds
           
           do p=1,nTimeSteps             
              strTimep(firstPar:lastPar)=time                        !StrTimep tells the starting time for each particle
              firstPar=lastPar+1
              lastPar=lastPar+nParTstep
              time=time+timeStep              
           enddo
        endif
         
  
  
        strPar=1
        endPar=nParTstep 
        !Moving Particles.
        !Notes: strTime is of particles, i.e. strTime(particles) NOT strTime(time)



lpTS:   do t=1,nTimeSteps 
           
           timeForAvging=timeForAvging+timeStep 
           
           do pIni=strPar,endPar
              
              xPosCell=int(xPos(pIni)/dx) + 1                        ! getting nearest cell location CHECK THIS**
              yPosCell=int(yPos(pIni)/dy) + 1
              zPosCell=int((zPos(pIni) + dz)/dz) + 1
              
              call random_number(x)
              uPrimeCurr(pIni)= sigU(xPosCell,yPosCell,zPosCell)*x
              call random_number(x)
              vPrimeCurr(pIni)= sigV(xPosCell,yPosCell,zPosCell)*x
              call random_number(x)
              wPrimeCurr(pIni)= sigW(xPosCell,yPosCell,zPosCell)*x 
           enddo
           
           ! Here we will  move the previously initialized particles.
           do pMov=1,endPar
              
              if(any(pOutDomain == pMov)) then
                 cycle
              else                 
                 call MoveParticles
              endif
           enddo

           do pConc=1,endPar
              if(any(pOutDomain == pConc)) then
                 cycle
              else
              
                 xConcBoxInd=int( (xPos(pConc)-XConcBoxLowLmt)/(xConcBoxSize) ) + 1
                 xConcBoxInd=sign(float(xConcBoxInd),(xPos(pConc)-XConcBoxLowLmt))   !To check if position of the particle is before the collecting box start, xConcBoxInd will be -ve for that case and will not enter the if condition below
                 yConcBoxInd=int( (yPos(pConc)-YConcBoxLowLmt)/(yConcBoxSize) ) + 1
                 yConcBoxInd=sign(float(yConcBoxInd),(yPos(pConc)-YConcBoxLowLmt))
                 zConcBoxInd=int( (zPos(pConc)-ZConcBoxLowLmt)/(zConcBoxSize) ) + 1
                 zConcBoxInd=sign(float(zConcBoxInd),(zPos(pConc)-ZConcBoxLowLmt))                 

                 if(xConcBoxInd>=1 .and. xConcBoxInd<=numXConcBox .and. yConcBoxInd>=1 .and. yConcBoxInd<=numYConcBox .and. zConcBoxInd>=1 .and. zConcBoxInd<=numZConcBox) then
                    Conc(xConcBoxInd,yConcBoxInd,zConcBoxInd)=Conc(xConcBoxInd,yConcBoxInd,zConcBoxInd) + ConcConst 
                 endif  
              endif
           enddo

           if(timeForAvging >= ConAvgTime) then
             call WriteConc
             Conc=0.0
             timeForAvging=timeForAvging - ConAvgTime
           endif

           !the below two line MUST go after the above do loop as "endPar" changes in the lines below.
           strpar=endPar+1                                              ! This is true ONLY for continous release
           endPar=endPar+nParTstep                                      ! This is true ONLY for continous release

              
        enddo lpTS
        end
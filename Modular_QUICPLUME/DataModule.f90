       module datamodule
          Implicit None
          
          integer i,j,k,nx,ny,nz
          integer windTimeSteps,windTimeIncr,windTimeStamp,roofFlag,upwind,streetCanyon
          integer, dimension(:,:,:),allocatable:: icellflag

          real dx,dy,dz,vnKarman
          real tm(2)                                                          ! for time elapsed          
          real znautDiss,wallRough
          real Tau11,Tau12,Tau13,Tau21,Tau22,Tau23,Tau31,Tau32,Tau33,Lam11,Lam12,Lam13,Lam21,Lam22,Lam23,Lam31,Lam32,Lam33
          real DLam11Dx,DLam12Dx,DLam13Dx,DLam21Dx,DLam22Dx,DLam23Dx,DLam31Dx,DLam32Dx,DLam33Dx
          real DLam11Dy,DLam12Dy,DLam13Dy,DLam21Dy,DLam22Dy,DLam23Dy,DLam31Dy,DLam32Dy,DLam33Dy
          real DLam11Dz,DLam12Dz,DLam13Dz,DLam21Dz,DLam22Dz,DLam23Dz,DLam31Dz,DLam32Dz,DLam33Dz          
          real DTau11Dx,DTau12Dx,DTau13Dx,DTau21Dx,DTau22Dx,DTau23Dx,DTau31Dx,DTau32Dx,DTau33Dx
          real DTau11Dy,DTau12Dy,DTau13Dy,DTau21Dy,DTau22Dy,DTau23Dy,DTau31Dy,DTau32Dy,DTau33Dy
          real DTau11Dz,DTau12Dz,DTau13Dz,DTau21Dz,DTau22Dz,DTau23Dz,DTau31Dz,DTau32Dz,DTau33Dz

          real, dimension(:,:,:),allocatable:: DuDx,DuDy,DuDz,DvDx,DvDy,DvDz,DwDx,DwDy,DwDz          
          real, dimension(:,:,:),allocatable:: u,v,w,disBotSurf,disTopSurf,disNegXSurf,disPosXSurf,disNegYSurf,disPosYSurf
          real, dimension(:,:,:),allocatable:: sigU,sigV,sigW

          
          
       end module



       module disperseModule
           
          integer pMov,xPosCell,yPosCell,zPosCell,t,numPar,RelType,SrcType 
          integer numXConcBox,numYConcBox,numZConcBox
          integer, dimension(:),allocatable:: pOutDomain

          real Dur,ConAvgTime,strTime,TotRel,xCordSrc,yCordSrc,zCordSrc,radSrc,timeStep
          real dissXCenter,dissYCenter,dissZCenter,uMeanInterp,vMeanInterp,wMeanInterp,xDomainLen,yDomainLen,zDomainLen
          real XConcBoxLowLmt,XConcBoxUpLmt,YConcBoxLowLmt,YConcBoxUpLmt,ZConcBoxLowLmt,ZConcBoxUpLmt
          real, dimension(:),allocatable:: xPos,yPos,zPos,uPrimeCurr,vPrimeCurr,wPrimeCurr,uPrimePrev,vPrimePrev,wPrimePrev
          real, dimension(:),allocatable:: strTimep,xConcBoxCen,yConcBoxCen,zConcBoxCen
          real, dimension(:,:,:),allocatable:: Conc
          

       end module
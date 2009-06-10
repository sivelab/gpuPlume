        subroutine InitVars
        use DataModule
        Implicit None
        
        vnKarman=0.4
        allocate(disBotSurf(nx,ny,nz),disNegXSurf(nx,ny,nz),disPosXSurf(nx,ny,nz),disNegYSurf(nx,ny,nz),disPosYSurf(nx,ny,nz))
        allocate(disTopSurf(nx,ny,nz),sigU(nx,ny,nz),sigV(nx,ny,nz),sigW(nx,ny,nz))
        allocate(DuDx(nx,ny,nz),DuDy(nx,ny,nz),DuDz(nx,ny,nz),DvDx(nx,ny,nz),DvDy(nx,ny,nz),DvDz(nx,ny,nz),DwDx(nx,ny,nz),DwDy(nx,ny,nz),DwDz(nx,ny,nz))
        
        disTopSurf=0
        disBotSurf=0
        disNegXSurf=0
        disPosXSurf=0
        disNegYSurf=0
        disPosYSurf=0
        
        
        DuDx=0
        DuDy=0
        DuDz=0
        DvDx=0
        DvDy=0
        DvDz=0
        DwDx=0
        DwDy=0
        DwDz=0    
                


        DLam11Dx=0
        DLam12Dx=0
        DLam13Dx=0
        DLam21Dx=0
        DLam22Dx=0
        DLam23Dx=0
        DLam31Dx=0
        DLam32Dx=0
        DLam33Dx=0
        DLam11Dy=0
        DLam12Dy=0
        DLam13Dy=0
        DLam21Dy=0
        DLam22Dy=0
        DLam23Dy=0
        DLam31Dy=0
        DLam32Dy=0
        DLam33Dy=0
        DLam11Dz=0
        DLam12Dz=0
        DLam13Dz=0
        DLam21Dz=0
        DLam22Dz=0
        DLam23Dz=0
        DLam31Dz=0
        DLam32Dz=0
        DLam33Dz=0
          
        DTau11Dx=0
        DTau12Dx=0
        DTau13Dx=0
        DTau21Dx=0
        DTau22Dx=0
        DTau23Dx=0
        DTau31Dx=0
        DTau32Dx=0
        DTau33Dx=0
        DTau11Dy=0
        DTau12Dy=0
        DTau13Dy=0
        DTau21Dy=0
        DTau22Dy=0
        DTau23Dy=0
        DTau31Dy=0
        DTau32Dy=0
        DTau33Dy=0
        DTau11Dz=0
        DTau12Dz=0
        DTau13Dz=0
        DTau21Dz=0
        DTau22Dz=0
        DTau23Dz=0
        DTau31Dz=0
        DTau32Dz=0
        DTau33Dz=0
        
    

        end

        subroutine WriteFiles
        use DataModule
        Implicit None 
        !---------------------------
        !Write all FIles
        !called by Main.f90
        !---------------------------

    
        
        open(50,file='test.dat',status='unknown')
        open(52,file='OutputFiles/Turbulence.dat',status='unknown')
    
        write(50,'(4i6,2x,e10.3)')(((i,j,k,icellflag(i,j,k),i=1,nx),j=1,ny),k=1,nz)  
        write(52,'(3i6,2x,3e123.3)')(((i,j,k,sigU(i,j,k),sigV(i,j,k),sigW(i,j,k),i=1,nx),j=1,ny),k=1,nz)          
    
        return
        end


        subroutine WriteConc
        use DataModule
        use disperseModule
        Implicit None 
        



        write(51,102)float(t)*timeStep
        write(51,103),numXConcBox,numYConcBox,numZConcBox
        write(51,104)           
        write(51,'(4e10.3)')((((xConcBoxCen(i),yConcBoxCen(j),yConcBoxCen(k),Conc(i,j,k)), i=1,numXConcBox),j=1,numYConcBox),k=1,numZConcBox)
        
102     format('ZONE T="',f8.1,' seconds"')   
103     format('I=',i6,', J=',i6,',  K=',i6,', F=POINT') 
104     format('DT=(SINGLE SINGLE SINGLE SINGLE)')    
    
        end
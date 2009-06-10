        subroutine ReadFiles
        use DataModule
        Implicit None

        real dum
        
        open(1,file='Modular_QUICPLUME/InputFiles/QU_Velocity.dat',status='old')  
        open(2,file='Modular_QUICPLUME/InputFiles/QU_Celltype.dat',status='old')
        open(3,file='Modular_QUICPLUME/InputFiles/QU_buildings.inp',status='old')
        open(4,file='Modular_QUICPLUME/InputFiles/QU_simparams.inp',status='old')

        !For writing

        !Reading Simparams
        read(4,*)nx
        read(4,*)ny
        read(4,*)nz
        read(4,*)dx
        read(4,*)dy
        read(4,*)dz
        read(4,*)
        read(4,*)
        read(4,*)
        read(4,*)roofFlag
        read(4,*)upwind
        read(4,*)streetCanyon
            
        !Reading QU_Velocity
        read(1,*)           !Reads the char string, not stored in a var as it was not imp
        read(1,*),windTimeSteps
        read(1,*)           !Reads the char string, not stored in a var as it was not imp
        read(1,*),windTimeIncr
        read(1,*),windTimeStamp

        allocate(u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz),icellflag(nx,ny,nz))

        do k=1,nz
           do j=1,ny
                do i=1,nx
                   read(1,*)dum,dum,dum,u(i,j,k),v(i,j,k),w(i,j,k)  !dum reads the value of x,y and z location of QU cell
                   read(2,*)dum,dum,dum,icellflag(i,j,k)            !dum reads the value of x,y and z location of QU cell
                   icellflag(i,j,k)=min(1,icellflag(i,j,k))         ! 0 for surface and 1 for free atomsphere
              enddo
           enddo
        enddo

        !input from QU_buildings.inp
        read(3,*)x_subdomain_start	! Subdomain coordinates
        read(3,*)y_subdomain_start	!
        read(3,*)x_subdomain_end	!
        read(3,*)y_subdomain_end	!
        read(3,*)zo	! MAN 8-19-2005 Updated input output file structures

        ! MAN 7/27/2005 var dz subdomain changed into meters
        x_subdomain_start = x_subdomain_start*dx
        y_subdomain_start = y_subdomain_start*dy
        x_subdomain_end = x_subdomain_end*dx
        y_subdomain_end = y_subdomain_end*dy
        !end MAN 7/27/2005

        read(3,*)inumbuild	!number of buildings
        read(3,*)

        allocate(xfo(inumbuild),yfo(inumbuild),zfo(inumbuild))
        allocate(gamma(inumbuild))
        allocate(aa(inumbuild),bb(inumbuild))
        allocate(Ht(inumbuild),Wti(inumbuild),Lti(inumbuild))
        allocate(bldnum(inumbuild),bldtype(inumbuild),group_id(inumbuild))
        allocate(atten(inumbuild))

        ! Read in the building number, type, height, width, length ,xfo,yfo,zfo,gamma and atten
        ! atten - attenuation coefficient for vegetation
        !erp 1/31/2003
        ! note that for now if the building is cylindrical, enter Lti = 0.

        do i=1,inumbuild
           read(3,*)bldnum(i),group_id(i),bldtype(i),Ht(i),Wti(i),Lti(i),xfo(i),yfo(i),zfo(i),gamma(i),atten(i)

           Ht(i)=Ht(i)*dz	!convert back to real world unit, TZ 10/29/04
           Wti(i)=Wti(i)*dy	!convert back to real world unit, TZ 10/29/04
           Lti(i)=Lti(i)*dx	!convert back to real world unit, TZ 10/29/04
           xfo(i)=xfo(i)*dx	!convert back to real world unit, TZ 10/29/04
           yfo(i)=yfo(i)*dy	!convert back to real world unit, TZ 10/29/04
           zfo(i)=zfo(i)*dz	!convert back to real world unit, TZ 10/29/04

           Ht(i)=Ht(i)+zfo(i)			!erp 7/25/03
           !erp	    zfo(i)=zfo(i)+2				!erp 7/23/03
           if(bldtype(i).eq.2)then		!if the building is a cylinder/ellipse
              bb(i)=Wti(i)/2.			!set minor axis to input Width
              aa(i)=Lti(i)/2.			!set major axis to input Lenth
           endif

           !erp 10/22/04 Pentagon addition
           if(bldtype(i).eq.3)then		!if the building is a Pentagon
              bb(i)=Wti(i)/2.				!Radius Pentagon is inscribed in
              xfo(i)=xfo(i)-bb(i)
           endif
        enddo
        
        close(1)
        close(2)
        close(3)
        close(4)        !closed Celltype.dat file as we are not going to use it again
        return
        end











        Subroutine ReadFilesDispersion
        use DataModule
        use disperseModule
        Implicit None
        
        !For Reading
        open(10,file='Modular_QUICPLUME/InputFiles/QP_Params.inp', status='old')
        open(11,file='Modular_QUICPLUME/InputFiles/QP_Source.inp', status='old')
        !For writing
        open(51,file='Modular_QUICPLUME/OutputFiles/QP_Confield.dat', status='unknown')
        
        
        write(51,100)
        write(51,101)
100     format('TITLE     = "Concentrations in 3D (grams per cubic meter)"')
101     format('VARIABLES = "X" "Y" "Z" "C"') 





        read(10,*);read(10,*);read(10,*);read(10,*);read(10,*);read(10,*);read(10,*)

        read(10,*)znautDiss
        read(10,*)wallRough

        read(10,*);read(10,*);read(10,*);read(10,*)

        read(10,*)numPar
        read(10,*)timeStep
        read(10,*)Dur
        read(10,*)ConAvgTime
        read(10,*)strTime
        
        read(10,*);read(10,*)

        read(10,*)numXConcBox
        read(10,*)numYConcBox
        read(10,*)numZConcBox
        read(10,*)XConcBoxLowLmt
        read(10,*)XConcBoxUpLmt
        read(10,*)YConcBoxLowLmt
        read(10,*)YConcBoxUpLmt
        read(10,*)ZConcBoxLowLmt
        read(10,*)ZConcBoxUpLmt


        read(11,*);read(11,*)

        read(11,*)TotRel
        read(11,*)RelType
        
        read(11,*);read(11,*)

        read(11,*)SrcType
        read(11,*)xCordSrc
        read(11,*)yCordSrc
        read(11,*)zCordSrc
        read(11,*)radSrc

        !Whatif if source is in the surface cell or building cell
       
        return
        end



           


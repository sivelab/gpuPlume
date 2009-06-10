        program main
        use datamodule
        
        
        call ReadFiles
        call InitVars
        call Turbulence
        call ReadFilesDispersion
        call Dispersion
        call WriteFiles
        call DeallocateVars
        print*,'Waheguru'
        
        print *,etime(tm)
        
        print *,tm(1)
        print *,tm(2)
        


        stop
        end
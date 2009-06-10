      
      subroutine RanGauss(x)
         real xGauss
         call random_number(xGauss)
         u=xGauss
         if (u.le..8638) then
            call random_number(xGauss)
            v=2.0*xGauss-1.0
            call random_number(xGauss)
            w=2.0*xGauss-1.0
            x=2.3153508*u-1.0+v+w  
            return       
         endif
         if (u.le..9745) then
            call random_number(xGauss)
            v=xGauss
            x=1.5*(v-1.0+9.0334237*(u-.8638))   
            return
         endif
         if (u.gt..9973002) then 
            v = 4.5
            X = 1. 
            do while((x*v*v).gt.4.5)
               call random_number(xGauss)
               v=xGauss
               call random_number(xGauss)
               w=xGauss
               w=max(1.e-07,w)
               x=4.5-log(w)
            enddo
            x=sign(sqrt(2.0*x),(u-.9986501))
         else
            u = 50.
            v = 0.
            sum = 49.0024445
            w = 0.
            do while(u.gt.(49.0024445*exp(-v*v/2.0)-sum-w))
               call random_number(xGauss)
               x=6.0*xGauss-3.0
               call random_number(xGauss)
               u=xGauss
               v=abs(x)
               w=6.6313339*(3.0-v)**2.0
               sum=0.
               if (v.lt.1.5) sum=6.0432809*(1.5-v)
               if (v.lt.1.0)  sum=sum+13.2626678*(3.0-v*v)-w  
            enddo
         endif
         return
      end
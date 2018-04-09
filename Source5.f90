!  Console1.f90 
!
!  FUNCTIONS:
!  Console1 - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Console1
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program Console1

    implicit none

        ! 变量声明
        external aW1, aE1, aW2, aE2, aW3, aE3
        REAL(8) :: aE, aW, aP, De, Fe, Dw, Fw, dif, den, u, L, dx, fai0, faiL, aE1, aW1, aW2, aE2, aW3, aE3
        INTEGER :: I,J,K,n
        !已知参数
        dif=1.0e-1
        den=1.0e0
        L=1.0e0
        n=5
        fai0=1
        faiL=0
        dx=L/n
        !输入速度
        print *,'请输入速度u='
        READ(*,*)u
        !计算D和F
        De=dif/dx
        Dw=dif/dx
        Fe=den*u
        Fw=den*u
        !计算aE, aW, aP
        DO K=1,3
            IF(K==1) THEN
                aE=aE1(De,Fe)
                aW=aW1(Dw,Fw) 
            ELSEIF(K==2) THEN
                aE=aE2(De,Fe)
                aW=aW2(Dw,Fw) 
            ELSE
                aE=aE3(De,Fe)
                aW=aW3(Dw,Fw) 
            ENDIF
            aP=aE+aW+Fe-Fw
        !得到系数矩阵，求解方程组
        call solve(aE,aW,aP,n,fai0,faiL,K)
        ENDDO
        !计算精确解
        call exact(faiL,fai0,Fe,dx,L,n,dif)
    end program Console1
    
    
    
    
    !使用不同方法求aW，aE
    function aE1(De,Fe)
        implicit none
        real*8::aE1,De,Fe
        aE1=De-Fe/2
    end function
    
    function aW1(Dw,Fw)
        implicit none
        real*8::aW1,Dw,Fw
        aW1=Dw+Fw/2
    end function
    
    function aE2(De,Fe)
        implicit none
        real*8::aE2,De,Fe
        IF(Fe>0) THEN
            aE2=De
        ELSE
            aE2=De-Fe
        ENDIF
    end function  
    
    function aW2(Dw,Fw)
        implicit none
        real*8::aW2,Dw,Fw
        IF(Fw>0) THEN
            aW2=Dw+Fw
        ELSE
            aW2=Dw
        ENDIF
    end function
    
    function aE3(De,Fe)
        implicit none
        real*8::aE3,De,Fe
        IF((De-Fe/2)>0) THEN
            IF((De-Fe/2)>(-Fe)) THEN
                aE3=(De-Fe/2)
            ELse
                aE3=-Fe
            ENDIF 
        ELSE
            IF((-Fe)>0) THEN
                aE3=-Fe
            ELSE
                aE3=0
            ENDIF 
        ENDIF
    end function
    
     function aW3(Dw,Fw)
        implicit none
        real*8::aW3,Dw,Fw
        IF(Fw>0) THEN
            IF(Fw>(Dw+Fw/2)) THEN
                aW3=Fw
            ELSE
                aW3=Dw+Fw/2
            ENDIF
        ELSE
            IF((Dw+Fw/2)>0) THEN
                aW3=Dw+Fw/2
            ELSE
                aW3=0
            ENDIF
        ENDIF
    end function
     
    subroutine solve(aE,aW,aP,n,fai0,faiL,K)
        implicit none
        real*8::aE,aW,aP,fai0,faiL,A(6,6), b(6), c(5), d(6), x(6)
        INTEGER :: I,J,K,n
        !系数矩阵初始化
        !首先全部初始化为0
        DO I=1,(n+1)
            DO J=1,(n+1)
            A(I,J)=0
            ENDDO
            b(I)=0
        ENDDO
        !初始化三对角矩阵
        A(1,1)=1
        A(1,2)=0
        b(1)=fai0
    
        DO I=2,n
            A(I,I-1)=-aW
            A(I,I)=aP
            A(I,I+1)=-aE
            b(I)=0
        ENDDO
    
        
        A(n+1,n)=0
        A(n+1,n+1)=1
        b(n+1)=faiL 
    
        !TDMA求解方程
        c(1)=A(1,2)/A(1,1)
        d(1)=b(1)/A(1,1)
        DO I=2,(n)
            c(I)=A(I,I+1)/(A(I,I)-c(I-1)*A(I,I-1))
            d(I)=(b(I)-d(I-1)*A(I,I-1))/(A(I,I)-c(I-1)*A(I,I-1))
        ENDDO
        d(n+1)=(b(n+1)-d(n)*A(n+1,n))/(A(n+1,n+1)-c(n)*A(n+1,n))
        x(n+1)=d(n+1)
        !输出结果
        DO I=1,(n)
            x(n+1-I)=d(n+1-I)-c(n+1-I)*x(n+1-I+1)
        ENDDO
        IF(K==1) THEN
            print *, '中心查分格式fai=',x 
        ELSEIF(K==2) THEN
            print *, '迎风格式fai=',x 
        ELSE
            print *, '混合格式fai=',x  
        ENDIF
        pause
    end subroutine solve

    subroutine exact(faiL,fai0,Fe,dx,L,n,dif)
        implicit none
        real*8::faiL,fai0,Fe,dx,L,dif,x(6)
        INTEGER ::I, n
        DO I=1,(n+1)
            x(I)=(faiL-fai0)*((EXP(Fe*dx*(I-1)/(dif))-1)/(EXP(Fe*L/dif)-1))+fai0
        ENDDO
        print *, '精确解fai=',x 
        pause
    end subroutine exact
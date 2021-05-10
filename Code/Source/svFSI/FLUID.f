!
! Copyright (c) Stanford University, The Regents of the University of
!               California, and others.
!
! All Rights Reserved.
!
! See Copyright-SimVascular.txt for additional details.
!
! Permission is hereby granted, free of charge, to any person obtaining
! a copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject
! to the following conditions:
!
! The above copyright notice and this permission notice shall be included
! in all copies or substantial portions of the Software.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!--------------------------------------------------------------------
!
!     This is for solving fluid transport equation solving Navier-Stokes
!     equations. Dirichlet boundary conditions are either treated
!     strongly or weakly.
!
!--------------------------------------------------------------------

      SUBROUTINE CONSTRUCT_FLUID(lM, Ag, Yg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo)

      INTEGER(KIND=IKIND) a, e, g, Ac, eNoN, cPhys
      REAL(KIND=RKIND) w, Jac, ksix(nsd,nsd)

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), al(:,:), yl(:,:),
     2   bfl(:,:), N(:), Nx(:,:), lR(:,:), lK(:,:,:)

      eNoN = lM%eNoN

!     FLUID: dof = nsd+1
      ALLOCATE(ptr(eNoN), xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN),
     2   bfl(nsd,eNoN), N(eNoN), Nx(nsd,eNoN), lR(dof,eNoN),
     3   lK(dof*dof,eNoN,eNoN))

!     Loop over all elements of mesh
      DO e=1, lM%nEl
!        Update domain and proceed if domain phys and eqn phys match
         cDmn  = DOMAIN(lM, cEq, e)
         cPhys = eq(cEq)%dmn(cDmn)%phys
         IF (cPhys .NE. phys_fluid) CYCLE

!        Update shape functions for NURBS
         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

!        Create local copies
         DO a=1, eNoN
            Ac = lM%IEN(a,e)
            ptr(a)   = Ac
            xl(:,a)  = x(:,Ac)
            al(:,a)  = Ag(:,Ac)
            yl(:,a)  = Yg(:,Ac)
            bfl(:,a) = Bf(:,Ac) ! JP 2021_05_04: "Bf" stands for "Body force" from MOD.f
         END DO

!        Gauss integration
         lR = 0._RKIND
         lK = 0._RKIND
         DO g=1, lM%nG
            IF (g.EQ.1 .OR. .NOT.lM%lShpF) THEN
               CALL GNN(eNoN, nsd, lM%Nx(:,:,g), xl, Nx, Jac, ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
            END IF
            w = lM%w(g) * Jac
            N = lM%N(:,g)

            IF (nsd .EQ. 3) THEN
               CALL FLUID3D(eNoN, w, N, Nx, al, yl, bfl, ksix, lR, lK)

            ELSE IF (nsd .EQ. 2) THEN
               CALL FLUID2D(eNoN, w, N, Nx, al, yl, bfl, ksix, lR, lK)

            END IF
         END DO ! g: loop

!        Assembly
#ifdef WITH_TRILINOS
         IF (eq(cEq)%assmTLS) THEN
            CALL TRILINOS_DOASSEM(eNoN, ptr, lK, lR)
         ELSE
#endif
            CALL DOASSEM(eNoN, ptr, lK, lR)
#ifdef WITH_TRILINOS
         END IF
#endif
      END DO ! e: loop

      DEALLOCATE(ptr, xl, al, yl, bfl, N, Nx, lR, lK)

      RETURN
      END SUBROUTINE CONSTRUCT_FLUID
!####################################################################
      SUBROUTINE FLUID3D(eNoN, w, N, Nx, al, yl, bfl, Kxi, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(3,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN), bfl(3,eNoN), Kxi(3,3)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      INTEGER(KIND=IKIND) a, b
      REAL(KIND=RKIND) ctM, ctC, tauM, tauC, tauB, kT, kS, kU, mu, rho,
     2   divU, amd, wl, wr, p, pa, u(3), ud(3), px(3), f(3), up(3),
     3   ua(3), ux(3,3), es(3,3), rV(3), rM(3,3), uNx(eNoN), upNx(eNoN),
     4   uaNx(eNoN), NxNx, gam, mu_s, mu_x, es_x(3,eNoN), T1, T2, T3

      ctM  = 1._RKIND
      ctC  = 36._RKIND

      rho  = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      f(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      f(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      f(3) = eq(cEq)%dmn(cDmn)%prop(f_z)

      T1   = eq(cEq)%af * eq(cEq)%gam * dt
      amd  = eq(cEq)%am/T1
      wl   = w*T1
      wr   = w*rho

!     Indices are not selected based on the equation only
!     because fluid equation always come first
      p  = 0._RKIND
      u  = 0._RKIND
      ud = -f
      px = 0._RKIND
      ux = 0._RKIND
      DO a=1, eNoN
         p  = p + N(a)*yl(4,a)

         ud(1) = ud(1) + N(a)*(al(1,a)-bfl(1,a))
         ud(2) = ud(2) + N(a)*(al(2,a)-bfl(2,a))
         ud(3) = ud(3) + N(a)*(al(3,a)-bfl(3,a))

         px(1) = px(1) + Nx(1,a)*yl(4,a)
         px(2) = px(2) + Nx(2,a)*yl(4,a)
         px(3) = px(3) + Nx(3,a)*yl(4,a)

         u(1) = u(1) + N(a)*yl(1,a)
         u(2) = u(2) + N(a)*yl(2,a)
         u(3) = u(3) + N(a)*yl(3,a)

         ux(1,1) = ux(1,1) + Nx(1,a)*yl(1,a)
         ux(2,1) = ux(2,1) + Nx(2,a)*yl(1,a)
         ux(3,1) = ux(3,1) + Nx(3,a)*yl(1,a)
         ux(1,2) = ux(1,2) + Nx(1,a)*yl(2,a)
         ux(2,2) = ux(2,2) + Nx(2,a)*yl(2,a)
         ux(3,2) = ux(3,2) + Nx(3,a)*yl(2,a)
         ux(1,3) = ux(1,3) + Nx(1,a)*yl(3,a)
         ux(2,3) = ux(2,3) + Nx(2,a)*yl(3,a)
         ux(3,3) = ux(3,3) + Nx(3,a)*yl(3,a)
      END DO
      divU = ux(1,1) + ux(2,2) + ux(3,3)

      IF (mvMsh) THEN
         DO a=1, eNoN
            u(1) = u(1) - N(a)*yl(5,a)
            u(2) = u(2) - N(a)*yl(6,a)
            u(3) = u(3) - N(a)*yl(7,a)
         END DO
      END IF

!     Strain rate tensor 2*e_ij := (u_ij + u_ji)
      es(1,1) = ux(1,1) + ux(1,1)
      es(2,1) = ux(2,1) + ux(1,2)
      es(3,1) = ux(3,1) + ux(1,3)
      es(1,2) = es(2,1)
      es(2,2) = ux(2,2) + ux(2,2)
      es(3,2) = ux(3,2) + ux(2,3)
      es(1,3) = es(3,1)
      es(2,3) = es(3,2)
      es(3,3) = ux(3,3) + ux(3,3)

      DO a=1, eNoN
        es_x(1,a) = es(1,1)*Nx(1,a) + es(2,1)*Nx(2,a) + es(3,1)*Nx(3,a)
        es_x(2,a) = es(1,2)*Nx(1,a) + es(2,2)*Nx(2,a) + es(3,2)*Nx(3,a)
        es_x(3,a) = es(1,3)*Nx(1,a) + es(2,3)*Nx(2,a) + es(3,3)*Nx(3,a)
      END DO

!     Shear-rate := (2*e_ij*e_ij)^.5
      gam = es(1,1)*es(1,1) + es(2,1)*es(2,1) + es(3,1)*es(3,1)
     2    + es(1,2)*es(1,2) + es(2,2)*es(2,2) + es(3,2)*es(3,2)
     3    + es(1,3)*es(1,3) + es(2,3)*es(2,3) + es(3,3)*es(3,3)
      gam = SQRT(0.5_RKIND*gam)

!     Compute viscosity based on shear-rate and chosen viscosity model
      CALL GETVISCOSITY(eq(cEq)%dmn(cDmn), gam, mu, mu_s, mu_x)
      IF (ISZERO(gam)) THEN
         mu_x = 0._RKIND
      ELSE
         mu_x = mu_x/gam
      END IF

      kT = 4._RKIND*(ctM/dt)**2_RKIND

      kU = u(1)*u(1)*Kxi(1,1) + u(2)*u(1)*Kxi(2,1) + u(3)*u(1)*Kxi(3,1)
     2   + u(1)*u(2)*Kxi(1,2) + u(2)*u(2)*Kxi(2,2) + u(3)*u(2)*Kxi(3,2)
     3   + u(1)*u(3)*Kxi(1,3) + u(2)*u(3)*Kxi(2,3) + u(3)*u(3)*Kxi(3,3)

      kS = Kxi(1,1)*Kxi(1,1) + Kxi(2,1)*Kxi(2,1) + Kxi(3,1)*Kxi(3,1)
     2   + Kxi(1,2)*Kxi(1,2) + Kxi(2,2)*Kxi(2,2) + Kxi(3,2)*Kxi(3,2)
     3   + Kxi(1,3)*Kxi(1,3) + Kxi(2,3)*Kxi(2,3) + Kxi(3,3)*Kxi(3,3)
      kS = ctC * kS * (mu/rho)**2._RKIND

      tauM = 1._RKIND / (rho * SQRT( kT + kU + kS ))
      tauC = 1._RKIND / (tauM * (Kxi(1,1) + Kxi(2,2) + Kxi(3,3)))

      rV(1) = ud(1) + u(1)*ux(1,1) + u(2)*ux(2,1) + u(3)*ux(3,1)
      rV(2) = ud(2) + u(1)*ux(1,2) + u(2)*ux(2,2) + u(3)*ux(3,2)
      rV(3) = ud(3) + u(1)*ux(1,3) + u(2)*ux(2,3) + u(3)*ux(3,3)

      up(1) = -tauM*(rho*rV(1) + px(1))
      up(2) = -tauM*(rho*rV(2) + px(2))
      up(3) = -tauM*(rho*rV(3) + px(3))

      tauB = up(1)*up(1)*Kxi(1,1) + up(2)*up(1)*Kxi(2,1)
     2     + up(3)*up(1)*Kxi(3,1) + up(1)*up(2)*Kxi(1,2)
     3     + up(2)*up(2)*Kxi(2,2) + up(3)*up(2)*Kxi(3,2)
     4     + up(1)*up(3)*Kxi(1,3) + up(2)*up(3)*Kxi(2,3)
     5     + up(3)*up(3)*Kxi(3,3)
      IF (ISZERO(tauB)) tauB = eps
      tauB = rho/SQRT(tauB)

      ua(1) = u(1) + up(1)
      ua(2) = u(2) + up(2)
      ua(3) = u(3) + up(3)
      pa    = p - tauC*divU

      rV(1) = tauB*(up(1)*ux(1,1) + up(2)*ux(2,1) + up(3)*ux(3,1))
      rV(2) = tauB*(up(1)*ux(1,2) + up(2)*ux(2,2) + up(3)*ux(3,2))
      rV(3) = tauB*(up(1)*ux(1,3) + up(2)*ux(2,3) + up(3)*ux(3,3))

      rM(1,1) = mu*es(1,1) - rho*up(1)*ua(1) + rV(1)*up(1) - pa
      rM(2,1) = mu*es(2,1) - rho*up(1)*ua(2) + rV(1)*up(2)
      rM(3,1) = mu*es(3,1) - rho*up(1)*ua(3) + rV(1)*up(3)

      rM(1,2) = mu*es(1,2) - rho*up(2)*ua(1) + rV(2)*up(1)
      rM(2,2) = mu*es(2,2) - rho*up(2)*ua(2) + rV(2)*up(2) - pa
      rM(3,2) = mu*es(3,2) - rho*up(2)*ua(3) + rV(2)*up(3)

      rM(1,3) = mu*es(1,3) - rho*up(3)*ua(1) + rV(3)*up(1)
      rM(2,3) = mu*es(2,3) - rho*up(3)*ua(2) + rV(3)*up(2)
      rM(3,3) = mu*es(3,3) - rho*up(3)*ua(3) + rV(3)*up(3) - pa

      rV(1) = ud(1) + ua(1)*ux(1,1) + ua(2)*ux(2,1) + ua(3)*ux(3,1)
      rV(2) = ud(2) + ua(1)*ux(1,2) + ua(2)*ux(2,2) + ua(3)*ux(3,2)
      rV(3) = ud(3) + ua(1)*ux(1,3) + ua(2)*ux(2,3) + ua(3)*ux(3,3)

      DO a=1, eNoN
         uNx(a)  = u(1)*Nx(1,a)  + u(2)*Nx(2,a)  + u(3)*Nx(3,a)
         upNx(a) = up(1)*Nx(1,a) + up(2)*Nx(2,a) + up(3)*Nx(3,a)
         uaNx(a) = uNx(a) + upNx(a)

         lR(1,a) = lR(1,a) + wr*N(a)*rV(1) + w*(Nx(1,a)*rM(1,1)
     2      + Nx(2,a)*rM(2,1) + Nx(3,a)*rM(3,1))

         lR(2,a) = lR(2,a) + wr*N(a)*rV(2) + w*(Nx(1,a)*rM(1,2)
     2      + Nx(2,a)*rM(2,2) + Nx(3,a)*rM(3,2))

         lR(3,a) = lR(3,a) + wr*N(a)*rV(3) + w*(Nx(1,a)*rM(1,3)
     2      + Nx(2,a)*rM(2,3) + Nx(3,a)*rM(3,3))

         lR(4,a) = lR(4,a) + w*(N(a)*divU - upNx(a))
      END DO

      DO a=1, eNoN
         DO b=1, eNoN
            rM(1,1) = Nx(1,a)*Nx(1,b)
            rM(2,1) = Nx(2,a)*Nx(1,b)
            rM(3,1) = Nx(3,a)*Nx(1,b)
            rM(1,2) = Nx(1,a)*Nx(2,b)
            rM(2,2) = Nx(2,a)*Nx(2,b)
            rM(3,2) = Nx(3,a)*Nx(2,b)
            rM(1,3) = Nx(1,a)*Nx(3,b)
            rM(2,3) = Nx(2,a)*Nx(3,b)
            rM(3,3) = Nx(3,a)*Nx(3,b)

            NxNx = Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b) + Nx(3,a)*Nx(3,b)

            T1 = mu*NxNx + tauB*upNx(a)*upNx(b)
     2         + rho*( N(a)*(amd*N(b) + uaNx(b))
     3         + rho*tauM*uaNx(a)*(uNx(b) + amd*N(b)) )

            T2 = rho*tauM*uaNx(a)

            T3 = rho*tauM*(amd*N(b) + uNx(b))

!           dM/dU
            lK(1,a,b)  = lK(1,a,b)  + wl*((mu + tauC)*rM(1,1) + T1
     2         + mu_x*es_x(1,a)*es_x(1,b))
            lK(2,a,b)  = lK(2,a,b)  + wl*(mu*rM(2,1) + tauC*rM(1,2)
     2         + mu_x*es_x(1,a)*es_x(2,b))
            lK(3,a,b)  = lK(3,a,b)  + wl*(mu*rM(3,1) + tauC*rM(1,3)
     2         + mu_x*es_x(1,a)*es_x(3,b))

            lK(5,a,b)  = lK(5,a,b)  + wl*(mu*rM(1,2) + tauC*rM(2,1)
     2         + mu_x*es_x(2,a)*es_x(1,b))
            lK(6,a,b)  = lK(6,a,b)  + wl*((mu + tauC)*rM(2,2) + T1
     2         + mu_x*es_x(2,a)*es_x(2,b))
            lK(7,a,b)  = lK(7,a,b)  + wl*(mu*rM(3,2) + tauC*rM(2,3)
     2         + mu_x*es_x(2,a)*es_x(3,b))

            lK(9,a,b)  = lK(9,a,b)  + wl*(mu*rM(1,3) + tauC*rM(3,1)
     2         + mu_x*es_x(3,a)*es_x(1,b))
            lK(10,a,b) = lK(10,a,b) + wl*(mu*rM(2,3) + tauC*rM(3,2)
     2         + mu_x*es_x(3,a)*es_x(2,b))
            lK(11,a,b) = lK(11,a,b) + wl*((mu + tauC)*rM(3,3) + T1
     2         + mu_x*es_x(3,a)*es_x(3,b))

!           dM/dP
            lK(4,a,b)  = lK(4,a,b)  - wl*(Nx(1,a)*N(b) - Nx(1,b)*T2)
            lK(8,a,b)  = lK(8,a,b)  - wl*(Nx(2,a)*N(b) - Nx(2,b)*T2)
            lK(12,a,b) = lK(12,a,b) - wl*(Nx(3,a)*N(b) - Nx(3,b)*T2)

!           dC/dU
            lK(13,a,b) = lK(13,a,b) + wl*(N(a)*Nx(1,b) + Nx(1,a)*T3)
            lK(14,a,b) = lK(14,a,b) + wl*(N(a)*Nx(2,b) + Nx(2,a)*T3)
            lK(15,a,b) = lK(15,a,b) + wl*(N(a)*Nx(3,b) + Nx(3,a)*T3)

!           dC/dP
            lK(16,a,b) = lK(16,a,b) + wl*(tauM*NxNx)
         END DO
      END DO

      RETURN
      END SUBROUTINE FLUID3D
!--------------------------------------------------------------------
      SUBROUTINE FLUID2D(eNoN, w, N, Nx, al, yl, bfl, Kxi, lR, lK)
        ! JP 2021_05_04: bfl stands for body force local (the local element body source)
        ! JP 2021_05_04: I think Kxi represents the greek letter, "xi", which Mahdi uses in his thesis (eqn 2.8, pg 15) to compute tau_M, tau_bar, and tau_C parameters. where xi is a matrix of size nsd x nsd, which is exactly what is coded below: "Kxi(2,2)", where here we are in 2D (FLUID2D), so nsd = 2.
        ! JP 2021_05_04: I think eNoN is "Number of nodes (control points) in a single element" from MOD.f (under section, TYPE mshType)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(2,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN), bfl(2,eNoN), Kxi(2,2)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      INTEGER(KIND=IKIND) a, b
      REAL(KIND=RKIND) ctM, ctC, tauM, tauC, tauB, kT, kS, kU, mu, rho,
     2   divU, amd, wl, wr, p, pa, u(2), ud(2), px(2), f(2), up(2),
     3   ua(2), ux(2,2), es(2,2), rV(2), rM(2,2), uNx(eNoN), upNx(eNoN),
     4   uaNx(eNoN), NxNx, gam, mu_s, mu_x, es_x(2,eNoN), T1, T2, T3

      ctM  = 1._RKIND ! JP 2021_05_04: ctM represents c1 in Madhi's thesis, eqn 2.8 pg 15
      ctC  = 36._RKIND ! JP 2021_05_04: ctC represents c2 in Madhi's thesis, eqn 2.8 pg 15

      rho  = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      f(1) = eq(cEq)%dmn(cDmn)%prop(f_x) ! JP 2021_05_04: what is the difference between "f" and bfl? Are they both body forces?
      f(2) = eq(cEq)%dmn(cDmn)%prop(f_y)

      T1   = eq(cEq)%af * eq(cEq)%gam * dt
      amd  = eq(cEq)%am/T1
      wl   = w*T1
      wr   = w*rho

!     Indices are not selected based on the equation only
!     because fluid equation always come first
      ! JP 2021_05_04: this section computes the pressure and the velocity, using the shape functions. E.g. recall that the pressure = sum (over the nodes in the local element) of the pressure degree-of-freedom times the shape function. The pressure and velocity here are the values at the t_(n + alpha_f) time step
      p  = 0._RKIND ! JP 2021_05_04: p = pressure?
      u  = 0._RKIND ! JP 2021_05_04: u = velocity?
      ud = -f
      px = 0._RKIND ! JP 2021_05_04: px = spatial derivative of pressure?
      ux = 0._RKIND ! JP 2021_05_04: ux = spatial derivative of velocity?
      DO a=1, eNoN
         p  = p + N(a)*yl(3,a) ! JP 2021_05_04: note that yl is an array of shape (tDof, eNoN), where tDof = "Total number of degrees of freedom per node" (from MOD.f). So that means that the pressure soltn is stored in the 3rd row of yl and the x and y velocities are stored in the 1st and 2nd rows of yl. and yl is obtained from "yl(:,a)  = Yg(:,Ac)", where Yg contains the value of the all soltns in the model evaluated at the t_(n + alpha_f) time step (as discussed in MAIN.f and PIC.f)

         ud(1) = ud(1) + N(a)*(al(1,a)-bfl(1,a)) ! JP 2021_05_04: al(1,a) is the time derivative of x-velocity
         ud(2) = ud(2) + N(a)*(al(2,a)-bfl(2,a)) ! JP 2021_05_04: al(2,a) is the time derivative of y-velocity

         px(1) = px(1) + Nx(1,a)*yl(3,a) ! JP 2021_05_05: this line computes the pressure gradient
         px(2) = px(2) + Nx(2,a)*yl(3,a)

         u(1) = u(1) + N(a)*yl(1,a)
         u(2) = u(2) + N(a)*yl(2,a)

         ux(1,1) = ux(1,1) + Nx(1,a)*yl(1,a) ! JP 2021_05_04: ux is an array with shape (nsd, nsd) (here nsd = 2, for 2D problem). ux(i, j) = the ith derivative of the jth component of u. So for example, ux(1, 2) represents the x-derivative of the y-velocity.
         ux(2,1) = ux(2,1) + Nx(2,a)*yl(1,a)
         ux(1,2) = ux(1,2) + Nx(1,a)*yl(2,a)
         ux(2,2) = ux(2,2) + Nx(2,a)*yl(2,a)
      END DO
      divU = ux(1,1) + ux(2,2) ! JP 2021_05_04: divU = divergence of velocity, evaluated at the t_(n + alpha_f) time step

      IF (mvMsh) THEN ! JP 2021_04_21: according to MOD.f, mvMsh is a boolean that indicates "Whether mesh is moving"; I think this parameter is used only in FSI cases (see READFILES.f with the line "CASE ('FSI')", where this section sets mvMsh to be true (mvMsh has a value of false by default))
         DO a=1, eNoN
            u(1) = u(1) - N(a)*yl(4,a)
            u(2) = u(2) - N(a)*yl(5,a)
         END DO
      END IF

!     Strain rate tensor 2*e_ij := (u_ij + u_ji) ! JP 2021_05_04: I think this section computes the 2 times the symmetric gradient of velocity, where recall that the symmetric gradient of velocity is = 0.5 * (grad(u) + (grad(u))^T)
      es(1,1) = ux(1,1) + ux(1,1)
      es(2,1) = ux(2,1) + ux(1,2)
      es(1,2) = es(2,1)
      es(2,2) = ux(2,2) + ux(2,2)

      DO a=1, eNoN
        es_x(1,a) = es(1,1)*Nx(1,a) + es(2,1)*Nx(2,a) ! the spatial derivatives of the "Strain rate tensor" or the divergence of the "Strain rate tensor"???
        es_x(2,a) = es(1,2)*Nx(1,a) + es(2,2)*Nx(2,a)
      END DO

!     Shear-rate := (2*e_ij*e_ij)^.5
      gam = es(1,1)*es(1,1) + es(2,1)*es(2,1)
     2    + es(1,2)*es(1,2) + es(2,2)*es(2,2)
      gam = SQRT(0.5_RKIND*gam)

!     Compute viscosity based on shear-rate and chosen viscosity model
      CALL GETVISCOSITY(eq(cEq)%dmn(cDmn), gam, mu, mu_s, mu_x)
      IF (ISZERO(gam)) THEN
         mu_x = 0._RKIND
      ELSE
         mu_x = mu_x/gam
      END IF

      kT = 4._RKIND*(ctM/dt)**2_RKIND

      kU = u(1)*u(1)*Kxi(1,1) + u(2)*u(1)*Kxi(2,1)
     2   + u(1)*u(2)*Kxi(1,2) + u(2)*u(2)*Kxi(2,2)

      kS = Kxi(1,1)*Kxi(1,1) + Kxi(2,1)*Kxi(2,1)
     2   + Kxi(1,2)*Kxi(1,2) + Kxi(2,2)*Kxi(2,2)
      kS = ctC * kS * (mu/rho)**2._RKIND

      tauM = 1._RKIND / (rho * SQRT( kT + kU + kS )) ! JP 2021_05_04: the kT term represents the first term ( (2 * c1 / dt ) ** 2.0 ) in the tau_M eqn in eqn 2.8 (pg 15) of Mahdi's thesis; the kU term represents the 2nd term (u dotted with xi times u) in that same eqn; the kS represents the 3rd term (c2 * (mu / rho) ** 2.0 * inner product of xi and xi) in that same equation
      tauC = 1._RKIND / (tauM * (Kxi(1,1) + Kxi(2,2))) ! JP 2021_05_04: this line matches the eqn for tau_C in eqn 2.8 (pg 15) of Mahdi's thesis; where "(Kxi(1,1) + Kxi(2,2))" represents the trace of the matrix, xi

      rV(1) = ud(1) + u(1)*ux(1,1) + u(2)*ux(2,1)
        ! JP 2021_05_04: I think "u(1)*ux(1,1) + u(2)*ux(2,1)" is the dot product between the velocity and the gradient of the velocity
        ! JP 2021_05_04: I think rV is the residual of the PDE. Note that here, rV contains 1) the source term/body force (inside ud) 2) the time derivative of velocity (inside ud), 3) the convective term, but 4) is missing the diffusion term
                ! JP 2021_05_04: I think the diffusion term is not present because is Mahdi is using linear velocity and pressure interpolations (linear shape functions), so that means the diffusion term is zero (b/c the 2nd derivative of linear functions are zero)
      rV(2) = ud(2) + u(1)*ux(1,2) + u(2)*ux(2,2)

      up(1) = -tauM*(rho*rV(1) + px(1)) ! JP 2021_05_05: this line represents the negative tau_M times the PDE residual; issue - matches eqn 2.8 pg 15 in Mahdi's thesis for up (with the extra multiplication of the density, rho)
      up(2) = -tauM*(rho*rV(2) + px(2))


      ! JP 2021_04_22: I think "tauB" stands for tau_bar, where tau_bar is defined in eqn 2.8 in mahdi's thesis
      tauB = up(1)*up(1)*Kxi(1,1) + up(2)*up(1)*Kxi(2,1)
     2     + up(1)*up(2)*Kxi(1,2) + up(2)*up(2)*Kxi(2,2)
      IF (ISZERO(tauB)) tauB = eps
      tauB = rho/SQRT(tauB) ! JP 2021_05_05: this line exactly matches eqn 2.8 pg 15 in Mahdi's thesis (for tau_bar)

      ua(1) = u(1) + up(1)
      ua(2) = u(2) + up(2)
      pa    = p - tauC*divU

      rV(1) = tauB*(up(1)*ux(1,1) + up(2)*ux(2,1))
      rV(2) = tauB*(up(1)*ux(1,2) + up(2)*ux(2,2))

      rM(1,1) = mu*es(1,1) - rho*up(1)*ua(1) + rV(1)*up(1) - pa ! JP 2021_05_05: issue - I think the term "rho*up(1)*ua(1)" should not be multiplied by ua(1) but by u(1) instead and also i think there is an extra unneeded multiplication by rho in this same term
                 ! JP 2021_05_05: from my GoodNotes (written in light red), "mu*es(1,1)" and "- pa" account for terms 3 and 7; "rho*up(1)*ua(1)" accounts for term 6 (but with the error/issue mentioned above); "rV(1)*up(1)" accounts for term 9
      rM(2,1) = mu*es(2,1) - rho*up(1)*ua(2) + rV(1)*up(2)

      rM(1,2) = mu*es(1,2) - rho*up(2)*ua(1) + rV(2)*up(1)
      rM(2,2) = mu*es(2,2) - rho*up(2)*ua(2) + rV(2)*up(2) - pa

      rV(1) = ud(1) + ua(1)*ux(1,1) + ua(2)*ux(2,1) ! JP 2021_05_05: from my GoodNotes (written in light red), "ud(1)" accounts for terms 1, 4; "ua(1)*ux(1,1) + ua(2)*ux(2,1)" accounts for terms 2, 8
      rV(2) = ud(2) + ua(1)*ux(1,2) + ua(2)*ux(2,2)

      ! JP 2021_05_04: I think this section computes the local element residual vector
      DO a=1, eNoN
         uNx(a)  = u(1)*Nx(1,a)  + u(2)*Nx(2,a)
         upNx(a) = up(1)*Nx(1,a) + up(2)*Nx(2,a)
         uaNx(a) = uNx(a) + upNx(a)

         lR(1,a) = lR(1,a) + wr*N(a)*rV(1) + w*(Nx(1,a)*rM(1,1)
     2      + Nx(2,a)*rM(2,1))
                ! JP 2021_05_05: issue - I think the "wr*N(a)*rV(1)" term looks good, except that I think it is multiplying the up term by an extra unneeded rho term (where that rho term lives inside the "wr")
                ! JP 2021_05_05: The "wr*N(a)*rV(1)" term includes terms 1, 2, 4, 8 from my GoodNotes (written in light red)
                ! JP 2021_05_05: The "w*(Nx(1,a)*rM(1,1) + Nx(2,a)*rM(2,1))" term includes terms 3, 6, 7, 9 from my GoodNotes (written in light red)

         lR(2,a) = lR(2,a) + wr*N(a)*rV(2) + w*(Nx(1,a)*rM(1,2)
     2      + Nx(2,a)*rM(2,2))

         lR(3,a) = lR(3,a) + w*(N(a)*divU - upNx(a)) ! JP 2021_05_05: issue - this line matches eqn 2.10 pg 16 of Mahdi's thesis, for the continuity residual (also matches my eqn for the continuity residual in my GoodNotes), but I think the upNx term here is missing a division by density (rho).
                ! JP 2021_05_05: from my GoodNotes (written in light red), "N(a)*divU" accounts for term 10; "upNx(a)" accounts for term 11
      END DO

      ! JP 2021_05_04: I think this section computes the local element tangent matrix
      DO a=1, eNoN
         DO b=1, eNoN
            rM(1,1) = Nx(1,a)*Nx(1,b)
            rM(2,1) = Nx(2,a)*Nx(1,b)
            rM(1,2) = Nx(1,a)*Nx(2,b)
            rM(2,2) = Nx(2,a)*Nx(2,b)
                ! JP 2021_05_10: Here, let's match terms in "rM" with red-font terms in GoodNotes: "rM" = term 5 + 7

            NxNx = Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b)

            T1 = mu*NxNx + tauB*upNx(a)*upNx(b)
     2         + rho*( N(a)*(amd*N(b) + uaNx(b))
     3         + rho*tauM*uaNx(a)*(uNx(b) + amd*N(b)) )
                    ! JP 2021_05_10: T1 contains all the terms with the kronecker delta in my GoodNotes
                        ! JP 2021_05_10: Here, let's match terms in "T1" with red-font terms in GoodNotes: "mu*NxNx" = term 4; "tauB*upNx(a)*upNx(b)" = term 9; "rho*( N(a)*(amd*N(b)" = term 1; "rho*( N(a)*(" and "uaNx(b))" both = term 3 + 8; "rho*tauM*uaNx(a)*" and "amd*N(b))" both = term 2; "rho*tauM*uaNx(a)*(uNx(b)" = term 6
                            ! JP 2021_05_10: issue - note that in the "rho*tauM*uaNx(a)" term, there is an extra "up" term in this svFSI code, even though that extra "up" is not present in Madhi's eqn 2.18 for K or my GoodNotes)

            T2 = rho*tauM*uaNx(a)

            T3 = rho*tauM*(amd*N(b) + uNx(b))

!           dM/dU ! JP 2021_04_22: I think "M" stands for Momentum and "U" stands for velocity
            lK(1,a,b) = lK(1,a,b) + wl*((mu + tauC)*rM(1,1) + T1
     2         + mu_x*es_x(1,a)*es_x(1,b)) ! JP 2021_05_10: this line corresponds to i = 1, j = 1 (using 1-based indexing) in my GoodNotes
                        ! JP 2021_05_10: matches eqn 2.18 pg 18 (eqn for K^{ab}) in Mahdi's thesis, except for one difference, as mentioned in the "issue" above for the "rho*tauM*uaNx(a)" term
                        ! JP 2021_05_10: I think the "mu_x" term is for non-newtonian viscosity or some other kind of fancy viscosity that I don't need for the simple NS eqns that I wrote in my GoodNotes or the ones that Mahdi has in his thesis
            lK(2,a,b) = lK(2,a,b) + wl*(mu*rM(2,1) + tauC*rM(1,2)
     2         + mu_x*es_x(1,a)*es_x(2,b)) ! JP 2021_05_10: this line corresponds to i = 1, j = 2 (using 1-based indexing) in my GoodNotes

            lK(4,a,b) = lK(4,a,b) + wl*(mu*rM(1,2) + tauC*rM(2,1)
     2         + mu_x*es_x(2,a)*es_x(1,b)) ! JP 2021_05_10: this line corresponds to i = 2, j = 1 (using 1-based indexing) in my GoodNotes
            lK(5,a,b) = lK(5,a,b) + wl*((mu + tauC)*rM(2,2) + T1
     2         + mu_x*es_x(2,a)*es_x(2,b)) ! JP 2021_05_10: this line corresponds to i = 2, j = 2 (using 1-based indexing) in my GoodNotes

!           dM/dP
            lK(3,a,b) = lK(3,a,b) - wl*(Nx(1,a)*N(b) - Nx(1,b)*T2) ! JP 2021_05_10: issue - matches eqn 2.18 pg 18 (eqn for G^{ab}) in Mahdi's thesis, except for one difference. In "T2", there is an extra "up" term in this svFSI code, even though that extra "up" is not present in Madhi's eqn 2.18 for G.
            lK(6,a,b) = lK(6,a,b) - wl*(Nx(2,a)*N(b) - Nx(2,b)*T2)

!           dC/dU ! JP 2021_04_22: I think "C" stands for Continuity and "P" stands for Pressure
            lK(7,a,b) = lK(7,a,b) + wl*(N(a)*Nx(1,b) + Nx(1,a)*T3) ! JP 2021_05_10: matches eqn 2.18 pg 18 (eqn for D^{ab}) in Mahdi's thesis
            lK(8,a,b) = lK(8,a,b) + wl*(N(a)*Nx(2,b) + Nx(2,a)*T3)

!           dC/dP
            lK(9,a,b) = lK(9,a,b) + wl*(tauM*NxNx) ! JP 2021_05_10: matches eqn 2.18 pg 18 (eqn for L^{ab}) in Mahdi's thesis
         END DO
      END DO

      RETURN
      END SUBROUTINE FLUID2D
!####################################################################
      PURE SUBROUTINE BFLUID(eNoN, w, N, y, h, nV, lR, lK)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), y(tDof), h, nV(nsd)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      INTEGER(KIND=IKIND) a, b, i, j
      REAL(KIND=RKIND) T1, wl, hc(nsd), udn, u(nsd)

      wl  = w*eq(cEq)%af*eq(cEq)%gam*dt
      udn = 0._RKIND
      IF (mvMsh) THEN
         DO i=1, nsd
            j    = i + nsd + 1
            u(i) = y(i) - y(j)
            udn  = udn + u(i)*nV(i)
         END DO
      ELSE
         DO i=1, nsd
            u(i) = y(i)
            udn  = udn + u(i)*nV(i)
         END DO
      END IF

      udn = 0.5_RKIND*eq(cEq)%dmn(cDmn)%prop(backflow_stab)*
     2   eq(cEq)%dmn(cDmn)%prop(fluid_density)*(udn - ABS(udn))
      hc  = h*nV + udn*u

!     Here the loop is started for constructing left and right hand side
      IF (nsd .EQ. 2) THEN
         DO a=1, eNoN
            lR(1,a) = lR(1,a) - w*N(a)*hc(1)
            lR(2,a) = lR(2,a) - w*N(a)*hc(2)
            DO b=1, eNoN
               T1        = wl*N(a)*N(b)*udn
               lK(1,a,b) = lK(1,a,b) - T1
               lK(5,a,b) = lK(5,a,b) - T1
            END DO
         END DO
      ELSE
         DO a=1, eNoN
            lR(1,a) = lR(1,a) - w*N(a)*hc(1)
            lR(2,a) = lR(2,a) - w*N(a)*hc(2)
            lR(3,a) = lR(3,a) - w*N(a)*hc(3)
            DO b=1, eNoN
               T1 = wl*N(a)*N(b)*udn
               lK(1,a,b)  = lK(1,a,b)  - T1
               lK(6,a,b)  = lK(6,a,b)  - T1
               lK(11,a,b) = lK(11,a,b) - T1
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE BFLUID
!####################################################################
      SUBROUTINE BWFLUID3D(eNoN, w, N, Nx, yl, ub, nV, tauB, lR, lK) ! JP 2021_05_04: I think the "W" in "BWFLUID3D" stands for "weak" and the "B" stands for "boudary", so I think this function is used to apply the BCs (in a weakly enforced sense)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(3,eNoN),
     2   yl(tDof,eNoN), ub(3), nV(3), tauB(2)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      INTEGER(KIND=IKIND) :: a, b
      REAL(KIND=RKIND) :: rho, nu, T1, wr, wl, wrl, tauT, tauN, p, uhn,
     2   ubn, u(3), uh(3), ux(3,3), sigman(3), Nxn(eNoN), rV(3),
     3   rM(3,3), nu_s, es(3,3), gam, nu_x

      rho  = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      T1   = eq(cEq)%af * eq(cEq)%gam * dt
      wr   = w * rho
      wl   = w * T1
      wrl  = wr * T1
      tauT = tauB(1) / rho
      tauN = tauB(2) / rho

      p    = 0._RKIND
      u    = 0._RKIND
      ux   = 0._RKIND
      DO a=1, eNoN
         p = p + N(a)*yl(4,a)

         u(1) = u(1) + N(a)*yl(1,a)
         u(2) = u(2) + N(a)*yl(2,a)
         u(3) = u(3) + N(a)*yl(3,a)

         ux(1,1) = ux(1,1) + Nx(1,a)*yl(1,a)
         ux(2,1) = ux(2,1) + Nx(2,a)*yl(1,a)
         ux(3,1) = ux(3,1) + Nx(3,a)*yl(1,a)
         ux(1,2) = ux(1,2) + Nx(1,a)*yl(2,a)
         ux(2,2) = ux(2,2) + Nx(2,a)*yl(2,a)
         ux(3,2) = ux(3,2) + Nx(3,a)*yl(2,a)
         ux(1,3) = ux(1,3) + Nx(1,a)*yl(3,a)
         ux(2,3) = ux(2,3) + Nx(2,a)*yl(3,a)
         ux(3,3) = ux(3,3) + Nx(3,a)*yl(3,a)

         Nxn(a)  = Nx(1,a)*nV(1) + Nx(2,a)*nV(2) + Nx(3,a)*nV(3)
      END DO

      uh = 0._RKIND
      IF (mvMsh) THEN
         DO a=1, eNoN
            uh(1) = uh(1) + N(a)*yl(5,a)
            uh(2) = uh(2) + N(a)*yl(6,a)
            uh(3) = uh(3) + N(a)*yl(7,a)
         END DO
      END IF
      ubn = (u(1)-ub(1))*nV(1) + (u(2)-ub(2))*nV(2) + (u(3)-ub(3))*nV(3)
      uhn = (u(1)-uh(1))*nV(1) + (u(2)-uh(2))*nV(2) + (u(3)-uh(3))*nV(3)
      uhn = (ABS(uhn) - uhn) * 0.5_RKIND

!     Strain rate tensor 2*e_ij := (u_ij + u_ji)
      es(1,1) = ux(1,1) + ux(1,1)
      es(2,2) = ux(2,2) + ux(2,2)
      es(3,3) = ux(3,3) + ux(3,3)

      es(1,2) = ux(1,2) + ux(2,1)
      es(1,3) = ux(1,3) + ux(3,1)
      es(2,3) = ux(2,3) + ux(3,2)

      es(2,1) = es(1,2)
      es(3,1) = es(1,3)
      es(3,2) = es(2,3)

!     Shear-rate := (2*e_ij*e_ij)^.5
      gam = es(1,1)*es(1,1) + es(1,2)*es(1,2) + es(1,3)*es(1,3) +
     2      es(2,1)*es(2,1) + es(2,2)*es(2,2) + es(2,3)*es(2,3) +
     3      es(3,1)*es(3,1) + es(3,2)*es(3,2) + es(3,3)*es(3,3)
      gam = SQRT(0.5_RKIND*gam)

!     Compute viscosity based on shear-rate and chosen viscosity model
      CALL GETVISCOSITY(eq(cEq)%dmn(cDmn), gam, nu, nu_s, nu_x)
      nu = nu/rho

      sigman(1) = nu*(es(1,1)*nV(1) + es(1,2)*nV(2) + es(1,3)*nV(3))
     2   - (p/rho)*nV(1)
      sigman(2) = nu*(es(2,1)*nV(1) + es(2,2)*nV(2) + es(2,3)*nV(3))
     2   - (p/rho)*nV(2)
      sigman(3) = nu*(es(3,1)*nV(1) + es(3,2)*nV(2) + es(3,3)*nV(3))
     2   - (p/rho)*nV(3)

      rV(1) = -sigman(1) + (tauT + uhn)*(u(1)-ub(1))
      rV(2) = -sigman(2) + (tauT + uhn)*(u(2)-ub(2))
      rV(3) = -sigman(3) + (tauT + uhn)*(u(3)-ub(3))

      DO a=1, eNoN
c         lR(1,a) = lR(1,a) + wr*( N(a)*(rV(1) + (tauT-tauN)*ubn*nV(1)) -
c     2      nu*(Nxn(a)*(u(1)-ub(1)) + Nx(1,a)*ubn) )
         lR(1,a) = lR(1,a) + wr*( N(a)*(rV(1) + (tauT-tauN)*ubn*nV(1)) -
     2      nu*(u(1)-ub(1))*(Nxn(a) + Nx(1,a)*nV(1)) )

c         lR(2,a) = lR(2,a) + wr*( N(a)*(rV(2) + (tauT-tauN)*ubn*nV(2)) -
c     2      nu*(Nxn(a)*(u(2)-ub(2)) + Nx(2,a)*ubn) )
         lR(2,a) = lR(2,a) + wr*( N(a)*(rV(2) + (tauT-tauN)*ubn*nV(2)) -
     2      nu*(u(2)-ub(2))*(Nxn(a) + Nx(2,a)*nV(2)) )

c         lR(3,a) = lR(3,a) + wr*( N(a)*(rV(3) + (tauT-tauN)*ubn*nV(3)) -
c     2      nu*(Nxn(a)*(u(3)-ub(3)) + Nx(3,a)*ubn) )
         lR(3,a) = lR(3,a) + wr*( N(a)*(rV(3) + (tauT-tauN)*ubn*nV(3)) -
     2      nu*(u(3)-ub(3))*(Nxn(a) + Nx(3,a)*nV(3)) )

         lR(4,a) = lR(4,a) - w*N(a)*ubn
      END DO

      DO a=1, eNoN
         DO b=1, eNoN
            T1 = (tauT + uhn)*N(a)*N(b) - nu*(N(a)*Nxn(b) + N(b)*Nxn(a))

            rM(1,1) = (tauT-tauN)*N(a)*N(b)*nV(1)*nV(1) -
     2         nu*(N(a)*Nx(1,b)*nV(1) + Nx(1,a)*N(b)*nV(1))
c            rM(1,2) = (tauT-tauN)*N(a)*N(b)*nV(1)*nV(2) -
c     2         nu*(N(a)*Nx(1,b)*nV(2) + Nx(1,a)*N(b)*nV(2))
c            rM(1,3) = (tauT-tauN)*N(a)*N(b)*nV(1)*nV(3) -
c     2         nu*(N(a)*Nx(1,b)*nV(3) + Nx(1,a)*N(b)*nV(3))
            rM(1,2) = (tauT-tauN)*N(a)*N(b)*nV(1)*nV(2) -
     2         nu*N(a)*Nx(1,b)*nV(2)
            rM(1,3) = (tauT-tauN)*N(a)*N(b)*nV(1)*nV(3) -
     2         nu*N(a)*Nx(1,b)*nV(3)

c            rM(2,1) = (tauT-tauN)*N(a)*N(b)*nV(2)*nV(1) -
c     2         nu*(N(a)*Nx(2,b)*nV(1) + Nx(2,a)*N(b)*nV(1))
            rM(2,1) = (tauT-tauN)*N(a)*N(b)*nV(2)*nV(1) -
     2         nu*N(a)*Nx(2,b)*nV(1)
            rM(2,2) = (tauT-tauN)*N(a)*N(b)*nV(2)*nV(2) -
     2         nu*(N(a)*Nx(2,b)*nV(2) + Nx(2,a)*N(b)*nV(2))
c            rM(2,3) = (tauT-tauN)*N(a)*N(b)*nV(1)*nV(3) -
c     2         nu*(N(a)*Nx(2,b)*nV(3) + Nx(2,a)*N(b)*nV(3))
            rM(2,3) = (tauT-tauN)*N(a)*N(b)*nV(1)*nV(3) -
     2         nu*N(a)*Nx(2,b)*nV(3)

c            rM(3,1) = (tauT-tauN)*N(a)*N(b)*nV(3)*nV(1) -
c     2         nu*(N(a)*Nx(3,b)*nV(1) + Nx(3,a)*N(b)*nV(1))
c            rM(3,2) = (tauT-tauN)*N(a)*N(b)*nV(3)*nV(2) -
c     2         nu*(N(a)*Nx(3,b)*nV(2) + Nx(3,a)*N(b)*nV(2))
            rM(3,1) = (tauT-tauN)*N(a)*N(b)*nV(3)*nV(1) -
     2         nu*N(a)*Nx(3,b)*nV(1)
            rM(3,2) = (tauT-tauN)*N(a)*N(b)*nV(3)*nV(2) -
     2         nu*N(a)*Nx(3,b)*nV(2)
            rM(3,3) = (tauT-tauN)*N(a)*N(b)*nV(3)*nV(3) -
     2         nu*(N(a)*Nx(3,b)*nV(3) + Nx(3,a)*N(b)*nV(3))

!           dM/dU
            lK(1,a,b) = lK(1,a,b) + wrl*(T1 + rM(1,1))
            lK(2,a,b) = lK(2,a,b) + wrl*rM(1,2)
            lK(3,a,b) = lK(3,a,b) + wrl*rM(1,3)

            lK(5,a,b) = lK(5,a,b) + wrl*rM(2,1)
            lK(6,a,b) = lK(6,a,b) + wrl*(T1 + rM(2,2))
            lK(7,a,b) = lK(7,a,b) + wrl*rM(2,3)

            lK(9,a,b)  = lK(9,a,b)  + wrl*rM(3,1)
            lK(10,a,b) = lK(10,a,b) + wrl*rM(3,2)
            lK(11,a,b) = lK(11,a,b) + wrl*(T1 + rM(3,3))

!           dM/dP
            lK(4,a,b)  = lK(4,a,b)  + wl*N(a)*N(b)*nV(1)
            lK(8,a,b)  = lK(8,a,b)  + wl*N(a)*N(b)*nV(2)
            lK(12,a,b) = lK(12,a,b) + wl*N(a)*N(b)*nV(3)

!           dC/dU
            lK(13,a,b) = lK(13,a,b) - wl*N(a)*N(b)*nV(1)
            lK(14,a,b) = lK(14,a,b) - wl*N(a)*N(b)*nV(2)
            lK(15,a,b) = lK(15,a,b) - wl*N(a)*N(b)*nV(3)
         END DO
      END DO

      RETURN
      END SUBROUTINE BWFLUID3D
!--------------------------------------------------------------------
      SUBROUTINE BWFLUID2D(eNoN, w, N, Nx, yl, ub, nV, tauB, lR, lK) ! JP 2021_05_04: I think the "W" in "BWFLUID3D" stands for "weak" and the "B" stands for "boudary", so I think this function is used to apply the BCs (in a weakly enforced sense)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(2,eNoN),
     2   yl(tDof,eNoN), ub(2), nV(2), tauB(2)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      INTEGER(KIND=IKIND) :: a, b
      REAL(KIND=RKIND) :: rho, nu, T1, wr, wl, wrl, tauT, tauN, p, uhn,
     2   ubn, u(2), uh(2), ux(2,2), sigman(2), Nxn(eNoN), rV(2),
     3   rM(2,2), nu_s, es(2,2), gam, nu_x

      rho  = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      T1   = eq(cEq)%af * eq(cEq)%gam * dt
      wr   = w*rho
      wl   = w  * T1
      wrl  = wr * T1
      tauT = tauB(1) / rho
      tauN = tauB(2) / rho

      p    = 0._RKIND
      u    = 0._RKIND
      ux   = 0._RKIND
      DO a=1, eNoN
         p = p + N(a)*yl(3,a)

         u(1) = u(1) + N(a)*yl(1,a)
         u(2) = u(2) + N(a)*yl(2,a)

         ux(1,1) = ux(1,1) + Nx(1,a)*yl(1,a)
         ux(2,1) = ux(2,1) + Nx(2,a)*yl(1,a)
         ux(1,2) = ux(1,2) + Nx(1,a)*yl(2,a)
         ux(2,2) = ux(2,2) + Nx(2,a)*yl(2,a)

         Nxn(a)  = Nx(1,a)*nV(1) + Nx(2,a)*nV(2)
      END DO

      uh = 0._RKIND
      IF (mvMsh) THEN
         DO a=1, eNoN
            uh(1) = uh(1) + N(a)*yl(4,a)
            uh(2) = uh(2) + N(a)*yl(5,a)
         END DO
      END IF
      uhn = (u(1)-uh(1))*nV(1) + (u(2)-uh(2))*nV(2)
      ubn = (u(1)-ub(1))*nV(1) + (u(2)-ub(2))*nV(2)
      uhn = (ABS(uhn) - uhn) * 0.5_RKIND

!     Strain rate tensor 2*e_ij := (u_ij + u_ji)
      es(1,1) = ux(1,1) + ux(1,1)
      es(2,2) = ux(2,2) + ux(2,2)

      es(1,2) = ux(1,2) + ux(2,1)
      es(2,1) = es(1,2)

!     Shear-rate := (2*e_ij*e_ij)^.5
      gam = es(1,1)*es(1,1) + es(1,2)*es(1,2) + es(2,1)*es(2,1)
     2    + es(2,2)*es(2,2)
      gam = SQRT(0.5_RKIND*gam)

!     Compute viscosity based on shear-rate and chosen viscosity model
      CALL GETVISCOSITY(eq(cEq)%dmn(cDmn), gam, nu, nu_s, nu_x)
      nu = nu/rho

      sigman(1) = -(p/rho)*nV(1) + nu*(es(1,1)*nV(1) + es(1,2)*nV(2))
      sigman(2) = -(p/rho)*nV(2) + nu*(es(2,1)*nV(1) + es(2,2)*nV(2))

      rV(1) = -sigman(1) + (tauT + uhn)*(u(1)-ub(1))
      rV(2) = -sigman(2) + (tauT + uhn)*(u(2)-ub(2))

      DO a=1, eNoN
c         lR(1,a) = lR(1,a) + wr*( N(a)*(rV(1) + (tauT-tauN)*ubn*nV(1)) -
c     2      nu*(Nxn(a)*(u(1)-ub(1)) + Nx(1,a)*ubn) )
         lR(1,a) = lR(1,a) + wr*( N(a)*(rV(1) + (tauT-tauN)*ubn*nV(1)) -
     2      nu*(u(1)-ub(1))*(Nxn(a) + Nx(1,a)*nV(1)) )

c         lR(2,a) = lR(2,a) + wr*( N(a)*(rV(2) + (tauT-tauN)*ubn*nV(2)) -
c     2      nu*(Nxn(a)*(u(2)-ub(2)) + Nx(2,a)*ubn) )
         lR(2,a) = lR(2,a) + wr*( N(a)*(rV(2) + (tauT-tauN)*ubn*nV(2)) -
     2      nu*(u(2)-ub(2))*(Nxn(a) + Nx(2,a)*nV(2)) )

         lR(3,a) = lR(3,a) - w*N(a)*ubn
      END DO

      DO a=1, eNoN
         DO b=1, eNoN
            T1 = (tauT + uhn)*N(a)*N(b) - nu*(N(a)*Nxn(b) + N(b)*Nxn(a))

            rM(1,1) = (tauT-tauN)*N(a)*N(b)*nV(1)*nV(1) -
     2         nu*(N(a)*Nx(1,b)*nV(1) + Nx(1,a)*N(b)*nV(1))
c            rM(1,2) = (tauT-tauN)*N(a)*N(b)*nV(1)*nV(2) -
c     2         nu*(N(a)*Nx(1,b)*nV(2) + Nx(1,a)*N(b)*nV(2))
            rM(1,2) = (tauT-tauN)*N(a)*N(b)*nV(1)*nV(2) -
     2         nu*N(a)*Nx(1,b)*nV(2)

c            rM(2,1) = (tauT-tauN)*N(a)*N(b)*nV(2)*nV(1) -
c     2         nu*(N(a)*Nx(2,b)*nV(1) + Nx(2,a)*N(b)*nV(1))
            rM(2,1) = (tauT-tauN)*N(a)*N(b)*nV(2)*nV(1) -
     2         nu*N(a)*Nx(2,b)*nV(1)
            rM(2,2) = (tauT-tauN)*N(a)*N(b)*nV(2)*nV(2) -
     2         nu*(N(a)*Nx(2,b)*nV(2) + Nx(2,a)*N(b)*nV(2))

!           dM/dU
            lK(1,a,b) = lK(1,a,b) + wrl*(T1 + rM(1,1))
            lK(2,a,b) = lK(2,a,b) + wrl*rM(1,2)

            lK(4,a,b) = lK(4,a,b) + wrl*rM(2,1)
            lK(5,a,b) = lK(5,a,b) + wrl*(T1 + rM(2,2))

!           dM/dP
            lK(3,a,b) = lK(4,a,b) + wl*N(a)*N(b)*nV(1)
            lK(6,a,b) = lK(8,a,b) + wl*N(a)*N(b)*nV(2)

!           dC/dU
            lK(7,a,b) = lK(7,a,b) - wl*N(a)*N(b)*nV(1)
            lK(8,a,b) = lK(8,a,b) - wl*N(a)*N(b)*nV(2)
         END DO
      END DO

      RETURN
      END SUBROUTINE BWFLUID2D
!####################################################################
      SUBROUTINE GETVISCOSITY(lDmn, gamma, mu, mu_s, mu_x)
      USE COMMOD
      IMPLICIT NONE
      TYPE(dmnType), INTENT(IN) :: lDmn
      REAL(KIND=RKIND), INTENT(INOUT)  :: gamma
      REAL(KIND=RKIND), INTENT(OUT) :: mu, mu_s, mu_x

      REAL(KIND=RKIND) :: mu_i, mu_o, lam, a, n, T1, T2

      SELECT CASE (lDmn%visc%viscType)
      CASE (viscType_Const)
         mu   = lDmn%visc%mu_i
         mu_s = mu
         mu_x = 0._RKIND

      CASE (viscType_CY)
         mu_i = lDmn%visc%mu_i
         mu_o = lDmn%visc%mu_o
         lam  = lDmn%visc%lam
         a    = lDmn%visc%a
         n    = lDmn%visc%n

         T1   = 1._RKIND + (lam*gamma)**a
         T2   = T1**((n-1._RKIND)/a)
         mu   = mu_i + (mu_o-mu_i)*T2
         mu_s = mu_i

         T1   = T2/T1
         T2   = lam**a * gamma**(a-1._RKIND) * T1
         mu_x = (mu_o-mu_i)*(n-1._RKIND)*T2

      CASE (viscType_Cass)
         mu_i = lDmn%visc%mu_i
         mu_o = lDmn%visc%mu_o
         lam  = lDmn%visc%lam

         IF (gamma .LT. lam) THEN
            mu_o  = mu_o/SQRT(lam)
            gamma = lam
         ELSE
            mu_o = mu_o/SQRT(gamma)
         END IF
         mu   = (mu_i + mu_o) * (mu_i + mu_o)
         mu_s = mu_i*mu_i
         mu_x = 2._RKIND*mu_o*(mu_o + mu_i)/gamma

      END SELECT

      RETURN
      END SUBROUTINE GETVISCOSITY
!####################################################################

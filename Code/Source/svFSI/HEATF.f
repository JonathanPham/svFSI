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
!     This is for solving heat equation in a fluid (stabilized
!     advection-diffusion equation)
!
!--------------------------------------------------------------------

      SUBROUTINE CONSTRUCT_HEATF(lM, Ag, Yg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo)

      INTEGER(KIND=IKIND) a, e, g, Ac, eNoN, cPhys
      REAL(KIND=RKIND) w, Jac, ksix(nsd,nsd)

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), al(:,:), yl(:,:),
     2   N(:), Nx(:,:), lR(:,:), lK(:,:,:)

      eNoN = lM%eNoN

!     HEATF: dof = 1
      ALLOCATE(ptr(eNoN), xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN),
     2   N(eNoN), Nx(nsd,eNoN), lR(dof,eNoN), lK(dof*dof,eNoN,eNoN))

!     Loop over all elements of mesh
      DO e=1, lM%nEl
!        Update domain and proceed if domain phys and eqn phys match
         cDmn  = DOMAIN(lM, cEq, e)
         cPhys = eq(cEq)%dmn(cDmn)%phys
         IF (cPhys .NE. phys_heatF) CYCLE

!        Update shape functions for NURBS
         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

!        Create local copies
         DO a=1, eNoN
            Ac = lM%IEN(a,e)
            ptr(a)  = Ac
            xl(:,a) = x(:,Ac)
            al(:,a) = Ag(:,Ac)
            yl(:,a) = Yg(:,Ac)
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
               CALL HEATF3D(eNoN, w, N, Nx, al, yl, ksix, lR, lK)

            ELSE IF (nsd .EQ. 2) THEN
               CALL HEATF2D(eNoN, w, N, Nx, al, yl, ksix, lR, lK)

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

      DEALLOCATE(ptr, xl, al, yl, N, Nx, lR, lK)

      RETURN
      END SUBROUTINE CONSTRUCT_HEATF
!####################################################################
!     This is for solving heat equation in a fluid.
      PURE SUBROUTINE HEATF3D (eNoN, w, N, Nx, al, yl, ksix, lR, lK) ! JP 2021_04_21: see my notes in the HEATF2D subroutine to understand what the HEATF3D subroutine does. These two subroutines basically compute the same thing, except one does it for 2d domain problems and the other for 3d.
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(nsd,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN), ksix(nsd,nsd)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(1,eNoN), lK(1,eNoN,eNoN)

      REAL(KIND=RKIND), PARAMETER :: ct(4) = (/4._RKIND, 1._RKIND,
     2   3._RKIND, 1._RKIND/)

      INTEGER(KIND=IKIND) i, a, b
      REAL(KIND=RKIND) nu, tauM, kU, kS, nTx, Tp, udTx, udNx(eNoN), T1,
     2   amd, wl, Td, Tx(nsd), u(nsd), s

      T1  = eq(cEq)%af*eq(cEq)%gam*dt
      amd = eq(cEq)%am/T1
      nu  = eq(cEq)%dmn(cDmn)%prop(conductivity)
      s   = eq(cEq)%dmn(cDmn)%prop(source_term)
      i   = eq(cEq)%s
      wl  = w*T1

      u  = 0._RKIND
      Td = -s
      Tx = 0._RKIND
      DO a=1, eNoN
         u(1) = u(1) + N(a)*yl(1,a)
         u(2) = u(2) + N(a)*yl(2,a)
         u(3) = u(3) + N(a)*yl(3,a)

         Td = Td + N(a)*al(i,a)

         Tx(1) = Tx(1) + Nx(1,a)*yl(i,a)
         Tx(2) = Tx(2) + Nx(2,a)*yl(i,a)
         Tx(3) = Tx(3) + Nx(3,a)*yl(i,a)
      END DO

      IF (mvMsh) THEN
         DO a=1, eNoN
            u(1) = u(1) - N(a)*yl(5,a)
            u(2) = u(2) - N(a)*yl(6,a)
            u(3) = u(3) - N(a)*yl(7,a)
         END DO
      END IF

      kU = u(1)*u(1)*ksix(1,1) + u(2)*u(1)*ksix(2,1)
     2   + u(3)*u(1)*ksix(3,1) + u(1)*u(2)*ksix(1,2)
     3   + u(2)*u(2)*ksix(2,2) + u(3)*u(2)*ksix(3,2)
     4   + u(1)*u(3)*ksix(1,3) + u(2)*u(3)*ksix(2,3)
     5   + u(3)*u(3)*ksix(3,3)

      kS = ksix(1,1)*ksix(1,1) + ksix(2,1)*ksix(2,1)
     2   + ksix(3,1)*ksix(3,1) + ksix(1,2)*ksix(1,2)
     3   + ksix(2,2)*ksix(2,2) + ksix(3,2)*ksix(3,2)
     4   + ksix(1,3)*ksix(1,3) + ksix(2,3)*ksix(2,3)
     5   + ksix(3,3)*ksix(3,3)

      nTx = ksix(1,1)*Tx(1)*Tx(1)
     2    + ksix(2,2)*Tx(2)*Tx(2)
     3    + ksix(3,3)*Tx(3)*Tx(3)
     4    + (ksix(1,2) + ksix(2,1))*Tx(1)*Tx(2)
     5    + (ksix(1,3) + ksix(3,1))*Tx(1)*Tx(3)
     6    + (ksix(2,3) + ksix(3,2))*Tx(2)*Tx(3)

      IF (ISZERO(nTx)) nTx = eps

      udTx = u(1)*Tx(1) + u(2)*Tx(2) + u(3)*Tx(3)
      Tp   = ABS(Td + udTx)
!      nu   = nu + kDC
      nu   = nu + 0.5_RKIND*Tp/SQRT(nTx)
      tauM = ct(4)/SQRT((ct(1)/(dt*dt)) + ct(2)*kU + ct(3)*nu*nu*kS)
      Tp   = -tauM*(Td + udTx)

      DO a=1,eNoN
         udNx(a) = u(1)*Nx(1,a) + u(2)*Nx(2,a) + u(3)*Nx(3,a)
      END DO

      DO a=1,eNoN
         lR(1,a) = lR(1,a) + w*(N(a)*(Td + udTx) + (Nx(1,a)*Tx(1)
     2      + Nx(2,a)*Tx(2) + Nx(3,a)*Tx(3))*nu - udNx(a)*Tp)

         DO b=1,eNoN
            lK(1,a,b) = lK(1,a,b) + wl*(nu*(Nx(1,a)*Nx(1,b)
     2         + Nx(2,a)*Nx(2,b) + Nx(3,a)*Nx(3,b))
     3         + (N(a) + tauM*udNx(a))*(N(b)*amd + udNx(b)))
         END DO
      END DO

      RETURN
      END SUBROUTINE HEATF3D
!--------------------------------------------------------------------
!     This is for solving heat equation in a fluid

! JP 2021_04_08: i think the original author of this advection-diffusion FEM solver (HEATF.f) was mahdi, so look at his paper to understand the weak form and FEM formulation etc. The relevant paper is: "A non-discrete method for computation of residence time in fluid mechanics simulations". Physics of Fluids. 2013. See eqn 5 in this paper for the weak form?
      PURE SUBROUTINE HEATF2D (eNoN, w, N, Nx, al, yl, ksix, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(nsd,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN), ksix(nsd,nsd)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(1,eNoN), lK(1,eNoN,eNoN)

      REAL(KIND=RKIND), PARAMETER :: ct(4) = (/4._RKIND, 1._RKIND,
     2   3._RKIND, 1._RKIND/)

      INTEGER(KIND=IKIND) i, a, b
      REAL(KIND=RKIND) nu, tauM, kU, kS, nTx, Tp, udTx, udNx(eNoN), T1,
     2   amd, wl, Td, Tx(nsd), u(nsd), s

      T1  = eq(cEq)%af*eq(cEq)%gam*dt
      amd = eq(cEq)%am/T1
      nu  = eq(cEq)%dmn(cDmn)%prop(conductivity)
      s   = eq(cEq)%dmn(cDmn)%prop(source_term)
      i   = eq(cEq)%s
      wl  = w*T1

      u  = 0._RKIND
      Td = -s
      Tx = 0._RKIND
      DO a=1, eNoN
         u(1) = u(1) + N(a)*yl(1,a) ! JP 2021_04_21: i think this line computes the x-component of the velocity vector
         u(2) = u(2) + N(a)*yl(2,a) ! JP 2021_04_21: i think this line computes the y-component of the velocity vector

         Td = Td + N(a)*al(i,a) ! JP 2021_04_21: shape function multiplied by the time derivative of the degree of freedom

         Tx(1) = Tx(1) + Nx(1,a)*yl(i,a) ! JP 2021_04_21: spatial derivative of the shape function multiplied by the the degree of freedom
         Tx(2) = Tx(2) + Nx(2,a)*yl(i,a)
      END DO

      IF (mvMsh) THEN ! JP 2021_04_21: according to MOD.f, mvMsh is a boolean that indicates "Whether mesh is moving"; I think this parameter is used only in FSI cases (see READFILES.f with the line "CASE ('FSI')", where this section sets mvMsh to be true (mvMsh has a value of false by default))
         DO a=1, eNoN
            u(1) = u(1) - N(a)*yl(5,a)
            u(2) = u(2) - N(a)*yl(6,a)
            u(3) = u(3) - N(a)*yl(7,a)
         END DO
      END IF

      ! JP 2021_04_21: I believe the below kU, kS, and nTx parameters are used to compute the kappa_tilda and tau_M terms
      kU = u(1)*u(1)*ksix(1,1) + u(2)*u(1)*ksix(2,1)
     2   + u(1)*u(2)*ksix(1,2) + u(2)*u(2)*ksix(2,2)

      kS = ksix(1,1)*ksix(1,1) + ksix(2,1)*ksix(2,1)
     2   + ksix(1,2)*ksix(1,2) + ksix(2,2)*ksix(2,2)

      nTx = ksix(1,1)*Tx(1)*Tx(1) + ksix(2,2)*Tx(2)*Tx(2)
     4    + (ksix(1,2) + ksix(2,1))*Tx(1)*Tx(2)

      IF (ISZERO(nTx)) nTx = eps

      udTx = u(1)*Tx(1) + u(2)*Tx(2) ! JP 2021_04_21: I think this line computes everything in term 2 (blue) of my GoodNotes derivation (the advection matrix) (but without the N_{a}^{e} term). This term is the dot product between the velocity vector, u, and the gradient of the -the shape function multiplied by the degree of freedom soltn-, a.k.a. grad(N_{b}^{e} * d_{b}^{e})
      Tp   = ABS(Td + udTx) ! JP 2021_04_21: Tp is the norm of the residual here
!      nu   = nu + kDC
      nu   = nu + 0.5_RKIND*Tp/SQRT(nTx) ! JP 2021_04_19: I think this line refers to eqn 7 of esmaily 2013, where here, Tp stands for the residual, L (from eqn 7).
      tauM = ct(4)/SQRT((ct(1)/(dt*dt)) + ct(2)*kU + ct(3)*nu*nu*kS) ! JP 2021_04_19: this term mostly matches eqn 6 of esmaily 2013
      Tp   = -tauM*(Td + udTx) ! JP 2021_04_19: I think this line represents the negative of tau_m multiplied by the residual, L ("Td + udTx"), where this entire term is eqn in the last integral term of eqn 5 of Mahdi's residence time paper (esmaily 2013). However, note that the residual, as represented by "Td + udTx", is missing the diffusion term tho (it only accounts for 1) the multiplication of the shape function by the acceleration term (the time derivative of the "d" degree of freedom) and 2) the convection term). So note that this "Tp" term, which is again missing the diffusion term is used to compute the residual that is solved by gen alpha (which is different from the PDE residual); so my question is: why are we missing the diffusion term here in the PDE residual, "Tp". Dont we need the diffusion term term in Tp to compute the correct gen alpha residual? Or maybe it is possible that we are using equal-order velocity and pressure interpolation for the FEM of the Navier-Stokes eqn (where I believe we use first-order interpolation (shape functions) for both velocity and pressure); in this case, if we are indeed using first-order shape functions, then the diffusion term is zero because the diffusion term is a 2nd derivative and the 2nd derivative of a first-order (linear) shape function is zero. In Mahdi's residence time paper, he does say that he is using "equal-order velocity and pressure interpolation"

      DO a=1, eNoN
         udNx(a) = u(1)*Nx(1,a) + u(2)*Nx(2,a) ! JP 2021: "udNx" is = the dot product between the velocity vector, u, and the grad(N_{a}^{e})
      END DO

      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*(N(a)*(Td + udTx)
     2      + (Nx(1,a)*Tx(1) + Nx(2,a)*Tx(2))*nu - udNx(a)*Tp) ! JP 2021_04_19: based on the above comment for Tp, I think udNx(a) represents the grad(q)*u term in eqn 5 of esmaily 2013

                    ! JP 2021_04_21: The term "N(a)*(Td + udTx)" computes term 1 and term 6 and term 2 in my GoodNotes derivation (in blue) (basically the mass matrix and source term contributions, and the advection martix contribution), where the "Td" term includes the contributions from the mass matrix and source term and the "udTx" term includes the contribution from the advection matrix

                    ! JP 2021_04_21: the "(Nx(1,a)*Tx(1) + Nx(2,a)*Tx(2))*nu" term represents the stiffness matrix (term 3 (blue) in my GoodNotes), where the "Tx" is = N_{b}^{e} * d_{b]^{e}, summed over all b's from 1 to eNoN (where eNoN = n_{en} I think)

                    ! JP 2021_04_21: "- udNx(a)*Tp" represents the SUPG stabilization integral term, where "Tp" is = - tau_m * L (where L is just the residual, not the norm of the residual) and "udNx" is = the dot product between the velocity vector, u, and the grad(N_{a}^{e})

         DO b=1, eNoN
            lK(1,a,b) = lK(1,a,b) + wl*(nu*(Nx(1,a)*Nx(1,b)
     2         + Nx(2,a)*Nx(2,b))
     3         + (N(a) + tauM*udNx(a))*(N(b)*amd + udNx(b)))

                    ! JP 2021_04_21: In the line, "tauM*udNx(a)", "udNx" is = the dot product between the velocity vector, u, and the grad( N_{a}^{e} ).

                    ! JP 2021_04_21: "(N(a) + tauM*udNx(a))*(N(b)*amd" represents the derivative of the residual (used in gen alpha, not the PDE residual) wrt the acceleration (the time derivative of the "d" degree of freedom) evaluated at n+alpha_m; this corresponds to term 1 and 2 in GoodNotes (in orange writing).

                    ! JP 2021_04_21: "(N(a) + tauM*udNx(a))" * "udNx(b))" represents the term 3 and term 5 (orange) in GoodNotes, except for one small difference which is this code is missing the diffusion term tho. I think this might be related to 1) the idea of using an inconsistent tangent matrix so we dont need to compute a fully correct/exact tangent matrix and we can still hit convergence in our newton solver or 2) the idea of using "equal-order velocity and pressure interpolation" (and first-order linear shape functions for velocity) (see above discussion for more details (locatied where we compute the term, "Tp   = -tauM*(Td + udTx)"))

                    ! JP 2021_04_21: the term, "nu*(Nx(1,a)*Nx(1,b)" + "Nx(2,a)*Nx(2,b))", represents term 4 (orange) in GoodNotes, where this term corresponds to the stiffness matrix contribution to the tangent matrix. Specifically, "Nx(1,a)*Nx(1,b)" + "Nx(2,a)*Nx(2,b)" is the dot product between the grad(N_{a}^{e}) and grad(N_{b}^{e})
         END DO
      END DO

      RETURN
      END SUBROUTINE HEATF2D
!####################################################################
      PURE SUBROUTINE BHEATF (eNoN, w, N, y, h, nV, lR, lK)
      USE COMMOD
      IMPLICIT NONE

      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), y(tDof), h, nV(nsd)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      INTEGER(KIND=IKIND) a, b, i
      REAL(KIND=RKIND) T1, wl, T, udn

      wl = w*eq(cEq)%af*eq(cEq)%gam*dt
      T  = y(eq(cEq)%s)

      udn = 0._RKIND
      DO i=1, nsd
         udn = udn + y(i)*nV(i)
      END DO
      udn = 0.5_RKIND*(udn - ABS(udn))
      T1  = h - udn*T

      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*N(a)*T1
         DO b=1, eNoN
            lK(1,a,b) = lK(1,a,b) - wl*N(a)*N(b)*udn
         END DO
      END DO

      RETURN
      END SUBROUTINE BHEATF
!####################################################################

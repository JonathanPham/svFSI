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
!     This is for solving heat equation in a solid (simple diffusion/
!     Laplace equation)
!
!--------------------------------------------------------------------

! JP: I think this code uses fortran 77 style. If you look at https://web.stanford.edu/class/me200c/tutorial_77/03_basics.html, for example, then you will see in the "Continuation" section that the "2" character is used a continuation character and the continued line begins on position 6 (column 6)

      SUBROUTINE CONSTRUCT_HEATS(lM, Ag, Yg) ! JP 2021_04_14: I think this subroutine computes local element tangent matrix and residual vector FOR ALL ELEMENTS (using Gauss quadrature) (in the subroutine HEATS3D or HEATS2D) and then automatically assembles those element contributions into the global tangent matrix and residual vector (in the subroutine DOASSEM from LHSA.f)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo) ! JP: from MOD.f, tDof is the " Total number of degrees of freedom per node"; tnNo is the "Total number of nodes"

      INTEGER(KIND=IKIND) a, e, g, Ac, eNoN, cPhys
      REAL(KIND=RKIND) w, Jac, ksix(nsd,nsd)

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), al(:,:), yl(:,:),
     2   N(:), Nx(:,:), lR(:,:), lK(:,:,:)

      eNoN = lM%eNoN

!     HEATS: dof = 1
      ALLOCATE(ptr(eNoN), xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN),
     2   N(eNoN), Nx(nsd,eNoN), lR(dof,eNoN), lK(dof*dof,eNoN,eNoN))

!     Loop over all elements of mesh
      DO e=1, lM%nEl ! JP: 3/25/21: e is element number or element ID?; nEl is number of elements?; lM is location array or location matrix (from pinsky notes)???
!        Update domain and proceed if domain phys and eqn phys match
         cDmn  = DOMAIN(lM, cEq, e)
         cPhys = eq(cEq)%dmn(cDmn)%phys
         IF (cPhys .NE. phys_heatS) CYCLE

!        Update shape functions for NURBS
         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

!        Create local copies
         DO a=1, eNoN ! JP: i think eNoN = nen = number of local nodes on each element (i think this because here, we are looping a over 1 to eNoN and a the first argument to the IEN array below and recall from ME335A that IEN(a,e) returns the global number number for local node a of element e)
            Ac = lM%IEN(a,e) ! JP: i think Ac is the global node number that local node a of element e corresponds to
            ptr(a)  = Ac
            xl(:,a) = x(:,Ac)
            al(:,a) = Ag(:,Ac) ! JP 2021_04_02: al = Ag = (1 - alpha_m) * Ao + alpha_m * An (from PIC.f; the exact line is: Ag(s:e,a) = Ao(s:e,a)*coef(1) + An(s:e,a)*coef(2) )
                ! JP 2021_04_02: this 'al' gets used in the below subroutine, HEATS2D, in the line "Td = Td + N(a)*al(i,a)"
                ! N gets allocated in the line "N = lM%N(:,g)" in this file
                ! JP 2021_04_14: I think "al" stands for "acceleration local"
            yl(:,a) = Yg(:,Ac) ! JP 2021_04_14: I think "yl" stands for "velocity local"; recall from PIC.f that Yg is the velocity at time n + alpha_f
         END DO

!        Gauss integration
         lR = 0._RKIND ! JP 2021_04_14: I think this line initializes the local element residual vector to be zeros
         lK = 0._RKIND ! JP 2021_04_14: I think this line initializes the local element tangent matrix to be zeros
         DO g=1, lM%nG ! JP 2021_04_02: is nG = number of gauss points? ! JP 2021_04_14: yes, from MOD.f under the section "TYPE mshType", "nG" is the "Number of Gauss points for integration"
            IF (g.EQ.1 .OR. .NOT.lM%lShpF) THEN
               CALL GNN(eNoN, nsd, lM%Nx(:,:,g), xl, Nx, Jac, ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
            END IF
            w = lM%w(g) * Jac ! JP 2021_04_14: from NN.f in the subroutine, GNN, I think "Jac" is the jacobian of the parent element (recall from ME335A that the jacobian of the parent element is = the area of the parent element in 2D or something like that) ! JP 2021_04_14: lM%w(g) is the "Gauss weights" (from MOD.f under section "TYPE mshType")
            N = lM%N(:,g) ! JP 2021_04_02: lM = msh(iM) b/c lM is the first argument of CONSTRUCT_HEATS and CONSTRUCT_HEATS is called in EQASSEM.f, where lM is the first argument of GLOBALEQASSEM which is called in MAIN.f (in the line "CALL GLOBALEQASSEM(msh(iM), Ag, Yg, Dg)"). and msh is defined in MOD.f (in the line "TYPE(mshType), ALLOCATABLE :: msh(:)" under the section "DERIVED TYPE VARIABLES") so that means that msh is of type mshType. And one of the "fields" or "properties" of mshType is "N", which is supposed to be the "Parent shape function". I think lM%N is of size (eNoN, nG), where eNoN is the number nodes per element and nG is the number of Gauss points (for integration) but I am just speculating here and I dont know for sure (I think I am right now) since later in this code, N (not lM%N) is indexed using the element local node number

            IF (nsd .EQ. 3) THEN
               CALL HEATS3D(eNoN, w, N, Nx, al, yl, lR, lK)

            ELSE IF (nsd .EQ. 2) THEN
               CALL HEATS2D(eNoN, w, N, Nx, al, yl, lR, lK)

            END IF
         END DO ! g: loop

!        Assembly
#ifdef WITH_TRILINOS
         IF (eq(cEq)%assmTLS) THEN
            CALL TRILINOS_DOASSEM(eNoN, ptr, lK, lR)
         ELSE
#endif
            CALL DOASSEM(eNoN, ptr, lK, lR) ! JP 2021_04_14: I think line assembles element e's local tangent matrix and residual vector into the global tangent matrix and residual vector
#ifdef WITH_TRILINOS
         END IF
#endif
      END DO ! e: loop

      DEALLOCATE(ptr, xl, al, yl, N, Nx, lR, lK)

      RETURN
      END SUBROUTINE CONSTRUCT_HEATS
!####################################################################
!     This is for solving the heat equation in a solid
      PURE SUBROUTINE HEATS3D (eNoN, w, N, Nx, al, yl, lR, lK) ! JP 2021_04_08: weiguang said that lR = local residual vector; lK = local tangent matrix (Jacobian) (where the residual vector and the tangent matrix are used in the newton solver in the generalized alpha method) ! JP 2021_04_14: this subroutine is the same as HEATS2D, except this subroutine is for 3D problems instead of 2D problems.
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(nsd,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(1,eNoN), lK(1,eNoN,eNoN)

      INTEGER(KIND=IKIND) i, a, b
      REAL(KIND=RKIND) nu, T1, amd, wl, Td, Tx(nsd), s, rho

      nu  = eq(cEq)%dmn(cDmn)%prop(conductivity)
      s   = eq(cEq)%dmn(cDmn)%prop(source_term)
      rho = eq(cEq)%dmn(cDmn)%prop(solid_density)
            ! JP: why do we need density here? In ME335A course reader, there is no density term involved in the heat equation. maybe there is "another" heat equation that is more general than the one in pinsky's notes? see references below:
                    ! https://www.google.com/search?q=heat+equatin&sxsrf=ALeKk03wL9g4J1LfoXJRIzGWt8ZpX4K8dA%3A1616898296176&ei=-OhfYNGWCo7K0PEPzPON4AM&oq=heat+equatin&gs_lcp=Cgdnd3Mtd2l6EAMyBwgAELEDEAoyBAgAEAoyBAgAEAoyBAgAEEMyBAgAEAoyBAgAEAoyBAgAEAoyBAgAEAoyBAgAEAoyBAgAEAo6BwgAEEcQsAM6CAgAELEDEJECOgoIABCHAhCxAxAUOgUIABCxAzoCCAA6BQgAEJECUPYqWMcvYK0yaAFwAngAgAH8AYgB1AeSAQU0LjMuMZgBAKABAaoBB2d3cy13aXrIAQjAAQE&sclient=gws-wiz&ved=0ahUKEwjR8tLE99HvAhUOJTQIHcx5AzwQ4dUDCA0&uact=5
                    ! https://en.wikipedia.org/wiki/Thermal_conduction
                    ! https://en.wikipedia.org/wiki/Heat_equation
                    ! https://tutorial.math.lamar.edu/classes/de/theheatequation.aspx

      T1  = eq(cEq)%af*eq(cEq)%gam*dt
      amd = eq(cEq)%am * rho/T1
      i   = eq(cEq)%s

      wl = w*T1

      Td = -s
      Tx = 0._RKIND
      DO a=1, eNoN
         Td = Td + N(a)*al(i,a)

         Tx(1) = Tx(1) + Nx(1,a)*yl(i,a)
         Tx(2) = Tx(2) + Nx(2,a)*yl(i,a)
         Tx(3) = Tx(3) + Nx(3,a)*yl(i,a)
      END DO
      Td = Td * rho

      DO a=1,eNoN
         lR(1,a) = lR(1,a) + w*(N(a)*Td
     2      + (Nx(1,a)*Tx(1) + Nx(2,a)*Tx(2) + Nx(3,a)*Tx(3))*nu)

         DO b=1,eNoN
            lK(1,a,b) = lK(1,a,b) + wl*(N(a)*N(b)*amd
     2         + nu*(Nx(1,a)*Nx(1,b) +Nx(2,a)*Nx(2,b) +Nx(3,a)*Nx(3,b)))
         END DO
      END DO

      RETURN
      END SUBROUTINE HEATS3D
!--------------------------------------------------------------------
!     This is for solving the heat equation in a solid

! JP: i think this subroutines constructs and computes and returns the local element stiffness matrix (lK) and the local element force vector (lR) (which includes contributions from the source term; not sure if the dirichlet BC is accounted for here; the neumann BC is accounted for in BHEATS); The reason that I think this code computes and returns the local element stiffness matrix (lK) and the local element force vector (lR) is because below, it says "REAL(KIND=RKIND), INTENT(INOUT) :: lR(1,eNoN), lK(1,eNoN,eNoN)", which means that this subroutines can both read in the values of lR and lK and modify them (return them) as well (because INTENT is INOUT) ! reference: https://fortran-lang.org/learn/quickstart/organising_code
    ! JP 2021_04_02: wait is lR the element force vector or is it the element residual vector?? (maybe the 'R' stands for "residual"??)

! JP 2021_04_08: I think mahdi might have written the svFSI code, so maybe refer to his thesis to see which technique/form/scheme he used to implement the generalized alpha method. and for more details on the weak forms, etc for the fluid Navier-stokes equations, etc

      PURE SUBROUTINE HEATS2D (eNoN, w, N, Nx, al, yl, lR, lK) ! JP 2021_04_14: I think this subroutine computes for a given element and a given Gauss point, its local tangent matrix (lK) and residual vector (lR) where w, N, and Nx are the weight used in Gauss integration, shape function, and first derivative of the shape function, respectively (all already evaluated at the given Gauss point). Note that this subroutine does not account for the contributions of the EBC and the NBC (dirichlet and neumann BCs) to the force vector, which is used to compute the residual vector, lR; it only accounts for the contribution of the source term to the force vector
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(nsd,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(1,eNoN), lK(1,eNoN,eNoN)

      INTEGER(KIND=IKIND) i, a, b
      REAL(KIND=RKIND) nu, T1, amd, wl, Td, Tx(nsd), s, rho

      nu  = eq(cEq)%dmn(cDmn)%prop(conductivity) ! JP 2021_04_14: I think is nu is constant here (not a function of space) ! JP 2021_04_14: also I think we are assuming a thermally isotropic material here, so the thermal "conductivity" tensor is a diagonal matrix (eqn 6.7 in ME335A Ch6 Couse reader)
      s   = eq(cEq)%dmn(cDmn)%prop(source_term) ! JP 2021_04_14: I think the source term, s, is constant here (not a function of space)
      rho = eq(cEq)%dmn(cDmn)%prop(solid_density)
      ! JP 2021_04_14: note that in ME335A Ch6 of the course reader, there is a heat capacity term in the strong form (and thus the mass matrix). This HEATS.f code does not explicitly use a heat capacity so i think that the heat capacity is built-in directly into the value of "rho"

      T1  = eq(cEq)%af*eq(cEq)%gam*dt
      amd = eq(cEq)%am * rho/T1
      i   = eq(cEq)%s

      wl = w*T1

      Td = -s
      Tx = 0._RKIND
      DO a=1, eNoN ! JP: is eNoN = nen = number of element nodes??
         Td = Td + N(a)*al(i,a) ! JP: note that the "al" used here is different from the "al" defined in MOD.f (which is "Min norm of face normals in contact"). The "al" used here is defined in this file in the section of "ALLOCATE(ptr(eNoN), xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN)"
         ! JP 2021_04_02: al is set in the line "al(:,a) = Ag(:,Ac)", with al = Ag = (1 - alpha_m) * Ao + alpha_m * An (from PIC.f; the exact line is: Ag(s:e,a) = Ao(s:e,a)*coef(1) + An(s:e,a)*coef(2) )
         ! JP 2021_04_02: I think al acts like the time derivative of "g" (g_dot) in the galerkin form for the unsteady heat equation in ME335B (see my GoodNotes notes on the unsteady heat eqn and FEM). However, al is not necessarily g_dot because we are using the generalized alpha method here. Thus, I think in order to better understand what this line is doing, I need to try to using the generalized alpha method to solve unsteady heat equation with FEM
                ! JP 2021_04_14: this line is computing the sum of 1) the product of the mass matrix and the acceleration vector (well, more specifically the product of the shape function and the acceleration vector) and 2) the source term (where this sum is used in the residual vector)

         Tx(1) = Tx(1) + Nx(1,a)*yl(i,a) ! JP 2021_04_14: I think this line computes the product of the multiplication between the stiffness matrix and the velocity (well more specifically, the product of the multiplication between derivative of the shape function and the velocity) (which is used in the residual)
         Tx(2) = Tx(2) + Nx(2,a)*yl(i,a) ! JP 2021_04_14: I think this line computes the product of the multiplication between the stiffness matrix and the velocity (well more specifically, the product of the multiplication between derivative of the shape function and the velocity) (which is used in the residual)
      END DO
      Td = Td * rho ! JP 2021_04_14: Note that "Td" contains the source term. So, why are we multiplying the source term "s" by the density here? Maybe this is a bug or maybe mahdi or whoever wrote this code is solving a slightly different heat eqn from the one that I learned in ME335B

      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*(N(a)*Td + (Nx(1,a)*Tx(1)
     2      + Nx(2,a)*Tx(2))*nu) ! JP: I think this line might be computing the RHS forcing term, since Td looks like it might be the source term here
                    ! JP 2021_04_14: in the above calculation of lR(1,a), "N(a)*Td" represents the -- mass matrix multiplied by the acceleration vector plus the source term (the force vector with contribution due to only the source term); the term "(Nx(1,a)*Tx(1) + Nx(2,a)*Tx(2))*nu" represents the stiffness matrix multiplied by the velocity vector
                    ! JP 2021_04_14: I finally understand the computation of lR(1,a) here. The code here "lR(1,a) = lR(1,a) + w*(N(a)*Td + (Nx(1,a)*Tx(1)" matches my derivation of what the residual vector for the unsteady heat eqn looks like (see my notes titled "Unsteady Heat Eqn w/ generalized-alpha (using an accleration update/correction in the Newton solver)" in the "FEM" notebook in GoodNotes)

         DO b=1, eNoN
            lK(1,a,b) = lK(1,a,b) + wl*(N(a)*N(b)*amd
     2         + nu*(Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b))) ! JP 2021_03_27: i think this line might be computing the stiffness matrix? maybe lK is the local element stiffness matrix, Ke? And maybe Ke(a, b) = lK(1, a, b)? Idk what the first dimension of lK corresponds to though (what that "1" corresponds to). And it looks like maybe Nx is the element shape function, where the first parameter corresponds to the dimension?? -- compare these thoughts against what is in Pinsky's ME335A notes (recall that vijay mentioned that he uses pinsky's continuum mechanics (maybe he meant FEM notes) for voigh notation in elasticity at least, but maybe he uses it for heat eqn too?)
                    ! JP 2021_03_27: i think the the "l" in front of "lK" stands for "local", so maybe here "lK" is the local stiffness matrix (ie the local stiffness matrix for element e?)
                    ! JP 2021_03_27: as such, maybe "lR" is the "local" force vector (ie the local RHS) for element e, where lR should include all contributions from the source term, the neumann BC, etc
                    ! JP 2021_04_02: actually, from MOD.f Nx is the derivative (spatial) of the element shape function and I think N is the element shape function
                    ! JP 2021_04_14: in the above computation of lK(1,a,b), the term "N(a)*N(b)*amd" represents the mass matrix (m^e(a, b)); and the term "nu*(Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b))" represents the stiffness matrix
                    ! JP 2021_04_14: I finally understand the computation of lK(1,a,b) here. The code here "lK(1,a,b) = lK(1,a,b) + wl*(N(a)*N(b)*amd ..." matches my derivation of what the tangent matrix for the unsteady heat eqn looks like (see my notes titled "Unsteady Heat Eqn w/ generalized-alpha (using an accleration update/correction in the Newton solver)" in the "FEM" notebook in GoodNotes)
         END DO
      END DO

      RETURN
      END SUBROUTINE HEATS2D
!--------------------------------------------------------------------
      PURE SUBROUTINE BHEATS (eNoN, w, N, h, lR)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), h
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN)

      INTEGER(KIND=IKIND) a

      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*N(a)*h ! JP: is this computing the neumann BC here?? yes, i think so because BHEATS is called from EQASSEM.f under the section with comment "Construct Neumann BCs"
                ! JP 2021_06_03: issue: I think there is an error here. I think the "h" (the Neumann BC value) should be multiplied by a negative sign (or it should be subtracted from lR instead of added to lR), b/c according to my GoodNotes for the Transient Heat Equation, we have to bring the (known) Neumann BC over to the LHS which requires us to subtract it.
                    ! JP 2021_06_03: other the the above issue, this line "lR(1,a) = lR(1,a) + w*N(a)*h" makes sense to me and agrees with my FEM GoodNotes notes for the Transient Heat Equation
      END DO

      RETURN
      END SUBROUTINE BHEATS
!####################################################################

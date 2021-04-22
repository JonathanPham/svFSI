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
!     This routine contains predictor/initiator/corrector routines, as
!     part of time integration scheme. ! JP 2021_04_02: PIC (the file name) stands for Predictor/Initiator/Corrector
!
!--------------------------------------------------------------------

!     This is the predictor
      SUBROUTINE PICP
      USE COMMOD
      IMPLICIT NONE

      INTEGER(KIND=IKIND) iEq, s, e
      REAL(KIND=RKIND) coef

!     Prestress initialization
      IF (pstEq) THEN
         pS0 = pS0 + pSn
         Ao = 0._RKIND
         Yo = 0._RKIND
         Do = 0._RKIND
      END IF

!     IB treatment: Set dirichlet BC and update traces. For explicit
!     coupling, compute FSI forcing and freeze it for the time step.
!     For implicit coupling, project IB displacement on background
!     mesh and predict quantities at next time step
      IF (ibFlag) THEN
!        Set IB Dirichlet BCs
         CALL IB_SETBCDIR(ib%Yb, ib%Ubo)

!        Update IB location and tracers
         CALL IB_UPDATE(Do)

         IF (ib%cpld .EQ. ibCpld_E) THEN
!           FSI forcing for immersed bodies (explicit coupling)
            CALL IB_CALCFFSI(Ao, Yo, Do, ib%Auo, ib%Ubo)

         ELSE IF (ib%cpld .EQ. ibCpld_I) THEN
!           Project IB displacement (Ubn) to background mesh
            CALL IB_PRJCTU(Do)

!           Predictor step for implicit coupling
            CALL IB_PICP()

         END IF
      END IF

      DO iEq=1, nEq
         s = eq(iEq)%s
         e = eq(iEq)%e
         coef = (eq(iEq)%gam - 1._RKIND)/eq(iEq)%gam
         An(s:e,:) = Ao(s:e,:)*coef ! JP 2021_04_14: I think this the same predictor used for the acceleration in Bazilevs 2007 eqn 87 ("Variational multiscale residual-based turbulence modeling for large eddy simulation of incompressible flows", CMAME)

!        electrophysiology
         IF (eq(iEq)%phys .EQ. phys_CEP) THEN
            CALL CEPINTEG(iEq, e, Do)
         END IF

         Yn(s:e,:) = Yo(s:e,:) ! JP 2021_04_14: I think this the same predictor used for the velocity in Bazilevs 2007 eqn 86 ("Variational multiscale residual-based turbulence modeling for large eddy simulation of incompressible flows", CMAME)

         IF (dFlag) THEN ! JP 2021_04_14: I think dFlag is set in INITIALIZE.f. (e.g. on line "dFlag = .TRUE." under the section "CASE (phys_lElas)"), but dFlag has a default value of ".FALSE." heatF and heatS in that file. I think based on this section, "dFlag = .TRUE." means that we want to provide some kind of special initialization for the displacement term (or maybe it means that the PDE of interest has a 2nd order derivative in time, such as the eqns for elastodynamics, so there, we would need some kind of initial guess for the displacement, similar to what is done in the Newmark method (Lecture 7, ME335B)); I also wrote a note about dFlag in MAIN.f (see that note for more details)
            IF (.NOT.sstEq) THEN
!              struct, lElas, FSI (struct, mesh)
               coef = dt*dt*(0.5_RKIND*eq(iEq)%gam - eq(iEq)%beta)
     2            /(eq(iEq)%gam - 1._RKIND)
               Dn(s:e,:) = Do(s:e,:) + Yn(s:e,:)*dt + An(s:e,:)*coef
            ELSE
!              ustruct, FSI
               IF (eq(iEq)%phys .EQ. phys_ustruct .OR.
     2             eq(iEq)%phys .EQ. phys_FSI) THEN
                  coef = (eq(iEq)%gam - 1._RKIND)/eq(iEq)%gam
                  Ad(:,:)   = Ad(:,:)*coef
                  Dn(s:e,:) = Do(s:e,:)
               ELSE IF (eq(iEq)%phys .EQ. phys_mesh) THEN
!              mesh
                  coef = dt*dt*(0.5_RKIND*eq(iEq)%gam - eq(iEq)%beta)
     2               /(eq(iEq)%gam - 1._RKIND)
                  Dn(s:e,:) = Do(s:e,:) + Yn(s:e,:)*dt + An(s:e,:)*coef
               END IF
            END IF
         ELSE
            Dn(s:e,:) = Do(s:e,:)
         END IF
      END DO

      RETURN
      END SUBROUTINE PICP
!====================================================================
!     This is the initiator
! JP 2021_04_02: what the heck is an initiator step?? does the initiator step come before or after the predictor step? from the MAIN.f file, it looks like the order of calls is: PICP, PICI, and then PICC (predictor, initiator and then corrector). So maybe the initiator step computes the initial guess for the Newton solver??
! JP 2021_04_14: according to MAIN.f, in the location where PICI is called, it says "(quantities at n+am, n+af)", so i think that maybe PICI is computing what initial guesses of what the soltns at the n + alpha_m and n + alpha_f steps are, e.g. eqns 89 and 90 in Bazilevs 2007 ! JP 2021_04_14: after looking through this subroutine in more detail, this is indeed what the code is doing. It is indeed computing the A(n + alpha_m), Y(n + alpha_f)
      SUBROUTINE PICI(Ag, Yg, Dg)
      USE COMMOD
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(INOUT) :: Ag(tDof,tnNo), Yg(tDof,tnNo),
     2   Dg(tDof,tnNo) ! JP: from MOD.f, tDof is the " Total number of degrees of freedom per node"; tnNo is the "Total number of nodes"

      INTEGER(KIND=IKIND) s, e, i, a
      REAL(KIND=RKIND) coef(4)

      dof         = eq(cEq)%dof ! JP: cEq = "Current equation" in MOD.f; eq is a 1d array of type eqType that stores all data related to equations (MOD.f)
      eq(cEq)%itr = eq(cEq)%itr + 1

      DO i=1, nEq
         s       = eq(i)%s ! JP: from MOD.f in the eqType type, i think 's' stands for 'start'
         e       = eq(i)%e ! JP: from MOD.f in the eqType type, i think 'e' stands for 'end'
         coef(1) = 1._RKIND - eq(i)%am ! JP: from MOD.f, I think 'am' is one of the generalized alpha time integrator parameters
         coef(2) = eq(i)%am
         coef(3) = 1._RKIND - eq(i)%af ! JP: from MOD.f, I think 'af' is one of the generalized alpha time integrator parameters
         coef(4) = eq(i)%af

         ! JP: I think this below for / DO loop is computing some generalized alpha method variables / parameters? however, observe below that there are quantities for acceleration, velocity, and displacement. As such, I think this method is not using the Jansen 2000 generalized alpha method, but i think it is using the original generalized alpha method from the paper, '''J. Chung and G. M. Hulbert, “A time integration algorithm for structural dynamics with improved numerical dissipation: The generalized-α method”, Journal of Applied Mechanics, 60 (1993) 371–75.'''.
         DO a=1, tnNo ! JP: from MOD.f, tDof is the " Total number of degrees of freedom per node"; tnNo is the "Total number of nodes"
            Ag(s:e,a) = Ao(s:e,a)*coef(1) + An(s:e,a)*coef(2) ! JP: from MOD.f in the eqType type, I think 'A' stands for acceleration; this eqn comes from eqn 12 in the Chung 1993 paper. ! JP 2021_04_14: this line matches eqn 89 in Bazilevs 2007
                        ! however, i have an important question: from this eqn, coef(1) is = 1 - alpha_m (am) and coef(2) is = alpha_m (am) and I thought that Ao is the old acceleration and An is the new acceleration, so why is Ao being multipled by 1 - alpha_m here? Shouldnt it be multipled by alpha_m instead, and shouldnt An be multiplied by An (at least, this is how it is done in the Chung 1993 paper).
                            ! from my notes in INITIALIZE.f, which is where the am and af terms are initialized (from roInf), note that there is a difference between the generalized alpha methods between the Jansen 2000 and Chung 1993 papers. Maybe this difference is why the eqn listed here in PIC.f for Ag and Yg and Dg dont match what is shown in the Chung 1993 paper. So my question is: which method (Chung or Jansen) is correct? Maybe see Ju and Ingrid's paper on generalized alpha for idea of which is the correct version?
                            ! JP 2021_04_02: i checked and yes, the Chung 1993 and Jansen 2000 paper have different defitions of alpha, so thats why their definitions of alpha_m and alpha_f differ. I think this svFSI code is using the Jansen 2000 definition of alpha_m and alpha_f.
            Yg(s:e,a) = Yo(s:e,a)*coef(3) + Yn(s:e,a)*coef(4) ! JP: from MOD.f in the eqType type, I think 'Y' stands for velocity; this eqn comes from eqn 11 in the Chung 1993 paper. ! JP 2021_04_14: this line matches eqn 90 in Bazilevs 2007
            Dg(s:e,a) = Do(s:e,a)*coef(3) + Dn(s:e,a)*coef(4) ! JP: from MOD.f in the eqType type, I think 'D' stands for displacement; this eqn comes from eqn 10 in the Chung 1993 paper.
         END DO
      END DO

      IF (pstEq) THEN
         pSn(:,:) = 0._RKIND
         pSa(:)   = 0._RKIND
      END IF

      RETURN
      END SUBROUTINE PICI
!====================================================================
!     This is the corrector. Decision for next eqn is also made here
      SUBROUTINE PICC
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      LOGICAL :: l1, l2, l3, l4
      INTEGER(KIND=IKIND) :: s, e, a, Ac
      REAL(KIND=RKIND) :: coef(4), r1, dUl(nsd)

      s       = eq(cEq)%s
      e       = eq(cEq)%e
      coef(1) = eq(cEq)%gam*dt
      coef(2) = eq(cEq)%beta*dt*dt
      coef(3) = 1._RKIND / eq(cEq)%am
      coef(4) = eq(cEq)%af*coef(1)*coef(3)

      IF (sstEq) THEN ! JP 2021_04_14: sstEq is a boolean than says "Whether velocity-pressure based structural dynamics solver is used" according to MOD.f; this variable is set in READFILES.f
!        ustruct, FSI (ustruct)
         IF (eq(cEq)%phys .EQ. phys_ustruct .OR.
     2       eq(cEq)%phys .EQ. phys_FSI) THEN
            DO a=1, tnNo
               An(s:e,a)   = An(s:e,a)   - R(:,a)
               Yn(s:e,a)   = Yn(s:e,a)   - R(:,a)*coef(1)
               dUl(:)      = Rd(:,a)*coef(3) + R(1:dof-1,a)*coef(4)
               Ad(:,a)     = Ad(:,a)     - dUl(:)
               Dn(s:e-1,a) = Dn(s:e-1,a) - dUl(:)*coef(1)
            END DO
         ELSE IF (eq(cEq)%phys .EQ. phys_mesh) THEN
            DO a=1, tnNo
               An(s:e,a)   = An(s:e,a) - R(:,a)
               Yn(s:e,a)   = Yn(s:e,a) - R(:,a)*coef(1)
               Dn(s:e,a)   = Dn(s:e,a) - R(:,a)*coef(2)
            END DO
         END IF
      ELSE
         DO a=1, tnNo
            ! JP 2021_04_14: '-R(:,a)' represents the update/correction for the acceleration, where this correction/update was obtained in the newton solver by solving "K*delta = -R" (or in other words, "K*delta = R" where after we solve for delta, we must apply the minus sign in the correction step in gen alpha), where K represents the tangent matrix and R represents the residual and delta represents the update/correction in the Newton solver in gen alpha ! actually, is R(:,a) the acceleration update (the 'delta' in the above eqn) or is it the residual?? for now, I will assume that it is the 'delta' term, since it does not say anywhere in this entire svFSI if R is the residual or the delta... I think from SOLVE.f, in the call to NSSOLVER or GMRESS, R is the residual initially but then after NSSOLVER or GMRESS solve the linear system with that residual R, as a last step they overwrite R with the delta update / correction instead (maybe this is done to save memory??), but I am not sure if this is actually the case or not... I need to ask someone, maybe weiguang, for confirmation

            An(s:e,a) = An(s:e,a) - R(:,a) ! JP 2021_04_14: this line matches eqn 94 in Bazilevs 2007, assuming that "R(:,a)" is the update/correction, as mentioned in the above comment...

            Yn(s:e,a) = Yn(s:e,a) - R(:,a)*coef(1) ! JP 2021_04_14: this line matches eqn 95 in Bazilevs 2007, assuming that "R(:,a)" is the update/correction, as mentioned in the above comment...

            Dn(s:e,a) = Dn(s:e,a) - R(:,a)*coef(2)
         END DO
      END IF

      IF ((eq(cEq)%phys .EQ. phys_ustruct) .OR.
     2    (eq(cEq)%phys .EQ. phys_stokes)) THEN
         CALL PICETH()
      END IF

      IF (eq(cEq)%phys .EQ. phys_FSI) THEN
         s = eq(2)%s
         e = eq(2)%e
         DO Ac=1, tnNo
            IF (ISDOMAIN(cEq, Ac, phys_struct) .OR.
     2          ISDOMAIN(cEq, Ac, phys_ustruct) .OR.
     3          ISDOMAIN(cEq, Ac, phys_lElas)) THEN
               An(s:e,Ac) = An(1:nsd,Ac)
               Yn(s:e,Ac) = Yn(1:nsd,Ac)
               Dn(s:e,Ac) = Dn(1:nsd,Ac)
            END IF
         END DO
      END IF

!     Update Xion for cardiac electrophysiology
      IF (eq(cEq)%phys .EQ. phys_CEP) THEN
         s = eq(cEq)%s
         DO a=1, tnNo
            Xion(1,a) = Yn(s,a)
         END DO
      END IF

!     Update prestress at the nodes and re-initialize
      IF (pstEq) THEN
         CALL COMMU(pSn)
         CALL COMMU(pSa)
         DO a=1, tnNo
            IF (.NOT.ISZERO(pSa(a))) THEN
               pSn(:,a) = pSn(:,a) / pSa(a)
            END IF
         END DO
         pSa = 0._RKIND
      END IF

!     Filter out the non-wall displacements for CMM equation
      IF (eq(cEq)%phys.EQ.phys_CMM .AND. .NOT.cmmInit) THEN
         DO a=1, tnNo
            r1 = REAL(cmmBdry(a), KIND=RKIND)
            Dn(s:e-1,a) = Dn(s:e-1,a)*r1
         END DO
      END IF

!     IB treatment
      IF (ibFlag) CALL IB_PICC()

      IF (ISZERO(eq(cEq)%FSILS%RI%iNorm)) eq(cEq)%FSILS%RI%iNorm = eps
      IF (ISZERO(eq(cEq)%iNorm)) eq(cEq)%iNorm = eq(cEq)%FSILS%RI%iNorm
      IF (eq(cEq)%itr .EQ. 1) THEN
         eq(cEq)%pNorm = eq(cEq)%FSILS%RI%iNorm/eq(cEq)%iNorm
      END IF
      r1 = eq(cEq)%FSILS%RI%iNorm/eq(cEq)%iNorm

      l1 = eq(cEq)%itr .GE. eq(cEq)%maxItr
      l2 = r1 .LE. eq(cEq)%tol
      l3 = r1 .LE. eq(cEq)%tol*eq(cEq)%pNorm
      l4 = eq(cEq)%itr .GE. eq(cEq)%minItr
      IF (l1 .OR. ((l2.OR.l3).AND.l4)) eq(cEq)%ok = .TRUE.
      IF (ALL(eq%ok)) RETURN

      IF (eq(cEq)%coupled) THEN
         cEq = cEq + 1
         IF (ALL(.NOT.eq%coupled .OR. eq%ok)) THEN
            DO WHILE (cEq .LE. nEq)
               IF (.NOT.eq(cEq)%coupled) EXIT
               cEq = cEq + 1
            END DO
         ELSE
            IF (cEq .GT. nEq) cEq = 1
            DO WHILE (.NOT.eq(cEq)%coupled)
               cEq = cEq + 1
               IF (cEq .GT. nEq) cEq = 1
            END DO
         END IF
      ELSE
         IF (eq(cEq)%ok) cEq = cEq + 1
      END IF

      RETURN
      END SUBROUTINE PICC
!====================================================================
!     Pressure correction at edge nodes for Taylor-Hood type element
!     via interpolation
      SUBROUTINE PICETH()
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      LOGICAL THflag
      INTEGER(KIND=IKIND) a, b, e, g, s, iM, Ac, eType, eNoN, eNoNq
      REAL(KIND=RKIND) Jac, eVol, p, xp(nsd), xi0(nsd), xi(nsd),
     2   ksix(nsd,nsd)

      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), xql(:,:), pl(:), Nq(:),
     2   Nqx(:,:), sA(:), sF(:)

      THflag = .FALSE.
      DO iM=1, nMsh
         IF (msh(iM)%nFs .EQ. 2) THEN
            THflag = .TRUE.
            EXIT
         END IF
      END DO
      IF (.NOT.THflag) RETURN

      ALLOCATE(sA(tnNo), sF(tnNo))
      sF(:) = 0._RKIND
      sA(:) = 0._RKIND

      s = eq(cEq)%s
      DO iM=1, nMsh
         IF (msh(iM)%nFs .EQ. 1) CYCLE

         eType = msh(iM)%fs(2)%eType

         eNoN  = msh(iM)%fs(1)%eNoN
         eNoNq = msh(iM)%fs(2)%eNoN
         ALLOCATE(xl(nsd,eNoN), xql(nsd,eNoNq), pl(eNoNq), Nq(eNoNq),
     2      Nqx(nsd,eNoNq))

         xi0 = 0._RKIND
         DO g=1, msh(iM)%fs(2)%nG
            xi0 = xi0 + msh(iM)%fs(2)%xi(:,g)
         END DO
         xi0 = xi0 / REAL(msh(iM)%fs(2)%nG, KIND=RKIND)

         DO e=1, msh(iM)%nEl
            cDmn = DOMAIN(msh(iM), cEq, e)
            IF ((eq(cEq)%dmn(cDmn)%phys .NE. phys_ustruct) .AND.
     2          (eq(cEq)%dmn(cDmn)%phys .NE. phys_stokes)) CYCLE

            DO a=1, eNoN
               Ac = msh(iM)%IEN(a,e)
               xl(:,a) = x(:,Ac)
            END DO

            DO a=1, eNoNq
               Ac = msh(iM)%IEN(a,e)
               pl(a)    = Yn(s+nsd,Ac)
               xql(:,a) = xl(:,a)
            END DO

            eVol = 0._RKIND
            DO g=1, msh(iM)%fs(2)%nG
               IF (g.EQ.1 .OR. .NOT.msh(iM)%fs(2)%lShpF) THEN
                  CALL GNN(eNoNq, nsd, msh(iM)%fs(2)%Nx(:,:,g), xql,
     2               Nqx, Jac, ksix)
                  IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
               END IF
               eVol = eVol + msh(iM)%fs(2)%w(g)*Jac
            END DO

            DO a=eNoNq+1, eNoN
               Ac = msh(iM)%IEN(a,e)
               xp = xl(:,a)

               xi = xi0
               CALL GETNNX(eType, eNoNq, xql, msh(iM)%fs(2)%xib,
     2            msh(iM)%fs(2)%Nb, xp, xi, Nq, Nqx)

               p = 0._RKIND
               DO b=1, eNoNq
                  p = p + pl(b)*Nq(b)
               END DO

               sF(Ac) = sF(Ac) + p*eVol
               sA(Ac) = sA(Ac) + eVol
            END DO

         END DO ! e-loop
         DEALLOCATE(xl, xql, pl, Nq, Nqx)
      END DO ! iM-loop

      CALL COMMU(sA)
      CALL COMMU(sF)

      DO a=1, tnNo
         IF (.NOT.ISZERO(sA(a))) Yn(s+nsd,a) = sF(a)/sA(a)
      END DO

      DEALLOCATE(sA, sF)

      RETURN
      END SUBROUTINE PICETH
!====================================================================

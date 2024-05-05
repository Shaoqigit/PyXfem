import numpy as np


class RBSolver_fromResidual:

    def __init__(self, maxmodes, operatorModesSpace, operatorModesParam, rhsModesSpace, rhsModesParam,
                 parameterTrainingSpace, dtype_=float):

        self.maxModes = maxmodes

        # Affine decomposition : op(x,mu) = sum_i op_i(x) * theta_i(mu)
        self.operatorModesSpace = operatorModesSpace  # op_i
        self.operatorModesParam = operatorModesParam  # theta_i

        # Affine decomposition : rhs(x,mu) = sum_i rhs_i(x) * theta_i(mu)
        self.rhsModesSpace = rhsModesSpace  # rhs_i
        self.rhsModesParam = rhsModesParam  # theta_i

        self.parameterTrainingSpace = parameterTrainingSpace

        self.RB_Basis = np.zeros((self.operatorModesSpace[0].shape[0], self.maxModes), dtype=dtype_)
        self.RB_Coeffs = np.zeros((len(self.parameterTrainingSpace), self.maxModes), dtype=dtype_)
        self.nbModes = 0
        self.worstParamHistory = []

        self.TildeComputedFor = -1  # The number of mode for which Tilde quantities have been precomputed
        self.Atilde_k = []  # Precomputed projected affine decomposition (operator space part)
        self.Btilde_k = []  # Precomputed projected affine decomposition (rhs space part)

        self.dtype = dtype_

    # From the current state of the greedy strategy :
    # Computes:
    # * the worst snapshot and associated residual
    # * the reduced problem for this particular snapshot

    def computeWorstParameterAndAssociatedData(self, basis, coeffs):

        maxResidNorm = -1.
        imuTilde = 0

        # Loop over mu_i
        for imu, mui in enumerate(self.parameterTrainingSpace):
            # mui = self.parameterSnapshots[imu]

            # Residual_mui = np.zeros(Bi[0].shape)

            Bcurr = np.zeros(self.rhsModesSpace[0].shape, dtype=self.dtype)
            Acurr = np.zeros(self.operatorModesSpace[0].shape, dtype=self.dtype)

            for k in range(len(self.rhsModesSpace)):
                Bcurr += self.rhsModesSpace[k] * self.rhsModesParam[k](mui)

            for k in range(len(self.operatorModesSpace)):
                Acurr += self.operatorModesSpace[k] * self.operatorModesParam[k](mui)

            Residual_mui = Bcurr - Acurr @ (basis @ coeffs[imu])

            norm_Rmu = np.linalg.norm(Residual_mui)

            # BIDOUILLE : CHECK
            # if imu in self.worstParamHistory:
            #     norm_Rmu *= 0.
            #     print('Bidouille !')


            # print('mu=',mui,'Norm resid',norm_Rmu)

            if norm_Rmu > maxResidNorm:
                maxResidNorm = norm_Rmu
                maxResid = Residual_mui
                imuTilde = imu

                opMaxResid = Acurr
                rhsMaxResid = Bcurr
                # print('max')

        print('Max residual norm', maxResidNorm, 'for snapshot', imuTilde)

        return maxResid, opMaxResid, rhsMaxResid, imuTilde



    # For the current reduced basis, compute the Galerkin coefficients
    def solveRomEqn(self, basis):

        m = basis.shape[1]
        coeffs = np.zeros((self.parameterTrainingSpace.shape[0], basis.shape[1]), dtype=self.dtype)

        # print(basis, basis.shape)

        ###### Construction du problème réduit (partie deterministe)

        # Construct Atilde_k and Btilde_k ONLY IF NEEDED
        if m != self.TildeComputedFor:
            # Construction des Atilde_k
            self.Atilde_k = []

            for k in range(len(self.operatorModesSpace)):
                Atilde = basis[:, :m].T @ (self.operatorModesSpace[k] @ basis[:, :m])
                self.Atilde_k.append(Atilde)

            # Construction des Btilde_k
            self.Btilde_k = []

            for k in range(len(self.rhsModesSpace)):
                Btilde = basis[:, :m].T @self.rhsModesSpace[k]
                self.Btilde_k.append(Btilde)

            self.TildeComputedFor = m
        else :
            print('Atilde and Btilde do already precomputed')



        ###### Resolution du problème réduit pour les snapshots
        # Loop over mu_i
        for imu, mui in enumerate(self.parameterTrainingSpace):


            # Particularisation du problème réduit pour mu
            RB_Matrix = np.zeros((m, m), dtype=self.dtype)
            RB_RHS = np.zeros(m, dtype=self.dtype)

            for k in range(len(self.operatorModesParam)):
                RB_Matrix += self.operatorModesParam[k](mui) * self.Atilde_k[k]

            for k in range(len(self.rhsModesParam)):
                RB_RHS += self.rhsModesParam[k](mui) * self.Btilde_k[k]
            # aaaaa
            # Solve
            coeffs[imu] = np.linalg.solve(RB_Matrix, RB_RHS)

        return coeffs



    def updateBasis(self, op, rhs):
        # Nouvelle fonction de base
        newBasis = np.linalg.solve(op, rhs)
        # newBasis /= np.linalg.norm(newBasis)
        return newBasis



    # Solve Greedy
    def solveGreedyLeblond(self, tol):
        m = 0
        err = 1.
        errorHistory = [];

        while err > tol and m < self.maxModes:

            if m == 0:
                # Calcul du résidu
                maxResid, Acurr, Bcurr, imuTilde = self.computeWorstParameterAndAssociatedData(self.RB_Basis[:, :m + 1],
                                                                                               self.RB_Coeffs[:, :m + 1])

            self.worstParamHistory.append(imuTilde)
            newBasis = self.updateBasis(Acurr, maxResid)

            # MAJ des bases reduites
            # gramm schmit orthogonalization
            othoNewBasis = newBasis - np.sum(np.dot(self.RB_Basis[:, j], newBasis)/np.dot(self.RB_Basis[:, j], self.RB_Basis[:, j]) *self.RB_Basis[:, j]  for j in range(m))
            othoNewBasis /= np.linalg.norm(othoNewBasis)
            
            self.RB_Basis[:, m] = othoNewBasis
            # self.RB_Basis[:, m], R = np.linalg.qr(self.RB_Basis[:, m])
            m += 1
            self.nbModes = m

            # Projection de la solution sur la base reduite
            self.RB_Coeffs[:, :m] = self.solveRomEqn(self.RB_Basis[:, :m])

            # Calcul du résidu
            maxResid, Acurr, Bcurr, imuTilde = self.computeWorstParameterAndAssociatedData(self.RB_Basis[:, :m],
                                                                                               self.RB_Coeffs[:, :m])
            # maxResid_ = RB_Basis[:,:m]@RB_Coeffs[:,:m].T

            # Erreur de la boucle while
            err = np.max(np.abs(maxResid))
            errorHistory.append(err)

            print('Mode', m, 'Error', err, 'Worst snapshot', imuTilde, '\n')

        return errorHistory

    # Solve Greedy
    def solveGreedy(self, tol):
        m = 0
        err = 1.
        errorHistory = [];

        while err > tol and m < self.maxModes:

            if m == 0:
                # Calcul du résidu
                maxResid, Acurr, Bcurr, imuTilde = self.computeWorstParameterAndAssociatedData(self.RB_Basis[:, :m + 1],
                                                                                               self.RB_Coeffs[:,
                                                                                               :m + 1])

            self.worstParamHistory.append(imuTilde)
            newBasis = self.updateBasis(Acurr, Bcurr)

            # orthonormalization
            othoNewBasis = newBasis - np.sum(np.dot(self.RB_Basis[:, j], newBasis)/np.dot(self.RB_Basis[:, j], self.RB_Basis[:, j]) *self.RB_Basis[:, j]  for j in range(m))
            othoNewBasis /= np.linalg.norm(othoNewBasis)
            # MAJ des bases reduites
            self.RB_Basis[:, m] = othoNewBasis
            m += 1
            self.nbModes = m

            # Projection de la solution sur la base reduite
            self.RB_Coeffs[:, :m] = self.solveRomEqn(self.RB_Basis[:, :m])

            # aaaa

            # Calcul du résidu
            maxResid, Acurr, Bcurr, imuTilde = self.computeWorstParameterAndAssociatedData(self.RB_Basis[:, :m],
                                                                                           self.RB_Coeffs[:, :m])
            # maxResid_ = RB_Basis[:,:m]@RB_Coeffs[:,:m].T

            # Erreur de la boucle while
            err = np.max(np.abs(maxResid))
            errorHistory.append(err)

            print('Mode', m, 'Error', err, 'Worst snapshot', imuTilde, '\n')

        return errorHistory


    # mu is part of the training space: we just take its id
    def reconstructHiFiApproximatedSolution(self, imu):
        return self.RB_Basis[:, :self.nbModes] @ self.RB_Coeffs[imu, :self.nbModes]

    # mu was not part of the training space here: we take its value
    def reconstructHiFiApproximatedSolutionArbitrary(self, mu):

        ###### Resolution du problème réduit pour le parametre donne
        # Particularisation du problème réduit pour mu
        m = self.nbModes

        # Check if the pre-computation of the reduced operators is consistent
        if m != self.TildeComputedFor:
            raise

        RB_Matrix = np.zeros((m, m), dtype=self.dtype)
        RB_RHS = np.zeros(m, dtype=self.dtype)

        for k in range(len(self.operatorModesParam)):
            RB_Matrix += self.operatorModesParam[k](mu) * self.Atilde_k[k]

        for k in range(len(self.rhsModesParam)):
            RB_RHS += self.rhsModesParam[k](mu) * self.Btilde_k[k]

        # Solve
        coeff = np.linalg.solve(RB_Matrix, RB_RHS)

        return self.RB_Basis[:, :self.nbModes] @ coeff


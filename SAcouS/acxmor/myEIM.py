
import numpy as np
from numpy import unravel_index
from numpy.linalg import solve
from numpy.linalg import inv


class nonIntrusiveEIMV2:
    def __init__(self, snapshots, maxRank, tolerance=1.e-8, orthogonalizeBasis=False, dtype_=float):
        self.snapshots = snapshots
        self.maxRank = maxRank
        self.tolerance = tolerance
        self.orthogonalizeBasis = orthogonalizeBasis
        self.display = True
        self.metric = lambda x: abs(x).max(axis=0)
        self.indXi = 0
        self.indX = 0
        self.rank = 0
        self.dtype = dtype_

    def vecMax(self, arr):
        return arr.max(), arr.argmax()

    def interpolate(self, intrusive=False):
        V = np.zeros((self.snapshots.shape[0], self.maxRank), dtype=self.dtype)
        indXi = np.zeros(self.maxRank, int)
        indX = np.zeros(self.maxRank, int)
        fullApprox = np.zeros(self.snapshots.shape, dtype=self.dtype)
        residual = self.snapshots - fullApprox
        normResidual, indXi[0] = self.vecMax(self.metric(residual))
        trash, indX[0] = self.vecMax(abs(residual[:, indXi[0]]))
        l = 0

        if self.display:
            print('Iteration %3d - Residual %e\n' % (l, normResidual))
            print(indXi)

        while (self.maxRank > l) and (normResidual > self.tolerance):
            v = residual[:, indXi[l]]
            v = v / v[indX[l]]

            V[:, l] = v[:]
            Vl = V[:, :l + 1]

            lambda_ = solve(Vl[indX[:l + 1], :], self.snapshots[indX[:l + 1], :])  # gamma

            fullApprox = np.dot(Vl, lambda_)
            lastResidual = residual
            residual = self.snapshots - fullApprox

            l += 1

            # Pour le dernier mode, on ne fait que calculer le residu car in n'ira pas plus loin dans l'algo
            if l < self.maxRank:
                normResidual, indXi[l] = self.vecMax(self.metric(residual))
                trash, indX[l] = self.vecMax(abs(residual[:, indXi[l]]))
            else:
                normResidual, trash = self.vecMax(self.metric(residual))
            #                aaa
            #            print(indXi)

            if self.display:
                print('Iteration %3d - Residual %e\n' % (l, normResidual))
                # print(indXi)

        if intrusive:
            self.indXi = indXi
            self.indX = indX
            self.rank = l
            return V[:, :l], lambda_.T
        else:
            self.indXi = indXi[:l]
            self.indX = indX[:l]
            self.rank = l
            lambda_ = self.snapshots[:, indXi[:l]];
            R = self.snapshots[:, indXi[:l]];
            lambda_ = solve(lambda_[indX[:l], :], self.snapshots[indX[:l], :]);
	    # R basis, lambda is the parametric dependence function
            return R, lambda_.T





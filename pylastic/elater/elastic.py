import json
from pylastic.utils.vectors import *

import numpy as np

from scipy import optimize


class Elastic:
    """An elastic tensor, along with methods to access it"""

    def __init__(self, s):
        """Initialize the elastic tensor from a string"""

        if not s:
            raise ValueError("no matrix was provided")

        # Argument can be a 6-line string, a list of list, or a string representation of the list of list
        try:
            if type(json.loads(s)) == list: s = json.loads(s)
        except:
            pass

        if type(s) == str:
            # Remove braces and pipes
            s = s.replace("|", " ").replace("(", " ").replace(")", " ")

            # Remove empty lines
            lines = [line for line in s.split('\n') if line.strip()]
            if len(lines) != 6:
                raise ValueError("should have six rows")

            # Convert to float
            try:
                mat = [list(map(float, line.split())) for line in lines]
            except:
                raise ValueError("not all entries are numbers")
        elif type(s) == list:
            # If we already have a list, simply use it
            mat = s
        else:
            raise ValueError("invalid argument as matrix")

        # Make it into a square matrix
        try:
            mat = np.array(mat)
        except:
            # Is it upper triangular?
            if list(map(len, mat)) == [6, 5, 4, 3, 2, 1]:
                mat = [[0] * i + mat[i] for i in range(6)]
                mat = np.array(mat)

            # Is it lower triangular?
            if list(map(len, mat)) == [1, 2, 3, 4, 5, 6]:
                mat = [mat[i] + [0] * (5 - i) for i in range(6)]
                mat = np.array(mat)

        if not isinstance(mat, np.ndarray):
            raise ValueError("should be a square or triangular matrix")
        if mat.shape != (6, 6):
            raise ValueError("should be a square or triangular matrix")

        # Check that is is symmetric, or make it symmetric
        if np.linalg.norm(np.tril(mat, -1)) == 0:
            mat = mat + np.triu(mat, 1).transpose()
        if np.linalg.norm(np.triu(mat, 1)) == 0:
            mat = mat + np.tril(mat, -1).transpose()
        if np.linalg.norm(mat - mat.transpose()) > 1e-3:
            raise ValueError("should be symmetric, or triangular")
        elif np.linalg.norm(mat - mat.transpose()) > 0:
            mat = 0.5 * (mat + mat.transpose())

        # Store it
        self.CVoigt = mat

        # Put it in a more useful representation
        try:
            self.SVoigt = np.linalg.inv(self.CVoigt)
        except:
            raise ValueError("matrix is singular")

        VoigtMat = [[0, 5, 4], [5, 1, 3], [4, 3, 2]]

        def SVoigtCoeff(p, q):
            return 1. / ((1 + p // 3) * (1 + q // 3))

        self.Smat = [[[[SVoigtCoeff(VoigtMat[i][j], VoigtMat[k][l]) * self.SVoigt[VoigtMat[i][j]][VoigtMat[k][l]]
                        for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]
        return

    def isOrthorhombic(self):
        def iszero(x): return (abs(x) < 1.e-3)

        return (iszero(self.CVoigt[0][3]) and iszero(self.CVoigt[0][4]) and iszero(self.CVoigt[0][5])
                and iszero(self.CVoigt[1][3]) and iszero(self.CVoigt[1][4]) and iszero(self.CVoigt[1][5])
                and iszero(self.CVoigt[2][3]) and iszero(self.CVoigt[2][4]) and iszero(self.CVoigt[2][5])
                and iszero(self.CVoigt[3][4]) and iszero(self.CVoigt[3][5]) and iszero(self.CVoigt[4][5]))

    def isCubic(self):
        def iszero(x): return (abs(x) < 1.e-3)

        return (iszero(self.CVoigt[0][3]) and iszero(self.CVoigt[0][4]) and iszero(self.CVoigt[0][5])
                and iszero(self.CVoigt[1][3]) and iszero(self.CVoigt[1][4]) and iszero(self.CVoigt[1][5])
                and iszero(self.CVoigt[2][3]) and iszero(self.CVoigt[2][4]) and iszero(self.CVoigt[2][5])
                and iszero(self.CVoigt[3][4]) and iszero(self.CVoigt[3][5]) and iszero(self.CVoigt[4][5])
                and iszero(self.CVoigt[0][0] - self.CVoigt[1][1]) and iszero(self.CVoigt[0][0] - self.CVoigt[2][2])
                and iszero(self.CVoigt[0][0] - self.CVoigt[1][1]) and iszero(self.CVoigt[0][0] - self.CVoigt[2][2])
                and iszero(self.CVoigt[3][3] - self.CVoigt[4][4]) and iszero(self.CVoigt[3][3] - self.CVoigt[5][5])
                and iszero(self.CVoigt[0][1] - self.CVoigt[0][2]) and iszero(self.CVoigt[0][1] - self.CVoigt[1][2]))

    def Young(self, x):
        a = dirVec(x[0], x[1])
        r = sum([a[i] * a[j] * a[k] * a[l] * self.Smat[i][j][k][l]
                 for i in range(3) for j in range(3) for k in range(3) for l in range(3)])
        return 1 / r

    def Young_2(self, x, y):
        a = dirVec(x, y)
        r = sum([a[i] * a[j] * a[k] * a[l] * self.Smat[i][j][k][l]
                 for i in range(3) for j in range(3) for k in range(3) for l in range(3)])
        return 1 / r

    def LC(self, x):
        a = dirVec(x[0], x[1])
        r = sum([a[i] * a[j] * self.Smat[i][j][k][k]
                 for i in range(3) for j in range(3) for k in range(3)])
        return 1000 * r

    def LC_2(self, x, y):
        a = dirVec(x, y)
        r = sum([a[i] * a[j] * self.Smat[i][j][k][k]
                 for i in range(3) for j in range(3) for k in range(3)])
        return 1000 * r

    def shear(self, x):
        a = dirVec(x[0], x[1])
        b = dirVec2(x[0], x[1], x[2])
        r = sum([a[i] * b[j] * a[k] * b[l] * self.Smat[i][j][k][l]
                 for i in range(3) for j in range(3) for k in range(3) for l in range(3)])
        return 1 / (4 * r)

    def Poisson(self, x):
        a = dirVec(x[0], x[1])
        b = dirVec2(x[0], x[1], x[2])
        r1 = sum([a[i] * a[j] * b[k] * b[l] * self.Smat[i][j][k][l]
                  for i in range(3) for j in range(3) for k in range(3) for l in range(3)])
        r2 = sum([a[i] * a[j] * a[k] * a[l] * self.Smat[i][j][k][l]
                  for i in range(3) for j in range(3) for k in range(3) for l in range(3)])
        return -r1 / r2

    def averages(self):
        A = (self.CVoigt[0][0] + self.CVoigt[1][1] + self.CVoigt[2][2]) / 3
        B = (self.CVoigt[1][2] + self.CVoigt[0][2] + self.CVoigt[0][1]) / 3
        C = (self.CVoigt[3][3] + self.CVoigt[4][4] + self.CVoigt[5][5]) / 3
        a = (self.SVoigt[0][0] + self.SVoigt[1][1] + self.SVoigt[2][2]) / 3
        b = (self.SVoigt[1][2] + self.SVoigt[0][2] + self.SVoigt[0][1]) / 3
        c = (self.SVoigt[3][3] + self.SVoigt[4][4] + self.SVoigt[5][5]) / 3

        KV = (A + 2 * B) / 3
        GV = (A - B + 3 * C) / 5

        KR = 1 / (3 * a + 6 * b)
        GR = 5 / (4 * a - 4 * b + 3 * c)

        KH = (KV + KR) / 2
        GH = (GV + GR) / 2

        return [[KV, 1 / (1 / (3 * GV) + 1 / (9 * KV)), GV, (1 - 3 * GV / (3 * KV + GV)) / 2],
                [KR, 1 / (1 / (3 * GR) + 1 / (9 * KR)), GR, (1 - 3 * GR / (3 * KR + GR)) / 2],
                [KH, 1 / (1 / (3 * GH) + 1 / (9 * KH)), GH, (1 - 3 * GH / (3 * KH + GH)) / 2]]

    def shear2D(self, x):
        ftol = 0.001
        xtol = 0.01

        def func1(z): return self.shear([x[0], x[1], z])

        r1 = optimize.minimize(func1, np.pi / 2.0, args=(), method='Powell',
                               options={"xtol": xtol, "ftol": ftol})  # , bounds=[(0.0,np.pi)])

        def func2(z): return -self.shear([x[0], x[1], z])

        r2 = optimize.minimize(func2, np.pi / 2.0, args=(), method='Powell',
                               options={"xtol": xtol, "ftol": ftol})  # , bounds=[(0.0,np.pi)])
        return (float(r1.fun), -float(r2.fun))

    def shear3D(self, x, y, guess1=np.pi / 2.0, guess2=np.pi / 2.0):
        tol = 0.005

        def func1(z): return self.shear([x, y, z])

        r1 = optimize.minimize(func1, guess1, args=(), method='COBYLA', options={"tol": tol})  # , bounds=[(0.0,np.pi)])

        def func2(z): return -self.shear([x, y, z])

        r2 = optimize.minimize(func2, guess2, args=(), method='COBYLA', options={"tol": tol})  # , bounds=[(0.0,np.pi)])
        return (float(r1.fun), -float(r2.fun), float(r1.x), float(r2.x))

    def Poisson2D(self, x):
        ftol = 0.001
        xtol = 0.01

        def func1(z): return self.Poisson([x[0], x[1], z])

        r1 = optimize.minimize(func1, np.pi / 2.0, args=(), method='Powell',
                               options={"xtol": xtol, "ftol": ftol})  # , bounds=[(0.0,np.pi)])

        def func2(z): return -self.Poisson([x[0], x[1], z])

        r2 = optimize.minimize(func2, np.pi / 2.0, args=(), method='Powell',
                               options={"xtol": xtol, "ftol": ftol})  # , bounds=[(0.0,np.pi)])
        return (min(0, float(r1.fun)), max(0, float(r1.fun)), -float(r2.fun))

    def poisson3D(self, x, y, guess1=np.pi / 2.0, guess2=np.pi / 2.0):
        tol = 0.005

        def func1(z): return self.Poisson([x, y, z])

        r1 = optimize.minimize(func1, guess1, args=(), method='COBYLA', options={"tol": tol})  # , bounds=[(0.0,np.pi)])

        def func2(z): return -self.Poisson([x, y, z])

        r2 = optimize.minimize(func2, guess2, args=(), method='COBYLA', options={"tol": tol})  # , bounds=[(0.0,np.pi)])
        return (min(0, float(r1.fun)), max(0, float(r1.fun)), -float(r2.fun), float(r1.x), float(r2.x))


import sys
import io
from sympy import sqrt, exp, N, S, sin, cos
from sympy.matrices import Matrix
from numpy import float32

canonicalShapes = {"coplanes": Matrix([[1.0, 0.0, 0.0, 0.0],
                                       [0.0, 0.0, 0.0, 0.0],
                                       [0.0, 0.0, 0.0, 0.0],
                                       [0.0, 0.0, 0.0, 0.0]]),
                   "parplanes": Matrix([[1.0, 0.0, 0.0, 0.0],
                                        [0.0, 0.0, 0.0, 0.0],
                                        [0.0, 0.0, 0.0, 0.0],
                                        [0.0, 0.0, 0.0, 1.0]]),
                   "interplanes": Matrix([[1.0, 0.0, 0.0, 0.0],
                                          [0.0, -1.0, 0.0, 0.0],
                                          [0.0, 0.0, 0.0, 0.0],
                                          [0.0, 0.0, 0.0, 0.0]]),
                   "ellcylinder": Matrix([[1.0, 0.0, 0.0, 0.0],
                                          [0.0, 1.0, 0.0, 0.0],
                                          [0.0, 0.0, 0.0, 0.0],
                                          [0.0, 0.0, 0.0, -1.0]]),
                   "hypcylinder": Matrix([[1.0, 0.0, 0.0, 0.0],
                                          [0.0, -1.0, 0.0, 0.0],
                                          [0.0, 0.0, 0.0, 0.0],
                                          [0.0, 0.0, 0.0, -1.0]]),
                   "ellcone": Matrix([[1.0, 0.0, 0.0, 0.0],
                                      [0.0, 1.0, 0.0, 0.0],
                                      [0.0, 0.0, -1.0, 0.0],
                                      [0.0, 0.0, 0.0, 0.0]]),
                   "ellipsoid": Matrix([[1.0, 0.0, 0.0, 0.0],
                                        [0.0, 1.0, 0.0, 0.0],
                                        [0.0, 0.0, 1.0, 0.0],
                                        [0.0, 0.0, 0.0, -1.0]]),
                   "hyperboloid1": Matrix([[1.0, 0.0, 0.0, 0.0],
                                           [0.0, 1.0, 0.0, 0.0],
                                           [0.0, 0.0, -1.0, 0.0],
                                           [0.0, 0.0, 0.0, 1.0]]),
                   "hyperboloid2": Matrix([[1.0, 0.0, 0.0, 0.0],
                                           [0.0, 1.0, 0.0, 0.0],
                                           [0.0, 0.0, -1.0, 0.0],
                                           [0.0, 0.0, 0.0, -1.0]]),
                   "parcylinder": Matrix([[1.0, 0.0, 0.0, 0.0],
                                          [0.0, 0.0, 0.0, 0.0],
                                          [0.0, 0.0, 0.0, 0.5],
                                          [0.0, 0.0, 0.5, 0.0]]),
                   "ellparaboloid": Matrix([[1.0, 0.0, 0.0, 0.0],
                                            [0.0, 1.0, 0.0, 0.0],
                                            [0.0, 0.0, 0.0, 0.5],
                                            [0.0, 0.0, 0.5, 0.0]]),
                   "hypparaboloid": Matrix([[-1.0, 0.0, 0.0, 0.0],
                                            [0.0, 1.0, 0.0, 0.0],
                                            [0.0, 0.0, 0.0, 0.5],
                                            [0.0, 0.0, 0.5, 0.0]]),
}

def applyTForm(quad, tform):
    return tform.T * quad * tform

def translateMtx(x, y, z):
    return Matrix([[1.0, 0.0, 0.0, x],
                   [0.0, 1.0, 0.0, y],
                   [0.0, 0.0, 1.0, z],
                   [0.0, 0.0, 0.0, 1.0]])

def scaleXMtx(scale):
    return Matrix([[scale, 0.0, 0.0, 0.0],
                   [0.0, 1.0, 0.0, 0.0],
                   [0.0, 0.0, 1.0, 0.0],
                   [0.0, 0.0, 0.0, 1.0]])

def scaleYMtx(scale):
    return Matrix([[1.0, 0.0, 0.0, 0.0],
                   [0.0, scale, 0.0, 0.0],
                   [0.0, 0.0, 1.0, 0.0],
                   [0.0, 0.0, 0.0, 1.0]])

def scaleZMtx(scale):
    return Matrix([[1.0, 0.0, 0.0, 0.0],
                   [0.0, 1.0, 0.0, 0.0],
                   [0.0, 0.0, scale, 0.0],
                   [0.0, 0.0, 0.0, 1.0]])

def scaleMtx(scaleX, scaleY, scaleZ):
    return Matrix([[scaleX, 0.0, 0.0, 0.0],
                   [0.0, scaleY, 0.0, 0.0],
                   [0.0, 0.0, scaleZ, 0.0],
                   [0.0, 0.0, 0.0, 1.0]])

def scaleAllMtx(scale):
    return Matrix([[1.0, 0.0, 0.0, 0.0],
                   [0.0, 1.0, 0.0, 0.0],
                   [0.0, 0.0, 1.0, 0.0],
                   [0.0, 0.0, 0.0, 1.0 / S(scale)]])

def rotateMtx(axis, angle):
    ct = cos(angle)
    st = sin(angle)
    tpMtx = Matrix([[S(axis[0]) ** 2, S(axis[0]) * S(axis[1]), S(axis[0]) * S(axis[2]), S(0)],
                    [S(axis[0]) * S(axis[1]), S(axis[1]) ** 2, S(axis[1]) * S(axis[2]), S(0)],
                    [S(axis[0]) * S(axis[2]), S(axis[1]) * S(axis[2]), S(axis[2]) ** 2, S(0)],
                    [S(0), S(0), S(0), S(0)]]) * S(1 - ct)
    return Matrix([[ct + S(axis[0]) ** 2 * S(1 - ct),
                    S(axis[0]) * S(axis[1]) * S(1 - ct) - S(axis[2]) * st,
                    S(axis[0]) * S(axis[2]) * S(1 - ct) + S(axis[1]) * st,
                    0]
                   [S(axis[0]) * S(axis[1]) * S(1 - ct) - S(axis[2]) * st,
                    ct + S(axis[1]) ** 2 * S(1 - ct),
                    S(axis[1]) * S(axis[2]) * S(1 - ct) - S(axis[0]) * st,
                    0]
                   [S(axis[0]) * S(axis[2]) * S(1 - ct) - S(axis[1]) * st,
                    S(axis[1]) * S(axis[2]) * S(1 - ct) + S(axis[0]) * st,
                    ct + S(axis[2]) ** 2 * S(1 - ct),
                    S(0)]
                   [S(0), S(0), S(0), S(1)]])

def writeScene(fname, quads):
    f = io.open(fname, "w")
    f.write("-20 1\n")
    i = 0
    for q in quads:
        if i > 0:
            f.write("\n")
        f.write("m")
        for y in range(q.shape[0]):
            val = float32(q[y * q.shape[1]])
            vstr = "{:.30f}".format(val).rstrip("0").rstrip(".")
            f.write("\n" + vstr)
            for x in range(1, q.shape[1]):
                val = float32(q[y * q.shape[1] + x])
                vstr = "{:.30f}".format(val).rstrip("0").rstrip(".")
                f.write(" " + vstr)
        i += 1

def genAxialCylinders(numCyls, eps = (2 ** -8)):
    defCyl = canonicalShapes["ellcylinder"]
    scene = [applyTForm(applyTForm(defCyl,
                                   scaleMtx(2.0 * (N(1.0) + N(i) * N(eps) - N(eps)),
                                            2.0 * (N(1.0) + N(i) * N(i) * N(eps) - N(eps)),
                                            1.0)),
                        translateMtx(-0.5, -0.5, 0.0))
             for i in range(1, numCyls + 1)]
    return scene

def genSingleAmbiguousEllipsoids(numEllipsoids):
    x0, y0, z0 = 0.5, 0.5, 0.5
    radius0 = 0.5
    epsilon = 2 ** -19
    delta = (radius0 - (numEllipsoids + 8) * epsilon) / numEllipsoids
    defEll = canonicalShapes["ellipsoid"]
    scene = []
    for i in range(numEllipsoids):
        x = x0 + i * delta
        y = y0
        z = z0
        radius = radius0 - (i * delta + i * epsilon)
        scene.append(applyTForm(applyTForm(defEll,
                                           scaleAllMtx(1 / S(radius))),
                                translateMtx(-x, -y, -z)))
    return scene

def genCenteredEllipsoids(numElls, eps = (2 ** -8)):
    defEll = canonicalShapes["ellipsoid"]
    scene = [applyTForm(applyTForm(defEll,
                                   scaleMtx(2.0 * (N(1.0) +
                                                   N(i) * N(eps) -
                                                   N(eps)),
                                            2.0 * (N(1.0) +
                                                   sqrt(N(i) * N(i)) * N(eps) -
                                                   N(eps)),
                                            2.0 * (N(1.0)
                                                   + N(i) * N(i) *
                                                   N(eps) - N(eps)))),
                        translateMtx(-0.5, -0.5, -0.5))
             for i in range(1, numElls + 1)]
    return scene

if __name__ == "__main__":
    print("Generating " + sys.argv[2] +
          " centered ellipsoids in " +
          sys.argv[1])
    cyls = genSingleAmbiguousEllipsoids(int(sys.argv[2]))
    writeScene(sys.argv[1], cyls)

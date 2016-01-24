
import sys
import io
from sympy import exp, N, S
from sympy.matrices import Matrix

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
                                            [0.0, 0.0, 0.0, 1.0]]),
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
                                          [0.0, 0.0, 0.0, 1.0]]),
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

def rotateMtx(axis, angle):
    ct = np.cos(angle)
    st = np.sin(angle)
    tpMtx = Matrix([[axis[0] ** 2, axis[0] * axis[1], axis[0] * axis[2], 0],
                    [axis[0] * axis[1], axis[1] ** 2, axis[1] * axis[2], 0],
                    [axis[0] * axis[2], axis[1] * axis[2], axis[2] ** 2, 0],
                    [0, 0, 0, 0]]) * (1 - ct)
    return Matrix([[ct + axis[0] ** 2 * (1 - ct),
                    axis[0] * axis[1] * (1 - ct) - axis[2] * st,
                    axis[0] * axis[2] * (1 - ct) + axis[1] * st,
                    0]
                   [axis[0] * axis[1] * (1 - ct) - axis[2] * st,
                    ct + axis[1] ** 2 * (1 - ct),
                    axis[1] * axis[2] * (1 - ct) - axis[0] * st,
                    0]
                   [axis[0] * axis[2] * (1 - ct) - axis[1] * st,
                    axis[1] * axis[2] * (1 - ct) + axis[0] * st,
                    ct + axis[2] ** 2 * (1 - ct),
                    0]
                   [0, 0, 0, 1]])

def writeScene(fname, quads):
    f = io.open(fname, "w")
    i = 0
    for q in quads:
        if i > 0:
            f.write("\n")
        f.write("m")
        for y in range(q.shape[0]):
            val = float(q[y * q.shape[1]])
            vstr = "{:.30f}".format(val)
            f.write("\n" + vstr)
            for x in range(1, q.shape[1]):
                val = float(q[y * q.shape[1] + x])
                vstr = "{:.30f}".format(val)
                f.write(" " + vstr)
        i += 1

def genAxialCylinders(numCyls, eps = (2 ** -20)):
    defCyl = canonicalShapes["ellcylinder"]
    scene = [defCyl * scaleXMtx(1 + i * eps) for i in range(numCyls)]
    return scene

if __name__ == "__main__":
    cyls = genAxialCylinders(int(sys.argv[2]))
    writeScene(sys.argv[1], cyls)

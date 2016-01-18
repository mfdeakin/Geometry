
import sys
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
    return Matrix([[S("1.0 / scale"), 0.0, 0.0, 0.0],
                     [0.0, 1.0, 0.0, 0.0],
                     [0.0, 0.0, 1.0, 0.0],
                     [0.0, 0.0, 0.0, 1.0]])

def scaleYMtx(scale):
    return Matrix([[1.0, 0.0, 0.0, 0.0],
                     [0.0, 1.0 / scale, 0.0, 0.0],
                     [0.0, 0.0, 1.0, 0.0],
                     [0.0, 0.0, 0.0, 1.0]])

def scaleZMtx(scale):
    return Matrix([[1.0, 0.0, 0.0, 0.0],
                     [0.0, 1.0, 0.0, 0.0],
                     [0.0, 0.0, 1.0 / scale, 0.0],
                     [0.0, 0.0, 0.0, 1.0]])

def scaleMtx(scaleX, scaleY, scaleZ):
    return Matrix([[1.0 / scaleX, 0.0, 0.0, 0.0],
                     [0.0, 1.0 / scaleY, 0.0, 0.0],
                     [0.0, 0.0, 1.0 / scaleZ, 0.0],
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


if __name__ == "__main__":
    print(canonicalShapes[sys.argv[1]])

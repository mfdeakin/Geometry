import csv
import os
import numpy as np
import matplotlib.pyplot as plt

files = os.listdir()

resLSQ = {}
resLSQAvg = {}
resLSQVar = {}
res5th = {}
res95th = {}

incLSQ = {}
incLSQAvg = {}
incLSQVar = {}
inc5th = {}
inc95th = {}

for fname in files:
    testInfo = fname.split('.')
    if len(testInfo) != 5:
        continue
    testName, testSize, _, testMachine, fileType = fname.split('.')
    if fileType != 'csv':
        continue
    testSize = int(testSize)
    if not testMachine in resLSQ.keys():
        resLSQ[testMachine] = {}
        incLSQ[testMachine] = {}
        resLSQAvg[testMachine] = {}
        incLSQAvg[testMachine] = {}
        resLSQVar[testMachine] = {}
        incLSQVar[testMachine] = {}
    if not testName in resLSQ[testMachine].keys():
        resLSQ[testMachine][testName] = {}
        incLSQ[testMachine][testName] = {}
        resLSQAvg[testMachine][testName] = np.array([0.0, 0.0])
        incLSQAvg[testMachine][testName] = np.array([0.0, 0.0])
        resLSQVar[testMachine][testName] = np.array([0.0, 0.0])
        incLSQVar[testMachine][testName] = np.array([0.0, 0.0])
    if not testSize in resLSQ[testMachine][testName].keys():
        resLSQ[testMachine][testName][testSize] = []
        incLSQ[testMachine][testName][testSize] = []
    csvFile = open(fname, 'r')
    csvRead = csv.reader(csvFile)
    resTimes = []
    incTimes = []
    precIncreases = []
    numCorrect = 0
    numLines = 0
    next(csvRead)
    for row in csvRead:
        if len(row) != 5:
            break
        resTimes.append(np.float64(row[1]))
        incTimes.append(np.float64(row[2]))
        numCorrect += int(row[3])
        precIncreases.append(int(row[4]))
        numLines += 1
    resLSqr = np.linalg.lstsq(np.vstack([precIncreases,
                                         np.ones(len(precIncreases))]).T,
                              resTimes)
    mRes, cRes = resLSqr[0]
    resLSQ[testMachine][testName][testSize].append(np.array([mRes, cRes]))
    incLSqr = np.linalg.lstsq(np.vstack([precIncreases,
                                         np.ones(len(precIncreases))]).T,
                              incTimes)
    mInc, cInc = incLSqr[0]
    incLSQ[testMachine][testName][testSize].append(np.array([mInc, cInc]))

for machine in resLSQ:
    for test in resLSQ[machine]:
        mTot = 0.0
        cTot = 0.0
        for size in resLSQ[machine][test]:
            print(test, machine, resLSQ[machine][test][size])
            mTot += resLSQ[machine][test][size][0][0]
            cTot += resLSQ[machine][test][size][0][1]
        resLSQAvg[machine][test][0] = mTot / len(resLSQ[machine][test])
        resLSQAvg[machine][test][1] = cTot / len(resLSQ[machine][test])
        print(resLSQAvg[machine][test])
        exit()

#    spRes = plt.scatter(precIncreases, resTimes,
#                        c = "blue", label = "Resultant Method")
#    spInc = plt.scatter(precIncreases, incTimes,
#                        c = "red", label = "Naive Method")
#    plt.legend(bbox_to_anchor=(0.95, 0.25))
#    print(fname)
#    print("% Correct: " + str(numCorrect) + " / " + str(numLines))
#    print(m, c)
#    plt.show()

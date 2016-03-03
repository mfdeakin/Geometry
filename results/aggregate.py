import csv
import os
import numpy as np
import matplotlib.pyplot as plt

files = os.listdir()

resAvg = {}
resVar = {}
res5th = {}
res95th = {}

for fname in files:
    testInfo = fname.split('.')
    if len(testInfo) != 5:
        continue
    testName, testSize, _, testMachine, fileType = fname.split('.')
    if fileType != 'csv':
        continue
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
        resTimes.append(np.float128(row[1]))
        incTimes.append(np.float128(row[2]))
        numCorrect += int(row[3])
        precIncreases.append(int(row[4]))
        numLines += 1
    spRes = plt.scatter(precIncreases, resTimes,
                        c = "blue", label = "Resultant Method")
    spInc = plt.scatter(precIncreases, incTimes,
                        c = "red", label = "Naive Method")
    plt.legend(bbox_to_anchor=(0.95, 0.25))
    print(fname)
    print("% Correct: " + str(numCorrect) + " / " + str(numLines))
    plt.show()

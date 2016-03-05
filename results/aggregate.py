import csv
import os
import numpy as np
import matplotlib.pyplot as plt
from pylab import savefig

def readData(fname, size):
    fd = open(fname, 'r')
    reader = csv.reader(fd)
    next(reader)
    data = []
    for row in reader:
        if len(row) != 5:
            continue
        number, resTime, mpTime, correct, precIncreases = map(int, row)
        data.append([size, precIncreases, resTime, mpTime, correct, number])
    return np.array(data)

def validateFName(fname):
    testInfo = fname.split('.')
    if len(testInfo) != 5:
        return None
    testName, testSize, _, testMachine, fileType = fname.split('.')
    if fileType != 'csv':
        return None
    try:
        testSize = int(testSize)
    except ValueError:
        return None
    return [testMachine, testName, testSize]

def aggregateData():
    files = os.listdir()
    datum = {}
    for fname in files:
        dataInfo = validateFName(fname)
        if dataInfo == None:
            continue
        print(fname)
        machine, test, size = dataInfo
        if not machine in datum:
            datum[machine] = {}
        if not test in datum[machine]:
            datum[machine][test] = []
        testData = datum[machine][test]
        data = readData(fname, size)
        precIncreases = data.T[1]
        resTimes = data.T[2]
        mpTimes = data.T[3]
        resLSqr = np.linalg.lstsq(np.vstack([precIncreases,
                                             np.ones(len(precIncreases))]).T,
                                  resTimes)
        mpLSqr = np.linalg.lstsq(np.vstack([precIncreases,
                                            np.ones(len(precIncreases))]).T,
                                  mpTimes)
        plt.scatter(precIncreases, resTimes, c = "blue",
                    label = "Resultant Method")
        plt.scatter(precIncreases, mpTimes, c = "red",
                    label = "Resultant Method")
        testData.append((size, data, resLSqr, mpLSqr))
        plt.show()

if __name__ == "__main__":
    aggregateData()

import csv
import os
import numpy as np
import matplotlib.pyplot as plt
import sklearn.decomposition

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

def createPlot(data, resLSqr, mpLSqr, fname):
    precIncreases = data.T[1]
    resTimes = data.T[2]
    mpTimes = data.T[3]
    plt.scatter(precIncreases, resTimes, c = "blue",
                label = "Resultant Method")
    plt.scatter(precIncreases, mpTimes, c = "red",
                label = "Resultant Method")
    plt.plot(precIncreases,
             precIncreases * resLSqr[0][0] + resLSqr[0][1],
             c = "cyan", linewidth = 4)
    plt.plot(precIncreases,
             precIncreases * mpLSqr[0][0] + mpLSqr[0][1],
             c = "purple", linewidth = 4)
    plt.axes().set_xlim(left = 0, right = max(precIncreases) * 1.0625)
    maxValue = max(max(precIncreases) * mpLSqr[0][0] + mpLSqr[0][1],
                   max(precIncreases) * resLSqr[0][0] + resLSqr[0][1],
                   max(mpTimes), max(resTimes))
    plt.axes().set_ylim(bottom = 0,
                        top =  maxValue * 1.0625)
    plt.yticks(np.arange(0, maxValue + 1, maxValue / 10.0))
    plt.axes().yaxis.get_major_formatter().set_powerlimits((0, 3))
    print("Saving '" + fname + "'")
    plt.savefig(fname, format = "png", dpi = 300)
    plt.close()

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
            datum[machine][test] = np.array([])
        testData = datum[machine][test]
        data = readData(fname, size)
        newShape = (testData.shape[0] + data.shape[0],
                    data.shape[1])
        datum[machine][test] = np.append(testData, data).reshape(newShape)
    return datum

def analyzeData(datum):
    for machine in datum:
        for test in datum[machine]:
            testData = datum[machine][test]
            precIncreases = testData.T[1]
            resTimes = testData.T[2]
            resLSqr = np.linalg.lstsq(np.vstack([precIncreases,
                                                 np.ones(len(precIncreases))]).T,
                                      resTimes)
            incTimes = testData.T[3]
            incLSqr = np.linalg.lstsq(np.vstack([precIncreases,
                                                 np.ones(len(precIncreases))]).T,
                                      incTimes)
            plotFName = test + "." + machine + ".results.png"
            createPlot(testData, resLSqr, incLSqr, plotFName)
            pca = sklearn.decomposition.PCA(6).fit(testData)
            print("Resultant Least Square Coefficients:", resLSqr[0])
            print("Increased Precision Coefficients:", incLSqr[0])
            print("PCA Components of " +
                  "[numQuadrics, numPrecIncreases, " +
                  "resultantTime_ns, incTime_ns, correct?, counter\n",
                  pca.components_[0])

if __name__ == "__main__":
    datum = aggregateData()
    analyzeData(datum)
    

import csv
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sklearn.decomposition

def readData(fname, size):
    fd = open(fname, 'r')
    reader = csv.reader(fd)
    next(reader)
    data = []
    for row in reader:
        if len(row) != 6:
            continue
        number, resNum, resTime, mpNum, mpTime, correct = map(int, row)
        data.append([size, resNum, resTime, mpNum, mpTime, correct, number])
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

def createPlot(data, fname):
    numQuads = data.T[0]
    resNum = data.T[1]
    resTimes = data.T[2]
    mpNum = data.T[3]
    mpTimes = data.T[4]
    plt.suptitle(fname)
    plt.scatter(resNum, resTimes, c = "blue",
               label = "Resultant Method")
    plt.scatter(mpNum, mpTimes, c = "red",
                label = "Increased Precision Method")
    maxValue = max(max(mpTimes), max(resTimes))
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
            datum[machine][test] = {}
        if not size in datum[machine][test]:
            datum[machine][test][size] = np.array([])
        testData = datum[machine][test][size]
        data = readData(fname, size)
        newShape = (testData.shape[0] + data.shape[0],
                    data.shape[1])
        datum[machine][test][size] = np.append(testData, data).reshape(newShape)
    return datum

def analyzeData(datum):
    for machine in datum:
        for test in datum[machine]:
            mtData = np.array([])
            for size in datum[machine][test]:
                testData = datum[machine][test][size]
                mtShape = (mtData.shape[0] + testData.shape[0],
                           testData.shape[1])
                mtData = np.append(mtData, testData)
                mtData = mtData.reshape(mtShape)
                resNum = testData.T[1]
                resTimes = testData.T[2]
                resIndep = np.vstack([resNum,
                                      np.ones(len(resNum))]).T
                resLSqr = np.linalg.lstsq(resIndep, resTimes)
                mpNum = testData.T[3]
                mpTimes = testData.T[4]
                mpIndep = np.vstack([mpNum,
                                     np.ones(len(mpNum))]).T
                mpLSqr = np.linalg.lstsq(mpIndep, mpTimes)
                plotFName = (test + "." + machine + "." +
                             str(size) + ".results.png")
                if len(datum[machine][test].keys()) == 1:
                    print("Test Results: " + machine + ", " + test)
                    createPlot(testData, plotFName)
                    try:
                        pca = sklearn.decomposition.PCA(6).fit(testData)
                        print("PCA Components of " +
                              "[numQuadrics, numResultants, " +
                              "resultantTime_ns, numMPComparisons, " +
                              "mpTime_ns, correct?, counter\n",
                              pca.components_)
                    except np.linalg.linalg.LinAlgError:
                        print("Could not compute the PCA of " + plotFName)
                    print("Resultant Least Square:", resLSqr)
                    print("Increased Precision Least Square:", mpLSqr)
                    print()
            if len(datum[machine][test].keys()) > 1:
                numQuads = mtData.T[0]
                resNum = mtData.T[1]
                resTimes = mtData.T[2]
                resIndep = np.vstack([numQuads, resNum,
                                   np.ones(len(resNum))]).T
                resLSqr = np.linalg.lstsq(resIndep, resTimes)
                mpNum = mtData.T[3]
                mpTimes = mtData.T[4]
                mpIndep = np.vstack([numQuads, mpNum,
                                     np.ones(len(mpNum))]).T
                mpLSqr = np.linalg.lstsq(mpIndep, mpTimes)
                numIncorrect = len(resNum) - np.sum(mtData.T[5])
                adjuster = np.identity(mtData.shape[1])
                adjuster[2][0] = -resLSqr[0][0]
                adjuster[4][0] = -mpLSqr[0][0]
                sizeAdjusted = adjuster.dot(mtData.T).T
                plotFName = test + "." + machine + ".all.png"
                createPlot(mtData, plotFName)
                plotFName = test + "." + machine + ".adjusted.png"
                createPlot(sizeAdjusted, plotFName)
                print("Overall Test Results: " + machine + ", " + test)
                createPlot(testData, plotFName)
                try:
                    pca = sklearn.decomposition.PCA(6).fit(testData)
                    print("PCA Components of " +
                          "[numQuadrics, numResultants, " +
                          "resultantTime_ns, numMPComparisons, " +
                          "mpTime_ns, correct?, counter\n",
                          pca.components_)
                except np.linalg.linalg.LinAlgError:
                    print("Could not compute the PCA of " + plotFName)
                print("Resultant Least Square:", resLSqr)
                print("Increased Precision Least Square:", mpLSqr)
                print()
                

if __name__ == "__main__":
    datum = aggregateData()
    analyzeData(datum)
    

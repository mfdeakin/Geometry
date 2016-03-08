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

def createPlot(data, resLSqr, mpLSqr, fname):
    numQuads = data.T[0]
    resNum = data.T[1]
    resTimes = data.T[2]
    mpNum = data.T[3]
    mpTimes = data.T[4]
    fig = plt.figure()
    fig.suptitle(fname)
    plt.scatter(resNum, resTimes, c = "blue",
               label = "Resultant Method")
    plt.scatter(mpNum, mpTimes, c = "red",
                label = "Increased Precision Method")
    maxValue = max(max(mpNum) * mpLSqr[0][0] + mpLSqr[0][1],
                   max(resNum) * resLSqr[0][0] + resLSqr[0][1],
                   max(mpTimes), max(resTimes))
    plt.axes().set_ylim(bottom = 0,
                        top =  maxValue * 1.0625)
    plt.yticks(np.arange(0, maxValue + 1, maxValue / 10.0))
    plt.axes().yaxis.get_major_formatter().set_powerlimits((0, 3))
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(numQuads, resNum, resTimes, c = "blue",
    #            label = "Resultant Method")
    # ax.scatter(numQuads, mpNum, mpTimes, c = "red",
    #             label = "Increased Precision Method")
    # ax.plot_wireframe(numQuads, resNum,
    #                   numQuads * resLSqr[0][0] +
    #                   resNum * resLSqr[0][1] + resLSqr[0][2])
    # ax.plot_wireframe(numQuads, mpNum,
    #                   numQuads * mpLSqr[0][0] +
    #                   mpNum * mpLSqr[0][1] + mpLSqr[0][2])
    # ax.axes().set_xlim(left = 0, right = max(max(mpNum), max(resNum)) * 1.0625)
    plt.show()
    print("Saving '" + fname + "'")
    fig.savefig(fname, format = "png", dpi = 300)
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
                resNum = testData.T[1]
                resTimes = testData.T[2]
                resLSqr = np.linalg.lstsq(np.vstack([resNum,
                                                     np.ones(len(resNum))]).T,
                                          resTimes)
                mpNum = testData.T[3]
                mpTimes = testData.T[4]
                mpLSqr = np.linalg.lstsq(np.vstack([mpNum,
                                                    np.ones(len(mpNum))]).T,
                                         mpTimes)
                print(resLSqr[0])
                plotFName = (test + "." + machine + "." +
                             str(size) + ".results.png")
                createPlot(testData, resLSqr, mpLSqr, plotFName)
                try:
                    pca = sklearn.decomposition.PCA(6).fit(testData)
                    print("PCA Components of " +
                          "[numQuadrics, numResultants, " +
                          "resultantTime_ns, numMPComparisons, " +
                          "mpTime_ns, correct?, counter\n",
                          pca.components_)
                except np.linalg.linalg.LinAlgError:
                    print("Could not compute the PCA of " + plotFName)
                print("Resultant Least Square Coefficients:", resLSqr[0])
                print("Average Residual:", resLSqr[1][0] / len(resNum))
                print("Increased Precision Coefficients:", mpLSqr[0])
                print("Average Residual:", mpLSqr[1][0] / len(mpNum))

if __name__ == "__main__":
    datum = aggregateData()
    analyzeData(datum)
    

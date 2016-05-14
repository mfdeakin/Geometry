import csv
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sklearn.decomposition

testNameMap = {"hardEllipsoidsSingle": ("Shifted", "Ellipsoids"),
               "packedEll": ("Packed", "Ellipsoids")}

def readData(fname, size):
    fd = open(fname, 'r')
    reader = csv.reader(fd)
    next(reader)
    data = []
    for row in reader:
        if len(row) != 8:
            continue
        number, fpTime, fpCorrect, resNum, resTime, resCorrect, mpNum, mpTime = map(int, row)
        data.append([size,
                     resNum, resTime / 1e6, resCorrect,
                     resNum, fpTime / 1e6, fpCorrect,
                     mpNum, mpTime / 1e6,
                     number])
    return np.array(data)

def validateFName(fname):
    testInfo = fname.split('.')
    if len(testInfo) != 5:
        return None
    testName, testSize, _, testMachine, fileType = fname.split('.')
    if fileType != "csv" or testName == "cylinders":
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
            datum[machine][test] = {}
        if not size in datum[machine][test]:
            datum[machine][test][size] = np.array([])
        testData = datum[machine][test][size]
        data = readData(fname, size)
        newShape = (testData.shape[0] + data.shape[0],
                    data.shape[1])
        datum[machine][test][size] = np.append(testData,
                                               data).reshape(newShape)
    return datum

def createPlot(data, fname, maxPrec,
               resCoeff, resConst,
               fpCoeff, fpConst,
               mpCoeff, mpConst):
    (numQuads,
     resNum, resTimes, resCorrect,
     fpNum, fpTimes, fpCorrect,
     mpNum, mpTimes, _) = data.T
    #plt.suptitle(fname)
    plt.axes().set_color_cycle(["cyan", "orange", "chartreuse"])
    plt.scatter(resNum, resTimes, c = "blue",
               label = "Resultant Method")
    plt.scatter(mpNum, mpTimes, c = "red",
                label = "Increased Precision Method")
    plt.scatter(fpNum, fpTimes, c = "green",
                label = "Machine Precision Method")
    if not (np.isnan(resCoeff) or np.isnan(resConst) or
            np.isnan(mpCoeff) or np.isnan(mpConst)):
        lines = plt.plot([0, maxPrec],
                         [resConst,
                          maxPrec * resCoeff + resConst], '-',
                         [0, maxPrec],
                         [mpConst,
                          maxPrec * mpCoeff + mpConst], '-',
                         [0, maxPrec],
                         [fpConst,
                          maxPrec * fpCoeff + fpConst], '-')
        plt.setp(lines, linewidth=3)
    plt.axes().set_xlim(left = 0,
                        right = maxPrec)
    maxTime = max(max(mpTimes), max(resTimes)) * 1.0625
    plt.axes().set_ylim(bottom = 0,
                        top = maxTime)
    plt.yticks(np.arange(0, maxTime + 1, maxTime / 10.0))
    plt.axes().yaxis.get_major_formatter().set_powerlimits((0, 3))
    print("Saving '" + fname + "'")
    plt.savefig(fname, format = "png", dpi = 300)
    plt.close()

def analyzeData(datum, plotData = False):
    analysis = {}
    for machine in datum:
        analysis[machine] = {}
        for test in datum[machine]:
            mtData = np.array([])
            for size in datum[machine][test]:
                testData = datum[machine][test][size]
                mtShape = (mtData.shape[0] + testData.shape[0],
                           testData.shape[1])
                mtData = np.append(mtData, testData)
                mtData = mtData.reshape(mtShape)
            (numQuads,
             resNum, resTimes, resCorrect,
             fpNum, fpTimes, fpCorrect,
             mpNum, mpTimes, _) = mtData.T
            if max(numQuads) == min(numQuads):
                quadComp = False
                resIndep = np.vstack([resNum,
                                      np.ones(len(resNum))]).T
                fpIndep = np.vstack([fpNum,
                                     np.ones(len(fpNum))]).T
                mpIndep = np.vstack([mpNum,
                                     np.ones(len(mpNum))]).T
            else:
                quadComp = True
                resIndep = np.vstack([numQuads, resNum,
                                  np.ones(len(resNum))]).T
                fpIndep = np.vstack([numQuads, fpNum,
                                     np.ones(len(fpNum))]).T
                mpIndep = np.vstack([numQuads, mpNum,
                                     np.ones(len(mpNum))]).T
            resLSqr = np.linalg.lstsq(resIndep, resTimes)
            fpLSqr = np.linalg.lstsq(fpIndep, fpTimes)
            mpLSqr = np.linalg.lstsq(mpIndep, mpTimes)
            fpDisagree = len(resNum) - np.sum(mtData.T[6])
            resDisagree = len(resNum) - np.sum(mtData.T[3])
            analysis[machine][test] = (resLSqr, resDisagree,
                                       fpLSqr, fpDisagree,
                                       mpLSqr)
            print(machine, test)
            print("Resultant Least Square ((quadrics, precIncreases, constant), (residual), rank, singular values):\n", resLSqr)
            print("Machine Precision Least Square ((quadrics, precIncreases, constant), (residual), rank, singular values):\n", fpLSqr)
            print("Increased Precision Least Square ((quadrics, precIncreases, constant), (residual), rank, singular values):\n", mpLSqr)
            print()
            if plotData:
                plotFName = test + "." + machine + ".all.png"
                createPlot(mtData, plotFName,
                           max(max(resNum), max(mpNum)),
                           resLSqr[0][-2], resLSqr[0][-1],
                           fpLSqr[0][-2], fpLSqr[0][-1],
                           mpLSqr[0][-2], mpLSqr[0][-1])
                if quadComp:
                    adjuster = np.identity(mtData.shape[1])
                    adjuster[2][0] = -resLSqr[0][0]
                    adjuster[5][0] = -fpLSqr[0][0]
                    adjuster[8][0] = -mpLSqr[0][0]
                    sizeAdjusted = adjuster.dot(mtData.T).T
                    plotFName = (test + "." + machine +
                                 ".adjusted.png")
                    createPlot(sizeAdjusted, plotFName,
                               max(max(resNum), max(mpNum)),
                               resLSqr[0][-2], resLSqr[0][-1],
                               fpLSqr[0][-2], fpLSqr[0][-1],
                               mpLSqr[0][-2], mpLSqr[0][-1])
    return analysis

def buildTable(analysis, fname):
    tableOut = open(fname, 'w')
    header = ("\\begin{tabular}{|l|l|ll|lll|l|}\n" +
              "\\hline\n" +
              "Machine & Scene & Method & Disagreements & " +
              "ms/Quadric & ms/Comp & Const ms & " +
              "Residual ($err^2$)\\\\\n" +
              "\\hline\n")
    tableOut.write(header)
    for m in analysis:
        oneMachine = False
        for t in analysis[m]:
            (resLSqr, resDisagree,
             fpLSqr, fpDisagree, mpLSqr) = analysis[m][t]
            rowStr = "\\hline\n"
            if oneMachine == False:
                rowStr += m.capitalize() + " "
                oneMachine = True
            rowStr += ("& {:s} & Increased Prec. & -- & " +
                       "{:s} & {:0.06f} & {:0.06f} & " +
                       "{:0.06f}\\\\\n")
            rSqr = mpLSqr[1][0]
            if len(fpLSqr[0]) == 3:
                quadCoeffStr = "{:0.06f}".format(mpLSqr[0][0])
            else:
                quadCoeffStr = "--"
            rowStr = rowStr.format(testNameMap[t][0],
                                   quadCoeffStr, mpLSqr[0][-2],
                                   mpLSqr[0][-1], rSqr)
            tableOut.write(rowStr)
            
            rowStr = ("& {:s} & Approximate & {:d} & " +
                      "{:s} & {:0.06f} & {:0.06f} & " +
                      "{:0.06f}\\\\\n")
            rSqr = fpLSqr[1][0]
            if len(fpLSqr[0]) == 3:
                quadCoeffStr = "{:0.06f}".format(fpLSqr[0][0])
            else:
                quadCoeffStr = "--"
            rowStr = rowStr.format(testNameMap[t][1],
                                   int(fpDisagree),
                                   quadCoeffStr, fpLSqr[0][-2],
                                   fpLSqr[0][-1], rSqr)
            tableOut.write(rowStr)
            
            rowStr = ("&& Resultant & {:d} & " +
                       "{:s} & {:0.06f} & {:0.06f} & " +
                       "{:0.06f}\\\\\n")
            rSqr = resLSqr[1][0]
            if len(resLSqr[0]) == 3:
                quadCoeffStr = "{:0.06f}".format(resLSqr[0][0])
            else:
                quadCoeffStr = "--"
            rowStr = rowStr.format(int(resDisagree),
                                   quadCoeffStr, resLSqr[0][-2],
                                   resLSqr[0][-1], rSqr)
            tableOut.write(rowStr)
    footer = ("\\hline\n" +
              "\\end{tabular}\n")
    tableOut.write(footer)

if __name__ == "__main__":
    datum = aggregateData()
    analysis = analyzeData(datum, False)
    buildTable(analysis, "comparison_table.tex")
    

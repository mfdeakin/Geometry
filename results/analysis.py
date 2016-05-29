import csv
import os
import sys
import math
import random
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

testNameMap = {"hardEllipsoidsSingle": ["Nested", "Spheres"],
               "hardEllipsoidsSingle_fixedlines": ["Nested", "Spheres", "Fixed"],
               "packedEll_fixedlines": ["Packed", "Spheres"]}

machineNameMap = {"gentoo": ["Gentoo"], "arch": ["Arch"]}

def readData(fname, size):
    fd = open(fname, 'r')
    reader = csv.reader(fd)
    data = []
    for row in reader:
        if len(row) != 19:
            continue
        try:
            (number,
             singleNum, singleTime, singleCorrect,
             doubleNum, doubleTime, doubleCorrect,
             mpNum, mpTime, mpCorrect,
             resNum, resTime) = map(int, row[:12])
            (lDirX, lDirY, lDirZ, lIntX, lIntY, lIntZ) = map(float, row[13:])
        except ValueError:
            continue
        data.append([size,
                     mpNum, mpTime / 1e6, mpCorrect,
                     singleNum, singleTime / 1e6, singleCorrect,
                     doubleNum, doubleTime / 1e6, doubleCorrect,
                     resNum, resTime / 1e6,
                     lDirX, lDirY, lDirZ,
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
               mpCoeff, mpConst,
               singleCoeff, singleConst,
               doubleCoeff, doubleConst,
               resCoeff, resConst):
    (numQuads,
     mpNum, mpTimes, mpCorrect,
     singleNum, singleTimes, singleCorrect,
     doubleNum, doubleTimes, doubleCorrect,
     resNum, resTimes) = data.T[:12]
    plt.axes().set_color_cycle(["cyan", "orange"])
    plt.scatter(resNum, resTimes, c = "blue",
                label = "Resultant Method", marker = "+")
    plt.scatter(mpNum, mpTimes, c = "red",
                label = "Increased Precision Method", marker = "x")
    if not (np.isnan(resCoeff) or np.isnan(resConst) or
            np.isnan(mpCoeff) or np.isnan(mpConst)):
        lines = plt.plot([0, maxPrec],
                         [resConst,
                          maxPrec * resCoeff + resConst], '-',
                         [0, maxPrec],
                         [mpConst,
                          maxPrec * mpCoeff + mpConst], '-',
                         [0, maxPrec])
        plt.setp(lines, linewidth=2)
    plt.axes().set_xlim(left = 0,
                        right = maxPrec)
    maxTime = max(max(mpTimes), max(resTimes)) * 1.0625
    plt.axes().set_ylim(bottom = 0,
                        top = maxTime)
    plt.xlabel("Accurate Comparisons")
    plt.ylabel("Time (ms)")
    plt.yticks(np.arange(0, 400, 50.0))
    plt.axes().yaxis.get_major_formatter().set_powerlimits((0, 3))
    print("Saving '" + fname + "'")
    plt.savefig(fname, format = "png", dpi = 300)
    plt.close()

def cycleVals(vals):
    return np.vstack((np.zeros(len(vals)),
                      vals)).T.reshape(2 * len(vals))

def cycleShuffleCoords(cX, cY, cZ):
    if len(cX) == 0:
        return ([], [], [])
    c = np.vstack((cX, cY, cZ)).T.tolist()
    random.shuffle(c)
    c = np.array(c).T
    cSX = cycleVals(c[0])
    cSY = cycleVals(c[1])
    cSZ = cycleVals(c[2])
    return (cSX, cSY, cSZ)

def removeLabeledVals(vals, labels, keepLabel):
    return [vals[i] for i in range(len(vals))
            if labels[i] == keepLabel]

def plotVectors(fname, correct, dirX, dirY, dirZ):
    print("Plotting Vectors in " + fname)
    maxLinesPlotted = 10000
    f = plt.figure()
    ax = f.add_subplot(111, projection="3d")
    cX = removeLabeledVals(dirX, correct, 1)
    cY = removeLabeledVals(dirY, correct, 1)
    cZ = removeLabeledVals(dirZ, correct, 1)
    cSX, cSY, cSZ = cycleShuffleCoords(cX, cY, cZ)
    ax.plot(cSX[:2 * maxLinesPlotted],
            cSY[:2 * maxLinesPlotted],
            cSZ[:2 * maxLinesPlotted], c = "gray")
    iX = removeLabeledVals(dirX, correct, 0)
    iY = removeLabeledVals(dirY, correct, 0)
    iZ = removeLabeledVals(dirZ, correct, 0)
    iSX, iSY, iSZ = cycleShuffleCoords(iX, iY, iZ)
    ax.plot(iSX[:2 * maxLinesPlotted],
            iSY[:2 * maxLinesPlotted],
            iSZ[:2 * maxLinesPlotted], c = "chartreuse")
    plt.show()
    #plt.savefig(fname, format = "png", dpi = 300)
    plt.close()

def saveData(data, machine, test):
    fname = "dataAggregate.{:s}.{:s}.csv".format(test, machine)
    print("Saving aggregated data in " + fname)
    f = open(fname, 'w')
    writer = csv.writer(f, delimiter = ',')
    writer.writerow(("Number of Quadrics",
                     "Number of Accurate Computations " +
                     "Performed with Increased Precision (IPAC)",
                     "Time of IPAC in ns",
                     "Correct Result from IPAC?",
                     "Number of Accurate Computations " +
                     "Performed with Single Precision (SPAC)",
                     "Time of SPAC in ns",
                     "Correct Result from SPAC?",
                     "Number of Accurate Computations " +
                     "Performed with Double Precision (DPAC)",
                     "Time of DPAC in ns",
                     "Correct Result from DPAC?",
                     "Number of Accurate Computations " +
                     "Performed with Resultants (RAC)",
                     "Time of RAC in ns",
                     "",
                     "Line Direction X Component",
                     "Line Direction Y Component",
                     "Line Direction Z Component"))
    for i in range(len(data)):
        rowData = data[i].T.tolist()
        row = rowData[:12] + [""] + rowData[12:]
        writer.writerow(row)
    return None

def analyzeData(datum, plotData = False):
    analysis = {}
    for machine in datum:
        for test in datum[machine]:
            mtData = np.array([])
            for size in datum[machine][test]:
                testData = datum[machine][test][size]
                mtShape = (mtData.shape[0] + testData.shape[0],
                           testData.shape[1])
                mtData = np.append(mtData, testData)
                mtData = mtData.reshape(mtShape)
            (numQuads,
             mpNum, mpTimes, mpCorrect,
             singleNum, singleTimes, singleCorrect,
             doubleNum, doubleTimes, doubleCorrect,
             resNum, resTimes,
             lDirX, lDirY, lDirZ) = mtData.T[:15]
            if max(numQuads) == min(numQuads):
                quadComp = False
                resIndep = np.vstack([resNum,
                                      np.ones(len(resNum))]).T
                singleIndep = np.vstack([singleNum,
                                     np.ones(len(singleNum))]).T
                doubleIndep = np.vstack([doubleNum,
                                     np.ones(len(doubleNum))]).T
                mpIndep = np.vstack([mpNum,
                                     np.ones(len(mpNum))]).T
            else:
                saveData(mtData, machine, test)
                quadComp = True
                resIndep = np.vstack([numQuads, resNum,
                                  np.ones(len(resNum))]).T
                singleIndep = np.vstack([numQuads, singleNum,
                                         np.ones(len(singleNum))]).T
                doubleIndep = np.vstack([numQuads, doubleNum,
                                         np.ones(len(doubleNum))]).T
                mpIndep = np.vstack([numQuads, mpNum,
                                     np.ones(len(mpNum))]).T
            resLSqr = np.linalg.lstsq(resIndep, resTimes)
            singleLSqr = np.linalg.lstsq(singleIndep, singleTimes)
            doubleLSqr = np.linalg.lstsq(doubleIndep, doubleTimes)
            mpLSqr = np.linalg.lstsq(mpIndep, mpTimes)
            singleDisagree = len(singleNum) - np.sum(singleCorrect)
            doubleDisagree = len(doubleNum) - np.sum(doubleCorrect)
            mpDisagree = len(mpNum) - np.sum(mpCorrect)
            
            if not test in analysis:
                analysis[test] = {}
            analysis[test][machine] = (mpLSqr, mpDisagree,
                                       singleLSqr, singleDisagree,
                                       doubleLSqr, doubleDisagree,
                                       resLSqr)
            print(machine, test)
            print("Single Precision Least Square " +
                  "((quadrics, precIncreases, constant), " +
                  "(residual), rank, singular values):\n",
                  singleLSqr)
            print("Double Precision Least Square " +
                  "((quadrics, precIncreases, constant), " +
                  "(residual), rank, singular values):\n",
                  doubleLSqr)
            print("Increased Precision Least Square "
                  + "((quadrics, precIncreases, constant), " +
                  "(residual), rank, singular values):\n",
                  mpLSqr)
            print("Resultant Least Square " +
                  "((quadrics, precIncreases, constant), " +
                  "(residual), rank, singular values):\n",
                  resLSqr)
            print()
            
            if plotData:
                plotFName = test + "." + machine + ".all.png"
                createPlot(mtData, plotFName,
                           max(max(resNum), max(mpNum)),
                           mpLSqr[0][-2], mpLSqr[0][-1],
                           singleLSqr[0][-2], singleLSqr[0][-1],
                           doubleLSqr[0][-2], doubleLSqr[0][-1],
                           resLSqr[0][-2], resLSqr[0][-1])
                if quadComp:
                    adjuster = np.identity(mtData.shape[1])
                    adjuster[2][0] = -resLSqr[0][0]
                    adjuster[5][0] = -singleLSqr[0][0]
                    adjuster[8][0] = -mpLSqr[0][0]
                    sizeAdjusted = adjuster.dot(mtData.T).T
                    plotFName = (test + "." + machine +
                                 ".adjusted.png")
                    createPlot(sizeAdjusted, plotFName,
                               max(max(resNum), max(mpNum)),
                               mpLSqr[0][-2], mpLSqr[0][-1],
                               singleLSqr[0][-2], singleLSqr[0][-1],
                               doubleLSqr[0][-2], doubleLSqr[0][-1],
                               resLSqr[0][-2], resLSqr[0][-1])
                for type in [("mp", mpCorrect),
                             ("single", singleCorrect),
                             ("double", doubleCorrect)]:
                    plotFName = (test + "." + machine + "." +
                                 type[0] + ".vectors.png")
                    plotVectors(plotFName, type[1],
                                lDirX, lDirY, lDirZ)
    return analysis

def formatFloat(fpval, sigDigits, shiftRight, negativeSpace = True):
    digitsLeft = int(np.floor(np.log10(abs(fpval))))
    printable = round(fpval, sigDigits - 1 - digitsLeft)
    if printable == int(printable):
        printable = int(printable)
    digitsPrec = sigDigits - digitsLeft - 1
    printed = ("{:." + str(digitsPrec) + "f}").format(printable)
    if printed[0] != '-' and negativeSpace:
        printed = "\\hphantom{-}" + printed
    if digitsLeft >= 0:
        shiftRight -= digitsLeft + 1
    if shiftRight > 0:
        printed = "\\hphantom{" + "0" * shiftRight + "}" + printed
    return printed

def formatInt(ival, shiftRight, negativeSpace = False):
    printed = str(int(ival))
    digits = len(printed)
    spacer = ""
    if printed[0] == '-':
        digits -= 1
    elif negativeSpace:
        spacer = "\\hphantom{-}"
    shiftRight -= digits
    if shiftRight > 0:
        spacer = "\\hphantom{" + "0" * shiftRight + "}" + spacer
    return spacer + printed

def printRow(test, testIdx, machine, machIdx, method,
             lsqrVals, errors = None):
    if testIdx < len(testNameMap[test]):
        testPrint = testNameMap[test][testIdx]
    else:
        testPrint = ""
    if machIdx < len(machineNameMap[machine]):
        machPrint = machineNameMap[machine][machIdx]
    else:
        machPrint = ""
    if errors == None:
        errorPrint = "\\hphantom{---}---"
    else:
        errorPrint = formatInt(errors, 4)
    if len(lsqrVals[0]) == 3:
        mspquad = formatFloat(lsqrVals[0][0], 3, 0, False)
    else:
        mspquad = "--"
    mspcomp = formatFloat(lsqrVals[0][-2], 3, 0)
    msconst = formatFloat(lsqrVals[0][-1], 3, 0)
    residual = formatFloat(lsqrVals[1][0], 6, 6)
    fmtString = ("{:s} & {:s} & {:s} & {:s} & {:s} & {:s} & " +
                 "{:s} & {:s}\\\\\n")
    rowStr = fmtString.format(testPrint, machPrint, method,
                              errorPrint, mspquad, mspcomp, msconst, residual)
    return rowStr

def buildTable(analysis, fname):
    tableOut = open(fname, 'w')
    header = ("\\begin{tabular}{|l|l|ll|lll|l|}\n" +
              "\\hline\n" +
              "Scene & Machine & Method & Errors & " +
              "ms/Quadric & ms/Comp & Const ms & " +
              "$\\sum$ Residual ($\\text{ms}^2$)\\\\\n" +
              "\\hhline{|=|=|==|===|=|}\n")
    tableOut.write(header)
    printHLine = False
    for t in analysis:
        if printHLine:
            tableOut.write("\\hhline{|-|-|--|---|-|}\n")
        printHLine = False
        testPrinted = 0
        for m in analysis[t]:
            if printHLine:
                tableOut.write("\\hhline{|~|-|--|---|-|}\n")
            printHLine = True
            (mpLSqr, mpDisagree,
             singleLSqr, singleDisagree,
             doubleLSqr, doubleDisagree,
             resLSqr) = analysis[t][m]
            rowStr = printRow(t, testPrinted, m, 0,
                              "Single Prec", singleLSqr,
                              singleDisagree)
            tableOut.write(rowStr)
            testPrinted += 1
            
            rowStr = printRow(t, testPrinted, m, 1,
                              "Double Prec", doubleLSqr,
                              doubleDisagree)
            tableOut.write(rowStr)
            testPrinted += 1
            
            rowStr = printRow(t, testPrinted, m, 2,
                              "Increased Prec.", mpLSqr,
                              mpDisagree)
            tableOut.write(rowStr)
            testPrinted += 1
            
            rowStr = printRow(t, testPrinted, m, 3,
                              "Resultant", resLSqr,
                              None)
            tableOut.write(rowStr)
            testPrinted += 1
    footer = ("\\hline\n" +
              "\\end{tabular}\n")
    tableOut.write(footer)
    return None

if __name__ == "__main__":
    datum = aggregateData()
    genPlots = True
    if len(sys.argv) > 1:
        genPlots = not (sys.argv[1].lower() == "false")
    analysis = analyzeData(datum, genPlots)
    buildTable(analysis, "comparison_table.tex")
    

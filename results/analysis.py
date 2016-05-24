import csv
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sklearn.decomposition

testNameMap = {"hardEllipsoidsSingle": ["Nested", "Spheres"],
               "packedEll": ["Packed", "Spheres"]}

machineNameMap = {"gentoo": ["Gentoo"], "arch": ["Arch"]}

def readData(fname, size):
    fd = open(fname, 'r')
    reader = csv.reader(fd)
    next(reader)
    data = []
    for row in reader:
        if len(row) != 8:
            continue
        number, fpTime, fpCorrect, mpNum, mpTime, mpCorrect, resNum, resTime = map(int, row)
        data.append([size,
                     mpNum, mpTime / 1e6, mpCorrect,
                     resNum, fpTime / 1e6, fpCorrect,
                     resNum, resTime / 1e6,
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
               mpCoeff, mpConst,
               fpCoeff, fpConst):
    (numQuads,
     fpNum, fpTimes, mpCorrect,
     mpNum, mpTimes, fpCorrect,
     resNum, resTimes, _) = data.T
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
    plt.yticks(np.arange(0, 750, 50.0))
    plt.axes().yaxis.get_major_formatter().set_powerlimits((0, 3))
    print("Saving '" + fname + "'")
    plt.savefig(fname, format = "png", dpi = 300)
    plt.close()

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
             mpNum, mpTimes, fpCorrect,
             fpNum, fpTimes, mpCorrect,
             resNum, resTimes, _) = mtData.T
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
            fpDisagree = len(resNum) - np.sum(fpCorrect)
            mpDisagree = len(resNum) - np.sum(mpCorrect)
            if not test in analysis:
                analysis[test] = {}
            analysis[test][machine] = (mpLSqr, mpDisagree,
                                       fpLSqr, fpDisagree,
                                       resLSqr)
            print(machine, test)
            print("Machine Precision Least Square ((quadrics, precIncreases, constant), (residual), rank, singular values):\n", fpLSqr)
            print("Increased Precision Least Square ((quadrics, precIncreases, constant), (residual), rank, singular values):\n", mpLSqr)
            print("Resultant Least Square ((quadrics, precIncreases, constant), (residual), rank, singular values):\n", resLSqr)
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
            (fpLSqr, fpDisagree,
             mpLSqr, mpDisagree,
             resLSqr) = analysis[t][m]
            rowStr = printRow(t, testPrinted, m, 0,
                              "Approximate", fpLSqr,
                              fpDisagree)
            tableOut.write(rowStr)
            testPrinted += 1
            
            rowStr = printRow(t, testPrinted, m, 1,
                              "Increased Prec.", mpLSqr,
                              mpDisagree)
            tableOut.write(rowStr)
            testPrinted += 1
            
            rowStr = printRow(t, testPrinted, m, 2,
                              "Resultant", resLSqr,
                              None)
            tableOut.write(rowStr)
            testPrinted += 1
    footer = ("\\hline\n" +
              "\\end{tabular}\n")
    tableOut.write(footer)

if __name__ == "__main__":
    datum = aggregateData()
    genPlots = False
    if len(sys.argv) > 1:
        genPlots = (sys.argv[1].lower() == "true")
    analysis = analyzeData(datum, genPlots)
    buildTable(analysis, "comparison_table.tex")
    

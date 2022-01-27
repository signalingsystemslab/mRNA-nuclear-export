#!/usr/bin/env python3
class CheckArgs():  #class that checks arguments and ultimately returns a validated set of arguments to the main program
    
    def __init__(self):
        import argparse
        import os
        parser = argparse.ArgumentParser()
        parser.add_argument("-p", "--pe1", help = "First paired end reads", required = True)
        parser.add_argument("-q", "--pe2", help = "Second paired end reads")
        parser.add_argument("-b", "--barcodes", help = "Barcode reads", required = True)
        parser.add_argument("-c", "--barcodes2", help = "Second barcode reads")
        parser.add_argument("-d", "--outputDirectory", help = "Directory for demultiplexed partial outputs")
        parser.add_argument("-l", "--validBarcodeList", help = "Comma-separated list of valid barcodes", required = True)
        parser.add_argument("-a", "--useSampleNames", help = "Use sample names instead of barcode for output files", action = "store_true")
        parser.add_argument("-v", "--verbose", help = "Run in verbose mode", action = "store_true")
        parser.add_argument("-z", "--gzip", help = "Force gzip output", action = "store_true")
        parser.add_argument("-g", "--noGzip", help = "Disallow gzip output", action = "store_true")
        parser.add_argument("-s", "--qualityScore33Conversion", help = "Convert base 64 quality scores to base 33", action = "store_true")
        rawArgs = parser.parse_args()
        pe1 = rawArgs.pe1
        if not os.path.isfile(pe1):
            raise FileNotFoundError("Paired end 1 file not found at %s" %pe1)
        self.pe1 = pe1
        pe2 = rawArgs.pe2
        if pe2 and not os.path.isfile(pe2):
            raise FileNotFoundError("Paired end 2 file not found at %s" %pe2)
        self.pe2 = pe2
        barcodes = rawArgs.barcodes
        if not os.path.isfile(barcodes):
            raise FileNotFoundError("Barcode file not found at %s" %barcodes)
        self.barcodes = barcodes
        barcodes2 = rawArgs.barcodes2
        if barcodes2 and not os.path.isfile(barcodes2):
            raise FileNotFoundError("Second barcode file not found at %s" %barcodes2)
        self.barcodes2 = barcodes2
        self.outputDirectory = rawArgs.outputDirectory
        self.barcodeTable = processValidBarcodeList(rawArgs.validBarcodeList)
        self.useSampleNames = rawArgs.useSampleNames
        self.verbose = rawArgs.verbose
        self.gzip = rawArgs.gzip
        self.noGzip = rawArgs.noGzip
        if self.gzip and self.noGzip:
            raise RuntimeError("Force gzip and disallow gzip cannot both be set to true.")
        self.qualityScore33Conversion = rawArgs.qualityScore33Conversion
        
def processValidBarcodeList(barcodeList):
    barcodeList = barcodeList.split(",")
    try:
        barcodeLength = int(barcodeList[0])
    except ValueError:
        newBarcodeList = []
        barcodeLength = len(barcodeList[0])
        for barcode in barcodeList:
            validatedBarcode = validateBarcodeSeq(barcode, barcodeLength)
            if not validatedBarcode:
                raise RuntimeError("Invalid entry found in barcodes: %s\n%s" %(barcode, barcodeList))
            newBarcodeList.append([validatedBarcode, validatedBarcode])
        return newBarcodeList
    else:
        del barcodeList[0]
        newBarcodeList = []
        for item in barcodeList:
            barcode = item[:barcodeLength]
            name = item[barcodeLength:]
            validatedBarcode = validateBarcodeSeq(barcode, barcodeLength)
            if not validatedBarcode or not name:
                raise RuntimeError("Invalid entry found in barcodes: %s\n%s" %(barcode, barcodeList))
            newBarcodeList.append([validatedBarcode, name])
        return newBarcodeList
    
def validateBarcodeSeq(seq, lengthRequirement = False):
    seq = seq.upper()
    for letter in seq:
        if not letter in "ATGC":
            return False
    if lengthRequirement and not len(seq) == lengthRequirement:
        return False
    return seq
    
        
class QseqLine(object):

    def __init__(self, rawLine, barcodeSeq = False, readNumber = 1, runIDConversion = False, flowCellID = "FC000", skipFailedFilterLines = True): ####
        #### def __init__(self, rawLine, barcodeSeq = False, runIDConversion = False, flowCellID = "FC000", skipFailedFilterLines = True):
        rawLine = rawLine.strip()
        #print(rawLine)
        if rawLine.endswith("0") and skipFailedFilterLines:
            self.passedFilter = False  #only setting this value and skipping all others to save time.  This will be the only attribute of the object.
        else:
            self.barcodeSeq = barcodeSeq
            if self.barcodeSeq:
                self.barcodeSeq = self.barcodeSeq.strip()
            self.rawLine = rawLine
            self.flowCellID = flowCellID
            self.readNumber = readNumber ####
            self.rawLineList = rawLine.split("\t")
            self.machine = self.rawLineList[0]
            try:
                self.run = int(self.rawLineList[1])
            except ValueError:
                if runIDConversion:
                    self.run = runIDConversion
                else:
                    self.run = self.runIDConversion(self.rawLineList[1])
                if not self.run:
                    raise RuntimeError("Unexpected format for run ID on line:\n %s\nRun ID %s was not an integer and could not be converted with existing conversion scheme." %(rawLine, self.run))
                else:
                    self.flowCellID = self.rawLineList[1]
            self.lane = int(self.rawLineList[2])
            self.tile = int(self.rawLineList[3])
            self.xCoord = float(self.rawLineList[4])
            self.yCoord = float(self.rawLineList[5])
            self.index = int(self.rawLineList[6])
            self.readNumber = int(self.readNumber) ####
            #### self.readNumber = int(self.rawLineList[7])
            self.sequence = self.rawLineList[8]
            self.quality = self.rawLineList[9]
            self.filter = int(self.rawLineList[10])
            self.passedFilter = bool(self.filter)
            self.fastqConversionString = None
            self.calledBarcode = None  #holder for a corrected barcode
    
    def __str__(self):
        return self.rawLine
    
    def runIDConversion(self, runID):
        def bscrc(runID):
            import re
            result = re.match("(\D+)(\d+)", runID)
            if not result:
                return False
            letters = result.group(1).upper()
            numbers = result.group(2)
            if len(letters) > 3:
                return False
            machine = str(ord(letters[0]))
            try:
                side = letters[1]
            except IndexError:
                return False
            try:
                side = {"A":"0", "B":"1"}[side]  #using an ephemeral dictionary
            except KeyError:
                return False
            try:
                if letters[2] == "P":
                    pairCode = "1"
                else:
                    pairCode = "0"
            except IndexError:
                pairCode = "0"
            return int(machine + side + pairCode + numbers)
        conversion = bscrc(runID)
        if conversion:
            return conversion
        #room for elif statements for future conversion schemes if needed
        else:
            return False
        
    def convertQuality64to33(self):
        test = True
        self.quality = "".join([chr(ord(character)-31) for character in self.quality])
        for character in self.quality:
            if ord(character) < 33:
                raise InvalidQualityConversion("Invalid quality score conversion returned for line %s" %self.rawLine)
        
    def fastqConversion(self, convertQualityScore):
        if not self.fastqConversionString:
            if convertQualityScore:
                self.convertQuality64to33()
            self.fastqConversionString = self.createFastqLines()
        return self.fastqConversionString  #quick return if we have already made this value
    
    def coordinateConversion(self, coordinate):
        return int((coordinate * 10) + 1000.5) 
    
    def createFastqLines(self):
        cluster = []
        cluster.append("@" + self.machine)
        cluster.append(self.run)
        cluster.append(self.flowCellID)
        cluster.append(self.lane)
        cluster.append(self.tile)
        cluster.append(self.coordinateConversion(self.xCoord))
        cluster.append(self.coordinateConversion(self.yCoord))
        cluster = [str(item) for item in cluster]
        firstCluster = ":".join(cluster)
        cluster = []
        cluster.append(self.readNumber)
        cluster.append({0:"Y", 1:"N"}[self.filter]) ####
        #### cluster.append({1:"Y", 0:"N"}[self.filter])
        cluster.append("0")
        cluster.append(self.barcodeSeq)
        cluster = [str(item) for item in cluster]
        secondCluster = ":".join(cluster)
        lines = []
        lines.append(firstCluster + " " + secondCluster)
        lines.append(self.sequence.replace(".","N"))
        lines.append("+")
        lines.append(self.quality)
        return "\n".join(lines)
    
class InvalidQualityConversion(Exception):
    
    def __init__(self, message):
        self.message = message

def matchingReadCounts(barcodes, sequence):
    if not len(barcodes) == len(sequence):
        return False
    return True

def stripList(inputList):
    if inputList: # check that input list is not empty
        while not inputList[0].strip():
            del inputList[0]
        while not inputList[-1].strip():
            del inputList[-1]
        
def getFirstLineSeqLength(qseqFile):
    if qseqFile.endswith(".gz"):
        import gzip
        file = gzip.open(qseqFile, 'rt')
    else:
        file = open(qseqFile, 'r')
    line = False
    eof = False
    while not line and not eof:
        line = file.readline()
        eof = not line
        line = line.strip()
    line  = QseqLine(line, skipFailedFilterLines = False)
    return len(line.sequence)
     
def initialBarcodeFileProcessing(barcodeDataLines, validBarcodeTable, verbose = False, removeFailedFilter = True):
    readBarcodeCounts = {}
    calledBarcodeCounts = {}
    validBarcodeList = [item[0] for item in validBarcodeTable]
    emptyBarcodeCount = 0
    unmatchedBarcodeCount = 0
    qualityBase = 0
    if verbose:
        print("Analyzing first barcode line for run data.")
    firstLine = QseqLine(barcodeDataLines[0], skipFailedFilterLines = False)
    barcodeLength = len(firstLine.sequence)
    flowCellID = firstLine.flowCellID
    runIDConversion = firstLine.run
    tile = firstLine.tile
    lane = firstLine.lane
    if not "N" * barcodeLength in validBarcodeList:
        validBarcodeList.append("N" * barcodeLength)
    for i in range(0, len(barcodeDataLines)):
        if verbose and i % 10000 == 0:
            print("Parsing barcode line %s" %i, end = '\r')
        barcodeDataLines[i] = QseqLine(barcodeDataLines[i], runIDConversion = runIDConversion, flowCellID = flowCellID, skipFailedFilterLines = False)
    if verbose:
        print("Parsed barcode line %s    " %len(barcodeDataLines))
    if verbose:
        print("Checking quality score encoding and calling barcodes.")
    for i in range(0, len(barcodeDataLines)):
        if verbose and i % 10000 == 0:
            print("Calling barcode line %s" %(i), end = "\r")
        readBarcode = barcodeDataLines[i].sequence.replace(".", "N")
        if not qualityBase and "N" in readBarcode:
            locationOfN = readBarcode.index("N")
            qualityOfN = barcodeDataLines[i].quality[locationOfN]
            if qualityOfN == "@":
                if verbose:
                    print("Base 64 quality score encoding found.                    ")
                qualityBase = 64
            if qualityOfN == "!":
                if verbose:
                    print("Base 33 quality score encoding found.                    ")
                qualityBase = 33
        if readBarcode == "N" * barcodeLength:
            emptyBarcodeCount += 1
        calledBarcode = callBarcode(readBarcode, validBarcodeList, -(-barcodeLength//8))  #will allow for up to one mismatch per 8 bases of barcode length.  Can make this adjustable later if needed.
        if not calledBarcode:
            calledBarcode = "Unmatched"
            unmatchedBarcodeCount += 1
        if not readBarcode in readBarcodeCounts:
            readBarcodeCounts[readBarcode] = 1
        else:
            readBarcodeCounts[readBarcode] += 1
        if not calledBarcode in calledBarcodeCounts:
            calledBarcodeCounts[calledBarcode] = 1
        else:
            calledBarcodeCounts[calledBarcode] += 1
        barcodeDataLines[i].calledBarcode = calledBarcode
    if verbose:
        print("Called barcode line %s      " %len(barcodeDataLines))
        print("Lines before filtering: %s" %len(barcodeDataLines))
    barcodeList = []
    for line in barcodeDataLines:
        if removeFailedFilter and not line.passedFilter:
            continue
        barcodeList.append(line.calledBarcode)
    if verbose:
        print("Lines after filtering: %s" %len(barcodeList))
    return (barcodeList, barcodeLength, tile, lane, runIDConversion, flowCellID, readBarcodeCounts, calledBarcodeCounts, emptyBarcodeCount, unmatchedBarcodeCount, qualityBase)
    
def readInFullQseqFile(fileHandle, skipFailedFiltering = True):
    lineList = []
    line = fileHandle.readline()
    while line:
        if skipFailedFiltering:
            if line.strip().endswith("0"):
                line = fileHandle.readline()
                continue
        if line.strip():
            lineList.append(line.strip())
        line = fileHandle.readline()
    return lineList

def isUniformLength(inputList):
    firstLength = len(inputList[0])
    for element in inputList:
        if not len(element) == firstLength:
            return False
    return True

def convertDualBarcodes(validBarcodeTable, barcode1Length, barcode2Length):
    barcode1Table = []
    barcode2Table = []
    for line in validBarcodeTable:
        combinedBarcode, sample = line
        if not len(combinedBarcode) == barcode1Length + barcode2Length:
            raise RuntimeError("Combined barcode %s is not made up of a %s base sequence and a %s base sequence as expected from the individual barcode lengths." %(combinedBarcode, barcode1Length, barcode2Length))
        barcode1 = combinedBarcode[:barcode1Length]
        barcode2 = combinedBarcode[barcode1Length:]
        barcode1Table.append([barcode1, sample])
        barcode2Table.append([barcode2, sample])
    return (barcode1Table, barcode2Table)

def validateDualBarcodes(validBarcodeTable, calledCombinedBarcodeList):   #can do this in place instead of as a return
    validBarcodeList = [item[0] for item in validBarcodeTable]
    for i in range(0, len(calledCombinedBarcodeList)):
        if "Unmatched" in calledCombinedBarcodeList[i]:  #marks out any barcode where one of the members could not be matched
            calledCombinedBarcodeList[i] = "Unmatched"
            continue
        if not calledCombinedBarcodeList[i] in validBarcodeList:  #marks out any invalid barcode combinations
            calledCombinedBarcodeList[i] = "Unmatched"
            
def combineBarcodeLists(first, second):
    assert len(first) == len(second), "First and second barcode lists are not the same number of elements (%s and %s, respectively)." %(len(first), len(second))
    combined = []
    for i in range(0, len(first)):
        combined.append(first[i] + second[i])
    return combined

def callBarcode(readBarcode, validBarcodeList, maxDifference = 1):
    if readBarcode in validBarcodeList:
        return readBarcode
    if not len(readBarcode) == len(validBarcodeList[0]):
        raise RuntimeError("Read barcode %s had different length than the valid barcodes\n%s" %(readBarcode, validBarcodeList))
    possibleBarcodes = []
    for validBarcode in validBarcodeList:
        if calculateMismatches(readBarcode, validBarcode, maxDifference) <= maxDifference:
            possibleBarcodes.append(validBarcode)
    if len(possibleBarcodes) == 1:
        return possibleBarcodes[0]
    else:
        return False

def calculateMismatches(string1, string2, maxDifference = False):
    if not len(string1) == len(string2):
        raise RuntimeError("Cannot calculate hamming distances on strings of different lengths")
    mismatches = 0
    for i in range(0,len(string1)): 
        if not string1[i] == string2[i]:
            mismatches += 1
            if maxDifference:
                if mismatches > maxDifference:
                    return len(string1) + 1
    return mismatches

def convertData(data, barcodeList, readNumber, flowCellID, runIDConversion, verbose = False): ####
    #### def convertData(data, barcodeList, flowCellID, runIDConversion, verbose = False):
    for i in range(0, len(barcodeList)):
        if verbose and i % 10000 == 0:
            print("Converting sequence line %s" %(i), end = "\r")
        data[i] = QseqLine(data[i], barcodeList[i], readNumber, runIDConversion, flowCellID) ####
        #### data[i] = QseqLine(data[i], barcodeList[i], runIDConversion, flowCellID)
    if verbose:
        print("Converted sequence line %s" %len(barcodeList))
        
def demultiplex(data, qualityConvert, verbose = False):
    for i in range(0, len(data)):
        if verbose and i % 10000 == 0:
            print("Demultiplexing line %s" %(i), end = "\r")
        if i == 0:
            barcode = data[i].barcodeSeq
            seqHolder = data[i].fastqConversion(qualityConvert)
            data[i] = {barcode:[seqHolder]}
        else:
            if not data[i].barcodeSeq in data[0]:
                data[0][data[i].barcodeSeq] = []            
            data[0][data[i].barcodeSeq].append(data[i].fastqConversion(qualityConvert))
            data[i] = []
    if verbose:
        print("Demultiplexed line %s    " %len(data))
    if data:
        del data[1:]
        for barcode in list(data[0].keys()):
            if verbose:
                print("Joining data for barcode: %s          " %barcode, end = '\r')
            data[0][barcode] = "\n".join(data[0][barcode])
    if verbose:
        print()

def writeData(data, outputDir, tile, lane, readNumber, barcodeTable, useSampleName, gzipped = False, verbose = False, skipUnmatched = True):
    import os
    if gzipped:
        import gzip
    if useSampleName:
        nameHash = {}
        for item in barcodeTable:
            nameHash[item[0]] = item[1]
    
    if not outputDir:
        outputDir = os.getcwd()
    if not os.path.isdir(outputDir):
        try:
            os.mkdir(outputDir)
        except FileExistsError:
            pass
        if not os.path.isdir(outputDir):
            raise RuntimeError("Unable to find or create the needed output directory.")
    readNumber = str(readNumber)
    tile = str(tile)
    lane = str(lane)
    
    ### create filename first and creat file if it doesn't exist:
    for item in barcodeTable:
        if useSampleName: ###
            fileName = ".".join([nameHash[item[0]], lane, tile, readNumber, "part"]) ###
        else: ### 
            fileName = ".".join([item[0], lane, tile, readNumber, "part"]) ###
        if gzipped: ###
            if not os.path.exists(outputDir + os.sep + fileName + "gz" ): ###
                file = gzip.open(outputDir + os.sep + fileName + "gz", "wt") ###
                file.close() ###
        else: ###
            if not os.path.exists(outputDir + os.sep + fileName): ###
                file = open(outputDir + os.sep + fileName, "w") ###
                file.close() ### 
    
    if data:
        data = data[0]
        for barcode in list(data.keys()):
            if useSampleName: ###
                fileName = ".".join([nameHash[barcode], lane, tile, readNumber, "part"]) ###
            else: ### 
                fileName = ".".join([barcode, lane, tile, readNumber, "part"]) ###
            ### open and write only if barcode fund
            if skipUnmatched:
                if barcode == "N" * len(barcode) or barcode == "Unmatched":
                    continue
            ### if useSampleName:
                ### fileName = ".".join([nameHash[barcode], lane, tile, readNumber, "part"])
            ### else:
                ### fileName = ".".join([barcode, lane, tile, readNumber, "part"])
            if gzipped:
                file = gzip.open(outputDir + os.sep + fileName + ".gz", 'wt')
            else:
                file = open(outputDir + os.sep + fileName, 'w')
            if verbose:
                print("Writing %s" %fileName, end = "\r")
            print(data[barcode], file = file)
            file.close()
        if verbose:
            print()

def dualBarcodes(args):
    import datetime
    starttime = datetime.datetime.now()
    verbose = args.verbose
    barcode1Length = getFirstLineSeqLength(args.barcodes)
    barcode2Length = getFirstLineSeqLength(args.barcodes2)
    barcode1Table, barcode2Table = convertDualBarcodes(args.barcodeTable, barcode1Length, barcode2Length)
    if args.barcodes.endswith(".gz"):
        import gzip
        barcodeFile = gzip.open(args.barcodes, 'rt')
        gzipFile = True
    else:
        barcodeFile = open(args.barcodes, 'r')
        gzipFile = False
    if verbose:
        print("Loading first set of barcodes from file...", end = "", flush = True)
    barcodes = readInFullQseqFile(barcodeFile, skipFailedFiltering = False)
    if verbose:
        print("DONE")
    barcodeFile.close()
    stripList(barcodes)
    barcodeList1, barcodeLength1, tile, lane, runIDConversion, flowCellID, readBarcodeCounts, calledBarcodeCounts, emptyBarcodeCount, unmatchedBarcodeCount, qualityBase = initialBarcodeFileProcessing(barcodes, barcode1Table, verbose)
    import operator
    print("Most common read barcodes:")
    for line in (sorted(readBarcodeCounts.items(), key=operator.itemgetter(1))[-15:]):
        print(line)
    print("Most common called barcodes:")
    for line in (sorted(calledBarcodeCounts.items(), key=operator.itemgetter(1))):
        print(line)
    if args.barcodes2.endswith(".gz"):
        import gzip
        barcodeFile = gzip.open(args.barcodes2, 'rt')
        gzipFile = True
    else:
        barcodeFile = open(args.barcodes2, 'r')
        gzipFile = False
    if verbose:
        print("Loading second set of barcodes from file...", end = "", flush = True)
    barcodes = readInFullQseqFile(barcodeFile, skipFailedFiltering = False)
    if verbose:
        print("DONE")
    barcodeFile.close()
    stripList(barcodes)
    barcodeList2, barcodeLength2, tile, lane, runIDConversion, flowCellID, readBarcodeCounts, calledBarcodeCounts, emptyBarcodeCount, unmatchedBarcodeCount, qualityBase = initialBarcodeFileProcessing(barcodes, barcode2Table, verbose)
    del barcodes
    if qualityBase == 64 and args.qualityScore33Conversion:
        qualityConvert = True
    else:
        qualityConvert = False
    barcodeList = combineBarcodeLists(barcodeList1, barcodeList2)
    del barcodeList1
    del barcodeList2
    validateDualBarcodes(args.barcodeTable, barcodeList)
    print("Most common read barcodes:")
    for line in (sorted(readBarcodeCounts.items(), key=operator.itemgetter(1))[-15:]):
        print(line)
    print("Most common called barcodes:")
    for line in (sorted(calledBarcodeCounts.items(), key=operator.itemgetter(1))):
        print(line)
    if args.pe1.endswith(".gz"):
        import gzip
        pe1File = gzip.open(args.pe1, 'rt')
        gzipFile = True
    else:
        pe1File = open(args.pe1, 'r')
        gzipFile = False
    if verbose:
        print("Loading paired end 1 from file...", end = "", flush = True)
    pe = readInFullQseqFile(pe1File)
    if verbose:
        print("DONE")
    pe1File.close()
    stripList(pe)
    if not matchingReadCounts(barcodeList, pe):
        raise RuntimeError("Lengths of barcode and pe1 files do not match.  Please determine why these files have different line counts.\nBarcode: %s\nPE1: %s" %(len(barcodeList), len(pe)))
    convertData(pe, barcodeList, 1, flowCellID, runIDConversion, verbose) ####
    #### convertData(pe, barcodeList, flowCellID, runIDConversion, verbose)
    demultiplex(pe, qualityConvert, verbose)
    writeData(pe[0], args.outputDirectory, tile, lane, 1, args.barcodeTable, args.useSampleNames, (gzipFile or args.gzip) and not args.noGzip, verbose)
    if args.pe2:
        if args.pe2.endswith(".gz"):
            import gzip
            pe2File = gzip.open(args.pe2, 'rt')
            gzipFile = True
        else:
            pe2File = open(args.pe2, 'r')
            gzipFile = False
        if verbose:
            print("Loading paired end 2 from file...", end = "", flush = True)
        pe = readInFullQseqFile(pe2File)
        if verbose:
            print("DONE")
        pe2File.close()
        stripList(pe)
        if not matchingReadCounts(barcodeList, pe):
            raise RuntimeError("Lengths of barcode and pe2 files do not match.  Please determine why these files have different line counts.\nBarcode: %s\nPE1: %s" %(len(barcodeList), len(pe)))
        convertData(pe, barcodeList, 2, flowCellID, runIDConversion, verbose) ####
        #### convertData(pe, barcodeList, flowCellID, runIDConversion, verbose)
        demultiplex(pe, qualityConvert, verbose)
        writeData(pe[0], args.outputDirectory, tile, lane, 2, args.barcodeTable, args.useSampleNames, (gzipFile or args.gzip) and not args.noGzip, verbose)
    print("Runtime: %s" %(datetime.datetime.now() - starttime))

def main(args):
    import datetime
    starttime = datetime.datetime.now()
    verbose = args.verbose
    if args.barcodes.endswith(".gz"):
        import gzip
        barcodeFile = gzip.open(args.barcodes, 'rt')
        gzipFile = True
    else:
        barcodeFile = open(args.barcodes, 'r')
        gzipFile = False
    if verbose:
        print("Loading barcodes from file...", end = "", flush = True)
    barcodes = readInFullQseqFile(barcodeFile, skipFailedFiltering = False)
    if verbose:
        print("DONE")
    barcodeFile.close()
    stripList(barcodes)
    barcodeList, barcodeLength, tile, lane, runIDConversion, flowCellID, readBarcodeCounts, calledBarcodeCounts, emptyBarcodeCount, unmatchedBarcodeCount, qualityBase = initialBarcodeFileProcessing(barcodes, args.barcodeTable, verbose)
    del barcodes
    if qualityBase == 64 and args.qualityScore33Conversion:
        qualityConvert = True
    else:
        qualityConvert = False
    if verbose:
        import operator
        print("Most common read barcodes:")
        for line in (sorted(readBarcodeCounts.items(), key=operator.itemgetter(1))[-15:]):
            print(line)
        print("Most common called barcodes:")
        for line in (sorted(calledBarcodeCounts.items(), key=operator.itemgetter(1))):
            print(line)
    if args.pe1.endswith(".gz"):
        import gzip
        pe1File = gzip.open(args.pe1, 'rt')
        gzipFile = True
    else:
        pe1File = open(args.pe1, 'r')
        gzipFile = False
    if verbose:
        print("Loading paired end 1 from file...", end = "", flush = True)
    pe = readInFullQseqFile(pe1File)
    if verbose:
        print("DONE")
    pe1File.close()
    stripList(pe)
    if not matchingReadCounts(barcodeList, pe):
        raise RuntimeError("Lengths of barcode and pe1 files do not match.  Please determine why these files have different line counts.\nBarcode: %s\nPE1: %s" %(len(barcodeList), len(pe)))
    convertData(pe, barcodeList, 1, flowCellID, runIDConversion, verbose) ####
    #### convertData(pe, barcodeList, flowCellID, runIDConversion, verbose)
    demultiplex(pe, qualityConvert, verbose)
    writeData(pe, args.outputDirectory, tile, lane, 1, args.barcodeTable, args.useSampleNames, (gzipFile or args.gzip) and not args.noGzip, verbose)
    if args.pe2:
        if args.pe2.endswith(".gz"):
            import gzip
            pe2File = gzip.open(args.pe2, 'rt')
            gzipFile = True
        else:
            pe2File = open(args.pe2, 'r')
            gzipFile = False
        if verbose:
            print("Loading paired end 2 from file...", end = "", flush = True)
        pe = readInFullQseqFile(pe2File)
        if verbose:
            print("DONE")
        pe2File.close()
        stripList(pe)
        if not matchingReadCounts(barcodeList, pe):
            raise RuntimeError("Lengths of barcode and pe2 files do not match.  Please determine why these files have different line counts.\nBarcode: %s\nPE1: %s" %(len(barcodeList), len(pe)))
        convertData(pe, barcodeList, 2, flowCellID, runIDConversion, verbose) ####
        #### convertData(pe, barcodeList, flowCellID, runIDConversion, verbose)
        demultiplex(pe, qualityConvert, verbose)
        writeData(pe, args.outputDirectory, tile, lane, 2, args.barcodeTable, args.useSampleNames, (gzipFile or args.gzip) and not args.noGzip, verbose)
    print("Runtime: %s" %(datetime.datetime.now() - starttime))
    
if __name__ == '__main__':
    args = CheckArgs()
    if args.barcodes2:
        dualBarcodes(args)
    else:
        main(args)
    
'''  
['CAGATCA',
'ATGTCAG',
'CTTGTAA',
'GGCTACA',
'TGACCAA',
'AGTCAAC',
'CGATGTA',
'CCGTCCC',
'AGTTCCG',
'GCCAATA',
'GTGGCCT',
'ACAGTGA',
'GTCCGCA',
'GTGAAAC']

CAGATCA,ATGTCAG,CTTGTAA,GGCTACA,TGACCAA,AGTCAAC,CGATGTA,CCGTCCC,AGTTCCG,GCCAATA,GTGGCCT,ACAGTGA,GTCCGCA,GTGAAAC
'''

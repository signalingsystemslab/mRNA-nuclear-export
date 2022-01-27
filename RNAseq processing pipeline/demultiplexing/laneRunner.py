#!/usr/bin/env python3

class CheckArgs():  #class that checks arguments and ultimately returns a validated set of arguments to the main program
    
    def __init__(self):
        import argparse
        import os
        parser = argparse.ArgumentParser()
        parser.add_argument("-i", "--inputDirectory", help = "Input Directory", required = True)
        parser.add_argument("-c", "--noConcatenation", help = "Don't concatenate tiles into one fastq", action = 'store_true')
        parser.add_argument("-b", "--barcodes", help = "Barcode file", required = True)
        parser.add_argument("-x", "--nextera", help = "Using Nextera dual barcodes.", action = 'store_true')
        parser.add_argument("-o", "--outputDirectory", help = "Directory for demultiplexed partial outputs", required = True)
        parser.add_argument("-n", "--noCleanup", help = "Do not remove partial files after concatenation", action = 'store_true')
        parser.add_argument("-m", "--forceMemory", help = "Force a specific RAM allocation for each demultiplex node", type = int)
        parser.add_argument("--mock", help = "Mock run, do not submit to scheduler.", action = "store_true")
        parser.add_argument("-p", "--multiprocessing", help = "Run with multiprocessing instead of through SGE scheduler.  User must pass the number of parallel processes to allow.", type = int)
        parser.add_argument("-s", "--qualityScore33Conversion", help = "Convert base 64 quality scores to base 33", action = "store_true")
        rawArgs = parser.parse_args()
        inputDirectory = rawArgs.inputDirectory
        if not os.path.isdir(inputDirectory):
            raise FileNotFoundError("Unable to find input directory at %s" %inputDirectory)
        self.inputDirectory = inputDirectory
        outputDirectory = rawArgs.outputDirectory
        if not os.path.isdir(outputDirectory):
            try:
                os.mkdir(outputDirectory)
                os.mkdir(outputDirectory + os.sep + "parts")
            except FileExistsError:
                pass
        if not os.path.isdir(outputDirectory):
            raise RuntimeError("Unable to find output directory %s and unable to create it." %outputDirectory)
        self.outputDirectory = rawArgs.outputDirectory
        self.noConcatenation = rawArgs.noConcatenation
        self.noCleanup = rawArgs.noCleanup or rawArgs.noConcatenation #if no concatenation is set, this will set noCleanup to true regardless to avoid automatically deleting the data
        barcodes = rawArgs.barcodes
        if not os.path.isfile(barcodes):
            raise FileNotFoundError("Unable to find barcode file at %s" %barcodes)
        self.barcodes = barcodes
        self.forceRam = rawArgs.forceMemory
        if self.forceRam and self.forceRam < 1:
            raise RuntimeError("Cannot set a memory value less than 1.  Passed value was %s" %self.forceRam)
        self.mock = rawArgs.mock
        self.nextera = rawArgs.nextera
        self.multiprocessing = rawArgs.multiprocessing
        if self.multiprocessing and self.multiprocessing < 1:
            raise RuntimeError("Cannot set a parallel process value less than 1.  Passed value was %s" %self.multiprocessing)
        self.qualityScore33Conversion = rawArgs.qualityScore33Conversion

class QseqFile(object):
    
    def __init__(self, file, directory):
        import os
        self.fileName = file
        nameSplit = file.split(".")
        if nameSplit[-1] == "gz":
            self.gzipped = True
        else:
            self.gzipped = False
        first, lane, read, tile, format = nameSplit[0].split("_")
        self.lane = int(lane)
        self.read = int(read)
        self.tile = int(tile)
        self.fullPath = directory + os.sep + file
        self.uncompressedSize = self.calculateSize()
        
    def __str__(self):
        return self.fullPath
        
    def calculateSize(self):
        import os
        if not self.gzipped:
            return os.path.getsize(self.fullPath)
        else:
            zipFile = open(self.fullPath, 'rb')
            zipFile.seek(-4, os.SEEK_END)
            sizeBytes = zipFile.read()
            zipFile.close()
            return int.from_bytes(sizeBytes, 'little')
        
def createFileTree(directory):
    import os
    lane = None
    largestFileSize = 0
    dirList = os.listdir(directory)
    fileTree = {}  #by tile and then by read
    for fileName in dirList:
        if "QSEQ.TXT" in fileName.upper() and os.path.isfile(os.path.join(directory, fileName)):
        #    print("Using %s" %fileName)
            file = QseqFile(fileName, directory)
            if type(lane) == type(None):
                lane = file.lane
            if not file.tile in fileTree:
                fileTree[file.tile] = {}
            if file.read in fileTree[file.tile]:
                raise RuntimeError("Name/read file collision between \n%s and \n%s" %(fileTree[file.tile][file.read], file))
            fileTree[file.tile][file.read] = file
            largestFileSize = max([largestFileSize, file.uncompressedSize])
        #else:
        #    print("Rejecting %s" %fileName)
    return (fileTree, largestFileSize, lane)

def getBarcodeDataFromFile(barcodeFileName):
    barcodeTable = []
    barcodeFile = open(barcodeFileName, 'r')
    line = barcodeFile.readline()
    while line:
        line = line.strip()
        if not line:
            line = barcodeFile.readline()
            continue
        line = line.split()
        if len(line) == 1:
            barcodeTable.append([line[0].strip(), line[0].strip()])
        elif len(line) == 2:
            barcodeTable.append([line[0].strip(), line[1].strip()])
        else:
            raise RuntimeError("Invalid line elements in barcode line %s" %line)
        line = barcodeFile.readline()
    barcodeLength = len(barcodeTable[0][0])
    for i in range(0, len(barcodeTable)):
        if not len(barcodeTable[i][0]) == barcodeLength:
            raise RuntimeError("Detected barcodes of different lengths in barcode table:\n%s" %barcodeTable)  
        barcodeTable[i][0] = barcodeTable[i][0].upper()
        for character in barcodeTable[i][0]:
            if not character in "ATGC":
                raise RuntimeError("Invalid DNA characters in barcode table:\n%s" %barcodeTable)
        if not barcodeTable[i][1].strip():
            barcodeTable[i][1] = barcodeTable[i][0]
    return (barcodeTable, barcodeLength)

def createBarcodeListArg(barcodeTable, barcodeLength):
    joinedPairs = ["".join(pair) for pair in barcodeTable]
    joinedArg = ",".join(joinedPairs)
    return str(barcodeLength) + "," + joinedArg

def dualBarcode(args):
    import arrayWrapper
    import os
    args = CheckArgs()
    forceRam = args.forceRam
    inputDir = args.inputDirectory
    outputDir = args.outputDirectory
    emailAddress = ""  #for debugging purposes
    emailCondition = "" #for debugging purposes
    barcodeFile = args.barcodes
    if not args.noConcatenation:
        partDir = outputDir + os.sep + "parts"
    else:
        partDir = outputDir
    fileTree, largestFile, lane = createFileTree(inputDir)
    demuxTempDir = arrayWrapper.createTempDir("demux%s" %lane, outputDir)
    if forceRam:
        ram = forceRam
    else:
        ram = max([int(largestFile * 10 //1000000000), 4])
    barcodeTable, barcodeLength = getBarcodeDataFromFile(barcodeFile)
    barcodeHash = {}
    for pair in barcodeTable:
        barcodeHash[pair[0]] = pair[1]
    expectedOutputs = {}  #barcode/read/tile
    for pair in barcodeTable:
        expectedOutputs[pair[0]] = {}
    demuxJobList = arrayWrapper.JobList(demuxTempDir)
    for tile in list(fileTree.keys()):
        job = ["filterAndConvert.py"]
        jobArgs = ["--pe1 %s" %fileTree[tile][1], "--barcodes %s" %fileTree[tile][2], "--barcodes2 %s" %fileTree[tile][3], "--outputDirectory %s" %partDir, "--validBarcodeList %s" %createBarcodeListArg(barcodeTable, barcodeLength)]
        for barcode in list(barcodeHash.keys()):
            if not 1 in expectedOutputs[barcode]:
                expectedOutputs[barcode][1] = []
            expectedOutputs[barcode][1].append(tile)
        if 4 in fileTree[tile]:
            jobArgs.append("--pe2 %s" %fileTree[tile][4])
            for barcode in list(barcodeHash.keys()):
                if not 2 in expectedOutputs[barcode]:
                    expectedOutputs[barcode][2] = []
                expectedOutputs[barcode][2].append(tile)
        if args.noConcatenation:
            jobArgs.append("--useSampleNames")
        if args.qualityScore33Conversion:
            jobArgs.append("--qualityScore33Conversion")
        fullJob = job + jobArgs
        demuxJobList.addJob(fullJob)
    demuxJobList.finish()
    filterAndConvert = arrayWrapper.HoffmanArrayJob([], demuxJobList.length, "demux%s" %lane, demuxTempDir, emailAddress, emailCondition, cores = 1, memory = ram, maxRetries = 10, mock = args.mock)
    if args.noConcatenation:
        quit("Not starting a concatenation process.")
    concatTempDir = arrayWrapper.createTempDir("concat%s" %lane, outputDir)
    concatJobList = arrayWrapper.JobMatrix(concatTempDir, 3)
    for barcode in list(barcodeHash.keys()):
        for read in list(expectedOutputs[barcode].keys()):
            if read == 2:
                originalRead = 4
            elif read == 1:
                originalRead = read
            else:
                raise RuntimeError("Got an inappropriate paired end read number.  Not sure how this happened.")
            concatenateTiles = sorted(expectedOutputs[barcode][read])
            partialFilesList = []
            for tile in concatenateTiles:
                if fileTree[tile][originalRead].gzipped:
                    gzipOut = True
                    partialFileName = ".".join([barcode, str(lane), str(tile), str(read), "part", "gz"])
                    partialFileName = partDir + os.sep + partialFileName
                else:
                    gzipOut = False
                    partialFileName = ".".join([barcode, str(lane), str(tile), str(read), "part"])
                    partialFileName = partDir + os.sep + partialFileName
                partialFilesList.append(partialFileName)
            outputFileName = outputDir + os.sep + barcodeHash[barcode] + "." + str(lane) + "." + str(read) + ".fastq"
            if gzipOut:
                outputFileName += ".gz"
            if args.noCleanup:
                concatJobList.addJob(["cat", " ".join(partialFilesList), ">", outputFileName])
            else:
                concatJobList.addJobList([["cat", " ".join(partialFilesList), ">", outputFileName], ["rm -rf", " ".join(partialFilesList)]])
    concatJobList.finish()
    concatenate = arrayWrapper.HoffmanArrayJob(filterAndConvert.jobID, concatJobList.length, "concat%s" %lane, concatTempDir, emailAddress, emailCondition, cores = 1, memory = 2, maxRetries = 10, mock = args.mock)
            
def main(args):
    import arrayWrapper
    import os
    args = CheckArgs()
    forceRam = args.forceRam
    inputDir = args.inputDirectory
    outputDir = args.outputDirectory
    emailAddress = ""  #for debugging purposes
    emailCondition = "" #for debugging purposes
    barcodeFile = args.barcodes
    if not args.noConcatenation:
        partDir = outputDir + os.sep + "parts"
    else:
        partDir = outputDir
    fileTree, largestFile, lane = createFileTree(inputDir)
    demuxTempDir = arrayWrapper.createTempDir("demux%s" %lane, outputDir)
    if forceRam:
        ram = forceRam
    else:
        ram = max([int(largestFile * 10 //1000000000), 4])
    barcodeTable, barcodeLength = getBarcodeDataFromFile(barcodeFile)
    barcodeHash = {}
    for pair in barcodeTable:
        barcodeHash[pair[0]] = pair[1]
    expectedOutputs = {}  #barcode/read/tile
    for pair in barcodeTable:
        expectedOutputs[pair[0]] = {}
    demuxJobList = arrayWrapper.JobList(demuxTempDir)
    for tile in list(fileTree.keys()):
        job = ["filterAndConvert.py"]
        jobArgs = ["--pe1 %s" %fileTree[tile][1], "--barcodes %s" %fileTree[tile][2], "--outputDirectory %s" %partDir, "--validBarcodeList %s" %createBarcodeListArg(barcodeTable, barcodeLength)]
        for barcode in list(barcodeHash.keys()):
            if not 1 in expectedOutputs[barcode]:
                expectedOutputs[barcode][1] = []
            expectedOutputs[barcode][1].append(tile)
        if 3 in fileTree[tile]:
            jobArgs.append("--pe2 %s" %fileTree[tile][3])
            for barcode in list(barcodeHash.keys()):
                if not 2 in expectedOutputs[barcode]:
                    expectedOutputs[barcode][2] = []
                expectedOutputs[barcode][2].append(tile)
        if args.noConcatenation:
            jobArgs.append("--useSampleNames")
        if args.qualityScore33Conversion:
            jobArgs.append("--qualityScore33Conversion")
        fullJob = job + jobArgs
        demuxJobList.addJob(fullJob)
    demuxJobList.finish()
    if not args.multiprocessing:
        filterAndConvert = arrayWrapper.HoffmanArrayJob([], demuxJobList.length, "demux%s" %lane, demuxTempDir, emailAddress, emailCondition, cores = 1, memory = ram, maxRetries = 10, mock = args.mock)
    else:
        arrayWrapper.multiprocessJob(demuxJobList.jobList, args.multiprocessing)
    if args.noConcatenation:
        quit("Not starting a concatenation process.")
    concatTempDir = arrayWrapper.createTempDir("concat%s" %lane, outputDir)
    concatJobList = arrayWrapper.JobMatrix(concatTempDir, 3)
    for barcode in list(barcodeHash.keys()):
        for read in list(expectedOutputs[barcode].keys()):
            if read == 2:
                originalRead = 3
            elif read == 1:
                originalRead = read
            else:
                raise RuntimeError("Got an inappropriate paired end read number.  Not sure how this happened.")
            concatenateTiles = sorted(expectedOutputs[barcode][read])
            partialFilesList = []
            for tile in concatenateTiles:
                if fileTree[tile][originalRead].gzipped:
                    gzipOut = True
                    partialFileName = ".".join([barcode, str(lane), str(tile), str(read), "part", "gz"])
                    partialFileName = partDir + os.sep + partialFileName
                else:
                    gzipOut = False
                    partialFileName = ".".join([barcode, str(lane), str(tile), str(read), "part"])
                    partialFileName = partDir + os.sep + partialFileName
                partialFilesList.append(partialFileName)
            outputFileName = outputDir + os.sep + barcodeHash[barcode] + "." + str(lane) + "." + str(read) + ".fastq"
            if gzipOut:
                outputFileName += ".gz"
            if args.noCleanup:
                concatJobList.addJob(["cat", " ".join(partialFilesList), ">", outputFileName])
            else:
                concatJobList.addJobList([["cat", " ".join(partialFilesList), ">", outputFileName], ["rm -rf", " ".join(partialFilesList)]])
    concatJobList.finish()
    if not args.multiprocessing:
        concatenate = arrayWrapper.HoffmanArrayJob(filterAndConvert.jobID, concatJobList.length, "concat%s" %lane, concatTempDir, emailAddress, emailCondition, cores = 1, memory = 2, maxRetries = 10, mock = args.mock)
    else:
        arrayWrapper.multiprocessJob(concatJobList.jobList, args.multiprocessing)

if __name__ == '__main__':
    args = CheckArgs()
    if args.nextera:
        dualBarcode(args)
    else:
        main(args)
        

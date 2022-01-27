#!/usr/bin/env python3

global genericRunnerPaths
genericRunnerPaths = {#"python3" : "/u/local/apps/python/3.4.3/bin/python3",  #hoffman
                      #"python3" : "/Library/Frameworks/Python.framework/Versions/3.4/bin/python3"  #local for debugging
                      "qsub" : "/u/systems/UGE8.0.1vm/bin/lx-amd64/qsub",
                      "arrayWrapper" : "./arrayWrapper.py"}

def multiprocessJob(jobList, maxSubprocesses):
    import functools
    import multiprocessing.dummy
    import subprocess
    commands = []
    for job in jobList:
        if type(job) == str:
            if job == "Zeroth job placeholder":
                continue
            commands.append(job)
        elif type(job) == list:
            for subjob in job:
                commands.append(subjob)
    pool = multiprocessing.dummy.Pool(maxSubprocesses) # two concurrent commands at a time
    for i, returncode in enumerate(pool.imap(functools.partial(subprocess.call, shell=True), commands)):
        print("Completed %s of %s" %(i+1, len(commands)))
        if returncode != 0:
            print("Command %s failed: Exit status %s on %s" % (i, returncode, " ".join(commands[i])))
        

class HoffmanArrayJob(object):
    
    def __init__(self, dependencies, jobRange, jobName, tempDir, emailAddress, emailConditions, cores = 1, memory = 6, maxRetries = 10, mock = False):
        import os
        self.dependencies = dependencies
        if type(self.dependencies) in (int, str):
            self.dependencies = [self.dependencies]
        self.jobName = jobName
        if self.jobName[0].isdigit():
            self.jobName = "j" + self.jobName
        if type(jobRange) in (list, tuple):
            assert len(jobRange) <= 2 and len(jobRange) > 0, "Job range must either be a single integer or a tuple/list containing two integers.  Passed: %s" %jobRange
            for item in jobRange:
                assert type(item) == int, "Job range values must be integers.  Value %s was a %s." %(item, type(item))
            if len(jobRange) == 1:
                self.endJob = jobRange[0]
                self.startJob = 1
            else:
                assert jobRange[0] <= jobRange[1], "The end of the job range must be greater than or equal to the start of the job range.  Start: %s   End: %s" %(jobRange[0], jobRange[1])
                self.startJob, self.endJob = jobRange
        elif type(jobRange) == int:
            self.endJob = jobRange
            self.startJob = 1
        assert self.startJob >= 1, "Start of job range must be greater than or equal to 1.  Cluster job arrays are indexed to 1 and job zero will return as invalid.  Start value was %s" %self.startJob
        self.emailAddress = emailAddress
        self.emailConditions = emailConditions
        self.mock = mock
        self.cores = cores
        if not type(cores) == int:
            raise RuntimeError("Core count must be an integer")
        self.memory = memory
        if not type(memory) == int or type(memory) == float:
            raise RuntimeError("Memory limit must be a number")
        self.memoryPerCore = self.calculateMemoryPerCore()
        self.tempDir = tempDir
        if not os.path.isdir(self.tempDir):
            raise RuntimeError("Unable to find working directory")
        self.jobID = self.submitJob(maxRetries)
        
    def calculateMemoryPerCore(self):  #because pe shared mode tends to give you your desired memory per core
        if self.cores == 1:
            return self.memory
        if self.memory % self.cores == 0:  #checking to see if it is evenly divisible
            return self.memory // self.cores
        return round(float(self.memory)/self.cores, 1)
        
    # def clearedDependencies(self, dependencies):
    #     import os
    #     for dependency in dependencies:
    #         if not os.path.isfile(dependency):
    #             return False
    #         if os.path.getsize(dependency) == 0:  #make sure it's not an empty file
    #             return False
    #     return True
    # 
    # def pulledFromJobBoard(self, jobNumber):  #job board will essentially have files as job tokens.  Only one process should be able to delete each file.
    #     import os
    #     jobBoardFile = self.tempDir + os.sep + "jobBoard" + os.sep + str(jobNumber) + ".job"
    #     try:
    #         os.remove(jobBoardFile)
    #     except FileNotFoundError:
    #         return False
    #     return True
    # 
    # def jobReady(self):
    #     if self.clearedDependencies(self.dependencies):
    #         if self.pulledFromJobBoard(self.jobNumber):  #very important they get tested in this order.  This line will take the token.
    #            return True
            
    def submitJob(self, maxRetries = 10):
        import os
        import subprocess
        peArg = []
        if self.cores > 1:
            peArg = ["-pe", "shared", str(self.cores)]
        limitsArgs = ["-l", "h_rt=23:59:59,h_data=" + str(self.memory) + "G"]
        jobRangeArg = ["-t", "%s-%s" %(self.startJob, self.endJob)]
        outputsDir = self.tempDir + os.sep + "outputs"
        outputsArg = ["-o", outputsDir, "-e", outputsDir]
        jobNameArg = ["-N", self.jobName]
        emailArg = []
        if self.emailAddress:
            emailArg = ["-M", self.emailAddress, "-m", self.emailConditions]
        holdArg = []
        if self.dependencies:
            dependencies = [str(dependency) for dependency in self.dependencies]
            holdArg = ["-hold_jid", ",".join(dependencies)]
        startingArgs = [genericRunnerPaths["qsub"], "-cwd"]
        schedulerCommandList = startingArgs + jobNameArg + limitsArgs + peArg + outputsArg + jobRangeArg + emailArg + holdArg
        commandToExecute = [genericRunnerPaths["python3"], genericRunnerPaths["arrayWrapper"], "-d", self.tempDir]
        commandToExecuteString = " ".join(commandToExecute)
        if self.mock:  #trap mock submissions here
            print("MOCK SUBMIT: " + " ".join(schedulerCommandList))
            print("MOCK COMMUNICATE: " + commandToExecuteString)
            print("Mock job, no number assigned")
            return 123456   #placeholder number
        successfullySubmitted = False
        retries = 0
        while not successfullySubmitted:
            print()
            print("JOB: %s" %(commandToExecuteString))
            print("SUB: %s" %(" ".join(schedulerCommandList)))
            child = subprocess.Popen(schedulerCommandList, stdout = subprocess.PIPE, stdin = subprocess.PIPE, stderr = subprocess.PIPE)
            childOut, childErr = child.communicate(input = commandToExecuteString.encode())
            childExitStatus = child.returncode
            qsubOut = childOut.decode().strip()
            qsubError = childErr.decode().strip()
            successfullySubmitted = childExitStatus == 0  #if it went through, it should exit status zero
            if successfullySubmitted:
                import re
                print("QSUB: %s" %(qsubOut))
                regex = re.search('Your job.* (\d+)', qsubOut)
                return int(regex.group(1))
            else:
                if retries >= maxRetries:
                    raise RuntimeError("Failed to qsub after %s attempts.\nQSUB OUT: %s\n QSUB ERR: %s" %(retries, qsubOut, qsubError))
                else:
                    print("Submission to qsub failed. Retrying.\nQSUB OUT: %s\n QSUB ERR: %s" %(qsubOut, qsubError))
                    retries += 1

class JobList(object):
    
    def __init__(self, tempDir):
        import os
        self.jobFile = open(os.path.join(tempDir, "jobs.pkl"), 'wb')
        self.jobList = ["Zeroth job placeholder"]
        self.length = 0
        
    def addJob(self, args, listJob = False):
        if not self.jobFile:
            raise RuntimeError("Can't add jobs after job file has been finished.")
        argList = args
        if argList[0].upper() == "PYTHON3":
            argList[0] = genericRunnerPaths["python3"]
        argList = [str(arg) for arg in argList]
        if listJob:
            return (" ".join(argList))
        self.jobList.append(" ".join(argList))
        self.length += 1
        
    def addJobList(self, jobList, dependent = True):
        processedJobList = []
        for job in jobList:
            processedJobList.append(self.addJob(job, True))
        if dependent:
            processedJobList = " && ".join(processedJobList)
        self.jobList.append(processedJobList)
        self.length += 1
        
    def finish(self):
        import pickle
        pickle.dump(self.jobList, self.jobFile)
        self.jobFile.close()
        self.jobFile = False

class JobMatrix(object):  #Just different enough that I don't want to inherit from JobList
    
    def __init__(self, tempDir, width):
        import os
        try:
            width = int(width)
        except:  #using a catch-all exception here because we will sort this out properly in a moment.
            pass
        assert type(width) == int, "Width value for job matrix must be an integer. A %s was passed instead." %type(width)
        assert width > 0, "Width must be a positive value.  Value given was %s" %width
        self.width = width
        self.jobFile = open(os.path.join(tempDir, "jobs.pkl"), 'wb')
        self.jobList = ["Zeroth job placeholder", []]
        self.length = 1
        
    def addJob(self, args, multiJob = False):
        if not self.jobFile:
            raise RuntimeError("Can't add jobs after job file has been finished.")
        argList = args
        if argList[0].upper() == "PYTHON3":
            argList[0] = genericRunnerPaths["python3"]
        argList = [str(arg) for arg in argList]
        if multiJob:
            return (" ".join(argList))
        if len(self.jobList[self.length]) >= self.width:
            self.jobList.append([])
            self.length += 1
        self.jobList[self.length].append(" ".join(argList))
        
    def addJobList(self, jobList, dependent = True):
        processedJobList = []
        for job in jobList:
            processedJobList.append(self.addJob(job, True))
        if dependent:
            joiner = " && "
        else:
            joiner = " ; "
        processedJobString = joiner.join(processedJobList)
        if len(self.jobList[self.length]) >= self.width:
            self.jobList.append([])
            self.length += 1
        self.jobList[self.length].append(processedJobString)
        
    def finish(self):
        import pickle
        pickle.dump(self.jobList, self.jobFile)
        self.jobFile.close()
        self.jobFile = False
                    
def createTempDir(name, useDir = False):
    import re
    import os
    import datetime
    successful = False
    if not useDir:
        useDir = os.getcwd()
    while not successful:
        currenttime = datetime.datetime.now()
        currenttime = str(currenttime)
        currenttime = re.sub(r'\W','',currenttime)
        tempDir = useDir + os.sep + "." + name + currenttime
        if os.path.isdir(tempDir):
            continue
        try:
            os.mkdir(tempDir)
        except OSError:
            continue
        successful = True
    os.mkdir(tempDir + "/outputs")
    return tempDir

class CheckArgs():  #class that checks arguments and ultimately returns a validated set of arguments to the main program
    
    def __init__(self):
        import argparse
        import os
        parser = argparse.ArgumentParser()
        parser.add_argument("-d", "--directory", help = "Directory with job data")
        parser.add_argument("-j", "--job", help = "Force this instance to use a job number without looking at env variables.", type = int)
        parser.add_argument("-s", "--subjob", help = "Force this instance to only run a specific subjob from the line")
        parser.add_argument("--noClockOut", help = "Do not clock out via this wrapper.", action = 'store_true')
        parser.add_argument("-i", "--internetRequired", help = "Reschedule this job if the node is unable to reach the internet.", action = "store_true")
        rawArgs = parser.parse_args()
        if rawArgs.directory:
            if os.path.isdir(rawArgs.directory):
                self.directory = rawArgs.directory
            else:
                raise RuntimeError("Temporary directory not found: " + rawArgs.directory)
        else:
            raise RunTimeError("No temporary directory specified.")
        if rawArgs.job:
            self.job = rawArgs.job
        else:
            self.job = False
        subjob = rawArgs.subjob
        if subjob and not self.job:
            raise RuntimeError("A job is needed if a subjob is specified.")
        if subjob:
            try:
                subjob = int(subjob)
            except ValueError:
                raise RuntimeError("Subjob must be specified as an integer.")
            self.subjob = subjob
        else:
            self.subjob = -1
        self.noClockOut = rawArgs.noClockOut
        self.internetRequired = rawArgs.internetRequired
        
def main():
    import datetime
    startTime = datetime.datetime.now()
    import os  #import the library for making os system calls
    import pickle
    import random
    import time
    import os
    os.environ["PATH"] += ":/u/local/apps/R/current/bin"  #adding R to the path for the subshell because some GATK functions need this
    global args  #declare args as a global
    args = CheckArgs()  #get an object containing validated arguments
    try:
        nodeNumber = os.environ["HOSTNAME"].replace("n","")
        nodeNumber = int(nodeNumber)
    except ValueError:
        nodeNumber = False
    if args.internetRequired:
        testInternetConnection = os.system("ping 8.8.8.8 -c 1 -W 2")  #fire off a ping to google's open DNS
        if not testInternetConnection == 0:  #ping operation will return a non-zero value if it fails
            import sys
            sys.exis(99)
    if not args.job:
        try:
            jobListNumber = int(os.environ["SGE_TASK_ID"])   #get the array job number from environmental variables
        except KeyError:  #if it cannot get that value
            raise RuntimeError("Unable to find a valid task ID in OS environment variables.")   #something is wrong so quit
    else:
        jobListNumber = args.job
    directory = args.directory  #get the tempdir from arguments
    if not directory.endswith(os.sep):  #if it does not end with a separator
        directory += os.sep  #add one
    jobListFile = open(directory + "jobs.pkl",'rb')
    jobLists = pickle.load(jobListFile)
    jobListFile.close()
    thisJobList = jobLists[jobListNumber]
    if type(thisJobList) == str:
        thisJobList = [thisJobList]
    # random.seed(jobListNumber)  #introducing some jitter here
    # delay = random.uniform(0,60)
    # time.sleep(delay)
    if args.subjob != -1:
        jobRange = [args.subjob]
    else:
        jobRange = range(0, len(thisJobList))
    if not args.noClockOut:
        startTouchFile = directory + str(jobListNumber) + ".runner.started"
        touchFile = open(startTouchFile, 'w')
        touchFile.close()
    for thisJobNumber in jobRange:
        if os.path.isfile(directory + str(jobListNumber) + "." + str(thisJobNumber) + ".done") or os.path.isfile(directory + str(jobListNumber) + "." + str(thisJobNumber) + ".already.done"):
            print("Job %s.%s has already been marked as done.  Skipping it." %(jobListNumber, thisJobNumber))
            continue
        print("Running job %s.%s" %(jobListNumber, thisJobNumber))
        jobStatus = os.system(thisJobList[thisJobNumber])  #run the bash file we just identified and set jobStatus to its exit status
        print("Job exit status: %s" %(jobStatus))
        if not args.noClockOut:
            if jobStatus == 0:  #if the job finished successfully
                touchFilePath = directory + str(jobListNumber) + "." + str(thisJobNumber) + ".done"  #define our clockout file
                touchFile = open(touchFilePath, 'w')  #create our clockout file
                touchFile.close()  #close it without writing anything
                failedFile = directory + str(jobListNumber) + "." + str(thisJobNumber) + ".failed"
                if os.path.isfile(failedFile):
                    os.remove(failedFile)
            elif jobStatus == 10752 or jobStatus == 42:  #I have no idea why I tell it to exit 42 and it exits 10752, but I can't be arsed to fix it right now.
                touchFilePath = directory + str(jobListNumber) + "." + str(thisJobNumber) + ".already.done"  #define our clockout file
                touchFile = open(touchFilePath, 'w')  #create our clockout file
                touchFile.close()  #close it without writing anything
                failedFile = directory + str(jobListNumber) + "." + str(thisJobNumber) + ".failed"
                if os.path.isfile(failedFile):
                    os.remove(failedFile)
            else:
                try:
                    nodeNumber = os.environ["HOSTNAME"]
                except ValueError:
                    nodeNumber = "Unable to get node number"
                touchFilePath = directory + str(jobListNumber) + "." + str(thisJobNumber) + ".failed"  #define our clockout file
                touchFile = open(touchFilePath, 'w')  #create our clockout file
                touchFile.write("Failed on node " + nodeNumber + "\n")
                touchFile.close()  #close it without writing anything
                import sys
                sys.exit(100)  #this exit status should tell the cluster that the job failed and not to proceed with the next steps.  They will have to be manually purged from teh queue.
        endTouchFile = directory + str(jobListNumber) + ".runner.ended"
        touchFile = open(endTouchFile, 'w')
        touchFile.close()
    while not (datetime.datetime.now() - startTime).seconds > 200:
        time.sleep(1)

if __name__ == '__main__':
    main()

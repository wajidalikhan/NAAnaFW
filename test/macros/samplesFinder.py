import os, commands


###This funcTion searches recursively for root files in a remote path, requesting the path contains the word "sampleToFind" and returning an output file list:
def findRootFileInPath(ls_command,srm, initialPath,sampleToFind,extraString="",date="", veto = "failed"):
       
    ls_command_exec = ls_command+srm+initialPath
    print "command is: "+ ls_command_exec
    status,ls_la = commands.getstatusoutput( ls_command_exec )
    list = []
    list = ls_la.split(os.linesep)
    
    outlist =[]
    fetchlist =[]
    pulledlist =[]
    #print list
        
    for l in list:
        if l == initialPath:
#            print " can't search further! "
            continue
#       print "list element is: "+ l
#        print "sampleToFind is: "+ sampleToFind
#       print l
#       print date
#       print extraString 
        if (sampleToFind in l and not "failed" in l and not veto in l and (((extraString == "") or (extraString in l))) and (((date == "") or (("/"+date) in l))) ): 
#            print date
#            print extraString 
#            print " there's  sample in this element! "
#            print l
#            print "veto is " + veto + " is it in l? " + str(veto in l)
        
            if ".root" in l:
                #       print " there's a rootfile in this element! "
                outlist.append(l)
            else:
                fetchlist.append(l)
        else:
            if ".root" in l:
                continue
            else: fetchlist.append(l)

    if(len(outlist)==0):
        for l in fetchlist:
            pulledlist = findRootFileInPath(ls_command=ls_command,srm=srm, initialPath=l,sampleToFind=sampleToFind,extraString = extraString, date = date, veto = veto)

            print " pulled list:"
            print pulledlist
            for f in pulledlist:
                outlist.append(f)

    print "path is "+ initialPath
    print " check outlist length "
    print len(outlist)
#    if status: raise RFIOErr('Dir not found')
    return outlist

###This one copies a files in an input list to another directory, returning a list wiht the copied file names:
def copyFilesInList(cp_command,srm, inputList, outputDir, maxSimultaneousCopies=1):
    
    m =0
    command_cp_exec = cp_command+srm
    for p in inputList:
        outName = outputDir
        words = p.split("/")
        for word in words:
            if(".root" in word):
                outName = outputDir+"/"+word
        command_ls_local = "ls "+ outName        

        check_doubles = commands.getstatusoutput(command_ls_local)
        
        print check_doubles
        if check_doubles[0]==512: 
            final=" &"
            m=m+1;
            if (m % maxSimultaneousCopies ==0):
                final = ""
#            print "final is: " + final
            command_cp_exec = cp_command+srm+p + " "+ outName +final
            print "command cp is :" +command_cp_exec
            os.system(command_cp_exec)

###This one copies a files in an input list to another directory, returning a list wiht the copied file names:
def printFilesInList(inputList, path, nameOut, outputDir, option = 'w'):
 
    filename = outputDir+"list"+nameOut+".py"
    if (option == 'a'):
        filename = outputDir+"listAllSamples.py"
        
    f = open(filename, option)
    

#    f.write("import FWCore.ParameterSet.Config as cms \n")
    if(option == 'a'): f.write("filelist"+nameOut.replace("-","_")+" = [")
    if(option == 'w'): f.write("filelist = [")
    for l in inputList:
        for lsplit in l.split(path):
            print lsplit
            if (".root" in lsplit): f.write("'root://xrootd.ba.infn.it//"+lsplit+"', \n")
    f.write("]\n")        
        
###MAIN###

cmdls = "lcg-ls"
cmdcp = "lcg-cp"
srm = "  -b -D srmv2 srm://storage01.lcg.cscs.ch:8443/srm/managerv2\?SFN="

#srm = " -b -D srmv2 srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN="


#samplesAndPaths = {"TTJets":"/pnfs/lcg.cscs.ch/cms/trivcat/store/user/decosa/ttDM/Phys14_Tree_v2/"}
samplesAndPaths = {
    "ST_t-channel":"/pnfs/lcg.cscs.ch/cms/trivcat/store/user/oiorio/ttDM/samples/2016/Oct/",
    }

def writeSamplesAndPaths(sap, outputDir="./", option = 'w'):
    filename = outputDir+"samplesAndPaths.py"
        
    f = open(filename, option)
    
    f.write("samplesAndPaths = {")

    for sample in sap:
        f.write("'"+str(sample)+"':'"+str(sap[sample])+"',")
    f.write("}")  
  
outDir = "/tmp/oiorio/"

writeSamplesAndPaths(samplesAndPaths)
os.system("mkdir samples")
os.system("rm samples/listAllSamples.py")
#os.system("rm jer/listAllSamples.py")

for sample in samplesAndPaths:
    path = samplesAndPaths[sample]

#    outputList = []
    outputListJES = []
    print path
    outputListJES= findRootFileInPath(cmdls,srm,path,sample,extraString ="",date="", veto = "failed")



    print "final output list for "+sample+" is : "
    print outputListJES

    printFilesInList(outputListJES, "/pnfs/lcg.cscs.ch/cms/trivcat/", sample,"./samples/","w")
    printFilesInList(outputListJES, "/pnfs/lcg.cscs.ch/cms/trivcat/", sample,"./samples/","a")

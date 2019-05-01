import os
import shutil
import glob
import subprocess
import fileinput
import pandas as pd 
import matplotlib.pyplot as plt
import math

class Storm:
  def __init__(self, timestep, cov_th, radius_th, area, speedx, speedy, i10, i25, i50, i75, i90, loc):
      self.timestep = timestep
      self.cov_th = cov_th
      self.radius_th=radius_th
      self.area=area
      self.speedx=speedx
      self.speedy=speedy
      self.i10=i10
      self.i25=i25
      self.i50=i50
      self.i75=i75
      self.i90=i90
      self.loc=loc

ref_dir= "/home/z5194283/hdrive/MET_Tutorial/MyData/RealData/FinalData/MRMS-Regridded/"
IMERG_dir= "/home/z5194283/hdrive/MET_Tutorial/MyData/RealData/FinalData/IMERG/"
MRMS_Remap_NN_dir = "/home/z5194283/hdrive/MET_Tutorial/MyData/RealData/FinalData/MRMS-RemapNN/"
MRMS_Fixed_dir = "/home/z5194283/hdrive/MET_Tutorial/MyData/RealData/FinalData/MRMS-Regridded-Fixed/"
dst_dir= "/home/z5194283/hdrive/MET_Tutorial/MyData/RealData/FinalData/Structured/"
SampleConfigAddress = "/home/z5194283/hdrive/MET_Tutorial/MyData/RealData/FinalData/SampleConfig"

# get the Full address of files in each data group
reference_files=glob.glob(ref_dir+"*.nc")
reference_files.sort()
IMERG_files = glob.glob(IMERG_dir+"*.nc")
IMERG_files.sort()
MRMS_Remap_NN_files = glob.glob(MRMS_Remap_NN_dir+"*.nc")
MRMS_Remap_NN_files.sort()
MRMS_Fixed_files = glob.glob(MRMS_Fixed_dir+"*.nc")
MRMS_Fixed_files.sort()

# get the names of files in reference data: here is MRMS-Regridded
dirs = os.listdir(ref_dir)

# Creating Folders based on the Observatory Data: Here the reference data is MRMS-Regridded
i=0
R=[0.5]
T=[2]
#CheckExistingFolders = glob.glob(dst_dir+"/*")
if len(os.listdir(dst_dir) ) == 0:
    for namefile in dirs:         
        try:
            os.makedirs(dst_dir+namefile[:-3])
            newdir = dst_dir+namefile[:-3]
            os.makedirs(newdir + "/IMERG")
            os.makedirs(newdir + "/MRMS")
            os.makedirs(newdir + "/MRMSNN")
            os.makedirs(newdir + "/MRMSFixed")
            
            #defining the number of thresholds and Radius

            
            #Creating output folders
            for r in R:
                for t in T:
                    os.makedirs(newdir + "/Output-IMERG-MRMS-R"+str(r)+"-T"+str(t))
            os.makedirs(newdir + "/Output-MRMSNN-MRMS")           
                      
            # Copying files into different folders
            shutil.copy(reference_files[i], newdir + "/MRMS"+"/")
            shutil.copy(reference_files[i+1], newdir + "/MRMS"+"/")
            shutil.copy(reference_files[i+2], newdir + "/MRMS"+"/")
            shutil.copy(IMERG_files[i], newdir + "/IMERG"+"/")
            shutil.copy(IMERG_files[i+1], newdir + "/IMERG"+"/")
            shutil.copy(IMERG_files[i+2], newdir + "/IMERG"+"/")
            shutil.copy(MRMS_Remap_NN_files[i], newdir + "/MRMSNN"+"/")
            shutil.copy(MRMS_Remap_NN_files[i+1], newdir + "/MRMSNN"+"/")
            shutil.copy(MRMS_Remap_NN_files[i+2], newdir + "/MRMSNN"+"/")
            shutil.copy(MRMS_Fixed_files[i], newdir + "/MRMSFixed"+"/")
            shutil.copy(MRMS_Fixed_files[i+1], newdir + "/MRMSFixed"+"/")
            shutil.copy(MRMS_Fixed_files[i+2], newdir + "/MRMSFixed"+"/")
            i=i+1

        except:
            
            print("last loops!")
else:
    print("All the input folders has been already created !")

# getting the complete address of all created folders
new_folders=glob.glob("/home/z5194283/hdrive/MET_Tutorial/MyData/RealData/FinalData/Structured/*")
new_folders.sort()

# creating config folders and files
ConfigFolderAddress = "/home/z5194283/hdrive/MET_Tutorial/MyData/RealData/FinalData/Structured-Configs"
for r in R:
    for t in T:
        if not os.path.exists(ConfigFolderAddress):
            os.makedirs(ConfigFolderAddress)  
        shutil.copy(SampleConfigAddress, ConfigFolderAddress + "/IM-R" + str(r) + "-T" +str(t))            

# modifying config file
for r in R:
    for t in T:
        filename = ConfigFolderAddress + "/IM-R" + str(r) + "-T" +str(t)
        
        file = open(filename,"r")
        text = file.readlines()
        file.close()
        oldR = text[64][3:26]
        oldT = text[65][3:28]
        
        newR="conv_radius       = " + str(r)
        newT="conv_thresh       = >=" + str(t)
        
        #replacing the thresholds in the config files
        with fileinput.FileInput(filename, inplace=True, backup='.bak') as file:
            for line in file:
                print(line.replace(oldR, newR), end='')            
        with fileinput.FileInput(filename, inplace=True, backup='.bak') as file:
            for line in file:        
                print(line.replace(oldT, newT), end='')    
        os.remove(filename+".bak")


 
# Running MTD for different parameters

CheckExistingFile2=0 #new_folders[0]+"/IMERG"

CheckExistingFolders1 = glob.glob(new_folders[0]+"/*")
for folderss in CheckExistingFolders1:
    if "Output" in folderss:
       CheckExistingFile2=len(glob.glob(folderss+"/*"))
#       CheckExistingFile2 = Vaset[0] 
       break
if(CheckExistingFile2 >= 1):
    print("The process has been already done !")

else:
    
    for folder in new_folders:
        for r in R:
            for t in T:
                LoadMTD =  "module load gsl/2.3" +" \n" + "module load netcdf-c/4.4.1.1" + " \n" + "module load netcdf-c++/4-4.3.0" +" \n" + "module load met/8.20" + " \n" + "mtd \\"
                ObsAddress =    folder +"/MRMS/*.nc \\"
                FcstAddress = folder +"/IMERG/*.nc \\"
                ConfigAddress = "/home/z5194283/hdrive/MET_Tutorial/MyData/RealData/FinalData/Structured-Configs/IM-R"+str(r)+"-T"+str(t)+" \\"
                OutputAddress = folder+"/Output-IMERG-MRMS-R"+str(r)+"-T"+str(t)+" \\"
                Script= LoadMTD + "\n" + "-fcst " + FcstAddress + "\n" + "-obs " + ObsAddress + "\n" + "-config " + ConfigAddress + "\n" + "-outdir " + OutputAddress + "\n" + "-v 2" + "\n"    
                subprocess.getstatusoutput(Script)

# modifying the txt tables (removing the additional spaces)
    oldtxt="0 m above mean sea level"
    newtxt="0_m_above_mean_sea_level"
    for folder in new_folders:
        for r in R:
            for t in T:
                txtOutputAddress = glob.glob(folder+"/Output-IMERG-MRMS-R"+str(r)+"-T"+str(t)+"/*.txt")
                for txtoutputfiles in txtOutputAddress:
                    with fileinput.FileInput(txtoutputfiles, inplace=True, backup='.bak') as file:
                        for line in file:
                            print(line.replace(oldtxt, newtxt), end='') 
                    os.remove(txtoutputfiles+".bak")



# Getting Outputs From the Created txt files
IMERG=list()
MRMS=list()
Time=0
try: 
    os.remove('/home/z5194283/hdrive/MET_Tutorial/MyData/RealData/FinalData/log_file.txt')
except:
    k=0
j=0

# if this index equals = 0 means the information of a correcct data is imported to the object if it is 1 therefore the file is corrupted and the information should be zero to add to the file
FileCorruptFinder = 0
for folder in new_folders:
    for r in R:
        for t in T:
            address = glob.glob(folder+"/Output-IMERG-MRMS-R"+str(r)+"-T"+str(t)+"/*.txt")
#            address = "/home/ho0man/Mygit/Structured/MRMS_PrecipRate_00.00_20180913-213000-compressed-l5-double/Output-IMERG-MRMS-R2-T2/mtd_20180913_212959V_3d_pair_simple.txt"
            for filez in address:
                if "3d_pair_simple" in filez:    
                    
                    table = pd.read_table(filez, delim_whitespace=True)
                    
                    #finding the cluster with highest interest
                    highestInterest = max(table["INTEREST"].values)
                    index1 = int(table[table["INTEREST"]==highestInterest].index[0])
                    PairedClusters=table["OBJECT_CAT"][index1]
                    ClusterNames = PairedClusters.split("_")
                    
            
            for filezz in address:
                try:
                    if "3d_single_cluster" in filezz:
                        
                        # finding the Characteristcs of Clusters in another txt Table                
                        table2 = pd.read_table(filezz, delim_whitespace=True)
                        
                        #Properties of Forecast and Observed Cluster
                        FcstIndex=int(table2[table2["OBJECT_CAT"]==ClusterNames[0]].index[0])
                        Fcsti10=table2["INTENSITY_10"].values[FcstIndex]
                        Fcsti25=table2["INTENSITY_25"].values[FcstIndex]
                        Fcsti50=table2["INTENSITY_50"].values[FcstIndex]
                        Fcsti75=table2["INTENSITY_75"].values[FcstIndex]
                        Fcsti90=table2["INTENSITY_90"].values[FcstIndex]
                        FcstLoc=list()
                        FcstLoc.append([table2["CENTROID_X"].values[FcstIndex],table2["CENTROID_Y"].values[FcstIndex]])
                        FcstSpeedX=table2["X_DOT"].values[FcstIndex]
                        FcstSpeedY=table2["Y_DOT"].values[FcstIndex]
                        NumberofFilesInaDataFolder =3
                        FcstArea=int(table2["VOLUME"].values[FcstIndex]/NumberofFilesInaDataFolder)
                        
                        ObsIndex=int(table2[table2["OBJECT_CAT"]==ClusterNames[1]].index[0])
                        Obsi10=table2["INTENSITY_10"].values[ObsIndex]
                        Obsi25=table2["INTENSITY_25"].values[ObsIndex]
                        Obsi50=table2["INTENSITY_50"].values[ObsIndex]
                        Obsi75=table2["INTENSITY_75"].values[ObsIndex]
                        Obsi90=table2["INTENSITY_90"].values[ObsIndex]
                        ObsLoc=list()
                        ObsLoc.append([table2["CENTROID_X"].values[ObsIndex],table2["CENTROID_Y"].values[ObsIndex]])
                        ObsSpeedX=table2["X_DOT"].values[ObsIndex]
                        ObsSpeedY=table2["Y_DOT"].values[ObsIndex]
                        NumberofFilesInaDataFolder =3
                        ObsArea=int(table2["VOLUME"].values[ObsIndex]/NumberofFilesInaDataFolder)
                        FileCorruptFinder=0
                        break
                       
                    else:                 
                        FileCorruptFinder=1
                        # it means that it might exist data from at least one of products, in the future I might handle it.   
                except:

                    f = open('/home/z5194283/hdrive/MET_Tutorial/MyData/RealData/FinalData/log_file.txt', 'a')
                    f.write(filezz+"\n")                    
                    f.close()
                    FileCorruptFinder=1
                    
            #Creating a Cluster Object based on Storm Class:
            if FileCorruptFinder == 0:
                IMERG.append(Storm(Time,t,r,FcstArea,FcstSpeedX,FcstSpeedY,Fcsti10,Fcsti25,Fcsti50,Fcsti75,Fcsti90,FcstLoc))
                MRMS.append(Storm(Time,t,r,ObsArea,ObsSpeedX,ObsSpeedY,Obsi10,Obsi25,Obsi50,Obsi75,Obsi90,ObsLoc))
                log="Timestep="+str(IMERG[j].timestep)+" I90 = "+str(IMERG[j].i90)+" Locx = " + str(IMERG[j].loc[0][0]) +" Locy = " + str(IMERG[j].loc[0][1])  
                j=j+1
                f = open('/home/z5194283/hdrive/MET_Tutorial/MyData/RealData/FinalData/log_file.txt', 'a')
                f.write(log+"\n")
                f.close()
            else:
                IMERG.append(Storm(Time,t,r,0,0,0,0,0,0,0,0,FcstLoc))
                MRMS.append(Storm(Time,t,r,0,0,0,0,0,0,0,0,ObsLoc))
                log="Timestep="+str(IMERG[j].timestep)+" I90 = "+str(IMERG[j].i90)+" Locx = " + str(IMERG[j].loc[0][0]) +" Locy = " + str(IMERG[j].loc[0][1])
                j=j+1
                f = open('/home/z5194283/hdrive/MET_Tutorial/MyData/RealData/FinalData/log_file.txt', 'a')
                f.write(log+"\n")
                f.close()
                    
    Time=Time+1
    
    
print("plotting ...")    
# plotting    
MI90=list()
II90 =list()
TT=list()
II75=list()
MI75=list()
II50=list()
MI50=list()
MA=list()
IA=list()
MV=list()
IV=list()

for jj in range(0, len(IMERG)-1):
    if IMERG[jj].cov_th == 4 and IMERG[jj].radius_th == 6:
        IV.append(math.sqrt(IMERG[jj].speedx**2+IMERG[jj].speedy**2)*20)
        MV.append(math.sqrt(MRMS[jj].speedx**2+MRMS[jj].speedy**2)*20)
        TT.append(IMERG[jj].timestep)
        MI75.append(MRMS[jj].i75)
        II75.append(IMERG[jj].i75)
        MI90.append(MRMS[jj].i90)
        II90.append(IMERG[jj].i90)
        MI50.append(MRMS[jj].i50)
        II50.append(IMERG[jj].i50)
        MA.append(MRMS[jj].area)
        IA.append(IMERG[jj].area)

#plt.subplot(3,1,1)
#plt.plot(TT, II90, 'C3', zorder=1, lw=2, label='IMERG')
#plt.plot(TT, MI90, 'C2', zorder=2, lw=2, label='MRMS')
#plt.title('Intensity 90')
#plt.legend(loc='lower right')
##plt.xlabel("Time Step (30 min)")  
#plt.ylabel("Intensity (mm/hr)")     
#plt.tight_layout()  
#
#
#plt.subplot(3,1,2)
#plt.plot(TT, II75, 'C3', zorder=1, lw=2, label='IMERG')
#plt.plot(TT, MI75, 'C2', zorder=2, lw=2, label='MRMS')
#plt.title('Intensity 75')
#plt.legend(loc='lower right')
##plt.xlabel("Time Step (30 min)")  
#plt.ylabel("Intensity (mm/hr)") 
#plt.tight_layout()  
#
#plt.subplot(3,1,3)
#plt.plot(TT, II50, 'C3', zorder=1, lw=2, label='IMERG')
#plt.plot(TT, MI50, 'C2', zorder=2, lw=2, label='MRMS')
#plt.title('Intensity 50')
#plt.legend(loc='lower right')
#plt.xlabel("Time Step (30 min)")  
#plt.ylabel("Intensity (mm/hr)") 
#plt.tight_layout() 


#plt.subplot(2,1,1)
#plt.plot(TT, IA, 'C3', zorder=1, lw=2, label='IMERG')
#plt.plot(TT, MA, 'C2', zorder=2, lw=2, label='MRMS')
#plt.title('area')
#plt.legend(loc='lower right')
#plt.xlabel("Time Step (30 min)")  
#plt.ylabel("Area (Pixel=10 Km))")
#plt.tight_layout() 
#
#plt.subplot(2,1,2)
#plt.plot(TT, IV, 'C3', zorder=1, lw=2, label='IMERG')
#plt.plot(TT, MV, 'C2', zorder=2, lw=2, label='MRMS')
#plt.title('Velocity')
#plt.legend(loc='lower right')
#plt.xlabel("Time Step (30 min)")  
#plt.ylabel("Velocity (Km/hr)")
#plt.tight_layout() 

plt.scatter(MA, IA, s=20)
plt.title('Radar Vs IMERG')
#plt.ylim((13,30))
#plt.xlim((13,30))
plt.legend(loc=2)
plt.xlabel("Radar (100 Km^2)")  
plt.ylabel("IMERG (100 Km^2)")
plt.tight_layout() 

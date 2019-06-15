import glob
import subprocess
import pandas as pd 
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import os
import shutil
import fileinput
import math
from netCDF4 import Dataset
from scipy import stats
from scipy.interpolate import interp1d
import numpy.ma as ma
from operator import truediv

#

class Pixel:
  def __init__(self, ID, Bias,Intensity,Quantile,QualityIndex,RBias):
      self.ID = ID
      self.Bias = Bias
      self.Intensity=Intensity
      self.Quantile=Quantile
      self.QualityIndex=QualityIndex
      self.RBias=RBias

# Netcdf file mode outputs for PC:
#outputfolder="/home/ho0man/Temp/modeoutputs/ModeOutPut-Modified/"
# Netcdf file mode outputs for Storms:
outputfolder="/srv/ccrc/data60/z5194283/OutPuts/IMERG-MRMS-Larger-Domain-Mode-Run/"
outpufiles=glob.glob(outputfolder+"*_obj.nc")
outpufiles.sort() 
# TXT mode outputs
txtoutputfiles=glob.glob(outputfolder+"*_obj.txt")
txtoutputfiles.sort()
## input IMERG and MRMS files for PC:
#folder="/home/ho0man/Temp/modeoutputs/Sample/"
# input IMERG and MRMS files for Storms:
folder="/srv/ccrc/data60/z5194283/Data/LargerDomain/"
IMERG_files=glob.glob(folder+"IMERG_HHFV06B/Selected/*.nc")
IMERG_files.sort()  
MRMS_files=glob.glob(folder+"MRMS-RegriddedLargeDomain-6Days-30min/Selected/*.nc")
MRMS_files.sort()  
########################
#Filtering the Hurricane by using Mask created by MODE

i=0
B=list()
S=list()
I=list()
rows=0
cols=0

# For Quantile Approach
Pixellist=list()
I90_MRMS=list()
I90_IMERG=list()
I99_IMERG=list()
I50_MRMS=list()
I50_IMERG=list()
I10_MRMS=list()
I10_IMERG=list()
I90_IMERG_SensorID=list()
I99_IMERG_SensorID=list()
I99_locations=list()
IQ_IMERG=list()
IQ_MRMS=list()
IQ_IMERG_SensorID=list()
MEQI99_MRMS=list()

IMERG=list()
MRMS=list()
Sensor=list()
IMERGCentroidX=list()
IMERGCentroidY=list()
MRMSCentroidX=list()
MRMSCentroidY=list()

ACC_IMERG=list()
ACC_MRMS=list()
ACC_Sensor=list()

# Defining a matrix for recording number of evets in each pixel in order to take the best pixels:
NoE=np.zeros((148, 183))

# For Pixel by Pixel Approach
Pixellist0=list()

for file in  outpufiles:
    print(str(i))
    
    # extracting centorid of each cluster
    filez=txtoutputfiles[i]    
    table = pd.read_table(filez, delim_whitespace=True)  
#    print(table["OBJECT_ID"].values[-2])
    ICX = table["CENTROID_X"].values[-3]
    MCX = table["CENTROID_X"].values[-2]
    MRMSCentroidX.append(MCX)
    IMERGCentroidX.append(ICX)
    ICY = table["CENTROID_Y"].values[-3]
    MCY = table["CENTROID_Y"].values[-2]
    MRMSCentroidY.append(MCY)
    IMERGCentroidY.append(ICY)    
    Modeoutput = Dataset(file, 'r')    
    ######Below is the mask of hurricane:    
    IMERGMask=Modeoutput.variables['fcst_clus_id'][:,:]
    
    # Isolating the main cluster 
    IMERGMask[IMERGMask<0]=0
    IMERGMask[IMERGMask>1]=0
        
#    ######Below is the Mask of other objects excluding Hurricane
#    IMERGMask[IMERGMask<1]=-1
#    IMERGMask[IMERGMask>0]=0
#    IMERGMask[IMERGMask<0]=1
#    
#    ######Below is the mask of all places
#    IMERGMask[IMERGMask<1]=1

    MRMSMask=Modeoutput.variables['obs_clus_id'][:,:]
    MRMSMask[MRMSMask<0]=0
    MRMSMask[MRMSMask>1]=0
    IMERGfile = Dataset(IMERG_files[i], 'r')
    SensorType=IMERGfile.variables['HQprecipSource'][0,:,:]
    SensorType[SensorType<0]=0
    IMERGData=IMERGfile.variables['precipitationCal'][0,:,:] 
    IMERGData[IMERGData<0]=0
    QIData=IMERGfile.variables['precipitationQualityIndex'][0,:,:] 
    QIData[QIData<0]=0
#    IRonly=IMERGfile.variables['IRprecipitation'][:,:] 
#    IRonly[IRonly<0]=0
    MRMSfile = Dataset(MRMS_files[i], 'r')
    MRMSData=MRMSfile.variables['PrecipRate_0mabovemeansealevel'][0,:,:]
    MRMSData[MRMSData<0]=0

    MRMSFiltered=np.multiply(MRMSData, MRMSMask)   
    IMERGFiltered=np.multiply(np.transpose(IMERGData), IMERGMask)     
    SensorTypeFiltered=np.multiply(np.transpose(SensorType), IMERGMask) 
    QIFiltered=np.multiply(np.transpose(QIData), IMERGMask)
    QI1=ma.getdata(QIFiltered)
    S1=ma.getdata(SensorTypeFiltered)
    Sensor.append(S1)
    I1=ma.getdata(IMERGFiltered) 
    IMERG.append(I1)
    M1=ma.getdata(MRMSFiltered)
    MRMS.append(M1)   
    
    
    ### Accumulated Precip. for IMERG and MRMS ... Sensor is the highest value 
    ACC_IMERG.append(np.sum(I1)/np.size(I1))
    ACC_MRMS.append(np.sum(M1)/np.size(M1))
    ACC_Sensor.append(np.max(S1))

    
    
               
    # reshaping 2d to 1d and filtering zero values
    M=list()
    rows=len(M1)
    cols=len(M1[0])
    for ii in range(rows):
        for jj in range(cols):
            if M1[ii][jj] > 0:
                M.append(M1[ii][jj])
      
    I=list()
    B=list()
    S=list()
    QI=list()
    rows=len(I1)
    cols=len(I1[0])
    Pixellist00=list()
    for ii in range(rows):
        for jj in range(cols):
            if I1[ii][jj] > 0:  
                
                

                I.append(I1[ii][jj])
                S.append(S1[ii][jj])
                B.append(M1[ii][jj]-I1[ii][jj])
                QI.append(QI1[ii][jj])
                Pixellist00.append(Pixel(S1[ii][jj],M1[ii][jj]-I1[ii][jj],I1[ii][jj],1,QI1[ii][jj],(M1[ii][jj]-I1[ii][jj])/I1[ii][jj]))
         
    ##### this matrix show the number of events in each pixel (NOE) 
#    for ii in range(rows):
#        for jj in range(cols):        
#            if I1[ii][jj] > np.quantile(I,0.5):
#                NoE[ii][jj]=NoE[ii][jj]+1        
    
    
    
    # Combining the IMERG data and Sensors ID together
#    IS_Combined=list()
#    IS_Combined= np.vstack((S, I))
#    
#    # The process of sorting and ranking of IMERG and MRMS product
#    I_df = pd.DataFrame(IS_Combined[1], columns=['a'])
#    
#    # sort it by the desired series and caculate the percentile
#    sI_df = I_df.sort_values('a').reset_index()
#    sI_df['b'] = sI_df.index / float(len(sI_df) - 1)
#    
#    # setup the interpolator using the value as the index
#    I_interp = interp1d(sI_df['a'], sI_df['b'])      
#    M_df = pd.DataFrame(M, columns=['a'])
#    
#    # sort it by the desired series and caculate the percentile
#    sM_df = M_df.sort_values('a').reset_index()
#    sM_df['b'] = sM_df.index / float(len(sM_df) - 1)
#    
#    # setup the interpolator using the value as the index
#    M_interp = interp1d(sM_df['a'], sM_df['b'])   
#    
#    # finding the third smallest and biggest value of MRMS to filter outliers
#    s = set(M) # used to convert any of the list/tuple to the distinct element and sorted sequence of elements
#    
#    # Note: above statement will convert list into sets 
#    M_3rdSmallest=sorted(s)[1]
#    M_3rdBiggest=sorted(s)[-2]
#    
#    # Matching Process:
#    j=0
#    stop=0
#    filefinder=list()
#    Pixellistt=list()
#    for j in range(len(IS_Combined[1])):
#        I_quantile=I_interp(IS_Combined[1][j])        
#        EQ_MRMS=np.quantile(M,I_quantile)
#        if I_quantile<=0.999 and I_quantile>=0.95:# and EQ_MRMS>=M_3rdSmallest and EQ_MRMS<=M_3rdBiggest:         
#            Biass=np.subtract(EQ_MRMS,IS_Combined[1][j])
#            RBiass=np.divide(Biass,IS_Combined[1][j])
##            if IS_Combined[0][j]==9:
##                stop=stop+1
##                filefinder.append(i)
#            # Adding all information of a pixel to Pixel class in order to find all related information to gether:
#            Pixellistt.append(Pixel(IS_Combined[0][j],Biass,IS_Combined[1][j],I_quantile,QI[j],RBiass))
#    print("Reading files No. "+str(i))

    
    #finding the 90th percentile of precipitation from IMERG, MRMS and Sensor ID:
   
    
    # Finding the closest IMERG value to 99th percentile:
#    IMERGCloseto99=min(I, key=lambda x:abs(x-np.quantile(I,1)))
#    I99_IMERG.append(IMERGCloseto99)
#    I99_IMERG_SensorID.append(max(IS_Combined[0][np.where(IS_Combined[1] == IMERGCloseto99)]))
#    I99_location=np.where(I1 == IMERGCloseto99)
#    I99_locations.append([I99_location[0][0],I99_location[1][0]])
#    
#    # Finding the MRMS value at the same location of I99
#    AllMEQI99_MRMS=M1[I99_location[0],I99_location[1]]
#    MEQI99_MRMS.append(min(AllMEQI99_MRMS))
    
    
#    IMERGCloseto50=min(I, key=lambda x:abs(x-np.quantile(I,0.50)))
#    I50_IMERG.append(IMERGCloseto50)
#    IMERGCloseto10=min(I, key=lambda x:abs(x-np.quantile(I,0.10)))
#    I10_IMERG.append(IMERGCloseto10)
#    # The same to MRMS
#    MRMSCloseto90=min(M, key=lambda x:abs(x-np.quantile(M,0.90)))
#    I90_MRMS.append(MRMSCloseto90)
#    MRMSCloseto50=min(I, key=lambda x:abs(x-np.quantile(M,0.50)))
#    I50_MRMS.append(MRMSCloseto50)
#    MRMSCloseto10=min(I, key=lambda x:abs(x-np.quantile(M,0.10)))
#    I10_MRMS.append(MRMSCloseto10)
    
   
    
    #######################################
#    percentile=[0.5,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,0.99]
##    percentile=[0.35,0.37,0.39,0.41,0.43,0.45,0.47,0.49,0.51,0.53,0.55,0.57,0.59,0.61,0.63,0.65,0.67,0.69,0.71,0.73]
##  
#    IQu_IMERG=list()
#    IQu_MRMS=list()
#    IQu_IMERG_SensorID=list()
#    for pe in percentile:
#        IMERGClosetoQu=min(I, key=lambda x:abs(x-np.quantile(I,pe)))
#        IQu_IMERG.append(IMERGClosetoQu)
#        MRMSClosetoQu=min(I, key=lambda x:abs(x-np.quantile(M,pe)))
#        IQu_MRMS.append(MRMSClosetoQu)
#        IQu_IMERG_SensorID.append(max(IS_Combined[0][np.where(IS_Combined[1] == IMERGClosetoQu)]))
#    IQ_IMERG.append(IQu_IMERG)
#    IQ_MRMS.append(IQu_MRMS)
#    IQ_IMERG_SensorID.append(IQu_IMERG_SensorID)
    ########################################    
    
    
    # Checking if there are two or more different PMW Sensors:
#    DifferentSensorChecker=list()
#    DifferentSensorChecker.append(np.unique(IS_Combined[0][np.where(IS_Combined[1] == IMERGCloseto99)]))
#
#    if len(DifferentSensorChecker)>2:
#        print("There are two sensors")
#                 
    i=i+1  
#    
#    if i>100:
#        break
#    Pixellist0.append(Pixellist00) 
#    Pixellist.append(Pixellistt) 
    
Timestep=np.arange(i) 
#######################
#
##Plotting the bar chart with percentiles: 
#PrecentileID=1   
#plt.subplot(2,1,1)
#Sen_ID=[QQ[15] for QQ in IQ_IMERG_SensorID]
#plt.bar(Timestep,Sen_ID)
#plt.title('PMW Sensors')
##plt.legend(loc='lower right') 
#plt.ylabel("Sensor ID")  
##plt.xlabel("Intensity (mm/hr)")  
#plt.yticks(np.arange(0, 12, step=1)) 
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  
#
#plt.subplot(2,1,2)
#IME=[QQ[15] for QQ in IQ_IMERG]
#plt.plot(Timestep,IME,label='IMERG-'+str(percentile[9]))#, lw=3)
##plt.plot(Timestep,I50_IMERG,label='IMERG-95')
##plt.plot(Timestep,I10_IMERG,label='IMERG-90')
#
#MRM=[QQ[15] for QQ in IQ_MRMS]
#plt.plot(Timestep,MRM,label='MRMS-'+str(percentile[9]))
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='lower right') 
#plt.ylabel("Intensity (mm/hr)")   
##plt.yticks(np.arange(-60, 100, step=20))  
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  
#
#plt.subplot(4,2,2)
#Sen_ID=[QQ[15] for QQ in IQ_IMERG_SensorID]
#plt.bar(Timestep,Sen_ID)
#plt.title('PMW Sensors')
##plt.legend(loc='lower right') 
#plt.ylabel("Sensor ID")  
##plt.xlabel("Intensity (mm/hr)")  
#plt.yticks(np.arange(0, 12, step=1)) 
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  
##
#plt.subplot(4,2,4)
#IME=[QQ[15] for QQ in IQ_IMERG]
#plt.plot(Timestep,IME,label='IMERG-'+str(percentile[15]))
##plt.plot(Timestep,I50_IMERG,label='IMERG-95')
##plt.plot(Timestep,I10_IMERG,label='IMERG-90')
#MRM=[QQ[15] for QQ in IQ_MRMS]
#plt.plot(Timestep,MRM,label='MRMS-'+str(percentile[15]))
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='lower right') 
#plt.ylabel("Intensity (mm/hr)")   
##plt.yticks(np.arange(-60, 100, step=20))  
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  
#
#plt.subplot(4,2,5)
#Sen_ID=[QQ[3] for QQ in IQ_IMERG_SensorID]
#plt.bar(Timestep,Sen_ID)
#plt.title('PMW Sensors')
##plt.legend(loc='lower right') 
#plt.ylabel("Sensor ID")  
##plt.xlabel("Intensity (mm/hr)")  
#plt.yticks(np.arange(0, 12, step=1)) 
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  
##
#plt.subplot(4,2,7)
#IME=[QQ[3] for QQ in IQ_IMERG]
#plt.plot(Timestep,IME,label='IMERG-'+str(percentile[3]))
##plt.plot(Timestep,I50_IMERG,label='IMERG-95')
##plt.plot(Timestep,I10_IMERG,label='IMERG-90')
#MRM=[QQ[3] for QQ in IQ_MRMS]
#plt.plot(Timestep,MRM,label='MRMS-'+str(percentile[3]))
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='lower right') 
#plt.ylabel("Intensity (mm/hr)")   
##plt.yticks(np.arange(-60, 100, step=20))  
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  
#
#plt.subplot(4,2,6)
#Sen_ID=[QQ[3] for QQ in IQ_IMERG_SensorID]
#plt.bar(Timestep,Sen_ID)
#plt.title('PMW Sensors')
##plt.legend(loc='lower right') 
#plt.ylabel("Sensor ID")  
##plt.xlabel("Intensity (mm/hr)")  
#plt.yticks(np.arange(0, 12, step=1)) 
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  
##
#plt.subplot(4,2,8)
#IME=[QQ[0] for QQ in IQ_IMERG]
#plt.plot(Timestep,IME,label='IMERG-'+str(percentile[0]))
##plt.plot(Timestep,I50_IMERG,label='IMERG-95')
##plt.plot(Timestep,I10_IMERG,label='IMERG-90')
#MRM=[QQ[0] for QQ in IQ_MRMS]
#plt.plot(Timestep,MRM,label='MRMS-'+str(percentile[0]))
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='lower right') 
#plt.ylabel("Intensity (mm/hr)")   
##plt.yticks(np.arange(-60, 100, step=20))  
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  
##

#
#listRBT=list()
#listRB0=list()
#listRB5=list()
#listRB7=list()
#listRB11=list()
#listRB9=list()
#listRB3=list()
#
#listBT=list()
#listB0=list()
#listB5=list()
#listB7=list()
#listB11=list()
#listB9=list()
#listB3=list()
#
#listIT=list()
#listI0=list()
#listI5=list()
#listI7=list()
#listI11=list()
#listI9=list()
#listI3=list()
#
#listQT=list()
#listQ5=list()
#listQ0=list()
#listQ7=list()
#listQ9=list()
#listQ11=list()
#listQ3=list()
#
#listQIT=list()
#listQI5=list()
#listQI0=list()
#listQI7=list()
#listQI9=list()
#listQI11=list()
#listQI3=list()

#i=0
#for member in B:
#    if S[i]==0:
#        B0.append(B[i])
#        I0.append(I[i])
#    if S[i]==5:
#        B5.append(B[i])
#        I5.append(I[i])
#    if S[i]==7:
#        B7.append(B[i]) 
#        I7.append(I[i])
#    if S[i]==11:
#        B11.append(B[i])
#        I11.append(I[i])
#    if S[i]==9:
#        B9.append(B[i])  
#        I9.append(I[i])
#    i=i+1
#    print("Creating Bias lists No. "+str(i))

#k=0
### For just Pixel by Pixel plotting ... it should be commented when you want to plot quantile outputs
#for eachmemebrofPixellist in Pixellist0:
### For just Quantile plotting ... it should be commented when you want to plot Pixel by Pixel outputs
##for eachmemebrofPixellist in Pixellist:
#    
#    RBT=list()
#    RB0=list()
#    RB5=list()
#    RB7=list()
#    RB11=list()
#    RB9=list()
#    RB3=list()
#    
#    BT=list()
#    B0=list()
#    B5=list()
#    B7=list()
#    B11=list()
#    B9=list()
#    B3=list()
#    
#    IT=list()
#    I0=list()
#    I5=list()
#    I7=list()
#    I11=list()
#    I9=list()
#    I3=list()
#    
#    QT=list()
#    Q5=list()
#    Q0=list()
#    Q7=list()
#    Q9=list()
#    Q11=list()
#    Q3=list()
#    
#    QIT=list()
#    QI5=list()
#    QI0=list()
#    QI7=list()
#    QI9=list()
#    QI11=list()
#    QI3=list()
#    
#
#    for member in eachmemebrofPixellist:
#        RBT.append(member.RBias)
#        BT.append(member.Bias)
#        IT.append(member.Intensity)
#        QT.append(member.Quantile)
#        QIT.append(member.QualityIndex)    
#        if member.ID==0:
#            RB0.append(member.RBias)
#            B0.append(member.Bias)
#            I0.append(member.Intensity)
#            Q0.append(member.Quantile)
#            QI0.append(member.QualityIndex)
#        if member.ID==5:
#            RB5.append(member.RBias)
#            B5.append(member.Bias)
#            I5.append(member.Intensity)
#            Q5.append(member.Quantile)
#            QI5.append(member.QualityIndex)
#        if member.ID==7:
#            RB7.append(member.RBias)
#            B7.append(member.Bias)  
#            I7.append(member.Intensity)
#            Q7.append(member.Quantile)
#            QI7.append(member.QualityIndex)
#        if member.ID==11:
#            RB11.append(member.RBias)
#            B11.append(member.Bias)
#            I11.append(member.Intensity)
#            Q11.append(member.Quantile)
#            QI11.append(member.QualityIndex)
#        if member.ID==9:
#            RB9.append(member.RBias)
#            B9.append(member.Bias)  
#            I9.append(member.Intensity)
#            Q9.append(member.Quantile)
#            QI9.append(member.QualityIndex)
#        if member.ID==3:
#            RB3.append(member.RBias)
#            B3.append(member.Bias)  
#            I3.append(member.Intensity)
#            Q3.append(member.Quantile)
#            QI3.append(member.QualityIndex)
#        k=k+1
##        print("Creating Bias lists No. "+str(k))
#    
#    listRBT.append(RBT)
#    listRB0.append(RB0)
#    listRB5.append(RB5)
#    listRB7.append(RB7)
#    listRB11.append(RB11)
#    listRB9.append(RB9)
#    listRB3.append(RB3)
#        
#    listBT.append(BT)
#    listB0.append(B0)
#    listB5.append(B5)
#    listB7.append(B7)
#    listB11.append(B11)
#    listB9.append(B9)
#    listB3.append(B3)
#    
#    listIT.append(IT)
#    listI0.append(I0)
#    listI5.append(I5)
#    listI7.append(I7)
#    listI11.append(I11)
#    listI9.append(I9)
#    listI3.append(I3)
#    
#    listQT.append(QT)
#    listQ5.append(Q5)
#    listQ0.append(Q0)
#    listQ7.append(Q7)
#    listQ9.append(Q9)
#    listQ11.append(Q11)
#    listQ3.append(Q3)
#    
#    listQIT.append(QIT)
#    listQI5.append(QI5)
#    listQI0.append(QI0)
#    listQI7.append(QI7)
#    listQI9.append(QI9)
#    listQI11.append(QI11)
#    listQI3.append(QI3)
    

####Plotting the bar chart with Shifted Random Pixels: 
################################################################## 
#
#XX=111
#YY=126  
#
#
#
#IME=[Id[XX][YY] for Id in IMERG]
#Sen_ID=[Id[XX][YY] for Id in Sensor]
#MRM=list()
#ShiftedBias=list()
#Dist=list()
#for kkk in range(len(IME)):
#    ICX=IMERGCentroidX[kkk]
#    ICY=IMERGCentroidY[kkk]
#    MCX=MRMSCentroidX[kkk]
#    MCY=MRMSCentroidY[kkk]
#    DeltaX=int(MCX-ICX)
#    DeltaY=int(MCY-ICY)
#    Dist.append((DeltaX**2+DeltaY**2)**0.5)
#    IMEE=[Id[XX][YY] for Id in IMERG][kkk] 
##            print(IMEE)
#    MRMM=[Id[XX+DeltaX][YY+DeltaY] for Id in MRMS][kkk]  
##            print(MRMM)
#    ShiftedBias.append(MRMM-IMEE)       
#    MRM.append(MRMM)     
#
#
#plt.subplot(4,1,1)
#plt.bar(Timestep,Sen_ID)
#plt.title('PMW Sensors')
##plt.legend(loc='lower right') 
#plt.ylabel("Sensor ID")  
##plt.xlabel("Intensity (mm/hr)")  
#plt.yticks(np.arange(0, 12, step=1)) 
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  
#
#plt.subplot(4,1,2)
#
#plt.plot(Timestep,IME,label='IMERG')#, lw=3)
#
#
#plt.plot(Timestep,MRM,label='MRMS')
#
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='lower right') 
#plt.ylabel("Intensity (mm/hr)")   
##plt.yticks(np.arange(-60, 100, step=20))  
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  
##
#
#plt.subplot(4,1,3)
#IMEnp=np.asarray(IME, dtype=np.float32)
#MRMnp=np.asarray(MRM, dtype=np.float32)
#plt.plot(Timestep,MRMnp-IMEnp,label='Bias')#, lw=3)
#
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='lower right') 
#plt.ylabel("Intensity (mm/hr)")   
##plt.yticks(np.arange(-60, 100, step=20))  
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  
#
#plt.subplot(4,1,4)
#IMEnp=np.asarray(IME, dtype=np.float32)
#MRMnp=np.asarray(MRM, dtype=np.float32)
#plt.plot(Timestep,Dist,label='Centroid Distance (Pixel)')#, lw=3)
#
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='lower right') 
#plt.ylabel("Distance (Pixel)")   
##plt.yticks(np.arange(-60, 100, step=20))  
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  



#####Plotting the bar chart with Window Random Pixels: 
################################################################### 
#
#XX=111
#YY=126   
#
#
#WinDist=list()
#IME=[Id[XX][YY] for Id in IMERG]
#Sen_ID=[Id[XX][YY] for Id in Sensor]
#MRM=list()
#window=4
#n=window*2+1
#for kkk in range(len(IME)):
#    Biiaass=np.zeros((n,n))
#    for iii in range(-int(n/2),int(n/2)+1):
#        for jjj in range(-int(n/2),int(n/2)+1):
#
#            IMEE=[Id[XX][YY] for Id in IMERG][kkk] 
##            print(IMEE)
#            MRMM=[Id[XX+iii][YY+jjj] for Id in MRMS][kkk]  
##            print(MRMM)
#            Biiaass[iii+int(n/2)][jjj+int(n/2)]=MRMM-IMEE
#    XXM=XX
#    YYM=YY         
#    for iii in range(len(Biiaass)):
#        for jjj in range(len(Biiaass[0])):
#            if np.min(np.absolute(Biiaass))==np.absolute(Biiaass[iii][jjj]) and np.absolute(Biiaass[iii][jjj])!=0:
#                XXM=XX+(iii-int(n/2))
#                YYM=YY+(jjj-int(n/2))                
#             
#    MRM.append([Id[XXM][YYM] for Id in MRMS][kkk])   
#    WinDist.append(((XX-XXM)**2+(YY-YYM)**2)**0.5)
#
#
#plt.subplot(3,1,1)
#plt.bar(Timestep,Sen_ID)
#plt.title('PMW Sensors, X='+str(XX)+' Y='+str(YY)+', Window = '+str(n)+' Pixels')
##plt.legend(loc='lower right') 
#plt.ylabel("Sensor ID")  
##plt.xlabel("Intensity (mm/hr)")  
#plt.yticks(np.arange(0, 12, step=1)) 
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  
#
#plt.subplot(3,1,2)
#
#plt.plot(Timestep,IME,label='IMERG')#, lw=3)
#
#
#plt.plot(Timestep,MRM,label='MRMS')
#
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='lower right') 
#plt.ylabel("Intensity (mm/hr)")   
##plt.yticks(np.arange(-60, 100, step=20))  
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  
##
#
#plt.subplot(3,1,3)
#IMEnp=np.asarray(IME, dtype=np.float32)
#MRMnp=np.asarray(MRM, dtype=np.float32)
#plt.plot(Timestep,WinDist,label='Matching Distance (Pixel)')#, lw=3)
#
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='lower right') 
#plt.ylabel("Distance (Pixel)")   
##plt.yticks(np.arange(-60, 100, step=20))  
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  



###########Plotting the I99 Sensor chart with MRMS at the same location
#
#
#plt.subplot(3,1,1)
#plt.bar(Timestep,I99_IMERG_SensorID)
#plt.title('PMW Sensors, I99')
##plt.legend(loc='lower right') 
#plt.ylabel("Sensor ID")  
##plt.xlabel("Intensity (mm/hr)")  
#plt.yticks(np.arange(0, 12, step=1)) 
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  
#
#plt.subplot(3,1,2)
#
#plt.plot(Timestep,I99_IMERG,label='IMERG')#, lw=3)
#
#
#plt.plot(Timestep,MEQI99_MRMS,label='MRMS')
#
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='lower right') 
#plt.ylabel("Intensity (mm/hr)")   
##plt.yticks(np.arange(-60, 100, step=20))  
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  
##
#
#plt.subplot(3,1,3)
#IMEnp=np.asarray(I99_IMERG, dtype=np.float32)
#MRMnp=np.asarray(MEQI99_MRMS, dtype=np.float32)
#plt.plot(Timestep,MRMnp-IMEnp,label='Bias (MRMS-IMERG)')#, lw=3)
#
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='lower right') 
#plt.ylabel("Distance (Pixel)")   
##plt.yticks(np.arange(-60, 100, step=20))  
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  



######Plotting I99 sensor bar chart with a window
####################################
#     
#
#
##XX=111
##YY=130
#Sen_ID=list()
#IME=list()
#WinDist=list()
#cntr=0
#for Id in IMERG:
#    Sen_ID.append(Sensor[cntr][I99_locations[cntr][0]][I99_locations[cntr][1]])   
#    IME.append(Id[I99_locations[cntr][0]][I99_locations[cntr][1]])
#    cntr=cntr+1
#MRM=list()
#window=4
#n=window*2+1
#for kkk in range(len(IME)):
#    Biiaass=np.zeros((n,n))
#    XX=I99_locations[kkk][0]
#    YY=I99_locations[kkk][1]
#    for iii in range(-int(n/2),int(n/2)+1):
#        for jjj in range(-int(n/2),int(n/2)+1):
#
#            IMEE=IME[kkk] 
##            print(IMEE)
#            MRMM=[Id[XX+iii][YY+jjj] for Id in MRMS][kkk]  
##            print(MRMM)
#            Biiaass[iii+int(n/2)][jjj+int(n/2)]=MRMM-IMEE
#    XXM=XX
#    YYM=YY         
#    for iii in range(len(Biiaass)):
#        for jjj in range(len(Biiaass[0])):
#            if np.min(np.absolute(Biiaass))==np.absolute(Biiaass[iii][jjj]) and np.absolute(Biiaass[iii][jjj])!=0:
#                XXM=XX+(iii-int(n/2))
#                YYM=YY+(jjj-int(n/2))                
#             
#    MRM.append([Id[XXM][YYM] for Id in MRMS][kkk])   
#    WinDist.append(((XX-XXM)**2+(YY-YYM)**2)**0.5)
#
#
#plt.subplot(3,1,1)
#plt.bar(Timestep,Sen_ID)
#plt.title('PMW Sensors, X='+str(XX)+' Y='+str(YY)+', Window = '+str(n)+' Pixels')
##plt.legend(loc='lower right') 
#plt.ylabel("Sensor ID")  
##plt.xlabel("Intensity (mm/hr)")  
#plt.yticks(np.arange(0, 12, step=1)) 
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  
#
#plt.subplot(3,1,2)
#
#plt.plot(Timestep,IME,label='IMERG')#, lw=3)
#
#
#plt.plot(Timestep,MRM,label='MRMS')
#
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='lower right') 
#plt.ylabel("Intensity (mm/hr)")   
##plt.yticks(np.arange(-60, 100, step=20))  
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  
##
#
#plt.subplot(3,1,3)
#IMEnp=np.asarray(IME, dtype=np.float32)
#MRMnp=np.asarray(MRM, dtype=np.float32)
#plt.plot(Timestep,WinDist,label='Matching Distance (Pixel)')#, lw=3)
#
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='lower right') 
#plt.ylabel("Distance (Pixel)")   
##plt.yticks(np.arange(-60, 100, step=20))  
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  
#




######Plotting the bar chart with Random Pixels with MRMS and IMERG IMAGES: 
################################################################### 

    

#i=0
#
#for file in outpufiles:
#    IMERGfile = Dataset(IMERG_files[i], 'r')
#    SensorType=IMERGfile.variables['HQprecipSource'][:,:]
#    SensorType[SensorType<0]=0
#    
#    Modeoutput = Dataset(file, 'r')
#    IMERGMask=Modeoutput.variables['fcst_obj_raw'][:,:]    
#    MRMSMask=Modeoutput.variables['obs_obj_raw'][:,:]
#    plt.subplot(2,2,1)
#    plt.title('I'+str(i))
#    plt.imshow(np.transpose(SensorType), cmap='Greys')
#    plt.imshow(IMERGMask, vmin=0, vmax=70)
#    #plt.scatter(126,111,color='red',s=2)
#    
#    plt.xticks([])
#    plt.yticks([])
#    plt.subplot(2,2,2)
#    plt.title('M'+str(i))
#    plt.imshow(MRMSMask, vmin=0, vmax=70)
#    #plt.scatter(126,111,color='red',s=2)
#    
#    plt.xticks([])
#    plt.yticks([])
#    
#    mng = plt.get_current_fig_manager()
#    mng.window.showMaximized()
#    plt.tight_layout()
#    plt.savefig('/srv/ccrc/data60/z5194283/OutPuts/2DSensorsImages/'+str(i)+'.png',dpi=600)
#    plt.close()
#    i=i+1




#plt.subplot(2,1,2)
#IME=[Id[XX][YY] for Id in IMERG]
#plt.plot(Timestep,IME,label='IMERG')#, lw=3)
#
#MRM=[Id[XX][YY] for Id in MRMS]
#plt.plot(Timestep,MRM,label='MRMS')
#
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='lower right') 
#plt.ylabel("Intensity (mm/hr)")   
##plt.yticks(np.arange(-60, 100, step=20))  
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  
#
    


#
    
##plotting bar chart with Accumulated precipitaiton:
################################    
#plt.subplot(3,1,1)
#
#
#plt.bar(Timestep,ACC_Sensor)
#plt.title('PMW Sensors Accumlated')
##plt.legend(loc='lower right') 
#plt.ylabel("Sensor ID")  
##plt.xlabel("Intensity (mm/hr)")  
#plt.yticks(np.arange(0, 12, step=1)) 
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  
#
#
#plt.subplot(3,1,2)
#
#plt.plot(Timestep,ACC_IMERG,label='IMERG')#, lw=3)
#
#
#plt.plot(Timestep,ACC_MRMS,label='MRMS')
#
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='lower right') 
#plt.ylabel("Intensity (mm/hr)")   
##plt.yticks(np.arange(-60, 100, step=20))  
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  
##
#
#plt.subplot(3,1,3)
#IMEnp=np.asarray(ACC_IMERG, dtype=np.float32)
#MRMnp=np.asarray(ACC_MRMS, dtype=np.float32)
#plt.plot(Timestep,MRMnp-IMEnp,label='Bias MRMS-IMERG')#, lw=3)
#
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='lower right') 
#plt.ylabel("Intensity (mm/hr)")   
##plt.yticks(np.arange(-60, 100, step=20))  
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  

    
    
#######Plotting the bar chart with Random Pixels: 
#################################################################### 
#plt.subplot(3,1,1)
#XX=113
#YY=127
#Sen_ID=[Id[XX][YY] for Id in Sensor]
#plt.bar(Timestep,Sen_ID)
#plt.title('PMW Sensors X='+str(XX)+' Y='+str(YY))
##plt.legend(loc='lower right') 
#plt.ylabel("Sensor ID")  
##plt.xlabel("Intensity (mm/hr)")  
#plt.yticks(np.arange(0, 12, step=1)) 
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  
#
#
#plt.subplot(3,1,2)
#IME=[Id[XX][YY] for Id in IMERG]
#plt.plot(Timestep,IME,label='IMERG')#, lw=3)
#
#MRM=[Id[XX][YY] for Id in MRMS]
#plt.plot(Timestep,MRM,label='MRMS')
#
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='lower right') 
#plt.ylabel("Intensity (mm/hr)")   
##plt.yticks(np.arange(-60, 100, step=20))  
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  
##
#
#plt.subplot(3,1,3)
#IMEnp=np.asarray(IME, dtype=np.float32)
#MRMnp=np.asarray(MRM, dtype=np.float32)
#plt.plot(Timestep,MRMnp-IMEnp,label='Bias')#, lw=3)
#
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='lower right') 
#plt.ylabel("Intensity (mm/hr)")   
##plt.yticks(np.arange(-60, 100, step=20))  
#plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  



##################################### 
####Scatter Plot of random pixels:
####################################################################
    
#pp=list()
#for ii in range(rows):
#    print('one row has completed')
#    for jj in range(cols):        
#        if NoE[ii][jj] > 0:
#            pp.append([ii,jj])
#            
            
    
    
    
    
#plt.subplot(2,2,1)
#XX=95
#YY=126
            
####### Plotting scatters for each with sensor classification:
#           
IMEnp=list()
MRMnp=list()
BiasIM=list()   

IMEnpS0=list()
MRMnpS0=list()
BiasIMS0=list()  
NEventS0=list()  
IMEnpS3=list()
MRMnpS3=list()
BiasIMS3=list()   
NEventS3=list()  
IMEnpS5=list()
MRMnpS5=list()
BiasIMS5=list() 
NEventS5=list()    
IMEnpS7=list()
MRMnpS7=list()
BiasIMS7=list()   
NEventS7=list()  
IMEnpS9=list()
MRMnpS9=list()
BiasIMS9=list()  
NEventS9=list()   
IMEnpS11=list()
MRMnpS11=list()
BiasIMS11=list()    
NEventS11=list()   
IMEnpS20=list()
MRMnpS20=list()
BiasIMS20=list() 
NEventS20=list()  
IMEnpS=list()
MRMnpS=list()
BiasIMS=list() 
NEventS=list()  

for nnn in range(len(IMERG)):
    print(str(nnn))
    IMEnp0=np.asarray(IMERG[nnn], dtype=np.float32)
    MRMnp0=np.asarray(MRMS[nnn], dtype=np.float32)
    SNnp0=np.asarray(Sensor[nnn], dtype=np.float32)
    if nnn!=0:
        SNnplast=np.asarray(Sensor[nnn-1], dtype=np.float32)
        if nnn!=1:
            SNnplastlast=np.asarray(Sensor[nnn-2], dtype=np.float32)
    if nnn<152:
        SNnpfuture=np.asarray(Sensor[nnn+1], dtype=np.float32)
        if nnn<151:
            SNnpfuturefuture=np.asarray(Sensor[nnn+2], dtype=np.float32)
            
            
    for XX in range(len(IMEnp0)):
        for YY in range(len(IMEnp0[0])):      

            #Filtering Zero-Value pixels in both IMERG and MRMS and Filtering Sensros
            if nnn!=0 and nnn!=152:
                if IMEnp0[XX][YY]!=0 and MRMnp0[XX][YY]!=0 and SNnp0[XX][YY]==0 and SNnplast[XX][YY]==0 and SNnpfuture[XX][YY]==0:
                    IMEnpS0.append(IMEnp0[XX][YY])
                    MRMnpS0.append(MRMnp0[XX][YY])
                    BiasIMS0.append(IMEnp0[XX][YY]-MRMnp0[XX][YY])  
                    NEventS0.append(nnn)  

                if nnn>1 and nnn<151:
                
                    #Filtering Zero-Value pixels in both IMERG and MRMS and Filtering Sensros
                    if IMEnp0[XX][YY] !=0 and MRMnp0[XX][YY]!=0 and SNnp0[XX][YY]==3:# and SNnplast[XX][YY]==0 and SNnpfuture[XX][YY]==0 and SNnplastlast[XX][YY]==0 and SNnpfuturefuture[XX][YY]==0:
                        IMEnpS3.append(IMEnp0[XX][YY])
                        MRMnpS3.append(MRMnp0[XX][YY])
                        BiasIMS3.append(IMEnp0[XX][YY]-MRMnp0[XX][YY])   
                        NEventS3.append(nnn) 

                    
            
                    #Filtering Zero-Value pixels in both IMERG and MRMS and Filtering Sensros
                    if IMEnp0[XX][YY] !=0 and MRMnp0[XX][YY]!=0 and SNnp0[XX][YY]==5:# and SNnplast[XX][YY]==0 and SNnpfuture[XX][YY]==0 and SNnplastlast[XX][YY]==0 and SNnpfuturefuture[XX][YY]==0:
                        IMEnpS5.append(IMEnp0[XX][YY])
                        MRMnpS5.append(MRMnp0[XX][YY])
                        BiasIMS5.append(IMEnp0[XX][YY]-MRMnp0[XX][YY])  
                        NEventS5.append(nnn) 

            
            
                    #Filtering Zero-Value pixels in both IMERG and MRMS and Filtering Sensros
                    if IMEnp0[XX][YY] !=0 and MRMnp0[XX][YY]!=0 and SNnp0[XX][YY]==7:# and SNnplast[XX][YY]==0 and SNnpfuture[XX][YY]==0 and SNnplastlast[XX][YY]==0 and SNnpfuturefuture[XX][YY]==0:
                        IMEnpS7.append(IMEnp0[XX][YY])
                        MRMnpS7.append(MRMnp0[XX][YY])
                        BiasIMS7.append(IMEnp0[XX][YY]-MRMnp0[XX][YY])   
                        NEventS7.append(nnn) 

            
            
            
                    #Filtering Zero-Value pixels in both IMERG and MRMS and Filtering Sensros
                    if IMEnp0[XX][YY] !=0 and MRMnp0[XX][YY]!=0 and SNnp0[XX][YY]==9:# and SNnplast[XX][YY]==0 and SNnpfuture[XX][YY]==0 and SNnplastlast[XX][YY]==0 and SNnpfuturefuture[XX][YY]==0:
                        IMEnpS9.append(IMEnp0[XX][YY])
                        MRMnpS9.append(MRMnp0[XX][YY])
                        BiasIMS9.append(IMEnp0[XX][YY]-MRMnp0[XX][YY])  
                        NEventS9.append(nnn) 
 
            
            
                    #Filtering Zero-Value pixels in both IMERG and MRMS and Filtering Sensros
                    if IMEnp0[XX][YY] !=0 and MRMnp0[XX][YY]!=0 and SNnp0[XX][YY]==11:# and SNnplast[XX][YY]==0 and SNnpfuture[XX][YY]==0 and SNnplastlast[XX][YY]==0 and SNnpfuturefuture[XX][YY]==0:
                        IMEnpS11.append(IMEnp0[XX][YY])
                        MRMnpS11.append(MRMnp0[XX][YY])
                        BiasIMS11.append(IMEnp0[XX][YY]-MRMnp0[XX][YY])  
                        NEventS11.append(nnn) 


                    #Filtering Zero-Value pixels in both IMERG and MRMS and Filtering Sensros
                    if IMEnp0[XX][YY] !=0 and MRMnp0[XX][YY]!=0 and SNnp0[XX][YY]==20:# and SNnplast[XX][YY]==0 and SNnpfuture[XX][YY]==0 and SNnplastlast[XX][YY]==0 and SNnpfuturefuture[XX][YY]==0:
                        IMEnpS20.append(IMEnp0[XX][YY])
                        MRMnpS20.append(MRMnp0[XX][YY])
                        BiasIMS20.append(IMEnp0[XX][YY]-MRMnp0[XX][YY])  
                        NEventS20.append(nnn) 
#
IMEnpS.append(IMEnpS0)
IMEnpS.append(IMEnpS3)
IMEnpS.append(IMEnpS5)
IMEnpS.append(IMEnpS7)
IMEnpS.append(IMEnpS9)
IMEnpS.append(IMEnpS11)
IMEnpS.append(IMEnpS20)


MRMnpS.append(MRMnpS0)
MRMnpS.append(MRMnpS3)
MRMnpS.append(MRMnpS5)
MRMnpS.append(MRMnpS7)
MRMnpS.append(MRMnpS9)
MRMnpS.append(MRMnpS11)
MRMnpS.append(MRMnpS20)


BiasIMS.append(BiasIMS0)
BiasIMS.append(BiasIMS3)
BiasIMS.append(BiasIMS5)
BiasIMS.append(BiasIMS7)
BiasIMS.append(BiasIMS9)
BiasIMS.append(BiasIMS11)
BiasIMS.append(BiasIMS20)


NEventS.append(list(set(NEventS0)))
NEventS.append(list(set(NEventS3))) 
NEventS.append(list(set(NEventS5))) 
NEventS.append(list(set(NEventS7))) 
NEventS.append(list(set(NEventS9))) 
NEventS.append(list(set(NEventS11))) 
NEventS.append(list(set(NEventS20))) 



#
##
#
#







#plt.scatter(IMEnpS0,MRMnpS0,label='Sensor 0')#, s=1)#, lw=3)
#plt.scatter(IMEnpS3,MRMnpS3,label='Sensor 3')#, s=1)
#plt.scatter(IMEnpS5,MRMnpS5,label='Sensor 5')#, s=1)
#plt.scatter(IMEnpS7,MRMnpS7,label='Sensor 7')#, s=1)
#plt.scatter(IMEnpS9,MRMnpS9,label='Sensor 9')#, s=1)
#plt.scatter(IMEnpS11,MRMnpS11,label='Sensor 11')#, s=1)
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='upper right') 
#plt.ylabel("MRMS Intensity (mm/hr)")# (IMERG-MRMS) (mm/hr)")   
#plt.xlabel("IMERG Intensity (mm/hr)")  
##plt.yticks(np.arange(-60, 100, step=20))  
##plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  

Sensorno=[0,3,5,7,9,11,20]
title='PMW Raw Data Hurricane Only Sensor No. '
############Scatterplotting for IMERG and MRMS
xmin=0
xmax=115
ymin=0
ymax=50
step=5
ylabel='IMERG Intensity (mm/hr)'
xlabel='MRMS Intensity (mm/hr)'
for pl in range(6):
    plt.subplot(2,3,pl+1)     
    plt.scatter(MRMnpS[pl],IMEnpS[pl],label='Sensor '+str(Sensorno[pl])+" NE="+str(len(NEventS[pl])))#, s=1)#, lw=3)
    plt.title(title+str(Sensorno[pl]))
    plt.legend(loc='upper right') 
    plt.ylabel(ylabel)# (IMERG-MRMS) (mm/hr)") 
    plt.xlabel(xlabel)  
    plt.axis('scaled')
    axes = plt.gca()
    axes.set_xlim([xmin,xmax])
    axes.set_ylim([ymin,ymax])
    plt.yticks(np.arange(ymin, ymax, step=step))  
    plt.xticks(np.arange(xmin, xmax, step=step))   
    plt.grid()
#    

###############Scatter Plotting MRMS and Bias
#xmin=0
#xmax=115
#ymin=-110
#ymax=50
#step=10
#ylabel='Bias IMERG-MRMS (mm/hr)'
#xlabel='MRMS Intensity (mm/hr)'
#for pl in range(6):
#    plt.subplot(2,3,pl+1)     
#    plt.scatter(MRMnpS[pl],BiasIMS[pl],label='Sensor '+str(Sensorno[pl])+" NE="+str(len(NEventS[pl])))#, s=1)#, lw=3)
#    plt.title(title+str(Sensorno[pl]))
#    plt.legend(loc='upper right') 
#    plt.ylabel(ylabel)# (IMERG-MRMS) (mm/hr)") 
#    plt.xlabel(xlabel)  
#    plt.axis('scaled')
#    axes = plt.gca()
#    axes.set_xlim([xmin,xmax])
#    axes.set_ylim([ymin,ymax])
#    plt.yticks(np.arange(ymin, ymax, step=step))  
#    plt.xticks(np.arange(xmin, xmax, step=step))   
#    plt.grid()




###############Scatter Plotting IMERG and Bias
#xmin=0
#xmax=115
#ymin=-110
#ymax=50
#step=10
#ylabel='Bias IMERG-MRMS (mm/hr)'
#xlabel='IMERG Intensity (mm/hr)'
#for pl in range(6):
#    plt.subplot(2,3,pl+1)     
#    plt.scatter(IMEnpS[pl],BiasIMS[pl],label='Sensor '+str(Sensorno[pl])+" NE="+str(len(NEventS[pl])))#, s=1)#, lw=3)
#    plt.title(title+str(Sensorno[pl]))
#    plt.legend(loc='upper right') 
#    plt.ylabel(ylabel)# (IMERG-MRMS) (mm/hr)") 
#    plt.xlabel(xlabel)  
#    plt.axis('scaled')
#    axes = plt.gca()
#    axes.set_xlim([xmin,xmax])
#    axes.set_ylim([ymin,ymax])
#    plt.yticks(np.arange(ymin, ymax, step=step))  
#    plt.xticks(np.arange(xmin, xmax, step=step))   
#    plt.grid()



#
#plt.subplot(2,3,1)     
#plt.scatter(MRMnpS0,IMEnpS0,label='Sensor 0')#, s=1)#, lw=3)
#plt.title('IR-Isolated 1 step')
#plt.legend(loc='upper right') 
#plt.ylabel(ylabel)# (IMERG-MRMS) (mm/hr)") 
#plt.xlabel(xlabel)  
#plt.axis('scaled')
#axes = plt.gca()
#axes.set_xlim([xmin,xmax])
#axes.set_ylim([ymin,ymax])
#plt.yticks(np.arange(ymin, ymax, step=step))  
#plt.xticks(np.arange(xmin, xmax, step=step))   
#plt.grid()
#
#plt.subplot(2,3,2)     
#plt.scatter(MRMnpS3,IMEnpS3,label='Sensor 3')#, s=1)
#plt.title('Sensor 3-Isolated 2 steps')
#plt.legend(loc='upper right') 
#plt.ylabel(ylabel)# (IMERG-MRMS) (mm/hr)") 
#plt.xlabel(xlabel)  
#plt.axis('scaled')
#axes = plt.gca()
#axes.set_xlim([xmin,xmax])
#axes.set_ylim([ymin,ymax])
#plt.yticks(np.arange(ymin, ymax, step=step))  
#plt.xticks(np.arange(xmin, xmax, step=step))   
#plt.grid()
#
#plt.subplot(2,3,3)     
#plt.scatter(MRMnpS5,IMEnpS5,label='Sensor 5')#, s=1)
#plt.title('Sensor 5-Isolated 2 steps')
#plt.legend(loc='upper right') 
#plt.ylabel(ylabel)# (IMERG-MRMS) (mm/hr)") 
#plt.xlabel(xlabel)  
#plt.axis('scaled')
#axes = plt.gca()
#axes.set_xlim([xmin,xmax])
#axes.set_ylim([ymin,ymax])
#plt.yticks(np.arange(ymin, ymax, step=step))  
#plt.xticks(np.arange(xmin, xmax, step=step))  
#plt.grid()
#
#plt.subplot(2,3,4)   
#plt.scatter(MRMnpS7,IMEnpS7,label='Sensor 7')#, s=1)  
#plt.title('Sensor 7-Isolated 2 steps')
#plt.legend(loc='upper right') 
#plt.ylabel(ylabel)# (IMERG-MRMS) (mm/hr)") 
#plt.xlabel(xlabel)  
#plt.axis('scaled')
#axes = plt.gca()
#axes.set_xlim([xmin,xmax])
#axes.set_ylim([ymin,ymax])
#plt.yticks(np.arange(ymin, ymax, step=step))  
#plt.xticks(np.arange(xmin, xmax, step=step))     
#plt.grid()
#
#plt.subplot(2,3,5)    
#plt.scatter(MRMnpS9,IMEnpS9,label='Sensor 9')#, s=1) 
#plt.title('Sensor 9-Isolated 2 steps')
#plt.legend(loc='upper right') 
#plt.ylabel(ylabel)# (IMERG-MRMS) (mm/hr)") 
#plt.xlabel(xlabel)  
#plt.axis('scaled')
#axes = plt.gca()
#axes.set_xlim([xmin,xmax])
#axes.set_ylim([ymin,ymax])
#plt.yticks(np.arange(ymin, ymax, step=step))  
#plt.xticks(np.arange(xmin, xmax, step=step))   
#plt.grid()
#
#plt.subplot(2,3,6)  
#plt.scatter(MRMnpS11,IMEnpS11,label='Sensor 11')#, s=1)   
#plt.title('Sensor 11-Isolated 2 steps')
#plt.legend(loc='upper right') 
#plt.ylabel(ylabel)# (IMERG-MRMS) (mm/hr)") 
#plt.xlabel(xlabel)  
#plt.axis('scaled')
#axes = plt.gca()
#axes.set_xlim([xmin,xmax])
#axes.set_ylim([ymin,ymax])
#plt.yticks(np.arange(ymin, ymax, step=step))  
#plt.xticks(np.arange(xmin, xmax, step=step))  
#
#plt.grid()
#
#plt.tight_layout()  
#
###################### histogram 1d plotting


##########Bias    
#for pl in range(6):
#    plt.subplot(2,3,pl+1) 
#    plt.hist(BiasIMS[0], bins=200, normed=True,label='IR, NE='+str(len(NEventS[0])),alpha=.4)
#    plt.hist(BiasIMS[pl], bins=200, normed=True,label='PMW, NE='+str(len(NEventS[pl])),alpha=.4)
#    plt.legend(loc='upper right')
#    plt.xlabel('Bias IMERG-MRMS (mm/hr)')
#    plt.ylabel('Percentile (0 to 1)')
#    plt.title(title+str(Sensorno[pl]))
#    plt.xticks(np.arange(-40, 60, step=10)) 
#    plt.xlim([-40,60])
    

##########MRMS and IMERG    
#for pl in range(6):
#    plt.subplot(2,3,pl+1)
#    plt.hist(IMEnpS[pl], bins=200, normed=True,label='IMERG',alpha=.4)
#    plt.hist(MRMnpS[pl], bins=200, normed=True,label='MRMS',alpha=.4)    
#    plt.legend(loc='upper right') 
#    plt.xlabel('Intensity (mm/hr)')
#    plt.ylabel('Percentile (0 to 1)')
#    plt.title(title+str(Sensorno[pl]))
#    plt.xticks(np.arange(0, 50, step=10)) 
#    plt.xlim([0,50])
    



#plt.subplot(2,3,1)         
#plt.hist(BiasIMS0, bins=200)
#plt.ylabel('No of times')
#plt.title('IR-Only-Isolated 1 step')
#plt.xticks(np.arange(-40, 60, step=10)) 
#plt.xlim([-40,60])
#
#plt.subplot(2,3,2)         
#plt.hist(BiasIMS3, bins=200)
#plt.ylabel('No of times')
#plt.title('Sensor 3')
#plt.xticks(np.arange(-40, 60, step=10)) 
#plt.xlim([-40,60])
#
#plt.subplot(2,3,3)         
#plt.hist(BiasIMS5, bins=200)
#plt.ylabel('No of times')
#plt.title('Sensor 5')
#plt.xticks(np.arange(-40, 60, step=10)) 
#plt.xlim([-40,60])
#
#plt.subplot(2,3,4)         
#plt.hist(BiasIMS7, bins=200)
#plt.ylabel('No of times')
#plt.title('Sensor 7')
#plt.xticks(np.arange(-40, 60, step=10)) 
#plt.xlim([-40,60])
#
#plt.subplot(2,3,5)         
#plt.hist(BiasIMS9, bins=200)
#plt.ylabel('No of times')
#plt.title('Sensor 9')
#plt.xticks(np.arange(-40, 60, step=10)) 
#plt.xlim([-40,60])
#plt.xlabel('Bias Distribution')
#
#plt.subplot(2,3,6)         
#plt.hist(BiasIMS11, bins=200)
#plt.ylabel('No of times')
#plt.title('Sensor 11')
#plt.xticks(np.arange(-40, 60, step=10)) 
#plt.xlim([-40,60])

#plt.tight_layout() 




######### histogram 2d plotting
    
    

################ MRMS IMERG
#xmin=0
#xmax=5
#ymin=0
#ymax=5
#
#ylabel='IMERG Intensity (mm/hr)'
#xlabel='MRMS Intensity (mm/hr)'
#
#for pl in range(6):
#    plt.subplot(2,3,pl+1)  
#    plt.hist2d(MRMnpS[pl], IMEnpS[pl], bins=int(np.max(MRMnpS[pl])/0.2), cmap='Reds', normed=True)
#
#    plt.title(title+str(Sensorno[pl]))
#    #plt.xticks(np.arange(-40, 60, step=10)) 
#    plt.axis('scaled')
#    plt.ylim([ymin,ymax])
#    plt.xlim([xmin,xmax])
#    plt.ylabel(ylabel)
#    plt.xlabel(xlabel)
#    cb = plt.colorbar()
#    cb.set_label('percentile 0 to 1, NE='+str(len(NEventS[pl])))
    
    
################# IMERG Bias
#xmin=0
#xmax=5
#ymin=-3
#ymax=5
#
#ylabel='Bias IMERG-MRMS (mm/hr)'
#xlabel='IMERG Intensity (mm/hr)'
#
#for pl in range(6):
#    plt.subplot(2,3,pl+1)  
#    plt.hist2d(IMEnpS[pl], BiasIMS[pl], bins=int(np.max(IMEnpS[pl])/0.2), cmap='Reds', normed=True)
#
#    plt.title(title+str(Sensorno[pl]))
#    #plt.xticks(np.arange(-40, 60, step=10)) 
#    plt.axis('scaled')
#    plt.ylim([ymin,ymax])
#    plt.xlim([xmin,xmax])
#    plt.ylabel(ylabel)
#    plt.xlabel(xlabel)
#    cb = plt.colorbar()
#    cb.set_label('percentile 0 to 1, NE='+str(len(NEventS[pl])))

################# MRMS Bias
#xmin=0
#xmax=40
#ymin=-15
#ymax=25
#
#ylabel='Bias IMERG-MRMS (mm/hr)'
#xlabel='MRMS Intensity (mm/hr)'

#for pl in range(6):
#    plt.subplot(2,3,pl+1)  
#    plt.hist2d(MRMnpS[pl], BiasIMS[pl], bins=int(np.max(MRMnpS[pl])/0.5), cmap='Reds', normed=True)
#
#    plt.title(title+str(Sensorno[pl]))
#    #plt.xticks(np.arange(-40, 60, step=10)) 
#    plt.axis('scaled')
#    plt.ylim([ymin,ymax])
#    plt.xlim([xmin,xmax])
#    plt.ylabel(ylabel)
#    plt.xlabel(xlabel)
#    cb = plt.colorbar()
#    cb.set_label('percentile 0 to 1')
#   
#xmin=0
#xmax=20
#ymin=0
#ymax=20   
#ylabel='IMERG Intensity (mm/hr)'
#    
#    
#plt.subplot(2,3,1)         
#plt.hist2d(MRMnpS0, IMEnpS0, bins=int(np.max(MRMnpS0)/0.5), cmap='Reds')
#plt.ylabel('No of times')
#plt.title('IR-Only-Isolated 1 step')
##plt.xticks(np.arange(-40, 60, step=10)) 
#plt.axis('scaled')
#plt.ylim([ymin,ymax])
#plt.xlim([xmin,xmax])
#plt.ylabel(ylabel)
#
#cb = plt.colorbar()
#cb.set_label('counts in bin')
#
#plt.subplot(2,3,2)         
#plt.hist2d(MRMnpS3, IMEnpS3, bins=int(np.max(MRMnpS3)/0.5), cmap='Reds')
#plt.ylabel('No of times')
#plt.title('Sensor 3-Isolated 2 steps')
##plt.xticks(np.arange(-40, 60, step=10))
#plt.axis('scaled')
#plt.ylim([ymin,ymax])
#plt.xlim([xmin,xmax])
#
#cb = plt.colorbar()
#cb.set_label('counts in bin')
#
#plt.subplot(2,3,3)         
#plt.hist2d(MRMnpS5, IMEnpS5, bins=int(np.max(MRMnpS5)/0.5), cmap='Reds')
#plt.ylabel('No of times')
#plt.title('Sensor 5-Isolated 2 steps')
##plt.xticks(np.arange(-40, 60, step=10)) 
#plt.axis('scaled')
#plt.ylim([ymin,ymax])
#plt.xlim([xmin,xmax])
#cb = plt.colorbar()
#cb.set_label('counts in bin')
#
#plt.subplot(2,3,4)         
#plt.hist2d(MRMnpS7, IMEnpS7, bins=int(np.max(MRMnpS7)/0.5), cmap='Reds')
#plt.ylabel('No of times')
#plt.title('Sensor 7-Isolated 2 steps')
##plt.xticks(np.arange(-40, 60, step=10)) 
#plt.axis('scaled')
#plt.ylim([ymin,ymax])
#plt.xlim([xmin,xmax])
#plt.ylabel(ylabel)
#cb = plt.colorbar()
#cb.set_label('counts in bin')
#
#plt.subplot(2,3,5)         
#plt.hist2d(MRMnpS9, IMEnpS9, bins=int(np.max(MRMnpS9)/0.5), cmap='Reds')
#plt.ylabel('No of times')
#plt.title('Sensor 9-Isolated 2 steps')
##plt.xticks(np.arange(-40, 60, step=10)) 
#plt.axis('scaled')
#plt.ylim([ymin,ymax])
#plt.xlim([xmin,xmax])
#plt.xlabel('MRMS-Intensity')
#
#cb = plt.colorbar()
#cb.set_label('counts in bin')
#
#plt.subplot(2,3,6)         
#plt.hist2d(MRMnpS11, IMEnpS11, bins=int(np.max(MRMnpS11)/0.5), cmap='Reds')
#plt.ylabel('No of times')
#plt.title('Sensor 11-Isolated 2 steps')
##plt.xticks(np.arange(-40, 60, step=10)) 
#plt.axis('scaled')
#plt.ylim([ymin,ymax])
#plt.xlim([xmin,xmax])
#cb = plt.colorbar()
#cb.set_label('counts in bin')
#
#plt.tight_layout() 





#######################################################################

#plt.subplot(2,2,2)
#XX=112
#YY=126
#
#IME=[Id[XX][YY] for Id in IMERG]
#MRM=[Id[XX][YY] for Id in MRMS]
#IMEnp0=np.asarray(IME, dtype=np.float32)
#MRMnp0=np.asarray(MRM, dtype=np.float32)
##Filtering Zero-Value pixels in both IMERG and MRMS 
#IMEnp=list()
#BiasIM=list()
#for i in range(len(IMEnp0)):
#    if IMEnp0[i] !=0 and MRMnp0[i]!=0:
#        IMEnp.append(IMEnp0[i])
#        BiasIM.append(MRMnp0[i]-IMEnp0[i])    
#plt.scatter(IMEnp,BiasIM,label='X='+str(XX)+', Y='+str(YY))#, lw=3)
#
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='upper right') 
#plt.ylabel("Bias (MRMS-IMERG) (mm/hr)")   
#plt.xlabel("Intensity (mm/hr)")  
##plt.yticks(np.arange(-60, 100, step=20))  
##plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
#plt.tight_layout()  
#
#plt.subplot(2,2,3)
#XX=110
#YY=132

#
#IME=[Id[XX][YY] for Id in IMERG]
#MRM=[Id[XX][YY] for Id in MRMS]
#IMEnp0=np.asarray(IME, dtype=np.float32)
#MRMnp0=np.asarray(MRM, dtype=np.float32)
##Filtering Zero-Value pixels in both IMERG and MRMS 
#IMEnp=list()
#BiasIM=list()
#for i in range(len(IMEnp0)):
#    if IMEnp0[i] !=0 and MRMnp0[i]!=0:
#        IMEnp.append(IMEnp0[i])
#        BiasIM.append(MRMnp0[i]-IMEnp0[i])    
#plt.scatter(IMEnp,BiasIM,label='X='+str(XX)+', Y='+str(YY))#, lw=3)
#
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='upper right') 
#plt.ylabel("Bias (MRMS-IMERG) (mm/hr)")   
#plt.xlabel("Intensity (mm/hr)")  
##plt.yticks(np.arange(-60, 100, step=20))  
##plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
##plt.tight_layout()  
#
#
#plt.subplot(2,2,4)
#XX=110
#YY=130
#
#IME=[Id[XX][YY] for Id in IMERG]
#MRM=[Id[XX][YY] for Id in MRMS]
#IMEnp0=np.asarray(IME, dtype=np.float32)
#MRMnp0=np.asarray(MRM, dtype=np.float32)
##Filtering Zero-Value pixels in both IMERG and MRMS 
#IMEnp=list()
#BiasIM=list()
#for i in range(len(IMEnp0)):
#    if IMEnp0[i] !=0 and MRMnp0[i]!=0:
#        IMEnp.append(IMEnp0[i])
#        BiasIM.append(MRMnp0[i]-IMEnp0[i])    
#plt.scatter(IMEnp,BiasIM,label='X='+str(XX)+', Y='+str(YY))#, lw=3)
#
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='upper right') 
#plt.ylabel("Bias (MRMS-IMERG) (mm/hr)")   
#plt.xlabel("Intensity (mm/hr)")  
##plt.yticks(np.arange(-60, 100, step=20))  
##plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
##plt.tight_layout()  
#

#
#plt.subplot(2,2,1)
#XX=111
#YY=125
#
#IME=[Id[XX][YY] for Id in IMERG]
#MRM=[Id[XX][YY] for Id in MRMS]
#IMEnp0=np.asarray(IME, dtype=np.float32)
#MRMnp0=np.asarray(MRM, dtype=np.float32)
##Filtering Zero-Value pixels in both IMERG and MRMS 
#IMEnp=list()
#BiasIM=list()
#for i in range(len(IMEnp0)):
#    if IMEnp0[i] !=0 and MRMnp0[i]!=0:
#        IMEnp.append(IMEnp0[i])
#        BiasIM.append(MRMnp0[i]-IMEnp0[i])    
#plt.scatter(IMEnp,BiasIM,label='X='+str(XX)+', Y='+str(YY))#, lw=3)
#
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='upper right') 
#plt.ylabel("Bias (MRMS-IMERG) (mm/hr)")   
#plt.xlabel("Intensity (mm/hr)")  
##plt.yticks(np.arange(-60, 100, step=20))  
##plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
##plt.tight_layout()  
#
#
#plt.subplot(2,2,2)
#XX=111
#YY=128
#
#IME=[Id[XX][YY] for Id in IMERG]
#MRM=[Id[XX][YY] for Id in MRMS]
#IMEnp0=np.asarray(IME, dtype=np.float32)
#MRMnp0=np.asarray(MRM, dtype=np.float32)
##Filtering Zero-Value pixels in both IMERG and MRMS 
#IMEnp=list()
#BiasIM=list()
#for i in range(len(IMEnp0)):
#    if IMEnp0[i] !=0 and MRMnp0[i]!=0:
#        IMEnp.append(IMEnp0[i])
#        BiasIM.append(MRMnp0[i]-IMEnp0[i])    
#plt.scatter(IMEnp,BiasIM,label='X='+str(XX)+', Y='+str(YY))#, lw=3)
#
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='upper right') 
#plt.ylabel("Bias (MRMS-IMERG) (mm/hr)")   
#plt.xlabel("Intensity (mm/hr)")  
##plt.yticks(np.arange(-60, 100, step=20))  
##plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
##plt.tight_layout()  
#
#
#
#plt.subplot(2,2,3)
#XX=112
#YY=128
#
#IME=[Id[XX][YY] for Id in IMERG]
#MRM=[Id[XX][YY] for Id in MRMS]
#IMEnp0=np.asarray(IME, dtype=np.float32)
#MRMnp0=np.asarray(MRM, dtype=np.float32)
##Filtering Zero-Value pixels in both IMERG and MRMS 
#IMEnp=list()
#BiasIM=list()
#for i in range(len(IMEnp0)):
#    if IMEnp0[i] !=0 and MRMnp0[i]!=0:
#        IMEnp.append(IMEnp0[i])
#        BiasIM.append(MRMnp0[i]-IMEnp0[i])    
#plt.scatter(IMEnp,BiasIM,label='X='+str(XX)+', Y='+str(YY))#, lw=3)
#
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='upper right') 
#plt.ylabel("Bias (MRMS-IMERG) (mm/hr)")   
#plt.xlabel("Intensity (mm/hr)")  
##plt.yticks(np.arange(-60, 100, step=20))  
##plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
##plt.tight_layout()  
#
#
#
#plt.subplot(2,2,4)
#XX=109
#YY=130
#
#IME=[Id[XX][YY] for Id in IMERG]
#MRM=[Id[XX][YY] for Id in MRMS]
#IMEnp0=np.asarray(IME, dtype=np.float32)
#MRMnp0=np.asarray(MRM, dtype=np.float32)
##Filtering Zero-Value pixels in both IMERG and MRMS 
#IMEnp=list()
#BiasIM=list()
#for i in range(len(IMEnp0)):
#    if IMEnp0[i] !=0 and MRMnp0[i]!=0:
#        IMEnp.append(IMEnp0[i])
#        BiasIM.append(MRMnp0[i]-IMEnp0[i])    
#plt.scatter(IMEnp,BiasIM,label='X='+str(XX)+', Y='+str(YY))#, lw=3)
#
##plt.title('Intensity 90th Precentile')
#plt.legend(loc='upper right') 
#plt.ylabel("Bias (MRMS-IMERG) (mm/hr)")   
#plt.xlabel("Intensity (mm/hr)")  
##plt.yticks(np.arange(-60, 100, step=20))  
##plt.xticks(np.arange(0, 74, step=1))  
#plt.grid()
##plt.tight_layout()  
####################################################################
#######Scatter Plotting
#plt.subplot(3,2,1)
#
#for n in range(len(listRB0)):    
#    plt.scatter(listI0[n],listRB0[n])    
#plt.title('RBias Sensor IR-Only')
#plt.legend(loc='lower right') 
#plt.ylabel("Bias mm/hr")  
#plt.xlabel("Intensity (mm/hr)")  
##plt.yticks(np.arange(-60, 100, step=20)) 
##plt.yticks(np.arange(-1, 3, step=0.5))  
#plt.grid()
#plt.tight_layout()  
#
#plt.subplot(3,2,2)   
#for n in range(len(listRB9)):    
#    plt.scatter(listI9[n],listRB9[n])    
#plt.title('RBias Sensor No. 9')
#plt.legend(loc='lower right') 
#plt.ylabel("Bias mm/hr")  
#plt.xlabel("Intensity (mm/hr)")  
##plt.yticks(np.arange(-60, 100, step=20)) 
##plt.yticks(np.arange(-1, 3, step=0.5))     
##plt.xticks(np.arange(0, 940, step=1))    
#plt.grid()
#plt.tight_layout()  
#
#plt.subplot(3,2,3)
#for n in range(len(listRB7)):    
#    plt.scatter(listI7[n],listRB7[n])    
#plt.title('RBias Sensor No. 7')
#plt.legend(loc='lower right') 
#plt.ylabel("Bias mm/hr") 
#plt.xlabel("Intensity (mm/hr)")  
##plt.yticks(np.arange(-60, 100, step=20)) 
##plt.yticks(np.arange(-1, 3, step=0.5))     
#plt.grid()
#plt.tight_layout()  
#
#plt.subplot(3,2,4)
#for n in range(len(listRB5)):    
#    plt.scatter(listI5[n],listRB5[n])    
#plt.title('RBias Sensor No. 5')
#plt.legend(loc='lower right') 
#plt.ylabel("Bias mm/hr") 
#plt.xlabel("Intensity (mm/hr)")  
##plt.yticks(np.arange(-60, 100, step=20)) 
##plt.yticks(np.arange(-1, 3, step=0.5))    
#plt.grid()
#plt.tight_layout()  
#
#plt.subplot(3,2,5)
#for n in range(len(listRB11)):    
#    plt.scatter(listI11[n],listRB11[n])    
#plt.title('RBias Sensor No. 11')
#plt.legend(loc='lower right') 
#plt.ylabel("Bias mm/hr") 
#plt.xlabel("Intensity (mm/hr)")  
##plt.yticks(np.arange(-60, 100, step=20)) 
##plt.yticks(np.arange(-1, 3, step=0.5))  
#plt.grid()
#plt.tight_layout()  
#
#plt.subplot(3,2,6)
#for n in range(len(listRB3)):    
#    plt.scatter(listI3[n],listRB3[n])    
#plt.title('RBias Sensor No. 3')
#plt.legend(loc='lower right') 
#plt.ylabel("Bias mm/hr") 
#plt.xlabel("Intensity (mm/hr)") 
#plt.ylim(-5, 10) 
##plt.yticks(np.arange(-10, 20, step=5)) 
##plt.yticks(np.arange(-1, 3, step=0.5))   
#plt.grid()
#plt.tight_layout()  
#


### Plotting Sensors fluctuations
#
#plt.subplot(3,2,1)
#TStep=0
#for n in range(len(listRB0)):    
#    plt.plot(np.arange(TStep,TStep+len(listRB0[n])),listRB0[n])
#    TStep=TStep+len(listRB0[n])
#plt.title('Bias Sensor IR-Only')
#plt.legend(loc='lower right') 
#plt.ylabel("RBias (%)")   
##plt.xlabel("Intensity (mm/hr)")  
##plt.yticks(np.arange(-1, 3, step=0.5)) 
#plt.grid()
#plt.tight_layout()  
#
#plt.subplot(3,2,2)   
##plt.subplot(2,1,1) 
##plt.plot(Q9)
##plt.xticks(np.arange(0, 940, step=1))  
##plt.grid()
##plt.tight_layout()  
##
##plt.subplot(2,1,2) 
#TStep=0
#for n in range(len(listRB9)):    
#    plt.plot(np.arange(TStep,TStep+len(listRB9[n])),listRB9[n])
#    TStep=TStep+len(listRB9[n])  
#plt.title('RBias Sensor No. 9')
#plt.legend(loc='lower right') 
#plt.ylabel("RBias (%)")     
##plt.yticks(np.arange(-1, 3, step=0.5)) 
##plt.xticks(np.arange(0, 940, step=1))    
#plt.grid()
#plt.tight_layout()  
#
#plt.subplot(3,2,3)
#TStep=0
#for n in range(len(listRB7)):    
#    plt.plot(np.arange(TStep,TStep+len(listRB7[n])),listRB7[n])
#    TStep=TStep+len(listRB7[n])
#plt.title('RBias Sensor No. 7')
#plt.legend(loc='lower right') 
#plt.ylabel("RBias (%)")    
##plt.yticks(np.arange(-1, 3, step=0.5))    
#plt.grid()
#plt.tight_layout()  
#
#plt.subplot(3,2,4)
#TStep=0
#for n in range(len(listRB5)):    
#    plt.plot(np.arange(TStep,TStep+len(listRB5[n])),listRB5[n])
#    TStep=TStep+len(listRB5[n])
#plt.title('RBias Sensor No. 5')
#plt.legend(loc='lower right') 
#plt.ylabel("RBias (%)")     
##plt.yticks(np.arange(-1, 3, step=0.5)) 
#plt.grid()
#plt.tight_layout()  
#
#plt.subplot(3,2,5)
#TStep=0
#for n in range(len(listRB11)):    
#    plt.plot(np.arange(TStep,TStep+len(listRB11[n])),listRB11[n])
#    TStep=TStep+len(listRB11[n])
#plt.title('RBias Sensor No.5')
#plt.legend(loc='lower right') 
#plt.ylabel("RBias (%)")    
##plt.yticks(np.arange(-1, 3, step=0.5)) 
#plt.grid()
#plt.tight_layout()  
#
#plt.subplot(3,2,6)
#TStep=0
#for n in range(len(listRB3)):    
#    plt.plot(np.arange(TStep,TStep+len(listRB3[n])),listRB3[n])
#    TStep=TStep+len(listRB3[n])
#plt.title('RBias Sensor No.3')
#plt.legend(loc='lower right') 
#plt.ylabel("RBias (%)")   
##plt.yticks(np.arange(-1, 3, step=0.5)) 
#plt.grid()
#plt.tight_layout()  
#      
 

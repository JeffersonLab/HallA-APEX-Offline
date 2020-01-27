#!/usr/bin/python3.6

#######################################################
# A python script to loop tree2ascii over dp_cuts intervals
# Author: John Williamson
# Created: 11 April 2019
#######################################################


import os


#set baseRootFileName = "/scratch/johnw/apex/Rootfiles/apex_online_"

baseRootFileName = "/w/work3/home/johnw/Rootfiles/apex_%d_1_7_2019.root"
baseRootFileName_part = "/w/work3/home/johnw/Rootfiles/apex_%d_1_7_2019_%d.root"
dpRootFileName = "/w/work3/home/johnw/Rootfiles/apex_%d_1_7_2019_dp_%d.root"

# #set RootFileName_1 = "${baseRootFileName}${Run_number}.root"
# set RootFileName_1 = "${baseRootFileName}${Run_number}_1_7_2019.root"


# echo "$RootFileName_1"

# set DataFileName = "text_cuts/SieveData.f51"
# set DataFileName = "text_cuts/Vars_with_cuts_${Run_number}_1_7_2019.dat"

dataName = "text_cuts/Vars_with_cuts_%d_1_7_2019_dp_all_rast.dat"
#dataName_part = "text_cuts/Vars_with_cuts_%d_1_7_2019_dp_%d.dat"



dp_lim_1 = 1
dp_lim_2 = 7
#dp_lim_2 = 1

# run_list = [[4179,3],[4180,1]]

# JW: Each entry in runlist is as followd: run_number, file_splits, dp_range_lowest, dp_range_highest
# where dp_range is if the cutfiles for that run are seperated into ranges in dp, parameters should be 0,0
# for runs that have not dp split cutfiles

run_list = [[4179,3,1,7],[4180,1,0,0],[4181,3,0,0]]
#run_list = [[4180,1,0,0],[4181,3,0,0]]

# run_list = [[4179,0]]


for i in run_list[:]:    
    print(i[0])

for i in run_list[:]:
    for k in range(0,i[1]+1):
#        print(baseRootFileName % (i,j))
        fullName = baseRootFileName % (i[0])
        outDataName = dataName % (i[0])        

        if(k==0):
            dpRootFileName_ = fullName
        if(k>0):
            dpRootFileName_ = baseRootFileName_part % (i[0],k)

            
        if(i[2]==0 & i[3]==0):
                
            cutName = baseRootFileName % i[0]

            print(" tree2ascii -ap -d SieveVarsBend_rast.def \n -c  %s.VertexCut.cut \n -g %s.FullCut.root \n -o %s \n %s" % (cutName,cutName,outDataName,dpRootFileName_) )

            mycmd = " tree2ascii -ap -d SieveVarsBend_rast.def  -c  %s.VertexCut.cut  -g %s.FullCut.root -o %s %s" % (cutName,cutName,outDataName,dpRootFileName_)

            os.system(mycmd)

            
        else:
        
            for j in range(i[2],i[3]+1 ):
                
                cutName = dpRootFileName % (i[0],j)


                print(" tree2ascii -ap -d SieveVarsBend_rast.def \n -c  %s.SieveCut.cut \n -g %s.FullCut.root \n -o %s \n %s" % (cutName,cutName,outDataName,dpRootFileName_) )
        

                mycmd = " tree2ascii -ap -d SieveVarsBend_rast.def  -c  %s.SieveCut.cut  -g %s.FullCut.root -o %s %s" % (cutName,cutName,outDataName,dpRootFileName_)

                os.system(mycmd)
                


###########################################################################
# LHRS central foil : SieveCentralLHRS
###########################################################################

# echo "# $Run_number"
# tree2ascii -p -N $OutN -d SieveVarsBend.def \
#     -c $RootFileName_1.VertexCut.cut \
#     -g $RootFileName_1.FullCut.root \
#     -o $DataFileName \
#     $RootFileName_1


# echo "# $Run_number  "
# tree2ascii -p -N $OutN -d SieveVarsBend.def \
#     -c $RootFileName_1.ColCut.cut \
#     -g $RootFileName_1.FullCut.root \
#     -o $DataFileName \
#     $RootFileName_1 


# set  RootFileName_1_1 = "${baseRootFileName}${Run_number}_1_7_2019_1.root" 

# tree2ascii -ap -N $OutN -d SieveVarsBend.def \
#     -c $RootFileName_1.ColCut.cut \
#     -g $RootFileName_1.FullCut.root \
#     -o $DataFileName \
#     $RootFileName_1_1 

# set  RootFileName_1_2 = "${baseRootFileName}${Run_number}_1_7_2019_2.root"

# tree2ascii -ap -N $OutN -d SieveVarsBend.def \
#     -c $RootFileName_1.ColCut.cut \
#     -g $RootFileName_1.FullCut.root \
#     -o $DataFileName \
#     $RootFileName_1_2 

 
# set  RootFileName_1_3 = "${baseRootFileName}${Run_number}_1_7_2019_3.root"

# tree2ascii -ap -N $OutN -d SieveVarsBend.def \
#     -c $RootFileName_1.ColCut.cut \
#     -g $RootFileName_1.FullCut.root \
#     -o $DataFileName \
#     $RootFileName_1_3 

  
# set Run_number = "4180"

# set  RootFileName_2 = "${baseRootFileName}${Run_number}_1_7_2019.root"

# echo "# $Run_number "
# tree2ascii -ap -N $OutN -d SieveVarsBend.def \
#     -c $RootFileName_2.ColCut.cut \
#     -g $RootFileName_2.FullCut.root \
#     -o $DataFileName \
#     $RootFileName_2 


# set  RootFileName_2_1 = "${baseRootFileName}${Run_number}_1_7_2019_1.root" 

# tree2ascii -ap -N $OutN -d SieveVarsBend.def \
#     -c $RootFileName_2.ColCut.cut \
#     -g $RootFileName_2.FullCut.root \
#     -o $DataFileName \
#     $RootFileName_2_1 



# set Run_number = "4181"

# set  RootFileName_3 = "${baseRootFileName}${Run_number}_1_7_2019.root"

# echo "# $Run_number "

# tree2ascii -ap -N $OutN -d SieveVarsBend.def \
#     -c $RootFileName_3.ColCut.cut \
#     -g $RootFileName_3.FullCut.root \
#     -o $DataFileName \
#     $RootFileName_3


# set  RootFileName_3_1 = "${baseRootFileName}${Run_number}_1_7_2019_1.root"

# echo "# $Run_number "
# tree2ascii -ap -N $OutN -d SieveVarsBend.def \
#     -c $RootFileName_3.ColCut.cut \
#     -g $RootFileName_3.FullCut.root \
#     -o $DataFileName \
#     $RootFileName_3_1 


# set  RootFileName_3_2 = "${baseRootFileName}${Run_number}_1_7_2019_2.root"

# echo "# $Run_number "
# tree2ascii -ap -N $OutN -d SieveVarsBend.def \
#     -c $RootFileName_3.ColCut.cut \
#     -g $RootFileName_3.FullCut.root \
#     -o $DataFileName \
#     $RootFileName_3_2 


# echo "# 1170 : 0% @ $Id"
# tree2ascii -p -O $Id -N $OutN -d SieveVarsBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1170.root.SieveCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1170.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/${DataFileName} \
# 	ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1170.root

# @ Id += $IdOffSet
# echo "# 1172 : +1% @ $Id"
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1172.rot.SieveCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1172.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/${DataFileName} \
# 	ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1172.root

# @ Id += $IdOffSet
# echo "# 1181 : +2% @ $Id"
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1181.root.SieveCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1181.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/${DataFileName} \
# 	ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1181.root

# @ Id += $IdOffSet
# echo "# 1183 : +3% @ $Id"
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1183.root.SieveCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1183.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/${DataFileName} \
# 	ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1183.root

# @ Id += $IdOffSet
# echo "# 1190 : +4% @ $Id"
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1190.root.SieveCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1190.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/${DataFileName} \
# 	ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1190.root  

# @ Id += $IdOffSet
# echo "# 1169-1170 : Elastic Tail for -3%~-1% @ $Id"
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1169.1170.Kine1.root.SieveCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1169.1170.Kine1.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/${DataFileName} \
# 	ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1169.1170.Kine1.root  
	
# @ Id += $IdOffSet
# echo "# 1169.1170.Kine2 : -4% @ $Id"
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1169.1170.Kine2.root.SieveCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1169.1170.Kine2.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/${DataFileName} \
# 	ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1169.1170.Kine2.root  


# echo "July Data ... "

# @ Id += $IdOffSet
# echo "# 1898-1899 : +0% @ $Id"
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1898.root.SieveCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1898.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/${DataFileName} \
# 	ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1898.root
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1898.root.SieveCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1898.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/${DataFileName} \
# 	ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1899.root
	
# @ Id += $IdOffSet
# echo "# 1911-1912 : +2% @ $Id"
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1911.root.SieveCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1911.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/${DataFileName} \
# 	ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1911.root
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1911.root.SieveCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1911.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/${DataFileName} \
# 	ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1912.root
	
# @ Id += $IdOffSet
# echo "# 1919 : +4% @ $Id"
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1919.root.SieveCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1919.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/${DataFileName} \
# 	ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1919.root
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1919.root.SieveCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1919.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/${DataFileName} \
# 	ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1920.root

# @ Id += $IdOffSet
# echo "# 1898-1899 : Elastic Tail for -3%~-1% @ $Id"
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1898.1899.Kine1.root.SieveCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1898.1899.Kine1.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/${DataFileName} \
# 	ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1898.1899.Kine1.root  
	
# @ Id += $IdOffSet
# echo "# 1898-1899.Kine2 : -4% @ $Id"
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1898.1899.Kine2.root.SieveCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1898.1899.Kine2.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/${DataFileName} \
# 	ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1898.1899.Kine2.root  

###########################################################################
# Run Tests : 2 Pass
###########################################################################

# echo "# 1814 : 2Pass @ $Id"
# tree2ascii -p -O $Id -N $OutN -d SieveVarsBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1814.root.VertexCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1814.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/${DataFileName}_2Pass \
# 	ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1814.root
# 	
# echo "# 1815 : 2Pass @ $Id"
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1815.root.VertexCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1815.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/${DataFileName}_2Pass \
# 	ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1815.root
# 
# 
# echo "# 1819 : 2Pass @ $Id"
# tree2ascii -p -O $Id -N $OutN -d SieveVarsBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1819.root.VertexCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1819.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/${DataFileName}_2Pass_sieveout \
# 	ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1819.root
	
###########################################################################
# Run Tests : 1 Pass
###########################################################################

# echo "# 1898 : 0% @ $Id"
# tree2ascii -p -O $Id -N $OutN -d SieveVarsBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1898.root.VertexCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1898.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/${DataFileName}_1Pass \
# 	ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1898.root
# 	
# echo "# 1899 : 0% @ $Id"
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1898.root.VertexCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1898.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/${DataFileName}_1Pass \
# 	ArchivedRootFiles/SieveCentralLHRS/apex_coinc_1899.root


###########################################################################
# Run Tests : 1 Pass, Iter 2
###########################################################################

# echo "# 1898 : 0% @ $Id"
# tree2ascii -p -O $Id -N $OutN -d SieveVarsRasterBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1898.root.VertexCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1898.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/${DataFileName}_1Pass \
# 	ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1898.root
# 
# echo "# 1899 : 0% @ $Id"
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsRasterBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1898.root.VertexCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1898.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/${DataFileName}_1Pass \
# 	ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1899.root
	
# # Rastered
# 	
# echo "# 1897 : 0% @ $Id"
# tree2ascii -p -O $Id -N $OutN -d SieveVarsRasterBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1897.root.VertexCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1897.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/${DataFileName}_1PassRaster \
# 	ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1897.root
# 
# echo "# 1899 : 0% @ $Id"
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsRasterBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1897.root.VertexCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1897.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/${DataFileName}_1PassRaster \
# 	ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1896.root


# echo "# 1898 from Eric @ $Id"
# tree2ascii -p -O $Id -N $OutN -d SieveVarsBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1898.root.VertexCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1898.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/${DataFileName}_1PassEric \
# 	/w/halla-1/apex/ejensen/rootfiles/apex_1898.0.root
# 	
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1898.root.VertexCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1898.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/${DataFileName}_1PassEric \
# 	/w/halla-1/apex/ejensen/rootfiles/apex_1898.1.root
# 
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBend.def \
# 	-c ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1898.root.VertexCut.cut \
# 	-g ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/apex_coinc_1898.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralLHRS/db_L.vdc.dat.21.Ph.PlusJulyData.replay/${DataFileName}_1PassEric \
# 	/w/halla-1/apex/ejensen/rootfiles/apex_1898.2.root





############## RHRS ##############

# echo "# 1898 : 0% @ $Id"
# tree2ascii -p -O $Id -N $OutN -d SieveVarsRasterBend.RHRS.def \
# 	-c ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1898.root.VertexCut.cut \
# 	-g ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1898.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/${DataFileName}_1Pass \
# 	ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1898.root
# 
# echo "# 1899 : 0% @ $Id"
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsRasterBend.RHRS.def \
# 	-c ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1898.root.VertexCut.cut \
# 	-g ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1898.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/${DataFileName}_1Pass \
# 	ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1899.root

# # unratered, unrastered beam
# 
# echo "# 1898 : 0% @ $Id"
# tree2ascii -p -O $Id -N $OutN -d SieveVarsBended.RHRS.def \
# 	-c ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1898.root.VertexCut.cut \
# 	-g ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1898.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/${DataFileName}_1Pass_urb \
# 	ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1898.root
# 
# echo "# 1899 : 0% @ $Id"
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBended.RHRS.def \
# 	-c ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1898.root.VertexCut.cut \
# 	-g ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1898.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/${DataFileName}_1Pass_urb \
# 	ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1899.root
	
	
# echo "# 1898 : 0% @ $Id"
# tree2ascii -p -O $Id -N $OutN -d SieveVarsBended.RHRS.def \
# 	-c ArchivedRootFiles/SieveCentralRHRS.Iter2/apex_coinc_1898.root.VertexCut.cut \
# 	-g ArchivedRootFiles/SieveCentralRHRS.Iter2/apex_coinc_1898.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralRHRS.Iter2/${DataFileName}_1Pass_urb \
# 	ArchivedRootFiles/SieveCentralRHRS.Iter2/apex_coinc_1898.root
# 
# echo "# 1899 : 0% @ $Id"
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBended.RHRS.def \
# 	-c ArchivedRootFiles/SieveCentralRHRS.Iter2/apex_coinc_1898.root.VertexCut.cut \
# 	-g ArchivedRootFiles/SieveCentralRHRS.Iter2/apex_coinc_1898.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralRHRS.Iter2/${DataFileName}_1Pass_urb \
# 	ArchivedRootFiles/SieveCentralRHRS.Iter2/apex_coinc_1899.root
	
# Rastered
	
# echo "# 1897 : 0% @ $Id"
# tree2ascii -p -O $Id -N $OutN -d SieveVarsRasterBend.RHRS.def \
# 	-c ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1897.root.VertexCut.cut \
# 	-g ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1897.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/${DataFileName}_1PassRaster \
# 	ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1897.root
# 
# echo "# 1899 : 0% @ $Id"
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsRasterBend.RHRS.def \
# 	-c ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1897.root.VertexCut.cut \
# 	-g ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1897.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/${DataFileName}_1PassRaster \
# 	ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1896.root


###########################################################################
# RHRS 2pass sieve in 
###########################################################################


# echo "# 1812 : Sieve in @ $Id"
# tree2ascii -p -O $Id -N $OutN -d SieveVarsBended.RHRS.def \
# 	-c ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1812.root.VertexCut.cut \
# 	-g ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1812.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/${DataFileName}_2Pass_sievein \
# 	ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1812.root
# 
# echo "# 1813 : Sieve in @ $Id"
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBended.RHRS.def \
# 	-c ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1813.root.VertexCut.cut \
# 	-g ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1813.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/${DataFileName}_2Pass_sievein \
# 	ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1813.root


# echo "# 1816 : Sieve half-in @ $Id"
# tree2ascii -p -O $Id -N $OutN -d SieveVarsBended.RHRS.def \
# 	-c ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1816.root.VertexCut.cut \
# 	-g ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1816.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/${DataFileName}_2Pass_sievehalfin \
# 	ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1816.root
	
# echo "# 1818 : Sieve out @ $Id"
# tree2ascii -p -O $Id -N $OutN -d SieveVarsBended.RHRS.def \
# 	-c ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1818.root.VertexCut.cut \
# 	-g ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1818.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/${DataFileName}_2Pass_sieveout \
# 	ArchivedRootFiles/SieveCentralRHRS.Iter2/db_R.vdc.dat.8.Ph.SortTerms.replay/apex_coinc_1818.root

###########################################################################
# RHRS central foil : SieveCentralRHRS
###########################################################################

# echo "# 1898-9 : 0% @ $Id"
# tree2ascii -p -O $Id -N $OutN -d SieveVarsBended.RHRS.def \
# 	-c ArchivedRootFiles/SieveCentralRHRS/apex_coinc_1898.root.SieveCut.cut \
# 	-g ArchivedRootFiles/SieveCentralRHRS/apex_coinc_1898.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralRHRS/${DataFileName} \
# 	ArchivedRootFiles/SieveCentralRHRS/apex_coinc_1898.root
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBended.RHRS.def \
# 	-c ArchivedRootFiles/SieveCentralRHRS/apex_coinc_1898.root.SieveCut.cut \
# 	-g ArchivedRootFiles/SieveCentralRHRS/apex_coinc_1898.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralRHRS/${DataFileName} \
# 	ArchivedRootFiles/SieveCentralRHRS/apex_coinc_1899.root
# 
# @ Id += $IdOffSet
# echo "# 1911-2 : +2% @ $Id"
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBended.RHRS.def \
# 	-c ArchivedRootFiles/SieveCentralRHRS/apex_coinc_1911.root.SieveCut.cut \
# 	-g ArchivedRootFiles/SieveCentralRHRS/apex_coinc_1911.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralRHRS/${DataFileName} \
# 	ArchivedRootFiles/SieveCentralRHRS/apex_coinc_1911.root
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBended.RHRS.def \
# 	-c ArchivedRootFiles/SieveCentralRHRS/apex_coinc_1911.root.SieveCut.cut \
# 	-g ArchivedRootFiles/SieveCentralRHRS/apex_coinc_1911.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralRHRS/${DataFileName} \
# 	ArchivedRootFiles/SieveCentralRHRS/apex_coinc_1912.root
# 
# 
# @ Id += $IdOffSet
# echo "# 1919-20 : +4% @ $Id"
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBended.RHRS.def \
# 	-c ArchivedRootFiles/SieveCentralRHRS/apex_coinc_1920.root.SieveCut.cut \
# 	-g ArchivedRootFiles/SieveCentralRHRS/apex_coinc_1920.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralRHRS/${DataFileName} \
# 	ArchivedRootFiles/SieveCentralRHRS/apex_coinc_1920.root
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBended.RHRS.def \
# 	-c ArchivedRootFiles/SieveCentralRHRS/apex_coinc_1920.root.SieveCut.cut \
# 	-g ArchivedRootFiles/SieveCentralRHRS/apex_coinc_1920.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralRHRS/${DataFileName} \
# 	ArchivedRootFiles/SieveCentralRHRS/apex_coinc_1919.root
# 
# @ Id += $IdOffSet
# echo "# 1898-1899 : -2.5% ~ -0.5% @ $Id"
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBended.RHRS.def \
# 	-c ArchivedRootFiles/SieveCentralRHRS/apex_coinc_1898-9.kine1.root.SieveCut.cut \
# 	-g ArchivedRootFiles/SieveCentralRHRS/apex_coinc_1898-9.kine1.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralRHRS/${DataFileName} \
# 	ArchivedRootFiles/SieveCentralRHRS/apex_coinc_1898-9.kine1.root
# 
# @ Id += $IdOffSet
# echo "# 1898-1899 : -4% ~ -2.5% @ $Id"
# tree2ascii -ap -O $Id -N $OutN -d SieveVarsBended.RHRS.def \
# 	-c ArchivedRootFiles/SieveCentralRHRS/apex_coinc_1898-9.kine2.root.SieveCut.cut \
# 	-g ArchivedRootFiles/SieveCentralRHRS/apex_coinc_1898-9.kine2.root.FullCut.root \
# 	-o ArchivedRootFiles/SieveCentralRHRS/${DataFileName} \
# 	ArchivedRootFiles/SieveCentralRHRS/apex_coinc_1898-9.kine2.root
# 
# 
#echo "Done!"

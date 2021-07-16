#  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#  :::                                                                                           :::
#  :::  Set your input parameters here for MC-Glauber simulation.                                :::
#  :::  The code will not read any line of this file starting with symbol "#" or any empty line. :::
#  :::                                                                                           :::
#  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



#collision energy
#[Info] in GeV
SNN 2760.0


#collision species & impact parameter
#[Info] species can be (Au,U,Pb,p) 
projectile Pb
target Pb
bmin  0.0
bmax  17.0

npp 2.25
xhard 0.14

#output file name
root_output_file_name    Pb_Pb_2760_min_bias_xhard_0.14.root

# :: END :: #

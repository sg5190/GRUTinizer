# GRUTinizer
Install GRUTinizer
- Needs a root installation to work, at present I have tested up tp 6.22.08 on FRIB computers. 
  Note between version 6.18 and 6.22 the library lPyROOT seems to have been replaced my lROOTTPython
  on the FRIB machines, the makefile has been adjusted to use either library based on the root version
- For root you will need a variable $ROOTSYS defined in your bashrc pointing to your root installation
 
- Add these lines to your bashrc
export GRUTSYS=/path/to/GRUTinizer
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GRUTSYS/lib
source $GRUTSYS/thisgrut.sh

- Go the the main GRUTinizer directory and type
  make -j N
where N is the number of cores to use

Assuming no errors (there are still some unused class rule warnings I am working on), you have successfully installed GRUTinizer Woo!

Unpacking a file

- Raw files have the file extension .evt (NSCLDAQ) or .dat file (GRETINA DAQ), both files are unpacked into trees with the command
  grutinizer -lmtq run.evt/Global.dat channels.cal(Optional)
- channels.cal contains a file map and calibrations, it is not too important for GRETINA data at the moment and can be ignored.


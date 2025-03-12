#! /bin/bash

executable=sector_corner_modes
make $executable

# Stuff to move
important_files="$executable
                 $executable.cc
           	     $executable.bash
                 Modes_sector/"

# Setup directories YOU MUST PICK A NAME FOR YOUR OURPUT DIRECTORY.
main_dir=NEW_$executable   # Change 0 to be the name of your output directory. e.g. NEW_xyz
if [ main_dir == 0 ]; then
   echo " "
   echo "You need to give your new output a name"
fi

# Check if ther main directory already exists, if so ask whether you want to delete it
if [ -e $main_dir ]; then
    echo " "
    echo "WARNING: Directory " $main_dir " already exists!"
    read -p "         DELETE it and continue? [Y/n] " yn
    case $yn in
        ''|[Yy]* ) rm -rf $main_dir;;
        [Nn]* ) echo "Can't continue until you move $main_dir"; exit;;
    esac
#    echo "DELETING IT"
#    rm -rf $main_dir
    echo " "
fi
mkdir $main_dir

# Transfer the important files to the main directory and go there
cp -r $important_files $main_dir
cd $main_dir

reslt_dir=RESLT
mkdir $reslt_dir

#-------------------------------------------------------------
# Oomph-lib time!
#-------------------------------------------------------------
echo "Doing "$reslt_dir
./$executable \
    $@ \
    > $reslt_dir/OUTPUT

echo " "
echo " "
echo "Done!"
echo " "
echo " "


exit


# Install/unInstall package classes in LAMMPS

if (test $1 = 1) then

    cp -p bond_centforce.cpp ..
    cp -p atom_vec_bondxo.cpp ..
    cp -p atom_vec_angxo.cpp ..

    cp -p bond_centforce.h ..
    cp -p atom_vec_bondxo.h ..
    cp -p atom_vec_angxo.h ..
    
elif (test $1 = 0) then

    rm -f ../bond_centforce.cpp
    rm -f ../atom_vec_bondxo.cpp
    rm -f ../atom_vec_angxo.cpp
    
    rm -f ../bond_centforce.h
    rm -f ../atom_vec_bondxo.h
    rm -f ../atom_vec_angxo.h
    
fi

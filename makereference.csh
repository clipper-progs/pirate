#!/bin/tcsh -f
echo "This script will help you make reference structures to use with"
echo "pirate."
echo
echo "You will need the PDB code of a structure for which structure factors"
echo "have been deposited, to at least the same resolution as your work"
echo "structure. Acurate, complete data is preferred."
echo
echo "Due to variations in the deposited files, this may not work."
rm log

# loop over pdbs
foreach pdbcode ($*)

# do a pdb
set name = $pdbcode:l:l:l:l
echo "PDBCODE: $name" >> log

# download info
wget ftp://ftp.ebi.ac.uk/pub/databases/msd/pdb_uncompressed/pdb{$name}.ent
wget ftp://ftp.ebi.ac.uk/pub/databases/rcsb/pdb/data/structures/all/structure_factors/r{$name}sf.ent.Z

# extract downloaded info
zcat r{$name}sf.ent.Z | sed 's/^_refln.F_meas$/_refln.F_meas_au/' | sed 's/^_refln.F_meas_sigma$/_refln.F_meas_sigma_au/' > tmp.ent
set CELL = `awk '/CRYST1/{print substr($0,8,48)}' pdb{$name}.ent`
set SYMM = `awk '/CRYST1/{print substr($0,56,13)}' pdb{$name}.ent`

# calc solvent content
cat pdb{$name}.ent | grep -v " HOH " | grep -v " WAT " | grep -v " H20 " > tmp0.pdb
ncsmask xyzin tmp0.pdb mskout tmp0.msk << eof
CELL $CELL
SYMM "$SYMM"
GRID 256 256 256
RADI 2.2
eof
mapmask mskin tmp0.msk mskout tmp1.msk << eof | tee -a log
extend overlap
xyzlim cell
pad 0.0
eof

# make structure factor file
cif2mtz hklin tmp.ent hklout tmp0.mtz << eof
CELL $CELL
SYMM "$SYMM"
eof
echo "NOTE: THIS STEP WILL FAIL IF THE DATA HAVE ALREADY BEEN TRUNCATED I->F."
echo "AN ERROR MESSAGE AT THIS STAGE IS PERFECTLY NORMAL FOR MOST PDBS."
truncate hklin tmp0.mtz hklout tmp1.mtz << eof
labin IMEAN=I SIGIMEAN=SIGI
labout F=FP SIGF=SIGFP
eof
if ($status) mv tmp0.mtz tmp1.mtz
mtzdump hklin tmp1.mtz << eof > tmp.log
eof
set RESO = `awk '/ A )/{res=0.999*$6;if(res<1.9)res=1.9;print res}' tmp.log`
unique hklout tmp2.mtz << eof
CELL $CELL
SYMM "$SYMM"
RESO $RESO
LABOUT F=FP SIGF=SIGFP
eof
freerflag hklin tmp2.mtz hklout tmp3.mtz << eof
eof
cad hklin1 tmp1.mtz hklin2 tmp3.mtz hklout tmp4.mtz << eof
LABIN FILE 1 E1=FP E2=SIGFP
LABIN FILE 2 E1=FreeR_flag
eof

awk '/^CRYST/||/^SCALE/||/^ATOM/||/^HETAT/{print}' pdb{$name}.ent > tmp0.pdb
refmac5 XYZIN tmp0.pdb XYZOUT tmp1.pdb HKLIN tmp4.mtz HKLOUT tmp5.mtz << eof | tee -a log
make check NONE
refi -
    type UNRE -
    resi MLKF -
    meth CGMAT
ncyc 0
scal -
    type BULK -
    LSSC -
    ANISO
solvent YES
weight -
    NOEX -
    MATRIX 0.3
scal WORK
LABIN FP=FP SIGFP=SIGFP
PNAME unknown
DNAME unknown
RSIZE 80
END
eof
if ($status) exit
chltofom -mtzin tmp5.mtz -mtzout tmp6.mtz -colin-phifom "/*/*/[PHWT,FOM]" -colout FC
mv tmp6.mtz $name.mtz
rm -f tmp.log tmp.ent tmp0.msk tmp1.msk tmp0.pdb tmp1.pdb tmp0.mtz tmp1.mtz tmp2.mtz tmp3.mtz tmp4.mtz tmp5.mtz

echo
echo "Reference structure complete: $name.mtz"

end

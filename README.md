oeante
======

Antechamber (AMBER GAFF parameterization engine) implementation using OpenEye tools

Original code by Richard Dixon, released as GPLv2.

# Introduction

Antechamber is a set of programs designed to ease the pain associated with
setting up and running novel molecules with the Amber force field and
program.  Prior to programmatic support from something like antechamber,
the general procedure was to a) manually select parameters in an ad hoc
fashion from related molecules or b) spend a fair amount of time running QM
calculations.

The structural parameters derived via option b) were often
in agreement with those obtained from option a) suggesting that the overall
process could be automated.  Incorporation of a number of parameter
estimation schemes used in other force fields, mmff94, uff and mm2 in
particular, created a third type of parameter generation.  If option a) can
be described as "ad hoc" and option b) as "rigorous", option c) is
semi-rigorous.  The parameter generation schemes do take into account the
chemical environment of novel molecules, but ultimately invoke a rule-based
procedure to assign values.

To salve my conscience, please consider option b) on molecules that are
important to you.  If you are going to spend days or weeks of simulation
time, a few hours of quantum chemistry could save you some trouble.

Now that that's out of the way:  of course you are going to run molecules
without even looking at the parameters that come out of the calculation.
That's why we automate things.  Here's how to use these tools

# Program descriptions

## OEante.py

###  Requirements
*  Requires the OEchem toolkit from OpenEye.  As currently configured, also
  requires the OEquacpac and OEomega toolkits to generate charges and
  initial coordinates respectively.  If you can generate these data in some
  other way, you can avoid the dependence on these toolkits.  Also requires
  python 2.6 due to the use of sets and product functions.

###  Usage
*  OEante.py -h will generate a help page
*  OEante.py FILE.[smi,sdf,mol2] is all that is really needed.  This will
  generate a mol2 format file (look.mol2 by default) with the atom type[1] and
  charge fields set.
*  The file gaffsmarts.txt must be in the current directory or explicitly
  indicated by the -t flag.

## OEparmchk.py

###  Requirements
*  Requires the OEchem toolkit from OpenEye.

###  Usage
*  OEparmchk.py -h will generate a help page
*  OEparmchk.py look.mol2 will generate the files:
     look_pchk.mol2, ante_look.leaprc and ante_look.frcmod
*  These files are sufficient to run leap to generate topology and
  coordinate files and run a gas phase minimization using sander.  The
  bare-bones leaprc file generated here can be modified as needed for more
  sophisticated calculations.

###  Considerations
*  The example usage above embodies three important decisions that may not
  work for your cases.  
  *  Decision 1:  all structural parameters are generated using empirical
    rules.  If you wish the program to use an existing set of parameters,
    e.g. the current gaff set, this file must be indicated with the -p
    flag, e.g. -p gaff.dat.

  *  Decision 2:  improper torsions receive a force constant of 0.0 by
    default.  The original concept of improper torsions is that they would
    only be applied to atoms which actually pucker.  So a true "amber"
    molecule might not have an improper torsion applied to every sp2
    center.  My own preference is to apply a force to all planar centers
    and I could easily be persuaded to make this the default.  To invoke
    this, use the -f flag, e.g. -f 10.5.

  *  Decision 3:  angle bending force constants between heavy atoms are fixed
    at 70 and those involving hydrogen at 35.  The uff scheme provides a
    way to tune these parameters a bit more and can be invoked with the -u
    flag.

## gaffsmarts.txt

  These are the rules to determine atom types.  The current incarnation
  tries to generate exclusive rules that only match the correct chemical
  environment, i.e. an atom matches one and only one rule.  This has lead
  to complex rules and will most likely be the main source of concern when
  using these programs.  I can only see tweaking the rules a bit going
  forward given their current complexity.  If wholesale changes are needed,
  I anticipate changing to a "last rule wins" system so an atom could match
  multiple rules.  Such ambiguities would be resolved by the order in which
  the rules are listed in the input file.

#Usage examples

  The sequence of commands:
```
    OEante.py nma.sdf
    OEparmchk.py -f 10.5 -u look.mol2
    $AMBERHOME/exe/tleap -f ante_look.leaprc
    $AMBERHOME/exe/sander -O -i min.in -o min.out -p lig.top -c lig.crd
    $AMBERHOME/exe/ambpdb -aatm -p lig.top < restrt > min_out.pdb
```
  generates a minimized structure using estimated parameters and am1-bcc
  charges.  The final energy is -27.41 kcal.  Replacing the OEparmchk.py
  and tleap lines with:
```
    OEParmchk.py -f 10.5 -u -p gaff.dat look.mol2
    $AMBERHOME/exe/tleap -f gaff_look.leaprc
```
  generates a minimized structure using the gaff parameter set and am1-bcc
  charges (from oequacpac).  The final energy in this case is -38.83.

## Notes
[1]  OpenEye file readers/writers try hard to impose "proper" formatting
     for the file type in question.  Reading of Mol2 format files assumes
     that Sybyl atom types are used and assigns atom properties according
     to this assumption.  In order to use GAFF types and defeat this
     process, I prepend an "_" to each type.  This forces oechem to retype.
     The "_" is removed by OEparmchky.py.
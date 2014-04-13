#!/usr/bin/env python

""" OEante:  AMBER force field and program atom-typing and
    parameter generation.  Richard Dixon, Vertex Pharma

    Copyright (2010) Vertex Pharma.
    This software is distributed under the terms of the
    GNU General Public License
    """

import sys, os, string
import re
from openeye.oechem import *
from math import *
from optparse import OptionParser

def uniqueVdw(mol, List, debug):
    """ Assign missing nonbonded parameters based on atomic number"""

    _Rstar = {1 : 1.443, 6 : 1.908, 7 : 1.824, 8 : 1.75, 9 : 1.75,
              15 : 2.1, 16 : 2.0, 17 : 1.948, 35 : 2.22, 53 : 2.35}
    _Eps = {1 : 0.157, 6 : 0.086, 7 : 0.17, 8 : 0.21, 9 : 0.061,
            15 : 0.2, 16 : 0.25, 17 : 0.265, 35 : 0.32, 53 : 0.4}
    
    VdW = []

    for atom in mol.GetAtoms():
        if atom.GetType().replace("_","") not in List:
            List.append(atom.GetType().replace("_",""))
            try:
                VdW.append("%s %f %f" % \
                               (atom.GetType().replace("_",""), _Rstar[atom.GetAtomicNum()], \
                                    _Eps[atom.GetAtomicNum()]))
            except:
                continue
    if debug == True:
        print "NONBON"
        for i in xrange(len(VdW)):
            tok = string.split(VdW[i])
            print "%4s \t %7.4f \t %7.4f" % (tok[0], float(tok[1]), float(tok[2]))

        print "\n\n"
    return(VdW)
        
def uniqueImpropers(mol, Force, List, debug):
    """ find 3-coordinate sp2 centers and apply improper parameters """
    """ All impropers are generic, so add a specific if no match is found """
    Impropers = []
    I = []
    J = []
    K = []
    L = []
    
    for atom in mol.GetAtoms():
        if atom.GetDegree() == 3 and (atom.GetHyb() == 2 or atom.GetType().replace("_","") == "n"):
            Itemp = -1
            Jtemp = -1
            Ktemp = atom.GetIdx()
            Ltemp = -1
            for bond in atom.GetBonds():
                if bond.GetBgnIdx() == atom.GetIdx():
                    if bond.GetEnd().GetDegree() == 1 and \
                       Ltemp == -1:
                        Ltemp = bond.GetEndIdx()
                    elif bond.GetEnd().GetDegree() == 1 and \
                         Jtemp == -1:
                        Jtemp = bond.GetEndIdx()
                    elif Itemp == -1:
                        Itemp = bond.GetEndIdx()
                    elif Jtemp == -1:
                        Jtemp = bond.GetEndIdx()
                    elif Ltemp == -1:
                        Ltemp = bond.GetEndIdx()
                    else:
                        continue
                elif bond.GetEndIdx() == atom.GetIdx():
                    if bond.GetBgn().GetDegree() == 1 and \
                       Ltemp == -1:
                        Ltemp = bond.GetBgnIdx()
                    elif bond.GetBgn().GetDegree() == 1 and \
                         Jtemp == -1:
                        Jtemp = bond.GetBgnIdx()
                    elif Itemp == -1:
                        Itemp = bond.GetBgnIdx()
                    elif Jtemp == -1:
                        Jtemp = bond.GetBgnIdx()
                    elif Ltemp == -1:
                        Ltemp = bond.GetBgnIdx()
                    else:
                        continue
            I.append(Itemp)
            J.append(Jtemp)
            K.append(Ktemp)
            L.append(Ltemp)

    for ii in xrange(len(I)):
        i = mol.GetAtom(OEHasAtomIdx(I[ii])).GetType().replace("_","")
        j = mol.GetAtom(OEHasAtomIdx(J[ii])).GetType().replace("_","")
        k = mol.GetAtom(OEHasAtomIdx(K[ii])).GetType().replace("_","")
        l = mol.GetAtom(OEHasAtomIdx(L[ii])).GetType().replace("_","")
        Type = "%-2s-%-2s-%-2s-%-2s" % (i, j, k, l)
        TypeGen1 = "X -%-2s-%-2s-%-2s" % (j, k, l)
        TypeGen2 = "X -X -%-2s-%-2s" % (k, l)
        if Type not in List and TypeGen1 not in List and TypeGen2 not in List:
            Phase = 180.0
            xVal = 2.0
            Vx2 = Force
            List.append(Type)
            Impropers.append("%s %f %f %f" % (Type, Vx2, Phase, xVal))

    if debug == True:
        print "IMPROPER"
        for i in xrange(len(Impropers)):
            tok = string.split(Impropers[i])
            if len(tok) == 7:
                print "%-2s%-3s%-3s%-3s \t %6.2f \t %6.2f \t %6.2f" % \
                      (tok[0], tok[1], tok[2], tok[3], \
                       float(tok[4]), float(tok[5]), float(tok[6]))
            elif len(tok) == 6:
                if len(tok[0]) > 3:
                    print "%-5s%-3s%-3s \t %6.2f \t %6.2f \t %6.2f" % \
                          (tok[0], tok[1], tok[2], \
                           float(tok[3]), float(tok[4]), float(tok[5]))
                elif len(tok[1]) > 3:
                    print "%-2s%-6s%-3s \t %6.2f \t %6.2f \t %6.2f" % \
                          (tok[0], tok[1], tok[2], \
                           float(tok[3]), float(tok[4]), float(tok[5]))
                else:
                    print "%-2s%-3s%-6s \t %6.2f \t %6.2f \t %6.2f" % \
                          (tok[0], tok[1], tok[2], \
                           float(tok[3]), float(tok[4]), float(tok[5]))
            elif len(tok) == 5:
                if len(tok[0]) > 3:
                    print "%-8s%-3s \t %6.2f \t %6.2f \t %6.2f" % \
                          (tok[0], tok[1], \
                           float(tok[2]), float(tok[3]), float(tok[4]))
                else:
                    print "%-2s%-9s \t %6.2f \t %6.2f \t %6.2f" % \
                          (tok[0], tok[1], \
                           float(tok[2]), float(tok[3]), float(tok[4]))
            else:
                print "%-11s \t %6.2f \t %6.2f \t %6.2f" % \
                      (tok[0], \
                       float(tok[1]), float(tok[2]), float(tok[3]))
        print ""
    
    return Impropers
    
    
def uniqueTorsions(mol, Angles, Iang, Jang, Kang, List, debug):
    """ Find standard torsion and apply types """
    """ Adds generic type if none found in List """
    Torsions = []
    I = []
    J = []
    K = []
    L = []

    for i in xrange(len(Iang)-1):
        for j in xrange(i+1,len(Iang)):
            if Jang[i] == Iang[j] and Iang[i] == Jang[j]:
                I.append(Kang[i])
                J.append(Jang[i])
                K.append(Jang[j])
                L.append(Kang[j])
            elif Jang[i] == Kang[j] and Iang[i] == Jang[j]:
                I.append(Kang[i])
                J.append(Jang[i])
                K.append(Jang[j])
                L.append(Iang[j])
            elif Jang[i] == Iang[j] and Kang[i] == Jang[j]:
                I.append(Iang[i])
                J.append(Jang[i])
                K.append(Jang[j])
                L.append(Kang[j])
            elif Jang[i] == Kang[j] and Kang[i] == Jang[j]:
                I.append(Iang[i])
                J.append(Jang[i])
                K.append(Jang[j])
                L.append(Iang[j])
            elif Jang[j] == Iang[i] and Iang[j] == Jang[i]:
                I.append(Kang[j])
                J.append(Jang[j])
                K.append(Jang[i])
                L.append(Kang[i])
            elif Jang[j] == Kang[i] and Iang[j] == Jang[i]:
                I.append(Kang[j])
                J.append(Jang[j])
                K.append(Jang[i])
                L.append(Iang[i])
            elif Jang[j] == Iang[i] and Kang[j] == Jang[i]:
                I.append(Iang[j])
                J.append(Jang[j])
                K.append(Jang[i])
                L.append(Kang[i])
            elif Jang[j] == Kang[i] and Kang[j] == Jang[i]:
                I.append(Iang[j])
                J.append(Jang[j])
                K.append(Jang[i])
                L.append(Iang[i])
    
    for ii in xrange(len(I)):
        i = mol.GetAtom(OEHasAtomIdx(I[ii])).GetType().replace("_","")
        j = mol.GetAtom(OEHasAtomIdx(J[ii])).GetType().replace("_","")
        k = mol.GetAtom(OEHasAtomIdx(K[ii])).GetType().replace("_","")
        l = mol.GetAtom(OEHasAtomIdx(L[ii])).GetType().replace("_","")
        if j > k:
            TypeGen = "X -%-2s-%-2s-X " % (k,j)
        else:
            TypeGen = "X -%-2s-%-2s-X " % (j,k)
        if i > l:
            Type = "%-2s-%-2s-%-2s-%-2s" % (l,k,j,i)
        else:
            Type = "%-2s-%-2s-%-2s-%-2s" % (i,j,k,l)
        if Type not in List and TypeGen not in List:
            List.append(TypeGen)
            numBonds, Vx2, Phase, xVal = doTorsion(mol, I[ii], J[ii], K[ii],
                                                   L[ii], debug)
            Torsions.append("%s %d %f %f %f" % (TypeGen, numBonds, Vx2, Phase,
                                                xVal))
    if debug == True:
        print "DIHE"
        for tor in Torsions:
            tok = string.split(tor)
            if len(tok) == 8:
                print "%-2s%-3s%-3s%-3s \t %d \t %6.3f \t %6.3f \t %6.3f" % \
                      (tok[0], tok[1], tok[2], tok[3], int(tok[4]), \
                       float(tok[5]), float(tok[6]), float(tok[7]))
            elif len(tok) == 7:
                if len(tok[0]) > 3:
                    print "%-5s%-3s%-3s \t %d \t %6.3f \t %6.3f \t %6.3f" % \
                          (tok[0], tok[1], tok[2], int(tok[3]), \
                           float(tok[4]), float(tok[5]), float(tok[6]))
                elif len(tok[1]) > 3:
                    print "%-2s%-6s%-3s \t %d \t %6.3f \t %6.3f \t %6.3f" % \
                          (tok[0], tok[1], tok[2], int(tok[3]), \
                           float(tok[4]), float(tok[5]), float(tok[6]))
                else:
                    print "%-2s%-3s%-6s \t %d \t %6.3f \t %6.3f \t %6.3f" % \
                          (tok[0], tok[1], tok[2], int(tok[3]), \
                           float(tok[4]), float(tok[5]), float(tok[6]))
            elif len(tok) == 6:
                if len(tok[0]) > 3:
                    print "%-8s%-3s \t %d \t %6.3f \t %6.3f \t %6.3f" % \
                          (tok[0], tok[1], int(tok[2]), \
                           float(tok[3]), float(tok[4]), float(tok[5]))
                else:
                    print "%-2s%-9s \t %d \t %6.3f \t %8.3f \t %6.3f" % \
                          (tok[0], tok[1], int(tok[2]), \
                           float(tok[3]), float(tok[4]), float(tok[5]))
            else:
                print "%-11s \t %d \t %6.3f \t %8.3f \t %6.3f" % \
                      (tok[0], int(tok[1]), \
                       float(tok[2]), float(tok[3]), float(tok[4]))
        print""
        
    return Torsions

def doTorsion(mol, I, J, K, L, debug):
    """ Estimate force field parameters for input torsion
    This method follows Halgren MMFF V."""

    Ummff = {6 : 2.0, 7 : 2.0, 8 : 2.0,
             14: 1.25, 15 : 1.25, 16: 1.25}
    Vmmff = {6 : 2.12, 7 : 1.5, 8 : 0.2,
             14: 1.22, 15 : 2.40, 16: 0.49}
    
    numBonds = 0
    Vx2 = 0.0
    Phase = 0.0
    xVal = 0.0
    vScale = 1.0
    
    iAtom = mol.GetAtom(OEHasAtomIdx(I))
    jAtom = mol.GetAtom(OEHasAtomIdx(J))
    kAtom = mol.GetAtom(OEHasAtomIdx(K))
    lAtom = mol.GetAtom(OEHasAtomIdx(L))

    numBonds = (jAtom.GetDegree() - 1) * (kAtom.GetDegree() - 1)

    if jAtom.GetHyb() == 1 or kAtom.GetHyb() == 1:
        Vx2 = 0.0
        Phase = 0.0
        xVal = 2.0
    elif jAtom.GetHyb() == 3 and kAtom.GetHyb == 3:
        xVal = 3.0
        Phase = 0.0

        """ Turn off force for certain J,K pairs """
        if jAtom.GetDegree() == 4:
            if kAtom.GetDegree() == 3:
                if kAtom.GetValence() == 4:
                    vScale = 0.0
                elif kAtom.GetAtomicNum() == 6 or kAtom.GetAtomicNum() == 7 \
                     or kAtom.GetAtomicNum() == 8 or kAtom.GetAtomicNum() \
                     == 14 or kAtom.GetAtomicNum() == 15 or \
                     kAtom.GetAtomicNum() == 16:
                    vScale = 0.0
            elif kAtom.GetDegree() == 2:
                if kAtom.GetValence() == 3:
                    vScale = 0.0
                elif kAtom.GetAtomicNum() == 6 or kAtom.GetAtomicNum() == 7 \
                     or kAtom.GetAtomicNum() == 8 or kAtom.GetAtomicNum() \
                     == 14 or kAtom.GetAtomicNum() == 15 or \
                     kAtom.GetAtomicNum() == 16:
                    vScale = 0.0
        elif kAtom.GetDegree() == 4:
            if jAtom.GetDegree() == 3:
                if jAtom.GetValence() == 4:
                    vScale = 0.0
                elif jAtom.GetAtomicNum() == 6 or jAtom.GetAtomicNum() == 7 \
                     or jAtom.GetAtomicNum() == 8 or jAtom.GetAtomicNum() \
                     == 14 or jAtom.GetAtomicNum() == 15 or \
                     jAtom.GetAtomicNum() == 16:
                    vScale = 0.0
            elif jAtom.GetDegree() == 2:
                if jAtom.GetValence() == 3:
                    vScale = 0.0
                elif jAtom.GetAtomicNum() == 6 or jAtom.GetAtomicNum() == 7 \
                     or jAtom.GetAtomicNum() == 8 or jAtom.GetAtomicNum() \
                     == 14 or jAtom.GetAtomicNum() == 15 or \
                     jAtom.GetAtomicNum() == 16:
                    vScale = 0.0
                
        
        if jAtom.GetAtomicNum() == 8:
            if kAtom.GetAtomicNum() == 8:
                Vx2 = -2.0
            elif kAtom.GetAtomicNum() == 16:
                Vx2 = -4.0
        elif jAtom.GetAtomicNum() == 16:
            if kAtom.GetAtomicNum() == 8:
                Vx2 = -4.0
            elif kAtom.GetAtomicNum() == 16:
                Vx2 = -8.0
        else:
            try:
                Vj = Vmmff[jAtom.GetAtomicNum()]
                Vk = Vmmff[kAtom.GetAtomicNum()]
                Vx2 = Vscale * sqrt(Vj * Vk)/ float(numBonds)
            except:
                Vx2 = 0.0
    else:
        Beta = 0.0
        Pie = 0.0
        xVal = 2.0
        Phase = 180.0
        for bond in jAtom.GetBonds():
            if bond.GetEndIdx() == kAtom.GetIdx() or \
               bond.GetBgnIdx() == kAtom.GetIdx():
                theBond = bond
                
        if jAtom.IsAromatic() == True and kAtom.IsAromatic() == True:
            if theBond.IsAromatic() == True:
                if jAtom.GetAtomicNum() != 7 and \
                   jAtom.GetAtomicNum() != 8 and \
                   jAtom.GetAtomicNum() != 16 and \
                   kAtom.GetAtomicNum() != 7 and \
                   kAtom.GetAtomicNum() != 8 and \
                   kAtom.GetAtomicNum() != 16:
                    Pie = 0.5
                else:
                    Pie = 0.3

                if jAtom.GetDegree() == 4 and kAtom.GetDegree() == 3:
                    Beta = 3.0
                elif jAtom.GetDegree() == 3 and kAtom.GetDegree() == 4:
                    Beta = 3.0
                else:
                    Beta = 6.0
            else:
                Pie = 0.15
                Beta = 6.0
        elif theBond.GetOrder() == 2:
            Pie = 1.0
            Beta = 6.0
        elif theBond.GetOrder() == 1:
            if (jAtom.GetAtomicNum() == 7 or \
                jAtom.GetAtomicNum() == 8 or \
                jAtom.GetAtomicNum() == 16) and \
                (kAtom.GetAtomicNum() == 7 or \
                 kAtom.GetAtomicNum() == 8 or \
                 kAtom.GetAtomicNum() == 16):
                Pie = 0.0
                Beta = 0.0
            elif (jAtom.GetAtomicNum() == 7 or \
                  jAtom.GetAtomicNum() == 8 or \
                  jAtom.GetAtomicNum() == 16) and \
                  (kAtom.GetAtomicNum() == 6 or \
                   kAtom.GetAtomicNum() == 14 or \
                   kAtom.GetAtomicNum() == 15):
                Beta = 6.0
                Pie = 0.5
            elif (kAtom.GetAtomicNum() == 7 or \
                  kAtom.GetAtomicNum() == 8 or \
                  kAtom.GetAtomicNum() == 16) and \
                  (jAtom.GetAtomicNum() == 6 or \
                   jAtom.GetAtomicNum() == 14 or \
                   jAtom.GetAtomicNum() == 15):
                Beta = 6.0
                Pie = 0.5
            else:
                Beta = 6.0
                if kAtom.GetAtomicNum() <= 10:
                    Pie = 0.3
                else:
                    Pie = 0.15
        try:
            Vx2 = Beta * Pie * sqrt(Ummff[jAtom.GetAtomicNum()] * \
                                    Ummff[kAtom.GetAtomicNum()])
        except:
            Vx2 = 0.0
        
    return numBonds, Vx2, Phase, xVal   

def uniqueAngles(mol, Bonds, uff, List, debug):
    I = []
    J = []
    K = []
    aList = []
    aAngles = []
    Angles = []

    for bond1 in mol.GetBonds():
        for bond2 in mol.GetBonds():
            if bond1.GetBgnIdx() == bond2.GetBgnIdx() and \
                   bond1.GetEndIdx() != bond2.GetEndIdx():

                i = bond2.GetEnd().GetType().replace("_","")
                j = bond1.GetBgn().GetType().replace("_","")
                k = bond1.GetEnd().GetType().replace("_","")
                if i > k:
                    Type = "%-2s-%-2s-%-2s" % (k,j,i)
                else:
                    Type = "%-2s-%-2s-%-2s" % (i,j,k)
                if Type not in List or Type not in aList:
                    T0, Keq = doAngle(bond1, bond2, Bonds, uff, debug)
                    if Type not in List:
                        List.append(Type)
                        Angles.append("%s %f %f" % (Type, Keq, T0))
                    if Type not in aList:
                        aList.append(Type)
                        aAngles.append("%s %f %f" % (Type, Keq, T0))

                if bond1.GetEndIdx() > bond2.GetEndIdx():
                    I.append(bond2.GetEndIdx())
                    J.append(bond1.GetBgnIdx())
                    K.append(bond1.GetEndIdx())
                else:
                    I.append(bond1.GetEndIdx())
                    J.append(bond1.GetBgnIdx())
                    K.append(bond2.GetEndIdx())
    
            elif bond1.GetBgnIdx() == bond2.GetEndIdx() and \
                     bond1.GetEndIdx() != bond2.GetBgnIdx():

                i = bond2.GetBgn().GetType().replace("_","")
                j = bond1.GetBgn().GetType().replace("_","")
                k = bond1.GetEnd().GetType().replace("_","")
                if i > k:
                    Type = "%-2s-%-2s-%-2s" % (k,j,i)
                else:
                    Type = "%-2s-%-2s-%-2s" % (i,j,k)
                if Type not in List or Type not in aList:
                    T0, Keq = doAngle(bond1, bond2, Bonds, uff, debug)
                    if Type not in List:
                        List.append(Type)
                        Angles.append("%s %f %f" % (Type, Keq, T0))
                    if Type not in aList:
                        aList.append(Type)
                        aAngles.append("%s %f %f" % (Type, Keq, T0))

                if bond1.GetEndIdx() > bond2.GetBgnIdx():
                    I.append(bond2.GetBgnIdx())
                    J.append(bond1.GetBgnIdx())
                    K.append(bond1.GetEndIdx())
                else:
                    I.append(bond1.GetEndIdx())
                    J.append(bond1.GetBgnIdx())
                    K.append(bond2.GetBgnIdx())
    
            elif bond1.GetEndIdx() == bond2.GetBgnIdx() and \
                     bond1.GetBgnIdx() != bond2.GetEndIdx():

                i = bond2.GetEnd().GetType().replace("_","")
                j = bond1.GetEnd().GetType().replace("_","")
                k = bond1.GetBgn().GetType().replace("_","")
                if i > k:
                    Type = "%-2s-%-2s-%-2s" % (k,j,i)
                else:
                    Type = "%-2s-%-2s-%-2s" % (i,j,k)
                if Type not in List or Type not in aList:
                    T0, Keq = doAngle(bond1, bond2, Bonds, uff, debug)
                    if Type not in List:
                        List.append(Type)
                        Angles.append("%s %f %f" % (Type, Keq, T0))
                    if Type not in aList:
                        aList.append(Type)
                        aAngles.append("%s %f %f" % (Type, Keq, T0))

                if bond1.GetBgnIdx() > bond2.GetEndIdx():
                    I.append(bond2.GetEndIdx())
                    J.append(bond1.GetEndIdx())
                    K.append(bond1.GetBgnIdx())
                else:
                    I.append(bond1.GetBgnIdx())
                    J.append(bond1.GetEndIdx())
                    K.append(bond2.GetEndIdx())
    
            elif bond1.GetEndIdx() == bond2.GetEndIdx() and \
                     bond1.GetBgnIdx() != bond2.GetBgnIdx():

                i = bond2.GetBgn().GetType().replace("_","")
                j = bond1.GetEnd().GetType().replace("_","")
                k = bond1.GetBgn().GetType().replace("_","")
                if i > k:
                    Type = "%-2s-%-2s-%-2s" % (k,j,i)
                else:
                    Type = "%-2s-%-2s-%-2s" % (i,j,k)
                if Type not in List:
                    T0, Keq = doAngle(bond1, bond2, Bonds, uff, debug)
                    if Type not in List or Type not in aList:
                        List.append(Type)
                        Angles.append("%s %f %f" % (Type, Keq, T0))
                    if Type not in aList:
                        aList.append(Type)
                        aAngles.append("%s %f %f" % (Type, Keq, T0))

                if bond1.GetBgnIdx() > bond2.GetBgnIdx():
                    I.append(bond2.GetBgnIdx())
                    J.append(bond1.GetEndIdx())
                    K.append(bond1.GetBgnIdx())
                else:
                    I.append(bond1.GetBgnIdx())
                    J.append(bond1.GetEndIdx())
                    K.append(bond2.GetBgnIdx())

    if debug == True:
        print "ANGLE"
        for i in xrange(len(Angles)):
            tok = string.split(Angles[i])
            if len(tok) == 5:
                print "%-2s%-3s%-3s \t %6.1f \t %8.3f" % (tok[0], tok[1],
        tok[2], float(tok[3]), float(tok[4]))
            elif len(tok) == 4:
                if len(tok[0]) > 2:
                    print "%-5s%-3s \t %6.1f \t %8.3f" % (tok[0], tok[1],
        float(tok[2]), float(tok[3]))
                else:
                    print "%-2s%-6s \t %6.1f \t %8.3f" % (tok[0], tok[1],
        float(tok[2]), float(tok[3]))
            else:
                print "%-8s \t %6.1f \t %8.3f" % (tok[0], float(tok[1]),
        float(tok[2]))
        print ""

    return Angles, aAngles, I, J, K
    
def doAngle(bond1, bond2, Bonds, uff, debug):
    """ Estimate equilibrium value and force constant for angle bending
    across input bonds.  Angle values will be generic.  Force constants can
    be either generic:  35/70 for H/noH angles, or estimated to be half
    those determined using the UFF formula
    Rappe et. al. JACS 1992, 114(25), 10024-10035, p 10028 center
    """

    Zeff = {0 : 0.000, 1 : 0.712, 2 : 0.098, 3 : 1.026, 4 : 1.565,
            5 : 1.755, 6 : 1.912, 7 : 2.544, 8 : 2.300, 9 : 1.735,
            10 : 0.194, 11 : 1.081, 12 : 1.787, 13 : 1.792, 14 : 2.323,
            15 : 2.863, 16 : 2.703, 17 : 2.348, 18 : 0.300, 19 : 1.165,
            20 : 2.141, 21 : 2.592, 22 : 2.659, 23 : 2.679, 24 : 2.463,
            25 : 2.430, 26 : 2.430, 27 : 2.430, 28 : 2.430, 29 : 1.756,
            30 : 1.308, 31 : 1.821, 32 : 2.789, 33 : 2.864, 34 : 2.764,
            35 : 2.519, 36 : 0.452, 37 : 1.592, 38 : 2.449, 39 : 3.257,
            40 : 3.667, 41 : 3.618, 42 : 3.400, 43 : 3.400, 44 : 3.400,
            45 : 3.508, 46 : 3.210, 47 : 1.956, 48 : 1.650, 49 : 2.070,
            50 : 2.961, 51 : 2.704, 52 : 2.882, 53 : 2.650, 54 : 0.556,
            55 : 1.573, 56 : 2.727, 57 : 3.300, 58 : 3.300, 59 : 3.300,
            60 : 3.300, 61 : 3.300, 62 : 3.300, 63 : 3.300, 64 : 3.300,
            65 : 3.300, 66 : 3.300, 67 : 3.416, 68 : 3.300, 69 : 3.300,
            70 : 2.618, 71 : 3.271, 72 : 3.921, 73 : 4.075, 74 : 3.700,
            75 : 3.700, 76 : 3.700, 77 : 3.731, 78 : 3.382, 79 : 2.625,
            80 : 1.750, 81 : 2.068, 82 : 2.846, 83 : 2.470, 84 : 2.330,
            85 : 2.240, 86 : 0.583, 87 : 1.847, 88 : 2.920, 89 : 3.900,
            90 : 4.202, 91 : 3.900, 92 : 3.900, 93 : 3.900, 94 : 3.900,
            95 : 3.900, 96 : 3.900, 97 : 3.900, 98 : 3.900, 99 : 3.900,
            100 : 3.900, 101 : 3.900, 102 : 3.900, 103 : 3.900, 104 : 3.900,
            105 : 3.900, 106 : 3.900, 107 : 3.900, 108 : 3.900, 109 : 3.900}
    
    Keq = 0.0
    T0 = 0.0

    if bond1.GetBgnIdx() == bond2.GetBgnIdx() and \
           bond1.GetEndIdx() != bond2.GetEndIdx():
        i = bond2.GetEnd().GetType().replace("_","")
        j = bond1.GetBgn().GetType().replace("_","")
        k = bond1.GetEnd().GetType().replace("_","")
        Iatom = bond2.GetEnd()
        Jatom = bond1.GetBgn()
        Katom = bond1.GetEnd()
    elif bond1.GetBgnIdx() == bond2.GetEndIdx() and \
             bond1.GetEndIdx() != bond2.GetBgnIdx():

        i = bond2.GetBgn().GetType().replace("_","")
        j = bond1.GetBgn().GetType().replace("_","")
        k = bond1.GetEnd().GetType().replace("_","")
        Iatom = bond2.GetBgn()
        Jatom = bond1.GetBgn()
        Katom = bond1.GetEnd()
    elif bond1.GetEndIdx() == bond2.GetBgnIdx() and \
             bond1.GetBgnIdx() != bond2.GetEndIdx():

        i = bond2.GetEnd().GetType().replace("_","")
        j = bond1.GetEnd().GetType().replace("_","")
        k = bond1.GetBgn().GetType().replace("_","")
        Iatom = bond2.GetEnd()
        Jatom = bond1.GetEnd()
        Katom = bond1.GetBgn()
    elif bond1.GetEndIdx() == bond2.GetEndIdx() and \
             bond1.GetBgnIdx() != bond2.GetBgnIdx():

        i = bond2.GetBgn().GetType().replace("_","")
        j = bond1.GetEnd().GetType().replace("_","")
        k = bond1.GetBgn().GetType().replace("_","")
        Iatom = bond2.GetBgn()
        Jatom = bond1.GetEnd()
        Katom = bond1.GetBgn()

    jDeg = Jatom.GetDegree()

    if jDeg == 2:
        if Jatom.GetAtomicNum() == 8:
            T0 = 105.0
        elif Jatom.GetHyb() == 1:
            T0 = 180.0
        elif Jatom.GetAtomicNum() >= 10:
            T0 = 95.0
        else:
            T0 = 120.0
    elif jDeg == 3:
        if Jatom.GetHyb() == 2:
            T0 = 120.0
        else:
            if Jatom.GetAtomicNum() == 7:
                T0 = 107.0
            else:
                T0 = 92.0
    elif jDeg == 4:
        T0 = 109.45
    else:
        T0 = 120.0

    """ Separate check for small rings
    OEDetermineRingSystems must have been called """
    if OEAtomGetSmallestRingSize(Jatom) == 3:
        T0 = 60.0
    elif OEAtomGetSmallestRingSize(Jatom) == 4:
        T0 = 90.0

    """ Generic Force constants """
    if uff == False:
        if Iatom.GetAtomicNum() == 1 or Katom.GetAtomicNum() == 1:
            Keq = 35.0
        else:
            Keq = 70.0

            """ UFF force constants """
    else:
        Zi = Zeff[Iatom.GetAtomicNum()]
        Zk = Zeff[Katom.GetAtomicNum()]
        Theta = T0 * pi/180.0

        if i > j:
            BTij = "%-2s-%-2s" % (j,i)
        else:
            BTij = "%-2s-%-2s" % (i,j)
        
        if j > k:
            BTjk = "%-2s-%-2s" % (k,j)
        else:
            BTjk = "%-2s-%-2s" % (j,k)

        rij = 0.0
        rjk = 0.0
        for i in xrange(len(Bonds)):
            temp = Bonds[i].split("-")
            temp2 = temp[1].split()
            Comp = "%-2s-%-2s" % (temp[0],temp2[0])
            if Comp == BTij:
                rij = float(temp2[2])
            if Comp == BTjk:
                rjk = float(temp2[2])
        Beta = 332.0/(rij *rjk)
        rik2 = (rij * rij) + (rjk * rjk) - (2.0 * rij * rjk * cos(Theta))
        rik = sqrt(rik2)
        rik5 = rik2 * rik2 * rik
        fac1 = Beta * Zi * Zk * rij * rjk / rik5
        fac2 = rij * rij * rij * rjk * (1.0 - cos(Theta) * cos(Theta))
        fac3 = rik2 * cos(Theta)
        Keq = fac1*(fac2-fac3)
            
    return T0, Keq

    
def uniqueBonds(mol, resList, List, debug):
    Bonds = []
    aList = []
    aBonds = []
    for bond in mol.GetBonds():
        i = bond.GetBgn().GetType().replace("_","")
        j = bond.GetEnd().GetType().replace("_","")
        if i > j:
            Type = "%-2s-%-2s" % (j,i)
        else:
            Type = "%-2s-%-2s" % (i,j)
        if Type not in List or Type not in aList:
            req, rk = doBond(mol, bond, resList, debug)
            if Type not in List:
                List.append(Type)
                Bonds.append("%s %f %f" % (Type, rk, req))
            if Type not in aList:
                aList.append(Type)
                aBonds.append("%s %f %f" % (Type, rk, req))
    if debug == True:
        print "BOND"
        for i in xrange(len(Bonds)):
            tok = string.split(Bonds[i])
            if len(tok) == 4:
                print "%-2s%-3s \t %7.2f \t %6.3f" % (tok[0], tok[1], float(tok[2]),
                                                      float(tok[3]))
            else:
                print "%-5s \t %7.2f \t %6.3f" % (tok[0], float(tok[1]),
                                                  float(tok[2]))
        print ""
                
    return Bonds, aBonds

def doBond(mol, bond, resList, debug):
    """ Generate bond lengths and force constants for input bond in the
    molecule. 

    Assumes explicit hydrogens have been added, hybridization
    determined, aromaticty model selected and applied, etc.

    Parameters take from: 
      Allinger, Norman L., Zhue, Xuefeng and Bergsma, John
      "Molecular mechanics parameters"
      J. Mol. Struct. (Theochem), 1994, 312, 69-83.
    """
    Cov = {0 : 0.000, 1 : 0.351, 2 : 0.000, 3 : 1.230, 4 : 0.890,
           5 : 0.880, 6 : 0.760, 7 : 0.710, 8 : 0.720, 9 : 0.720,
           10 : 1.310, 11 : 1.570, 12 : 1.170, 13 : 1.250, 14 : 1.170,
           15 : 1.100, 16 : 1.040, 17 : 0.990, 18 : 1.740, 19 : 2.030,
           20 : 1.740, 21 : 1.440, 22 : 1.320, 23 : 1.220, 24 : 1.170,
           25 : 1.170, 26 : 1.160, 27 : 1.160, 28 : 1.150, 29 : 1.170,
           30 : 1.250, 31 : 1.250, 32 : 1.220, 33 : 1.210, 34 : 1.440,
           35 : 1.140, 36 : 1.890, 37 : 2.160, 38 : 1.920, 39 : 1.620,
           40 : 1.450, 41 : 1.340, 42 : 1.290, 43 : 0.000, 44 : 1.240,
           45 : 1.250, 46 : 1.280, 47 : 1.340, 48 : 1.410, 49 : 1.500,
           50 : 1.400, 51 : 1.410, 52 : 1.370, 53 : 1.330, 54 : 2.090,
           55 : 2.350, 56 : 1.980, 57 : 1.690, 58 : 1.650, 59 : 1.650,
           60 : 1.640, 61 : 0.000, 62 : 1.660, 63 : 1.850, 64 : 1.610,
           65 : 1.590, 66 : 1.590, 67 : 1.580, 68 : 1.570, 69 : 1.560,
           70 : 1.700, 71 : 1.560, 72 : 1.440, 73 : 1.340, 74 : 1.300,
           75 : 1.280, 76 : 1.260, 77 : 1.260, 78 : 1.290, 79 : 1.340,
           80 : 1.440, 81 : 1.550, 82 : 1.540, 83 : 1.520, 84 : 1.530,
           85 : 0.000, 86 : 0.000, 87 : 0.000, 88 : 0.000, 89 : 0.000,
           90 : 0.000, 91 : 0.000, 92 : 0.000, 93 : 0.000, 94 : 0.000,
           95 : 0.000, 96 : 0.000, 97 : 0.000, 98 : 0.000, 99 : 0.000,
           100 : 0.000, 101 : 0.000, 102 : 0.000, 103 : 0.000, 104 : 0.000,
           105 : 0.000, 106 : 0.000, 107 : 0.000, 108 : 0.000, 109 : 0.000}

    Eneg = {0 : 0.000,  1 : 2.200,  2 : 0.000,  3 : 1.000,  4 : 1.500,
            5 : 2.000,  6 : 2.600,  7 : 3.050,  8 : 3.500,  9 : 3.900,
            10 : 0.000,  11 : 0.900,  12 : 1.200,  13 : 1.500,  14 : 1.900,
            15 : 2.150,  16 : 2.600,  17 : 3.150,  18 : 0.000,  19 : 0.800,
            20 : 1.000,  21 : 1.300,  22 : 1.500,  23 : 1.600,  24 : 1.600,
            25 : 1.500,  26 : 1.800,  27 : 1.800,  28 : 1.800,  29 : 1.900,
            30 : 1.600,  31 : 1.600,  32 : 1.900,  33 : 2.000,  34 : 2.450,
            35 : 2.850,  36 : 0.000,  37 : 0.800,  38 : 1.000,  39 : 1.300,
            40 : 1.600,  41 : 1.600,  42 : 1.800,  43 : 1.900,  44 : 2.200,
            45 : 2.200,  46 : 2.200,  47 : 1.900,  48 : 1.700,  49 : 1.700,
            50 : 1.800,  51 : 2.050,  52 : 2.300,  53 : 2.650,  54 : 0.000,
            55 : 0.700,  56 : 0.900,  57 : 1.100,  58 : 1.100,  59 : 1.100,
            60 : 1.100,  61 : 1.100,  62 : 1.100,  63 : 1.100,  64 : 1.100,
            65 : 1.100,  66 : 1.100,  67 : 1.100,  68 : 1.100,  69 : 1.100,
            70 : 1.100,  71 : 1.200,  72 : 1.300,  73 : 1.300,  74 : 1.700,
            75 : 1.900,  76 : 2.200,  77 : 2.200,  78 : 2.400,  79 : 1.900,
            80 : 1.800,  81 : 1.800,  82 : 1.900,  83 : 2.000,  84 : 2.200,
            85 : 0.000,  86 : 0.650,  87 : 0.900,  88 : 1.100,  89 : 1.300,
            90 : 1.500,  91 : 1.700,  92 : 1.300,  93 : 1.300,  94 : 1.300,
            95 : 1.300,  96 : 1.300,  97 : 1.300,  98 : 1.300,  99 : 1.300,
            100 : 1.300,  101 : 1.300,  102 : 1.300,  103 : 0.000,  104 : 0.000,
            105 : 0.000,  106 : 0.000,  107 : 2.200,  108 : 0.000,  109 : 0.000}

    cVal = 0.05
    bExp = 1.4
    bCorr = 0.008
    dbFix = 0.1
    tbFix = 0.17
    abFix = 0.04
    abpiFix = 0.075
    sp2fix = 0.03
    spfix = 0.08
    fExp = 0.75
    preFac = 1.0
    conVal = 0.3
    md2Kcal = 143.93

    chi1 = Eneg[bond.GetBgn().GetAtomicNum()]
    chi2 = Eneg[bond.GetEnd().GetAtomicNum()]
    e1 = pow(fabs(chi1-chi2),bExp)
    R1 = Cov[bond.GetBgn().GetAtomicNum()]
    R2 = Cov[bond.GetEnd().GetAtomicNum()]
    if bond.GetBgn().GetAtomicNum() == 7 or \
           bond.GetBgn().GetAtomicNum() == 8 or \
           bond.GetBgn().GetAtomicNum() == 16:
        R1lp = True
    else:
        R1lp = False
    if bond.GetEnd().GetAtomicNum() == 7 or \
           bond.GetEnd().GetAtomicNum() == 8 or \
           bond.GetEnd().GetAtomicNum() == 16:
        R2lp = True
    else:
        R2lp = False

    if bond.GetOrder() == 3:
        R1 -= tbFix
        R2 -= tbFix
    elif bond.GetOrder() == 2:
        if bond.IsAromatic() == True or bond.GetIdx() in resList:
            if R1lp == True:
                R1 -= abpiFix
            else:
                R1 -= abFix

            if R2lp == True:
                R2 -= abpiFix
            else:
                R2 -= abFix
        else:
            R1 -= dbFix
            R2 -= dbFix

    elif bond.GetOrder() == 1:
        if bond.IsAromatic() == True or bond.GetIdx() in resList:
            if R1lp == True:
                R1 -= abpiFix
            else:
                R1 -= abFix

            if R2lp == True:
                R2 -= abpiFix
            else:
                R2 -= abFix
        else:
            if bond.GetBgn().GetAtomicNum() == 7 and \
               bond.GetBgn().GetHyb() == 2 and \
               bond.GetBgn().GetDegree() == 3 and \
               bond.GetBgn().GetValence() == 3:
                funkyNb = True
            else:
                funkyNb = False
            if bond.GetEnd().GetAtomicNum() == 7 and \
               bond.GetEnd().GetHyb() == 2 and \
               bond.GetEnd().GetDegree() == 3 and \
               bond.GetEnd().GetValence() == 3:
                funkyNe = True
            else:
                funkyNe = False
                
            if bond.GetBgn().GetHyb() == 2 and funkyNb == False:
                if (bond.GetEnd().GetHyb() == 3 or funkyNe == True) and \
                   R2lp == True:
                    R1 -= abFix
                    R2 -= abFix
                else:
                    R1 -= sp2fix
            elif bond.GetEnd().GetHyb() == 2 and funkyNe == False:
                if (bond.GetBgn().GetHyb() == 3 or funkyNb == True) and \
                   R1lp == True:
                    R1 -= abFix
                    R2 -= abFix
                else:
                    R2 -= sp2fix
            if bond.GetBgn().GetHyb() == 1:
                R1 -= spfix
            if bond.GetEnd().GetHyb() == 1:
                R2 -= spfix

    Req = R1 + R2 - cVal * e1 - bCorr
    e2 = pow((chi1*chi2)/(Req*Req),fExp)
    Force = md2Kcal * preFac * e2 + conVal
    """
    validate = True
    if validate == True:
        try:
            Rexp = OEGetDistance(mol, bond.GetBgn(), bond.GetEnd())
            Err = fabs(Rexp-Req)/Rexp
        except:
            Rexp = 0.0
            Err = 100.0

        print "%s %s %d %s %.3f %.3f %.2f" % (bond.GetBgn().GetType().replace("_",""), \
                                              bond.GetEnd().GetType().replace("_",""), \
                                              bond.GetOrder(), \
                                              bond.IsAromatic(), \
                                              Req, Rexp, Err)
    """
    return Req, Force

def uniqueMasses(mol, List, debug):
    """ Atom tye list with atomic masses """
    """ Polarizations not estimated """
    Masses = []
    for atom in mol.GetAtoms():
        if atom.GetType().replace("_","") not in List:
            List.append(atom.GetType().replace("_",""))
            Masses.append("%s %f %f" % (atom.GetType().replace("_",""),
                                          OEGetAverageWeight(atom.GetAtomicNum()),0.0))
    if debug == True:
        print mol.GetTitle()
        print "MASS"
        for mass in Masses:
            tok = string.split(mass)
            print "%-2s \t %7.2f \t %7.3f" % (tok[0], float(tok[1]), float(tok[2]))
        print ""

    return Masses

def doResonant(mol):
    """ capture nitro and carboxylate """
    List = []

    for atom in mol.GetAtoms():
        oCount = 0
        for bond in atom.GetBonds():
            if atom.GetIdx() == bond.GetBgnIdx():
                if bond.GetEnd().GetAtomicNum() == 8:
                    oCount += 1
            else:
                if bond.GetBgn().GetAtomicNum() == 8:
                    oCount += 1

        if (atom.GetAtomicNum() == 7 or atom.GetAtomicNum() == 6) and \
                atom.GetDegree() == 3 and atom.GetHyb() == 2 and oCount == 2:
            for bond in atom.GetBonds():
                if atom.GetIdx() == bond.GetBgnIdx():
                    if bond.GetEnd().GetAtomicNum() == 8:
                        List.append(bond.GetIdx())
                else:
                    if bond.GetBgn().GetAtomicNum() == 8:
                        List.append(bond.GetIdx())

    return List

def fixType(inType):
    """ sort type string and return """
    tok = string.split(inType,"-")
    if len(tok) == 2:
        if tok[0].strip() > tok[1].strip():
            Type = "%-2s-%-2s" % (tok[1].strip(),tok[0].strip())
        else:
            Type = inType
    elif len(tok) == 3:
        if tok[0].strip() > tok[2].strip():
            Type = "%-2s-%-2s-%-2s" % (tok[2].strip(),tok[1].strip(),tok[0].strip())
        else:
            Type = inType
    elif len(tok) == 4:
        if tok[0].strip() == "X" and tok[3].strip() == "X":
            if tok[1].strip() > tok[2].strip():
                Type = "%-2s-%-2s-%-2s-%-2s" % (tok[3].strip(),tok[2].strip(),tok[1].strip(),tok[0].strip())
            else:
                Type = inType
        else:
            if tok[0].strip() > tok[3].strip():
                Type = "%-2s-%-2s-%-2s-%-2s" % (tok[3].strip(),tok[2].strip(),tok[1].strip(),tok[0].strip())
            else:
                Type = inType
        
    return Type

def main():

    """ Command line parsing information """
    desc = "Estimate bond lengths and force constants"
    vers = "OEparmchk v0.2 - 11/2010, Richard Dixon, Vertex"
    use = "%prog [options] FILE.mol2 (generated by OEante.py)"
    parser = OptionParser(description=desc, version=vers, usage=use)

    parser.add_option("-d", "--debug", action="store_true", dest="debug",
                      default=False,
                      help="Turn on debug printing [default: %default]")
    parser.add_option("-u", "--uff", action="store_true", dest="uff",
                      default=False,
                      help="Use UFF empirical angle bending force constants\
                      [default: %default], i.e. constant force value applied")
    parser.add_option("-f", "--force", action="store", dest="force",
                      default=0.0,
                      help="force constant for out of plane bending \
                      [default: %default, i.e. no oop force applied]")
    parser.add_option("-p", "--param", action="store", dest="param",
                      help="INPUT parameter file [optional]")
    parser.add_option("-o", "--outprefix", action="store", dest="out",
                      help="OUTPUT file prefix (for leaprc, frcmod and modified mol2 files)")

    (options, args) = parser.parse_args()

    if(len(args) != 1):
        parser.error("incorrect number of arguments")

    ifs = oemolistream(args[0])
    if options.out == None:
        Prefix,ext = os.path.splitext(args[0])
    else:
        Prefix = options.out
        
    ofs = oemolostream()
    oName = Prefix+"_pchk.mol2"
    ofs.open(oName)

    if options.param == None:
        frcName = "ante_"+Prefix+".frcmod"
        leapName = "ante_"+Prefix+".leaprc"
    else:
        Pre2,ext = os.path.splitext(options.param)
        frcName = Pre2+"_"+Prefix+".frcmod"
        leapName = Pre2+"_"+Prefix+".leaprc"

    Frc = open(frcName,'w')
    Leap = open(leapName,'w')

    """ Get force field parameters if offered """
    Mass = {}
    Polar = {}
    Req = {}
    Rk = {}
    Teq = {}
    Tk = {}
    atHybrid = {}
    atSym = {}
    Tor_numB = {}
    Tor_Vx2 = {}
    Tor_Phase = {}
    Tor_xVal = {}
    Imp_Vx2 = {}
    Imp_Phase = {}
    Imp_xVal = {}
    Rstar = {}
    Eps = {}
    Header=""
    MassPrint = []
    MassList = []
    BondPrint = []
    BondList = []
    AnglePrint = []
    AngleList = []
    TorPrint = []
    TorList = []
    ImpPrint = []
    ImpList = []
    VdWPrint = []
    VdWList = []
    Footer = []
    
    if options.param != None:
        ff = open(options.param, 'r')
        First = True
        for line in ff:
            if First == True:
                First = False
                Header = line
            else:
                if line[2:3] == " " and line[3:4] != " ":
                    tok = string.split(line)
                    try:
                        Mass[tok[0]] = float(tok[1])
                        Polar[tok[0]] = float(tok[2])
                        if tok[0] not in MassList:
                            MassList.append(tok[0])
                    except:
                        continue
                elif line[2:3] == "-" and line[5:6] == " ":
                    Type = fixType(line[0:5])
                    if Type not in BondList:
                        BondList.append(Type)
                    ExtraSpace = len(string.split(Type))-1
                    tok = string.split(line)
                    try:
                        Rk[Type] = float(tok[1+ExtraSpace])
                        Req[Type] = float(tok[2+ExtraSpace])
                    except:
                        continue
                elif line[2:3] == "-" and line[5:6] == "-" and \
                     line[8:9] == " ":
                    Type = fixType(line[0:8])
                    if Type not in AngleList:
                        AngleList.append(Type)
                    ExtraSpace = len(string.split(Type))-1
                    tok = string.split(line)
                    try:
                        Tk[Type] = float(tok[1+ExtraSpace])
                        Teq[Type] = float(tok[2+ExtraSpace])
                    except:
                        continue
                elif line[2:3] == "-" and line[5:6] == "-" and \
                     line[8:9] == "-" and line[11:12] == " ":

                    try:
                        int(line[12:17].strip())
                        Type = fixType(line[0:11])
                        if Type not in TorList:
                            TorList.append(Type)
                        ExtraSpace = len(string.split(Type))-1
                        tok = string.split(line)
                        try:
                            Tor_numB[Type] = int(tok[1+ExtraSpace])
                            Tor_Vx2[Type] = float(tok[2+ExtraSpace])
                            Tor_Phase[Type] = float(tok[3+ExtraSpace])
                            Tor_xVal[Type] = float(tok[4+ExtraSpace])
                        except:
                            continue
                    except:
                        Type = line[0:11] # not fixing type name, order dependent
                        if Type not in ImpList:
                            ImpList.append(Type)
                        ExtraSpace = len(string.split(Type))-1
                        tok = string.split(line)
                        try:
                            Imp_Vx2[Type] = float(tok[1+ExtraSpace])
                            Imp_Phase[Type] = float(tok[2+ExtraSpace])
                            Imp_xVal[Type] = float(tok[3+ExtraSpace])
                        except:
                            continue
                        
        ff.seek(0)
        vStart = re.compile("^MOD4")
        vEnd = re.compile("^END")
        Open = False
        Foot = False
        for line in ff:
            if Foot == True:
                Footer.append(line)     # In case I want to reconstruct the input file
            if vStart.search(line):
                Open = True
            if vEnd.search(line):
                Open = False
                Foot = True
            if Open == True and line[2:3] != " ":
                tok = string.split(line)
                try:
                    if tok[0] not in VdWList:
                        VdWList.append(tok[0])
                    Rstar[tok[0]] = float(tok[1])
                    Eps[tok[0]] = float(tok[2])
                except:
                    continue
        ff.close()

    """ Construct printing lists """
    for k,v in Mass.iteritems():
        MassPrint.append("%-2s \t %7.2f \t %7.3f \t %s\n" % (k, v, Polar[k], options.param))
    for k,v, in Rk.iteritems():
        BondPrint.append("%-5s \t %7.2f \t %6.3f \t %s\n" % (k, v, Req[k], options.param))
    for k,v, in Tk.iteritems():
        AnglePrint.append("%-8s \t %6.1f \t %8.3f \t %s\n" % (k, v, Teq[k], options.param))
    for k,v, in Tor_numB.iteritems():
        TorPrint.append("%-11s \t %2d \t %6.3f \t %8.3f \t %6.3f \t %s\n" %
                        (k, v, Tor_Vx2[k], Tor_Phase[k], Tor_xVal[k], options.param))
    for k,v, in Imp_Vx2.iteritems():
        ImpPrint.append("%-11s \t %6.3f \t %8.3f \t %6.3f \t %s\n" %
                        (k, v, Imp_Phase[k], Imp_xVal[k], options.param))
    for k,v, in Rstar.iteritems():
        VdWPrint.append("%4s \t %7.4f \t %7.4f \t %s\n" % (k, v, Eps[k], options.param))

    for molIndex,mol in enumerate(ifs.GetOEMols()):
        
        OEAssignAromaticFlags(mol)
        OEAssignHybridization(mol)
        OEAddExplicitHydrogens(mol,False,True)
        Nring = OEDetermineRingSystems(mol)

        resList = doResonant(mol)

        Masses = uniqueMasses(mol, MassList, options.debug)
        Bonds,aBonds = uniqueBonds(mol, resList, BondList, options.debug)
        Angles, aAngles, I, J, K = uniqueAngles(mol, aBonds, options.uff,
                                                AngleList, options.debug)
        Torsions = uniqueTorsions(mol, aAngles, I, J, K, TorList, options.debug)
        Impropers = uniqueImpropers(mol, float(options.force), ImpList, options.debug)
        VdW = uniqueVdw(mol, VdWList, options.debug)

        """ Harvest atom types, symbols and hybridizations for leaprc """
        for atom in mol.GetAtoms():
            iTemp = atom.GetHyb()
            if iTemp == 1 or iTemp == 2 or atom.GetType().replace("_","") == "n":
                atHybrid[atom.GetType().replace("_","")] = "sp2"
            else:
                atHybrid[atom.GetType().replace("_","")] = "sp3"

            atSym[atom.GetType().replace("_","")] = OEGetAtomicSymbol(atom.GetAtomicNum())

        """ Add newly found parameters to existing set """
        for mass in Masses:
            tok = string.split(mass)
            MassPrint.append("%-2s \t %7.2f \t %7.3f \t *Antechamber_generated*\n" % \
                             (tok[0], float(tok[1]), float(tok[2])))
        del Masses[:]
        for bond in Bonds:
            tok = string.split(bond)
            BondPrint.append("%-5s \t %7.2f \t %6.3f \t *Antechamber_generated*\n" % \
                             (bond[0:6], float(tok[-2]), float(tok[-1])))
        del Bonds[:]
        for angle in Angles:
            tok = string.split(angle)
            AnglePrint.append("%-8s \t %6.1f \t %8.3f \t *Antechamber_generated*\n" % \
                              (angle[0:9], float(tok[-2]), float(tok[-1])))
        del Angles[:]
        for tor in Torsions:
            tok = string.split(tor)
            TorPrint.append("%-11s \t %2d \t %6.3f \t %8.3f \t %6.3f \t *Antechamber_generated*\n" % \
                            (tor[0:12], int(tok[-4]), float(tok[-3]), float(tok[-2]), float(tok[-1])))
        del Torsions[:]
        for imp in Impropers:
            tok = string.split(imp)
            ImpPrint.append("%-11s \t %6.3f \t %8.3f \t %6.3f \t *Antechamber_generated*\n" % \
                            (imp[0:12], float(tok[-3]), float(tok[-2]), float(tok[-1])))
        del Impropers[:]
        for vdw in VdW:
            tok = string.split(vdw)
            VdWPrint.append("%4s \t %7.4f \t %7.4f \t *Antechamber_generated*\n" % \
                            (tok[0], float(tok[1]), float(tok[2])))
        del VdW[:]

        """ remove _ from atom type names """
        for atom in mol.GetAtoms():
            tTemp = atom.GetType()
            tNew = tTemp.replace("_","")
            atom.SetType(tNew)
        OEWriteMol2File(ofs,mol,True)

    """ Construct and write frcmod """
    Frc.write("Antechamber frcmod file generated for %s\n" % (args[0]))
    Frc.write("MASS\n")
    for line in MassPrint:
        Frc.write(line)
    Frc.write("\nBOND\n")
    for line in BondPrint:
        Frc.write(line)
    Frc.write("\nANGLE\n")
    for line in AnglePrint:
        Frc.write(line)
    Frc.write("\nDIHE\n")
    for line in TorPrint:
        Frc.write(line)
    Frc.write("\nIMPROPERS\n")
    for line in ImpPrint:
        Frc.write(line)
    Frc.write("\nNONBON\n")
    for line in VdWPrint:
        Frc.write(line)
    Frc.write("\n\n")

    """ Construct and write leaprc file """
    Leap.write("logfile leapante.log\n")
    Leap.write("addAtomTypes {\n")
    for k,v in atSym.iteritems():
        Leap.write("       { \"%s\" \t \"%s\" \t \"%s\" }\n" % (k,v,atHybrid[k]))
    Leap.write("}\n")
    Leap.write("ligfrc = loadAmberParams %s\n" % (frcName))
    Leap.write("lig = loadMol2 %s\n" % (oName))
    Leap.write("saveAmberParm lig lig.top lig.crd\n")
    Leap.write("quit\n")
    
    ifs.close()
    ofs.close()
    Frc.close()
    Leap.close()

if __name__ == "__main__":
    main()

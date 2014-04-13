#!/usr/bin/env python

""" OEante:  AMBER force field and program atom-typing and
    parameter generation.  Richard Dixon, Vertex Pharma

    Copyright (2010) Vertex Pharma.
    This software is distributed under the terms of the
    GNU General Public License
    """

import sys, os, string
from openeye.oechem import *
from openeye.oeomega import *
from openeye.oequacpac import *
import math
import re
from math import sqrt
from optparse import OptionParser
from itertools import *

class MatchSet:
    def __init__(self,matchSet):
        self.matchSet = matchSet
        self.valid = True

    def isValid(self):
        return self.valid

    def compare(self,other):
        if self.matchSet.issubset(other.matchSet):
            self.valid = False

    def __repr__(self):
        o=""
        for key in self.matchSet:
           o+= str(key)+" "
        o+= str(self.valid)
        return(o) 

class Matcher:
    def __init__(self,smartsLines):
        self.patList = []
        for line in smartsLines:
            toks = string.split(line)
            if len(toks) == 1:
                smarts = toks[0]
                pat=OESubSearch()
                if not pat.Init(smarts):
                    OEThrow.Fatal('Bad smarts: '+smarts)
                else:
                    self.patList.append([smarts,pat])

    def match(self,mol):
        matchList = []
        for smarts,pat in self.patList:
            hitList = []
            for match in pat.Match(mol,1):
                hitSet = set()
                for matchpair in match.GetAtoms():
                    hitSet.add(matchpair.target.GetIdx())
                hitList.append(MatchSet(hitSet))
            matchList.append([smarts,hitList])
        return matchList

    def refineSets(self,matchList):
        for a in matchList:
            smartsA,listA = a
            for b in matchList:
                smartsB,listB = b
                if a != b:
                    for setA in listA:
                        for setB in listB:
                            setA.compare(setB)
        return matchList

def RingParse(mol):

    Rchange = 0

    Nring = OEDetermineRingSystems(mol)
    if Nring[0] == 0:
        return Rchange

    for atom in mol.GetAtoms():

        if atom.IsInRing():
            if OEAtomGetSmallestRingSize(atom) == 3:
                if atom.GetType() == "c3":
                    atom.SetType("cx")
                    Rchange += 1
                if atom.GetType() == "c2":
                    atom.SetType("cu")
                    Rchange += 1
            if OEAtomGetSmallestRingSize(atom) == 4:
                if atom.GetType() == "c3":
                    atom.SetType("cy")
                    Rchange += 1
                if atom.GetType() == "c2":
                    atom.SetType("cv")
                    Rchange += 1
    return Rchange

def SetMulti(SDict, smarts, lst, mol):

    IDlist = []
    del IDlist[:]
    id = string.split(str(lst))
    for i in xrange(len(id)-1):
        IDlist.append(int(id[i]))
    pat = OESubSearch()
    pat.Init(smarts)
    for match in pat.Match(mol,1):
        count = 0
        for mp in match.GetAtoms():
            if mp.target.GetIdx() in IDlist:
                count += 1
        if count == match.NumAtoms():
            labs = string.split(SDict[smarts],";")
            i = -1
            for mp in match.GetAtoms():
                i += 1
                atom = mol.GetAtom(OEHasAtomIdx(mp.target.GetIdx()))
                atom.SetType(labs[i])

def combTypes(mol):
    """ enumerate combinations of multitype atoms """
    masterTypes = []
    aList = []

    for atom in mol.GetAtoms():
        Type = atom.GetType()

        if Type == "cp" or Type == "cq":
            masterTypes.append(["cp","cq"])
            aList.append(atom.GetIdx())

        if Type == "cc" or Type == "cd":
            masterTypes.append(["cc","cd"])
            aList.append(atom.GetIdx())

        if Type == "ce" or Type == "cf":
            masterTypes.append(["ce","cf"])
            aList.append(atom.GetIdx())

        if Type == "cg" or Type == "ch":
            masterTypes.append(["cg","ch"])
            aList.append(atom.GetIdx())

        if Type == "nc" or Type == "nd":
            masterTypes.append(["nc","nd"])
            aList.append(atom.GetIdx())

        if Type == "ne" or Type == "nf":
            masterTypes.append(["ne","nf"])
            aList.append(atom.GetIdx())

        if Type == "pc" or Type == "pd":
            masterTypes.append(["pc","pd"])
            aList.append(atom.GetIdx())

        if Type == "pe" or Type == "pf":
            masterTypes.append(["pe","pf"])
            aList.append(atom.GetIdx())

    for x in product(*masterTypes):
        for i in xrange(len(aList)):
            fix = mol.GetAtom(OEHasAtomIdx(aList[i]))
            fix.SetType(x[i])
        OK = checkTypes(mol)
        if OK == 0:
            Score = scoreTypes(mol)
            if Score == 0:
                break

def checkTypes(mol):
    """ Identify cases of mismatched bond types """
    nCon = 0
    tDict = {}

    for k,v in tDict.iteritems():
        del tDict[k]
    for bond in mol.GetBonds():
        t1 = bond.GetBgn().GetType()
        t2 = bond.GetEnd().GetType()

        if t1 > t2:
            tt = t1
            t1 = t2
            t2 = tt

        if tDict.get(t1) == None:
            tDict[t1] = {}
        currentBO = tDict[t1].get(t2)
        if currentBO == None or currentBO == bond.GetType():
            tDict[t1][t2] = bond.GetType()
        else:
            i = bond.GetBgnIdx()
            j = bond.GetEndIdx()
            bo = bond.GetOrder()
            nCon += 1

    return nCon

def scoreTypes(mol):
    """ Score atom type combinations for adherence to AMBER rules.
    Rules are like types across single bonds, unlike types across
    higher order bonds. """

    Score = 0
    for bond in mol.GetBonds():
        t1 = bond.GetBgn().GetType()
        t2 = bond.GetEnd().GetType()
        bo = bond.GetType()

        if t1 > t2:
            tt = t1
            t1 = t2
            t2 = tt

        if t1 == "cp" and t2 == "cp":
            if bo != "1":
                Score += 1
        elif t1 == "cq" and t2 == "cq":
            if bo != "1":
                Score += 1
        elif t1 == "cp" and t2 == "cq":
            if bo == "1":
                Score += 1

        elif t1 == "cc" and t2 == "cc":
            if bo != "1":
                Score += 1
        elif t1 == "cd" and t2 == "cd":
            if bo != "1":
                Score += 1
        elif t1 == "cc" and t2 == "cd":
            if bo == "1":
                Score += 1

        elif t1 == "ce" and t2 == "ce":
            if bo != "1":
                Score += 1
        elif t1 == "cf" and t2 == "cf":
            if bo != "1":
                Score += 1
        elif t1 == "ce" and t2 == "cf":
            if bo == "1":
                Score += 1

        elif t1 == "cg" and t2 == "cg":
            if bo != "1":
                Score += 1
        elif t1 == "ch" and t2 == "ch":
            if bo != "1":
                Score += 1
        elif t1 == "cg" and t2 == "ch":
            if bo == "1":
                Score += 1

        elif t1 == "nc" and t2 == "nc":
            if bo != "1":
                Score += 1
        elif t1 == "nd" and t2 == "nd":
            if bo != "1":
                Score += 1
        elif t1 == "nc" and t2 == "nd":
            if bo == "1":
                Score += 1

        elif t1 == "ne" and t2 == "ne":
            if bo != "1":
                Score += 1
        elif t1 == "nf" and t2 == "nf":
            if bo != "1":
                Score += 1
        elif t1 == "ne" and t2 == "nf":
            if bo == "1":
                Score += 1

        elif t1 == "pc" and t2 == "pc":
            if bo != "1":
                Score += 1
        elif t1 == "pd" and t2 == "pd":
            if bo != "1":
                Score += 1
        elif t1 == "pc" and t2 == "pd":
            if bo == "1":
                Score += 1

        elif t1 == "pe" and t2 == "pe":
            if bo != "1":
                Score += 1
        elif t1 == "pf" and t2 == "pf":
            if bo != "1":
                Score += 1
        elif t1 == "pe" and t2 == "pf":
            if bo == "1":
                Score += 1
    return Score

def fixUnknown(mol):
    """ Assign generic types to remaining UNK atoms.  The scheme is based
    on number of connections and charge for heavy atoms and connected atom
    for hydrogen """

    unkCount = 0
    
    for atom in mol.GetAtoms():
        if atom.GetType() == "UNK":
            unkCount += 1
            aNum = atom.GetAtomicNum()
            aSym = OEGetAtomicSymbol(aNum).lower()
            aConn = atom.GetDegree()
            aChg = atom.GetFormalCharge()
            if aConn == 0:
                if aChg == 0:
                    aPost = "A"
                elif aChg < 0:
                    aPost = "B"
                else:
                    aPost = "C"
            elif aConn == 1:
                if aChg == 0:
                    aPost = "D"
                elif aChg < 0:
                    aPost = "E"
                else:
                    aPost = "F"
            elif aConn == 2:
                if aChg == 0:
                    aPost = "G"
                elif aChg < 0:
                    aPost = "H"
                else:
                    aPost = "I"
            elif aConn == 3:
                if aChg == 0:
                    aPost = "J"
                elif aChg < 0:
                    aPost = "K"
                else:
                    aPost = "L"
            elif aConn == 4:
                if aChg == 0:
                    aPost = "M"
                elif aChg < 0:
                    aPost = "N"
                else:
                    aPost = "O"
            elif aConn == 5:
                if aChg == 0:
                    aPost = "P"
                elif aChg < 0:
                    aPost = "Q"
                else:
                    aPost = "R"
            elif aConn == 6:
                if aChg == 0:
                    aPost = "S"
                elif aChg < 0:
                    aPost = "T"
                else:
                    aPost = "U"
            else:
                aPost = "Z"
                
            newType = aSym+aPost
            atom.SetType(newType)

    return(unkCount)

def oeante(input_filename, typefile='gaffsmarts.txt', nocharges=False, debug=False, mol2file='look.mol2'):
    """
    Parse the specified input file (smi, sdf, mol2) for molecules to parameterize.

    Parameters
    ----------
    input_filename : string
        An OpenEye-readable file (e.g. smi, sdf, mol2) containing one or more molecules to parse.
    typefile : string, optional, default='gaffsmarts.txt'
        Text file containing SMARTS definitions and corresponding GAFF atom types.
    nocharges : bool, optional, default=False
        If False, no charges will be created.
    debug : bool, optional, default=False
        If True, debug output will be printed.
    mol2file : string, optional, default='look.mol2'
        Output file?

    """


    Smarts = {}
    ft = open(options.typefile, 'r')
    Comment = re.compile("^\"")
    for line in ft:
        if not Comment.search(line):
            fields = string.split(line)
            Smarts[fields[0]] = fields[1]
    ft.close()
    SList = []
    for k,v in Smarts.iteritems():
        SList.append(k)
    matcher = Matcher(SList)

    ifs = oemolistream(input_filename)
    ofs = oemolostream()
    ofs.open(options.mol2file)

    omega = OEOmega()
    omega.SetMaxConfs(1)

    for molIndex,mol in enumerate(ifs.GetOEMols()):
        OEAssignAromaticFlags(mol, OEAroModelTripos)
        OEAssignHybridization(mol)
        OEAddExplicitHydrogens(mol,False,True)

        """ check for presence of coordinates and generate if charges
        requested and coordinates missing """

        if options.nocharges == False:

            if mol.GetDimension() != 3:
                log = oeostream()
                OEThrow.SetOutputStream(log)
                quiet= omega(mol)

            OEAssignPartialCharges(mol,OECharges_AM1BCC)

        matchList = matcher.match(mol)
        matchList = matcher.refineSets(matchList)

        OETriposAtomTypeNames(mol)
        OETriposBondTypeNames(mol)
        OETriposAtomNames(mol)
        for atom in mol.GetAtoms():
            atom.SetType("UNK")

        """ set intitial types based on Tripos aromatization rules """
        for smarts,lst in matchList:
            if options.debug == True:
                print smarts,lst
            for x in lst:
                if x.isValid():
                    id = string.split(str(x))
                    if len(id) == 2:
                        atom = mol.GetAtom(OEHasAtomIdx(int(id[0])))
                        atom.SetType(Smarts[smarts])
                    else:
                        SetMulti(Smarts, smarts, x, mol)

        """ set OE aromatization rules and correct types as needed """
        OEAssignAromaticFlags(mol, OEAroModelOpenEye)
        matchList = matcher.match(mol)
        matchList = matcher.refineSets(matchList)
        if options.debug == True:
            print "Round 2"
        for smarts,lst in matchList:
            if options.debug == True:
                print smarts,lst
            for x in lst:
                if x.isValid():
                    id = string.split(str(x))
                    if len(id) == 2:
                        atom = mol.GetAtom(OEHasAtomIdx(int(id[0])))
                        if atom.GetType() != Smarts[smarts]:
                            if (atom.GetType() == "c2" or atom.GetType() == "UNK") and \
                                    (Smarts[smarts] == "ca" or Smarts[smarts] == "cp"):
                                atom.SetType("cc")
                            elif atom.GetType() == "c2" and Smarts[smarts] == "cp":
                                atom.SetType("cc")
                            elif Smarts[smarts] == "na":
                                atom.SetType(Smarts[smarts])
                            elif atom.GetType() == "n2" and atom.IsAromatic():
                                atom.SetType("nc")
                            elif atom.GetType() == "ne" and atom.IsAromatic():
                                atom.SetType("nc")
                            elif atom.GetType() == "ce" and atom.IsAromatic():
                                atom.SetType("cc")
                            elif atom.GetType() == "n3" and Smarts[smarts] == "nh":
                                atom.SetType("nh")
                            elif (atom.GetType() == "p2" or atom.GetType() == "pe") \
                                     and Smarts[smarts] == "pb":
                                atom.SetType("pc")

        """ Fix cc/ce inconsistency """
        for atom in mol.GetAtoms():
            if atom.GetType() == "cc" and not atom.IsAromatic():
                atom.SetType("ce")

        """ ID 3 and 4 membered rings and change types accordingly """
        numRing = RingParse(mol)

        """ Generate generic name for any untyped atoms.  Used mixed case for these names """
        numUnknown = fixUnknown(mol)
        if options.debug == True:
            print "Number of unknown atoms detected = ", numUnknown

        """ Enumerate combinations of equivalent atom types and select set that obeys rules """
        combTypes(mol)

        """ Save output for viewing and OEparmchk.  Prepend _ to type names
        to avoid ca == calcium problem """

        for atom in mol.GetAtoms():
            tTemp = atom.GetType()
            tNew = "_"+tTemp
            atom.SetType(tNew)

        OEWriteMol2File(ofs,mol,True)

    return

def main():

    """ Command line parsing information """
    desc = "Set GAFF atom types via SMARTS pattern matching"
    vers = "OEante v0.2 - 11/2010, Richard Dixon, Vertex"
    use = "%prog [options] FILE.[smi,sdf,mol2]"
    parser = OptionParser(description=desc, version=vers, usage=use)

    parser.add_option("-t", "--typefile", action="store", dest="typefile",
                      default="gaffsmarts.txt",
                      help="input atom type data file in format \
                      SMARTS Name [default: %default]")
    parser.add_option("-n", "--nocharges", action="store_true", dest="nocharges",
                      default=False,
                      help="Turn off charge calculation [default: %default]")
    parser.add_option("-d", "--debug", action="store_true", dest="debug",
                      default=False,
                      help="Turn on debug printing [default: %default]")
    parser.add_option("-m", "--mol2file", action="store", dest="mol2file",
                      default="look.mol2",
                      help="output mol2 format file for visualization \
                      and input to OEparmchk [default: %default]")

    (options, args) = parser.parse_args()

    if(len(args) != 1):
        parser.error("incorrect number of arguments")

    oeante(args[0], options.typefile, options.nocharges, options.debug, options.mol2file)

if __name__ == "__main__":
    main()

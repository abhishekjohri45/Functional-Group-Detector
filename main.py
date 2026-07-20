"""
=========================================================
Functional Group Detector

A beginner-friendly RDKit project that detects
common functional groups from a SMILES string.

=========================================================
"""

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


# SMARTS Patterns

alcohol = Chem.MolFromSmarts("[OX2H]")

aldehyde = Chem.MolFromSmarts("[CX3H1](=O)")

ketone = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")

carboxylic_acid = Chem.MolFromSmarts("C(=O)[OX2H1]")

ester = Chem.MolFromSmarts("C(=O)O[#6]")

ether = Chem.MolFromSmarts("[#6]-O-[#6]")

amine = Chem.MolFromSmarts("[NX3;H2,H1,H0]")

amide = Chem.MolFromSmarts("C(=O)N")

benzene = Chem.MolFromSmarts("c1ccccc1")


def detect_functional_groups(molecule):

    found = False

    if molecule.HasSubstructMatch(alcohol):

        print("Alcohol")

        found = True

    if molecule.HasSubstructMatch(aldehyde):

        print("Aldehyde")

        found = True

    if molecule.HasSubstructMatch(ketone):

        print("Ketone")

        found = True

    if molecule.HasSubstructMatch(carboxylic_acid):

        print("Carboxylic Acid")

        found = True

    if molecule.HasSubstructMatch(ester):

        print("Ester")

        found = True

    if molecule.HasSubstructMatch(ether):

        print("Ether")

        found = True

    if molecule.HasSubstructMatch(amine):

        print("Amine")

        found = True

    if molecule.HasSubstructMatch(amide):

        print("Amide")

        found = True

    if molecule.HasSubstructMatch(benzene):

        print("Benzene Ring")

        found = True

    if found == False:

        print("No functional groups detected.")


while True:

    print()

    print("=" * 55)
    print("         FUNCTIONAL GROUP DETECTOR")
    print("=" * 55)

    print("1. Detect Functional Groups")
    print("2. Molecular Formula")
    print("3. Molecular Weight")
    print("4. Show All Information")
    print("5. Exit")

    print()

    choice = input("Enter your choice : ")

    if choice == "5":

        print()

        print("Thank you for using Functional Group Detector.")

        break

    print()

    smiles = input("Enter SMILES : ")

    molecule = Chem.MolFromSmiles(smiles)

    if molecule == None:

        print()

        print("Invalid SMILES.")

        continue

    # ----------------------------------------------------

    if choice == "1":

        print()

        print("Functional Groups Found")

        print("-" * 30)

        detect_functional_groups(molecule)

    # ----------------------------------------------------

    elif choice == "2":

        formula = rdMolDescriptors.CalcMolFormula(molecule)

        print()

        print("Molecular Formula")

        print("-" * 30)

        print(formula)

    # ----------------------------------------------------

    elif choice == "3":

        molecular_weight = Descriptors.MolWt(molecule)

        exact_weight = Descriptors.ExactMolWt(molecule)

        print()

        print("Molecular Weight")

        print("-" * 30)

        print("Average Molecular Weight :", round(molecular_weight, 2), "g/mol")

        print("Exact Molecular Weight :", round(exact_weight, 4))

    # ----------------------------------------------------

    elif choice == "4":

        formula = rdMolDescriptors.CalcMolFormula(molecule)

        molecular_weight = Descriptors.MolWt(molecule)

        exact_weight = Descriptors.ExactMolWt(molecule)

        print()

        print("=" * 40)

        print("MOLECULAR INFORMATION")

        print("=" * 40)

        print("Formula :", formula)

        print("Average Molecular Weight :", round(molecular_weight, 2), "g/mol")

        print("Exact Molecular Weight :", round(exact_weight, 4))

        print()

        print("Functional Groups")

        print("-" * 30)

        detect_functional_groups(molecule)

    # ----------------------------------------------------

    else:

        print()

        print("Invalid Choice.")

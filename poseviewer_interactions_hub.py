"""
Detect interactions between a receptor and a set of ligands in a poseviewer
or a complex file. The script lists each H-bond and contact (good, bad, or ugly)
along with a summary reporting the number of ligands forming each kind of
interaction with each receptor residue. The script also detects hydrophobic,
salt-bridge, pi-cation, pi-pi interactions, water-mediated h-bonds (single H2O
molecule only), and metal contacts.

<file> and <file2> are structure files in Maestro, SD, Mol2, or PDB formats. If
only a single file is specified (e.g., a pose-viewer file), interactions are
identified between the first structure and each subsequent structure. If two
files are specified, interactions are identified between the (first)
structure in the first file and each structure in the second file.

Copyright Schrodinger, LLC. All rights reserved.
"""

import argparse
import csv
import itertools
import os
import sys
from collections import OrderedDict
from collections import defaultdict

from schrodinger import adapter
from schrodinger import structure
from schrodinger.infra import mm
from schrodinger.structutils import analyze
from schrodinger.structutils import interactions
from schrodinger.structutils import measure
from schrodinger.structutils.interactions import hbond as hbond_module
from schrodinger.structutils.interactions.salt_bridge import get_salt_bridges
from schrodinger.utils import cmdline
from schrodinger.utils import csv_unicode
from schrodinger.utils import fileutils
from schrodinger.utils import log

logger = log.get_output_logger('poserviewer_interactions.py')

HBONDS = 'Hydrogen Bonds'
XBONDS = 'Halogen Bonds'
GOOD = 'Good Contacts'
BAD = 'Bad Contacts'
UGLY = 'Ugly Contacts'
SALT = 'Salt Bridges'
PICAT = 'Pi-Cation'
PIPI_FACE = 'Pi-pi Face-to-Face'
PIPI_EDGE = 'Pi-pi Edge-to-Face'
HPHOB = 'Hydrophobic'
WATER_HBONDS = 'Water-mediated H-Bonds'
AROMATIC_HBONDS = 'Aromatic H-Bonds'
METAL = 'Metal Contacts'

CONTACT_PROPS = {
    HBONDS: 'i_pvi_total_Hydrogen_Bonds',
    XBONDS: 'i_pvi_total_Halogen_Bonds',
    GOOD: 'i_pvi_total_Good_Contacts',
    BAD: 'i_pvi_total_Bad_Contacts',
    UGLY: 'i_pvi_total_Ugly_Contacts',
    SALT: 'i_pvi_total_Salt_Bridges',
    PICAT: 'i_pvi_total_Pi-Cation',
    PIPI_EDGE: 'i_pvi_total_Pi-Pi_Edge',
    PIPI_FACE: 'i_pvi_total_Pi-Pi_Face',
    HPHOB: 'i_pvi_total_Hydrophobic',
    WATER_HBONDS: 'i_pvi_total_Water_Mediated_HBonds',
    AROMATIC_HBONDS: 'i_pvi_total_Water_Aromatic_HBonds',
    METAL: 'i_pvi_total_Metal_Contacts',
}

# Short version of the contact type strings, for display in the Type column.
SHORTFORM = {
    HBONDS: "HBond",
    XBONDS: "XBond",
    UGLY: "Ugly",
    BAD: "Bad",
    GOOD: "Good",
    SALT: "Salt",
    HPHOB: "HPhob",
    PICAT: "PiCat",
    PIPI_FACE: "PiFace",
    PIPI_EDGE: "PiEdge",
    WATER_HBONDS: "Wat-HBond",
    AROMATIC_HBONDS: "Ar-HBond",
    METAL: "Metal",
}

# Values accepted by the -contacts option:
OPTION_VALUES = {
    HBONDS: 'hbond',
    XBONDS: 'xbond',
    GOOD: 'good',
    BAD: 'bad',
    UGLY: 'ugly',
    SALT: 'salt',
    PICAT: 'picat',
    PIPI_FACE: 'pipi',
    PIPI_EDGE: 'pipi',
    HPHOB: 'hphob',
    WATER_HBONDS: 'water_hbonds',
    AROMATIC_HBONDS: 'arom_hbonds',
    METAL: 'metal',
}

HPHOBIC_ASL = "((SMARTS. \"[#6;+0;!$([#6]~[#7,#8,#9])]\") \
              OR (SMARTS. \"[s,S;v2X2;H0;+0]\") \
              OR (SMARTS. \"[$([SX3v3;H0;+0;!$([SX3v3;H0;+0]~[#7,#8,#9])])]\") \
              OR (SMARTS. \"[$([Cl,Br,I;v1X1;H0;+0]([!Cl&!Br&!I]))]\") \
              OR (SMARTS. \"[B]\"))"

# Adding default settings in LID as of 20-3 (PYAPP-8188)
# The metal interaction is based only on distance and is different
# from the implementation of this script from PYAPP-4797
LID_HBOND_DIST = 2.8
LID_SALT_DIST_CUTOFF = 5.0
LID_METAL_CUTOFF = 2.5

HPHOB_VDW_TOL = 1.2
REC_LIG_CUTOFF_DIST = 5.0
SALT_DIST_CUTOFF = 3.5
METAL_DIST_CUTOFF = 0.5

# NOTE: We are using 2.5A here to match behavior of previous
# version of this script.
HBOND_MAX_DIST_DEFAULT = 2.5
# For other H-bond parameters, default is to use Maestro default values, which
# is accomplished by passing in <None> for these parameters to the hbond API.

XBOND_MAX_DIST = None
XBOND_MIN_DONOR_ANGLE = None
XBOND_MIN_ACCEPTOR_ANGLE = None
XBOND_MAX_ACCEPTOR_ANGLE = None

ABOND_MAX_DIST = None
ABOND_MIN_DONOR_ANGLE = None
ABOND_MIN_ACCEPTOR_ANGLE = None
ABOND_MAX_ACCEPTOR_ANGLE = None

PI_PI_EDGE_TO_FACE_MAX_DIST = None
PI_PI_FACE_TO_FACE_MAX_DIST = None
PI_CAT_MAX_DIST = None

CONTACT_TYPES_ERROR = """
At least one contact type must be specified if the -contacts option
is used. All contacts can be specified using "all" or specific
contacts can be specified using a comma-separated list of the
desired contacts selected from: hbond, xbond, good, bad, ugly, salt,
picat, pipi, hphob, water_hbonds, arom_hbonds, and metal.
For example: -contacts all or -contacts good,bad
"""


def parse_args():

    parser = cmdline.create_argument_parser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(
        "infiles",
        nargs="+",
        help="Input files. If one file is given (i.e. a "
        "pose-viewer or a complex file), interactions are identified "
        "between the first structure and each subsequent "
        "structure. If two files are specified, "
        "interactions are identified between the (first)"
        "structure in the first file and each structure "
        "in the second file..")

    parser.add_argument("-csv",
                        action='store_true',
                        help="Write the output to a CSV file.")

    parser.add_argument('-good_cutoff',
                        action='store',
                        type=float,
                        default=1.3,
                        metavar='<cutoff for good contacts>',
                        help='Maximum ratio of VDW sum allowed for good '
                        'contacts (def. 1.3)')

    parser.add_argument('-bad_cutoff',
                        action='store',
                        type=float,
                        default=0.89,
                        metavar='<cutoff for bad contacts>',
                        help='Maximum ratio of VDW sum allowed for bad '
                        'contacts (def. 0.89)')

    parser.add_argument('-ugly_cutoff',
                        action='store',
                        type=float,
                        default=0.75,
                        metavar='<cutoff for ugly contacts>',
                        help='Maximum ratio of VDW sum allowed for ugly '
                        'contacts (def. 0.75)')

    parser.add_argument(
        '-metal_cutoff',
        type=float,
        default=None,
        metavar='<cutoff for metal contacts>',
        help='Detect metal contacts that are within this distance of each '
        'of each other (+ vdW radii of each atom). (def. 0.5) not '
        'compatible to use with -lid_setting')

    parser.add_argument(
        '-hbond_max_dist',
        type=float,
        default=None,
        help='Maximum H-bond distance. NOTE: Default value is 2.5A, for '
        'backwards compatibility. Default value in Maestro is 2.8A, '
        '-lid-setting is only compatible with 2.8A for this setting. '
        'Does not apply to halogen and aromatic H-bonds.')

    parser.add_argument(
        '-hbond_min_d_angle',
        required=False,
        type=float,
        help='The minimum hydrogen bond donor angle. If not specified, '
        'Maestro preference default value will be used. '
        'Does not apply to halogen and aromatic H-bonds.')

    parser.add_argument(
        '-hbond_min_a_angle',
        required=False,
        type=float,
        help='The minimum hydrogen bond acceptor angle. If not specified, '
        'Maestro preference default value will be used. '
        'Does not apply to halogen and aromatic H-bonds.')

    parser.add_argument(
        '-hbond_max_a_angle',
        required=False,
        type=float,
        help='The maximum hydrogen bond acceptor angle. If not specified, '
        'Maestro preference default value will be used. '
        'Does not apply to halogen and aromatic H-bonds.')

    parser.add_argument(
        '-lid_setting',
        action='store_true',
        default=False,
        help='Use the same settings with the Maestro Ligand Interaction Diagram. '
        'hydrogen bond distance 2.8A, salt bridge distance 5.0A and metal cutoff '
        'of 2.5A regardless of vDW radius')

    parser.add_argument('-contacts',
                        action='store',
                        default='all',
                        metavar='<contacts to calculate>',
                        help='Should be "all" or a comma-separated list of the '
                        'desired contacts selected from: hbond, xbond, good, '
                        'bad, ugly, salt, picat, pipi, and hphob. For '
                        'example: -contacts good,bad (def. all)')

    # TODO: Remove completely in the future:
    parser.add_argument('-residues',
                        action='store_true',
                        help=argparse.SUPPRESS)

    parser.add_argument('-smarts',
                        action='store',
                        default=None,
                        help='Report only interactions to the ligand atoms '
                        'matching this SMARTS pattern (including '
                        'hydrogens bound to them)')

    parser.add_argument(
        '-lig_asl',
        action='store',
        default=None,
        help=
        'ASL pattern for ligand atom types to be considered for interactions.'
        'If input is a complex file, this option can also be used for'
        'selecting which ligand to analyze'
        'example: -lig_asl "ptype CA"; -lig_asl "mol.num 1"; or'
        '-lig_asl "mol.num 1 AND ptype CA"')

    parser.add_argument(
        '-rec_asl',
        action='store',
        default=None,
        help=
        'ASL pattern for receptor atom types to be considered for interactions')

    parser.add_argument('-hbond-heavy',
                        action='store_true',
                        help='For H-bond interactions, report the heavy atoms '
                        'instead of the hydrogens')

    parser.add_argument(
        '-per_ligand',
        action='store_true',
        default=False,
        help=
        'Write out interactions grouped by ligands in a separate output file ')

    parser.add_argument('-opv',
                        action='store',
                        help='Write annotated PoseViewer file to this path. '
                        'Total interaction counts will be added as '
                        'properties to each ligand.')
    return parser


def get_resid(at):
    """
    Return a residue identifier string for an atom.
    Example: "A:100i(ALA)"
    """
    res_str = str(at.getResidue())
    return "%s(%s)" % (res_str, at.pdbres.strip())


def get_residue_key(res_str):
    """
    Return a key for the given residue ID - for sorting.
    """
    res_str = res_str.split('(')[0]
    chain, resnum_ins = res_str.split(':')
    try:
        resnum = int(resnum_ins)
        inscode = ''
    except ValueError:
        resnum = int(resnum_ins[:-1])
        inscode = resnum_ins[-1]

    return (chain, resnum, inscode)


def get_lig_atoms(rec_st, rec_atoms_matching_lig_asl):
    """
    Return list of ligand atom indices that have at least
    one atom matching the given asl, if input is a complex file.

    :return: List of Atom indices
    :rtype: List
    """
    lig_atom_indices = []

    lig_mol_nums = list(set([rec_st.atom[at].molecule_number for at in rec_atoms_matching_lig_asl]))

    for m_idx in lig_mol_nums:
        lig_atom_indices.extend(rec_st.molecule[m_idx].getAtomIndices())
        sr2.append(rec_st.molecule[m_idx].extractStructure())

    return lig_atom_indices

class Contact:
    """
    Data object for storing information about each contact found

    ctype = contact type
    ctindex = structure2 index (ligand number)
    title2 = structure2 title
    atom1 = structure1 (receptor) atom object
    atom2 = structure1 (ligand) atom object
    residue1 = structure1 (receptor) residue label
    residue2 = structure2 (ligand) residue label
    distance = contact distance
    angle = contact angle. For H-bonds, acceptor angle.
    """

    def __init__(self, ctype, at1, at2, distance, angle=None):
        """
        ctype = contact type
        at1 = structure1 atom object
        at2 = structure2 atom object
        distance = distance between the 2 atoms.
        angle = contact angle. For H-bonds, acceptor angle.
        """
        self.ctype = ctype
        self.ctindex = None  # structure2 index (ligand number)
        self.title2 = at2._ct.title
        self.atom1 = at1
        self.atom2 = at2
        self.residue1 = get_resid(at1)
        self.residue2 = get_resid(at2)
        self.distance = distance
        self.angle = angle


def report_contacts(contacts_by_type, hbond_heavy):
    """
    Generate report from the given contacts. To be written to `*.txt` and
    `*.csv` files by calling code.
    """

    header = [
        "%9s" % "Type",
        "%6s" % "Ligand",
        "%-20s" % "Title",
        "%12s" % "LigResidue",
        "%11s" % "LigAtom",
        "%12s" % "RecResidue",
        "%11s" % "RecAtom",
        "%6s" % "Dist",
        "%6s" % "Angle",  # For HBond, PiPi, and PiCat contacts only
    ]
    rows_by_ctype = OrderedDict()

    for ctype, contacts in contacts_by_type.items():
        rows = []
        for contact in contacts:
            # atom1 is from the receptor, atom2 is from the ligand.
            atom1 = contact.atom1
            atom2 = contact.atom2
            type_str = SHORTFORM[ctype]

            # If -hbond-heavy option was used, report heavy atom instead
            # of hydrogen. Also report acceptor/donor sub-type.

            if type_str == "HBond":
                # Report whether it's a donor or an acceptor (from ligand perspective):
                # NOTE: For halogen bonds, we do not currently report this
                if atom2.element == 'H':
                    type_str = "HDonor"
                    atom1_heavy = atom1
                    atom2_heavy = atom2.bond[1].atom2
                    if hbond_heavy:
                        # Report the heavy atom instead:
                        atom2 = atom2_heavy
                elif atom2.element in ('O', 'N', 'S', 'F'):
                    type_str = "HAccep"
                    atom2_heavy = atom2
                    atom1_heavy = atom1.bond[1].atom2
                    if hbond_heavy:
                        # Report the heavy atom instead:
                        atom1 = atom1_heavy
                else:
                    raise ValueError("Invalid element for H-bond: %s" %
                                     atom2.element)

                # PYAPP-6518 Append either "nn", "cn", "nc", or "cc" to the
                # contact type label, depending on the charge of heavy atoms.
                ch1_label = 'n' if atom1_heavy.formal_charge == 0 else 'c'
                ch2_label = 'n' if atom2_heavy.formal_charge == 0 else 'c'
                type_str += ' ' + ch1_label + ch2_label

            row = [
                "%9s" % type_str,
                "%6d" % contact.ctindex,
                "%-20s" % contact.title2,
                "%12s" % contact.residue2,
                "%5s(%4s)" % (atom2.index, atom2.pdbname),
                "%12s" % contact.residue1,
                "%5s(%4s)" % (atom1.index, atom1.pdbname),
                "%6.3f" % contact.distance,
                "%6.1f" %
                contact.angle if contact.angle is not None else '   N/A',
            ]
            rows.append(row)
        rows_by_ctype[ctype] = rows

    return header, rows_by_ctype


def write_contact_txt(contact_types, all_contacts, write_per_residue=False):
    """
    :param contact_types: contact types requested by user
    :param all_contacts: all contacts found
    :param write_per_residue: a flag to write out per-residue contact information
    """
    # Group contacts by type (requested types only):
    contacts_by_type = {}
    for ctype in contact_types:
        contacts_by_type[ctype] = [c for c in all_contacts if c.ctype == ctype]
    header, reports_by_ctype = report_contacts(contacts_by_type,
                                               cmd_args.hbond_heavy)
    header_str = ' '.join(header)
    report_text = ""
    res_report_by_ctype = report_by_residue(contacts_by_type)
    for ctype, rows in reports_by_ctype.items():
        if len(rows) == 0:
            continue
        report_text += "\n========== %s " % ctype + "=" * (70 - len(ctype) + 19)
        report_text += "\n%s\n" % header_str
        for row in rows:
            report_text += ' '.join(row) + '\n'
        if write_per_residue:
            report_text += "\nInteractions grouped by receptor residue:"
            report_text += res_report_by_ctype[ctype]
    return report_text, header, reports_by_ctype


def report_by_residue(contacts_by_type):
    """
    For each contact type, report each protein residue that interacts with
    any ligand, and show all ligands it is interacting with.
    """

    lines_by_ctype = OrderedDict()
    for ctype, contacts in contacts_by_type.items():
        type_str = SHORTFORM[ctype]
        lines = [""]
        lines.append("%8s %-14s %8s %8s" %
                     ("Type", "Residue", "# of ligands", "# of contacts"))

        # Group contacts of this type by residue:
        residue_contacts = defaultdict(list)
        for contact in contacts:
            residue_contacts[contact.residue1].append(contact)

        reslist = sorted(list(residue_contacts), key=get_residue_key)

        total_contacts = 0
        for res_str in reslist:
            contacts = residue_contacts[res_str]
            contacts_by_ligand = defaultdict(list)
            count = len(contacts)
            for contact in contacts:
                contacts_by_ligand[contact.ctindex].append(contact)

            total_contacts += count
            # Header row for each residue lists residue name, number of
            # interacting ligands, and total number of interactions:
            # NOTE: We can't report more specific type of HBond (donor/acceptor)
            # because we are grouping all interactions per ligand together.
            lines.append("%8s %-14s %8d %8s" %
                         (type_str, res_str, len(contacts_by_ligand), count))
            # Then for each residue, list each interacting ligand:
            for lig_index, contacts in sorted(contacts_by_ligand.items()):
                lig_title = contacts[0].title2
                num_inter = len(contacts)
                lines.append("    ligand title: %20s    interactions: %3s" %
                             (lig_title, num_inter))

        lines.append("Total %s interactions: %i" % (ctype, total_contacts))
        lines.append("")
        lines_by_ctype[ctype] = '\n'.join(lines)
    return lines_by_ctype


def get_vdw_radius(atom):
    """
    Finds the vdw radius of atom

    Maestro replaces atom type 18 with atom type 15 when it calculates the vdw
    radius for contacts, so we will do the same with this routine. This is done
    to keep carboxylate atoms equivalent. But this means that we can't use the
    simpler atom.vdw_radius property.

    :type atom: schrodinger.structure._StructureAtom object
    :param atom: the atom to calculate the vdw radius for

    :rtype: float
    :return: the vdw radius of atom
    """
    # FIXME starting with 17-2 release this is no longer necessary.

    mytype = atom.atom_type
    if mytype == 18:
        mytype = 15
    radius = mm.mmat_get_vdw_radius(mytype)
    return radius


def get_smallest_acceptor_angle(hydrogen_atom, acceptor_atom):
    # Ported to Python from structureinteraction.cpp
    smallest_angle = 360.0
    for bonded_to_acceptor in acceptor_atom.bonded_atoms:
        if hydrogen_atom == bonded_to_acceptor:
            # in case of metal interactions ignore the zero-order bond
            # between metal and acceptor
            continue
        bond_angle = measure.measure_bond_angle(hydrogen_atom, acceptor_atom,
                                                bonded_to_acceptor)
        if bond_angle < smallest_angle:
            smallest_angle = bond_angle
    return smallest_angle


def measure_acceptor_angle(atom1, atom2):
    """
    Measure the acceptor angle for the given H-bond. Either atom can be a
    donor (interacting hydrogen) or an acceptor.

    :param atom1: One of the atoms in the interaction.
    :type atom1: structure._StructureAtom

    :param atom2: The other atom in the interaction.
    :type atom2: structure._StructureAtom
    """
    if atom1.element == 'H':
        hydrogen = atom1
        acceptor = atom2
    elif atom2.element == 'H':
        hydrogen = atom2
        acceptor = atom1
    else:
        return None

    assert acceptor.element != 'H'  # sanity check

    return get_smallest_acceptor_angle(hydrogen, acceptor)


def find_matching_lig_atoms(lig_st, smarts, lig_asl=None):
    """
    Return a list of atoms matching the given SMARTS or ASL pattern, including
    any hydrogens bound to them. If SMARTS or ASL is None, returns a list of all
    atoms in the structure. If both SMARTS and ASL is given, returns no match,
    only one can be provided.
    """
    matching_lig_atoms = set()

    if not smarts and not lig_asl:
        return lig_st.getAtomIndices()

    elif smarts:
        for aset in adapter.evaluate_smarts(lig_st, smarts):
            matching_lig_atoms.update(aset)

    elif lig_asl:
        matching_lig_atoms = set(analyze.evaluate_asl(lig_st, lig_asl))

    if matching_lig_atoms:
        # Add any hydrogens bound to those atoms:
        for a in lig_st.atom:
            if a.element == "H" and a.index not in matching_lig_atoms:
                if a.bond[1].atom2.index in matching_lig_atoms:
                    matching_lig_atoms.add(a.index)
    return sorted(matching_lig_atoms)


def assemble_matching_pairs(rec_dist_cell,
                            lig_atoms,
                            rec_st,
                            lig_st,
                            rec_atoms=None):
    """
    :param rec_dist_cell: iterator that uses a distance cell to iterate through
        neighbors of the specified atoms
    :type rec_dist_cell: iter

    :param lig_atoms: List of ligand atom indices matching the ligand ASL
    :type lig_atoms: list

    :param rec_st: Receptor structure object
    :type rec_st: `Structure<structure.Structure>`

    :param lig_st: Ligand structure object
    :type lig_st: `Structure<structure.Structure>`

    :param rec_atoms: List of receptor atom indices matching the ligand ASL
    :type rec_atoms: list
    """

    matched_atom_pairs = []

    for lig_anum, close_rec_atoms in rec_dist_cell.iterateNeighboringAtoms(
            lig_st, lig_atoms):
        lig_atom = lig_st.atom[lig_anum]

        if rec_atoms:
            close_rec_atoms = set(close_rec_atoms) & set(rec_atoms)

        for rec_anum in sorted(close_rec_atoms):
            rec_atom = rec_st.atom[rec_anum]
            dist = measure.measure_distance(lig_atom, rec_atom)
            matched_atom_pairs.append((lig_atom, rec_atom, dist))

    return matched_atom_pairs


def categorize_contacts(contact_type, bonds, matched_atom_pairs):
    """
    Returns a subset of matched_atom_pairs that match the bonds
    list - where each item is a (donor atom object, acceptor atom object).
    """
    # Set of {rec atom object, lig atom object} sets.
    contacts = []
    bond_set = set(map(frozenset, bonds))
    for lig_at, rec_at, dist in matched_atom_pairs:
        if {rec_at, lig_at} in bond_set:
            contacts.append(Contact(contact_type, rec_at, lig_at, dist))

    return contacts


def find_hbond_contacts(rec_st, lig_st, matched_atom_pairs, hb_params):
    """
    Find HBond contacts between the receptor and the
    given subset of the ligand atoms.

    :return: Generator of contacts, each being type `Contact`
    :rtype: generator
    """
    hbonds = hbond_module.get_hydrogen_bonds(st=lig_st, st2=rec_st, **hb_params)
    return categorize_contacts(HBONDS, hbonds, matched_atom_pairs)


def find_xbond_contacts(rec_st, lig_st, matched_atom_pairs):
    """
    Find XBond contacts between the receptor and the
    given subset of the ligand atoms.

    :return: Generator of contacts, each being type `Contact`
    :rtype: generator
    """

    kwargs = {
        'max_dist': XBOND_MAX_DIST,
        'min_donor_angle': XBOND_MIN_DONOR_ANGLE,
        'min_acceptor_angle': XBOND_MIN_ACCEPTOR_ANGLE,
        'max_acceptor_angle': XBOND_MAX_ACCEPTOR_ANGLE,
    }
    xbonds = hbond_module.get_halogen_bonds(st=lig_st, st2=rec_st, **kwargs)

    return categorize_contacts(XBONDS, xbonds, matched_atom_pairs)


def find_aromatic_hbond_contacts(rec_st, lig_st, matched_atom_pairs):
    """
    Find Aromatic HBond contacts between the receptor and the
    given subset of the ligand atoms.

    :return: Generator of contacts, each being type `Contact`
    :rtype: generator
    """
    kwargs = {
        'max_dist': ABOND_MAX_DIST,
        'min_donor_angle': ABOND_MIN_DONOR_ANGLE,
        'min_acceptor_angle': ABOND_MIN_ACCEPTOR_ANGLE,
        'max_acceptor_angle': ABOND_MAX_ACCEPTOR_ANGLE,
    }
    abonds = hbond_module.get_aromatic_hydrogen_bonds(st=lig_st,
                                                      st2=rec_st,
                                                      **kwargs)
    return categorize_contacts(AROMATIC_HBONDS, abonds, matched_atom_pairs)


def find_water_mediated_hbond_contacts(rec_st, lig_st, hb_params):
    """
    Find receptor-ligand atom pairs which form H-bonds to the same water
    molecule.
    """

    num_rec_atoms = len(rec_st.atom)
    st = rec_st.merge(lig_st)
    all_water_atoms = analyze.evaluate_asl(st, 'water')
    if not all_water_atoms:
        return
    bonds_to_water = defaultdict(list)
    for atom1, atom2 in hbond_module.get_hydrogen_bonds(st, all_water_atoms,
                                                        **hb_params):
        if atom1.index in all_water_atoms:
            bonds_to_water[atom1.molecule_number].append(atom2)
        elif atom2.index in all_water_atoms:
            bonds_to_water[atom2.molecule_number].append(atom1)

    for wat_molnum, bound_atoms in bonds_to_water.items():
        for atom1, atom2 in itertools.combinations(bound_atoms, 2):
            if atom1.index <= num_rec_atoms and atom2.index > num_rec_atoms:
                rec_atom = atom1
                lig_atom = atom2
            elif atom2.index <= num_rec_atoms and atom1.index > num_rec_atoms:
                rec_atom = atom2
                lig_atom = atom1
            else:
                # Otherwise either both atoms are in ligand or in receptor
                continue
            lig_atom = lig_st.atom[lig_atom.index - num_rec_atoms]
            dist = measure.measure_distance(rec_atom, lig_atom)
            yield Contact(WATER_HBONDS, rec_atom, lig_atom, dist)


def categorize_contact_by_dist(at2, at1, dist):
    sum_vdw = get_vdw_radius(at2) + get_vdw_radius(at1)
    if dist <= cmd_args.ugly_cutoff * sum_vdw:
        return UGLY
    elif dist <= cmd_args.bad_cutoff * sum_vdw:
        return BAD
    elif dist < cmd_args.good_cutoff * sum_vdw:
        return GOOD
    return None


def find_salt_contacts(rec_st, lig_st, matched_atom_pairs, salt_dist_cutoff):
    """
    Updated in PYAPP-8188 to use the get_salt_bridge from the infrastructure

    :return: Generator of contacts, each being type `Contact`
    :rtype: generator
    """
    salt_bridges = get_salt_bridges(struc1=lig_st,
                                    struc2=rec_st,
                                    cutoff=salt_dist_cutoff)
    return categorize_contacts(SALT, salt_bridges, matched_atom_pairs)


def find_hphob_contacts(lig_st, matched_atom_pairs):
    """
    Find HPhob contacts between the receptor and the
    given subset of the ligand atoms.

    :return: Generator of contacts, each being type `Contact`
    :rtype: generator
    """

    hphobic_lig_atoms = analyze.evaluate_asl(lig_st, HPHOBIC_ASL)

    for lig_at, rec_at, dist in matched_atom_pairs:
        sum_vdw = get_vdw_radius(lig_at) + get_vdw_radius(rec_at)
        if rec_at.index in hphobic_rec_atoms and \
                        lig_at.index in hphobic_lig_atoms:
            if dist <= HPHOB_VDW_TOL * sum_vdw:
                yield Contact(HPHOB, rec_at, lig_at, dist)


def find_metal_contacts(rec_st, matched_atom_pairs, metal_cutoff,
                        metal_include_vdw):

    metal_atoms = analyze.evaluate_asl(rec_st, 'metals')

    for lig_atom, rec_atom, dist in matched_atom_pairs:
        if rec_atom.index in metal_atoms:
            # PYAPP-8188, LID uses only distance as the criterion for metal contacts
            if metal_include_vdw:
                final_metal_cutoff = get_vdw_radius(lig_atom) + get_vdw_radius(
                    rec_atom) + metal_cutoff
            else:
                final_metal_cutoff = LID_METAL_CUTOFF
            if dist <= final_metal_cutoff:
                yield Contact(METAL, rec_atom, lig_atom, dist)


def find_picat_contacts(rec_st, lig_st, rec_atoms, lig_atoms):
    """
    Find PICAT contacts between the receptor and the given subset
    of the ligand atoms.

    :param rec_atoms: List of receptor atom indices matching the ligand ASL
    :type rec_atoms: list

    :param lig_atoms: List of ligand atom indices matching the ligand ASL
    :type lig_atoms: list

    :return: Generator of contacts, each being type `Contact`
    :rtype: generator
    """

    params = interactions.infrastructure.PiCationParams()
    if PI_CAT_MAX_DIST is not None:
        params.setMaximumDistance(PI_CAT_MAX_DIST)

    picats = interactions.find_pi_cation_interactions(rec_st, lig_st)
    for picat in picats:
        pi_atom = picat.pi_structure.atom[picat.pi_centroid.atoms[0]]
        cat_atom = picat.cation_structure.atom[picat.cation_centroid.atoms[0]]
        if picat.pi_structure == lig_st:
            lig_atom = pi_atom
            rec_atom = cat_atom
        else:
            assert picat.pi_structure == rec_st
            lig_atom = cat_atom
            rec_atom = pi_atom

        if (lig_atom.index in lig_atoms) and rec_atom.index in rec_atoms:
            yield Contact(PICAT, rec_atom, lig_atom, picat.distance,
                          picat.angle)


def find_pipi_contacts(rec_st, lig_st, matching_rec_atoms, matching_lig_atoms):
    """
    Find PIPI_FACE and PIPI_EDGE contacts between the receptor and
    the given subset of the ligand atoms.

    :return: Generator of contacts, of type `Contact`
    :rtype: generator
    """

    params = interactions.infrastructure.PiPiParams()
    if PI_PI_EDGE_TO_FACE_MAX_DIST is not None:
        params.setEdgeToFaceMaximumDistance(PI_PI_EDGE_TO_FACE_MAX_DIST)
    if PI_PI_FACE_TO_FACE_MAX_DIST is not None:
        params.getFaceToFaceMaximumAngle(PI_PI_FACE_TO_FACE_MAX_DIST)

    pipis = interactions.find_pi_pi_interactions(rec_st, lig_st, params=params)
    for pipi in pipis:
        if pipi.face_to_face:
            pipitype = PIPI_FACE
        else:
            pipitype = PIPI_EDGE

        ring1a = pipi.struct1.atom[pipi.ring1.atoms[0]]
        ring2a = pipi.struct2.atom[pipi.ring2.atoms[0]]
        # In practice struct1 is always rec_st and struct2 is lig_st;
        # but verify them:
        assert ring1a._ct == rec_st and ring2a._ct == lig_st
        if (ring2a.index
                in matching_lig_atoms) and ring1a.index in matching_rec_atoms:
            yield Contact(pipitype, ring1a, ring2a, pipi.distance, pipi.angle)


def find_ligand_contacts(lig_index, lig_st, rec_st, rec_dist_cell, opts):
    """
    Returns a list of contacts for all interactions between the receptor
    and the ligand.

    :return: List of all contacts between the ligand and receptor.
    :rtype: list of `Contact`
    """

    # Make a list of the ligand atoms we are interested in:
    lig_atoms = find_matching_lig_atoms(lig_st, cmd_args.smarts,
                                        cmd_args.lig_asl)

    if not lig_atoms:
        logger.warning(
            "WARNING: None of the ligand atoms matched the ASL or SMARTS pattern provided, Skipping ligand %i"
            % lig_index)
        return

    if cmd_args.rec_asl:
        rec_atoms = analyze.evaluate_asl(rec_st, cmd_args.rec_asl)
        if not rec_atoms:
            logger.warning(
                "WARNING: Receptor ASL expression matches no receptor atoms!")
            return
    else:
        rec_atoms = [rec_at.index for rec_at in rec_st.atom]

    # Get a list of receptor atoms within 5A of the matching ligand atoms:
    # - list of (lig_atom_obj, rec_atom_obj, distance)
    matched_atom_pairs = assemble_matching_pairs(rec_dist_cell, lig_atoms,
                                                 rec_st, lig_st, rec_atoms)

    if not matched_atom_pairs:
        logger.warning(
            "WARNING: None of the receptor atoms matching the ASL were close to the ligand; skipping ligand %i"
            % lig_index)
        return

    # Custom parameters for H-bonds.  Does not apply for halogen h-bonds and
    # aromatic-bonds.
    hb_params = {
        'max_dist': opts.hbond_max_dist,
        'min_donor_angle': opts.hbond_min_d_angle,
        'min_acceptor_angle': opts.hbond_min_a_angle,
        'max_acceptor_angle': opts.hbond_max_a_angle,
    }

    contacts = []
    if HBONDS in contact_types:
        contacts += find_hbond_contacts(rec_st, lig_st, matched_atom_pairs,
                                        hb_params)

    if XBONDS in contact_types:
        contacts += find_xbond_contacts(rec_st, lig_st, matched_atom_pairs)

    if WATER_HBONDS in contact_types:
        contacts += find_water_mediated_hbond_contacts(rec_st, lig_st,
                                                       hb_params)

    if AROMATIC_HBONDS in contact_types:
        contacts += find_aromatic_hbond_contacts(rec_st, lig_st,
                                                 matched_atom_pairs)

    if GOOD in contact_types or BAD in contact_types or UGLY in contact_types:
        for lig_at, rec_at, dist in matched_atom_pairs:
            ctype = categorize_contact_by_dist(lig_at, rec_at, dist)
            if ctype and ctype in contact_types:
                contacts.append(Contact(ctype, rec_at, lig_at, dist))

    if SALT in contact_types:
        contacts += find_salt_contacts(rec_st, lig_st, matched_atom_pairs,
                                       opts.salt_dist_cutoff)

    if HPHOB in contact_types:
        contacts += find_hphob_contacts(lig_st, matched_atom_pairs)

    if PICAT in contact_types:
        contacts += find_picat_contacts(rec_st, lig_st, rec_atoms, lig_atoms)

    if PIPI_FACE in contact_types or PIPI_EDGE in contact_types:
        for contact in find_pipi_contacts(rec_st, lig_st, rec_atoms, lig_atoms):
            if contact.ctype in contact_types:
                contacts.append(contact)

    if METAL in contact_types:
        contacts += find_metal_contacts(rec_st, matched_atom_pairs,
                                        opts.metal_cutoff,
                                        opts.metal_include_vdw)

    for contact in contacts:
        contact.ctindex = lig_index
        # For Hydrogen bonds, calculate acceptor angle:
        if contact.ctype == HBONDS:
            # Currently halogen bond angles are not supported
            contact.angle = measure_acceptor_angle(contact.atom1, contact.atom2)

    return contacts


if __name__ == '__main__':

    parser = parse_args()
    cmd_args = parser.parse_args()

    if not len(cmd_args.infiles) or len(cmd_args.infiles) > 2:
        parser.error(
            "Either a single PV, complex file or two structures files are "
            "required.")
        sys.exit(1)

    if len(cmd_args.infiles) == 1:
        file1 = cmd_args.infiles[0]
        if fileutils.is_poseviewer_file(file1):
            file_type = "pv_file"
        else:
            file_type = "complex_file"
    else:
        file_type = "two_files"
        file1 = cmd_args.infiles[0]
        file2 = cmd_args.infiles[1]

    if not os.path.isfile(file1):
        parser.error("Error: The file '%s' does not exist." % file1)

    if not file_type and not os.path.isfile(file2):
        parser.error("Error: The file '%s' does not exist." % file2)

    if cmd_args.lid_setting:
        if cmd_args.hbond_max_dist is not None and cmd_args.hbond_max_dist != LID_HBOND_DIST:
            parser.error(
                "Error: The specified hbond cutoff is not consistent with LID.")
        if cmd_args.metal_cutoff is not None:
            # LID uses only distance (2.5A) as the cutoff regardless of vDW radius
            parser.error("Error: Cannot specify metal cutoff with lid_setting.")
        cmd_args.hbond_max_dist = LID_HBOND_DIST
        cmd_args.salt_dist_cutoff = LID_SALT_DIST_CUTOFF
        cmd_args.metal_cutoff = LID_METAL_CUTOFF
        cmd_args.metal_include_vdw = False
    else:
        if cmd_args.hbond_max_dist is None:
            cmd_args.hbond_max_dist = HBOND_MAX_DIST_DEFAULT
        cmd_args.salt_dist_cutoff = SALT_DIST_CUTOFF
        if cmd_args.metal_cutoff is None:
            cmd_args.metal_cutoff = METAL_DIST_CUTOFF
        cmd_args.metal_include_vdw = True

    if cmd_args.smarts and cmd_args.lig_asl:
        logger.warning(
            "WARNING: Provide either SMARTS or ASL pattern for matching ligand atoms!"
        )
        sys.exit(1)

    logger.info("Initializing")

    requested_contacts = cmd_args.contacts.lower().split(',')
    contact_types = []

    if 'all' in requested_contacts:
        contact_types = list(OPTION_VALUES)
    else:
        contact_types = [
            key for key, val in OPTION_VALUES.items()
            if val in requested_contacts
        ]
        if not contact_types:
            parser.error(CONTACT_TYPES_ERROR)

    sr1 = structure.StructureReader(file1)
    rec_st = next(sr1)

    if cmd_args.lig_asl:
        rec_atoms_matching_lig_asl = analyze.evaluate_asl(
            rec_st, cmd_args.lig_asl)
    else:
        rec_atoms_matching_lig_asl = None

    if file_type == "pv_file":
        # Read the ligands from the PV file:
        sr2 = sr1
    elif file_type == "complex_file":
        sr2 = []
        # If -lig_asl is specified, pull out only ligands that have at
        # least one atom matching that ASL.  Leave the rest in the receptor:
        ligand_atom_indices = get_lig_atoms(rec_st, rec_atoms_matching_lig_asl)

        rec_st.deleteAtoms(ligand_atom_indices)
        # rec_st.deleteAtoms(rec_atoms_matching_lig_asl)
    else:
        # Read the ligands from the second structure file:
        sr2 = structure.StructureReader(file2)

    hphobic_rec_atoms = None
    if HPHOB in contact_types:
        hphobic_rec_atoms = analyze.evaluate_asl(rec_st, HPHOBIC_ASL)

    rec_dist_cell = measure.DistanceCellIterator(rec_st, REC_LIG_CUTOFF_DIST)

    if cmd_args.opv:
        out_pv_writer = structure.StructureWriter(cmd_args.opv)
        out_pv_writer.append(rec_st)

    contacts_by_ligand = []
    all_contacts = []

    for index, lig_st in enumerate(sr2, start=1):
        print(index)
        logger.info("Processing ligand %i" % index)
        lig_contacts = find_ligand_contacts(index, lig_st, rec_st,
                                            rec_dist_cell, cmd_args)
        all_contacts.extend(lig_contacts)
        contacts_by_ligand.append(lig_contacts)
        if cmd_args.opv:
            # PYAPP-6517 Add count properties, and write to annotated PV file:
            for ctype in contact_types:
                count = len([c for c in lig_contacts if c.ctype == ctype])
                contact_prop = CONTACT_PROPS[ctype]
                lig_st.property[contact_prop] = count
            out_pv_writer.append(lig_st)

    if cmd_args.per_ligand:
        # PYAPP-8160, write output grouped by ligands
        out_text_file = fileutils.get_basename(
            file1) + "_pv_perligand_interactions.txt"
        with open(out_text_file, 'w') as fh:
            for index, lig_contacts in enumerate(contacts_by_ligand):
                lig_report_header = "Ligand # %s\n\n" % str(index + 1)
                fh.write(lig_report_header)
                lig_report_text, _, _ = write_contact_txt(
                    contact_types, lig_contacts, False)
                fh.write(lig_report_text)
        logger.info('Output per-ligand interaction file: %s' % out_text_file)

    out_text_file = fileutils.get_basename(file1) + "_pv_interactions.txt"
    with open(out_text_file, 'w') as fh:
        report_text_header = "Number of ligands: %s\n\n" % len(
            contacts_by_ligand)
        report_text, header, reports_by_ctype = write_contact_txt(
            contact_types, all_contacts, True)
        fh.write(report_text_header)
        fh.write(report_text)

    if cmd_args.csv:
        # If -csv option was used, save the main report to CSV file too:
        csv_file = fileutils.get_basename(file1) + "_pv_interactions.csv"
        with csv_unicode.writer_open(csv_file) as csv_fh:
            writer = csv.writer(csv_fh)
            writer.writerow((list(map(str.strip, header))))
            for _, rows in reports_by_ctype.items():
                for row in rows:
                    row = list(map(str.strip, row))
                    writer.writerow(row)
        logger.info('Output CSV file: %s' % csv_file)

    logger.info('Output text file: %s' % out_text_file)

    if cmd_args.opv:
        out_pv_writer.close()
        logger.info("Annotated PV file: %s" % cmd_args.opv)

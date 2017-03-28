# ModeMap.py by J.M. Skelton


import argparse;
import csv;
import math;
import os;
import tarfile;

import numpy as np;

# This script makes use of a number of routines from the Phonopy Python API.

from phonopy import Phonopy;
from phonopy.file_IO import parse_FORCE_SETS, parse_FORCE_CONSTANTS;
from phonopy.interface.vasp import read_vasp;


# Parse command-line arguments.

parser = argparse.ArgumentParser(description = "Use the Phonopy Python API to map phonon mode eigenvectors as a function of normal-mode coordinate, Q");

# Defailts: read the structure from "POSCAR", read "FORCE_SETS" and calculate the force constants, and generate a 1D map.

parser.set_defaults(
    CellFile = "POSCAR",
    ReadForceConstants = False,
    MapMode = '1D',
    ScaleQ = True
    );

# The following four options mirror relevant command-line arguments that can be passed to Phonopy.

group1 = parser.add_argument_group("Phonopy settings");

group1.add_argument(
    "-c", "--cell",
    metavar = "cell_file",
    type = str, dest = 'CellFile',
    help = "POSCAR file to read (default: POSCAR)"
    );

group1.add_argument(
    "--dim",
    metavar = "'x y z'",
    type = str, dest = 'SupercellMatrix',
    required = True,
    help = "Supercell used to calculated the forces (same as Phonopy --dim option)"
    );

group1.add_argument(
    "--pa", "--primitive_axis",
    metavar = "'xx xy xz yx yy yz zx zy zz'",
    type = str, dest = 'PrimitiveMatrix',
    help = "Transformation matrix to be used in the calculations (same as Phonopy --pa/--primitive_axis option)"
    );

group1.add_argument(
    "--readfc",
    dest = 'ReadForceConstants',
    action = 'store_true',
    help = "Read force constants from a FORCE_CONSTANTS file (same as Phonopy --readfc option)"
    );

# The following variables control the map mode (1D/2D), the mode(s) and amplitude range(s) to be mapped, and the size of the supercell in which to generate the modulated structures.

group2 = parser.add_argument_group("Map settings");

group2.add_argument(
    "--map_2d",
    dest = 'MapMode',
    action = 'store_const', const = '2D',
    help = "Map two modes as a function of Q, specified via the --mode_1 and --mode_2 arguments"
    );

group2.add_argument(
    "--mode", "--mode_1",
    metavar = "'q_x q_y q_z band_index'",
    type = str, dest = 'Mode1',
    required = True,
    help = "Mode to map (first mode if --map_2d is specified), in the form of \"q_x q_y q_z band_index\""
    );

group2.add_argument(
    "--mode_2",
    metavar = "'q_x q_y q_z band_index'",
    type = str, dest = 'Mode2',
    help = "Second mode to map for --map_2d"
    );

group2.add_argument(
    "--q_range", "--q_range_1",
    metavar = "'start end step'",
    type = str, dest = 'QRange1',
    required = True,
    help = "Range of mode coordinates to map, in the form of \"start, stop, step\""
    );

group2.add_argument(
    "--q_range_2",
    metavar = "'start end step'",
    type = str, dest = 'QRange2',
    help = "Optional different range of mode coordinates to map the second mode when using --map_2d"
    );

group2.add_argument(
    "--supercell",
    metavar = "'x y z'",
    type = str, dest = 'ModulationSupercellMatrix',
    required = True,
    help = "Supercell expansion in which to generate modulated structures."
    );

group2.add_argument(
    "--no_q_scale",
    action = 'store_true',
    help = "Do not scale normal-mode coordinates fed to Phonopy by sqrt(N_a) (this is required for correct dU(Q) curves; only turn this off if you know what you're doing...!)"
    );

args = parser.parse_args();

# Convert arguments and perform some basic validation.

elements = [int(element) for element in args.SupercellMatrix.strip().split()];

if len(elements) == 3:
    args.SupercellMatrix = np.diag(elements);
elif len(elements) == 9:
    args.SupercellMatrix = np.array(elements).reshape(3, 3);
else:
    raise Exception("Error: The supercell matrix passed via --dim should be a set of three of nine integers (see corresponding Phonopy option).");

if args.PrimitiveMatrix != None:
    elements = [float(element) for element in args.PrimitiveMatrix.strip().split()];

    if len(elements) != 9:
        raise Exception("Error: The primitive matrix passed via --pa/--primitive_axis, if supplied, must be a set of nine numbers (see correspondiong Phonopy option).");

    args.PrimitiveMatrix = np.array(elements).reshape(3, 3);

if args.MapMode == '2D' and args.Mode2 == None:
    raise Exception("Error: If --map_2d is specified, a second mode to follow must be set via --mode_2.");

modes = [args.Mode1, args.Mode2];

for i, mode in enumerate(modes):
    if mode != None:
        elements = [element for element in mode.strip().split()];

        if len(elements) != 4:
            raise Exception("Error: Modes to follow should be specified by four numbers (q_z, q_y, q_z, bandIndex); supplying a phase factor is not currently supported.");

        qx, qy, qz = [float(element) for element in elements[:3]];
        bandIndex = int(elements[3]);

        if bandIndex <= 0:
            raise Exception("Error: The band indices of modes to follow must be greater than zero.");

        modes[i] = ((qx, qy, qz), bandIndex);

args.Mode1, args.Mode2 = modes;

qRanges = [args.QRange1, args.QRange2];

for i, qRange in enumerate(qRanges):
    if qRange != None:
        elements = [float(element) for element in qRange.strip().split()];

        if len(elements) != 3:
            raise Exception("Error: Mode-coordinate ranges should be specified as three numbers (start, stop, step).");

        start, stop, step = elements;
        qRanges[i] = np.arange(start, stop + step, step);

args.QRange1, args.QRange2 = qRanges;

if args.MapMode == '2D' and args.QRange2 == None:
    args.QRange2 = args.QRange1;

elements = [int(element) for element in args.ModulationSupercellMatrix.strip().split()];

if len(elements) == 3 or len(elements) == 9:
    args.ModulationSupercellMatrix = np.array(elements);
else:
    raise Exception("Error: The supercell matrix passed via --supercell should be a set of three of nine integers.");

# Read the structure file.

structure = read_vasp(args.CellFile);

# Set up a Phonopy object to do the "heavy lifting" of generating the modulated structures.

phonon = Phonopy(
    structure,
    args.SupercellMatrix,
    primitive_matrix = args.PrimitiveMatrix
    );

# Set up the force constants.

if args.ReadForceConstants:
    # If the --readfc flag was set, read the force constants directly from a FORCE_CONSTANTS file.

    phonon.set_force_constants(
        parse_FORCE_CONSTANTS()
        );
else:
    # If not, read the FORCE_SETS file, pass the dataset to the Phonopy object, and calculate the force constants.

    phonon.set_displacement_dataset(
        parse_FORCE_SETS()
        );

    phonon.produce_force_constants();

# Set a scaling factor for the user-input normal-mode coordinate ranges.
# The conversion from normal-mode coordinate amplitude to cartesian displacements in the Phonopy modulation routine uses a definition that leads to fitted harmonic frequencies (for modes with harmonic potential-energy surfaces) being a factor of \sqrt(N_a) too small.
# By default, we therefore apply a scaling factor before passing the amplitudes to the Phonopy routines (unless overridden by the --no__q_scale argument).

qScale = None;

if args.ScaleQ:
    qScale = math.sqrt(len(structure.get_scaled_positions()));

if args.MapMode == '1D':
    # 1D mapping mode.

    qPoint, index = args.Mode1;
    index = index - 1;

    modulations = [];

    fileCounter = 1;

    # For detailed maps, a large number of modulated structures will be generated.
    # Generally, what we will want to do is to transfer the structures to a workstation/HPC service to perform the force calculations.
    # We therefore take advantage of the tarfile module to bundle the structures into a .tar.gz archive, rather than leaving them "loose" in the working directory.

    with tarfile.open(r"ModeMap.tar.gz", 'w:gz') as modulationFilesArchive:
        # For each modulation, use the Phonopy object to generate the modulated structures.
        # The modulated structure is saved to a file named "MPOSCAR", which is added to the archive with a sequentially-numbered file name.

        for q in args.QRange1:
            phonon.set_modulations(
                args.ModulationSupercellMatrix,

                # Scale the normal-mode coordinate before passing to Phonopy, if required.

                [[qPoint, index, q * qScale if qScale != None else q, 0.0]]
                );

            phonon.write_modulations();

            modulationFilesArchive.add("MPOSCAR", arcname = "MPOSCAR-{0:0>3}".format(fileCounter));

            modulations.append(q);

            fileCounter = fileCounter + 1;

    # The details of the modulated structures - in particular the mapping of file numbers to modulation amplitudes - is stored in a CSV-format file.
    # This file is used by ModeMap_PostProcess.py.

    with open("ModeMap.csv", 'w') as outputWriter:
        outputWriterCSV = csv.writer(outputWriter, delimiter = ',', quotechar = '\"', quoting = csv.QUOTE_ALL);

        qx, qy, qz = qPoint;

        outputWriterCSV.writerow(["q = ({0:.3f}, {1:.3f}, {2:.3f}), band = {3}".format(qx, qy, qz, index + 1)]);

        outputWriterCSV.writerow([]);

        outputWriterCSV.writerow(["Modulation #", "Q [amu^1/2 A"]);

        for i, dispAmplitude in enumerate(modulations):
            outputWriterCSV.writerow([i + 1, dispAmplitude]);

        outputWriter.close();

elif args.MapMode == '2D':
    # 2D mapping mode.

    qPoint1, index1 = args.Mode1;
    qPoint2, index2 = args.Mode2;

    index1, index2 = index1 - 1, index2 - 1;

    modulations = [];

    fileCounter = 1;

    # As in the 1D mapping mode, the modulated structures are bundled into a .tar.gz file for convenience.

    with tarfile.open(r"ModeMap.tar.gz", 'w:gz') as modulationFilesArchive:
        for q1 in args.QRange1:
            for q2 in args.QRange2:
                phonon.set_modulations(
                    args.ModulationSupercellMatrix,
                    [
                        # Again if needed, scale the normal-mode coordinates before passing to Phonopy.

                        [qPoint1, index1, q1 * qScale if qScale != None else q1, 0.0],
                        [qPoint2, index2, q2 * qScale if qScale != None else q2, 0.0]]
                    );

                phonon.write_modulations();

                modulationFilesArchive.add("MPOSCAR", arcname = "MPOSCAR-{0:0>3}".format(fileCounter));

                modulations.append((q1, q2));

                fileCounter = fileCounter + 1;

    # Again, the details of the modulated structures are written to a CSV-format file.

    with open("ModeMap.csv", 'w') as outputWriter:
        outputWriterCSV = csv.writer(outputWriter, delimiter = ',', quotechar = '\"', quoting = csv.QUOTE_ALL);

        qx1, qy1, qz1 = qPoint1;
        qx2, qy2, qz2 = qPoint2;

        outputWriterCSV.writerow(["Mode 1", "q = ({0:.3f}, {1:.3f}, {2:.3f}), band = {3}".format(qx1, qy1, qz1, index1 + 1)]);
        outputWriterCSV.writerow(["Mode 2", "q = ({0:.3f}, {1:.3f}, {2:.3f}), band = {3}".format(qx2, qy2, qz2, index2 + 1)]);

        outputWriterCSV.writerow([]);

        outputWriterCSV.writerow(["Modulation #", "Q_1 [amu^1/2 A]", "Q_2 [amu^1/2 A]"]);

        for i, (dispAmplitude1, dispAmplitude2) in enumerate(modulations):
            outputWriterCSV.writerow([i + 1, dispAmplitude1, dispAmplitude2]);

# Depending on the modulation parameters, the Phonopy routines write a number of files to the working directory.
# After generating the required modulated structures, any files present should be cleaned up.

for tempFile in "MPOSCAR", "MPOSCAR-001", "MPOSCAR-002", "MPOSCAR-orig":
    if os.path.isfile(tempFile):
        os.remove(tempFile);

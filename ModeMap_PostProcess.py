# ModeMap_PostProcess.py by J. M. Skelton


import argparse;
import csv;

import numpy as np;

from scipy.interpolate import interp1d, interp2d;


# Function to read the "ModeMap.csv" files produced by ModeMap.py.

def _ReadModeMapCSV(filePath = "ModeMap.csv", mapMode = '1D'):
    modulationModeCoordinates = [];

    skipRows = None;

    if mapMode == '1D':
        skipRows = 3;
    elif mapMode == '2D':
        skipRows = 4;
    else:
        raise Exception("Error: Unknown mapMode '{0}'.".format(mapMode));

    with open(filePath, 'r') as inputReader:
        inputReaderCSV = csv.reader(inputReader);

        for i in range(0, skipRows):
            next(inputReaderCSV);

        for row in inputReaderCSV:
            modulationModeCoordinates.append(
                tuple([int(row[0])] + [float(item) for item in row[1:]])
                );

    itemLength = len(modulationModeCoordinates[0]);

    for item in modulationModeCoordinates[1:]:
        if len(item) != itemLength:
            raise Exception("Error: Inconsistent row length in mode map CSV file \"{0}\".".format(filePath));

    return sorted(modulationModeCoordinates, key = lambda item : item[0]);

# Function to read the "ExtractTotalEnergies.csv" files produced by the ExtractTotalEnergies.py utility script.

def _ReadExtractTotalEnergiesCSV(filePath = "ExtractTotalEnergies.csv", headerRows = 1):
    modulationTotalEnergies = [];

    with open("ExtractTotalEnergies.csv", 'r') as inputReader:
        inputReaderCSV = csv.reader(inputReader);

        for i in range(0, headerRows):
            next(inputReaderCSV);

        for modulation, totalEnergy in inputReaderCSV:
            modulationTotalEnergies.append(
                (int(modulation), float(totalEnergy))
                );

    return sorted(modulationTotalEnergies, key = lambda item : item[0]);


# Parse command-line arguments.

parser = argparse.ArgumentParser(description = "Post process calculations set up using ModeMap.py");

# By default, the script expects to process a 1D map.

parser.set_defaults(
    MapMode = '1D'
    );

parser.add_argument(
    "--map_2d",
    dest = 'MapMode',
    action = 'store_const', const = '2D',
    help = "Post process data from a ModeMap.py calculation with the --map_2d option."
    );

args = parser.parse_args();

# Read the input data from "ModeMap.csv" and "ExtractTotalEnergies.csv" and perform some basic validation.

modulationModeCoordinates = _ReadModeMapCSV(filePath = "ModeMap.csv", mapMode = args.MapMode);

modulationTotalEnergies = _ReadExtractTotalEnergiesCSV(filePath = "ExtractTotalEnergies.csv");

if len(modulationModeCoordinates) != len(modulationTotalEnergies):
    raise Exception("Error: Number of entries in \"ModeMap.csv\" and \"ExtractTotalEnergies.csv\" do not match.");

for mod1, (mod2, _) in zip([item[0] for item in modulationModeCoordinates], modulationTotalEnergies):
    if mod1 != mod2:
        raise Exception("Error: Modulation numbers in \"ModeMap.csv\" and \"ExtractTotalEnergies.csv\" are inconsistent.");

# Post process.

if args.MapMode == '1D':
    # For 1D maps, join the data in the "ModeMap.csv" and "ExtractTotalEnergies.csv" files to output energy as a function of normal-mode coordinate.

    # To output relative energies (\DeltaU(Q)), we need an energy zero.

    uZero = None;

    # First, check whether q = 0 is included in the mapping points.

    for i, (_, q) in enumerate(modulationModeCoordinates):
        if q == 0.0:
            _, uZero = modulationTotalEnergies[i];
            break;

    # If not, estimate it using 1D cubic interpolation.

    if uZero == None:
        interpFunction = interp1d([q for _, q in modulationModeCoordinates], [u for _, u in modulationTotalEnergies], kind = 'cubic');
        uZero = interpFunction(0.0);

    # Output the map to a CSV-format file.

    with open("ModeMap_PostProcess.csv", 'w') as outputWriter:
        outputWriterCSV = csv.writer(outputWriter, delimiter = ',', quotechar = '\"', quoting = csv.QUOTE_ALL);

        outputWriterCSV.writerow(["U_0 [eV]", uZero]);
        outputWriterCSV.writerow([]);

        outputWriterCSV.writerow(["Q [amu^1/2 A]", "U(Q) [eV]", "dU(Q) [meV]"]);

        for (_, q), (_, u) in zip(modulationModeCoordinates, modulationTotalEnergies):
            outputWriterCSV.writerow([q, u, 1000.0 * (u - uZero)]);

elif args.MapMode == '2D':
    # For 2D maps, arrange the data as a 2D array.

    q1Values = np.sort(
        np.unique([q1 for _, q1, _ in modulationModeCoordinates])
        );

    q2Values = np.sort(
        np.unique([q2 for _, _, q2 in modulationModeCoordinates])
        );

    # The data is arranged so that q1 -> columns and q2 -> rows.

    totalEnergy2D = np.zeros((len(q2Values), len(q1Values)), dtype = np.float64);

    # If we made the assumption of equally-spaced mode coordinates along each axis, the array-index calculations could probably be done a lot faster.
    # However, this is more general (e.g. would support non-uniform mode coordinates), and probably not perceptibly slower.

    for (_, q1, q2), (_, u) in zip(modulationModeCoordinates, modulationTotalEnergies):
        totalEnergy2D[np.argwhere(q2Values == q2), np.argwhere(q1Values == q1)] = u;

    # See whether Q1/Q2 = 0 are included in the maps - needed below.

    q1ZeroIndex, q2ZeroIndex = np.argwhere(q1Values == 0.0), np.argwhere(q2Values == 0.0);

    # As for 1D maps, to output relative energies, we need an energy zero.

    uZero = None;

    if len(q1ZeroIndex) == 1 and len(q2ZeroIndex) == 1:
        # Check whether Q1 = 0 and Q2 = 0 are included in the mapping points.

        uZero = totalEnergy2D[q2ZeroIndex, q1ZeroIndex];
    else:
        # If not, estimate the energy zero using bicubic interpolation.

        interpFunction = interp2d(q2Values, q1Values, totalEnergy2D, kind = 'cubic');
        uZero = (interpFunction(0.0, 0.0));

    # Both the array lookup and the interpolation return an array rather than a scalar.
    # This causes some of the array manipulations below to return results with counter-intuitive dimensions.
    # To avoid this, we explicitly cast to a (scalar) float.

    uZero = float(uZero);

    # Output the 2D map to a CSV-format file.

    with open("ModeMap_PostProcess_2DMap.csv", 'w') as outputWriter:
        outputWriterCSV = csv.writer(outputWriter, delimiter = ',', quotechar = '\"', quoting = csv.QUOTE_ALL);

        outputWriterCSV.writerow(["Q1 [amu^1/2 A]", "Columns"]);
        outputWriterCSV.writerow(["Q2 [amu^1/2 A]", "Rows"]);
        outputWriterCSV.writerow([]);

        outputWriterCSV.writerow(["U_0 [eV]", uZero]);
        outputWriterCSV.writerow([]);

        outputWriterCSV.writerow(["U(Q_1,Q_2) [eV]"]);
        outputWriterCSV.writerow([]);

        outputWriterCSV.writerow([""] + [q1 for q1 in q1Values]);

        for i, q2 in enumerate(q2Values):
            outputWriterCSV.writerow([q2] + [u for u in totalEnergy2D[i, :]]);

        outputWriterCSV.writerow([]);

        outputWriterCSV.writerow(["dU(Q_1,Q_2) [meV]"]);
        outputWriterCSV.writerow([]);

        outputWriterCSV.writerow([""] + [q1 for q1 in q1Values]);

        for i, q2 in enumerate(q2Values):
            outputWriterCSV.writerow([q2] + [1000.0 * (u - uZero) for u in totalEnergy2D[i, :]]);

    # If Q1/Q2 = 0 are included in the mapping points, we output either or both of these "slices" of the surface to a separate CSV file.

    headers, dataColumns = [], [];

    # q1ZeroIndex/q2ZeroIndex are 1x1 2D arrays, and attempting to extract a 1D slice of totalEnergy2D with them produces a 1x1xN 3D array.
    # To pass 1D arrays to the CSV output routine, we convert them to integers before using them as indices.

    if len(q2ZeroIndex) == 1:
        totalEnergySlice = totalEnergy2D[int(q2ZeroIndex), :];

        headers = headers + ["Q_1 [amu^1/2 A]", "U(Q_1) [eV]", "dU(Q_1) [meV]"];
        dataColumns = dataColumns + [q1Values, totalEnergySlice, 1000.0 * (totalEnergySlice - uZero)];

    if len(q1ZeroIndex) == 1:
        totalEnergySlice = totalEnergy2D[:, int(q1ZeroIndex)];

        headers = headers + ["Q_2 [amu^1/2 A]", "U(Q_2) [eV]", "dU(Q_2) [meV]"];
        dataColumns = dataColumns + [q2Values, totalEnergySlice, 1000.0 * (totalEnergySlice - uZero)];

    if len(headers) > 0:
        with open("ModeMap_PostProcess_1DProfiles.csv", 'w') as outputWriter:
            outputWriterCSV = csv.writer(outputWriter, delimiter = ',', quotechar = '\"', quoting = csv.QUOTE_ALL);

            outputWriterCSV.writerow(["U_0 [eV]", uZero]);
            outputWriterCSV.writerow([]);

            outputWriterCSV.writerow(headers);

            for i in range(0, max(len(column) for column in dataColumns)):
                outputWriterCSV.writerow(
                    [column[i] if i < len(column) else ""
                        for column in dataColumns
                        ]
                    );

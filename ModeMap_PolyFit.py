# ModeMap_PolyFit.py by J. M. Skelton


import argparse;
import csv;
import os;

import numpy as np;

import matplotlib as mpl;
import matplotlib.pyplot as plt;


# Parse command-line arguments.

parser = argparse.ArgumentParser(description = "Fit potential-energy surfaces output by ModeMap_PostProcess.py to polynomial functions");

# Defaults: read in data from a 1D map, fit to an 8-power polynomial, and output plots with the default x and y axis ranges.

parser.set_defaults(
    MapMode = '1D',
    PolynomialDegree = 8,
    PlotXLimits = None,
    PlotYLimits = None
    );

parser.add_argument(
    "--map_2d",
    dest = 'MapMode',
    action = 'store_const', const = '2D',
    help = "Fit the 1D surface profiles produced by post-processing a 2D mode map"
    );

parser.add_argument(
    "--degree", "--polynomial_degree",
    type = int, dest = 'PolynomialDegree',
    help = "Degree of the polynomial used to fit the surface profile(s) (default: 8)"
    );

parser.add_argument(
    "--plot_x",
    type = str, dest = 'PlotXLimits',
    help = "Optional x-axis range for the plot(s) in amu^1/2 Angstroms (e.g. \"-10 10\") (default: automatically determined)"
    );

parser.add_argument(
    "--plot_y",
    type = str, dest = 'PlotYLimits',
    help = "Optional y-axis range for the plot(s) in meV (e.g. \"-40 160\") (default: automatically determined)"
    );

args = parser.parse_args();

# Convert arguments and perform validation.

if args.PolynomialDegree <= 0:
    raise Exception("Error: The polynomial degree passed using --degree/--polynomial_degree must be greater than zero.");

if args.PlotXLimits != None:
    elements = [float(element) for element in args.PlotXLimits.strip().split()];

    if len(elements) != 2:
        raise Exception("Error: If supplied, the x-axis range must be specified by a pair of numbers.");

    args.PlotXLimits = elements;

if args.PlotYLimits != None:
    elements = [float(element) for element in args.PlotYLimits.strip().split()];

    if len(elements) != 2:
        raise Exception("Error: If supplied, the y-axis range must be specified by a pair of numbers.");

    args.PlotYLimits = elements;

# The input file to read data from depends on the mapping mode.

inputFile = None;

if args.MapMode == '1D':
    inputFile = "ModeMap_PostProcess.csv";
elif args.MapMode == '2D':
    inputFile = "ModeMap_PostProcess_1DProfiles.csv";

    if not os.path.isfile(inputFile):
        raise Exception("Error: If the 2D mapping mode is selected, the file \"ModeMap_PostProcess_1DProfiles.csv\" must be present in the working directory; this is only created by ModeMap_PostProcess.py if Q_1/Q_2 = 0 were mapped.");

# Read input data.

dataSets = [];

with open(inputFile, 'r') as inputReader:
    inputReaderCSV = csv.reader(inputReader);

    uZero = float(next(inputReaderCSV)[1]);

    next(inputReaderCSV);

    headerRow = next(inputReaderCSV);

    for i in range(0, len(headerRow), 3):
        dataSets.append(([], []));

    for row in inputReaderCSV:
        for i, (x, y) in enumerate(dataSets):
            baseIndex = i * 3;

            if baseIndex < len(row) and row[baseIndex] != "":
                x.append(float(row[baseIndex]));
                y.append(float(row[baseIndex + 2]));

# Convert the data for each potential surface to NumPy arrays.
# The dU(Q) values (energies with respect to the zero amplitude) read from the files are in units of meV, so we convert to eV before fitting, because our subsequent post processing accepts the fit coefficients in eV.

for i, (x, y) in enumerate(dataSets):
    dataSets[i] = (np.array(x, dtype = np.float64), np.array(y, dtype = np.float64) / 1000.0);

# Initialise Matplotlib.

fontSize = 8;
lineWidth = 0.5;

mpl.rc('font', **{ 'family' : 'serif', 'size' : fontSize, 'serif' : 'Times New Roman' });
mpl.rc('mathtext', **{ 'fontset' : 'custom', 'rm' : 'Times New Roman', 'it' : 'Times New Roman:italic', 'bf' : 'Times New Roman:bold' });

mpl.rc('axes', **{ 'labelsize' : fontSize, 'linewidth' : lineWidth });
mpl.rc('lines', **{ 'linewidth' : lineWidth, 'markeredgewidth' : lineWidth });

tickParams = { 'major.width' : lineWidth, 'minor.width' : lineWidth, 'direction' : 'in' };

mpl.rc('xtick', **tickParams);
mpl.rc('ytick', **tickParams);

# Process the data sets.

for i, (x, y) in enumerate(dataSets):
    print("Data set {0}:".format(i + 1));
    print("");

    # Fit the data.

    p = np.polyfit(x, y, deg = args.PolynomialDegree);

    # Print the coefficients.
    # They are printed one per line, and in ascening order, i.e. Q^0, Q^1, ... Q^N.
    # This output can be copy/pasted into the input file for the Schrodinger-solver post-processing code.

    for j, coefficient in enumerate(p[::-1]):
        print("{0: >23.16e}".format(coefficient));

    print("");

    # Write the fit coefficients and the original and fitted values plus the difference, in eV and meV, to a CSV-format file.

    fileNameStem = "ModeMap_PolyFit" if len(dataSets) == 1 else "ModeMap_PolyFit_{0}".format(i + 1);

    fileName = "{0}.csv".format(fileNameStem);

    print("  -> Outputting data to \"{0}\"".format(fileName));

    with open(fileName, 'w') as outputWriter:
        outputWriterCSV = csv.writer(outputWriter, delimiter = ',', quotechar = '\"', quoting = csv.QUOTE_ALL);

        outputWriterCSV.writerow(["Power", "Coefficient"]);

        for j, coefficient in enumerate(p[::-1]):
            outputWriterCSV.writerow([j, coefficient]);

        outputWriterCSV.writerow([]);

        outputWriterCSV.writerow(["Q [amu^1/2 A]", "dU(Q) [eV]", "dU(Q)_Fit [eV]", "d_Fit [eV]", "dU(Q) [meV]", "dU(Q)_Fit [meV]", "d_Fit [meV]"]);

        for q, u, uFit in zip(x, y, np.polyval(p, x)):
            dUFit = uFit - u;
            outputWriterCSV.writerow([q, u, uFit, dUFit, 1000.0 * u, 1000.0 * uFit, 1000.0 * dUFit]);

    # Plot a graph showing the fit overlaid on the original data points.

    fileName = "{0}.png".format(fileNameStem);

    print("  -> Outputting plot to \"{0}\"".format(fileName));

    plt.figure(figsize = (8.6 / 2.54, 6.0 / 2.54));

    plt.axhline(0.0, color = 'k');
    plt.axvline(0.0, color = 'k');

    plt.scatter(x, y * 1000.0, s = 15.0, marker = '^', edgecolor = 'b', facecolor = 'none', linewidth = 0.5);

    xValues = np.linspace(x[0], x[-1], 1000);
    plt.plot(xValues, np.polyval(p, xValues) * 1000.0, color = 'k');

    plt.xlabel(r"$Q$ [$\mathrm{amu}^{\frac{1}{2}} \AA$]");
    plt.ylabel(r"$\Delta U$($Q$) [meV]");

    if args.PlotXLimits != None:
        xMin, xMax = args.PlotXLimits;
        plt.xlim(xMin, xMax);
    else:
        plt.xlim(x[0], x[-1]);

    if args.PlotYLimits != None:
        yMin, yMax = args.PlotYLimits;
        plt.ylim(yMin, yMax);

    axes = plt.gca();

    axes.xaxis.set_ticks_position('both');
    axes.yaxis.set_ticks_position('both');

    axes.xaxis.grid(color = (211 / 255.0, 211 / 255.0, 211 / 255.0), dashes = (2.0, 1.0), linewidth = lineWidth);
    axes.yaxis.grid(color = (211 / 255.0, 211 / 255.0, 211 / 255.0), dashes = (2.0, 1.0), linewidth = lineWidth);

    axes.set_axisbelow(True);

    plt.tight_layout();

    plt.savefig(fileName, format = 'png', dpi = 300);

    print("");

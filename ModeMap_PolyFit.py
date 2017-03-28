# ModeMap_PolyFit.py by J. M. Skelton


import argparse;
import csv;
import math;
import os;

import numpy as np;

import matplotlib as mpl;
import matplotlib.pyplot as plt;

from mpl_toolkits.axes_grid.anchored_artists import AnchoredText;

from scipy.optimize import curve_fit;


# Conversion factors.

_EVToJoules = 1.60218e-19;
_AMUToKg = 1.66054e-27;
_AngstromToMetre = 1e-10;

_THzToInvCm = 33.35641;


# Harmonic function U(Q) = 1/2 * omega ** 2 * Q ** 2.

def _HarmonicFunction(q, omega):
    return 0.5 * omega ** 2 * q ** 2;


# Parse command-line arguments.

parser = argparse.ArgumentParser(description = "Fit potential-energy surfaces output by ModeMap_PostProcess.py to polynomial or harmonic functions");

# Defaults: read in data from a 1D map, fit to an 8-power polynomial, and output plots with the default x and y axis ranges.

parser.set_defaults(
    MapMode = '1D',
    PolynomialDegree = 8,
    HarmonicFit = False,
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
    "--harmonic_fit",
    dest = 'HarmonicFit',
    action = 'store_true',
    help = "Fit to a harmonic function U(Q) = 1/2 \omega ^ 2 Q ^ 2 (overrides --degree/--polynomial_degree)"
    );

parser.add_argument(
    "--degree", "--polynomial_degree",
    metavar = "n",
    type = int, dest = 'PolynomialDegree',
    help = "Degree of the polynomial used to fit the surface profile(s) (default: 8)"
    );

parser.add_argument(
    "--plot_x",
    metavar = "'x_min x_max'",
    type = str, dest = 'PlotXLimits',
    help = "Optional x-axis range for the plot(s) in amu^1/2 Angstroms (e.g. \"-10 10\") (default: automatically determined)"
    );

parser.add_argument(
    "--plot_y",
    metavar = "'y_min y_max'",
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

    # Fit the data.

    fitData = None;

    if args.HarmonicFit:
        # Perform a fit to the harmonic function.

        (omega), _ = curve_fit(_HarmonicFunction, x, y);

        # For single-parameter functions, curve_fit returns a 1D array with one element; this confuses some library functions, so it's easier to convert it to a float.

        omega = float(omega);

        # Convert to THz and inverse cm.

        omegaTHz = omega * math.sqrt(_EVToJoules / (_AMUToKg * _AngstromToMetre ** 2)) / (1.0e12 * 2.0 * math.pi);
        omegaInvCm = omegaTHz * _THzToInvCm;

        fitData = (omega, omegaTHz, omegaInvCm);
    else:
        fitData = np.polyfit(x, y, deg = args.PolynomialDegree);

    # If performing a polynomial fit, print the coefficients.

    if not args.HarmonicFit:
        # Coefficients are printed one per line, and in ascening order, i.e. Q^0, Q^1, ... Q^N.
        # This output can be copy/pasted into the input file for the Schrodinger-solver post-processing code.

        print("  -> Coefficients (ascending order, Q^0 ... Q^N):");
        print("");

        for j, coefficient in enumerate(fitData[::-1]):
            print("{0: >23.16e}".format(coefficient));

        print("");

    # Write the fit coefficients/fitted harmonic freqency and the original and fitted values plus the difference, in eV and meV, to a CSV-format file.

    fileNameStem = "ModeMap_PolyFit" if len(dataSets) == 1 else "ModeMap_PolyFit_{0}".format(i + 1);

    fileName = "{0}.csv".format(fileNameStem);

    print("  -> Outputting data to \"{0}\"".format(fileName));

    with open(fileName, 'w') as outputWriter:
        outputWriterCSV = csv.writer(outputWriter, delimiter = ',', quotechar = '\"', quoting = csv.QUOTE_ALL);

        if args.HarmonicFit:
            _, omegaTHz, omegaInvCm = fitData;

            outputWriterCSV.writerow(["\omega [THz]", omegaTHz]);
            outputWriterCSV.writerow(["\omega [cm^-1]", omegaInvCm]);
        else:
            outputWriterCSV.writerow(["Power", "Coefficient"]);

            for j, coefficient in enumerate(fitData[::-1]):
                outputWriterCSV.writerow([j, coefficient]);

        outputWriterCSV.writerow([]);

        outputWriterCSV.writerow(["Q [amu^1/2 A]", "dU(Q) [eV]", "dU(Q)_Fit [eV]", "d_Fit [eV]", "dU(Q) [meV]", "dU(Q)_Fit [meV]", "d_Fit [meV]"]);

        yFit = None;

        if args.HarmonicFit:
            omega, _, _ = fitData;
            yFit = _HarmonicFunction(x, omega);
        else:
            yFit = np.polyval(fitData, x);

        for q, u, uFit in zip(x, y, yFit):
            dUFit = uFit - u;
            outputWriterCSV.writerow([q, u, uFit, dUFit, 1000.0 * u, 1000.0 * uFit, 1000.0 * dUFit]);

    # Plot a graph showing the fit overlaid on the original data points.

    fileName = "{0}.png".format(fileNameStem);

    print("  -> Outputting plot to \"{0}\"".format(fileName));

    plt.figure(figsize = (8.6 / 2.54, 6.0 / 2.54));

    plt.axhline(0.0, color = 'k');
    plt.axvline(0.0, color = 'k');

    plt.scatter(x, y * 1000.0, s = 15.0, marker = '^', edgecolor = 'b', facecolor = 'none', linewidth = 0.5);

    xFit = np.linspace(x[0], x[-1], 1000);

    if args.HarmonicFit:
        omega, _, _ = fitData;
        plt.plot(xFit, _HarmonicFunction(xFit, omega) * 1000.0, color = 'k');
    else:
        plt.plot(xFit, np.polyval(fitData, xFit) * 1000.0, color = 'k');

    plt.xlabel(r"$Q$ [$\mathrm{amu}^{\frac{1}{2}} \mathrm{\AA}$]");
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

    # TODO: It _must_ be possible to set this using mpl.rc().

    axes.xaxis.set_ticks_position('both');
    axes.yaxis.set_ticks_position('both');

    axes.xaxis.grid(color = (211 / 255.0, 211 / 255.0, 211 / 255.0), dashes = (2.0, 1.0), linewidth = lineWidth);
    axes.yaxis.grid(color = (211 / 255.0, 211 / 255.0, 211 / 255.0), dashes = (2.0, 1.0), linewidth = lineWidth);

    axes.set_axisbelow(True);

    # If performing a harmonic fit, display the fitted frequency in a text box on the plot.

    if args.HarmonicFit:
        _, omegaTHz, omegaInvCm = fitData;

        anchoredText = AnchoredText(
            "Fit w/ $\omega$ = {0:.3f} THz ({1:.0f} cm$^{{-1}}$)".format(omegaTHz, omegaInvCm),
            loc = 9, frameon = True, prop = { 'size' : fontSize }
            );

        anchoredText.patch.set_linewidth(lineWidth);

        axes.add_artist(anchoredText);

    plt.tight_layout();

    plt.savefig(fileName, format = 'png', dpi = 300);

    print("");

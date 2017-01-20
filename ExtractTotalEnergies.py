# ExtractTotalEnergies.py by J. M. Skelton


import csv;
import os;
import re;


# This script assumes that the single-point calculations on the modulated structures were run in numbered directories, e.g. 001, 002, ..., NNN.
# It expects each directory matching this pattern to contain a VASP OSSICAR file with a total energy.

dirRegex = re.compile("^\d+$");

oszicarEnergyRegex = re.compile("E0= (?P<e0>[+-]?\.\d+E[+-]?\d+)");

# Extract total energies from the OSZICAR files in each numbered directory.

data = [];

for entry in os.listdir("./"):
	if os.path.isdir(entry) and dirRegex.match(entry):
		totalEnergy = None;

		inputReader = open(os.path.join(entry, "OSZICAR"), 'r');

		for line in inputReader:
			match = oszicarEnergyRegex.search(line);

			if match:
				# The calculations on the modulated structures should be single-point energy calculations.
				# If multiple total energies are found in an OSZICAR file, this is probably a sign that something is wrong, so we print a warning.

				if totalEnergy != None:
					print("WARNING: Multiple E0 values found in \"{0}\"".format(os.path.join(entry, "OSZICAR")));

				totalEnergy = float(match.group('e0'));

		inputReader.close();

		data.append((int(entry), totalEnergy));

# Write a CSV-format file listing the directory (moduation) numbers and extracted total energies.

outputWriter = open(r"ExtractTotalEnergies.csv", 'w');
outputWriterCSV = csv.writer(outputWriter, delimiter = ',', quotechar = '\"', quoting = csv.QUOTE_ALL);

outputWriterCSV.writerow(["Modulation #", "E0 [eV]"]);

for row in sorted(data):
	outputWriterCSV.writerow(row);

outputWriter.close();

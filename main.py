from src import readSequenceFile, getParameterDetails, pcaRegressionAlgorithm, motifsAlgorithm, processResults

import sys

filepath = "Test_Seq.txt"
try:
    f =open(filepath)
except NameError:
    filepath = str(raw_input("Please enter the input sequence file."))

# try:
#     f = open(writeFilePath)
# except NameError:
#     writeFilePath = str(raw_input("Please enter the output file path."))

sequence_map = readSequenceFile.readSequenceFile(filepath)
# print (sequence_map.keys())

parameter_map = {}
try:
    map_pca = getParameterDetails.iterateSequences(sequence_map)
    # print (parameter_map.keys())
except:
    print("Could not load parameters")

# map_pca = pcaRegressionAlgorithm.iterateSequences(parameter_map)

map_motif = motifsAlgorithm.iterateSequences(map_pca)

final_map = {}
for seq in map_pca:
    map_pca_1 = map_pca[seq]
    map_motif_1 = map_motif[seq]
    final_map[seq] = {}
    for start in map_pca_1:
        final_map[seq][start] = []
        map_pca_2 = map_pca_1[start]
        map_motif_2 = map_motif_1[start]
        map_pca_2 = map_pca_1[start]
        map_motif_2 = map_motif_1[start]

for seq in final_map:
    for start in final_map[seq]:
        if (final_map[seq][start][0] == 1):
            final_map[seq][start] = 1
        else:
            final_map[seq][start] = 0

final_map = processResults.predictSequencewiseTss(final_map)
pos_map = processResults.getTssPositions(final_map)



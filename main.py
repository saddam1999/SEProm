from src import readSequenceFile, getParameterDetails, pcaRegressionAlgorithm, motifsAlgorithm
import sys

filepath = "Test_Seq.txt"
try:
    f =open(filepath)
except NameError:
    filepath = str(raw_input("Please enter the input sequence file."))

try:
    f = open(writeFilePath)
except NameError:
    writeFilePath = str(raw_input("Please enter the output file path."))

sequence_map = readSequenceFile.readSequenceFile(filepath)
# print (sequence_map.keys())

parameter_map = {}
try:
    parameter_map = getParameterDetails.iterateSequences(sequence_map)
    # print (parameter_map.keys())
except:
    print("Could not load parameters")

map_pca = pcaRegressionAlgorithm.iterateSequences(parameter_map)
# print(map_pca)
map_motif = motifsAlgorithm.iterateSequences()
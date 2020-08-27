from src import readSequenceFile, getParameterDetails, pcaRegressionAlgorithm
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
# print (sequence_map)
keys = sequence_map.keys()
try:
    parameter_map = getParameterDetails.iterateSequences(keys)
except:
    print("Could not load parameters")

map_pca = pcaRegressionAlgorithm.iterateSequences(parameter_map)
# try:
#     getParameterDetails.iterateSequences(keys)
# except:
#     print("Error in calculating parameters")

import principalComponentAnalysis as pca
from lib import pcaEquations as pca_equations
from lib import regEquations as reg_equations
import logisticRegression as lr

MOVING_AVG_WINDOW_SIZE = 25
NO_TSS_WINDOW_LENGTH = 200
MOTIFS_NO_TSS_WINDOW_LENGTH = 200
ITR_WINDOW_SIZE = 100
SKIP_WINDOW = 1 #TODO 8 Previous 25
SKIP_WINDOW_SEQUENCE = 1
SKIP_WINDOW_RANGE = 40

WINDOW_40_PROB = 0.50
WINDOW_80_PROB = 0.50
WINDOW_100_PROB = 0.50
MOTIF_PROB = 0.50

RESULT_ITR_WINDOW = 5
RESULT_ITR_WINDOW_THRESH = 3
IGNORE_TSS_SEQ_THRESH = 35


def iterateSequences(param_map):
    final_result = {}
    normalized_map = param_map['normalized_params_map']
    # print(len(normalized_map[0]['aa']))
    for seq in normalized_map:
        final_result[seq] = iterate(seq,normalized_map)
    return final_result

def iterate(seq,normalized_map):
    print("Running PCA_Reg algo on sequence = ",seq)
    # normalized_map = param_map['normalized_params_map']
    params_map = normalized_map[seq]
    params = params_map.keys()
    # print(params)
    length = len(params_map[params[0]])
    # print (length)
    i = 0
    seq_40_map = {}
    seq_80_map = {}
    seq_100_map = {}

    for i in range(0,length-ITR_WINDOW_SIZE-ITR_WINDOW_SIZE-NO_TSS_WINDOW_LENGTH, SKIP_WINDOW):
        tss_motif_start = i
        tss_motif_stop = tss_motif_start + ITR_WINDOW_SIZE
        no_tss_motif_start = tss_motif_stop + NO_TSS_WINDOW_LENGTH
        no_tss_motif_stop = no_tss_motif_start + ITR_WINDOW_SIZE
        tss_window_40_arr = extractWindowx(tss_motif_start, tss_motif_stop, params, params_map,40)
        no_tss_window_40_arr = extractWindowx(no_tss_motif_start, no_tss_motif_stop, params, params_map,40)
        # print("-------------------DONE--------------------")

        tss_window_80_arr = extractWindowx(tss_motif_start, tss_motif_stop, params, params_map,80)
        no_tss_window_80_arr = extractWindowx(no_tss_motif_start, no_tss_motif_stop, params, params_map,80)

        tss_window_100_arr = extractWindowx(tss_motif_start, tss_motif_stop, params, params_map,80)
        no_tss_window_100_arr = extractWindowx(no_tss_motif_start, no_tss_motif_stop, params, params_map,80)

        seq_40_map[tss_motif_start] = []
        seq_40_map[tss_motif_start].append(tss_window_40_arr)
        seq_40_map[tss_motif_start].append(no_tss_window_40_arr)

        seq_80_map[tss_motif_start] = []
        seq_80_map[tss_motif_start].append(tss_window_80_arr)
        seq_80_map[tss_motif_start].append(no_tss_window_80_arr)

        seq_100_map[tss_motif_start] = []
        seq_100_map[tss_motif_start].append(tss_window_100_arr)
        seq_100_map[tss_motif_start].append(no_tss_window_100_arr)
    print (seq_40_map[tss_motif_start][0])
    import sys
    sys.exit()
    # print (no_tss_motif_stop)
    return predictPCA(seq, seq_40_map, seq_80_map, seq_100_map)

def extractWindow(start, stop, params, params_map):
    arr = []
    for p in params:
        param_arr = params_map[p]
        # print (param_arr)
        sum = 0.0
        length = 0
        itr_undefined = 0
        for i in range(start, stop):
            # print (i)
            try:
                param_arr[i]
            except:
                print("index error at", i)
            if param_arr[i]==None:
                itr_undefined += 1
            sum = sum+ param_arr[i]
            length+=1
        # print (float(sum)/float(length))
        number = round(sum / length,6)
        # print (number)
        arr.append(round(number / length,6))
    # print (arr)
    return arr

def extractWindowx(motif_start, motif_stop, params, params_map,x):
    start = motif_stop - x
    stop = motif_stop
    return extractWindow(start, stop, params, params_map)

def predictPCA(seq, seq_40_map, seq_80_map, seq_100_map):
    for start in seq_40_map:
        seq_40_map[start][0] = pca.getPCAs(seq_40_map[start][0], pca_equations.window_40)
        seq_40_map[start][1] = pca.getPCAs(seq_40_map[start][1], pca_equations.window_40)
    for start in seq_80_map:
        seq_80_map[start][0] = pca.getPCAs(seq_80_map[start][0], pca_equations.window_80)
        seq_80_map[start][1] = pca.getPCAs(seq_80_map[start][1], pca_equations.window_80)
    for start in seq_100_map:
        seq_100_map[start][0] = pca.getPCAs(seq_100_map[start][0], pca_equations.window_100)
        seq_100_map[start][1] = pca.getPCAs(seq_100_map[start][1], pca_equations.window_100)
    return predictRegression(seq, seq_40_map, seq_80_map, seq_100_map)

def predictRegression(seq, seq_40_map, seq_80_map, seq_100_map):
    for start in seq_40_map:
        seq_40_map[start][0] = lr.predict(seq_40_map[start][0], reg_equations.window_40, WINDOW_40_PROB)
        seq_40_map[start][1] = lr.predict(seq_40_map[start][1], reg_equations.window_40, WINDOW_40_PROB)
    for start in seq_80_map:
        seq_80_map[start][0] = lr.predict(seq_80_map[start][0], reg_equations.window_80, WINDOW_80_PROB)
        seq_80_map[start][1] = lr.predict(seq_80_map[start][1], reg_equations.window_80, WINDOW_80_PROB)
    for start in seq_100_map:
        seq_100_map[start][0] = lr.predict(seq_100_map[start][0], reg_equations.window_100, WINDOW_100_PROB)
        seq_100_map[start][1] = lr.predict(seq_100_map[start][1], reg_equations.window_100, WINDOW_100_PROB)

    return processResults(seq, seq_40_map, seq_80_map, seq_100_map)

def processResults(seq, seq_40_map, seq_80_map, seq_100_map):
    combined_result_map = {}
    for start in seq_40_map:
        res_40_tss = seq_40_map[start][0]
        res_40_notss = seq_40_map[start][1]
        res_80_tss = seq_80_map[start][0]
        res_80_notss = seq_80_map[start][1]
        res_100_tss = seq_100_map[start][0]
        res_100_notss = seq_100_map[start][1]
        combined_result_map[start] = []
        sum_res_tss = res_40_tss + res_80_tss + res_100_tss
        sum_res_notss = res_40_notss + res_80_notss + res_100_notss
        combined_result_map[start].append(1 if sum_res_tss >= 2 else 0)
        combined_result_map[start].append(1 if sum_res_notss >= 2 else 0)
    return combined_result_map



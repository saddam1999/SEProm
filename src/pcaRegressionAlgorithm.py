# sequence_map
# import main
import sys

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
    for seq in normalized_map:
        final_result[seq] = iterate(seq,param_map)
    return final_result

def iterate(seq,param_map):
    print("Running PCA_Reg algo on sequence = ",seq)
    normalized_map = param_map['normalized_params_map']
    params_map = normalized_map[seq]
    params = params_map.keys()
    length = len(params_map[params[0]])
    i = 0
    seq_40_map = {}
    seq_80_map = {}
    seq_100_map = {}

    while (i + ITR_WINDOW_SIZE) <= (length + ITR_WINDOW_SIZE - 1):
        tss_motif_start = i
        tss_motif_stop = tss_motif_start + ITR_WINDOW_SIZE
        no_tss_motif_start = tss_motif_stop + NO_TSS_WINDOW_LENGTH
        no_tss_motif_stop = no_tss_motif_start + ITR_WINDOW_SIZE
        tss_window_40_arr = extractWindow40(tss_motif_start, tss_motif_stop, params, params_map)
        no_tss_window_40_arr = extractWindow40(no_tss_motif_start, no_tss_motif_stop, params, params_map)

        tss_window_80_arr = extractWindow80(tss_motif_start, tss_motif_stop, params, params_map)
        no_tss_window_80_arr = extractWindow80(no_tss_motif_start, no_tss_motif_stop, params, params_map)

        tss_window_100_arr = extractWindow100(tss_motif_start, tss_motif_stop, params, params_map)
        no_tss_window_100_arr = extractWindow100(no_tss_motif_start, no_tss_motif_stop, params, params_map)

        seq_40_map[tss_motif_start] = []
        seq_40_map[tss_motif_start].append(tss_window_40_arr)
        seq_40_map[tss_motif_start].append(no_tss_window_40_arr)

        seq_80_map[tss_motif_start] = []
        seq_80_map[tss_motif_start].append(tss_window_80_arr)
        seq_80_map[tss_motif_start].append(no_tss_window_80_arr)

        seq_100_map[tss_motif_start] = []
        seq_100_map[tss_motif_start].append(tss_window_100_arr)
        seq_100_map[tss_motif_start].append(no_tss_window_100_arr)

        i += SKIP_WINDOW
    return predictPCA(seq, seq_40_map, seq_80_map, seq_100_map)

def extractWindow(start, stop, params, params_map):
    arr = []
    for p in params:
        param_arr = params_map[p]
        sum = 0
        length = 0
        itr_undefined = 0
        for i in range(start, stop+1):
            if param_arr[i]:
                continue
            else:
                itr_undefined+=1
                sum += param_arr[i]
            length+=1

        number = int((sum / length).toFixed(6))
        arr.append(int((sum / length).toFixed(6)))
    return arr

def extractWindow80(motif_start, motif_stop, params, params_map):
    start = motif_stop - 80
    stop = motif_stop
    return extractWindow(start, stop, params, params_map)




from lib import constants

def predictSequencewiseTss(final_map):
    result_map = {}
    for seq in final_map:
        start_arr = (final_map[seq]).keys()
        result_map[seq] = {}
        for i in range(len(start_arr)):
            tss = 0
            for j in range(i, i+constants.RESULT_ITR_WINDOW):
                if start_arr[j] and final_map[seq][start_arr[j]]:
                    tss += final_map[seq][start_arr[j]]
            if (tss >= constants.RESULT_ITR_WINDOW_THRESH):
                result_map[seq][start_arr[i]] = 1
            else:
                result_map[seq][start_arr[i]] = 0
    return result_map

def getTssPositions():

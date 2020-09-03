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

motif_positions = [['486-493','463-472','420-437','720-737','763-772','786-793'],
				['487-494','464-471','453-462','440-452','787-794','764-771','753-762','740-752'],
				['483-494','463-472','450-459','434-445','783-794','763-772','750-759','734-745'],
				['487-498','463-472','450-459','432-445','787-798','763-772','750-759','732-745']]

motifs_size = [[8,10,18],[8,8,10,13],[12,10,10,12],[12,10,10,14]]

motifs_dist_end = [[19,42,85],[18,41,52,65],[22,42,55,71],[18,42,55,73]]

def extractWindow(tss_start, tss_stop, motif, dist_end, struct_ener_map):
    motif_start = tss_stop - dist_end
    motif_stop = motif_start + motif - 1
    params_arr = []
    for k in struct_ener_map:
        obj = struct_ener_map[k]
        sum = 0
        length = 0
        itr_undefined = 0
        for i in range(motif_start, motif_stop+1):
            if obj[i] == None:
                # obj[i]
                itr_undefined += 1
            sum+=float(obj[i])
            length+=1
        params_arr.append(sum/length)
    return params_arr

def transformStructMap(struct_map):
    new_map = {}
    for seq in struct_map:
        map = struct_map[seq]
        for k in map:
            if new_map[k] == None:
                new_map[k] =[]
                for obj in map[k]:
                    new_map[k]=(obj)
    return new_map

def iterateSequences(parameters_map):
    final_result = {}
    combined_params_map = parameters_map['combined_params_map']
    for seq in combined_params_map:
        # print(seq)
        final_result[seq] = iterate(seq, parameters_map)
    return final_result

def iterate(seq, parameters_map):
    print("Running Motifs algorithm on sequence: ",seq)
    motif_map = {}
    struct_ener_map = parameters_map['combined_params_map'][seq]
    # print(seq)
    params = struct_ener_map.keys()
    print (params)

    length = len(struct_ener_map[params[0]].keys())
    print(length)
    for i in range(0,length):
        tss_motif_start = i
        tss_motif_stop = tss_motif_start + ITR_WINDOW_SIZE
        no_tss_motif_start = tss_motif_stop + MOTIFS_NO_TSS_WINDOW_LENGTH
        no_tss_motif_stop = no_tss_motif_start + ITR_WINDOW_SIZE
        # print(tss_motif_start, tss_motif_stop, no_tss_motif_stop, no_tss_motif_start)
        if motif_map[tss_motif_start] == None:
            motif_map[tss_motif_start] = {}

        # len = Object.keys(struct_ener_map[params[0]]).length;


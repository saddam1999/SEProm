import principalComponentAnalysis as pca
from lib import pcaEquations as pca_equations
from lib import regEquations as reg_equations
import src.logisticRegression as lr

MOVING_AVG_WINDOW_SIZE = 25
NO_TSS_WINDOW_LENGTH = 200
MOTIFS_NO_TSS_WINDOW_LENGTH = 200
ITR_WINDOW_SIZE = 100
SKIP_WINDOW = 100 #TODO 8 Previous 25
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
    # print (struct_ener_map)
    # print(tss_stop, dist_end,motif)
    motif_start = tss_stop - dist_end
    motif_stop = motif_start + motif
    # print(motif_start,motif_stop)
    params_arr = []
    for k in struct_ener_map.keys():
        obj = struct_ener_map[k]
        sum = 0
        length = 0
        itr_undefined = 0
        for i in range(motif_start, motif_stop):
            try:
                obj[i]
                sum += float(obj[i])
                length += 1
            except:
                itr_undefined += 1
        if length!=0:
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
    for seq in parameters_map['combined_params_map']:
        final_result[seq] = iterate(seq, parameters_map)
        import sys
        sys.exit()
    return final_result

def iterate(seq, parameters_map):
    print("Running Motifs algorithm on sequence: ",seq)
    motif_map = {}
    struct_ener_map = parameters_map['combined_params_map'][seq]
    # print(struct_ener_map)
    params = struct_ener_map.keys()
    # print(params)
    length = len(struct_ener_map[params[0]].keys())
    i=0
    # arr = [i for i in range(length)]
    # print(arr)
    # import sys
    # sys.exit()
    no_tss_motif_stop = 0
    # while(no_tss_motif_stop<length-1):

    for i in range(0, length - ITR_WINDOW_SIZE - ITR_WINDOW_SIZE - NO_TSS_WINDOW_LENGTH, SKIP_WINDOW):
        # print (i)
        tss_motif_start = i
        tss_motif_stop = tss_motif_start + ITR_WINDOW_SIZE
        no_tss_motif_start = tss_motif_stop + MOTIFS_NO_TSS_WINDOW_LENGTH
        no_tss_motif_stop = no_tss_motif_start + ITR_WINDOW_SIZE
        # print(tss_motif_start, tss_motif_stop, no_tss_motif_stop, no_tss_motif_start)
        try:
            (motif_map[tss_motif_start])
        except:
            motif_map[tss_motif_start] = {}
        for k in range(len(motifs_size)):
            motif_size = motifs_size[k]
            motif_dist_end = motifs_dist_end[k]
            for j in range(len(motif_size)):
                motif = motif_size[j]
                dist_end = motif_dist_end[j]
                try:
                    motif_map[tss_motif_start]['m_'+str(k)]
                except:
                    motif_map[tss_motif_start]['m_' + str(k)] = {}
                try:
                    (motif_map[tss_motif_start]['m_' + str(k)][j])
                except:
                    motif_map[tss_motif_start]['m_' + str(k)][j] = {}

                motif_map[tss_motif_start]['m_' + str(k)][j]['tss'] = extractWindow(tss_motif_start, tss_motif_stop,motif, dist_end, struct_ener_map)
                motif_map[tss_motif_start]['m_' + str(k)][j]['no_tss'] = extractWindow(no_tss_motif_start,no_tss_motif_stop, motif, dist_end,struct_ener_map)

    predictPCA(seq,motif_map)

def predictPCA(seq, motif_map):

    for start in motif_map.keys():
        motif_map[start]['m_0'][0]['tss'] = pca.getPCAs(motif_map[start]['m_0'][0]['tss'], pca_equations.m_0_0)
        motif_map[start]['m_0'][0]['no_tss'] = pca.getPCAs(motif_map[start]['m_0'][0]['no_tss'], pca_equations.m_0_0)

        motif_map[start]['m_0'][1]['tss'] = pca.getPCAs(motif_map[start]['m_0'][1]['tss'], pca_equations.m_0_1)
        motif_map[start]['m_0'][1]['no_tss'] = pca.getPCAs(motif_map[start]['m_0'][1]['no_tss'], pca_equations.m_0_1)
        motif_map[start]['m_0'][2]['tss'] = pca.getPCAs(motif_map[start]['m_0'][2]['tss'], pca_equations.m_0_2)
        motif_map[start]['m_0'][2]['no_tss'] = pca.getPCAs(motif_map[start]['m_0'][2]['no_tss'], pca_equations.m_0_2)

        motif_map[start]['m_1'][0]['tss'] = pca.getPCAs(motif_map[start]['m_1'][0]['tss'], pca_equations.m_1_0)
        motif_map[start]['m_1'][0]['no_tss'] = pca.getPCAs(motif_map[start]['m_1'][0]['no_tss'], pca_equations.m_1_0)
        motif_map[start]['m_1'][1]['tss'] = pca.getPCAs(motif_map[start]['m_1'][1]['tss'], pca_equations.m_1_1)
        motif_map[start]['m_1'][1]['no_tss'] = pca.getPCAs(motif_map[start]['m_1'][1]['no_tss'], pca_equations.m_1_1)
        motif_map[start]['m_1'][2]['tss'] = pca.getPCAs(motif_map[start]['m_1'][2]['tss'], pca_equations.m_1_2)
        motif_map[start]['m_1'][2]['no_tss'] = pca.getPCAs(motif_map[start]['m_1'][2]['no_tss'], pca_equations.m_1_2)
        motif_map[start]['m_1'][3]['tss'] = pca.getPCAs(motif_map[start]['m_1'][3]['tss'], pca_equations.m_1_3)
        motif_map[start]['m_1'][3]['no_tss'] = pca.getPCAs(motif_map[start]['m_1'][3]['no_tss'], pca_equations.m_1_3)

        motif_map[start]['m_2'][0]['tss'] = pca.getPCAs(motif_map[start]['m_2'][0]['tss'], pca_equations.m_2_0)
        motif_map[start]['m_2'][0]['no_tss'] = pca.getPCAs(motif_map[start]['m_2'][0]['no_tss'], pca_equations.m_2_0)
        motif_map[start]['m_2'][1]['tss'] = pca.getPCAs(motif_map[start]['m_2'][1]['tss'], pca_equations.m_2_1)
        motif_map[start]['m_2'][1]['no_tss'] = pca.getPCAs(motif_map[start]['m_2'][1]['no_tss'], pca_equations.m_2_1)
        motif_map[start]['m_2'][2]['tss'] = pca.getPCAs(motif_map[start]['m_2'][2]['tss'], pca_equations.m_2_2)
        motif_map[start]['m_2'][2]['no_tss'] = pca.getPCAs(motif_map[start]['m_2'][2]['no_tss'], pca_equations.m_2_2)
        motif_map[start]['m_2'][3]['tss'] = pca.getPCAs(motif_map[start]['m_2'][3]['tss'], pca_equations.m_2_3)
        motif_map[start]['m_2'][3]['no_tss'] = pca.getPCAs(motif_map[start]['m_2'][3]['no_tss'], pca_equations.m_2_3)

        motif_map[start]['m_3'][0]['tss'] = pca.getPCAs(motif_map[start]['m_3'][0]['tss'], pca_equations.m_3_0)
        motif_map[start]['m_3'][0]['no_tss'] = pca.getPCAs(motif_map[start]['m_3'][0]['no_tss'], pca_equations.m_3_0)
        motif_map[start]['m_3'][1]['tss'] = pca.getPCAs(motif_map[start]['m_3'][1]['tss'], pca_equations.m_3_1)
        motif_map[start]['m_3'][1]['no_tss'] = pca.getPCAs(motif_map[start]['m_3'][1]['no_tss'], pca_equations.m_3_1)
        motif_map[start]['m_3'][2]['tss'] = pca.getPCAs(motif_map[start]['m_3'][2]['tss'], pca_equations.m_3_2)
        motif_map[start]['m_3'][2]['no_tss'] = pca.getPCAs(motif_map[start]['m_3'][2]['no_tss'], pca_equations.m_3_2)
        motif_map[start]['m_3'][3]['tss'] = pca.getPCAs(motif_map[start]['m_3'][3]['tss'], pca_equations.m_3_3)
        motif_map[start]['m_3'][3]['no_tss'] = pca.getPCAs(motif_map[start]['m_3'][3]['no_tss'], pca_equations.m_3_3)
    return predictRegression(seq,motif_map)

def predictRegression(seq,motif_map):
    # print (motif_map.keys())
    for start in motif_map.keys():
        motif_map[start]['m_0'][0]['tss'] = lr.predict(motif_map[start]['m_0'][0]['tss'], reg_equations.m_0_0,MOTIF_PROB)
        motif_map[start]['m_0'][0]['no_tss'] = lr.predict(motif_map[start]['m_0'][0]['no_tss'], reg_equations.m_0_0,MOTIF_PROB)

        motif_map[start]['m_0'][1]['tss'] = lr.predict(motif_map[start]['m_0'][1]['tss'], reg_equations.m_0_1,MOTIF_PROB)
        motif_map[start]['m_0'][1]['no_tss'] = lr.predict(motif_map[start]['m_0'][1]['no_tss'], reg_equations.m_0_1,MOTIF_PROB)

        motif_map[start]['m_0'][2]['tss'] = lr.predict(motif_map[start]['m_0'][2]['tss'], reg_equations.m_0_2,MOTIF_PROB)
        motif_map[start]['m_0'][2]['no_tss'] = lr.predict(motif_map[start]['m_0'][2]['no_tss'], reg_equations.m_0_2,MOTIF_PROB)

        motif_map[start]['m_1'][0]['tss'] = lr.predict(motif_map[start]['m_1'][0]['tss'], reg_equations.m_1_0,MOTIF_PROB)
        motif_map[start]['m_1'][0]['no_tss'] = lr.predict(motif_map[start]['m_1'][0]['no_tss'], reg_equations.m_1_0,MOTIF_PROB)

        motif_map[start]['m_1'][1]['tss'] = lr.predict(motif_map[start]['m_1'][1]['tss'], reg_equations.m_1_1,MOTIF_PROB)
        motif_map[start]['m_1'][1]['no_tss'] = lr.predict(motif_map[start]['m_1'][1]['no_tss'], reg_equations.m_1_1,MOTIF_PROB)

        motif_map[start]['m_1'][2]['tss'] = lr.predict(motif_map[start]['m_1'][2]['tss'], reg_equations.m_1_2,MOTIF_PROB)
        motif_map[start]['m_1'][2]['no_tss'] = lr.predict(motif_map[start]['m_1'][2]['no_tss'], reg_equations.m_1_2,MOTIF_PROB)

        motif_map[start]['m_1'][3]['tss'] = lr.predict(motif_map[start]['m_1'][3]['tss'], reg_equations.m_1_3,MOTIF_PROB)
        motif_map[start]['m_1'][3]['no_tss'] = lr.predict(motif_map[start]['m_1'][3]['no_tss'], reg_equations.m_1_3,MOTIF_PROB)

        motif_map[start]['m_2'][0]['tss'] = lr.predict(motif_map[start]['m_2'][0]['tss'], reg_equations.m_2_0,MOTIF_PROB)
        motif_map[start]['m_2'][0]['no_tss'] = lr.predict(motif_map[start]['m_2'][0]['no_tss'], reg_equations.m_2_0,MOTIF_PROB)

        motif_map[start]['m_2'][1]['tss'] = lr.predict(motif_map[start]['m_2'][1]['tss'], reg_equations.m_2_1,MOTIF_PROB)
        motif_map[start]['m_2'][1]['no_tss'] = lr.predict(motif_map[start]['m_2'][1]['no_tss'], reg_equations.m_2_1,MOTIF_PROB)

        motif_map[start]['m_2'][2]['tss'] = lr.predict(motif_map[start]['m_2'][2]['tss'], reg_equations.m_2_2,MOTIF_PROB)
        motif_map[start]['m_2'][2]['no_tss'] = lr.predict(motif_map[start]['m_2'][2]['no_tss'], reg_equations.m_2_2,MOTIF_PROB)

        motif_map[start]['m_2'][3]['tss'] = lr.predict(motif_map[start]['m_2'][3]['tss'], reg_equations.m_2_3,MOTIF_PROB)
        motif_map[start]['m_2'][3]['no_tss'] = lr.predict(motif_map[start]['m_2'][3]['no_tss'], reg_equations.m_2_3,MOTIF_PROB)

        motif_map[start]['m_3'][0]['tss'] = lr.predict(motif_map[start]['m_3'][0]['tss'], reg_equations.m_3_0,MOTIF_PROB)
        motif_map[start]['m_3'][0]['no_tss'] = lr.predict(motif_map[start]['m_3'][0]['no_tss'], reg_equations.m_3_0,MOTIF_PROB)

        motif_map[start]['m_3'][1]['tss'] = lr.predict(motif_map[start]['m_3'][1]['tss'], reg_equations.m_3_1,MOTIF_PROB)
        motif_map[start]['m_3'][1]['no_tss'] = lr.predict(motif_map[start]['m_3'][1]['no_tss'], reg_equations.m_3_1,MOTIF_PROB)

        motif_map[start]['m_3'][2]['tss'] = lr.predict(motif_map[start]['m_3'][2]['tss'], reg_equations.m_3_2,MOTIF_PROB)
        motif_map[start]['m_3'][2]['no_tss'] = lr.predict(motif_map[start]['m_3'][2]['no_tss'], reg_equations.m_3_2,MOTIF_PROB)

        motif_map[start]['m_3'][3]['tss'] = lr.predict(motif_map[start]['m_3'][3]['tss'], reg_equations.m_3_3,MOTIF_PROB)
        motif_map[start]['m_3'][3]['no_tss'] = lr.predict(motif_map[start]['m_3'][3]['no_tss'], reg_equations.m_3_3,MOTIF_PROB)
    # print (motif_map)
    return processResults(seq,motif_map)

def processResults(seq,motif_map):
    combined_result_map = {}
    for start in motif_map.keys():
        row_map = motif_map[start]
        combined_result_map[start] = {}
        # print(row_map.keys())
        for row in row_map.keys():
            m_map = row_map[row]
            # print(m_map)
            combined_result_map[start][row] = {}
            threshold = len(m_map.keys()) - 1
            tss = 0
            no_tss = 0
            for motif in m_map:
                tss += m_map[motif]['tss']
                no_tss += m_map[motif]['no_tss']
            combined_result_map[start][row]['tss'] = (1 if (tss >= threshold) else 0)
            combined_result_map[start][row]['no_tss'] = (1 if (no_tss >= threshold) else 0)
        # print(combined_result_map)
        # import sys
        # sys.exit()
    for start in combined_result_map:
        tss = 0
        no_tss = 0
        result = []
        for row in combined_result_map[start]:
            tss += combined_result_map[start][row]['tss']
            no_tss += combined_result_map[start][row]['no_tss']
        result.append(1 if (tss > 0) else 0)
        result.append(1 if (no_tss > 0) else 0)
        combined_result_map[start] = result
    print(combined_result_map)
    return combined_result_map

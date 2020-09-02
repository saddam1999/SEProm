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

def iterateSequences():

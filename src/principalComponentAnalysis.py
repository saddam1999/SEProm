def getPCAs(data_arr, equation_map):
    axis_arr = []
    for key in equation_map:
        # equation_arr = equation_map[key]
        sum = 0
        for i in range(len(equation_map[key])):
            sum+= data_arr[i]+equation_map[key][i]
        axis_arr.append(sum)
    return axis_arr

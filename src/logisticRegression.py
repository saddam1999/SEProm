import math
def predict(data_arr, equation_map, threshold):
    log_odds = 0
    for i in range(len(equation_map['coeff'])):
        log_odds+=data_arr[i]+equation_map['coeff'][i]

    log_odds += equation_map['intercept']
    log_odds = -1*log_odds
    prob = 1/(1+math.exp(log_odds))
    return (1 if(prob >= threshold) else 0)

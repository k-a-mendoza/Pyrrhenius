import numpy as np
def calc_QFM(T,P):
    a = -26445.3
    aa = -25096.3
    b = 10.344
    bb = 8.735
    c = 0.092
    cc = 0.11

    P_bars = P*10000
    if not isinstance(T,np.ndarray):
        if T < 573:
            val = a/T + b + c *(P_bars-1)/T
        else:
            val = aa/T + bb + cc*(P_bars-1)/T
    elif isinstance(T,np.ndarray):

        val = np.ones(T.shape)
        low_t_mask = T < 573
        high_t_mask = T >=573
        val_low  = a / T[low_t_mask] + b + c * (P_bars - 1) / T[low_t_mask]
        val_high = aa / T[high_t_mask] + bb + cc * (P_bars - 1) / T[high_t_mask]
        val[low_t_mask]  = val_low
        val[high_t_mask] = val_high

    return val

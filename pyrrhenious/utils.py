import numpy as np
def calc_QFM(T,P):
    P_bars = P*10000
    if isinstance(T,np.ndarray) and isinstance(P,np.ndarray):
        log_fO2 = np.ones(T.shape)
        log_fO2[T<573]  = (-26445.3/T[T<573] ) + 10.344 + 0.092 * (P_bars[T<573]-1)/T[T<573]
        log_fO2[T>=573] = (-25096.3/T[T>=573]) + 8.735  + 0.11  * (P_bars[T>=573]-1)/T[T>=573]
    elif isinstance(T,np.ndarray) and not isinstance(P,np.ndarray):
        log_fO2 = np.ones(T.shape)
        log_fO2[T<573]  = (-26445.3/T[T<573] ) + 10.344 + 0.092 * (P_bars-1)/T[T<573]
        log_fO2[T>=573] = (-25096.3/T[T>=573]) + 8.735  + 0.11  * (P_bars-1)/T[T>=573]
    elif not isinstance(T,np.ndarray) and isinstance(P,np.ndarray):
        log_fO2 = np.ones(P.shape)
        if T > 573:
            log_fO2 = (-26445.3/T ) + 10.344 + 0.092 * (P_bars-1)/T
        else:
            log_fO2 = (-25096.3/T) + 8.735  + 0.11  * (P_bars-1)/T
    else:
        if T > 573:
            log_fO2 = (-26445.3/T ) + 10.344 + 0.092 * (P_bars-1)/T
        else:
            log_fO2 = (-25096.3/T) + 8.735  + 0.11  * (P_bars-1)/T
    return log_fO2
import numpy as np
cimport numpy as np
import cython

cdef extern from "math.h":
    float fabs(float t)

@cython.wraparound(False)
@cython.boundscheck(False)
def medoid(np.ndarray[np.float32_t, ndim=4] x, float fill_val=255.):
    cdef int i, j, k, m, n
    cdef float dist
    cdef float min_dist
    cdef float dist_sum
    cdef int argmin_dist
    cdef int na_count

    cdef np.ndarray[np.float32_t, ndim=2] dist_matrix
    dist_matrix = np.zeros((x.shape[0], x.shape[0]), dtype=x.dtype)    
    
    cdef np.ndarray[np.float32_t, ndim=3] out
    out = np.empty((x.shape[1], x.shape[2], x.shape[3]), dtype=x.dtype)

    for m in range(x.shape[2]):
        for n in range(x.shape[3]):
            na_count = 0
            for i in range(x.shape[0]):
                for k in range(x.shape[1]):
                    if x[i, k, m, n] == fill_val:
                        na_count += 1
                        break

            if na_count == x.shape[0]:
                for k in range(x.shape[1]):
                    out[k, m, n] = fill_val
                continue
            elif x.shape[0] - na_count < 3:
                i_data = 0
                for i in range(x.shape[0]):
                    found_na = False
                    for k in range(x.shape[1]):
                        if x[i, k, m, n] == fill_val:
                            found_na = True
                            break

                    if not found_na:
                       i_data = i
                       break

                for k in range(x.shape[1]):
                    out[k, m, n] = x[i_data, k, m, n]
                continue

            argmin_dist = -1
            min_dist = 0.

            for i in range(x.shape[0]):
                for j in range(i+1, x.shape[0]):
                    dist = 0.
                    for k in range(x.shape[1]):
                        if x[i, k, m, n] == fill_val or x[j, k, m, n] == fill_val:
                            dist = 255.
                            break
                        
                        dist += abs(x[i, k, m, n] - x[j, k, m, n])
                    
                    dist_matrix[i, j] = dist
                    dist_matrix[j, i] = dist
 
            for i in range(x.shape[0]):
                if x[i, 0, m, n] == fill_val:
                  continue
                dist_sum = 0.
                for j in range(x.shape[0]):
                    dist_sum += dist_matrix[i, j]

                if argmin_dist < 0:
                    argmin_dist = i
                    min_dist = dist_sum
                else:
                    if dist_sum < min_dist:
                        min_dist = dist_sum
                        argmin_dist = i

            for k in range(x.shape[1]):
              out[k, m, n] = x[argmin_dist, k, m, n]    
            
    return out        







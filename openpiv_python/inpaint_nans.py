import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

def inpaint_nans(A, method=0):
    """
    Interpolate NaN values in an array

    Parameters:
    -----------
    A : ndarray
        Array with NaN values to be filled
    method : int, optional
        Interpolation method to use:
        0 - (DEFAULT) Uses del^2 operator, builds smaller system for few NaNs
        1 - Uses del^2 operator over entire array
        2 - Uses del^2 operator, direct solve
        3 - Uses del^4 operator
        4 - Uses spring metaphor
        5 - Uses average of 8 nearest neighbors

    Returns:
    --------
    B : ndarray
        Array with NaN values filled in

    Notes:
    ------
    Solves approximation to one of several PDEs to interpolate and
    extrapolate holes in an array. All methods are capable of extrapolation,
    some are better than others. There are also speed and accuracy differences.

    Original MATLAB implementation by John D'Errico
    e-mail: woodchips@rochester.rr.com
    Release: 2
    Release date: 4/15/06
    """
    # Get array size
    n, m = A.shape
    A = A.flatten()
    nm = n * m

    # Find NaN elements
    k = np.isnan(A)

    # List the nodes which are known, and which will be interpolated
    nan_list = np.where(k)[0]
    known_list = np.where(~k)[0]

    # How many nans overall
    nan_count = len(nan_list)

    # If there are no NaNs, just return the original array
    if nan_count == 0:
        return A.reshape(n, m)

    # Convert NaN indices to (r,c) form
    nr, nc = np.unravel_index(nan_list, (n, m))

    # Both forms of index in one array:
    # column 0 == unrolled index
    # column 1 == row index
    # column 2 == column index
    nan_list = np.column_stack((nan_list, nr, nc))

    # Different methods
    if method == 0:
        # The same as method == 1, except only work on those
        # elements which are NaN, or at least touch a NaN.

        # Is it 1-d or 2-d?
        if (m == 1) or (n == 1):
            # Really a 1-d case
            work_list = nan_list[:, 0]
            work_list = np.unique(np.concatenate((work_list, work_list - 1, work_list + 1)))
            work_list = work_list[(work_list > 0) & (work_list < nm)]
            nw = len(work_list)

            # Build sparse matrix for del^2 operator
            fda = sparse.lil_matrix((nw, nm))
            for i in range(nw):
                fda[i, work_list[i]-1] = 1
                fda[i, work_list[i]] = -2
                fda[i, work_list[i]+1] = 1
        else:
            # A 2-d case

            # Horizontal and vertical neighbors only
            talks_to = np.array([[-1, 0], [0, -1], [1, 0], [0, 1]])
            neighbors_list = identify_neighbors(n, m, nan_list, talks_to)

            # List of all nodes we have identified
            all_list = np.vstack((nan_list, neighbors_list))

            # Generate sparse array with second partials on row
            # variable for each element in either list, but only
            # for those nodes which have a row index > 1 or < n
            L = np.where((all_list[:, 1] > 0) & (all_list[:, 1] < n-1))[0]
            nl = len(L)
            if nl > 0:
                fda = sparse.lil_matrix((nm, nm))
                for i in range(nl):
                    idx = all_list[L[i], 0]
                    fda[idx, [idx-1, idx, idx+1]] = [1, -2, 1]
            else:
                fda = sparse.lil_matrix((nm, nm))

            # 2nd partials on column index
            L = np.where((all_list[:, 2] > 0) & (all_list[:, 2] < m-1))[0]
            nl = len(L)
            if nl > 0:
                for i in range(nl):
                    idx = all_list[L[i], 0]
                    if idx - n >= 0:
                        fda[idx, idx-n] += 1
                    fda[idx, idx] += -2
                    if idx + n < nm:
                        fda[idx, idx+n] += 1

        # Eliminate knowns
        rhs = -fda[:, known_list].dot(A[known_list])
        k = np.where(np.any(fda[:, nan_list[:, 0]], axis=1))[0]

        # And solve...
        B = A.copy()
        B[nan_list[:, 0]] = spsolve(fda[k, :][:, nan_list[:, 0]], rhs[k])

    elif method == 1:
        # Least squares approach with del^2

        # Is it 1-d or 2-d?
        if (m == 1) or (n == 1):
            # A 1-d case
            fda = sparse.lil_matrix((nm-2, nm))
            for i in range(nm-2):
                fda[i, i:i+3] = [1, -2, 1]
        else:
            # A 2-d case

            # Compute finite difference for second partials on row variable first
            fda = sparse.lil_matrix((nm, nm))
            for i in range(1, n-1):
                for j in range(m):
                    idx = i + j*n
                    fda[idx, [idx-1, idx, idx+1]] = [1, -2, 1]

            # Now second partials on column variable
            for i in range(n):
                for j in range(1, m-1):
                    idx = i + j*n
                    fda[idx, [idx-n, idx, idx+n]] += [1, -2, 1]

        # Eliminate knowns
        rhs = -fda[:, known_list].dot(A[known_list])
        k = np.where(np.any(fda[:, nan_list[:, 0]], axis=1))[0]

        # And solve...
        B = A.copy()
        B[nan_list[:, 0]] = spsolve(fda[k, :][:, nan_list[:, 0]], rhs[k])

    elif method == 2:
        # Direct solve for del^2 BVP across holes

        # Is it 1-d or 2-d?
        if (m == 1) or (n == 1):
            # Really just a 1-d case
            raise ValueError('Method 2 has problems for vector input. Please use another method.')
        else:
            # A 2-d case
            L = np.where((nan_list[:, 1] > 0) & (nan_list[:, 1] < n-1))[0]
            nl = len(L)
            if nl > 0:
                fda = sparse.lil_matrix((nm, nm))
                for i in range(nl):
                    idx = nan_list[L[i], 0]
                    fda[idx, [idx-1, idx, idx+1]] = [1, -2, 1]
            else:
                fda = sparse.lil_matrix((nm, nm))

            # 2nd partials on column index
            L = np.where((nan_list[:, 2] > 0) & (nan_list[:, 2] < m-1))[0]
            nl = len(L)
            if nl > 0:
                for i in range(nl):
                    idx = nan_list[L[i], 0]
                    fda[idx, [idx-n, idx, idx+n]] += [1, -2, 1]

            # Fix boundary conditions at extreme corners
            # of the array in case there were nans there
            if 0 in nan_list[:, 0]:
                fda[0, 0] = -2
                fda[0, 1] = 1
                fda[0, n] = 1
            if n-1 in nan_list[:, 0]:
                fda[n-1, n-1] = -2
                fda[n-1, n-2] = 1
                fda[n-1, 2*n-1] = 1
            if nm-n in nan_list[:, 0]:
                fda[nm-n, nm-n] = -2
                fda[nm-n, nm-n+1] = 1
                fda[nm-n, nm-2*n] = 1
            if nm-1 in nan_list[:, 0]:
                fda[nm-1, nm-1] = -2
                fda[nm-1, nm-2] = 1
                fda[nm-1, nm-n-1] = 1

            # Eliminate knowns
            rhs = -fda[:, known_list].dot(A[known_list])

            # And solve...
            B = A.copy()
            k = nan_list[:, 0]
            B[k] = spsolve(fda[k, :][:, k], rhs[k])

    elif method == 3:
        # The same as method == 0, except uses del^4 as the
        # interpolating operator.

        # Del^4 template of neighbors
        talks_to = np.array([[-2, 0], [-1, -1], [-1, 0], [-1, 1], [0, -2], [0, -1],
                             [0, 1], [0, 2], [1, -1], [1, 0], [1, 1], [2, 0]])
        neighbors_list = identify_neighbors(n, m, nan_list, talks_to)

        # List of all nodes we have identified
        all_list = np.vstack((nan_list, neighbors_list))

        # Generate sparse array with del^4, but only
        # for those nodes which have a row & column index
        # >= 3 or <= n-2
        L = np.where((all_list[:, 1] >= 2) &
                     (all_list[:, 1] <= (n-3)) &
                     (all_list[:, 2] >= 2) &
                     (all_list[:, 2] <= (m-3)))[0]
        nl = len(L)
        if nl > 0:
            # Do the entire template at once
            fda = sparse.lil_matrix((nm, nm))
            for i in range(nl):
                idx = all_list[L[i], 0]
                # Apply the stencil values individually
                if idx - 2*n >= 0:
                    fda[idx, idx - 2*n] = 1
                if idx - n - 1 >= 0:
                    fda[idx, idx - n - 1] = 2
                if idx - n >= 0:
                    fda[idx, idx - n] = -8
                if idx - n + 1 >= 0:
                    fda[idx, idx - n + 1] = 2
                if idx - 2 >= 0:
                    fda[idx, idx - 2] = 1
                if idx - 1 >= 0:
                    fda[idx, idx - 1] = -8
                fda[idx, idx] = 20
                if idx + 1 < nm:
                    fda[idx, idx + 1] = -8
                if idx + 2 < nm:
                    fda[idx, idx + 2] = 1
                if idx + n - 1 < nm:
                    fda[idx, idx + n - 1] = 2
                if idx + n < nm:
                    fda[idx, idx + n] = -8
                if idx + n + 1 < nm:
                    fda[idx, idx + n + 1] = 2
                if idx + 2*n < nm:
                    fda[idx, idx + 2*n] = 1
        else:
            fda = sparse.lil_matrix((nm, nm))

        # On the boundaries, reduce the order around the edges
        L = np.where(((all_list[:, 1] == 1) | (all_list[:, 1] == (n-2))) &
                     (all_list[:, 2] >= 1) & (all_list[:, 2] <= (m-2)) |
                     ((all_list[:, 2] == 1) | (all_list[:, 2] == (m-2))) &
                     (all_list[:, 1] >= 1) & (all_list[:, 1] <= (n-2)))[0]
        nl = len(L)
        if nl > 0:
            for i in range(nl):
                idx = all_list[L[i], 0]
                # Apply the stencil values individually
                if idx - n >= 0:
                    fda[idx, idx - n] += 1
                if idx - 1 >= 0:
                    fda[idx, idx - 1] += 1
                fda[idx, idx] += -4
                if idx + 1 < nm:
                    fda[idx, idx + 1] += 1
                if idx + n < nm:
                    fda[idx, idx + n] += 1

        L = np.where(((all_list[:, 1] == 0) | (all_list[:, 1] == n-1)) &
                     (all_list[:, 2] >= 1) & (all_list[:, 2] <= (m-2)))[0]
        nl = len(L)
        if nl > 0:
            for i in range(nl):
                idx = all_list[L[i], 0]
                # Apply the stencil values individually
                if idx - n >= 0:
                    fda[idx, idx - n] += 1
                fda[idx, idx] += -2
                if idx + n < nm:
                    fda[idx, idx + n] += 1

        L = np.where(((all_list[:, 2] == 0) | (all_list[:, 2] == m-1)) &
                     (all_list[:, 1] >= 1) & (all_list[:, 1] <= (n-2)))[0]
        nl = len(L)
        if nl > 0:
            for i in range(nl):
                idx = all_list[L[i], 0]
                # Apply the stencil values individually
                if idx - 1 >= 0:
                    fda[idx, idx - 1] += 1
                fda[idx, idx] += -2
                if idx + 1 < nm:
                    fda[idx, idx + 1] += 1

        # Eliminate knowns
        rhs = -fda[:, known_list].dot(A[known_list])
        k = np.where(np.any(fda[:, nan_list[:, 0]], axis=1))[0]

        # And solve...
        B = A.copy()
        B[nan_list[:, 0]] = spsolve(fda[k, :][:, nan_list[:, 0]], rhs[k])

    elif method == 4:
        # Spring analogy

        # List of all springs between a node and a horizontal
        # or vertical neighbor
        hv_list = np.array([[-1, -1, 0], [1, 1, 0], [-n, 0, -1], [n, 0, 1]])
        hv_springs = []
        for i in range(4):
            hvs = nan_list + np.tile(hv_list[i], (nan_count, 1))
            k = (hvs[:, 1] >= 0) & (hvs[:, 1] < n) & (hvs[:, 2] >= 0) & (hvs[:, 2] < m)
            hv_springs.append(np.column_stack((nan_list[k, 0], hvs[k, 0])))

        hv_springs = np.vstack(hv_springs)

        # Delete replicate springs
        hv_springs = np.unique(np.sort(hv_springs, axis=1), axis=0)

        # Build sparse matrix of connections, springs
        nhv = hv_springs.shape[0]
        springs = sparse.lil_matrix((nhv, nm))
        for i in range(nhv):
            springs[i, hv_springs[i, 0]] = 1
            springs[i, hv_springs[i, 1]] = -1

        # Eliminate knowns
        rhs = -springs[:, known_list].dot(A[known_list])

        # And solve...
        B = A.copy()
        B[nan_list[:, 0]] = spsolve(springs[:, nan_list[:, 0]], rhs)

    elif method == 5:
        # Average of 8 nearest neighbors

        # Generate sparse array to average 8 nearest neighbors
        # for each nan element, be careful around edges
        fda = sparse.lil_matrix((nm, nm))

        # -1,-1
        L = np.where((nan_list[:, 1] > 0) & (nan_list[:, 2] > 0))[0]
        nl = len(L)
        if nl > 0:
            for i in range(nl):
                idx = nan_list[L[i], 0]
                fda[idx, idx-n-1] = 1
                fda[idx, idx] = -1

        # 0,-1
        L = np.where(nan_list[:, 2] > 0)[0]
        nl = len(L)
        if nl > 0:
            for i in range(nl):
                idx = nan_list[L[i], 0]
                fda[idx, idx-n] = 1
                fda[idx, idx] = -1

        # +1,-1
        L = np.where((nan_list[:, 1] < n-1) & (nan_list[:, 2] > 0))[0]
        nl = len(L)
        if nl > 0:
            for i in range(nl):
                idx = nan_list[L[i], 0]
                fda[idx, idx-n+1] = 1
                fda[idx, idx] = -1

        # -1,0
        L = np.where(nan_list[:, 1] > 0)[0]
        nl = len(L)
        if nl > 0:
            for i in range(nl):
                idx = nan_list[L[i], 0]
                fda[idx, idx-1] = 1
                fda[idx, idx] = -1

        # +1,0
        L = np.where(nan_list[:, 1] < n-1)[0]
        nl = len(L)
        if nl > 0:
            for i in range(nl):
                idx = nan_list[L[i], 0]
                fda[idx, idx+1] = 1
                fda[idx, idx] = -1

        # -1,+1
        L = np.where((nan_list[:, 1] > 0) & (nan_list[:, 2] < m-1))[0]
        nl = len(L)
        if nl > 0:
            for i in range(nl):
                idx = nan_list[L[i], 0]
                fda[idx, idx+n-1] = 1
                fda[idx, idx] = -1

        # 0,+1
        L = np.where(nan_list[:, 2] < m-1)[0]
        nl = len(L)
        if nl > 0:
            for i in range(nl):
                idx = nan_list[L[i], 0]
                fda[idx, idx+n] = 1
                fda[idx, idx] = -1

        # +1,+1
        L = np.where((nan_list[:, 1] < n-1) & (nan_list[:, 2] < m-1))[0]
        nl = len(L)
        if nl > 0:
            for i in range(nl):
                idx = nan_list[L[i], 0]
                fda[idx, idx+n+1] = 1
                fda[idx, idx] = -1

        # Eliminate knowns
        rhs = -fda[:, known_list].dot(A[known_list])

        # And solve...
        B = A.copy()
        k = nan_list[:, 0]
        B[k] = spsolve(fda[k, :][:, k], rhs[k])

    else:
        raise ValueError('If supplied, method must be one of: {0,1,2,3,4,5}.')

    # All done, make sure that B is the same shape as A was when we came in
    B = B.reshape(n, m)

    return B

def identify_neighbors(n, m, nan_list, talks_to):
    """
    Identify all the neighbors of those nodes in nan_list, not including the nans themselves

    Parameters:
    -----------
    n, m : int
        Size of the array
    nan_list : ndarray
        List of every nan element in the array
        nan_list[i,0] == linear index of i'th nan element
        nan_list[i,1] == row index of i'th nan element
        nan_list[i,2] == column index of i'th nan element
    talks_to : ndarray
        Defines which nodes communicate with each other, i.e., which nodes are neighbors
        talks_to[i,0] - defines the offset in the row dimension of a neighbor
        talks_to[i,1] - defines the offset in the column dimension of a neighbor

    Returns:
    --------
    neighbors_list : ndarray
        List of all neighbors of all the nodes in nan_list
    """
    if len(nan_list) == 0:
        return np.array([])

    # Use the definition of a neighbor in talks_to
    nan_count = nan_list.shape[0]
    talk_count = talks_to.shape[0]

    nn = np.zeros((nan_count * talk_count, 2), dtype=int)
    for i in range(talk_count):
        nn[i*nan_count:(i+1)*nan_count, :] = nan_list[:, 1:3] + np.tile(talks_to[i], (nan_count, 1))

    # Drop those nodes which fall outside the bounds of the original array
    L = (nn[:, 0] < 0) | (nn[:, 0] >= n) | (nn[:, 1] < 0) | (nn[:, 1] >= m)
    nn = nn[~L]

    # Form the same format 3 column array as nan_list
    neighbors_list = np.column_stack((np.ravel_multi_index((nn[:, 0], nn[:, 1]), (n, m)), nn))

    # Delete replicates in the neighbors list
    neighbors_list = np.unique(neighbors_list, axis=0)

    # And delete those which are also in the list of NaNs
    if len(neighbors_list) > 0:
        neighbors_list = np.array([x for x in neighbors_list if not np.any(np.all(x == nan_list, axis=1))])

    return neighbors_list

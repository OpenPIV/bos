import numpy as np

def write_openpiv_vec(filename, data, xUnits, tUnits, numrows, numcols):
    """
    Write OpenPIV vector data to a file in VEC format
    
    Parameters:
    -----------
    filename : str
        Full path to the output file
    data : ndarray
        Data to write, shape (N, 5) or (5, N)
    xUnits : str
        Units for spatial coordinates
    tUnits : str
        Units for time
    numrows : int
        Number of rows in the vector field
    numcols : int
        Number of columns in the vector field
        
    Returns:
    --------
    success : bool
        True if the file was written successfully, False otherwise
        
    Notes:
    ------
    The VEC format is a text file with a header and data columns.
    The data columns are: X, Y, U, V, CHC (correlation height)
    """
    try:
        with open(filename, 'w') as fid:
            # Check data shape and transpose if necessary
            if data.shape[1] == 5:
                data = data.T
            elif data.shape[0] != 5:
                raise ValueError('Wrong number of columns')
            
            # Write header
            zone = f'ZONE I={numrows}, J={numcols}'
            header = f'VARIABLES= "X {xUnits}", "Y {xUnits}", "U {xUnits}/{tUnits}", "V {xUnits}/{tUnits}", "CHC", {zone}\n'
            fid.write(header)
            
            # Write data
            for i in range(data.shape[1]):
                fid.write(f'{data[0, i]:3d} {data[1, i]:3d} {data[2, i]:7.4f} {data[3, i]:7.4f} {data[4, i]:7.4f}\n')
        
        return True
    except:
        return False

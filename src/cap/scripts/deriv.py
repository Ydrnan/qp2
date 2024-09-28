import math
import numpy as np
import sys

def deriv(x,y):
    dx = x[1] - x[0]
    df = np.zeros((len(y)))

    df[0] = 1.0/(2.0*dx) * (-3.0 * y[0] + 4.0 * y[1] - y[2])
    for i in range(1,len(y)-1):
        df[i] = 1.0/(2.0 * dx) * (-y[i-1] + y[i+1])
    df[-1] = 1.0/(2.0*dx) * (y[-3] - 4.0 * y[-2] + 3.0 * y[-1])

    return df

def deriv_complex(x,y_re,y_im):
    d_re = deriv(x,y_re)
    d_im = deriv(x,y_im)

    for i in range(len(x)):
        print("{:8.3e} {:12.6f} {:12.6f} {:12.4f}".format(x[i], x[i] * d_re[i], x[i] * d_im[i], x[i] * abs(math.sqrt(d_re[i]**2 + d_im[i]**2))))


def read_text_file(file_path):
    with open(file_path,'r') as f:
        lines = f.readlines()
    data = []
    for line in lines:
        if len(line.strip()) != 0: 
            data.append([float(i) for i in line.split()])

    return data

def main():
    # Check if the correct number of command-line arguments is provided
    if len(sys.argv) != 2:
        print("Usage: python script.py <file_path>")
        sys.exit(1)

    # Get the file path from the command-line argument
    file_path = sys.argv[1]

    # Read the text file and create the NumPy array
    data = read_text_file(file_path)
    data = np.asarray(data, dtype=float)

    # Print the resulting array
    print(data)
    deriv_complex(data[:,0], data[:,1], data[:,2])

if __name__ == "__main__":
    main()

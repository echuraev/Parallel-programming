nx = 201 

def write_matrix(file_name, matrix):
    two_dim_array = False
    with open(file_name, 'wt') as f:
        for i in range(len(matrix)):
			if i == 0:
				if type(matrix[i]) is list:
					two_dim_array = True
			if two_dim_array is True:
				for j in range(len(matrix[i])):
					if j == len(matrix[i])-1:
						f.write(str(matrix[i][j]) + '\n')
					else:
						f.write(str(matrix[i][j]) + ';')
			else:
				if i == len(matrix)-1:
					f.write(str(matrix[i]) + '\n')
				else:
					f.write(str(matrix[i]) + ';')

def gen_ro():
    ro = [1 if i < nx/2.0 else 0 for i in range(0, nx)]
    return ro

if __name__ == '__main__':
    write_matrix('ro.init', gen_ro())
    

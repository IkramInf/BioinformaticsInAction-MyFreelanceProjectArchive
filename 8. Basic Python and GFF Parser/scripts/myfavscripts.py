# create a function named validateMatrix to confirm the size and the values of same type
def validateMatrix(matrix, ncol=6, nrow=4, num=True):
    # assign row and column value of the given matrix
    row, col = len(matrix), len(matrix[0])

    # check the size
    if row == nrow and col == ncol:  
        print("Matrix size is OK.")
    else:
        print("Matrix size is not OK.")

    # check the type of all values
    # if all values are numbers, all() will return True
    if all([isinstance(c, (int, float)) for r in matrix for c in r]) == num:
        print("All Matrix values are numeric.")
    # if all values are characters, all() will return True
    elif all([isinstance(c, str) for r in matrix for c in r]):
        print("All Matrix values are characters.")
    else:
        print("Matrix contains different types of values.")
            

# create a function to find out length of variable            
def len(data):
    count = 0
    # iterating over given data
    for i in data:
        count += 1
    return count


if __name__ == "__main__":
    matrix = [[1, 2, 3, 4, 5, 6],
              [11, 12, 13, 14, 15, 16],
              [21, 22, 23, 24, 25, 26],
              [31, 32, 33, 34, 35, 36]]
    
    # testing validateMatrix function
    validateMatrix(matrix, col, row)

    # testing len function
    print("Length of row:", len(matrix))
            
            

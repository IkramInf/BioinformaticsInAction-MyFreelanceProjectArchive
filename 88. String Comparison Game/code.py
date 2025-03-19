# function to compare between two strings
def strNcompare(s1, s2, n):
    s1_portion = s1[0:n].lower()
    s2_portion = s2[0:n].lower()

    if s1_portion == s2_portion:
        return 0    # s1 is equal to s2
    elif s1_portion > s2_portion:
        return -1   # s1 is greater than s2
    else:
        return 1    # s1 is less than s2

# keep running the loop until getting -1 from user
while True:
    flag = int(input("String comparison [1(play), -1(quit)]: ", ))
    
    if flag == -1:
        break
    else:
        s1 = input("Enter the first string: ", )
        s2 = input("Enter the second string: ", )
        n = int(input("Number of characters: ", ))
        
        # call the function strNcompare
        output = strNcompare(s1, s2, n)
        
        if output == 0:
            print(f"{s1} is equal to {s2}\n")
        elif output == -1:
            print(f"{s1} is greater than {s2}\n")
        else:
            print(f"{s1} is less than {s2}\n")
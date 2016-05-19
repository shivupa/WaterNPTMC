import numpy as np
a= 3
"""
0   1   2   3   4
    5   6   7   8
        9   10  11
            12  13
                14

0   1   2   3   4
5   6   7   8   9
10  11  12  13  14
15  16  17  18  19
20  21  22  23  24

0,0 0,1 0,2 0,3 0,4
1,0 1,1 1,2 1,3 1,4
2,0 2,1 2,2 2,3 2,4
3,0 3,1 3,2 3,3 3,4
4,0 4,1 4,2 4,3 4,4
"""
count = 3
print count
for i in range(4,2,-1):
    count += i
    print count
for i in range(5-a):
    print a*(a+1) + i

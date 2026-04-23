
def R(i,j,k):
    if k == 0:
        if i== 0:
            if j==0:
                return ""
            elif j==1:
                return "1"
            else:
                return "2"
        elif i==1:
            if j==0:
                return "2"
            elif j==1:
                return ""
            else:
                return "1"
        else:
            if j==0:
                return "1"
            elif j==1:
                return "2"
            else:
                return ""
    temp = "("+ R(k, k, k-1) + ")*" if R(k, k, k-1) != "" else ""
    temp1 = R(i,j,k-1) + "|" if R(i,j,k-1) != "" else ""
    return temp1+ "("+ R(i,k,k-1) + temp +R(k, j, k-1) + ")"

print(R(0,0,3))

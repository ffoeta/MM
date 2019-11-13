import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root
from scipy import linalg
from texttable import Texttable


p4_ = 2
p5_ = 2
p9_ = 22
p8_ = 1000

x2root = 0

columnp = 0
columnl1 = 5
columnl2 = 6
columnl3 = 7
columnl4 = 8

array = np.ones((1000))
res = np.ones((100,10))

def drawOne(n,x,y,l1,l2,l3,l4):

    plt.grid(True)
    plt.title('X' + str((n+1)))

    # for i in range(np.prod(x.shape)):
    #     if ( (l1[i].imag == 0) and (l2[i].imag == 0) and (l3[i].imag == 0) and (l4[i].imag == 0) ):
    #         if ( (l1[i] < 0) and (l2[i] < 0) and (l3[i] < 0) and (l4[i] < 0) ) :
    #             color = 'green' 
    #             marker = '.'
    #         else :
    #             color = 'red'
    #             marker = 'o'
    #     else :
    #         color = 'black'
    #         marker = 'x'

    for i in range(np.prod(x.shape)):
        if ( (l1[i].real < 0) and (l2[i].real < 0) and (l3[i].real < 0) and (l4[i].real < 0) ):
            color = 'green' 
            marker = '.'
            plt.plot(x[i], y[i], color = color, marker = marker, markersize = 2.0)
        else :
            if ( (l1[i].imag != 0) and (l2[i].imag != 0) or (l3[i].imag != 0) and (l4[i].imag != 0) ):
                color = 'red'
                marker = '+'
                plt.plot(x[i], y[i], color = color, marker = marker, markersize = 2.0, zorder=10)
            else :
                color = 'blue'
                marker = '_'
                plt.plot(x[i], y[i], color = color, marker = marker, markersize = 2.0, zorder=11 )


        

    

    return 0

def drawFour(array):

    axes = [221,222,223,224]

    p = array[:,columnp].real
    l1 = array[:,columnl1]
    l2 = array[:,columnl2]
    l3 = array[:,columnl3]
    l4 = array[:,columnl4]

    for i in range(4):
        x = array[:,i+1].real
        plt.subplot(axes[i])
        drawOne(i,p,x,l1,l2,l3,l4)


    # p = array[:,columnp].real
    # l1 = array[:,columnl1]
    # l2 = array[:,columnl2]
    # l3 = array[:,columnl3]
    # l4 = array[:,columnl4]
    # x = array[:,1].real
    # drawOne(1,p,x,l1,l2,l3,l4)

    return 0

def fprint(array):

    t = Texttable()
    t.add_rows([['P', 'X1', 'X2', 'X3', 'X4','L1','L2','L3','L4']])
    t.set_cols_align(["l","r","r","r","r","r","r","r","r"])
    t.set_cols_width([6,6,6,6,6,16,16,16,16])
    t.set_precision(5)
    file = open('text.txt', 'w')
    
    for row in array:
        t.add_row([str(row[0].real),str(row[1].real),str(row[2].real),str(row[3].real),str(row[4].real),str(row[5]),str(row[6]),str(row[7]),str(row[8])])
    file.write(t.draw())
    file.close()

    return 0

def x1(x2):
    return (3*x2/22)

def x3(x2,x4):
    return ((2*x2 + 3*x4)/22)

def p(x2):
    return 3*x2/((22 - 3*x2)*np.exp(x2/(1 + x2/1000)))

def f(x2,x4):
    return x2 - 3*x4 + 22*(3*x2/((22-3*x2)*np.exp(x2/(1+x2/1000))))*(1-(2*x2+3*x4)/22)*np.exp(x4/(1+x4/1000))

def froot(x4):
    return x2root - 3*x4 + 22*(3*x2root/((22-3*x2root)*np.exp(x2root/(1+x2root/1000))))*(1-(2*x2root+3*x4)/22)*np.exp(x4/(1+x4/1000))


def yacobi(values):

    def exp1(x):
        return np.exp(x/(1 + x/p8_))

    def exp2(x):
        return (np.exp((p8_*x)/(p8_+x))/(p8_+x)**2)

    x1_ = values[0]
    x2_ = values[1]
    x3_ = values[2]
    x4_ = values[3]
    p_  = values[4]

    return np.array([[  (-1-p_*exp1(x2_)),      (p_*(p8_**2)*(1-x1_)*exp2(x2_)),            
                        0,                      0                                           ],
                    [   (-p_*p9_*exp1(x2_)),    (-1-p4_+p_*(p8_**2)*p9_*(1-x1_)*exp2(x2_)),  
                        0,                      0                                           ],
                    [   0,                      0,                                          
                        (-1-p_*exp1(x4_)),      (p_*(p8_**2)*(1-x3_)*exp2(x4_))             ],
                    [   0,                      0,                                          
                        (-p_*p9_*exp1(x4_)),    (-1-p5_+p_*(p8_**2)*p9_*(1-x3_)*exp2(x4_)) ]])

x2res = 0
x4res = 0

results = np.array([[0,0,0,0,0,0,0,0,0]])

for j in range(800):
    a = 0
    b = 0
    for i in range(1000):

        a = b
        b = f(j/100,i/100)

        if (a*b < 0):
            x2root = j/100
            temp = root(froot, [(i-1)/100, i/100])

            _x1_ = x1(x2root)
            _x2_ = x2root
            _x3_ = x3(x2root,temp.x[0])
            _x4_ = temp.x[0]
            _p_  = p(x2root)

            yac = yacobi([  _x1_,
                            _x2_,
                            _x3_,
                            _x4_,
                            _p_ ])
            yac = np.linalg.eig(yac)[0]
            results = np.vstack((results, [[    _p_ ,
                                                _x1_,
                                                _x2_,
                                                _x3_,
                                                _x4_,
                                                yac[0],
                                                yac[1],
                                                yac[2],
                                                yac[3],
                                                ]]))

results = np.delete(results, 0, 0)


fprint(results)
drawFour(results)

plt.show()

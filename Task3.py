def leitor(inputFile='T3input.txt'):
    input = open(inputFile, encoding='utf-8').readlines()
    global ICOD, n, x, y, x_target
    x = []
    y = []
    try:
        ICOD = int(input[1])
        n = int(input[3])
        x_target = float(input[-1])
        print(f'ICOD = {ICOD}\nn = {n}\n')
        for linha in input[6:6+n]: # Lê os pontos do input
            xi, yi = linha.split(' ')
            x.append(float(xi)), y.append(float(yi))
    except:
        with open('T3output.txt', 'w', encoding='utf-8') as file:
            print('O arquivo de input não foi devidamente escrito. Para mais informações, leia as instruções de uso em: "linkGitHUb"')
            file.write('O arquivo de input não foi devidamente escrito. Para mais informações, leia as instruções de uso em: "linkGitHUb"')


def Interpolate(x, y, x_target: int):
    resultado = 0 # resultado final (soma dos y_i*phi_i) inicializado em 0
    for i in range(len(y)):
        parcial = y[i] # Valor do termo inicializado como y_i
        for k in range(len(x)):   #OBS len(x) = len(y)
            if k != i:
                parcial *= ((x_target - x[k])/(x[i]-x[k]))
        resultado += parcial # soma termo i ao resultado final
    return resultado


def Regression(x, y, x_target: int): 
    a11 = len(x) # a11 = N
    a12 = a21 = a22 = c11 = c21 = 0
    for i in range(len(x)):
        a12 += x[i]
        a22 += x[i]*x[i]
        c11 += y[i]
        c21 += x[i]*y[i]
    a21 = a12 # matriz A eh simetrica
    # A=[[a11, a12],[a21, a22]] e C=[[c11],[c21]]
    detA = a11*a22 - a12*a21 # Diagonal principal - Diagonal secundaria
    # OBS: Matriz adjunta 2x2: Troca posisao na diag. prin. e troca sinal diag. sec.
    invA = [[a22/detA, -a12/detA],[-a21/detA, a11/detA]] # Matriz adjunta / determinante
    b0= invA[0][0]*c11+invA[0][1]*c21 # B = A^-1 * C
    b1= invA[1][0]*c11+invA[1][1]*c21 # B = A^-1 * C
    print(f"y = {b0} + {b1}x", ) if b1 > 0 else print(f"y = {b0} {b1}x", )
    return b0 + b1*x_target


leitor()
with open('T3output.txt', 'w', encoding='utf-8') as file:
    if ICOD == 1:
        result = Interpolate(x, y, x_target)
        print(f"Pelo método de Lagrange de interpolação, temos:\nValor estimado de y para x={x_target} é {result}")
        file.write(f"Pelo método de Lagrange de interpolação, temos:\nValor estimado de y para x={x_target} é {result}")
    elif ICOD == 2:
        result = Regression(x, y, x_target)
        print(f"Por Regressão, temos:\nValor estimado de y para x={x_target} é {result}")
        file.write(f"Por Regressão, temos:\nValor estimado de y para x={x_target} é {result}")
    else:
        print("Por favor, configure ICOD=1 ou ICOD=2.\nPara mais informações, leia as instruções de uso em: https://github.com/FTPaiva/AlgebraLinearTrab1 ")
        file.write("Por favor, configure ICOD=1 ou ICOD=2.\nPara mais informações, leia as instruções de uso em: https://github.com/FTPaiva/AlgebraLinearTrab1 ")
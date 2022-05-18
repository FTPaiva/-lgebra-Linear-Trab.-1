def LU(a, b):
    for k in range(n-1): #  "n-1" pq a ultima coluna n precisa ser zerada
        
        # Loop abaixo define os elementos da matriz L (apenas os de baixo da diagonal principal)
        for i in range(k+1, len(a)): # "n" pq tem q passar por todas as linhas abaixo do pivo
            try:
                a[i][k] /= a[k][k] # Posições que seriam zeradas são utilizadas pra guardar valores da matriz L
            except:
                print("Não é possível realizar a decomposição LU para matrizes com zero na diagonal principal!\nPor favor utilize outro algoritmo.")
                return 1
            # No Loop acima, o M[i][k] sempre sera um elemento abaixo da diagonal principal, ja que i>k
            
        # Loop abaixo define os elementos da matriz U (da diagonal principal para cima)
        for j in range(k+1, len(a)):
            for i in range(k+1, len(a)):
                a[i][j] -= a[i][k]*a[k][j]
    return (BackwardSubstitution(a, LUForwardSubstitution(a, b)))


def LUForwardSubstitution(A, B): # Vamos usar matriz A para obter os dados da matriz L
    Y = [B[0]]
    for i in range(1, len(B)): # Itera por todas as linhas a partir da segunda
        temp = B[i] # Novo elemento que sera calculado
        for j in range(0, i): # Itera pelas colunas e para antes de "chegar" no elemento da diagonal principal
            temp -= A[i][j]*Y[j]
        Y.append(temp) # Novo elemento adicionado ao vetor Y
    return Y


def BackwardSubstitution(A, Y):
    size = len(A)
    X = [0 for i in range(size)] # Cria o vetor X cheio de zeros
    X[size-1] = Y[size-1]/A[size-1][size-1] # Define o último elemento de X
    for i in range(size-1, -1, -1): # Percorre as linhas de A de baixo para cima
        temp = Y[i] # Novo elemento que será calculado
        for j in range(size-1, i, -1): # Percorre as colunas de A da direita para a esquerda e para quando "chega" no elemento da diagonal principal
            temp -= A[i][j]*X[j] 
        temp /= A[i][i]
        X[i] = temp # Novo elemento adicionado ao vetor X
    return X


def leitor(inputFile='T3input.txt', funcFile=''):
    input = open(inputFile, encoding='utf-8').readlines()
    global ICOD, n, x, y, x_target, funcs
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
    if funcFile !='':
        expressionsFile=open(funcFile, encoding='utf-8').readlines()
        funcs = []
        for linha in expressionsFile:
            funcs.append(linha.strip('\n'))


def Interpolate(x, y, x_target):
    resultado = 0 # resultado final (soma dos y_i*phi_i) inicializado em 0
    for i in range(len(y)):
        parcial = y[i] # Valor do termo inicializado como y_i
        for k in range(len(x)):   #OBS len(x) = len(y)
            if k != i:
                parcial *= ((x_target - x[k])/(x[i]-x[k]))
        resultado += parcial # soma termo i ao resultado final
    return resultado


def Regression(x, y, x_target): 
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


def MultilinearRegression(x, y, x_target, funcs):
    P, Y = [], []
    for i in range(n): # For each point
        Y.append([y[i]]) # Create Y vector
        P.append([]) #Add the lines
        for func in funcs: # For each given expression
            P[i].append(CalculateExpression(x[i], func))
    A = multiplyMatrix(P, P, 1) # A = Pt * P
    C = multiplyMatrix(P, Y, 1) # C = Pt * Y
    C_alterado = []
    for line in C:
        C_alterado.append(line[0]) #Funcao LU requer outro formato
    B = LU(A, C_alterado) # original AX = B, nesse caso AB = C (descobrimos B)
    result = 0
    for i in range(len(funcs)):
        result += B[i]* CalculateExpression(x_target, funcs[i])
    return result


def CalculateExpression(val, func):
    tempstr = ''
    for char in func:
        if char == 'x':
            tempstr += f'{val}'
        elif char == 'e':
            tempstr += f'2.71828'
        elif char == 'π':
            tempstr += f'3.14159'
        else:
            tempstr += char
    try:
        return eval(tempstr)
    except:
        print("Expressao dada é inválida")


def multiplyMatrix(M1, M2, transp=0):
    result = []
    for i in range(len(M1[0])): # Iterando pelas linhas de M1 (n1)
        result.append([])
        for j in range(len(M2[0])): # Iterando pelas colunas de M2 (m2)
            result[i].append(0)
            for k in range(len(M2)): # Iterando pelas linhas de M2 (n2)
                if transp == 1: # Troquei k por i na M1 (considerar M1 transposta)
                    result[i][j] += M1[k][i] * M2[k][j]
                else: # Anda uma col na M1 e uma linha na M2
                    result[i][j] += M1[i][k] * M2[k][j]
    return result


leitor(funcFile='T3funcs.txt')
with open('T3output.txt', 'w', encoding='utf-8') as file:
    if ICOD == 1:
        result = Interpolate(x, y, x_target)
        print(f"Pelo método de Lagrange de interpolação, temos:\nValor estimado de y para x={x_target} é {result}")
        file.write(f"Pelo método de Lagrange de interpolação, temos:\nValor estimado de y para x={x_target} é {result}")
    elif ICOD == 2:
        result = MultilinearRegression(x, y, x_target, funcs)
        print(f"Por Regressão Multilinear, temos:\nValor estimado de y para x={x_target} é {result}")
        file.write(f"Por Regressão Multilinear, temos:\nValor estimado de y para x={x_target} é {result}")
    elif ICOD == 3:
        result = Regression(x, y, x_target)
        print(f"Por Regressão Linear, temos:\nValor estimado de y para x={x_target} é {result}")
        file.write(f"Por Regressão Linear, temos:\nValor estimado de y para x={x_target} é {result}")
    else:
        print("Por favor, configure ICOD=1 ou ICOD=2.\nPara mais informações, leia as instruções de uso em: https://github.com/FTPaiva/AlgebraLinearTrab1 ")
        file.write("Por favor, configure ICOD=1 ou ICOD=2.\nPara mais informações, leia as instruções de uso em: https://github.com/FTPaiva/AlgebraLinearTrab1 ")


''' 
def transpose(M):
    for i in range(len(M)):
        for j in range(len(M[0])):
            if i > j:
                temp = M[i][j]
                M[i][j] = M[j][i]
                M[j][i] = temp
    return M
'''
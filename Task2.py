import math # Biblioteca math utilizada para funções como arctg, cos e sen
def leitor(inputFile='T2input.txt',):
    input = open(inputFile, encoding='utf-8').readlines()
    inputA = open("T2A.txt", encoding='utf-8').readlines()
    global ICOD, n, IDET, TOLm, det, A, numIter
    A = []

    try:
        n = int(input[1]) # Lê os parâmetros de input do arquivo txt
        ICOD = int(input[3])
        IDET = int(input[5])
        TOLm = float(input[7])

        for linha in inputA: # Lê o arquivo de input da matriz A
            A.append([float(x) for x in linha.split(' ')])

        print(f'n = {n}\nICOD = {ICOD}')
    except:
        with open('T2output.txt', 'w', encoding='utf-8') as file:
            print('Os arquivos de input não foram devidamente escritos. Para mais informações, leia as instruções de uso em: https://github.com/FTPaiva/AlgebraLinearTrab1')
            file.write('Os arquivos de input não foram devidamente escritos. Para mais informações, leia as instruções de uso em: https://github.com/FTPaiva/AlgebraLinearTrab1')


def MultiplyAX(A, X): # Função acessória para o método de potência para aplicar a multiplicação de matrizes A*X
    global n
    result = []
    for i in range(n): # Iterando pelas linhas de A
        result.append(0)
        for j in range(n): # Iterando pelos elementos de X
            result[i] += A[i][j]*X[j]

    return result

def multiplyMatrix(M1, M2, transp=0): # Função que implementa a multiplicação de matrizes
    result = []
    for i in range(n): # Iterando pelas linhas de M1 (n1)
        result.append([])
        for j in range(n): # Iterando pelas colunas de M2 (m2)
            result[i].append(0)
            for k in range(n): # Iterando pelas linhas de M2 (n2)
                if transp == 1: # Troquei k por i na M1 (considerar M1 transposta)
                    result[i][j] += M1[k][i] * M2[k][j]
                else: # Anda uma col na M1 e uma linha na M2
                    result[i][j] += M1[i][k] * M2[k][j]
    return result

def PowerMethod():
    try:
        global TOLm, A, numIter
        X = [1 for i in range(n)] # Vetor X composto por elementos 1
        lambdaOld = 1 # lambda0 = 1 
        lambdaNew = 0
        residuo = float('inf')

        numIter = 0
        while abs(residuo) >= TOLm:

            X = MultiplyAX(A, X) # Calculando Y e armazenando em X

            lambdaNew = X[0]
            for i in range(n): # Atualizando lambda = maior elemento (valor absoluto) de Y
                if (abs(X[i]) > abs(lambdaNew)):
                    lambdaNew = abs(X[i])

            for k in range(n): # Dividindo os termos de X por lambda
                X[k] = X[k]/lambdaNew

            residuo = (lambdaNew - lambdaOld)/lambdaNew # Calculando resíduo
            lambdaOld = lambdaNew
            numIter += 1

        lambdaNew = round(lambdaNew, 3) # Aproximando os valores do autovalor e autovetor
        for i in range(n):
            X[i] = round(X[i], 3)
    except:
        return(0, 0, True)
    return (lambdaNew, X, False) 


def Jacobi():
    global TOLm, det, A, numIter
    I = []
    for i in range(n): # Criando template de matriz identidade para a criação da matriz P
        I.append([0 for x in range(n)])
        I[i][i] = 1
    X = [linha[:] for linha in I] # X inicialmente identidade, deepcopy de I


    for i in range(n): # Verificando se a matriz A é simétrica
        for j in range(i, n):
            if(A[i][j] != A[j][i]):
                return(0, 0, True) # Caso contrário, para o método e acusa erro

    
    
    maiorI = 0 # Índices do maior elemento fora da diagonal principal de A
    maiorJ = 1

    for j in range(1, n): # Procura o maior elemento na primeira linha 
        if (A[0][j] >  A[0][maiorJ]):
            maiorI = 0
            maiorJ = j

    
    phi = 0 # Inicializando parâmetro phi da matriz P
    numIter = 0 # Número de iterações
    while True:

        for i in range(n): # Encontrar maior elemento de A fora da diagonal principal (valor absoluto)
            for j in range(n):
                if (i != j and abs(A[i][j]) > abs(A[maiorI][maiorJ])):
                    maiorI = i
                    maiorJ = j

        if (A[maiorI][maiorI] != A[maiorJ][maiorJ]): # Define o valor do parâmetro phi
            phi = math.atan((2*A[maiorI][maiorJ]/(A[maiorI][maiorI] - A[maiorJ][maiorJ]))/2)
        else:
            phi = math.pi/4

        P = [linha[:] for linha in I] # Define a matriz P a partir da identidade
        P[maiorI][maiorI], P[maiorJ][maiorJ] = math.cos(phi), math.cos(phi) # Elementos Pii e Pjj iguais a cos(phi)
        P[maiorI][maiorJ], P[maiorJ][maiorI] = -math.sin(phi), math.sin(phi) # Elementos Pij e Pji iguais a
        # -sen(phi) e sen(phi), respectivamente

        A = multiplyMatrix(P, A, 1) # Calculando nova matriz A (P transposta * A * P)
        A = multiplyMatrix(A, P)
        
        X = multiplyMatrix(X, P) # Calculando nova matriz X (X * P)

        repeat = False
        for i in range(n): # Testando se todos os elementos fora da diagonal de A são menores que a tolerância
            for j in range(n):
                if (i != j and abs(A[i][j]) > TOLm): repeat = True
        
        numIter += 1

        if (not repeat): break

    det = 1 # Calculando o determinante da matriz A original
    for i in range(n): # Multiplicando todos os autovalores da matriz
        det *= A[i][i]
    det = round(det, 3)

    for i in range(n): # Aproximando os valores das matrizes A e X
        for j in range(n):
            if (i != j): A[i][j] = 0 # Tornando 0 os valores fora da diagonal principal de A
            A[i][j] = round(A[i][j], 3)
            X[i][j] = round(X[i][j], 3)

    return (A, X, False)

leitor()
with open('T2output.txt', 'w', encoding='utf-8') as file:
    if ICOD == 1:
        result = PowerMethod()
        if (result[2]):
            print("A matriz inserida não é válida para o método de potência!")
            file.write("A matriz inserida não é válida para o método de potência!")
        else:
            print("Pelo método de potência, temos:"+
            f"\nAutovalor estimado = {result[0]}"+
            f"\nAutovetor estimado para o autovalor encontrado = {result[1]}T"+
            f"\nO número de iterações necessárias para obter os resultados foi {numIter}")
            file.write("Pelo método de potência, temos:"+
            f"\nAutovalor estimado = {result[0]}"+
            f"\nAutovetor estimado para o autovalor encontrado = {result[1]}T"+
            f"\nO número de iterações necessárias para obter os resultados foi {numIter}")
    elif ICOD == 2:
        result = Jacobi()
        if (result[2]):
            print("A matriz inserida não é válida para o método de potência, pois não é simétrica!")
            file.write("A matriz inserida não é válida para o método de potência, pois não é simétrica!")
        else:
            print("Pelo método de Jacobi, temos:"+
            f"\nMatriz dos autovalores estimados = {result[0]}"+
            f"\nMatriz dos autovetores estimados para os autovalores encontrados = {result[1]}"+
            f"\nO número de iterações necessárias para obter os resultados foi {numIter}")
            file.write("Pelo método de Jacobi, temos:"+
            f"\nMatriz dos autovalores estimados = {result[0]}"+
            f"\nMatriz dos autovetores estimados para os autovalores encontrados = {result[1]}"+
            f"\nO número de iterações necessárias para obter os resultados foi {numIter}")
    else:
        print("Por favor, configure ICOD=1 ou ICOD=2.\nPara mais informações, leia as instruções de uso em: https://github.com/FTPaiva/AlgebraLinearTrab1")
        file.write("Por favor, configure ICOD=1 ou ICOD=2.\nPara mais informações, leia as instruções de uso em: https://github.com/FTPaiva/AlgebraLinearTrab1")
    if (IDET > 0 and not result[2]):
        if (ICOD == 2):
            print(f"O determinante encontrado para a matriz foi: {det}")
            file.write(f"\nO determinante encontrado para a matriz foi: {det}")
        else: 
            print("O método utilizado não criou caminhos para calcularmos o determinante! " +
            "Por favor, utilize outro algoritmo.")
            file.write("\nO método utilizado não criou caminhos para calcularmos o determinante! " +
            "Por favor, utilize outro algoritmo.")
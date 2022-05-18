def leitor(inputFile='T2input.txt',):
    input = open(inputFile, encoding='utf-8').readlines()
    global ICOD, n, IDET, TOLm

    try:
        n = int(input[1])
        ICOD = int(input[3])
        IDET = int(input[5])
        TOLm = float(input[7])

        print(f'n = {n}\nICOD = {ICOD}\n')
    except:
        with open('T2output.txt', 'w', encoding='utf-8') as file:
            print('O arquivo de input não foi devidamente escrito. Para mais informações, leia as instruções de uso em: https://github.com/FTPaiva/AlgebraLinearTrab1')
            file.write('O arquivo de input não foi devidamente escrito. Para mais informações, leia as instruções de uso em: https://github.com/FTPaiva/AlgebraLinearTrab1')


def MultiplyAX(A, X):
    global n
    result = []
    for i in range(n): # Iterando pelas linhas de A
        result.append(0)
        for j in range(n): # Iterando pelos elementos de X
            result[i] += A[i][j]*X[j]

    return result

def PowerMethod():
    global TOLm
    A = [[1.0, 0.2, 0.0], [0.2, 1.0, 0.5], [0.0, 0.5, 1.0]]
    X = [1, 1, 1]
    lambdaOld = 1
    lambdaNew = 0
    residuo = float('inf')

    while abs(residuo) >= TOLm:
        X = MultiplyAX(A, X) # Calculando Y e armazenando em X
        lambdaNew = X[0] # Atualizando lambda = Y1
        X[0] = 1 # Atualizando X1 = 1
        for i in range(1, len(X)): # Dividindo os termos de X a partir de X2 por Y1 (lambda)
            X[i] = X[i]/lambdaNew
        residuo = (lambdaNew - lambdaOld)/lambdaNew # Calculando resíduo
        lambdaOld = lambdaNew 
    
    # CALCULAR DET
    return (lambdaNew, X)


def Jacobi():
    1


leitor()
with open('T2output.txt', 'w', encoding='utf-8') as file:
    if ICOD == 1:
        result = PowerMethod()
        print("Pelo método de potência, temos:"+
        f"\nAutovalor estimado = {result[0]}"+
        f"\nAutovetor estimado para o autovalor encontrado = {result[1]}T")
        file.write("Pelo método de potência, temos:"+
        f"\nAutovalor estimado = {result[0]}"+
        f"\nAutovetor estimado para o autovalor encontrado = {result[1]}T")
    elif ICOD == 2:
        1
    else:
        print("Por favor, configure ICOD=1 ou ICOD=2.\nPara mais informações, leia as instruções de uso em: https://github.com/FTPaiva/AlgebraLinearTrab1")
        file.write("Por favor, configure ICOD=1 ou ICOD=2.\nPara mais informações, leia as instruções de uso em: https://github.com/FTPaiva/AlgebraLinearTrab1")




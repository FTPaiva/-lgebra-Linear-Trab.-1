det = 1

def leitor(inputFile='T1input.txt'):
    config = open(inputFile, encoding='utf-8').readlines()
    inputA = open("T1A.txt", encoding='utf-8').readlines()
    inputB = open("T1B.txt", encoding='utf-8').readlines()
    global A, B, n, ICOD, IDET, TOLm
    A = []
    B = []
    try:
        n = int(config[0][2:])
        ICOD = int(config[1][5:])
        IDET = int(config[2][5:])
        TOLm = float(config[3][5:])
        print(f'n = {n}\nICOD = {ICOD}\nIDET = {IDET}\nTOLm = {TOLm}\n')
        for linha in inputA: # Lê o arquivo de input da matriz A
            A.append([int(x) for x in linha.split(' ')])
        for linha in inputB: # Lê o arquivo de input do vetor B
            B.append(int(linha))
    except:
        print('O arquivo de input não foi devidamente escrito. Por favor, leia as instruções de uso em: "linkGitHUb" ')


def LUForwardSubstitution(A, B): # Vamos usar matriz A para obter os dados da matriz L
    Y = [B[0]]
    for i in range(1, n): # Itera por todas as linhas a partir da segunda
        temp = B[i] # Novo elemento que sera calculado
        for j in range(0, i): # Itera pelas colunas e para antes de "chegar" no elemento da diagonal principal
            temp -= A[i][j]*Y[j]
        Y.append(temp) # Novo elemento adicionado ao vetor Y
    return Y

def CholeskyForwardSubstitution(A, B): # Vamos usar matriz A para obter os dados da matriz L
    Y = [0 for i in range(n)]
    Y[0] = B[0]/A[0][0]
    for i in range(1, n): # Itera por todas as linhas a partir da segunda
        temp = B[i] # Novo elemento que sera calculado
        for j in range(0, i+1): # Itera pelas colunas e para antes de "chegar" no elemento da diagonal principal
            temp -= A[i][j]*Y[j]
        temp /= A[i][i]
        Y[i] = temp # Novo elemento adicionado ao vetor Y
    return Y
    
        
def BackwardSubstitution(A, Y):
    X = [0 for i in range(n)] # Cria o vetor X cheio de zeros
    X[n-1] = Y[n-1]/A[n-1][n-1] # Define o último elemento de X
    for i in range(n-1, -1, -1): # Percorre as linhas de A de baixo para cima
        temp = Y[i] # Novo elemento que será calculado
        for j in range(n-1, i, -1): # Percorre as colunas de A da direita para a esquerda e para quando "chega" no elemento da diagonal principal
            temp -= A[i][j]*X[j] 
        temp /= A[i][i]
        X[i] = temp # Novo elemento adicionado ao vetor X
    return X


def LU():
    
    for k in range(n-1): #  "n-1" pq a ultima coluna n precisa ser zerada
        
        # Loop abaixo define os elementos da matriz L (apenas os de baixo da diagonal principal)
        for i in range(k+1, n): # "n" pq tem q passar por todas as linhas abaixo do pivo
            try:
                A[i][k] /= A[k][k] # Posições que seriam zeradas são utilizadas pra guardar valores da matriz L
            except:
                print("Não é possível realizar a decomposição LU para matrizes com zero na diagonal principal!\nPor favor utilize outro algoritmo.")
                return 1
            # No Loop acima, o M[i][k] sempre sera um elemento abaixo da diagonal principal, ja que i>k
            
        # Loop abaixo define os elementos da matriz U (da diagonal principal para cima)
        for j in range(k+1, n):
            for i in range(k+1, n):
                A[i][j] -= A[i][k]*A[k][j]
    
    global det
    for k in range(n): # Itera por todos elementos da diagonal principal
        det *= A[k][k]
    if det == 0:
        print("Não é possível realizar a decomposição LU para matrizes singulares!\nPor favor utilize outro algoritmo.")
        return 2
    else:
        return (BackwardSubstitution(A, LUForwardSubstitution(A, B)))


def Cholesky():

    for i in range(n):
        temp = 0
        for k in range(i):
            temp += (A[i][k])**2
            if (A[i][k] != A[i][k]):
                print("Não é possível realizar a decomposição de Cholesky para matrizes assimétricas!\nPor favor utilize outro algoritmo.")
                return 1
        if (A[i][i] - temp > 0):
            A[i][i] = (A[i][i] - temp)**0.5
        else: 
            print("Só é possível realizar a decomposição de Cholesky para matrizes positivas definidas!\nPor favor utilize outro algoritmo.")
            return 2

        for j in range(i + 1, n):
            temp = 0
            for k in range(i-1):
                temp += A[i][k] * A[j][k]
            A[j][i] = (A[i][j] - temp)/A[i][i]
            A[i][j] = A[j][i]
        
        global det
        
    for k in range(n): # Itera por todos elementos da diagonal principal
        det *= A[k][k] * A[k][k]

    return (BackwardSubstitution(A, CholeskyForwardSubstitution(A, B)))

def Jacobi():
    OldX = [1 for i in range(n)] # Vetor solução inicial composto por valores unitários
    NewX = OldX.copy()
    residuo = float('inf') # resíduo inicial infinito para garantir primeira execução do while
    numIter = 0 # Contador de iteracoes necessarias
    HistResiduo = []
    while residuo >= TOLm:
        for i in range(n):
            temp = B[i]
            testeDiagonalDominante = 0
            for j in range(n):            
                if i != j:
                    temp -= A[i][j] * OldX[j]
                    testeDiagonalDominante += abs(A[i][j]) # soma os valores absolutos da linha, execto o pivo
            if (abs(A[i][i]) < testeDiagonalDominante): # Testa se matriz eh diagonal dominante
                print("Não é possível realizar o método de Jacobi para matrizes que não são diagonal dominantes!\nPor favor utilize outro algoritmo.")
                return 1, 1, 1
            try:
                temp /= A[i][i]
                NewX[i] = temp # Salva o termo calculado no vetor X mais recente
            except:
                print("A matriz inserida não é válida para o método de Jacobi, pois existe zero na diagonal principal.\nPor favor utilize outro algoritmo.")
                return 2, 2, 2

        normaDif, normaNovo = 0, 0 # Normas calculadas para obter o novo residuo
        for i in range(n):
            try:
                normaDif += (NewX[i] - OldX[i])**2
                normaNovo += NewX[i]**2
            except:
                print("A matriz inserida não é válida para o método de Jacobi!\nPor favor utilize outro algoritmo.")
                return 3, 3, 3

        normaDif = normaDif**(0.5)
        normaNovo = normaNovo**(0.5)
        residuo = normaDif/normaNovo # Atualiza o residuo
        HistResiduo.append(residuo)
        
        numIter += 1
        OldX = NewX.copy() # Atualiza o vetor x antigo para ser igual ao mais novo

    return NewX, HistResiduo, numIter

def GaussSeidel():
    OldX = [1 for i in range(n)] # Vetor solução inicial composto por valores unitários
    NewX = OldX.copy()
    residuo = float('inf') # resíduo inicial infinito para garantir primeira execução do while
    numIter = 0 # Contador de iteracoes necessarias
    HistResiduo = []
    while residuo >= TOLm:
        for i in range(n):
            temp = B[i]
            testeDiagonalDominante = 0
            for j in range(n):
                if i != j:
                    temp -= A[i][j] * NewX[j] # Valores do NewX ate i-1 estao atualizados, resto ainda ta igual ao OldX
            try:
                temp /= A[i][i]
                NewX[i] = temp # Salva o termo calculado no vetor X mais recente
            except:
                print("A matriz inserida não é válida para o método de Gauss-Seidel, pois existe zero na diagonal principal.\nPor favor utilize outro método.")
                return 1, 1, 1
        
        normaDif, normaNovo = 0, 0 # Normas calculadas para obter o novo residuo
        for i in range(n):
            try: # Se ocorrer erro, valores não convergem!
                normaDif += (NewX[i] - OldX[i])**2
                normaNovo += NewX[i]**2
            except:
                print("A matriz inserida não é válida para o método de GaussSeidel!")
                return 2, 2, 2

        normaDif = normaDif**(0.5)
        normaNovo = normaNovo**(0.5)
        residuo = normaDif/normaNovo # Atualiza o residuo
        HistResiduo.append(residuo)
        
        numIter += 1
        OldX = NewX.copy() # Atualiza o vetor x antigo para ser igual ao mais novo

    return NewX, HistResiduo, numIter

def PrintVector(X, n):
    for elemento in X[:-1]:
        print(f"%.{n}f"%elemento, end=", ")
    print(f"%.{n}f ]"%X[-1])

def WriteVector(X, n):
    for elemento in X[:-1]:
        string = f"%.{n}f, "%elemento
        file.write(string)
    string = f"%.{n}f]"%X[-1]
    file.write(string)
    

leitor()


with open('T1output.txt', 'w', encoding='utf-8') as file:

    if ICOD == 1:
        X = LU()
        if type(X)!= int:
            file.write("Solução do sistema: X = [")
            WriteVector(X, 3)
        elif X == 1:
            file.write("Não é possível realizar a decomposição LU para matrizes com zero na diagonal principal!" +
             "\nPor favor utilize outro algoritmo.")
        elif X == 2:
            file.write("Não é possível realizar a decomposição LU para matrizes singulares!" +
            "\nPor favor utilize outro algoritmo.")
    elif ICOD == 2:
        X = Cholesky()
        if type(X)!= int:
            file.write("Solução do sistema: X = [")
            WriteVector(X, 3)
        elif X == 1:
            file.write("Não é possível realizar a decomposição de Cholesky para matrizes assimétricas!" +
            "\nPor favor utilize outro algoritmo.")
        elif X == 2:
            file.write("Só é possível realizar a decomposição de Cholesky para matrizes positivas definidas!" + 
            "\nPor favor utilize outro algoritmo.")
    elif ICOD == 3:
        X, HistResiduo, numIter = Jacobi()
        if type(X)!= int:
            file.write("Solução do sistema: X = [")
            WriteVector(X, 3)
            str = f"\nNúmero de iterações para convergência: {numIter}"
            file.write(str)
            file.write("\nHistórico do resíduo: X = [")
            WriteVector(HistResiduo, 5)
        elif X == 1:
            file.write("Não é possível realizar o método de Jacobi para matrizes que não são diagonal dominantes!"
            + "\nPor favor utilize outro algoritmo.")
        elif X == 2:
            file.write("A matriz inserida não é válida para o método de Jacobi, pois existe zero na diagonal principal." +
            "\nPor favor utilize outro algoritmo.")
        elif X == 3:
            file.write("A matriz inserida não é válida para o método de Jacobi!" +
            "\nPor favor utilize outro algoritmo.")
    elif ICOD == 4:
        X, HistResiduo, numIter = GaussSeidel()
        if type(X)!= int:
            file.write("Solução do sistema: X = [")
            WriteVector(X, 3)
            str = f"\nNúmero de iterações para convergência: {numIter}"
            file.write(str)
            file.write("\nHistórico do resíduo: X = [")
            WriteVector(HistResiduo, 5)
        elif X == 1:
            file.write("A matriz inserida não é válida para o método de Gauss-Seidel, pois existe zero na diagonal principal." + 
            "\nPor favor utilize outro método.")
        elif X == 2:
            file.write("A matriz inserida não é válida para o método de GaussSeidel!")
    else:
        file.write("O valor de ICOD fornecido não é válido! Para mais informações, leia as instruções de uso em: https://github.com/FTPaiva/AlgebraLinearTrab1")
        raise Exception("O valor de ICOD fornecido não é válido! Para mais informações, leia as instruções de uso em: https://github.com/FTPaiva/AlgebraLinearTrab1")

    if (type(X)!= int) and (IDET > 0):
        if ICOD == 1 or ICOD == 2:
            str = f"\nDeterminante: %.4f"%det
            file.write(str)
        else:
            file.write("\nA técnica de obter a solução do sistema não criou caminhos para calcularmos o determinante! " +
             "Por favor, utilize outro algoritmo.")


# Álgebra Linear Computacional - Trabalho 1
Primeiro trabalho para avaliação da disciplina Álgebra Linear Computacional COC473.

Programa em Python para solução de problemas de Álgebra Linear.

## Task 01:

### O usuário pode escolher entre os métodos:

1. Decomposição LU (ICOD = 1);
2. Decomposição de Cholesky (ICOD = 2);
3. Procedimento iterativo Jacobi (ICOD = 3);
4. Procedimento iterativo Gauss-Seidel (ICOD = 4);

### INPUTS (arquivo de entrada):

1. Ordem n do sistema de equações;
2. ICOD relativo ao método de análise;
3. IDET;
4. Matriz A;
5. Vetor B;
6. TOLm (Tolerância máxima utilizada para os métodos iterativos);

### OUTPUTS (arquivo de saída):
1. Solução X do sistema;
2. Avisos sobre erros de uso (quando presentes);
3. Determinante (quando solicitado);
4. Número de iterações e histórico do TOL no caso dos métodos iterativos;

## Task 02:

### O usuário pode escolher entre os métodos:

1. Decomposição LU (ICOD = 1);
2. Decomposição de Cholesky (ICOD = 2);

### INPUTS (arquivo de entrada):

1. Ordem n do sistema de equações;
2. ICOD relativo ao método de análise;
3. IDET;
4. Matriz A;
5. Vetor B;
6. TOLm (Tolerância máxima utilizada para os métodos iterativos);

### OUTPUTS (arquivo de saída):
1. Solução X do sistema;
2. Avisos sobre erros de uso (quando presentes);
3. Determinante (quando solicitado);
4. Número de iterações e histórico do TOL no caso dos métodos iterativos;

## Task 03:

### O usuário pode escolher entre os métodos:

1. Interpolação (ICOD = 1);
2. Regressão Multilinear (ICOD = 2);
3. Regressão (ICOD = 3);

### INPUTS (arquivo de entrada):

1. ICOD;
2. Número de pontos conhecidos;
3. Coordenadas x e y dos pontos;
4. Coordenada x para estimar y;
5. Funções, caso utilize a Regressão Multilinear;

### OUTPUTS (arquivo de saída):
1. Valor estimado de y para o x dado;
2. Avisos sobre erros de uso (quando presentes);

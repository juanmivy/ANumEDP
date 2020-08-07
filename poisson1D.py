#El objetivo es desarrollar un programa que resuelva numericamente el siguiente problema de valores iniciales conocido como la ecuación de Poisson en 1D
#u''(x)=f(x) con x en [a,b] y u(a) , u(b) valores conocidos

import numpy as np
import math
import matplotlib.pyplot as plt

#Primero definimos la malla en el intervalo

a = float (input('Introduce el valor inicial del intervalo: '))
b = float (input('Introduce el valor final del intervalo: '))
nodos = int (input('Introduce el número de nodos deseado: '))
paso=(b-a)/((nodos-1) * 1.0)

x = np.linspace(a,b,nodos)

#Introducimos la función deseada
def f(x):
    return np.sin(x)

#Introducimos los valores de la función u en los extremos

u_a = float (input('Introduce el valor de u en a: '))
u_b = float (input('Introduce el valor de u en b: '))

#Hay que crear una matriz de tamaño nxn

#matriz = np.zeros((nodos-2,nodos-2))
#print(matriz)

#creamos un vector de unos de tamaño nodos - 2
vector_unos=[]
for i in range(nodos-3): #como este vector va en la diagonal superior e inferior ha de tener un elemento menos que la diagonal principal
    vector_unos = vector_unos +  [1]

vector_menos_dos = []
for i in range(nodos-2):
    vector_menos_dos += [-2]
    
matriz_metodo = np.diag(vector_menos_dos) + np.diag(vector_unos,1) + np.diag(vector_unos,-1)  
matriz_metodo = matriz_metodo*(1/paso**2)

valores_funcion = f(x[1:nodos-1])
#ahora hay que corregir el primer valor y el ultimo

valores_funcion[0] = valores_funcion[0] - u_a/(paso**2 * 1.0) 
valores_funcion[nodos-3] = valores_funcion[nodos-3] - u_b/(paso**2 * 1.0) 

#Ahora solo queda resolver el sistema de ecuaciones anterior y añadirle los valores de los extremos

u_solucion =  np.linalg.solve(matriz_metodo, valores_funcion) 

u_solucion =u_solucion.tolist() #pasamos a lista para que sea mas facil trabajar con el 

u_solucion = [u_a] + u_solucion + [u_b]
#Ahora podemos sacar los datos en una grafica


# Plot the data
plt.plot(x, u_solucion, label='u(x)')

# Add a legend
plt.legend()

# Show the plot
plt.show()




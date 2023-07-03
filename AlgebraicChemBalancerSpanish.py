#!/usr/bin/env python3
#Author: Arnoldo Fernando Chue Sánchez
#Contact: arnoldwork20@gmail.com
#License: GNU/GPL 

import numpy as np
import sympy

aclaracionParaInputs = '''\nIndicaciones para los inputs:
    Los subíndices de los elementos van en tamaño normal después del elemento en el que aparecen y deben ser estrictamente menores a 10.
    Dejar espacio entre cada compuesto y símbolo de la reacción.
    Poner el símbolo => para indicar que se realiza la reacción.
Ejemplo:
    HCl + Na3PO4 => H3PO4 + NaCl\n'''

class reaccion:
    def __init__(self, reaccion:str):
        self.reaccion = reaccion
        self.react = []
        self.prod = []
        self.elem = set()
        self.cantidad = {}

    def reactivos(self):
        self.prod = self.reaccion.strip().split("=>")
        self.prod = self.prod[0].split("+")
        self.prod = [x.strip() for x in self.prod]
        return self.prod
    
    def productos(self):
        self.react = self.reaccion.strip().split("=>")
        self.react = self.react[1].split("+")
        self.react = [x.strip() for x in self.react]
        return self.react
    
    def elementos(self):
        x = ""
        for e in self.reaccion:
            if e == "=" or e == ">": break
            elif e == "+": pass
            elif e == " " or e.isnumeric(): 
                self.elem.add(x)
                x = ""
            elif e.isupper():
                self.elem.add(x)
                x = e
            else: x += e
        self.elem.add(x)
        self.elem.discard("")
        return self.elem
    
    def poner_1_explicitos(self):
        cont = 0
        for i, e in enumerate(self.reaccion):
            i += cont
            if i == 0 or e == "+" or e == ">": continue
            if e.isupper() and (self.reaccion[i-1].isupper() or self.reaccion[i-1].islower()): 
                self.reaccion = self.reaccion[:i]+"1"+self.reaccion[i:]
                cont += 1
                continue
            if e == " " and (self.reaccion[i-1].isupper() or self.reaccion[i-1].islower()): 
                self.reaccion = self.reaccion[:i]+"1"+self.reaccion[i:]
                cont += 1
                continue
        if not(self.reaccion[-1].isnumeric()): self.reaccion += "1"
        return self.reaccion
    
    def cantidad_elementos_por_molecula(self, reaccion_1_explicitos:str):
        e, num, ladoizquierdo, moleculas = "", "", True, 1
        for x in reaccion_1_explicitos:
            if x.isnumeric(): num += x
            elif x.isalpha():
                if num != "": 
                    if ladoizquierdo: 
                        try: self.cantidad[e].append(int(num))
                        except: self.cantidad[e] = [int(num)]
                    else: 
                        try: self.cantidad[e].append(-int(num))
                        except: self.cantidad[e] = [-int(num)]
                    num, e = "", x
                else: e += x
            elif x == ">": ladoizquierdo = False
            elif x == "+" or x == "=":
                if ladoizquierdo: 
                    try: self.cantidad[e].append(int(num))
                    except: self.cantidad[e] = [int(num)]
                else: 
                    try: self.cantidad[e].append(-int(num))
                    except: self.cantidad[e] = [-int(num)]
                num, e = "", ""
                for y in self.elem:
                    try: 
                        if len(self.cantidad[y]) < moleculas: 
                            try: self.cantidad[y].append(0)
                            except: self.cantidad[y] = [0]
                    except: self.cantidad[y] = [0]
                moleculas += 1
        try: self.cantidad[e].append(-int(num))
        except: self.cantidad[e] = [-int(num)]
        for y in self.cantidad:
            if len(self.cantidad[y]) < moleculas: 
                try: self.cantidad[y].append(0)
                except: self.cantidad[y] = [0] 
        return self.cantidad
    
    def esta_balanceada_inicialmente(self):
        for x in self.cantidad:
            if sum(self.cantidad[x]) != 0: return False
        return True
    
if __name__ == '__main__':
    print(aclaracionParaInputs)
    reaccion = reaccion(input("Ingrese la reacción:\n"))

    reactivos = reaccion.reactivos()
    productos = reaccion.productos()
    elementos = reaccion.elementos()
    reaccion_explicita = reaccion.poner_1_explicitos()
    cantidad = reaccion.cantidad_elementos_por_molecula(reaccion_explicita)

    if reaccion.esta_balanceada_inicialmente(): print("\nEstá balanceada desde el inicio")
    else:
        M = np.matrix(len(elementos)*[([0]*(len(productos)+len(reactivos)))])
        for i, c in enumerate(cantidad): M[i] = cantidad[c].copy()
        
        f_e_r = sympy.Matrix(M)
        f_e_r, pivotes = f_e_r.rref()

        espacio_nulo_fraccionario = sympy.Matrix(f_e_r.nullspace())
        denominadores = espacio_nulo_fraccionario.applyfunc(lambda x: x.as_numer_denom()[1])
        denominadores = set(denominadores) 
        if len(denominadores) == 1: 
            denominadores = list(denominadores) 
            denominadores.append(denominadores[0])
            comun_denominador = sympy.lcm(denominadores)
        else: comun_denominador = sympy.lcm(*set(denominadores))
        coeficientes_enteros = (comun_denominador * espacio_nulo_fraccionario).applyfunc(lambda x: int(x))

        reaccion_balanceada = ""
        moleculas = reactivos + productos
        for i, m in enumerate(moleculas):
            reaccion_balanceada += str(coeficientes_enteros[i]) + m + " "
            if i == len(reactivos) - 1: reaccion_balanceada += "=> "
            else: reaccion_balanceada += "+ "
        print("\nLa ecuación balanceada es:")
        print(reaccion_balanceada[:-2])
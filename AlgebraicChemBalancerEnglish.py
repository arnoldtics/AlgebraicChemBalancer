#!/usr/bin/env python3
#Author: Arnoldo Fernando Chue Sánchez
#Contact: arnoldwork20@gmail.com
#License: GNU/GPL 

import numpy as np
import sympy

clarificationForInputs = '''\nInstructions for inputs:
    Subscripts for elements come after the element itself and must be strictly less than 10.
    Leave a space between each compound and the reaction symbol.
    Use the symbol => to indicate a reaction.
Example:
    HCl + Na3PO4 => H3PO4 + NaCl\n'''

class Reaction:
    def __init__(self, reaction: str):
        self.reaction = reaction
        self.reactants = []
        self.products = []
        self.elements = set()
        self.amounts = {}

    def get_reactants(self):
        self.products = self.reaction.strip().split("=>")
        self.products = self.products[0].split("+")
        self.products = [x.strip() for x in self.products]
        return self.products

    def get_products(self):
        self.reactants = self.reaction.strip().split("=>")
        self.reactants = self.reactants[1].split("+")
        self.reactants = [x.strip() for x in self.reactants]
        return self.reactants

    def get_elements(self):
        x = ""
        for e in self.reaction:
            if e == "=" or e == ">": break
            elif e == "+": pass
            elif e == " " or e.isnumeric():
                self.elements.add(x)
                x = ""
            elif e.isupper():
                self.elements.add(x)
                x = e
            else: x += e
        self.elements.add(x)
        self.elements.discard("")
        return self.elements

    def add_explicit_ones(self):
        count = 0
        for i, e in enumerate(self.reaction):
            i += count
            if i == 0 or e == "+" or e == ">": continue
            if e.isupper() and (self.reaction[i-1].isupper() or self.reaction[i-1].islower()):
                self.reaction = self.reaction[:i] + "1" + self.reaction[i:]
                count += 1
                continue
            if e == " " and (self.reaction[i-1].isupper() or self.reaction[i-1].islower()):
                self.reaction = self.reaction[:i] + "1" + self.reaction[i:]
                count += 1
                continue
        if not(self.reaction[-1].isnumeric()): self.reaction += "1"
        return self.reaction

    def get_element_amounts_per_molecule(self, explicit_reaction: str):
        e, num, left_side, molecules = "", "", True, 1
        for x in explicit_reaction:
            if x.isnumeric(): num += x
            elif x.isalpha():
                if num != "":
                    if left_side:
                        try: self.amounts[e].append(int(num))
                        except: self.amounts[e] = [int(num)]
                    else:
                        try: self.amounts[e].append(-int(num))
                        except: self.amounts[e] = [-int(num)]
                    num, e = "", x
                else: e += x
            elif x == ">": left_side = False
            elif x == "+" or x == "=":
                if left_side:
                    try: self.amounts[e].append(int(num))
                    except: self.amounts[e] = [int(num)]
                else:
                    try: self.amounts[e].append(-int(num))
                    except: self.amounts[e] = [-int(num)]
                num, e = "", ""
                for y in self.elements:
                    try:
                        if len(self.amounts[y]) < molecules:
                            try: self.amounts[y].append(0)
                            except: self.amounts[y] = [0]
                    except:self.amounts[y] = [0]
                molecules += 1
        try: self.amounts[e].append(-int(num))
        except: self.amounts[e] = [-int(num)]
        for y in self.amounts:
            if len(self.amounts[y]) < molecules:
                try: self.amounts[y].append(0)
                except: self.amounts[y] = [0]
        return self.amounts

    def is_initially_balanced(self):
        for x in self.amounts:
            if sum(self.amounts[x]) != 0: return False
        return True

if __name__ == '__main__':
    print(clarificationForInputs)
    input_reaction = Reaction(input("Enter the reaction:\n"))

    reactants = input_reaction.get_reactants()
    products = input_reaction.get_products()
    elements = input_reaction.get_elements()
    explicit_reaction = input_reaction.add_explicit_ones()
    amounts = input_reaction.get_element_amounts_per_molecule(explicit_reaction)

    if input_reaction.is_initially_balanced(): print("\nThe reaction is already balanced.")
    else:
        M = np.matrix(len(elements) * [([0] * (len(products) + len(reactants)))])
        for i, c in enumerate(amounts): M[i] = amounts[c].copy()

        f_e_r = sympy.Matrix(M)
        f_e_r, pivots = f_e_r.rref()

        null_space = sympy.Matrix(f_e_r.nullspace())
        denominators = null_space.applyfunc(lambda x: x.as_numer_denom()[1])
        denominators = set(denominators)
        if len(denominators) == 1:
            denominators = list(denominators)
            denominators.append(denominators[0])
            common_denominator = sympy.lcm(denominators)
        else: common_denominator = sympy.lcm(*set(denominators))
        integer_coefficients = (common_denominator * null_space).applyfunc(lambda x: int(x))

        balanced_reaction = ""
        molecules = reactants + products
        for i, m in enumerate(molecules):
            balanced_reaction += str(integer_coefficients[i]) + m + " "
            if i == len(reactants) - 1: balanced_reaction += "=> "
            else: balanced_reaction += "+ "
        print("\nThe balanced equation is:")
        print(balanced_reaction[:-2])
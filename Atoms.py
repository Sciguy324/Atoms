# To whom it may concern: avert your eyes!  Hide your children!  This code was written quite poorly...  Not my best script...
# However, if you have no idea what any of this means, or how to code for that matter, you should be fine.
import itertools

class Element:
    """Class for basic elements"""

    elements = {}

    def __init__(self, symbol, molar_mass, atomic_number, neutrons=None):
        self.symbol = symbol
        self.molar_mass = molar_mass
        self.atomic_number = atomic_number
        self.electrons = atomic_number
        if neutrons:
            self.neutrons = neutrons
        else:
            self.neutrons = round(molar_mass) - atomic_number
        Element.elements[symbol] = self

    def __repr__(self):
        """Return the element's symbol for string representation"""
        return str(self.symbol)

    def grams_to_moles(self, grams): # Tested logically
        """Converts grams of an element to moles"""
        return grams / self.molar_mass

    def g2mol(self, grams):
        """Alias for grams_to_moles"""
        return self.grams_to_moles(grams)

    def moles_to_grams(self, moles):
        """Converts moles of an element to grams"""
        return self.molar_mass * moles

    @classmethod
    def list(cls):     # This is organized/placed poorly
        """Print a list of all registered elements"""
        for i, j in cls.elements.items():
            print("{}{} - Atomic #: {},{} Molar Mass: {}".format(i,
                                                                 abs(3-len(i))*" ",
                                                                 j.atomic_number,
                                                                 (3 - len(str(j.atomic_number)))*" ",
                                                                 j.molar_mass))

    @staticmethod
    def moles_to_atoms(moles): # Tested implicitly, explicitly & logically
        """Converts moles of an element to grams"""
        return moles * 6.022e23

    @staticmethod
    def atoms_to_moles(atoms): # Tested logically
        """Converts atoms of an element to moles"""
        return atoms / 6.022e23
    
    def grams_to_atoms(self, grams): # Tested implicitly logically (Subfunctions should work, therefore this should too.)
        """Convert grams of an element to atoms"""
        return Element.moles_to_atoms(self.grams_to_moles(grams))

    def atoms_to_grams(self, atoms):
        """Converts atoms of an element to grams"""
        return self.moles_to_grams(Element.atoms_to_moles(atoms))

    def __mul__(self, x):
        if type(x) == Element:
            # Element * Element -> Molecule
            return Molecule([self, x])
        elif type(x) == Molecule:
            # Element * Molecule -> Molecule
            return Molecule([self] + x.atom_list * x.count)
        elif type(x) == int:
            # Element * n -> Molecule of 'n' Elements
            return Molecule([self for i in range(x)])
        else:
            # Invalid operation
            raise TypeError("Unsupported operation between {} and {}".format(type(self), type(x)))

    def __rmul__(self, x):
        if type(x) == int:
            # n * Element -> 'n' Molecules of monatomic Element
            return Molecule([self], count=x)
        else:
            # Support for multiplication... BUT IN REVERSE!!!
            return self.__mul__(x)

    def __add__(self, x):
        if type(x) == Element:
            # Element + Element -> Combination
            return Combination([Molecule([self]), Molecule([x])])
        elif type(x) == Molecule:
            # Element + Molecule -> Combination
            return Combination([Molecule([self]), x])
        elif type(x) == Combination:
            # Element + Combination -> Combination
            return Combination([self] + x.parts)
        else:
            # Invalid operation
            raise TypeError("Unsupported operation between {} and {}".format(type(self), type(x)))


class Molecule:
    """Class for molecules"""

    sub = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

    def __init__(self, atom_list, count=1):
        #self.atom_list = sorted(atom_list, key=lambda x: x.atomic_number)
        self.atom_list = atom_list
        self.count = count
        self.molar_mass = count * sum([i.molar_mass for i in atom_list])
        self.__repr__()

    def __repr__(self):
        """Return a string representation of the molecule.  Does not understand polyatomic ions."""
        symbol_dict = {}
        molecule_list = []
        for i in self.atom_list:
            if i.symbol not in symbol_dict:
                symbol_dict[i.symbol] = 1
                if type(i) == Molecule and i.count > 1:
                    molecule_list.append(i.symbol)
            else:
                symbol_dict[i.symbol] += 1
        if self.count == 1:
            symbol = ""
        else:
            symbol = str(self.count)
        for i, j in symbol_dict.items():
            if i in molecule_list:
                symbol += "("
            symbol += i
            if j > 1:
                symbol += str(j).translate(Molecule.sub)
            if i in molecule_list:
                symbol += ")"
        self.symbol = symbol
        return symbol

    def no_num(self):
        """Print the molecule's symbol without the leading coefficient"""
        symbol_dict = {}
        molecule_list = []
        symbol = ""
        for i in self.atom_list:
            if i.symbol not in symbol_dict:
                symbol_dict[i.symbol] = 1
                if type(i) == Molecule and i.count > 1:
                    molecule_list.append(i.symbol)
            else:
                symbol_dict[i.symbol] += 1
        for i, j in symbol_dict.items():
            if i in molecule_list:
                symbol += "("
            symbol += i
            if j > 1:
                symbol += str(j).translate(Molecule.sub)
            if i in molecule_list:
                symbol += ")"
        return symbol

    def copy(self):
        """Returns a new copy of the molecule"""
        return Molecule([i for i in self.atom_list], count=self.count)

    def set_count(self, new_count):
        """Sets the count of a molecule"""
        self.count = new_count
        self.molar_mass = new_count * sum([i.molar_mass for i in self.atom_list])
        return self

    def get_atom_count(self, element):
        """Count how many atoms of a particular element are in a molecule"""
        atom_dict = dict((i, 0) for i in self.atom_list)
        for i in self.atom_list:
            atom_dict[i] += 1
        if element not in atom_dict:
            return 0
        return atom_dict[element]

    def singular(self):
        """Returns a singular copy of a molecule"""
        return Molecule(self.atom_list, 1)

    def remove_poly(self, polyatomic, count=None):
        """Returns a copy of a molecule with a set of atoms removed"""
        if type(polyatomic) == Element:
            polyatomic = Molecule([polyatomic])
        elif type(polyatomic) == Molecule:
            pass
        else:
            raise TypeError("Invalid polyatomic type '{}'".format(type(polyatomic)))
        
        new_molecule = self.copy()

        if count:
            if count < 0:
                raise ValueError("'{}' is not a valid value for count".format(count))
        else:
            count = polyatomic.count
        
        if self.count != 1:
            raise ValueError("Function is not responsible for figuring out what to do when molecule count != 1")
        else:
            for i in polyatomic.atom_list * count:
                if i in new_molecule.atom_list:
                    new_molecule.atom_list.pop(new_molecule.atom_list.index(i))
                else:
                    raise ValueError("Not enough elements in molecule to complete this action")

        return new_molecule

    def grams_to_moles(self, grams, element=None): # Tested implicitly
        """Converts grams of a molecule to moles.  If an element is specified, the number of
        moles of a particular element will be isolated and returned"""
        if element:
            return len([i for i in self.atom_list if i == element]) * grams / self.molar_mass
        else:
            return grams / self.molar_mass

    def g2mol(self, grams):
        """Alias for grams_to_moles"""
        return self.grams_to_moles(grams)
    
    def grams_to_molecules(self, grams, element=None): # Tested
        """Converts grams of a molecule to the number of molecules."""
        return Element.moles_to_atoms(self.grams_to_moles(grams, element))

    def grams_to_atoms(self, grams, element):   # Tested
        """Computes how many atoms of a particular element are present in the given sample"""
        return self.grams_to_molecules(grams, element)

    def moles_to_grams(self, moles, element=None): # Tested implicitly
        """Converts the number of moles of a molecule to grams.  If an element is specified,
        the mass of that element will be isolated and returned instead."""
        if element:
            return moles * (Molecule([i for i in self.atom_list if i == element])).molar_mass
        else:
            return moles * self.molar_mass

    def molecules_to_moles(self, molecules, element=None): # Tested
        """Converts molecules of a sample to moles.  If an element is specified, the number of
        moles of that element will be isolated and returned instead."""
        if element:
            return len([i for i in self.atom_list if i == element]) * molecules / 6.022e23
        else:
            return molecules / 6.022e23

    def moles_to_molecules(self, moles):
        """Converts a specified number of moles to molecules"""
        return moles * 6.022e23

    def molecules_to_grams(self, molecules, element=None): # Tested
        """Converts molecules of a sample to grams.  If an element is specified, the mass of that
        element will be isolated and returned instead."""
        return self.moles_to_grams(self.molecules_to_moles(molecules), element)

    def _mass_to_empirical(self, mass_dict):
        """Takes the percent mass of a set of elements and determines the emprical
        formula."""
        new_mass_dict = {}
        for i, j in mass_dict.items():
            new_mass_dict[i] = i.grams_to_moles(j)

        mass1 = list(new_mass_dict.items())[0][1]
        #print(new_mass_dict)
        for i, j in new_mass_dict.items():
            new_mass_dict[i] /= mass1
            #print(i, j, j / mass1)

        for i in range(1, 20):
            guess = dict((j, k*i) for j, k in new_mass_dict.items())
            done_guessing = True
            for j, k in guess.items():
                if abs((round(k, 0) - k) / k) > 0.1:
                    done_guessing = False
                    break
            if done_guessing: break

        if done_guessing:
            formatted_guess = []
            for i, j in guess.items():
                formatted_guess += [i,] * round(j)
            print("Guess: ", Molecule(formatted_guess))

        return new_mass_dict

    def mass_to_empirical(self):
        """Wrapper for the _mass_to_empirical function"""
        mass_dict = {}
        for i in self.atom_list:
            mass_dict[i] = float(input("{}: ".format(i)))
        return self._mass_to_empirical(mass_dict)

    def _percents(self):
        """Get the percent mass of each element in a molecule"""
        atom_dict = {}
        for i in self.atom_list:
            if not i in atom_dict:
                atom_dict[i] = self.get_atom_count(i) * i.molar_mass / self.molar_mass * 100

        return atom_dict

    def percents(self):
        """Wrapper for _percents to output data in more readable format"""
        for i, j in self._percents().items():
            print("{}: {}%".format(i, j))

    def check_for_poly(self, polyatomic):
        """Checks if a given polyatomic is in the given molecule.  Returns a list of atom ratios if true."""        
        if type(polyatomic) == Element:
            polyatomic = Molecule([polyatomic])
        elif type(polyatomic) == Molecule:
            pass
        else:
            raise TypeError("{} is an invalid argument type, expected element or molecule".format(type(polyatomic)))

        #print("[Debug]: Checked for {} in {}:".format(polyatomic, self))
        
        # If the "polyatomic" is a single element, just check for that element
        if len(polyatomic.atom_list) == 1:
            if self.get_atom_count(polyatomic.atom_list[0]):
                #print("\tFound!".format(polyatomic, self))
                #print(self.get_atom_count(polyatomic.atom_list[0]))
                #return True
                return [self.get_atom_count(polyatomic.atom_list[0])]

        atom_ratios = []
        # Scan molecule and produce the atom ratios
        for i in set(polyatomic.atom_list):
            atom_ratios.append(self.get_atom_count(i) / polyatomic.get_atom_count(i))
            if atom_ratios[-1] < 1:
                # If not enough of an element in the polyatomic was found, polyatomic cannot be present
                #print("\tNot found!".format(polyatomic, self))
                return None

        #print("\tFound!".format(polyatomic, self))
        return atom_ratios
    
    def check_solubility(self):
        """Checks the molecule against a SIMPLIFIED set of rules to determine whether its soluble.
Either returns the piece that indicated solubility or None"""
        ## SOLUBILITY CHECK:
        
        # Check for group 1A and ammonium compounds
        for i in [Li, Na, K, Rb, Cs, Fr, N*H*H*H]:
            if self.check_for_poly(i):
                return i
        
        # Check for nitrate compound
        if self.check_for_poly(N*O*O*O):
            return N*O*O*O
        
        # Check for chlorides, bromides, and iodides
        for i in [Cl, Br, I]:
            if any([self.check_for_poly(i) for i in [Ca, Ag, Hg*Hg, Sr, Ba, Pb]]):
                return None
            if self.check_for_poly(i):
                return i

        # Check for sulfates
        if self.check_for_poly(S*O*O*O*O):
            # Check for exceptions to this rule
            
            if any([self.get_atom_count(i) for i in [Ca, Sr, Ba, Pb]]):
                return None
            return S*O*O*O*O
        
        # Check for chlorates
        if self.check_for_poly(Cl*O*O*O):
            return Cl*O*O*O

        # Check for perchlorates
        if self.check_for_poly(Cl*O*O*O*O):
            return Cl*O*O*O*O

        # Check for acetates
        if self.check_for_poly(C*H*H*H*C*O*O):
            return C*H*H*H*C*O*O

        ## INSOLUBILITY CHECK:

        # Check for phoshates
        if self.check_for_poly(P*O*O*O*O):
            return None

        # Check for oxalates
        if self.check_for_poly(C*C*O*O*O*O):
            return None

        # Check for carbonates
        if self.check_for_poly(C*O*O*O):
            return None

        # Check for chromates
        if self.check_for_poly(Cr*O*O*O*O):
            # Check for exceptions
            if self.get_atom_count(Mg):
                return Mg
            return None

        # Check for hydroxides
        if self.check_for_poly(O*H):
            return None

        # Check for oxides
        if self.get_atom_count(O):
            return None

        # Check for sulfides
        if self.get_atom_count(S):
            # Check for exceptions
            for i in [Mg, Ca, Ba]:
                if self.get_atom_count(i):
                    return i
            return None

        # No solubility rule availble, assume insoluble
        print("Warning: no applicable rule found for '{}', assuming insoluble".format(self))
        return None

    def __mul__(self, x):
        """Overload the multiplication operator for 'ease of use'"""
        if type(x) == Element:
            # Molecule * Element -> Molecule
            return Molecule(self.atom_list * self.count + [x])
        elif type(x) == Molecule:
            # Molecule * Molecule -> Molecule
            return Molecule(self.atom_list * self.count + x.atom_list * x.count)
        elif type(x) == int or type(x) == float:
            # n * Molecule -> n amount of Molecule
            return Molecule(self.atom_list, self.count * x)
        else:
            # Invalid operation
            raise TypeError("Unsupported operation between {} and {}".format(type(self), type(x)))

    def compare(self, x):
        """Check if two molecules are equivalent"""
        if type(x) == Molecule:
            list1 = sorted(self.atom_list, key=lambda x: x.atomic_number)
            list2 = sorted(x.atom_list, key=lambda x: x.atomic_number)
            return (list1 == list2) and (self.count == x.count)
        else:
            raise TypeError("Unsupported comparison between {} and {}".format(type(self), type(x)))

    def __rmul__(self, x):
        """Support for reverse-order multiplication"""
        return self.__mul__(x)

    def __add__(self, x):
        if type(x) == Molecule:
            # Molecule + Molecule -> Combination
            return Combination([self, x])
        elif type(x) == Element:
            # Molecule + Element -> Combination
            return Combination([self, Molecule([x])])
        elif type(x) == Combination:
            # Molecule + Combination -> Combination
            return Combination([self] + x.parts)
        else:
            # Invalid operation
            raise TypeError("Unsupported operation between {} and {}".format(type(self), type(x)))

    def __gt__(self, x):
        """Overloads the > operator to declare the reaction in a chemical reaction"""
        if type(x) == Combination:
            return Reaction(Combination([self]), x)
        elif type(x) == Molecule:
            return Reaction(Combination([self]), Combination([x]))
        else:
            raise TypeError("Unsupported operation between {} and {}".format(type(self), type(x)))

class Combination:
    """Class representing a series of elements/molecules"""

    def __init__(self, molecule_list):
        self.parts = molecule_list

    def __repr__(self):
        result = ""
        for i in self.parts:
            result += i.__repr__() + " + "
        return result[:-3]

    def __sub__(self, x):
        """For internal use only"""
        if type(x) == Combination:
            new_list = []
            for i in self.parts:
                for j in x.parts:
                    if i.no_num() == j.no_num():
                        new_list.append(i.copy().set_count(i.count - j.count))
            return Combination(new_list)
                
        else:
            raise TypeError("Unsupported operation between {} and {}".format(type(self), type(x)))

    def __add__(self, x):
        if type(x) == Element or type(x) == Molecule:
            return Combination(self.parts + [x])
        elif type(x) == Combination:
            return Combination(self.parts + x.parts)
        else:
            raise TypeError("Unsupported operation between {} and {}".format(type(self), type(x)))

    def __gt__(self, x):
        if type(x) == Combination:
            return Reaction(self, x)
        elif type(x) == Molecule:
            return Reaction(self, Combination([x]))
        else:
            raise TypeError("Unsupported operation between {} and {}".format(type(self), type(x)))

    def copy(self):
        """Return a copy of this combination"""
        return Combination([i.copy() for i in self.parts])

    def total_mass(self):
        """Get the total mass, in g/mol, of the combination"""
        result = 0
        for i in self.molecule_list:
            result += i.molar_mass
        return result

    def get_count(self, symbol):
        """Gets the number of a molecule within a combination"""
        result = 0
        for i in self.parts:
            if i.__repr__() == symbol:
                result += 1
        return result

    def count(self):
        """Counts the number of elements in a combination"""
        result = {}
        for i in self.parts:
            if type(i) == Element:
                if i.symbol in result:
                    result[i.symbol] += 1
                else:
                    result[i.symbol] = 1
                    
            if type(i) == Molecule:
                for j in i.atom_list:
                    if j.symbol in result:
                        result[j.symbol] += i.count
                    else:
                        result[j.symbol] = i.count
        return result

    def dissolve(self):
        """Dissolves the components in the combination and finds the net ionic equation"""
        print("Given:", self)
        ions = {}
        for i in self.parts:
            temp_molecule = i.singular()
            # Check component for solubility
            j = temp_molecule.check_solubility()
            if j:
                # If component is soluble, isolate the pieces it will dissolve into
                if not (temp_molecule.check_for_poly(j)):
                    print("Debug:")
                    print("j -", j)
                    print("temp_molecule -", temp_molecule)
                    print("temp_molecule.check_for_poly(j) -", temp_molecule.check_for_poly(j))
                    raise Exception("An error has occured")
                num = round(min(temp_molecule.check_for_poly(j)))
                ions[temp_molecule] = ((i.count * num) * j, i.count * temp_molecule.remove_poly(j, count=num))
                
        if len(ions) == 0:
            raise ValueError("Could not verify that any of the reactants were soluble... no net ionic equation?")

        precips = []
        # Run through possible combinations of these ions to find any insoluble precipitates
        for i, j in ions.items():
            for m, n in ions.items():
                if n != j:
                    if not (j[0]*n[0]).check_solubility():
                        precips.append(Reaction(Combination([j[0], n[0]]), Combination([j[0]*n[0]])))
                        print("Found:", j[0]*n[0])
                    if not (j[1]*n[1]).check_solubility():
                        precips.append(Reaction(Combination([j[1], n[1]]), Combination([j[1]*n[1]])))
                        print("Found:", j[1]*n[1])

        if len(precips) == 0:
            print("No net ionic equation found!  No precipitation reaction occurs?")
        else:
            new_precips = []
            for i in precips:
                present = False
                if not any([i.right.parts[0].compare(j.right.parts[0]) for j in new_precips]):
                    new_precips.append(i)

            precips = list(new_precips)
            del new_precips

            print("Form of net ionic equation:")
            for i in precips:
                print(i)

        #print(precips)
        #print(ions)
            
        
class Reaction:
    """Class representing chemical reactions, composed of combinations"""

    def __init__(self, LHS, RHS):
        self.left = LHS
        self.right = RHS    

    def __repr__(self):
        return self.left.__repr__() + " -> " + self.right.__repr__()

    def copy(self):
        """Returns a new, separate copy of a reaction"""
        return Reaction(self.left.copy(), self.right.copy())

    def verify(self):
        """Checks whether a reaction is balanced"""
        return self.left.count() == self.right.count()

    def solve(self, limit=10):
        """Auto balances a chemcial reaction.
        Will continue until reaching the maximum number of a molecule allowed."""
        # First check if the two sides have the same elements
        element_list1 = [i for i, j in self.left.count().items()]
        element_list2 = [i for i, j in self.right.count().items()]
        for i in element_list1:
            if i not in element_list2:
                raise ValueError("Reaction impossible to balance, elements do not match:\nLeft: {}\nRight:{}".format(element_list1, element_list2))

        # Iterate through possible combinations, assuming the reaction is solvable
        if limit > 20:
            l_coeffs = sorted(itertools.product(range(1, limit), repeat=len(self.left.parts)), key=sum)
            r_coeffs = sorted(itertools.product(range(1, limit), repeat=len(self.right.parts)), key=sum)
        else:
            l_coeffs = itertools.product(range(1, limit), repeat=len(self.left.parts))
            r_coeffs = list(itertools.product(range(1, limit), repeat=len(self.right.parts)))

        for i in l_coeffs:
            # First note whether all coefficients are even
            even_left = True
            for j in i:
                if j % 2:
                    even_left = False
                    break
            
            #print("Attempting on left:", i)
            new_left = self.left.copy()
            # Place test values into left-hand side of reaction
            for j in range(len(new_left.parts)):
                new_left.parts[j].set_count(i[j])
            
            for j in r_coeffs:
                # First check if all coefficients are even, and if so, discard (won't be in simplest form)
                even_right = True
                for k in j:
                    if k % 2:
                        even_right = False
                        break
                if even_left and even_right: continue

                # Place test values into right-hand side of reaction and check validity
                new_right = self.right.copy()
                for k in range(len(new_right.parts)):
                    new_right.parts[k].set_count(j[k])
                    
                if Reaction(new_left, new_right).verify():
                    return Reaction(new_left, new_right)

        raise TimeoutError("""Unable to solve the reaction within the specified limit ({}),
Consider raising the limit by using .solve(limit=NEW_LIMIT_HERE)""".format(limit))

    def scale(self, factor):
        """Scales a reaction by an amount"""
        new_left = self.left.copy()
        for i in range(len(new_left.parts)):
            new_left.parts[i].set_count(new_left.parts[i].count * factor)

        new_right = self.right.copy()
        for i in range(len(new_right.parts)):
            new_right.parts[i].set_count(new_right.parts[i].count * factor)

        return Reaction(new_left, new_right)

    def _simulate(self, reactants_dict):
        """Internal function for simulating chemical reactions, takes a
        dictionary of elements matched with moles."""
        product = self.right.parts[0]
        for i, j in zip(self.left.parts, reactants_dict.items()):
            #print("{} moles of {} produces {}".format(j[1], j[0], j[1] * product.count / i.count))
            product_moles = j[1] * product.count / i.count
            try:
                if product_moles < limiting_moles:
                    limiting_moles = product_moles
                    limitor = i.copy()
            except NameError:
                limiting_moles = product_moles
                limitor = i.copy()
                
        #print("Limiting moles:", limiting_moles)
        #print("Limitor:", limitor)

        best = (self.copy()).scale(limiting_moles / product.count)
        #print("Best-case reaction:", best)

        excess = Combination([])
        for i, j, k in zip(best.left.parts, self.left.parts, reactants_dict.items()):
            excess += (i.copy()).set_count(k[1] - i.count)

        return best, excess
    
    def simulate(self):
        """Simulates a chemical reaction with a set of starting conditions.
        Conditions are provided through a user interface"""
        if not self.verify():
            if input("WARNING: Reaction is not balanced, run auto-balance? (Y/N): ").lower() == "n":
                raise ValueError("Reaction '{}' is not balanced!", self)
            else:
                self = self.solve()
        reactors = {}
        print("Reacting:", self)
        print("Enter starting conditions (append 'grams' to denote a value in grams, otherwise moles are assumed):")
        for i in self.left.parts:
            new_reactor = input(i.no_num() + ": ")
            if "grams" in new_reactor:
                new_reactor = abs(float(new_reactor.replace(" grams", "")))
                new_reactor = (i.singular()).grams_to_moles(new_reactor)
                reactors[i.singular()] = new_reactor
            else:
                reactors[i.singular()] = abs(float(new_reactor))

        best, excess = self._simulate(reactors)
        print("Best Possible:", best)
        print("Excess: ", excess)
        
        print("Products (In grams):")
        for i in best.right.parts:
            print(i.singular(), " - ",(i.singular()).moles_to_grams(i.count), "grams")

        print("Excess (In grams):")
        for i in excess.parts:
            if i.count != 0.0: print(i.singular(), " - ",(i.singular()).moles_to_grams(i.count), "grams")

    def sim(self):
        """Alias for simulate"""
        self.simulate()
    

print("Things to change:")
print("\t-Add concept of electrons for elements, as well as ion calculation")
print("\t-Add algorithm for solving net ionic equations given a set of reactants")
print("(Cd*S*O*O*O*O + K*K*S).dissolve()")

# Declare the elements of the periodic table
H = Element("H", 1.008, 1)
He = Element("He", 4.003, 2)
Li = Element("Li", 6.94, 3)
Be = Element("Be", 9.012, 4)
B = Element("B", 10.81, 5)
C = Element("C", 12.01, 6)
N = Element("N", 14.01, 7)
O = Element("O", 16.00, 8)
F = Element("F", 19.00, 9)
Ne = Element("Ne", 20.18, 10)
Na = Element("Na", 22.99, 11)
Mg = Element("Mg", 24.31, 12)
Al = Element("Al", 26.98, 13)
Si = Element("Si", 28.09, 14)
P = Element("P", 30.97, 15)
S = Element("S", 32.06, 16)
Cl = Element("Cl", 35.45, 17)
Ar = Element("Ar", 39.95, 18)
K = Element("K", 39.10, 19)
Ca = Element("Ca", 40.08, 20)
Sc = Element("Sc", 44.96, 21)
Ti = Element("Ti", 47.87, 22)
V = Element("V", 50.94, 23)
Cr = Element("Cr", 52.00, 24)
Mn = Element("Mn", 54.94, 25)
Fe = Element("Fe", 55.85, 26)
Co = Element("Co", 58.93, 27)
Ni = Element("Ni", 58.69, 28)
Cu = Element("Cu", 63.55, 29)
Zn = Element("Zn", 65.38, 30)
Ga = Element("Ga", 69.72, 31)
Ge = Element("Ge", 72.63, 32)
As = Element("As", 74.92, 33)
Se = Element("Se", 78.97, 34)
Br = Element("Br", 79.90, 35)
Kr = Element("Kr", 83.80, 36)
Rb = Element("Rb", 85.47, 37)
Sr = Element("Sr", 87.62, 38)
Y = Element("Y", 88.91, 39)
Zr = Element("Zr", 91.22, 40)
Nb = Element("Nb", 92.91, 41)
Mo = Element("Mo", 95.95, 42)
Tc = Element("Tc", 97, 43)
Ru = Element("Ru", 101.1, 44)
Rh = Element("Rh", 102.9, 45)
Pd = Element("Pd", 106.4, 46)
Ag = Element("Ag", 107.9, 47)
Cd = Element("Cd", 112.4, 48)
In = Element("In", 114.8, 49)
Sn = Element("Sn", 118.7, 50)
Sb = Element("Sb", 121.8, 51)
Te = Element("Te", 127.6, 52)
I = Element("I", 126.9, 53)
Xe = Element("Xe", 131.3, 54)
Cs = Element("Cs", 132.9, 55)
Ba = Element("Ba", 137.3, 56)
La = Element("La", 138.9, 57)
Ce = Element("Ce", 140.1, 58)
Pr = Element("Pr", 140.9, 59)
Nd = Element("Nd", 144.2, 60)
Pm = Element("Pm", 145, 61)
Sm = Element("Sm", 150.4, 62)
Eu = Element("Eu", 152.0, 63)
Gd = Element("Gd", 157.3, 64)
Tb = Element("Tb", 158.9, 65)
Dy = Element("Dy", 162.5, 66)
Ho = Element("Ho", 164.9, 67)
Er = Element("Er", 167.3, 68)
Tm = Element("Tm", 168.9, 69)
Yb = Element("Yb", 173.1, 70)
Lu = Element("Lu", 175.0, 71)
Hf = Element("Hf", 178.5, 72)
Ta = Element("Ta", 180.9, 73)
W = Element("W", 183.8, 74)
Re = Element("Re", 186.2, 75)
Os = Element("Os", 190.2, 76)
Ir = Element("Ir", 192.2, 77)
Pt = Element("Pt", 195.1, 78)
Au = Element("Au", 197.0, 79)
Hg = Element("Hg", 200.6, 80)
Tl = Element("Tl", 204.4, 81)
Pb = Element("Pb", 207.2, 82)
Bi = Element("Bi", 209.0, 83)
Po = Element("Po", 209, 84)
At = Element("At", 210, 85)
Rn = Element("Rn", 222, 86)
Fr = Element("Fr", 223, 87)
Ra = Element("Ra", 226, 88)
Ac = Element("Ac", 227, 89)
Th = Element("Th", 232.0, 90)
Pa = Element("Pa", 231.0, 91)
U = Element("U", 238.0, 92)
Np = Element("Np", 237, 93)
Pu = Element("Pu", 244, 94)
Am = Element("Am", 243, 95)
Cm = Element("Cm", 247, 96)
Bk = Element("Bk", 247, 97)
Cf = Element("Cf", 251, 98)
Es = Element("Es", 252, 99)
Fm = Element("Fm", 257, 100)
Md = Element("Md", 258, 101)
No = Element("No", 259, 102)
Lr = Element("Lr", 262, 103)
Rf = Element("Rf", 267, 104)
Db = Element("Db", 270, 105)
Sg = Element("Sg", 271, 106)
Bh = Element("Bh", 270, 107)
Hs = Element("Hs", 277, 108)
Mt = Element("Mt", 276, 109)
Ds = Element("Ds", 281, 110)
Rg = Element("Rg", 282, 111)
Cn = Element("Cn", 285, 112)
Uut = Element("Uut", 285, 113)
Fl = Element("Fl", 289, 114)
Uup = Element("Uup", 288, 115)
Lv = Element("Lv", 293, 116)
Uus = Element("Uus", 294, 117)
Uuo = Element("Uuo", 294, 118)

#Unob = Element("Unobtanium", 999, 999)

# Declare some common molecules and diatomics

H2O = H*H*O
water = H2O
CO2 = C*O*O
NH3 = N*H*H*H
ammonia = NH3
O2 = O*O
H2 = H*H
N2 = N*N
F2 = F*F
I2 = I*I
Cl2 = Cl*Cl
Br2 = Br*Br

# Declare some polyatomics
OH = O*H

import itertools
from math import gcd

DEBUG = False
do_debug = False

sub = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
sup = str.maketrans("0123456789+-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻")

def find_gcd(l):
    x = gcd(l[0], l[1])
    for i in range(2, len(l)):
        x = gcd(x, l[i])
    return x

def docs(search=None):
    import sys, inspect
    for i, j in inspect.getmembers(sys.modules[__name__]):
        if inspect.isclass(j) and i != "__loader__":
            if search:
                if not(search == i or search == j):
                    continue
            print(i, "-", j.__doc__)
            for k, m in inspect.getmembers(j):
                if inspect.isfunction(m):
                    print("\t", k, ":", m.__doc__, "\n")
            print("")
    if search is None: print("\nExample of syntax:\nMg*ep*ep + Ba*(2*NO3) > Ba*ep*ep + Mg*(2*NO3)\nBecomes:\n[Mg]²⁺ + Ba(NO₃)₂ -> [Ba]²⁺ + Mg(NO₃)₂")

class Charge:
    """Utility class for ionizing molecules."""

    def __init__(self, charge):
        self.charge = charge

    def __repr__(self):
        result = "e"
        if self.charge > 0:
            result += "+"
        result += str(self.charge).translate(sup)
        if abs(self.charge) == 1:
            result = result[:-1]
        return result

    def __mul__(self, x):
        """Multiply element/molecule by 'en' or 'ep' to ionize element/molecule.  Also accepts numbers for condensation purposes."""
        if type(x) in [float, int]:
            return Charge(self.charge * x)

        elif type(x) == Charge:
            return Charge(self.charge + x.charge)

        elif type(x) == Element:
            return Molecule({x: 1}, count=1, charge=self.charge)

        elif type(x) == Molecule:
            result = x.copy()
            result.charge += self.charge
            return result

        else:
            raise TypeError("Unsupported operation between {} and {}".format(type(self), type(x)))
    
    def __rmul__(self, x):
        """Reverse-order multiplication compatability that typically does nothing"""
        return self.__mul__(x)

class Element:
    """Class for basic elements"""

    elements = {}

    def __init__(self, symbol, molar_mass, atomic_number, oxidation):
        """Elements are referenced by their atomic symbol.  They may be multiplied together (or with numbers) to create molecules"""
        self.sym = symbol
        self.mol_mass = molar_mass
        self.atomic_number = atomic_number
        #self.electrons = atomic_number
        self.charge = 0
        self.neutrons = round(molar_mass) - atomic_number
        self.ox = oxidation  
        Element.elements[symbol] = self

    def molar_mass(self):
        """Returns the molar mass of the element"""
        return self.mol_mass
    
    def oxidation(self):
        """Get the oxidation state of the thing"""
        return self.ox
    
    def get_atoms(self):
        """Returns a dictionary of element-count pairs in the polyatomic"""
        return {self: 1}

    def extract_count(self):
        """Returns 1.  Purely for compatibility purposes"""
        return 1

    def copy(self):
        """Returns self.  Purely for compatibility purposes"""
        return self
    
    def __repr__(self):
        """Return the element's symbol for string representation"""
        return str(self.sym)

    def symbol(self):
        """Return the element's symbol for string representation"""
        return str(self.sym)

    def grams_to_moles(self, grams):
        """Converts grams of an element to moles"""
        return grams / self.molar_mass()

    def g2mol(self, grams):
        """Alias for grams_to_moles"""
        return self.grams_to_moles(grams)

    def moles_to_grams(self, moles):
        """Converts moles of an element to grams"""
        return self.molar_mass() * moles

    @classmethod
    def list(cls):     # This is organized/placed poorly
        """Print a list of all registered elements"""
        for i, j in cls.elements.items():
            print("{}{} - Atomic #: {},{} Molar Mass: {}".format(i,
                                                                 abs(3-len(i))*" ",
                                                                 j.atomic_number,
                                                                 (3 - len(str(j.atomic_number)))*" ",
                                                                 j.molar_mass()))

    @staticmethod
    def moles_to_atoms(moles):
        """Converts moles of an element to grams"""
        return moles * 6.022e23

    @staticmethod
    def atoms_to_moles(atoms):
        """Converts atoms of an element to moles"""
        return atoms / 6.022e23
    
    def grams_to_atoms(self, grams):
        """Convert grams of an element to atoms"""
        return Element.moles_to_atoms(self.grams_to_moles(grams))

    def atoms_to_grams(self, atoms):
        """Converts atoms of an element to grams"""
        return self.moles_to_grams(Element.atoms_to_moles(atoms))

    def __mul__(self, x):
        if type(x) == Element:
            # Element * Element -> Molecule
            result = Molecule({self: 1})
            result.append(x)
            return result
        
        elif type(x) == Molecule:
            # Element * Molecule -> Molecule
            count = x.count
            result = Molecule(dict((i, j*count) for i, j in x.parts.items()))
            result.append(self)
            return result

        elif type(x) == Polyatomic:
            # Element * Polyatomic -> Molecule
            return Molecule({self: 1, x: x.extract_count()})
            
        elif type(x) == int:
            # Element * n -> Molecule of 'n' Elements
            return Molecule({self: x})

        elif type(x) == Charge:
            return Molecule({self: 1}, count=1, charge=x.charge)
        
        else:
            # Invalid operation
            raise TypeError("Unsupported operation between {} and {}".format(type(self), type(x)))

    def __rmul__(self, x):
        if type(x) == int:
            # n * Element -> 'n' Molecules of monatomic Element
            return Molecule({self: 1}, count=x)
        else:
            # Support for multiplication... BUT IN REVERSE!!!
            return self.__mul__(x)

    def __add__(self, x):
        result = Combination({Molecule({self: 1}): 1})
        if type(x) == Element:
            # Element + Element -> Combination
            result += Molecule({x: 1})
            return result
        elif type(x) == Molecule or type(x) == Polyatomic:
            # Element + Molecule -> Combination
            result += x.copy()
            return result
        elif type(x) == Combination:
            # Element + Combination -> Combination
            result += x.copy()
            return result
        else:
            # Invalid operation
            raise TypeError("Unsupported operation between {} and {}".format(type(self), type(x)))

class Molecule:
    """Class for molecules"""
    
    def __init__(self, parts_dict, count=1, charge=0):
        self.parts = parts_dict
        self.count = count
        self.mol_mass = count * sum([j*(i.molar_mass()) for i, j in parts_dict.items()])
        self.count = count
        self.charge = charge

    def __hash__(self):
        return hash((tuple(self.get_atoms().items()), self.count))

    def __eq__(self, x):
        return self.compare(x, ignore_num=False)

    def __repr__(self):
        return self.symbol()

    def extract_count(self):
        """Sets the count to 1 and returns what it was before"""
        count = self.count
        self.count = 1
        return count

    def append(self, x):
        """Appends something onto the molecule"""
        amount = x.extract_count()

        if x in self.parts:
            self.parts[x] += amount
        else:
            self.parts[x] = amount

    def molar_mass(self):
        """Calculate the molar mass of the element"""
        self.mol_mass = self.count * sum([j *(i.molar_mass()) for i, j in self.parts.items()])
        return self.mol_mass

    def oxidation(self):
        """Get the oxidation state/charge of the thing"""
        return self.charge

    def symbol(self):
        """Returns the symbol of the element"""
        result = ""
        if self.charge not in [0, None]:
            charge = self.charge
            result += "["
        if self.count > 1:
            result += str(self.count)
        for i, j in self.parts.items():
            if type(i) == Element:
                result += i.symbol()
                if j > 1:
                    result += str(j).translate(sub)
            else:
                if j > 1:
                    result += "(" + i.symbol() + ")" + str(j).translate(sub)
                else:
                    result += i.symbol()
        if self.charge not in [0, None]:
            result += "]"
            if charge > 0:
                if charge > 1:
                    result += str(charge).translate(sup)
                result += "+".translate(sup)
            elif charge < 0:
                if charge < -1:
                    result += str(abs(charge)).translate(sup)
                result += "-".translate(sup)
        return result

    def copy(self):
        return Molecule(dict(self.parts), self.count, self.charge)

    def grams_to_moles(self, grams, element=None):
        """Converts grams of a molecule to moles.  If an element is specified, the number of
        moles of a particular element will be isolated and returned"""
        if element:
            return self.parts[element] * grams / self.molar_mass()
        else:
            return grams / self.molar_mass()

    def g2mol(self, grams):
        """Alias for grams_to_moles"""
        return self.grams_to_moles(grams)
    
    def grams_to_molecules(self, grams, element=None):
        """Converts grams of a molecule to the number of molecules."""
        return Element.moles_to_atoms(self.grams_to_moles(grams, element))

    def grams_to_atoms(self, grams, element):
        """Computes how many atoms of a particular element are present in the given sample"""
        return self.grams_to_molecules(grams, element)

    def moles_to_grams(self, moles, element=None):
        """Converts the number of moles of a molecule to grams.  If an element is specified,
        the mass of that element will be isolated and returned instead."""
        if element:
            return moles * (Molecule({element: self.parts[element]})).molar_mass()
        else:
            return moles * self.molar_mass()

    def molecules_to_moles(self, molecules, element=None):
        """Converts molecules of a sample to moles.  If an element is specified, the number of
        moles of that element will be isolated and returned instead."""
        if element:
            return self.parts[element] * molecules / 6.022e23
        else:
            return molecules / 6.022e23

    def moles_to_molecules(self, moles):
        """Converts a specified number of moles to molecules"""
        return moles * 6.022e23

    def molecules_to_grams(self, molecules, element=None):
        """Converts molecules of a sample to grams.  If an element is specified, the mass of that
        element will be isolated and returned instead."""
        return self.moles_to_grams(self.molecules_to_moles(molecules), element)

    def get_atoms(self):
        """Returns a dictionary of element-count pairs in the polyatomic"""
        result = {}
        # Scan the parts dictionary of the molecule
        for i, j in self.parts.items():
            # Get the elements composing the individual ion
            part = i.get_atoms()
            # Add the elements to the resultant dictionary
            for k, m in part.items():
                if k in result:
                    result[k] += m*j
                else:
                    result[k] = m*j
        return result

    def compare(self, x, ignore_num=True):
        """Checks if two molecules are describing the same thing"""
        if ignore_num:
            return (self.get_atoms()) == (x.get_atoms())
        else:
            return (self.get_atoms()) == (x.get_atoms()) and (self.count == x.count)

    def contains(self, element):
        """Checks if the molecule contains a specified element"""
        atoms_dict = self.get_atoms()

        for i, j in element.get_atoms().items():
            if not (i in atoms_dict and atoms_dict[i] >= j):
                return False
        return True

    def check_solubility(self, ignore_ph=True):
        """Checks if a molecule will dissociate in water"""

        if not ignore_ph:
            if self == Ca*(2*OH):
                return False
            if self in acids:
                if not self in strong_acids:
                    return False
                return True
            elif self in bases:
                if not self in strong_bases:
                    return False
                return True

        # Apply solubility rules:
        if any([self.contains(i) for i in [Li, Na, K, Cs, Rb, Fr, NH4,
                                           NO3, CH3CO2, ClO3, ClO4]]):
            return True
        
        if any([self.contains(i) for i in [Cl, Br, I]]):
            if any([self.contains(i) for i in [Pb, Ag, Hg*Hg]]):
                return False
            return True
        
        if self.contains(SO4):
            if any([self.contains(i) for i in [Ca, Ag, Hg*Hg, Sr, Ba, Pb]]):
                return False
            return True

        if any([self.contains(i) for i in [OH, O*O, S*S, CO3, PO4, C2O4]]):
            return False

        if self.contains(CrO4):
            if self.contains(Mg):
                return True
            return False

    def atomize(self):
        """Converts a monotomic bonding coefficient to a quantity coefficient"""
        copy = self.copy()
        if len(self.parts) == 1:
            count = copy.extract_count()
            element = list(copy.parts.items())[0][0]
            count *= list(copy.parts.items())[0][1]
            return Molecule({element: 1}, count=count)
        else:
            return copy

    def deatomize(self):
        """Transfers the quantity coefficient into a bond coefficient"""
        copy = self.copy()
        count = copy.extract_count()
        new_parts = dict((i, j*count) for i, j in copy.parts.items())
        return Molecule(new_parts)

    def singular(self):
        """Returns a singular molecule of a given molecule"""
        copy = self.copy()
        copy.extract_count()
        return copy
    
    def dissolve(self, ignore_ph=True):
        """Dissociates the ions of a molecule, if it's soluble"""
        if not self.check_solubility(ignore_ph):
            return self.copy()

        count = self.extract_count()
        if len(self.parts.items()) != 2:
            raise ValueError("There is currently no system in place to handle this argument")

        copy = self.copy()

        # Assign oxidation states
        result = Combination({})
        for i, j in self._assign_oxidation().items():
            result += Molecule({i: 1}, count=copy.parts[i]*count, charge=j)

        # Check for potentially unregistered acid
        H_ion = H*ep
        if H_ion in result.parts:
            if self not in acids:
                print("Warning: Is {} an acid that is not currently registered in the lookup table?".format(self))
        
        return result

    def _percents(self):
        """Get the percent mass of each element in a molecule"""
        factor = 100 / self.molar_mass()
        percents_dict = {}

        for i, j in self.parts.items():
            if i not in percents_dict:
                percents_dict[i] = self.parts[i] * i.molar_mass() * factor

        return percents_dict

    def percents(self):
        """Wrapper for the _percents function"""
        for i, j in self._percents().items():
            print("{}: {}%".format(i, j))

    def _assign_oxidation(self):
        """Determines the oxidation states of the ions within a molecule and returns them as a dictionary"""
        result = {}
        current_charge = self.charge

        # Apply rules for elemental atoms and singular ions
        ion_list = [i for i, j in self.parts.items()]
        #print(ion_list)
        if len(ion_list) == 1:
            if self.charge == 0:
                return {ion_list[0]: 0}
            else:
                return {ion_list[0]: self.charge / self.parts[ion_list[0]]}

        # Assume polyatomic ions are correct
        for i in ion_list:
            if type(i) == Polyatomic:
                result[i] = i.oxidation()
                current_charge -= i.oxidation() * self.parts[i]
                ion_list.remove(i)
            if len(ion_list) == 1:
                result[ion_list[0]] = current_charge / self.parts[ion_list[0]]
                return result

        # Apply rules for group 1A compounds
        for i in ion_list:
            if i in [Li, Na, K, Rb, Cs, Fr]:
                result[i] = 1
                current_charge -= 1 * self.parts[i]
                ion_list.remove(i)
            if len(ion_list) == 1:
                result[ion_list[0]] = current_charge / self.parts[ion_list[0]]
                return result

        # Apply rules for group 2A compounds
        for i in ion_list:
            if i in [Be, Mg, Ca, Sr, Ba, Ra]:
                result[i] = 2
                current_charge -= 2 * self.parts[i]
                ion_list.remove(i)
            if len(ion_list) == 1:
                result[ion_list[0]] = current_charge / self.parts[ion_list[0]]
                return result

        # Apply rules for fluorine, hydrogen, and oxygen
        if F in ion_list:
            result[F] = -1
            current_charge += 1 * self.parts[F]
            ion_list.remove(F)

        if len(ion_list) == 1:
            result[ion_list[0]] = current_charge / self.parts[ion_list[0]]
            return result
        
        if H in ion_list:
            result[H] = 1
            current_charge -= 1 * self.parts[H]
            ion_list.remove(H)

        if len(ion_list) == 1:
            result[ion_list[0]] = current_charge / self.parts[ion_list[0]]
            return result

        if O in ion_list:
            result[O] = -2
            current_charge += 2 * self.parts[O]
            ion_list.remove(O)

        if len(ion_list) == 1:
            result[ion_list[0]] = current_charge / self.parts[ion_list[0]]
            return result

        # Apply rules for group 7A ions
        for i in ion_list:
            if i in [F, Cl, Br, I, At, Ts]:
                result[i] = -1
                current_charge += 1 * self.parts[i]
                ion_list.remove(i)
            if len(ion_list) == 1:
                result[ion_list[0]] = current_charge / self.parts[ion_list[0]]
                return result

        # Apply rules for group 6A ions
        for i in ion_list:
            if i in [F, Cl, Br, I, At, Ts]:
                result[i] = -2
                current_charge += 2 * self.parts[i]
                ion_list.remove(i)
            if len(ion_list) == 1:
                result[ion_list[0]] = current_charge / self.parts[ion_list[0]]
                return result

        if len(ion_list) != 1:
            raise ValueError("Error: unsufficient information to assign oxidation states")

    def assign_oxidation(self):
        """Wrapper for _assign_oxidation to print results in a more readable format"""
        for i, j in self._assign_oxidation().items():
            print(i, ":", j)

    def __mul__(self, x):
        """Overload the multiplication operator for 'ease of use'"""
        if type(x) == Element:
            # Molecule * Element -> Molecule
            result = self.copy()
            result.append(x)
            return result
        elif type(x) == Molecule:
            # Molecule * Molecule -> Molecule
            result = self.copy()
            copy2 = x.copy()
            copy2 = copy2.atomize()
            result.append(copy2)
            return result
        elif type(x) == Polyatomic:
            # Molecule * polyatomic
            result = self.copy()
            result = result.deatomize()
            copy2 = x.copy()
            result.append(copy2)
            return result
        elif type(x) == int or type(x) == float:
            # n * Molecule -> n amount of Molecule
            return Molecule(dict(self.parts), self.count * x, self.charge)
        elif type(x) == Charge:
            result = self.copy()
            result.charge += x.charge
            return result
        else:
            # Invalid operation
            raise TypeError("Unsupported operation between {} and {}".format(type(self), type(x)))

    def __rmul__(self, x):
        """Support for reverse-order multiplication"""
        return self.__mul__(x)

    def __add__(self, x):
        count = self.extract_count()
        result = Combination({self: count})
        if type(x) in [Molecule, Element, Polyatomic]:
            # Molecule + Molecule -> Combination OR
            # Molecule + Element -> Combination 
            result.append(x)
            return result
        elif type(x) == Combination:
            # Molecule + Combination -> Combination
            result += x
            return result
        else:
            # Invalid operation
            raise TypeError("Unsupported operation between {} and {}".format(type(self), type(x)))

    def __gt__(self, x):
        """Overloads the > operator to declare the reaction in a chemical reaction"""
        count = self.extract_count()
        comb = Combination({self: count})
        if type(x) == Combination:
            return Reaction(comb, x)
        elif type(x) == Molecule:
            count2 = x.extract_count()
            return Reaction(comb, Combination({x: count2}))
        else:
            raise TypeError("Unsupported operation between {} and {}".format(type(self), type(x)))
        
class Polyatomic:
    """Class for polyatomic ions"""

    polyatomic_list = []

    def __init__(self, molec, oxidation, count=1):
        """Polyatomics are called by entering a single term.  For example, C₂H₃O₂⁻ is referenced as C2H3O2"""
        if type(molec) == dict:
            self.parts = dict(molec)
        else:
            self.parts = dict(molec.parts)
        self.mol_mass = sum([j*i.molar_mass() for i, j in self.parts.items()])
        self.ox = oxidation
        self.count = count
        self.sym = self.symbol()
        self.charge = oxidation
        if self not in Polyatomic.polyatomic_list:
            Polyatomic.polyatomic_list.append(self)
        else:
            if DEBUG: print("WARNING: Conlict detected in defining", self, "({})".format(oxidation))
            pass

    def __hash__(self):
        return hash((tuple(self.get_atoms().items()), self.count))

    def __eq__(self, x):
        return self.compare(x, ignore_num=False)

    def __repr__(self):
        return self.symbol()

    def grams_to_moles(self, grams, element=None):
        """Converts grams of a molecule to moles.  If an element is specified, the number of
        moles of a particular element will be isolated and returned"""
        if element:
            return self.parts[element] * grams / self.molar_mass()
        else:
            return grams / self.molar_mass()

    def g2mol(self, grams):
        """Alias for grams_to_moles"""
        return self.grams_to_moles(grams)

    def moles_to_grams(self, moles, element=None):
        """Converts the number of moles of a molecule to grams.  If an element is specified,
        the mass of that element will be isolated and returned instead."""
        if element:
            return moles * (Molecule({element: self.parts[element]})).molar_mass()
        else:
            return moles * self.molar_mass()
    
    def molar_mass(self):
        """Gets the molar mass of the polyatomic ion"""
        self.mol_mass = sum([j*i.molar_mass() for i, j in self.parts.items()])
        return self.mol_mass

    def oxidation(self):
        """Get the oxidation state of the thing"""
        return self.ox
    
    def symbol(self):
        """Returns the symbol of the element"""
        result = ""
        if self.count > 1:
            result += str(self.count)
        for i, j in self.parts.items():
            if type(i) == Element:
                result += i.symbol()
                if j > 1:
                    result += str(j).translate(sub)
            else:
                if j > 1:
                    result += "(" + i.symbol() + ")" + str(j).translate(sub)
                else:
                    result += i.symbol()
        return result

    def atomize(self):
        """Converts a monotomic bonding coefficient to a quantity coefficient"""
        if len(self.parts) == 1:
            copy = self.copy()
            count = copy.extract_count()
            ion = list(copy.parts.items())[0][0]
            count *= list(copy.parts.items())[0][1]
            return Polyatomic({ion: 1}, count=count)
        else:
            return self.copy()

    def deatomize(self):
        """Transfers the quantity coefficient into a bond coefficient"""
        copy = self.copy()
        count = copy.extract_count()
        count *= list(copy.parts.items())[0][1]
        return Molecule({copy: count})

    def append(self, x):
        """Appends something onto the molecule"""
        amount = x.extract_count()

        if self == x:
            return 2*self
        
        if x in self.parts:
            self.parts[x] += amount
        else:
            self.parts[x] = amount

    def extract_count(self):
        """Removes the quantity part of the polyatomic and returns that value"""
        count = self.count
        self.count = 1
        return count

    def copy(self):
        """Returns a new copy of the polyatomic"""
        return Polyatomic(dict(self.parts), self.oxidation(), self.count)

    def get_atoms(self):
        """Returns a dictionary of element-count pairs in the polyatomic"""
        result = {}
        # Scan the parts dictionary of the polyatomic
        for i, j in self.parts.items():
            # Get the elements composing the individual ion
            part = i.get_atoms()
            # Add the elements to the resultant dictionary
            for k, m in part.items():
                if k in result:
                    result[k] += m*j
                else:
                    result[k] = m*j
        return result
    
    def compare(self, x, ignore_num=True):
        """Checks if two molecules are describing the same thing"""
        if ignore_num:
            return (self.get_atoms()) == (x.get_atoms())
        else:
            return (self.get_atoms()) == (x.get_atoms()) and (self.count == x.count)

    def _assign_oxidation(self):
        """Returns the oxidation state of the ion"""
        return {self: self.oxidation()}

    def assign_oxidation(self):
        """Wrapper for _assign_oxidation"""
        print(self, "-", self.oxidation())
    
    def __mul__(self, x):
        """Overload the multiplication operator for 'ease of use'"""
        if type(x) == Element:
            # Polyatomic * Element -> Molecule
            copy1 = self.copy()
            copy1 = copy1.atomize()
            count1 = copy1.extract_count()
            return Molecule({copy1: count1, x: 1})
        elif type(x) == Molecule:
            # Polyatomic * Molecule -> Molecule
            # Molecule * polyatomic
            result = self.copy()
            result.append(x)
            return result
        elif type(x) == Polyatomic:
            # Polyatomic * polyatomic
            if self == x:
                copy = self.copy()
                count = copy.extract_count()
                return Molecule({copy: 2})
            
            copy1 = self.copy()
            count1 = copy1.extract_count()
            copy2 = x.copy()
            count2 = copy2.extract_count()
            return Molecule({copy1: count1, copy2: count2})
        elif type(x) == int or type(x) == float:
            # n * Polyatomic -> n amount of Polyatomic
            result = self.copy()
            result.count *= x
            return result
            
        else:
            # Invalid operation
            raise TypeError("Unsupported operation between {} and {}".format(type(self), type(x)))

    def __rmul__(self, x):
        """Support for reverse-order multiplication"""
        return self.__mul__(x)

    def __add__(self, x):
        if type(x) == Molecule:
            # Molecule + Molecule -> Combination
            count = self.extract_count()
            result = Combination({self: count})
            result.append(x)
            return result
        elif type(x) == Element:
            # Molecule + Element -> Combination
            count = self.extract_count()
            return Combination({self: count, Molecule({x: 1}): 1})
        elif type(x) == Combination:
            # Molecule + Combination -> Combination
            count = self.extract_count()
            result = Combination({self: count})
            result += x
            return result
        else:
            # Invalid operation
            raise TypeError("Unsupported operation between {} and {}".format(type(self), type(x)))

    @classmethod
    def list(cls):
        """Lists all polyatomic ions currently available"""
        for i in cls.polyatomic_list:
            print("{} ({})".format(i, i.charge))

class Combination:
    """Class representing a combination of elements/molecules"""

    def __init__(self, molecule_dict):
        """Combinations are declared through the use of '+' between components"""
        self.parts = molecule_dict

    def __hash__(self):
        result = {}
        # Scan the parts dictionary of the combination
        for i, j in self.parts.items():
            # Get the elements composing the individual molecules
            part = i.get_atoms()
            # Add the elements to the resultant dictionary
            for k, m in part.items():
                if k in result:
                    result[k] += m*j
                else:
                    result[k] = m*j

        print(result)
        return hash(tuple(result.items()))

    def __eq__(self, x):
        return hash(self) == hash(x)

    def __repr__(self):
        return self.symbol()

    def __add__(self, x):
        result = self.copy()
        if type(x) == Element or type(x) == Molecule or type(x) == Polyatomic:
            result.append(x)
            return result
        elif type(x) == Combination:
            for i, j in x.parts.items():
                result.append(j*i)
            return result
        else:
            raise TypeError("Unsupported operation between {} and {}".format(type(self), type(x)))

    def __gt__(self, x):
        if type(x) == Combination:
            return Reaction(self, x)
        elif type(x) == Molecule:
            count = x.extract_count()
            return Reaction(self, Combination({x: count}))
        else:
            raise TypeError("Unsupported operation between {} and {}".format(type(self), type(x)))

    def symbol(self):
        """Returns the symbolic representation of the combination"""
        result = ""
        for i, j in self.parts.items():
            if j != 1:
                result += str(j)
            result += i.symbol() + " + "
        return result[:-3]

    def append(self, x):
        if type(x) not in [Element, Molecule, Polyatomic]:
            raise TypeError("Unable to append '{}' to 'Combination'".format(type(x)))
        
        count = x.extract_count()
        
        if x in self.parts:
            self.parts[x] += count
        else:
            self.parts[x] = count

    def get_atoms(self):
        """Returns a dictionary of element-count pairs in the entire combination"""
        result = {}
        # Scan the parts dictionary of the polyatomic
        for i, j in self.parts.items():
            # Get the elements composing the individual ion
            part = i.get_atoms()
            # Add the elements to the resultant dictionary
            for k, m in part.items():
                if k in result:
                    result[k] += m*j
                else:
                    result[k] = m*j
        return result

    def copy(self):
        """Return a copy of this combination"""
        return Combination(dict((i.copy(), j) for i, j in self.parts.items()))

    def compare(self, x):
        """Check if two combinations are balanced"""
        elements1 = [i for i, j in self.get_atoms().items()]
        elements2 = [i for i, j in x.get_atoms().items()]
        for i in elements1:
            if i not in elements2:
                return False
        return True

    def total_mass(self):
        """Get the total mass, in g/mol, of the combination"""
        result = 0
        for i, j in self.parts.items():
            result += i.molar_mass() * j
        return result

    def dissolve(self, ignore_ph=True):
        """Dissociates the components of the combination"""
        result = Combination({})
        for i, j in self.parts.items():
            result += (j*i).dissolve(ignore_ph)
        return result

    def _assign_oxidation(self):
        """Assigns oxidation numbers to the ions in the combination"""
        stage1 = []
        parts = list(dict(self.parts).items())
        parts.sort(key=lambda x: x[0].molar_mass())
        #print(parts)
        
        for i, j in parts:
            copy = i.copy()
            copy = j * copy
            stage1.append(copy._assign_oxidation())

        stage2 = []
        stage2.append([])
        
        for i in stage1:
            for j, k in i.items():
                for m, n in enumerate(stage2):
                    if [j, k] not in n:
                        stage2[m].append([j, k])

        stage3 = []
        for i, j in enumerate(stage2):
            stage3.append({})
            stage3[i] = dict(j)

        return stage3
                
        
    def assign_oxidation(self):
        """Wrapper for the _assign_oxidation function to make the output more readable"""
        result = self._assign_oxidation()

        for i in result:
            for j, k in i.items():
                print(j, "-", k)

class Reaction:
    """Class representing chemical reactions, composed of combinations"""

    def __init__(self, LHS, RHS):
        """Reactions are declared through the use of '>' between two combinations"""
        self.left = LHS
        self.right = RHS    

    def __repr__(self):
        return self.left.symbol() + " -> " + self.right.symbol()

    def copy(self):
        """Returns a new, separate copy of a reaction"""
        return Reaction(self.left.copy(), self.right.copy())

    def verify(self):
        """Checks whether a reaction is balanced"""
        left_atoms = self.left.get_atoms()
        right_atoms = self.right.get_atoms()
        for i, j in left_atoms.items():
            if i in right_atoms:
                if j != right_atoms[i]:
                    return False
            else:
                return False
        return True

    def verify_charge(self):
        """Checks whether the charges in a reaction are conserved, let alone possible"""
        # Determine net charges on both sides of equation
        left_charge = 0
        for i, j in self.left.parts.items():
            left_charge += j * i.charge

        right_charge = 0
        for i, j in self.right.parts.items():
            right_charge += j * i.charge

        # Check for noninteger charges
        if round(left_charge) != left_charge:
            return False
        if round(right_charge) != right_charge:
            return False

        for i in self.left._assign_oxidation():
            for j, k in i.items():
                if round(k) != k:
                    return False

        for i in self.right._assign_oxidation():
            for j, k in i.items():
                if round(k) != k:
                    return False

        return left_charge == right_charge

    def solve(self, limit=10):
        """Auto balances a chemcial reaction.  Will continue until reaching the maximum number of a molecule allowed."""
        # First check if the two sides have the same elements
        if not self.left.compare(self.right):
            raise ValueError("Reaction impossible to balance, elements do not match:\nLeft: {}\nRight:{}".format([i for i, j in self.left.get_atoms().items()],
                                                                                                                 [i for i, j in self.right.get_atoms().items()]))

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
            new_left.parts = dict(zip([j[0] for j in new_left.parts.items()], i))
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
                new_right.parts = dict(zip([k[0] for k in new_right.parts.items()], j))

                possible = Reaction(new_left, new_right)
                if possible.verify() and possible.verify_charge():
                    return possible

        raise TimeoutError("""Unable to solve the reaction within the specified limit ({}),
Consider raising the limit by using .solve(limit=NEW_LIMIT_HERE)""".format(limit))

    def scale(self, factor):
        """Scales a reaction by an amount"""
        new_left = self.left.copy()
        new_left.parts = dict((i, j*factor) for i, j in new_left.parts.items())

        new_right = self.right.copy()
        new_right.parts = dict((i, j*factor) for i, j in new_right.parts.items())
        return Reaction(new_left, new_right)

    def _simulate(self, reactants_dict):
        """Internal function for simulating chemical reactions, takes a dictionary of elements matched with moles."""
        product = list(self.right.parts.items())[0][1]

        # Determine the most amount of moles of a product that can be produced, and scale the equation accordingly
        for i, j in reactants_dict.items():
            product_moles = j / self.left.parts[i] * product
            try:
                if product_moles < limiting_moles:
                    limiting_moles = product_moles
            except NameError:
                limiting_moles = product_moles

        best = self.scale(limiting_moles / product)

        # Determine how much of the starting reactants are left over
        excess = Combination({})
        for i, j in best.left.parts.items():
            if (reactants_dict[i] - j) != 0.0:
                excess += (reactants_dict[i] - j) * i

        return best, excess
    
    def simulate(self):
        """Simulates a chemical reaction with a set of starting conditions.  Conditions are provided through a user interface"""
        #raise NotImplemented
        if not self.verify():
            if input("WARNING: Reaction is not balanced, run auto-balance? (Y/N): ").lower() == "n":
                raise ValueError("Reaction '{}' is not balanced!", self)
            else:
                self = self.solve()
        reactors = {}
        print("Reacting:", self)
        print("Enter starting conditions (append 'grams' to denote a value in grams, otherwise moles are assumed):")
        for i, j in self.left.parts.items():
            new_reactor = input(i.symbol() + ": ")
            if "grams" in new_reactor:
                new_reactor = abs(float(new_reactor.replace(" grams", "")))
                new_reactor = i.grams_to_moles(new_reactor)
                reactors[i] = new_reactor
            else:
                reactors[i] = abs(float(new_reactor))

        best, excess = self._simulate(reactors)
        print("Best Possible:", best)
        print("Excess: ", excess)
        
        print("Products (In grams):")
        for i, j in best.right.parts.items():
            print(i, " - ", i.moles_to_grams(j), "grams")

        print("Excess (In grams):")
        for i, j in excess.parts.items():
            if j != 0.0: print(i, " - ", i.moles_to_grams(j), "grams")

    def sim(self):
        """Alias for simulate"""
        self.simulate()

    def dissolve(self, ignore_ph=True):
        """Dissolves both sides of the chemical equation"""
        return Reaction(self.left.dissolve(not ignore_ph), self.right.dissolve())

    def net_ionic(self, auto_balance=False):
        """Generate the net ionic equation of a reaction thing"""
        copy = self.copy()

        if not auto_balance:
            if not copy.verify():
                user = input("Warning: reaction is not balanced, run auto-balancer? (Y/N): ")
                if user.lower() == "y":
                    copy = copy.solve()
        else:
            copy = copy.solve()
        
        dissolved_left = copy.left.dissolve(False)
        dissolved_right = copy.right.dissolve(True)
        dissolved = copy.dissolve()
        new_left = Combination({})
        new_right = Combination({})

        # Cancel out matching components
        for i, j in dissolved.left.parts.items():
            if i not in dissolved.right.parts:
                new_left += j * i.copy()
        for i, j in dissolved.right.parts.items():
            if i not in dissolved.left.parts:
                new_right += j * i.copy()

        # Check if the reaction is reducable
        factor = find_gcd([j for i, j in new_left.parts.items()] + [j for i, j in new_right.parts.items()])
        new_left.parts = dict((i, int(j/factor)) for i, j in new_left.parts.items())
        new_right.parts = dict((i, int(j/factor)) for i, j in new_right.parts.items())
        result = Reaction(new_left, new_right)

        if len(result.left.parts) == 0 and len(result.right.parts) == 0:
            print("No net ionic equation could be found")
            return None

        if not result.verify():
            raise ValueError("An error occured: resulting ionic equation cannot be balanced\n\t-{}".format(result))

        #print("Note: remember to check if one of the molecules is a weak acid.  If it is, it is written as though it is insoluble.")
        
        return result

    def assign_oxidation(self):
        """Assigns oxidation states to all ions in the reaction"""
        print("Left:")
        for i in self.left._assign_oxidation():
            for j, k in i.items():
                print(j, "-", k)
        
        print("\nRight:")
        for i in self.right._assign_oxidation():
            for j, k in i.items():
                print(j, "-", k)
        

    def solve_redox(self):
        """"""
        pass

# Declare elements
H = Element("H", 1.008, 1, 1)
He = Element("He", 4.003, 2, 0)
Li = Element("Li", 6.94, 3, 1)
Be = Element("Be", 9.012, 4, 2)
B = Element("B", 10.81, 5, None)
C = Element("C", 12.01, 6, -4)
N = Element("N", 14.01, 7, -3)
O = Element("O", 16.00, 8, -2)
F = Element("F", 19.00, 9, -1)
Ne = Element("Ne", 20.18, 10, 0)
Na = Element("Na", 22.99, 11, 1)
Mg = Element("Mg", 24.31, 12, 2)
Al = Element("Al", 26.98, 13, None)
Si = Element("Si", 28.09, 14, None)
P = Element("P", 30.97, 15, -3)
S = Element("S", 32.06, 16, -2)
Cl = Element("Cl", 35.45, 17, -1)
Ar = Element("Ar", 39.95, 18, 0)
K = Element("K", 39.10, 19, 1)
Ca = Element("Ca", 40.08, 20, 2)
Sc = Element("Sc", 44.96, 21, None)
Ti = Element("Ti", 47.87, 22, None)
V = Element("V", 50.94, 23, None)
Cr = Element("Cr", 52.00, 24, 3)
Mn = Element("Mn", 54.94, 25, None)
Fe = Element("Fe", 55.85, 26, None)
Co = Element("Co", 58.93, 27, None)
Ni = Element("Ni", 58.69, 28, None)
Cu = Element("Cu", 63.55, 29, None)
Zn = Element("Zn", 65.38, 30, 2)
Ga = Element("Ga", 69.72, 31, None)
Ge = Element("Ge", 72.63, 32, None)
As = Element("As", 74.92, 33, -3)
Se = Element("Se", 78.97, 34, -2)
Br = Element("Br", 79.90, 35, -1)
Kr = Element("Kr", 83.80, 36, 0)
Rb = Element("Rb", 85.47, 37, 1)
Sr = Element("Sr", 87.62, 38, 2)
Y = Element("Y", 88.91, 39, None)
Zr = Element("Zr", 91.22, 40, None)
Nb = Element("Nb", 92.91, 41, None)
Mo = Element("Mo", 95.95, 42, None)
Tc = Element("Tc", 97, 43, None)
Ru = Element("Ru", 101.1, 44, None)
Rh = Element("Rh", 102.9, 45, None)
Pd = Element("Pd", 106.4, 46, None)
Ag = Element("Ag", 107.9, 47, 1)
Cd = Element("Cd", 112.4, 48, 2)
In = Element("In", 114.8, 49, None)
Sn = Element("Sn", 118.7, 50, None)
Sb = Element("Sb", 121.8, 51, -3)
Te = Element("Te", 127.6, 52, -2)
I = Element("I", 126.9, 53, -1)
Xe = Element("Xe", 131.3, 54, 0)
Cs = Element("Cs", 132.9, 55, 1)
Ba = Element("Ba", 137.3, 56, 2)
La = Element("La", 138.9, 57, None)
Ce = Element("Ce", 140.1, 58, None)
Pr = Element("Pr", 140.9, 59, None)
Nd = Element("Nd", 144.2, 60, None)
Pm = Element("Pm", 145, 61, None)
Sm = Element("Sm", 150.4, 62, None)
Eu = Element("Eu", 152.0, 63, None)
Gd = Element("Gd", 157.3, 64, None)
Tb = Element("Tb", 158.9, 65, None)
Dy = Element("Dy", 162.5, 66, None)
Ho = Element("Ho", 164.9, 67, None)
Er = Element("Er", 167.3, 68, None)
Tm = Element("Tm", 168.9, 69, None)
Yb = Element("Yb", 173.1, 70, None)
Lu = Element("Lu", 175.0, 71, None)
Hf = Element("Hf", 178.5, 72, None)
Ta = Element("Ta", 180.9, 73, None)
W = Element("W", 183.8, 74, None)
Re = Element("Re", 186.2, 75, None)
Os = Element("Os", 190.2, 76, None)
Ir = Element("Ir", 192.2, 77, None)
Pt = Element("Pt", 195.1, 78, None)
Au = Element("Au", 197.0, 79, None)
Hg = Element("Hg", 200.6, 80, None)
Tl = Element("Tl", 204.4, 81, None)
Pb = Element("Pb", 207.2, 82, None)
Bi = Element("Bi", 209.0, 83, -3)
Po = Element("Po", 209, 84, -2)
At = Element("At", 210, 85, -1)
Rn = Element("Rn", 222, 86, 0)
Fr = Element("Fr", 223, 87, 1)
Ra = Element("Ra", 226, 88, 2)
Ac = Element("Ac", 227, 89, None)
Th = Element("Th", 232.0, 90, None)
Pa = Element("Pa", 231.0, 91, None)
U = Element("U", 238.0, 92, None)
Np = Element("Np", 237, 93, None)
Pu = Element("Pu", 244, 94, None)
Am = Element("Am", 243, 95, None)
Cm = Element("Cm", 247, 96, None)
Bk = Element("Bk", 247, 97, None)
Cf = Element("Cf", 251, 98, None)
Es = Element("Es", 252, 99, None)
Fm = Element("Fm", 257, 100, None)
Md = Element("Md", 258, 101, None)
No = Element("No", 259, 102, None)
Lr = Element("Lr", 262, 103, None)
Rf = Element("Rf", 267, 104, None)
Db = Element("Db", 270, 105, None)
Sg = Element("Sg", 271, 106, None)
Bh = Element("Bh", 270, 107, None)
Hs = Element("Hs", 277, 108, None)
Mt = Element("Mt", 276, 109, None)
Ds = Element("Ds", 281, 110, None)
Rg = Element("Rg", 282, 111, None)
Cn = Element("Cn", 285, 112, None)
Nh = Element("Nh", 285, 113, None)
Fl = Element("Fl", 289, 114, None)
Mc = Element("Mc", 288, 115, None)
Lv = Element("Lv", 293, 116, None)
Ts = Element("Ts", 294, 117, -1)
Og = Element("Og", 294, 118, 0)

# Declare polyatomics
CH3CO2 = Polyatomic(C*H*H*H*C*O*O, -1)
C2H3O2 = Polyatomic(C*C*H*H*H*O*O, -1)
C9H7O4 = Polyatomic((C*9)*(H*7)*(O*4), -1)
N3 = Polyatomic(N*N*N, -1)
NH2 = Polyatomic(N*H*H, -1)
AtO3 = Polyatomic(At*O*O*O, -1)
N3 = Polyatomic(N*N*N, -1)
BiO3 = Polyatomic(Bi*O*O*O, -1)
HCO3 = Polyatomic(H*C*O*O*O, -1)#
HSO4 = Polyatomic(H*S*O*O*O*O, -1)
BrO3 = Polyatomic(Br*O*O*O*O, -1)
ClO4 = Polyatomic(Cl*O*O*O, -1)
ClO3 = Polyatomic(Cl*O*O*O, -1)
ClO2 = Polyatomic(Cl*O*O, -1)
ClO = Polyatomic(Cl*O, -1)
NO2 = Polyatomic(N*O*O, -1)
NO3 = Polyatomic(N*O*O*O, -1)
OCN = Polyatomic(O*C*N, -1)
CN = Polyatomic(C*N, -1)
#CHO3 = Polyatomic(C*H*O*O*O, -1)
OH = Polyatomic(O*H, -1)
IO3 = Polyatomic(I*O*O*O, -1)
C3H5O3 = Polyatomic((C*3)*(H*5)*(O*3), -1)
MnO3 = Polyatomic(Mn*O*O*O, -1)
MnO4 = Polyatomic(Mn*O*O*O, -1)
PO3 = Polyatomic(P*O*O*O, -1)
C7H5O3 = Polyatomic((C*7)*(H*5)*(O*3), -1)
C18H35O2 = Polyatomic((C*18)*(H*35)*(O*2), -1)
NH2SO3 = Polyatomic(N*H*H*S*O*O*O, -1)
#SCN = Polyatomic(S*C*N, -1)
VO3 = Polyatomic(V*O*O*O, -1)
H2PO4 = Polyatomic(H*H*P*O*O*O*O, -1)

#B4O7 = Polyatomic((B*4)*(O*7), -2)
CO3 = Polyatomic(C*O*O*O, -2)
CrO4 = Polyatomic(Cr*O*O*O*O, -2)
Cr2O7 = Polyatomic((Cr*2)*(O*7), -2)
C4H2O4 = Polyatomic((C*4)*(H*2)*(O*4), -2)
SiO3 = Polyatomic(Si*O*O*O, -2)
MoO4 = Polyatomic(Mo*O*O*O*O, -2)
C2O4 = Polyatomic(C*C*O*O*O*O, -2)
S2O8 = Polyatomic((S*2)*(O*8), -2)
SnO3 = Polyatomic(Sn*O*O*O, -2)
SeO4 = Polyatomic(Se*O*O*O*O, -2)
C4H6O4 = Polyatomic((C*4)*(H*6)*(O*4), -2)
SO4 = Polyatomic(S*O*O*O*O, -2)
SO3 = Polyatomic(S*O*O*O*O, -2)
C4H4O6 = Polyatomic((C*4)*(H*4)*(O*6), -2)
TeO4 = Polyatomic(Te*O*O*O*O, -2)
#SCN = Polyatomic(S*C*N, -2)
S2O3 = Polyatomic(S*S*O*O*O, -2)
WO4 = Polyatomic(W*O*O*O*O, -2)
HPO4 = Polyatomic(H*P*O*O*O*O, -2)

#AsO4 = Polyatomic(As*O*O*O*O, -3)
C6H8O7 = Polyatomic((C*6)*(H*8)*(O*7), -3)
PO4 = Polyatomic(P*O*O*O*O, -3)
#B4O7 = Polyatomic((B*4)*(O*7), -3)

SiO4 = Polyatomic(Si*O*O*O*O, -4)
#AsO4 = Polyatomic(As*O*O*O*O, -4)
P2O7 = Polyatomic((P*2)*(O*7), -4)

NH4 = Polyatomic(N*H*H*H*H, 1)
PH4 = Polyatomic(P*H*H*H*H, 1)
CHO = Polyatomic(C*H*O, 1)
C6H5 = Polyatomic((C*6)*(H*5), 1)

acids = [H*Cl, H*Br, H*I, H*ClO4, H*NO3, H*H*SO4,
         H*C2H3O2, H*C9H7O4, H*H*CO3, CHO*OH, H*CN, H*F, H*H*S,
         H*ClO, H*NO2, H*H*C2O4, C6H5*OH, H*H*H*PO4, H*H*SO3]

bases = [Li*OH, Na*OH, K*OH, Ca*(2*OH), Ba*(2*OH), Sr*(2*OH),
         N*H*H*H, C6H5*NH2, NH2*OH, NH2]

strong_acids = [H*Cl, H*Br, H*I, H*ClO4, H*NO3, H*H*SO4]
strong_bases = [Li*OH, Na*OH, K*OH, Ca*(2*OH), Ba*(2*OH), Sr*(2*OH)]

en = Charge(-1)
ep = Charge(1)

print("Use the function 'docs()' to print out all docstrings associated with all classes and their member functions")
print("An optional argument may be entered to narrow the results to a single class")
if DEBUG: print("[Debug]: {} Polyatomic ions loaded".format(len(Polyatomic.polyatomic_list)))
do_debug = True

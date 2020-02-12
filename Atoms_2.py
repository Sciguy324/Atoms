import itertools
from math import gcd
from math import log
from math import e
from collections import Counter
try:
    # Import parsing-related functions/classes
    from sympy.parsing.sympy_parser import parse_expr
    from sympy.core.numbers import Float as sympyFloat
    from sympy import solve as solve_expr
    from sympy import nsolve as nsolve_expr
    from sympy.physics import units as u
    #from sympy.abc import x
    from sympy import symbols
    INITIALIZED = True
except ModuleNotFoundError:
    print("WARNING: Unable to import 'SymPy' module.  Most algebra-solving functions will be unavailable.")
    INITIALIZED = False

DEBUG = False
do_debug = False

sub = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
sup = str.maketrans("0123456789+-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻")

def find_gcd(l):
    x = gcd(l[0], l[1])
    for i in range(2, len(l)):
        x = gcd(x, l[i])
    return x

def rate(molarity_a:list, molarity_b:list, rates:list):
    """Assist in guessing a rate law.  Takes two lists of concentration data and one list of rate data."""
    order_a = float(input("Order A guess: "))
    order_b = float(input("Order B guess: "))
    print("Resulting rate constants:")
    for i, (j, k, l) in enumerate(zip(molarity_a, molarity_b, rates)):
        print(i, "-", l / (j ** order_a * k ** order_b))

def conv_unit(x, unit1, unit2):
    """Utility function for pressure unit conversion"""
    units = {"atm": 1,
             "kPa": 101.3,
             "in Hg": 29.92,
             "mm Hg": 760,
             "torr": 760,
             "psi": 14.69,
             "Pa": 1.013e5
             }
    if unit1 not in units:
        raise ValueError("'{}' not defined in unit conversion dictionary.  Available units are: {}".format(unit1, ("".join([i+", " for i, j in units.items()])[:-2])))
    elif unit2 not in units:
        raise ValueError("'{}' not defined in unit conversion dictionary.  Available units are: {}".format(unit2, ("".join([i+", " for i, j in units.items()])[:-2])))
    else:
        return x / units[unit1] * units[unit2]

def unit_all(x, unit):
    """Utility function for performing all pressure unit conversions"""
    for i in ["atm", "kPa", "in Hg", "mm Hg", "torr", "psi", "Pa"]:
        print(conv_unit(x, unit, i))

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
                    print("\t{}{}: {}\n".format(k, inspect.signature(m), m.__doc__))
            print("")
    if search is None: print("\nExample of syntax:\n  Mg*ep*ep*aq + Ba*(2*NO3)*aq > Ba*ep*ep*aq + Mg*(2*NO3)*aq\nBecomes:\n", Mg*ep*ep*aq + Ba*(2*NO3)*aq > Ba*ep*ep*aq + Mg*(2*NO3)*aq)

class Solver:
    """Utility/wrapping class used when it is necessary to solve for a specific variable.
See bottom of file for declared solvers"""

    def __init__(self, variables: list, equation: str, units=None, hints=None):
        self.variables = variables
        
        # Check equation usability
        if '=' in equation:
            raise ValueError("Equation must be set up such that it equals zero!  No equal sign should be necessary")
        else:
            self.equation = equation

        # Load units property
        if units is None:
            # Default unit: None
            self.units = [1] * len(variables)
        else:
            # Hint exists, check for validity
            if len(units) == len(variables):
                # Hint list is the correct length
                self.units = units
            else:
                # Hint list is the incorrect length
                raise ValueError("Dimension mismatch between 'variables' list and 'units' list.")

        # Load hints property
        if hints is None:
            # Default hint: Nothing
            if units is None:
                self.hints = [""] * len(variables)
            else:
                self.hints = [" ({})".format(i) for i in units]
        else:
            # Hint exists, check for validity
            if len(hints) == len(variables):
                # Hint list is the correct length
                self.hints = [" ({})".format(i) for i in hints]
            else:
                # Hint list is the incorrect length
                raise ValueError("Dimension mismatch between 'variables' list and 'hints' list.")

    def __repr__(self):
        """Return a string representation of the equation"""
        return self.equation + ' = 0'

    def _run(self, values: dict):
        """Compute the equation using a set of inputs"""
        # TODO: Look into using .format instead of .replace
        if not INITIALIZED:
            raise ModuleNotFoundError("Cannot run solver without SymPy module!")
        # Substitute the values into the equation
        eq = self.equation
        for i, j in values.items():
            eq = eq.replace(i, str(j))
        print(eq)

        # Parse and solve equation
        result = solve_expr(parse_expr(eq))[0]
        if type(result) in [int, float, sympyFloat]:
            return result
        elif type(result) == dict:
            a = list(result.items())[0]
            return a[0] * solve_expr(a[1] - 1)[0]
        else:
            print("WARNING: Result of type '{}' was unexpected".format(type(result)))
            return result
    
    def run(self):
        """Wrapper for the '_run' function"""
        values = []
        unknown = ''
        unknown_hint = ''
        unknown_exists = False
        print("Use 'x' to indicate the SINGULAR unknown variable")
        # Ask user to input values for each variable
        for i, j, k in zip(self.variables, self.hints, self.units):
            user_input = eval(input("{}{}: ".format(i, j)))
            if type(user_input) in [int, float]:
                user_input *= k
            values.append("({})".format(user_input))
            if 'x' in str(user_input):
                # Register the unknown variable
                unknown = i
                unknown_hint = j
                unknown_exists = True
                
        if unknown_exists:
            # There was an unknown, proceed to solving
            print("{}{} = {}".format(unknown,
                                     unknown_hint,
                                     self._run(dict(zip(["__" + i for i in self.variables],
                                                        values)))))
        else:
            # I'm sorry Dave, but I'm afraid I can't do that
            raise ValueError("No unknown variable provided")
    

class State:
    """Utility class for setting the state of a molecule"""

    def __init__(self, symbol):
        self.symbol = symbol

    def __repr__(self):
        return "({})".format(self.symbol)

    def __mul__(self, x):
        """Multiply by element/molecule/polyatomic to set state"""
        if type(x) == Element:
            result = (1*x)
            result.state = self
            return result
        elif type(x) in [Molecule, Polyatomic]:
            result = x.copy()
            result.state = self
            return result
        else:
            raise TypeError("Unsupported operation between {} and {}".format(type(self), type(x)))

    def __rmul__(self, x):
        """Reverse multiplication"""
        return self.__mul__(x)
        

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
            # Invalid operation
            raise TypeError("Unsupported operation between {} and {}".format(type(self), type(x)))
    
    def __rmul__(self, x):
        """Reverse-order multiplication compatability that typically does nothing"""
        return self.__mul__(x)

class Element:
    """Class for basic elements"""

    elements = {}

    def __init__(self, symbol, molar_mass, atomic_number, oxidation, specific_heat=None):
        """Elements are referenced by their atomic symbol.  They may be multiplied together (or with numbers) to create molecules"""
        self.sym = symbol
        self.mol_mass = molar_mass
        self.atomic_number = atomic_number
        self.charge = 0
        self.neutrons = round(molar_mass) - atomic_number
        self.ox = oxidation
        self.specific_heat = specific_heat
        self.state = None
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

    def volume_to_moles(self, volume, temperature=273.15, pressure=1):
        """Converts volume to moles, assuming that the substance is a gas"""
        return volume * pressure / (temperature * R)

    def vol2mol(self, volume, temperature=273.15, pressure=1):
        """Alias for volume_to_moles"""
        return self.volume_to_moles(volume, temperature, pressure)

    def moles_to_volume(self, moles, temperature=273.15, pressure=1):
        """Converts moles to volume, assuming that the substancee is a gas.  Default temperature and pressure are at gas STP"""
        return moles * R * temperature / pressure

    def mol2vol(self, moles, temperature=273.15, pressure=1):
        """Alias for moles_to_volume"""
        return self.moles_to_volume(moles, temperature, pressure)

    # PV=nRT functions.  Not intended for user use
    @staticmethod
    def vnt_to_pressure(volume, moles, temperature):
        """Applies PV=nRT to obtain pressure"""
        return moles * R * temperature / volume
    
    @staticmethod
    def pnt_to_volume(pressure, moles, temperature):
        """Applies PV=nRT to obtain volume"""
        return moles * R * temperature / pressure
    
    @staticmethod
    def pvt_to_moles(pressure, volume, temperature):
        """Applies PV=nRT to obtain moles"""
        return pressure * volume / (R * temperature)
    
    @staticmethod
    def pvn_to_temp(pressure, volume, moles):
        """Applies PV=nRT to obtain temperature"""
        return pressure * volume / (moles * R)

    def _pvnrt(self, pressure, volume, moles, temperature):
        """Applies PV=nRT to determine the unknown variable, given all other information is known"""
        # Check if there are too many unknowns
        if sum([i is None for i in [pressure, volume, moles, temperature]]) > 1:
            raise ValueError("Too many unknowns")

        # Check if all information is known, and apply PV=nRT as a verification
        if sum([i is None for i in [pressure, volume, moles, temperature]]) == 0:
            return round(pressure*volume, 4) == round(moles*R*temperature, 4)

        # Apply PV=nRT for each unknown case
        args = [pressure, volume, moles, temperature]
        function = ([Element.vnt_to_pressure,
                     Element.pnt_to_volume,
                     Element.pvt_to_moles,
                     Element.pvn_to_temp])[args.index(None)]

        args.remove(None)
        
        return function(args[0], args[1], args[2])
    
    def pvnrt(self):
        """Wrapper for _pvnrt"""            
        args = [0,0,0,0]
        
        print("Note: Enter 'None' without units to indicate the unknown.  Values are assumed to be in atm, L, moles, and K respectively.  For other units, append the unit onto the end of the entry, separated by a space.")
        for i, j in enumerate(["Pressure", "Volume", "Moles", "Temperature"]):
            user = (input(j + ": ")).split(" ")

            # Check for unknown indicator
            if user[0] == "None":
                args[i] = None
                unknown_units = ("atm", "L", "mol", "K")[i]
                continue

            # Check for mercury-based units and correct
            if len(user) == 3:
                if user[2] == "Hg":
                    user[1] += " Hg"
                    user.pop(2)
                
            if len(user) == 2:
                # Check for pressure units
                if i == 0:
                    args[i] = conv_unit(floatify(user[0]), user[1], "atm")
                # Check for volume units
                elif i == 1:
                    try:
                        args[i] = floatify(user[0]) * ({"mL": 0.001, "milliliters": 0.001, "Milliliters": 0.001, "L": 1, "liters": 1, "Liters": 1})[user[1]]
                    except KeyError:
                        raise ValueError("Unrecognized unit {}".format(user[1]))
                # Check for mass units
                elif i == 2:
                    if user[1] in ["m", "mols", "moles", "mol", "mole"]:
                        args[i] = floatify(user[0])
                    elif user[1] in ["g", "grams", "gram"]:
                        args[i] = self.grams_to_moles(floatify(user[0]))
                    else:
                        raise ValueError("Unrecognized unit {}".format(user[1]))
                # Check for temperature units
                elif i == 3:
                    if user[1] in ["K", "kelvin", "Kelvin"]:
                        args[i] = floatify(user[0])
                    elif user[1] in ["C", "celsius", "Celsius"]:
                        args[i] = floatify(user[0]) + 273.15
                    elif user[1] in ["F", "fahrenheit", "Fahrenheit"]:
                        args[i] = (floatify(user[0]) - 32) * 5 / 9
                    else:
                        raise ValueError("Unrecognized unit {}".format(user[1]))
            elif len(user) > 2:
                raise ValueError("Invalid argument {}".format(" ".join(user)))
            else:
                args[i] = float(user[0])

        print("Result:")
        print(self._pvnrt(args[0], args[1], args[2], args[3]), unknown_units)
    
    def gas_density(self, temperature=273.15, pressure=1):
        """Calculates the density of the substance, assuming that it is a gas.  Default temperature and pressure are at gas STP"""
        return self.molar_mass() * pressure / (R * temperature)

    def heat(self, mass, temperature_change):
        """Determines the energy change after enacting a temperature change on a mass of an element"""
        return mass * self.specific_heat * temperature_change

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

        elif type(x) == State:
            return x * self
        
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
    
    def __init__(self, parts_dict, count=1, charge=0, state=None):
        self.parts = parts_dict
        self.count = count
        self.mol_mass = count * sum([j*(i.molar_mass()) for i, j in parts_dict.items()])
        self.count = count
        self.charge = charge
        self.state = state

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
        
        # First check if the thing is an element that has been cast as a 'molecule'
        if type(x) == Molecule and len(x.parts) == 1 and tuple(x.parts.items())[0][1] == 1:
            x = tuple(x.parts.items())[0][0]        

        # Add thing to molecule dictionary by the amount previously extracted
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

        if self.state is not None:
            result += str(self.state)
            
        return result

    def copy(self):
        """Returns a new instance of the molecule identical to self."""
        return Molecule(dict(self.parts), count=self.count, charge=self.charge, state=self.state)

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

    def volume_to_moles(self, volume, temperature=273.15, pressure=1):
        """Converts volume to moles, assuming that the molecule is a gas"""
        return volume * pressure / (temperature * R)

    def vol2mol(self, volume, temperature=273.15, pressure=1):
        """Alias for volume_to_moles"""
        return self.volume_to_moles(volume, temperature, pressure)

    def moles_to_volume(self, moles, temperature=273.15, pressure=1):
        """Converts moles to volume, assuming that the molecule is a gas.  Default temperature and pressure are at gas STP"""
        return moles * R * temperature / pressure

    def mol2vol(self, moles, temperature=273.15, pressure=1):
        """Alias for moles_to_volume"""
        return self.moles_to_volume(moles, temperature, pressure)

    def _pvnrt(self, pressure, volume, moles, temperature):
        """Applies PV=nRT to solve for an unknown"""
        return Element._pvnrt(self, pressure, volume, moles, temperature)

    def pvnrt(self):
        """Wrapper for _pvnrt"""
        return Element.pvnrt(self)

    def gas_density(self, temperature=273.15, pressure=1):
        """Calculates the density of the molecule, assuming that it is a gas.  Default temperature and pressure are at gas STP"""
        return self.molar_mass() * pressure / (R * temperature)

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
            if round(j) != j:
                raise Exception("An unknown error occured while attempting to assign charges: non-integer charge")
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

        # Check for special cases otherwise not covered by the rules
        if self == ((C*3)*(H*8)):
            return {H: 0, C: 0}
        
        # Apply rules for elemental atoms and singular ions
        ion_list = [i for i, j in self.parts.items()]
        
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
            # Molecule * charge
            result = self.copy()
            result.charge += x.charge
            return result
        elif type(x) == State:
            # Molecule * State
            return x * self
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

    def __init__(self, molec, oxidation, count=1, state=None):
        """Polyatomics are called by entering a single term.  For example, C₂H₃O₂⁻ is referenced as C2H3O2"""
        if type(molec) == dict:
            self.parts = dict(molec)
        else:
            self.parts = dict(molec.parts)
        self.mol_mass = sum([j*i.molar_mass() for i, j in self.parts.items()])
        self.ox = oxidation
        self.count = count
        self.state = state
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

        if self.state is not None:
            result += str(self.state)
            
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
        return Polyatomic(dict(self.parts), self.oxidation(), count=self.count, state=self.state)

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
        elif type(x) == State:
            # Polyatomic * State
            return x * self
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
    raoult_solver = Solver(['Pₛₒₗₙ', 'Pₛₒₗᵥ', 'Molesₛₒₗᵥ', 'Molesₛₒₗᵤₜₑ'],
                           'Pₛₒₗᵥ * Molesₛₒₗᵥ / (Molesₛₒₗᵥ + Molesₛₒₗᵤₜₑ) - Pₛₒₗₙ',
                           hints=['atm', 'atm', '', ''])

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
        # Set up component-charge dictionary, sorted by molar mass
        parts = list(dict(self.parts).items())
        parts.sort(key=lambda x: x[0].molar_mass())

        # Assign oxidation states to the components
        for i, j in parts:
            copy = i.copy()
            copy = j * copy
            stage1.append(copy._assign_oxidation())
        
        stage2 = []
        stage2.append([])

        # Condense data by eliminating repeats
        for i in stage1:
            for j, k in i.items():
                for m, n in enumerate(stage2):
                    if [j, k] not in n:
                        stage2[m].append([j, k])

        # Convert to a list of dictionaries.  List format is present to allow same
        # atoms to have multiple oxiditation states on one side
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

    def equil_temp_f(self):
        """Calculates the final temperature given a series of initial conditions"""
        numerator = 0
        denominator = 0
        print("Please enter mass in grams")
        cal_const = float(input("Calorimeter Constant: "))
        cal_temp_i = float(input("Calorimeter Initial Temp - "))
        numerator += cal_const * cal_temp_i
        denominator += cal_const
        
        for i in self.parts:
            mass = float(input(str(i.singular()) + " Mass - "))
            sp_heat = float(input(str(i.singular()) + " Specific Heat - "))
            temp_i = float(input(str(i.singular()) + " Initial Temp - "))
            numerator += mass * sp_heat * temp_i
            denominator += mass * sp_heat

        return numerator / denominator

    def equil_sp_heat(self):
        result = 0
        print("Enter 'None' to indicate the unknown specific heat")
        final_t = float(input("Final Temp - "))
        
        cal_const = float(input("Calorimeter Constant: "))
        cal_temp_i = float(input("Calorimeter Initial Temp - "))
        result += cal_const * (final_t - cal_temp_i)
        
        for i in self.parts:
            mass = float(input(str(i.singular()) + " Mass - "))
            temp_i = float(input(str(i.singular()) + " Initial Temp - "))
            sp_heat = input(str(i.singular()) + " Specific Heat - ")
            if sp_heat == "None":
                mass_u = mass
                temp_i_u = temp_i
            else:
                sp_heat = float(sp_heat)
                result += mass * sp_heat * (final_t - temp_i)

        return result / (-1 * mass_u * (final_t - temp_i_u))

    def raoult(self):
        if len(self.parts) != 2:
            raise ValueError("Invalid number of molecules present in combination, this case has not been coded for!")
        # First molecule is the dissolved one
        if list(self.parts.items())[0][0].state == aq:
            Combination.raoult_solver.run()
        # Second molecule is the dissolved one
        elif list(self.parts.items())[1][0].state == aq:
            pass
        else:
            raise ValueError("Dissolved molecules must be denoted with an 'aq' state (multiply molecule by aq)")
        
            
    

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
            print("Warning: Noninteger charge detected")
            return False
        if round(right_charge) != right_charge:
            print("Warning: Noninteger charge detected")
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

    def solve(self, limit=10, ignore_charge=False):
        """Auto balances a chemcial reaction.  Will continue until reaching the maximum number of a molecule allowed."""
        # First check if the two sides have the same elements
        if not self.left.compare(self.right):
            raise ValueError("Reaction impossible to balance, elements do not match:\nLeft: {}\nRight:{}".format([i for i, j in self.left.get_atoms().items()],
                                                                                                                 [i for i, j in self.right.get_atoms().items()]))
        # Variables to track what type of verification cause problems
        element_imbalance = 0
        charge_imbalance = 0

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
                if possible.verify():
                    if (ignore_charge or possible.verify_charge()):
                        if DEBUG:
                            print("Element Check Fails:", element_imbalance)
                            print("Charge Check Fails:", charge_imbalance)
                        return possible
                    else:
                        # Save a backup in case the solution was actually valid, and no other solutions were found.
                        backup = possible.copy()
                        charge_imbalance += 1
                else:
                    element_imbalance += 1
        

        print("Element Check Fails:", element_imbalance)
        print("Charge Check Fails:", charge_imbalance)

        # If only one charge imbalance was detected, ask if that was the solution the user was looking for
        if charge_imbalance == 1:
            user = input("""No solutions were found within the limit, but a charge imbalance was detected in the closest solution: {}.
Take this as the solution? (Y/N): """.format(backup))
            if user.lower() == "y":
                return backup

        # If multiple charge imbalances were detected, ask the user if they want to re-solve, this time ignoring charges
        elif charge_imbalance > 1:
            user = input("""No solutions were found within the limit, but charge imbalances were detected.
Run the solver again without regard to charges? (Y/N): """)
            if user.lower() == "y":
                possible = self.solve(limit=limit, ignore_charge=True)
        
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

##    def equil(self):
##        """Calculated the equilibrium constant of a given reaction"""
##        print("Input equilibrium concentrations")
##        
##        denominator = 1
##        for i, j in self.left.parts.items():
##            if i.state not in [l, g]:
##                denominator *= eval(input("[{}]: ".format(i))) ** j
##                
##        numerator = 1
##        for i, j in self.right.parts.items():
##            if i.state not in [l, g]:
##                numerator *= eval(input("[{}]: ".format(i))) ** j
##
##        return numerator / denominator

    def equil(self):
        """General equation solver for equilibrium coefficients"""
        variables = ["k"]
        products = "1"
        for i, j in self.right.parts.items():
            if i.state not in [l, s]:
                variables.append("[{}]".format(i))
                products += "*(__[{}])**{}".format(i, j)

        reactants = "1"
        for i, j in self.left.parts.items():
            if i.state not in [l, s]:
                variables.append("[{}]".format(i))
                reactants += "*(__[{}])**{}".format(i, j)
            
        general = "({}) / ({}) - __k".format(products, reactants)

        return Solver(variables, general).run()
        

# Declare elements
H = Element("H", 1.008, 1, 1, None)
He = Element("He", 4.003, 2, 0, 5.19)
Li = Element("Li", 6.94, 3, 1, None)
Be = Element("Be", 9.012, 4, 2, None)
B = Element("B", 10.81, 5, None, None)
C = Element("C", 12.01, 6, -4, None)
N = Element("N", 14.01, 7, -3, None)
O = Element("O", 16.00, 8, -2, None)
F = Element("F", 19.00, 9, -1, None)
Ne = Element("Ne", 20.18, 10, 0, 1.03)
Na = Element("Na", 22.99, 11, 1, None)
Mg = Element("Mg", 24.31, 12, 2, 1.017)
Al = Element("Al", 26.98, 13, None, 0.9)
Si = Element("Si", 28.09, 14, None, 0.703)
P = Element("P", 30.97, 15, -3, None)
S = Element("S", 32.06, 16, -2, 0.732)
Cl = Element("Cl", 35.45, 17, -1, None)
Ar = Element("Ar", 39.95, 18, 0, None)
K = Element("K", 39.10, 19, 1, None)
Ca = Element("Ca", 40.08, 20, 2, None)
Sc = Element("Sc", 44.96, 21, None, None)
Ti = Element("Ti", 47.87, 22, None, 0.523)
V = Element("V", 50.94, 23, None, None)
Cr = Element("Cr", 52.00, 24, 3, 0.448)
Mn = Element("Mn", 54.94, 25, None, None)
Fe = Element("Fe", 55.85, 26, None, 0.444)
Co = Element("Co", 58.93, 27, None, None)
Ni = Element("Ni", 58.69, 28, None, 0.444)
Cu = Element("Cu", 63.55, 29, None, 0.385)
Zn = Element("Zn", 65.38, 30, 2, 0.388)
Ga = Element("Ga", 69.72, 31, None, None)
Ge = Element("Ge", 72.63, 32, None, None)
As = Element("As", 74.92, 33, -3, None)
Se = Element("Se", 78.97, 34, -2, None)
Br = Element("Br", 79.90, 35, -1, None)
Kr = Element("Kr", 83.80, 36, 0, 0.247)
Rb = Element("Rb", 85.47, 37, 1, None)
Sr = Element("Sr", 87.62, 38, 2, None)
Y = Element("Y", 88.91, 39, None, None)
Zr = Element("Zr", 91.22, 40, None, None)
Nb = Element("Nb", 92.91, 41, None, None)
Mo = Element("Mo", 95.95, 42, None, None)
Tc = Element("Tc", 97, 43, None, None)
Ru = Element("Ru", 101.1, 44, None, None)
Rh = Element("Rh", 102.9, 45, None, None)
Pd = Element("Pd", 106.4, 46, None, None)
Ag = Element("Ag", 107.9, 47, 1, 0.237)
Cd = Element("Cd", 112.4, 48, 2, 0.232)
In = Element("In", 114.8, 49, None, None)
Sn = Element("Sn", 118.7, 50, None, 0.213)
Sb = Element("Sb", 121.8, 51, -3, None)
Te = Element("Te", 127.6, 52, -2, None)
I = Element("I", 126.9, 53, -1, None)
Xe = Element("Xe", 131.3, 54, 0, 0.158)
Cs = Element("Cs", 132.9, 55, 1, None)
Ba = Element("Ba", 137.3, 56, 2, None)
La = Element("La", 138.9, 57, None, None)
Ce = Element("Ce", 140.1, 58, None, None)
Pr = Element("Pr", 140.9, 59, None, None)
Nd = Element("Nd", 144.2, 60, None, None)
Pm = Element("Pm", 145, 61, None, None)
Sm = Element("Sm", 150.4, 62, None, None)
Eu = Element("Eu", 152.0, 63, None, None)
Gd = Element("Gd", 157.3, 64, None, None)
Tb = Element("Tb", 158.9, 65, None, None)
Dy = Element("Dy", 162.5, 66, None, None)
Ho = Element("Ho", 164.9, 67, None, None)
Er = Element("Er", 167.3, 68, None, None)
Tm = Element("Tm", 168.9, 69, None, None)
Yb = Element("Yb", 173.1, 70, None, None)
Lu = Element("Lu", 175.0, 71, None, None)
Hf = Element("Hf", 178.5, 72, None, None)
Ta = Element("Ta", 180.9, 73, None, None)
W = Element("W", 183.8, 74, None, 0.133)
Re = Element("Re", 186.2, 75, None, None)
Os = Element("Os", 190.2, 76, None, None)
Ir = Element("Ir", 192.2, 77, None, None)
Pt = Element("Pt", 195.1, 78, None, 0.133)
Au = Element("Au", 197.0, 79, None, 0.129)
Hg = Element("Hg", 200.6, 80, None, 0.138)
Tl = Element("Tl", 204.4, 81, None, None)
Pb = Element("Pb", 207.2, 82, None, 0.159)
Bi = Element("Bi", 209.0, 83, -3, None)
Po = Element("Po", 209, 84, -2, None)
At = Element("At", 210, 85, -1, None)
Rn = Element("Rn", 222, 86, 0, None)
Fr = Element("Fr", 223, 87, 1, None)
Ra = Element("Ra", 226, 88, 2, None)
Ac = Element("Ac", 227, 89, None, None)
Th = Element("Th", 232.0, 90, None, None)
Pa = Element("Pa", 231.0, 91, None, None)
U = Element("U", 238.0, 92, None, 0.115)
Np = Element("Np", 237, 93, None, None)
Pu = Element("Pu", 244, 94, None, None)
Am = Element("Am", 243, 95, None, None)
Cm = Element("Cm", 247, 96, None, None)
Bk = Element("Bk", 247, 97, None, None)
Cf = Element("Cf", 251, 98, None, None)
Es = Element("Es", 252, 99, None, None)
Fm = Element("Fm", 257, 100, None, None)
Md = Element("Md", 258, 101, None, None)
No = Element("No", 259, 102, None, None)
Lr = Element("Lr", 262, 103, None, None)
Rf = Element("Rf", 267, 104, None, None)
Db = Element("Db", 270, 105, None, None)
Sg = Element("Sg", 271, 106, None, None)
Bh = Element("Bh", 270, 107, None, None)
Hs = Element("Hs", 277, 108, None, None)
Mt = Element("Mt", 276, 109, None, None)
Ds = Element("Ds", 281, 110, None, None)
Rg = Element("Rg", 282, 111, None, None)
Cn = Element("Cn", 285, 112, None, None)
Nh = Element("Nh", 285, 113, None, None)
Fl = Element("Fl", 289, 114, None, None)
Mc = Element("Mc", 288, 115, None, None)
Lv = Element("Lv", 293, 116, None, None)
Ts = Element("Ts", 294, 117, -1, None)
Og = Element("Og", 294, 118, 0, None)

# Declare unknown element/molecule
X = Element("???", 0, 0, None, None)
X2 = X*X

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
BrO3 = Polyatomic(Br*O*O*O, -1)
ClO4 = Polyatomic(Cl*O*O*O*O, -1)
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
MnO4 = Polyatomic(Mn*O*O*O*O, -1)
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
SO3 = Polyatomic(S*O*O*O, -2)
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

# Declare charges utility
en = Charge(-1)
ep = Charge(1)

# Declare states utility
s = State("s")
l = State("l")
g = State("g")
aq = State("aq")

# Declare catalogues of acids/bases
acids = [H*Cl, H*Br, H*I, H*ClO4, H*NO3, H*H*SO4,
         H*C2H3O2, H*C9H7O4, H*H*CO3, CHO*OH, H*CN, H*F, H*H*S,
         H*ClO, H*NO2, H*H*C2O4, C6H5*OH, H*H*H*PO4, H*H*SO3]

bases = [Li*OH, Na*OH, K*OH, Ca*(2*OH), Ba*(2*OH), Sr*(2*OH),
         N*H*H*H, C6H5*NH2, NH2*OH, NH2]

strong_acids = [H*Cl, H*Br, H*I, H*ClO4, H*NO3, H*H*SO4]
strong_bases = [Li*OH, Na*OH, K*OH, Ca*(2*OH), Ba*(2*OH), Sr*(2*OH)]

# Declare melting point data catalogues
# Heats/enthalpy of fusion are in kJ/mol, melting points are in °C
heat_of_fusion = {H*H*O: 6.007, (C*3)*(H*8)*O: 5.37, (C*3)*(H*6)*O: 5.69, (C*4)*(H*10)*O: 7.27,
                  N*(H*3): 5.65, H*Cl: 1.992, C*O: 0.836, C*(Cl*4): 2.5, Na*Cl: 28.8, Ne: 0.33,
                  O*2: 0.44, C*(H*4): 0.94, C*C*(H*6): 2.85, Cl*2: 6.4, C*(Cl*4): 2.67,
                  (C*9)*(H*20): 19.3, Hg: 2.3, Na: 2.6, Al: 10.9}
melting_point = {H*H*O: 0.0, (C*3)*(H*8)*O: -89.5, (C*3)*(H*6)*O: -94.8, (C*4)*(H*10)*O: -116.3,
                 Ne: -249.15, O*2: -219.15, C*(H*4): -182.45, C*C*(H*6): -183.15, Cl*2: -100.95,
                 C*(Cl*4): -23.15, (C*9)*(H*20): 79.85, Hg: -39.15, Na: 97.85, Al: 659.85}

# Declare boiling point data catalogues
heat_of_vaporization = {Ne: 1.8, O*2: 6.82, C*(H*4): 8.18, C*C*(H*6): 14.72, Cl*2: 20.41,
                        C*(Cl*4): 30.0, H*H*O: 40.7, (C*9)*(H*20): 40.5, Hg: 58.6, Na: 98, Al: 284}
boiling_point = {Ne: -246.15, O*2: -182.95, C*(H*4): -161.15, C*C*(H*6): -89.15, Cl*2: -34.15,
                 C*(Cl*4): 76.85, H*H*O: 99.95, (C*9)*(H*20): 217.85, Hg: 356.85, Na: 884.85, Al: 2326.85}

# Declare some constants
R = 0.082058 * u.l * u.atm / u.mol / u.K # Ideal gas constant in L atm / (mol K)
R2 = 8.314 * u.J / u.K / u.mol

# Declare some basic solvers
raoult_solver = Solver(['Pₛₒₗₙ', 'Pₛₒₗᵥ', 'Molesₛₒₗᵥ', 'Molesₛₒₗᵤₜₑ'],
                        '__Pₛₒₗᵥ * __Molesₛₒₗᵥ / (__Molesₛₒₗᵥ + __Molesₛₒₗᵤₜₑ) - __Pₛₒₗₙ',
                        units=[u.atm, u.atm, u.moles, u.moles],
                        hints=['atm', 'atm', '', ''])

freeze_depress = Solver(['ΔT_f', 'Molesₛₒₗᵤₜₑ', 'Massₛₒₗᵥ', 'i', 'K_f'],
                        '__Molesₛₒₗᵤₜₑ / __Massₛₒₗᵥ * __i * __K_f - __ΔT_f',
                        units=[u.kelvin, u.moles, u.kg, 1, u.kelvin * u.kg / u.moles],
                        hints=['K', 'moles', 'kg', '', 'K*kg/mol'])

boiling_elevate = Solver(['ΔT_b', 'Molesₛₒₗᵤₜₑ', 'Massₛₒₗᵥ', 'i', 'K_b'],
                         '__Molesₛₒₗᵤₜₑ / __Massₛₒₗᵥ * __i * __K_b - __ΔT_b',
                         units=[u.kelvin, u.moles, u.kg, 1, u.kelvin * u.kg / u.moles],
                         hints=['K', 'moles', 'kg', '', 'K*kg/mol'])

arrhenius = Solver(['ln(k)', 'ln(A)', 'Eₐ/R', 'T'],
                   '-__ln(k) + __ln(A) - __Eₐ/R / __T')

zero_order_rate = Solver(['[A]', '[A₀]', 'k', 't'],
                         '-__k * __t + __[A₀] - __[A]')

first_order_rate = Solver(['ln([A])', 'ln([A]₀)', 'k', 't'],
                          '-__ln([A]) + __ln([A]₀) - __k * __t')

second_order_rate = Solver(['1/[A]', '1/[A₀]', 'k', 't'],
                           '__k * __t + __1/[A₀] - __1/[A]')

equil_const_convert = Solver(['Kₚ', 'K_c', 'T', 'Δn'],
                             '__K_c * (0.082058 * __T) ** __Δn - __Kₚ')

# Declare global 'x' variable
x=symbols('x', real=True)

print("Use the function 'docs()' to print out all docstrings associated with all classes and their member functions")
print("An optional argument may be entered to narrow the results to a single class")
if DEBUG: print("[Debug]: {} Polyatomic ions loaded".format(len(Polyatomic.polyatomic_list)))
do_debug = True

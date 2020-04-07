# But it started out ok!

import itertools
from math import gcd
from math import log
from math import log10
from math import e
from collections import Counter
try:
    # Import parsing-related functions/classes
    from sympy.parsing.sympy_parser import parse_expr
    from sympy.core.numbers import Float as sympyFloat
    from sympy.core.add import Add as sympyAdd
    from sympy import solve as solve_expr
    from sympy import nsolve as nsolve_expr
    from sympy.physics import units as u
    from sympy import symbols
    INITIALIZED = True
except ModuleNotFoundError:
    print("WARNING: Unable to import 'SymPy' module.  Most algebra-solving functions will be unavailable.")
    INITIALIZED = False

DEBUG = False
do_debug = False

sub = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
sup = str.maketrans("0123456789+-.", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻ᐧ")

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

    def __call__(self):
        """Allows solvers to be called with the syntax foo()"""
        self.run()

    def _run(self, values: dict, return_all=False):
        """Compute the equation using a set of inputs"""
        if not INITIALIZED:
            raise ModuleNotFoundError("Cannot run solver without SymPy module!")
        # Substitute the values into the equation
        eq = self.equation
        for i, j in values.items():
            eq = eq.replace(i, str(j))
        print(eq)

        # Parse and solve equation
        if return_all:
            return solve_expr(parse_expr(eq))
        else:
            result = solve_expr(parse_expr(eq))
            if type(result) == list and len(result) == 1:
                result = result[0]
        if type(result) in [int, float, sympyFloat]:
            return result
        elif type(result) == dict:
            print(result)
            a = list(result.items())[0]
            return a[0] * solve_expr(a[1] - 1)[0]
        elif type(result) == list:
            for i in result:
                try:
                    a = list(i.items())[0]
                    return a[0] * solve_expr(a[1] - 1)[0]
                except IndexError:
                    pass
                except AttributeError:
                    return result
            print("An error occured while attempting to find the solution\n\t", result)
        else:
            print("WARNING: Result of type '{}' was unexpected".format(type(result)))
            return result
    
    def run(self, return_all=False, return_result=False):
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
            result = self._run(dict(zip(["__" + i for i in self.variables], values)), return_all)
            if return_result:
                return result
            else:                
                print("{}{} = {}".format(unknown,
                                         unknown_hint,
                                         result))
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

    @property
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

    @property
    def symbol(self):
        """Return the element's symbol for string representation"""
        return str(self.sym)

    @property
    def enthalpy(self):
        """Uses a lookup table to find the element's molar enthalpy"""
        if str(self) in molar_enthalpy:
            return molar_enthalpy[str(self)] * u.J * 1000 / u.moles
        else:
            if self.state is None:
                raise AttributeError("The ΔH⁰_f for '{}' could not be found in the molar_enthalpy lookup table (Hint: make sure to specify state!)".formate(self))
            else:
                raise AttributeError("The ΔH⁰_f for '{}' could not be found in the molar_enthalpy lookup table".format(self))

    @property
    def entropy(self):
        """Uses a lookup table to find the element's standard entropy"""
        if str(self) in standard_entropy:
            return standard_entropy[str(self)] * u.J / u.moles / u.K
        else:
            if self.state is None:
                raise AttributeError("The S⁰ for '{}' could not be found in the standard_entropy lookup table (Hint: make sure to specify state!)".format(self))
            else:
                raise AttributeError("The S⁰ for '{}' could not be found in the standard_entropy lookup table".format(self))

    def grams_to_moles(self, grams):
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
        return self.molar_mass * pressure / (R * temperature)

    def heat(self, mass, temperature_change):
        """Determines the energy change after enacting a temperature change on a mass of an element"""
        return mass * self.specific_heat * temperature_change

    def _remove_atom(self, x):
        """Handle case where an atom is attempted to be removed from an atom"""
        raise ValueError("It is not possible to remove an atom ({}) from an atom ({})!".format(x, self))

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
        self.mol_mass = count * sum([j*(i.molar_mass) for i, j in parts_dict.items()])
        self.count = count
        self.charge = charge
        self.state = state

    def __hash__(self):
        return hash((tuple(self.get_atoms().items()), self.count))

    def __eq__(self, x):
        return self.compare(x, ignore_num=False)

    def __repr__(self):
        return self.symbol

    def __contains__(self, key):
        return key in self.parts

    def extract_count(self):
        """Sets the count to 1 and returns what it was before"""
        count = self.count
        self.count = 1
        return count

    def append(self, x):
        """Appends something onto the molecule"""
        amount = x.extract_count()
        
        # First check if the 'x' is an element that has been cast as a 'molecule'
        if type(x) == Molecule and len(x.parts) == 1 and tuple(x.parts.items())[0][1] == 1:
            x = tuple(x.parts.items())[0][0]        

        # Add 'x' to molecule dictionary by the amount previously extracted
        if x in self:
            self.parts[x] += amount
        else:
            self.parts[x] = amount

        # Modify charge
        #self.charge += x.charge * amount

    def _remove_atom(self, x):
        """Removes something from the molecule.  Warning: for internal use only due too finiky-ness"""
        # TODO: I don't like how this handles removing the atom from components
        # First check if the 'x' is an Element that has been cast as a 'molecule'
        if type(x) == Molecule and len(x.parts) == 1 and tuple(x.parts.items())[0][1] == 1:
            x = tuple(x.parts.items())[0][0]

        if type(x) != Element:
            raise ValueError("Cannot remove non-element '{}' from '{}'".format(x, self))
        
        # Remove one 'x' from molecule dictionary, if possible
        if x in self:
            if self.parts[x] > 1:
                self.parts[x] -= 1
            elif self.parts[x] == 1:
                del self.parts[x]
            else:
                raise ValueError("An unexpected error occured while attempting to remove '{}' from '{}'".format(x, self))
        else:
            # 'x' was not found in root dictionary, check components
            new_parts = {}
            skip_checking = False
            for i, j in self.parts.items():
                # The component was already removed, skip checking
                if skip_checking:
                    new_parts[i.copy()] = j
                
                # Cancel scanning this component if there's more than one of it
                if j != 1:
                    new_parts[i.copy()] = j
                    continue
                try:
                    # Attempt to remove from component
                    copy = i.copy()
                    if type(copy) is Polyatomic:
                        copy = Molecule(copy.parts, 1, 0)
                    copy._remove_atom(x)
                    new_parts[copy] = 1
                    
                except ValueError:
                    # Not found in component
                    new_parts[i.copy()] = 1
                    continue
                else:
                    # Successly removed from a component
                    skip_checking = True

            if skip_checking:
                self.parts = new_parts
            else:
                # 'x' was not found in components
                raise ValueError("Could not remove '{}' from '{}'".format(x, self))

    @property
    def molar_mass(self):
        """Calculate the molar mass of the element"""
        self.mol_mass = self.count * sum([j *(i.molar_mass) for i, j in self.parts.items()])
        return self.mol_mass

    def oxidation(self):
        """Get the oxidation state/charge of the thing"""
        return self.charge

    @property
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
                result += i.symbol
                if j > 1:
                    result += str(j).translate(sub)
            else:
                if j > 1:
                    result += "(" + i.symbol + ")" + str(j).translate(sub)
                else:
                    result += i.symbol

        if self.charge not in [0, None]:
            result += "]"
            # Prevent unnecessary decimals in charges
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

    @property
    def enthalpy(self):
        """Uses a lookup table to find the molecule's molar enthalpy"""
        if str(self) in molar_enthalpy:
            return molar_enthalpy[str(self)] * u.J * 1000 / u.moles
        else:
            if self.state is None:
                raise AttributeError("The ΔH⁰_f for '{}' could not be found in the molar_enthalpy lookup table (Hint: make sure to specify state!)".formate(self))
            else:
                raise AttributeError("The ΔH⁰_f for '{}' could not be found in the molar_enthalpy lookup table".format(self))

    @property
    def entropy(self):
        """Uses a lookup table to find the molecule's standard entropy"""
        if str(self) in standard_entropy:
            return standard_entropy[str(self)] * u.J / u.moles / u.K
        else:
            if self.state is None:
                raise AttributeError("The S⁰ for '{}' could not be found in the standard_entropy lookup table (Hint: make sure to specify state!)".format(self))
            else:
                raise AttributeError("The S⁰ for '{}' could not be found in the standard_entropy lookup table".format(self))

    def copy(self):
        """Returns a new instance of the molecule identical to self."""
        return Molecule(dict(self.parts), count=self.count, charge=self.charge, state=self.state)

    def grams_to_moles(self, grams, element=None):
        """Converts grams of a molecule to moles.  If an element is specified, the number of
        moles of a particular element will be isolated and returned"""
        if element:
            return self.parts[element] * grams / self.molar_mass
        else:
            return grams / self.molar_mass

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
            return moles * (Molecule({element: self.parts[element]})).molar_mass
        else:
            return moles * self.molar_mass

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
        return self.molar_mass * pressure / (R * temperature)

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
        """Converts a monotomic bonding coefficient to a quantity coefficient: H₂→2H"""
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
    
    def dissolve(self, ignore_ph=True, override_check=False):
        """Dissociates the ions of a molecule, if it's soluble"""
        if not override_check:
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

    def dissolve_salt(self, ignore_ph=True, override_check=True):
        """Similar to .dissolve(), but returns the full dissolution reaction with the salt.
Why so salty?"""
        copy = self.copy()
        return copy*s > copy.dissolve(ignore_ph, override_check)

    def solubility_solver(self, dissolved:dict={}, approx=False):
        """General solver for the solubility of a salt"""
        if type(dissolved) is not dict:
            raise TypeError("Dictionary of dissolved molecules + molarity is not a dictionary")
        else:
            if len(dissolved) > 0 and not approx:
                user_input = input("Warning: There is a chance that the solver will not be able to solve this equation symbolically.  Would you like to approximate (y/n)?\n")
                if user_input.lower() in ["yes", "y", "t", "true"]:
                    approx = True
                else:
                    approx = False

            ions = {}
            for i, j in dissolved.items():
                # i: Molecule, j: Concentration
                for k, l in i.parts.items():
                    # k: Part, l: Part count
                    if k in ions:
                        ions[k.copy()] += l * j
                    else:
                        ions[k.copy()] = l * j
                    
        
        if len(self.parts) != 2:
            raise ValueError("Unable to dissolve molecule - This script doesn't know how")
        variables = ["Solubility", "kₛₚ"]
        equation = "1"
        for i, j in self.parts.items():
            if i in ions:
                if approx:
                    equation += "*({})**{}".format(ions[i], j)
                else:
                    equation += "*({}*__Solubility + {})**{}".format(j, ions[i], j)
            else:
                equation += "*({}*__Solubility)**{}".format(j, j)
        equation += "-__kₛₚ"
        result = Solver(variables, equation).run(return_result=True)
        if type(result) is list:
            for i in result.copy():
                if type(i) is sympyAdd or i < 0:
                    result.remove(i)
            return result
        else:
            return result
        

    def _percents(self):
        """Get the percent mass of each element in a molecule"""
        factor = 100 / self.molar_mass
        percents_dict = {}

        for i, j in self.parts.items():
            if i not in percents_dict:
                percents_dict[i] = self.parts[i] * i.molar_mass * factor

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
            raise ValueError("Error: insufficient information to assign oxidation states")

    def assign_oxidation(self):
        """Wrapper for _assign_oxidation to print results in a more readable format"""
        for i, j in self._assign_oxidation().items():
            print(i, ":", j)

    def conjugate_acid(self):
        """Returns the conjugate acid of the molecule"""
        return self * H * ep

    def conjugate_base(self):
        """Returns the conjugate base of the molecule"""
        result = self.copy()
        result._remove_atom(H)
        return result * en

    def acid_reaction(self):
        """Determines the reaction in which the molecule acts as an acid"""
        return (self + H*H*O > self.conjugate_base() + H*H*H*O*ep)

    def base_reaction(self):
        """Determines the  reaction in which the molecule acts as base"""
        return (self + H*H*O > self.conjugate_acid() + OH*en)

    def __mul__(self, x):
        """Overload the multiplication operator for 'ease of use'"""
        if type(x) == Element:
            # Molecule * Element -> Molecule
            if self.count == 1:
                result = self.copy()
                result.append(x)
                return result
            else:
                copy = self.copy()
                copy.extract_count()
                return Molecule({copy: self.count, x: 1})
        elif type(x) == Molecule:
            # Molecule * Molecule -> Molecule
            copy = self.copy()
            count = copy.extract_count()
            
            result = Molecule(dict([(i, j*count) for i, j in copy.parts.items()]), 1, self.charge * count)

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

    def __pow__(self, x):
        """Special case molecule construction, allows for greater control of structure"""
        if type(x) in [Molecule, Polyatomic, Element]:
            copy = self.copy()
            count = copy.extract_count()
            charge = copy.charge
            copy.charge = 0
            
            copy2 = x.copy()
            count2 = copy2.extract_count()
            charge2 = copy2.charge
            copy2.charge = 0
            
            return Molecule({copy: count, copy2: count2}, 1, charge + charge2)
        else:
            raise TypeError("Unsupported operation between {} and {}".format(type(self), type(x)))

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
        self.mol_mass = sum([j*i.molar_mass for i, j in self.parts.items()])
        self.ox = oxidation
        self.count = count
        self.state = state
        self.sym = self.symbol
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
        return self.symbol

    def __contains__(self, key):
        return key in self.parts

    def grams_to_moles(self, grams, element=None):
        """Converts grams of a molecule to moles.  If an element is specified, the number of
        moles of a particular element will be isolated and returned"""
        if element:
            return self.parts[element] * grams / self.molar_mass
        else:
            return grams / self.molar_mass

    def g2mol(self, grams):
        """Alias for grams_to_moles"""
        return self.grams_to_moles(grams)

    def moles_to_grams(self, moles, element=None):
        """Converts the number of moles of a molecule to grams.  If an element is specified,
        the mass of that element will be isolated and returned instead."""
        if element:
            return moles * (Molecule({element: self.parts[element]})).molar_mass
        else:
            return moles * self.molar_mass

    @property
    def molar_mass(self):
        """Gets the molar mass of the polyatomic ion"""
        self.mol_mass = sum([j*i.molar_mass for i, j in self.parts.items()])
        return self.mol_mass

    def oxidation(self):
        """Get the oxidation state of the thing"""
        return self.ox

    @property
    def symbol(self):
        """Returns the symbol of the element"""
        result = ""

        if self.count > 1:
            result += str(self.count)

        for i, j in self.parts.items():
            if type(i) == Element:
                result += i.symbol
                if j > 1:
                    result += str(j).translate(sub)
            else:
                if j > 1:
                    result += "(" + i.symbol + ")" + str(j).translate(sub)
                else:
                    result += i.symbol

        if self.state is not None:
            result += str(self.state)
            
        return result

    @property
    def enthalpy(self):
        """Uses a lookup table to find the polyatomic's molar enthalpy"""
        if str(self) in molar_enthalpy:
            return molar_enthalpy[str(self)] * u.J * 1000 / u.moles
        else:
            if self.state is None:
                raise AttributeError("The ΔH⁰_f for '{}' could not be found in the molar_enthalpy lookup table (Hint: make sure to specify state!)".formate(self))
            else:
                raise AttributeError("The ΔH⁰_f for '{}' could not be found in the molar_enthalpy lookup table".format(self))

    @property
    def entropy(self):
        """Uses a lookup table to find the polyatomic's standard entropy"""
        if str(self) in standard_entropy:
            return standard_entropy[str(self)] * u.J / u.moles / u.K
        else:
            if self.state is None:
                raise AttributeError("The S⁰ for '{}' could not be found in the standard_entropy lookup table (Hint: make sure to specify state!)".format(self))
            else:
                raise AttributeError("The S⁰ for '{}' could not be found in the standard_entropy lookup table".format(self))

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
        
        if x in self:
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

    def acid_reaction(self):
        """Determines the reaction in which the molecule acts as an acid"""
        return Molecule(self.parts, self.count, self.charge).acid_reaction()

    def base_reaction(self):
        """Determines the  reaction in which the molecule acts as base"""
        return Molecule(self.parts, self.count, self.charge).base_reaction()
    
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
        elif type(x) == Charge:
            # Polyatomic * charge
            copy = self.copy()
            count = copy.extract_count()
            return Molecule({copy: 1}, count) * x
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
        if type(x) in [Molecule, Polyatomic]:
            # Polyatomic + Molecule -> Combination
            count = self.extract_count()
            result = Combination({self: count})
            result.append(x)
            return result
        elif type(x) == Element:
            # Polyatomic + Element -> Combination
            count = self.extract_count()
            return Combination({self: count, Molecule({x: 1}): 1})
        elif type(x) == Combination:
            # Polyatomic + Combination -> Combination
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

        #print(result)
        return hash(tuple(result.items()))

    def __eq__(self, x):
        return hash(self) == hash(x)

    def __repr__(self):
        return self.symbol

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

    def __contains__(self, key):
        return key in self.parts

    @property
    def symbol(self):
        """Returns the symbolic representation of the combination"""
        result = ""
        for i, j in self.parts.items():
            if j != 1:
                result += str(j)
            result += i.symbol + " + "
        return result[:-3]

    @property
    def enthalpy(self):
        """Returns the total molar enthalpy of the combination"""
        result = 0
        for i, j in self.parts.items():
            result += j * i.enthalpy
        return result

    @property
    def entropy(self):
        """Returns the total standard entropy of the combination"""
        result = 0
        for i, j in self.parts.items():
            result += j * i.entropy
        return result

    def append(self, x):
        if type(x) not in [Element, Molecule, Polyatomic]:
            raise TypeError("Unable to append '{}' to 'Combination'".format(type(x)))
        
        count = x.extract_count()
        
        if x in self:
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
            result += i.molar_mass * j
        return result

    def dissolve(self, ignore_ph=True, ignore_check=False):
        """Dissociates the components of the combination"""
        result = Combination({})
        for i, j in self.parts.items():
            result += (j*i).dissolve(ignore_ph, ignore_check)
        return result

    def _assign_oxidation(self):
        """Assigns oxidation numbers to the ions in the combination"""
        stage1 = []
        # Set up component-charge dictionary, sorted by molar mass
        parts = list(dict(self.parts).items())
        parts.sort(key=lambda x: x[0].molar_mass)

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
        return self.left.symbol + " -> " + self.right.symbol

    def __mul__(self, x):
        if type(x) in [float, int]:
            result = self.copy()
            # Number * (Reaction)
            for j, k in self.left.parts.items():
                result.left.parts[j] *= x
            for j, k in self.right.parts.items():
                result.right.parts[j] *= x
            return result
        else:
            raise TypeError("Unsupported operation between {} and {}".format(type(self), type(x)))

    def __add__(self, x):
        if type(x) is Reaction:
            result = self.copy()
            copy = x.copy()
            result.left += copy.left
            result.right += copy.right
            return result
        else:
            raise TypeError("Unsupported operation between {} and {}".format(type(self), type(x)))

    def __rmul__(self, x):
        return self.__mul__(x)

    def copy(self):
        """Returns a new, separate copy of a reaction"""
        return Reaction(self.left.copy(), self.right.copy())

    def reverse(self):
        """Swaps both sides of the reaction"""
        return Reaction(self.right.copy(), self.left.copy())

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
                return self.solve(limit=limit, ignore_charge=True)
        
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
        if not self.verify():
            if input("WARNING: Reaction is not balanced, run auto-balance? (Y/N): ").lower() == "n":
                raise ValueError("Reaction '{}' is not balanced!", self)
            else:
                self = self.solve()
        reactors = {}
        print("Reacting:", self)
        print("Enter starting conditions (append 'grams' to denote a value in grams, otherwise moles are assumed):")
        for i, j in self.left.parts.items():
            new_reactor = input(i.symbol + ": ")
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

    def equil_const(self):
        """General equation solver for equilibrium coefficients"""
        variables = ["k"]
        products = "1"
        for i, j in self.right.parts.items():
            if i.state not in [l, s] and i != H*H*O:
                variables.append("[{}]".format(i))
                products += "*(__[{}])**{}".format(i, j)

        reactants = "1"
        for i, j in self.left.parts.items():
            if i.state not in [l, s] and i != H*H*O:
                variables.append("[{}]".format(i))
                reactants += "*(__[{}])**{}".format(i, j)
            
        general = "({}) / ({}) - __k".format(products, reactants)

        return Solver(variables, general).run()

    def init_to_equil(self):
        """Equation solver for equilibrium coefficients.  Applies for reactions that are some unit change, 'x,' away from equilibrium"""
    
        variables = ["k", "x_value"]
        products = "1"
        for i, j in self.right.parts.items():
            if i.state not in [l, s] and i != H*H*O:
                variables.append("[{}]".format(i))
                products += "*(__[{}] + {}*__x_value)**{}".format(i, j, j)

        reactants = "1"
        for i, j in self.left.parts.items():
            if i.state not in [l, s] and i != H*H*O:
                variables.append("[{}]".format(i))
                reactants += "*(__[{}] - {}*__x_value)**{}".format(i, j, j)
            
        general = "({}) / ({}) - __k".format(products, reactants)

        return Solver(variables, general).run()

    def reach_equil_detailed(self):
        """Calculates the concentrations/pressures of all species after a reaction that is not in equilibrium reaches equilibrium"""
        products = "1"
        reactants = "1"
        product_dict = {}
        product_counts = {}
        reactant_dict = {}
        reactant_counts = {}
        print(self)
        print("Enter initial concentrations/partial pressures")

        # Obtain products data
        for i, j in self.right.parts.items():
            if i != H*H*O:
                user_input = eval(input("[{}]: ".format(i)))
                product_dict[i] = user_input
                product_counts[i] = j
                if i.state not in [l, s]:
                    products += "*({} + {}*x)**{}".format(user_input, j, j)

        # Obtain reactants data
        for i, j in self.left.parts.items():
            if i != H*H*O:
                user_input = eval(input("[{}]: ".format(i)))
                reactant_dict[i] = user_input
                reactant_counts[i] = j
                if i.state not in [l, s]:
                    reactants += "*({} - {}*x)**{}".format(user_input, j, j)

        # Obtain equilibrium constant
        equil_constant = eval(input("k: "))

        # Solve the equation
        equation = "({}) / ({}) - {}".format(products, reactants, equil_constant)
        print(equation)
        solutions = solve_expr(parse_expr(equation))
        print("Solutions before checking:", solutions)

        from sympy import re, im, simplify
        
        # Determine which possible solutions are valid
        for i in solutions.copy():
            if type(i) is sympyAdd:
                solutions.remove(i)
                continue
            
            # Check if the resulting reactant concentrations are nonnegative
            for (k, m), (n, p) in zip(reactant_dict.items(), reactant_counts.items()):
                if m - p*i < 0:
                    solutions.remove(i)
                    break
            
            # Check if the resulting product concentrations are nonnegative
            for (k, m), (n, p) in zip(product_dict.items(), product_counts.items()):
                if m + p*i < 0:
                    solutions.remove(i)
                    break
            
        print("Valid solutions:", solutions)

        # Print the results of each valid solution
        for i, j in enumerate(solutions):
            print("Solution {}:".format(i))
            for (k, m), (n, p) in zip(reactant_dict.items(), reactant_counts.items()):
                print("{}: {}".format(k, m - p*j))
            for (k, m), (n, p) in zip(product_dict.items(), product_counts.items()):
                print("{}: {}".format(k, m + p*j))

        # There were no solutions
        if len(solutions) == 0:
            print("Unable to find a solution!")

    def enthalpy(self):
        """Returns the total net molar enthalpy change"""
        return self.right.enthalpy - self.left.enthalpy

    def entropy(self):
        """Returns the total molar standard entropy change"""
        return self.right.entropy - self.left.entropy

    def standard_free_energy(self, temperature, moles):
        """Returns the standard free energy of the reaction"""
        return (self.enthalpy() - self.entropy() * temperature * u.K) * moles * u.moles
                


class AcidBase:
    """Used to work with acid/base solutions"""

    def __init__(self, acid: Molecule, base: Molecule, ka: float=None, kb: float=None):
        # Check for presense of hydrogen in acid
        if H not in acid.get_atoms():
            raise ValueError("Acid should have at least one hydrogen!")
        # Check if hydrogen count in acid/base make sense
        elif H in base.get_atoms():
            if acid.get_atoms()[H] < base.get_atoms()[H]:
                raise ValueError("Acid has fewer hydrogen atoms than base?")
        self.acid = acid.copy()
        self.base = base.copy()
        self.__ka = ka
        self.__kb = kb

    def __repr__(self):
        return "<AcidBase: {} / {}>".format(self.acid, self.base)

    @property
    def ka(self):
        if self.__ka is None:
            if self.__kb is None:
                return None
            else:
                return (1e-14) / self.__kb
        else:
            return self.__ka

    @property
    def kb(self):
        if self.__kb is None:
            if self.__ka is None:
                return None
            else:
                return (1e-14) / self.__ka
        else:
            return self.__kb

    def acid_reach_equil(self):
        """Determine equilibrium concentrations given initial conditions.  Assumes acid is on the left."""
        self.acid_reaction().reach_equil_detailed()

    def base_reach_equil(self):
        """Determine equilibrium concentrations given initial conditions.  Assumes base is on the left."""
        self.base_reaction().reach_equil_detailed()

    def calc_pH(self):
        """Calculates the pH of the solution at equilibrium"""
        from sympy import re, im, simplify
        print(self)
        print("Enter initial concentrations:")

        # Obtain acid concentration data
        acid_conc = eval(input("[{}]: ".format(self.acid)))

        # Obtain base concentration data
        base_conc = eval(input("[{}]: ".format(self.base)))

        # Initialize H/H3O concentration
        x_conc = 0

        base_mode = None

        user_input = ""

        if acid_conc != 0 and base_conc != 0:
            user_input = input("Is the acid or base on the left? (acid/base): ").lower()

        # Guess the direction of the reaction given the concentrations
        if acid_conc == 0 or user_input == "base":
            # Assume base is on the left, since the alternative would yield no solution
            base_mode = True
            reaction = self.base_reaction()

            # Obtain concentration of OH-
            x_conc = eval(input("[OH⁻]: "))
            
            # Obtain equilibrium constant if none has been set
            if self.kb is None:
                self.set_kb(eval(input("k_b: ")))

            # Substitute values into the equation
            equation = "(({} + x) * ({} + x)) / (({} - x)) - {}".format(acid_conc, x_conc, base_conc, self.kb)
            print(equation)

        elif base_conc == 0 or user_input == "acid":
            # Assume acid is on the left, since the alternative would yield no solution
            base_mode = False
            reaction = self.acid_reaction()
            print("Assuming:", reaction)

            # Obtain concentration of H3O+
            x_conc = eval(input("[H₃O⁺]: "))

            # Obtain equilibrium constant if none has been set
            if self.ka is None:
                self.set_ka(eval(input("kₐ: ")))

            # Substitute values into the equation
            equation = "(({} + x) * ({} + x)) / (({} - x)) - {}".format(base_conc, x_conc, acid_conc, self.ka)
            print(equation)
        else:
            raise NotImplemented

        # Solve the equation
        solutions = solve_expr(parse_expr(equation))
        print("Solutions before checking:", solutions)

        # Remove invalid solutions
        for i in solutions.copy():
            # Remove imaginary solution
            if type(i) is sympyAdd:
                solutions.remove(i)
                continue

            if base_mode:
                # Check if resulting reactant concentrations are nonnegative
                if base_conc - i < 0:
                    solutions.remove(i)
                    continue
                # Check if resulting product concentrations are nonnegative
                if acid_conc + i < 0 or x_conc + i < 0:
                    solutions.remove(i)
                    continue
            else:
                # Check if resulting reactant concentrations are nonnegative
                if acid_conc - i < 0:
                    solutions.remove(i)
                    continue
                # Check if resulting product concentrations are nonnegative
                if base_conc + i < 0 or x_conc + i < 0:
                    solutions.remove(i)
                    continue
        
        print("Valid solutions:", solutions)

        # Print the results of each valid solution
        for i, j in enumerate(solutions):
            print("Solution {}:".format(i))
            # [X] + x gives the final [OH]
            if base_mode:
                print("pH:", 14 + log10(j + x_conc))
                
            # [X] + x gives the final [H]
            else:
                print("pH:", -log10(j + x_conc))

        if len(solutions) == 0:
            print("Unable to find a solution!")

    def equivalence_point(self, titrant=None, titrant_conc=None, conc=None, verbose=True):
        """Calculate the pH of the equivalence point in an acid/base titration"""
        # Obtain conditions
        if titrant is None:
            titrant = eval(input("Titrant Identity: "))
        if titrant in strong_acids:
            acid = True
        elif titrant in strong_bases:
            acid = False
        else:
            raise ValueError("Titrant was not a strong acid or base.  Make sure you are using polyatomics whenever possible!")
        if titrant_conc is None:
            titrant_conc = float(eval(input("[Titrant]: ")))
        if conc is None:
            if acid:
                conc = float(eval(input("[Base]₀: ")))
            else:
                conc = float(eval(input("[Acid]₀: ")))

        # Check for ka/kb
        if self.ka is None:
            raise ValueError("Please set either the ka or kb of the acid/base solution using [self].set_ka or [self].set_kb")

        # Compute the concentration of the conjugate solution
        conjugate_conc = conc / (conc/titrant_conc + 1)

        # Solve for the final concentration of H⁺/OH⁻
        if acid:
            # Titrant was a strong acid, use ka
            # Print information
            if verbose:
                print("Initial Reaction: ", self.add_H())
                print("Reaction at Equivalence Point: ", self.acid_reaction())
            # Solve for unit change
            solutions = solve_expr("x**2 / ({} - x) - {}".format(conjugate_conc, self.ka))
            if verbose: print("Potential solutions for [H⁺]: ", solutions)
        else:
            # Titrant was a strong base, use kb
            # Print information
            if verbose:
                print("Initial Reaction: ", self.add_OH())
                print("Reaction at Equivalence Point: ", self.base_reaction())
            # Solve for unit change
            solutions = solve_expr("x**2 / ({} - x) - {}".format(conjugate_conc, self.kb))
            if verbose: print("Potential solutions for [OH⁻]: ", solutions)

        # Eliminate invalid solutions
        for i, j in enumerate(solutions.copy()):
            # Solution should not be less than or equal to zero
            if j <= 0:
                del solutions[i]
            # Solution should not be imaginary
            elif type(j) is sympyAdd:
                del solutions[i]

        # First check if no valid solutions remain
        if len(solutions) == 0:
            raise Exception("An unknown error occured: all solutions were invalid")
        elif len(solutions) > 1:
            raise Exception("An unknown error occured: too many solutions were valid")

        # Print final output
        if verbose: 
            print("Valid Solution:", solutions[0])
            print("pH:")
        if acid:
            return -log10(solutions[0])
        else:
            return 14 + log10(solutions[0])
        
    def add_titrant(self, titrant=None, titrant_conc=None, titrant_volume=None, conc=None, volume=None):
        """Calculates the pH of an acid/base solution after adding a titrant.  General case."""
        # Obtain conditions
        # Titrant status as an acid or base
        if titrant is None:
            titrant = eval(input("Identity of Titrant: "))
        if titrant in strong_acids:
            acid = True
        elif titrant in strong_bases:
            acid = False
        else:
            raise ValueError("Titrant was not a strong acid or base.  Make sure you are using polyatomics whenever possible!")

        # Titrant concentration and volume added
        if titrant_conc is None:
            titrant_conc = float(eval(input("Titrant Concentration: ")))
        if titrant_volume is None:
            titrant_volume = float(eval(input("Volume of Titrant Added (L): ")))

        # In the event that there are two H⁺/OH⁻ ions
        if titrant in [H*H*SO4, Ca*(2*OH), Ba*(2*OH), Sr*(2*OH)]:
            titrant_conc *= 2

        # Acid/Base concentration and volume data
        if acid:
            if conc is None:
                conc = float(eval(input("[Base]₀: ")))
            if volume is None:
                volume = float(eval(input("Initial Volume of Base: ")))
        else:
            if conc is None:
                conc = float(eval(input("[Acid]₀: ")))
            if volume is None:
                volume = float(eval(input("Initial Volume of Acid: ")))

        if self.base in [Ca*(2*OH), Ba*(2*OH), Sr*(2*OH)]:
            conc *= 2
        elif self.acid in [H*H*SO4]:
            conc *= 2

        # Compute number of moles
        moles = conc * volume
        titrant_moles = titrant_conc * titrant_volume

        # Acid/Base happens to be a strong acid/base
        if self.acid in strong_acids or self.base in strong_bases:
            # Titrant does not completely react with acid/base
            if titrant_moles < moles:
                if acid:
                    return 14 + log10((moles - titrant_moles)/(volume + titrant_volume))
                else:
                    return -log10((moles - titrant_moles)/(volume + titrant_volume))
            # Titrant perfectly reacts with acid/base
            elif titrant_moles == moles:
                return 7
            # Titrant exceeds acid/base
            else:
                if acid:
                    return -log10((titrant_moles - moles)/(volume + titrant_volume))
                else:
                    return 14 + log10((titrant_moles-moles)/(volume + titrant_volume))

        # Acid/base is not a strong acid/base

        # Titrant does not completely react with acid/base
        if titrant_moles < moles:
            conjugate_moles = titrant_moles
            moles -= titrant_moles
            if acid:
                # Titrant is an acid, so 'moles' refers to moles of base
                return 14 + log10(self.kb) - log10(conjugate_moles / moles)
            else:
                # Titrant is an base, so 'moles' refers to moles of acid
                return -log10(self.ka) + log10(conjugate_moles / moles)

        # Titrant perfectly reacts with acid/base (equivalence point)
        elif titrant_moles == moles:
            return self.equivalence_point(titrant, titrant_conc, conc, False)

        # There were more moles of titrant than initial acid/base
        else:
            if acid:
                return -log10((titrant_moles - moles) / (volume + titrant_volume))
            else:
                return 14 + log10((titrant_moles - moles) / (volume + titrant_volume))
    

    def set_ka(self, value):
        """Manually set the Ka value of the solution"""
        self.__ka = value

    def set_kb(self, value):
        """Manually set the Kb value of the solution"""
        self.__kb = value

    def acid_reaction(self):
        """Returns the reaction where the acid is on the left"""
        return Reaction(self.acid.copy() + H*H*O, self.base.copy() + H*H*H*O*ep)

    def base_reaction(self):
        """Returns the reaction where the base is on the left"""
        return Reaction(self.base.copy() + H*H*O, self.acid.copy() + OH*en)

    def add_OH(self):
        """Returns the reaction that results from adding a strong base to the solution"""
        
        return Reaction(self.acid.copy() + OH*en,
                        self.acid.conjugate_base().copy() + H*H*O)

    def add_H(self):
        """Returns the reaction that results from adding a strong acid to the solution"""
        return Reaction(self.base.copy() + H*ep,
                        self.base.conjugate_acid().copy() + H*H*O)

    def add(self, molecule):
        """Returns the reaction that results from adding a strong acid/base to the solution"""
        # Molecule is a strong acid
        if molecule in strong_acids:
            return self.add_H()
        # Molecule is a strong base
        elif molecule in strong_bases:
            return self.add_OH()
        else:
            raise ValueError("Molecule '{}' is not a strong acid/base".format(molecule))


class Buffer(AcidBase):
    """Subclass of AcidBase for working with buffers"""

    def __init__(self, acid: Molecule, base: Molecule, ka: float=None, kb: float=None):
        super().__init__(acid, base)
        if acid in strong_acids:
            print("Warning: Solutions of strong acids and their conjugate base are generally not buffer solutions")
        elif base in strong_bases:
            print("Warning: Solutions of strong bases and their conjugate acid are generally not buffer solutions")

    def __repr__(self):
        return "<Buffer: {} / {}>".format(self.acid, self.base)

    def preparation_solver(self, approxilate=False):
        """General equation solver for preparing a buffer of a specific pH"""
        # Assume the acid is on the left
        # Note: Doesn't matter if acid/base is on the left, you get the same answer either way
        # variables: [BASE₀, ACID₀, pH, Ka]
        if self.ka is None:
            raise ValueError("Kₐ/K_b is not set.  Please define one or both using SELF.set_ka or/and SELF.set_kb")
        variables = ['[{}]₀'.format(self.base), '[{}]₀'.format(self.acid), '[H₃O⁺]₀', 'pH']
        equation = '(__[{}]₀ + 10**-__pH) * (__[H₃O⁺]₀ + 10**-__pH) / (__[{}]₀ - 10**-__pH) - {}'.format(self.base, self.acid, self.ka)
        Solver(variables, equation).run()

        if approxilate:
              approx = self.preparation_approxilator()
              print("'Incorrect' value:", approx)
              actual = eval(input("Paste true value here: "))
              print("Approxilation is {}% off".format(abs((actual - approx) / actual * 100)))

    def preparation_approxilator(self):
        """When [Homework System that Must Not be Named] marks you wrong because it approxilated while you didn't."""
        print("Approxilated Solver:")

        pH = eval(input("pH: "))
        
        base = eval(input("[{}]₀: ".format(self.base)))
        if base == 0:
            part1 = str(10**-pH)
        else:
            part1 = str(base)
        
        acid = eval(input("[{}]₀: ".format(self.acid)))
        if acid == 0:
            part2 = str(10**-pH)
        else:
            part2 = str(acid)
        
        H3O = eval(input("[H₃O⁺]₀: "))
        if H3O == 0:
            part3 = str(10**-pH)
        else:
            part3 = str(H3O)

        approxilated = part1 + ' * ' + part3 + ' / (' + part2 + ')' + ' - ' + str(self.ka)
        print(approxilated)
        return solve_expr(approxilated)[0]


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

# Declare data scraped from a table
molar_enthalpy = {'Al(s)': 0, 'AlCl₃(s)': -704.2, 'Al₂O₃(s)': -1675.7, 'Al(OH)₃(s)': -1277.0, 'Ba(s)': 0.0, 'BaCl₂(s)': -858.6, 'BaCO₃(s)': -1219.0, 'BaO(s)': -553.5, 'Ba(OH)₂(s)': -946.0, 'BaSO₄(s)': -1473.2, 'Be(s)': 0.0, 'BeO(s)': -599.0, 'Be(OH)₂(s)': -902.5, 'Br(g)': 111.9, 'Br₂(l)': 0.0, 'Br₂(g)': 30.9, 'Br₂(aq)': -3.0, '[Br]⁻(aq)': -121.0, 'BrF₃(g)': -255.6, 'HBr(g)': -36.3, 'Cd(s)': 0.0, 'CdO(s)': -258.0, 'Cd(OH)₂(s)': -561.0, 'CdS(s)': -162.0, 'CdSO₄(s)': -935.0, 'Ca(s)': 0.0, 'Ca(g)': 178.2, '[Ca]²⁺(g)': 1925.9, 'CaC₂(s)': -59.8, 'CaCO₃(s, calcite)': -1206.9, 'CaCl₂(s)': -795.8, 'CaF₂(s)': -1219.6, 'CaH₂(s)': -186.2, 'CaO(s)': -635.1, 'CaS(s)': -482.4, 'Ca(OH)₂(s)': -986.1, 'Ca(OH)₂(aq)': -1002.8, 'Ca₃(PO₄)₂(s)': -4126.0, 'CaSO₄(s)': -1434.1, 'CaSiO₃(s)': -1630.0, 'C (s, graphite)': 0.0, 'C (s, diamond)': 1.9, 'C(g)': 716.7, 'CCl₄(l)': -135.4, 'CCl₄(g)': -102.9, 'CHCl₃(l)': -134.5, 'CHCl₃(g)': -103.1, 'CH₄(g)': -74.8, 'CH₃OH(g)': -200.7, 'CH₃OH(l)': -238.7, 'H₂CO(g)': -116.0, 'HCOOH(g)': -363.0, 'HCN(g)': 135.1, 'C₂H₂(g)': 226.7, 'C₂H₄(g)': 52.3, 'CH₃CHO (g, acetaldehyde)': -166.0, 'C₂H₄O (g, ethylene oxide)': -53.0, 'CH₃CH₂OH(l)': -277.7, 'CH₃CH₂OH(g)': -235.1, 'CH₃COOH(l)': -484.0, 'C₂H₆(g)': -84.7, 'C₃H₆(g)': 20.9, 'C₃H₈(g)': -103.8, 'CH₂=CHCN(l)': 152.0, 'C₆H₆(l)': 49.0, 'C₆H₁₂O₆(s)': -1275.0, 'CO(g)': -110.5, 'CO₂(g)': -393.5, 'CS₂(g)': 117.4, 'COCl₂(g)': -218.8, 'Cl(g)': 121.7, 'Cl₂(g)': 0.0, 'Cl₂(aq)': -23.0, '[Cl]⁻(aq)': -167.0, '[Cl]⁻(g)': -233.1, 'HCl(g)': -92.3, 'HCl(aq)': -167.2, 'Cr(s)': 0.0, 'Cr₂O₃(s)': -1139.7, 'CrO₃(s)': -579.0, 'CrCl₃(s)': -556.5, 'Cu(s)': 0.0, 'CuCl₂(s)': -220.1, 'CuCO₃(s)': -595.0, 'Cu₂O(s)': -170.0, 'CuO(s)': -157.3, 'Cu(OH)₂(s)': -450.0, 'CuS(s)': -49.0, 'F₂(g)': 0.0, 'F(g)': 79.0, '[F]⁻(g)': -255.4, '[F]⁻(aq)': -332.6, 'HF(g)': -271.1, 'HF(aq)': -332.6, 'H₂(g)': 0.0, 'H(g)': 218.0, '[H]⁺(g)': 1536.2, '[H]⁺(aq)': 0.0, '[OH]⁻(aq)': -230.0, 'H₂O(l)': -285.8, 'H₂O(g)': -241.8, 'H₂O₂(l)': -187.8, 'I₂(s)': 0.0, 'I₂(g)': 62.4, 'I₂(aq)': 23.0, 'I(g)': 106.8, '[I]⁻(g)': -197.0, '[I]⁻(aq)': -55.0, 'ICl(g)': 17.8, 'Fe(s)': 0.0, 'Fe₃C(s)': 21.0, 'FeCl₂(s)': -341.8, 'FeCl₃(s)': -399.5, 'Fe₀.₉₅O(s) (wustite)': -264.0, 'FeO(s)': -272.0, 'Fe₃O₄(s)': -1118.4, 'Fe₂O₃(s)': -824.2, 'FeS(s)': -95.0, 'FeS₂(s)': -178.2, 'FeSO₄(s)': -929.0, 'Fe(CO)₅(l)': -774.0, 'Pb(s)': 0.0, 'PbCl₂(s)': -359.4, 'PbO (s)': -217.3, 'PbO₂(s)': -277.0, 'PbS(s)': -100.4, 'PbSO₄(s)': -920.0, 'Mg(s)': 0.0, 'MgCl₂(s)': -641.3, 'MgCO₃(s)': -1095.8, 'MgO(s)': -601.7, 'Mg(OH)₂(s)': -924.5, 'MgS(s)': -346.0, 'Mn(s)': 0.0, 'MnO(s)': -385.0, 'Mn₃O₄(s)': -1387.0, 'Mn₂O₃(s)': -971.0, 'MnO₂(s)': -521.0, '[MnO₄]⁻(aq)': -543.0, 'Hg(l)': 0.0, 'HgCl₂(s)': -224.3, 'Hg₂Cl₂(s)': -265.4, 'HgO (s)': -90.8, 'HgS (s)': -58.2, 'Ni(s)': 0.0, 'NiCl₂(s)': -305.3, 'NiO(s)': -239.7, 'Ni(OH)₂(s)': -538.0, 'NiS(s)': -93.0, 'N₂(g)': 0.0, 'N(g)': 472.7, 'NH₃(g)': -46.1, 'NH₃(aq)': -80.0, '[NH₄]⁺(aq)': -132.0, 'NO(g)': 90.3, 'NOCl(g)': 51.7, 'NO₂(g)': 33.2, 'N₂O(g)': 82.1, 'N₂O₄(g)': 9.2, 'N₂O₄(l)': -20.0, 'N₂O₅(s)': -42.0, 'N₂H₄(l)': 50.6, 'N₂H₃CH₃(l)': 54.0, 'HNO₃(aq)': -207.4, 'HNO₃(l)': -174.1, 'HNO₃(g)': -135.1, 'NH₄ClO₄(s)': -295.0, 'NH₄Cl(s)': -314.4, 'NH₄Cl(aq)': -299.7, 'NH₄NO₃(s)': -365.6, 'NH₄NO₃(aq)': -339.9, 'O₂(g)': 0.0, 'O(g)': 249.2, 'O₃(g)': 142.7, 'P (s, white)': 0.0, 'P (s, red)': -70.4, 'P (s, black)': -39.0, 'P(g)': 314.6, 'P₄(s, white)': 0.0, 'P₄(s, red)': -17.6, 'P₄(g)': 59.0, 'PF₅(g)': -1578.0, 'PH₃(g)': 5.4, 'PCl₃(g)': -287.0, 'H₃PO₄(l)': -1279.0, 'H₃PO₄(aq)': -1288.0, 'P₄O₁₀(s)': -2984.0, 'K(s)': 0.0, 'KCl(s)': -436.7, 'KClO₃(s)': -397.7, 'KClO₄(s)': -433.0, 'KI(s)': -327.9, 'K₂O(s)': -361.0, 'K₂O₂(s)': -496.0, 'KO₂(s)': -283.0, 'KOH(s)': -424.8, 'KOH(aq)': -482.4, 'Si(s)': 0.0, 'SiBr₄(l)': -457.3, 'SiC(s)': -65.3, 'SiCl₄(g)': -657.0, 'SiH₄(g)': 34.3, 'SiF₄(g)': -1614.9, 'SiO₂(s)': -910.9, 'Ag(s)': 0.0, '[Ag]⁺(aq)': 105.0, 'AgBr(s)': -100.0, 'AgCN(s)': 146.0, 'AgCl(s)': -127.1, 'Ag₂CrO₄(s)': -712.0, 'AgI(s)': -62.0, 'Ag₂O(s)': -31.1, 'AgNO₃(s)': -124.4, 'Ag₂S(s)': -32.0, 'Na(s)': 0.0, 'Na(g)': 107.3, '[Na]⁺(g)': 609.4, '[Na]⁺(aq)': -240.0, 'NaBr(s)': -361.0, 'Na₂CO₃(s)': -1130.7, 'NaHCO₃(s)': -948.0, 'NaCl(s)': -411.2, 'NaCl(g)': -176.7, 'NaCl(aq)': -407.3, 'NaH(s)': -56.0, 'NaI(s)': -288.0, 'NaNO₂(s)': -359.0, 'NaNO₃(s)': -467.0, 'Na₂O(s)': -416.0, 'Na₂O₂(s)': -515.0, 'NaOH(s)': -425.6, 'NaOH(aq)': -470.1, 'S (s, rhombic)': 0.0, 'S (s, monoclinic)': 0.3, 'S(g)': 278.8, '[S₂]⁻(aq)': 33.0, 'S₈(g)': 102.0, 'S₂Cl₂(g)': -18.4, 'SF₆(g)': -1209.0, 'H₂S(g)': -20.6, 'SO₂(g)': -296.8, 'SO₃(g)': -395.7, 'SOCl₂(g)': -212.5, '[SO₄]²⁻(aq)': -909.0, 'H₂SO₄(l)': -814.0, 'H₂SO₄(aq)': -909.3, 'Sn (s, white)': 0.0, 'Sn (s, gray)': -2.1, 'SnCl₄(l)': -511.3, 'SnCl₄(g)': -471.5, 'SnO(s)': -285.0, 'SnO₂(s)': -580.7, 'Sn(OH)₂(s)': -561.0, 'Ti(s)': 0.0, 'TiCl₄(l)': -804.2, 'TiCl₄(g)': -763.2, 'TiO₂(s)': -939.7, 'U(s)': 0.0, 'UF₆(s)': -2137.0, 'UF₆(g)': -2113.0, 'UO₂(s)': -1084.0, 'U₃O₈(s)': -3575.0, 'UO₃(s)': -1230.0}
standard_entropy = {'Al(s)': 28.3, 'AlCl₃(s)': 110.7, 'Al₂O₃(s)': 50.9, 'Ba(s)': 67.0, 'BaCl₂(s)': 123.7, 'BaCO₃(s)': 112.0, 'BaO(s)': 70.4, 'BaSO₄(s)': 132.2, 'Be(s)': 9.5, 'BeO(s)': 14.0, 'Be(OH)₂(s)': 51.9, 'Br(g)': 175.0, 'Br₂(l)': 152.2, 'Br₂(g)': 245.5, 'Br₂(aq)': 130.0, '[Br]⁻(aq)': 82.0, 'BrF₃(g)': 292.5, 'HBr(g)': 198.7, 'Cd(s)': 52.0, 'CdO(s)': 55.0, 'Cd(OH)₂(s)': 96.0, 'CdS(s)': 65.0, 'CdSO₄(s)': 123.0, 'Ca(s)': 41.4, 'Ca(g)': 158.9, 'CaC₂(s)': 70.0, 'CaCO₃(s)': 92.9, 'CaCl₂(s)': 104.6, 'CaF₂(s)': 68.9, 'CaH₂(s)': 42.0, 'CaO(s)': 39.8, 'CaS(s)': 56.5, 'Ca(OH)₂(s)': 83.4, 'Ca(OH)₂(aq)': -74.5, 'Ca₃(PO₄)₂(s)': 241.0, 'CaSO₄(s)': 106.7, 'CaSiO₃(s)': 84.0, 'C (s, graphite)': 5.7, 'C (s, diamond)': 2.4, 'C(g)': 158.1, 'CCl₄(l)': 216.4, 'CCl₄(g)': 309.9, 'CHCl₃(l)': 201.7, 'CHCl₃(g)': 295.7, 'CH₄(g)': 186.3, 'CH₃OH(g)': 239.8, 'CH₃OH(l)': 126.8, 'H₂CO(g)': 219.0, 'HCOOH(g)': 249.0, 'HCN(g)': 202.0, 'C₂H₂(g)': 200.9, 'C₂H₄(g)': 219.6, 'CH₃CHO (g, acetaldehyde)': 250.0, 'C₂H₄O (g, ethylene oxide)': 242.0, 'CH₃CH₂OH(l)': 160.7, 'CH₃CH₂OH(g)': 282.7, 'CH₃COOH(l)': 160.0, 'C₂H₆(g)': 229.6, 'C₃H₆(g)': 266.9, 'C₃H₈(g)': 269.9, 'CH₂=CHCN(l)': 274.0, 'C₆H₆(l)': 172.8, 'C₆H₁₂O₆(s)': 212.0, 'CO(g)': 197.7, 'CO₂(g)': 213.7, 'CS₂(g)': 237.8, 'COCl₂(g)': 283.5, 'Cl(g)': 165.2, 'Cl₂(g)': 223.1, 'Cl₂(aq)': 121.0, '[Cl]⁻(aq)': 57.0, 'HCl(g)': 186.9, 'HCl(aq)': 56.5, 'Cr(s)': 23.8, 'Cr₂O₃(s)': 81.2, 'CrO₃(s)': 72.0, 'CrCl₃(s)': 123.0, 'Cu(s)': 33.2, 'CuCl₂(s)': 108.1, 'CuCO₃(s)': 88.0, 'Cu₂O(s)': 93.0, 'CuO(s)': 42.6, 'Cu(OH)₂(s)': 108.0, 'CuS(s)': 67.0, 'F₂(g)': 202.8, 'F(g)': 158.8, '[F]⁻(aq)': -13.8, 'HF(g)': 173.8, 'HF(aq)': 88.7, 'H₂(g)': 130.7, 'H(g)': 114.7, '[OH]⁻(aq)': -11.0, 'H₂O(l)': 69.9, 'H₂O(g)': 188.8, 'H₂O₂(l)': 109.6, 'I₂(s)': 116.1, 'I₂(g)': 260.7, 'I₂(aq)': 137.0, 'I(g)': 180.8, '[I]⁻(aq)': 106.0, 'ICl(g)': 247.6, 'Fe(s)': 27.8, 'Fe₃C(s)': 108.0, 'FeCl₂(s)': 118.0, 'FeCl₃(s)': 142.3, 'Fe₀.₉₅O(s) (wustite)': 59.0, 'Fe₃O₄(s)': 146.4, 'Fe₂O₃(s)': 87.4, 'FeS(s)': 67.0, 'FeS₂(s)': 52.9, 'FeSO₄(s)': 121.0, 'Fe(CO)₅(l)': 338.1, 'Pb(s)': 64.8, 'PbCl₂(s)': 136.0, 'PbO (s)': 68.7, 'PbO₂(s)': 69.0, 'PbS(s)': 91.2, 'PbSO₄(s)': 149.0, 'Mg(s)': 32.7, 'MgCl₂(s)': 89.6, 'MgCO₃(s)': 65.7, 'MgO(s)': 26.9, 'Mg(OH)₂(s)': 63.2, 'MgS(s)': 50.3, 'Mn(s)': 32.0, 'MnO(s)': 60.0, 'Mn₃O₄(s)': 149.0, 'Mn₂O₃(s)': 110.0, 'MnO₂(s)': 53.0, '[MnO₄]⁻(aq)': 190.0, 'Hg(l)': 75.9, 'HgCl₂(s)': 146.0, 'HgO (s)': 70.3, 'HgS (s)': 82.4, 'Ni(s)': 29.9, 'NiCl₂(s)': 97.7, 'NiO(s)': 38.0, 'Ni(OH)₂(s)': 79.0, 'NiS(s)': 53.0, 'N₂(g)': 191.6, 'N(g)': 153.3, 'NH₃(g)': 192.5, 'NH₃(aq)': 111.0, '[NH₄]⁺(aq)': 113.0, 'NO(g)': 210.8, 'NOCl(g)': 261.8, 'NO₂(g)': 240.1, 'N₂O(g)': 219.9, 'N₂O₄(g)': 304.3, 'N₂O₄(l)': 209.0, 'N₂O₅(s)': 178.0, 'N₂H₄(l)': 121.2, 'N₂H₃CH₃(l)': 166.0, 'HNO₃(aq)': 146.4, 'HNO₃(l)': 155.6, 'HNO₃(g)': 266.4, 'NH₄ClO₄(s)': 186.0, 'NH₄Cl(s)': 94.6, 'NH₄Cl(aq)': 169.9, 'NH₄NO₃(s)': 151.1, 'NH₄NO₃(aq)': 259.8, 'O₂(g)': 205.1, 'O(g)': 161.1, 'O₃(g)': 238.9, 'P (s, white)': 164.4, 'P (s, red)': 91.2, 'P (s, black)': 23.0, 'P(g)': 163.2, 'P₄(s, white)': 41.1, 'P₄(s, red)': 22.8, 'P₄(g)': 280.0, 'PF₅(g)': 296.0, 'PH₃(g)': 210.2, 'PCl₃(g)': 311.8, 'H₃PO₄(l)': 110.5, 'H₃PO₄(aq)': 158.0, 'P₄O₁₀(s)': 228.9, 'K(s)': 64.2, 'KCl(s)': 82.6, 'KClO₃(s)': 143.1, 'KClO₄(s)': 151.0, 'KI(s)': 106.3, 'K₂O(s)': 98.0, 'K₂O₂(s)': 113.0, 'KO₂(s)': 117.0, 'KOH(s)': 78.9, 'KOH(aq)': 91.6, 'Si(s)': 18.3, 'SiBr₄(l)': 277.8, 'SiC(s)': 16.6, 'SiCl₄(g)': 330.7, 'SiH₄(g)': 204.6, 'SiF₄(g)': 282.5, 'SiO₂(s)': 41.8, 'Ag(s)': 42.6, '[Ag]⁺(aq)': 73.0, 'AgBr(s)': 107.0, 'AgCN(s)': 84.0, 'AgCl(s)': 96.2, 'Ag₂CrO₄(s)': 217.0, 'AgI(s)': 115.0, 'Ag₂O(s)': 121.3, 'AgNO₃(s)': 140.9, 'Ag₂S(s)': 146.0, 'Na(s)': 51.2, 'Na(g)': 153.7, '[Na]⁺(aq)': 59.0, 'NaBr(s)': 86.8, 'Na₂CO₃(s)': 135.0, 'NaHCO₃(s)': 102.0, 'NaCl(s)': 72.1, 'NaCl(g)': 229.8, 'NaCl(aq)': 115.5, 'NaH(s)': 40.0, 'NaI(s)': 91.0, 'NaNO₃(s)': 116.0, 'Na₂O(s)': 73.0, 'Na₂O₂(s)': 95.0, 'NaOH(s)': 64.5, 'NaOH(aq)': 48.1, 'S (s, rhombic)': 31.8, 'S (s, monoclinic)': 33.0, 'S(g)': 167.8, '[S₂]⁻(aq)': -15.0, 'S₈(g)': 431.0, 'S₂Cl₂(g)': 331.5, 'SF₆(g)': 291.8, 'H₂S(g)': 205.8, 'SO₂(g)': 248.2, 'SO₃(g)': 256.8, 'SOCl₂(g)': 309.8, '[SO₄]²⁻(aq)': 20.0, 'H₂SO₄(l)': 156.9, 'H₂SO₄(aq)': 20.1, 'Sn (s, white)': 51.6, 'Sn (s, gray)': 44.1, 'SnCl₄(l)': 258.6, 'SnCl₄(g)': 365.8, 'SnO(s)': 56.0, 'SnO₂(s)': 52.3, 'Sn(OH)₂(s)': 155.0, 'Ti(s)': 30.6, 'TiCl₄(l)': 252.3, 'TiCl₄(g)': 354.9, 'TiO₂(s)': 49.9, 'U(s)': 50.0, 'UF₆(s)': 228.0, 'UF₆(g)': 380.0, 'UO₂(s)': 78.0, 'U₃O₈(s)': 282.0, 'UO₃(s)': 99.0}


# Declare some constants
R = 0.082058 * u.l * u.atm / u.mol / u.K # Ideal gas constant in L atm / (mol K)
R2 = 8.314 * u.J / u.K / u.mol

# Declare some basic solvers

pvnrt_solver = Solver(['P', 'V', 'n', 'T'],
                      '__P * __V - __n * 0.082058*atmosphere*liter/(kelvin*mole) * __T',
                      units=[u.atm, u.l, u.moles, u.K],
                      hints=['atm', 'L', 'moles','K'])

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

gibbs_solver = Solver(['ΔG⁰', 'ΔH⁰', 'S⁰', 'T'],
                      '__ΔH⁰ - __S⁰ * __T - __ΔG⁰')

# Declare global 'x' variable
x=symbols('x', real=True)

print("Use the function 'docs()' to print out all docstrings associated with all classes and their member functions")
print("An optional argument may be entered to narrow the results to a single class")
if DEBUG: print("[Debug]: {} Polyatomic ions loaded".format(len(Polyatomic.polyatomic_list)))
do_debug = True

# EXPERIMENTAL FUNCTIONS
def clear(n=35):
    """Clears the console"""
    print("\n"*n)

def dissolve_in_ph(salt, k_sp, ka):
    """[NOT GUARANTEED TO BE STABLE, OR WORK, AT ALL]
Calculates the new equilibrium constant after dissolving a salt in an acid solution
and provides the underlying reaction"""
    dissolution = salt > salt.dissolve(True, True)
    base, n = list(salt.dissolve(True, True).parts.items())[0]
    base_reaction = n * AcidBase(base.conjugate_acid(), base).acid_reaction().reverse()
    net_reaction = (dissolution + base_reaction).net_ionic()
    print(net_reaction)
    if k_sp is not None and ka is not None:
            print("k = kₛₚ / kₐ{}:".format(str(n).translate(sup)))
            print("k = {:e}".format(k_sp / ka**n))


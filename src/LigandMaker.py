import numpy as np
import scipy 

import pybel
import rdkit
from rdkit import Chem
import MDAnalysis as mda

"""
Last Updated: 17/12/2020 

Author: Sang Young Noh 

TODO 

"""
class Converter():
    """
    Written Manual here 
    """
    def __init__(self, option):
        """
        What does this constructor do? 
        """
        if option not in ['mol', 'sdf']:
            except ValueError:
                print("Need a mol or sdf input") 
        else:
            self.option = option
            

    def _UniverseConverterFromMol():
        """
        What does this function do? 
        """
        if self.option == 'mol':
            # If the chosen file is a mol, then we need to convert it to smiles
            moleculeSmile = Chem.MolToSmiles( 
            
    

from src.ENSDF import ENSDF
from src.ENSDF_Record import Time_scale
import numpy
from os.path import isfile
from os import access, R_OK
import warnings
import csv
import sys



#class DecaySchemes contains all transitions from the specified ENSDF files
class DecaySchemes:
    def __init__(self, filenames):
        # initialization of parameters
        self.filenames = filenames
        self.decays = []
        self.transition_all = []
        self.parse(self.filenames)
        self.transition_setup()
        self.transition_fill0()

        # Get the starting points and end points of decay
        self.StartPoints = self.Search_StartPoint()
        self.FinalPoints = self.Search_FinalPoint()
        self.transition_all.sort(key = lambda x:float(x.get('energy')))
    # read filenames
    def parse(self, filenames):
            try:
                for filename in filenames:
                    # check if files exist and are readable
                    if isfile(filename) and access(filename, R_OK):
                        self.decays.append(Decay(filename))
                    elif filename:
                        raise Exception(filename)
            except Exception as inst:
                print("Error: " + inst.args[0]+ " does not exist or is not readable! ")
                sys.exit(1)



    # gather all transitions from decays
    def transition_setup(self):
        for decay in self.decays:
            self.transition_all.extend(decay.Transitions)

        # change energis into strings with 2 decimal digits
        for transition in self.transition_all:
            transition['energy'] = '{:.2f}'.format(transition['energy'])

    # fill the parameters which has no value with 0
    def transition_fill0(self):
        for transition in self.transition_all:
            if transition['Hf'] == None and transition['gamma'] == 1:
                transition['Hf'] = '0'
            elif transition['Hf'] == None and transition['gamma'] == 0:
                transition['Hf'] = '0'

            if transition['alpha'] == None:
                transition['alpha'] = 0

    # search the starting nuclidic level of the decays
    def Search_StartPoint(self):
        starts = [transition['InitialState'] for transition in self.transition_all if\
            transition['InitialState']   not in [transition['FinalState'] for \
            transition in self.transition_all]]
        return set([x for x in starts if x.split("_")[1] not\
             in [y.DaughterNuclide for y in self.decays]])

    # search the final nuclidic level of the decays
    def Search_FinalPoint(self):
        return set([transition['FinalState'] for transition in self.transition_all if\
            transition['FinalState']   not in [transition['InitialState'] for \
            transition in self.transition_all]])

    # write transition list to a csv file
    def OutPut(self, Ofile):
        with open(Ofile, 'w', encoding = 'utf8', newline = '') as OutCSV:
            fc = csv.DictWriter(OutCSV, fieldnames=self.transition_all[0].keys())
            fc.writeheader()
            fc.writerows(self.transition_all)


class Decay:
    def __init__(self, filename):
        # initialization of parameters
        self.filename = filename
        self.ParentNuclide = None
        self.DaughterNuclide = None
        self.HalfLife = None
        self.ParentEnergy = None
        self.Norm = None
        self.Level_energies = []
        self.ENSDF_lines = []
        self.Transitions = []

        self.open(self.filename)
        self.format()

        #parse ENSDF lines
        try:
            self.ENSDF = ENSDF(self.ENSDF_lines)
        except Exception as info:
            print("Error: Parsing failed at {}. \nReason: {}".format(self.filename, info.args[0]))

        #get necessary info from the ENSDF:
        self.DecayInfoSetup()

        #set up excited levels:
        self.LevelSetup()


        #set up transitions:
        self.TransSetup()


        #sort transitions according to their energies:
        self.sort_transition()

    def float(self, str):
        return None if str=='' else numpy.float64(str)

    def float0(self, str):
        return 0 if str=='' else numpy.float64(str)

    # find the index of the level which has specific energy
    def find_level(self, energy):
        l = [abs(x-energy) for x in self.Level_energies]
        return l.index(min(l))

    # sort the transition list according to their energies
    def sort_transition(self):
        self.Transitions.sort(key = lambda x:x.get('energy'))

    def open(self, filename):
        try:
            with open(filename, 'r') as inputfile:
                lines = inputfile.read().splitlines()
        except FileNotFoundError:
            print("Error: File \"{}\" is not found".format(filename))
        else:
            self.ENSDF_lines = lines

    def format(self):
        #get rid of the line that contains space or is empty
        self.ENSDF_lines = [line for line in self.ENSDF_lines if not(line.isspace()) and bool(line)]

        #set string length to 80
        self.ENSDF_lines = ["{:<80}".format(line) for line in self.ENSDF_lines]

    def DecayInfoSetup(self):
        #set up parent nuclide:
        self.ParentNuclide = self.ENSDF.Parent['NUCID'].strip()

        #set up daughter nuclide:
        self.DaughterNuclide = self.ENSDF.ID['NUCID'].strip()

        #set up half-life of parent nuclide:
        time_list = self.ENSDF.Parent['T'].strip().split()
        try:
            self.HalfLife = Time_scale[time_list[1]](self.float(time_list[0]))
        except Exception:
            print("Error occurs during time transformation!\n Time string is: ", self.ENSDF.Parent['T'],\
                "set half-life to 100 seconds by default")
            self.HalfLife = 100


        #set up initial energy of parent nuclide
        try:
            self.ParentEnergy =  self.float(self.ENSDF.Parent['E'].strip())
        except Exception:
            print("Error occurs during energy transfmation of parent nuclide initial energy: ",\
            self.ENSDF.Parent['E'])

        #set up normalization factor for relative intensity:
        try:
            self.Norm = self.float(self.ENSDF.Norm['NR'].strip())
        except Exception:
            print("Error occurs during nomalization setting: ",\
            self.ENSDF.Norm['NR'])

    # set up level information from ENSDF object
    def LevelSetup(self):
        for level in self.ENSDF.levels:
            try:
                level_energy = self.float(level.par['E'].strip())
            except Exception as err:
                print("Error occurs during level setup: {}".format(err))
                print(level.par)
            else:
                self.Level_energies.append(level_energy)

    # set up the transition information from ENSDF object
    def TransSetup(self):
        for level in self.ENSDF.levels:

            # get the excitation energy of initial level
            try:
                level_energy = self.float(level.par['E'].strip())
            except Exception as err:
                print("Error occurs at level_energy conversion during transition setup:\n {}"\
                    .format(err))
                print(level.par)

            if level_energy is None:
                break

            # get index of initial state
            init_index = self.find_level(level_energy)

            # get half life of the excited state
            try:
                level_hf = level.par['T'].strip()
                if not level_hf:
                    level_hf = 0
                elif level_hf == 'STABLE':
                    level_hf = 1e5
                else:
                    time_list = level_hf.split()
                    level_hf = Time_scale[time_list[1]](self.float0(time_list[0]))
            except Exception as err:
                print("Error occurs at level hf conversion during transition setup:\n {}"\
                    .format(err))
                print(level.par)

            # store the necessary information for each transition
            for transition in level.transitions:
                try:
                    # get the value of transition energy
                    tran_energy = self.float0(transition['E'].strip())

                    # get the final state of the transition
                    lower_level_energy = level_energy - tran_energy
                    final_index = self.find_level(lower_level_energy)

                    # start to add information to the list
                    if transition['Trans_type'] == 'Gamma':
                        # add information to the list if it's gamma transition
                        self.Transitions.append({
                        'energy': tran_energy,
                        'intensity': round(self.float0(transition['RI'].strip()) * self.Norm, 4),
                        'InitialState': '{:.3f}_{}'.format(self.Level_energies[init_index], self.DaughterNuclide),
                        'FinalState': '{:.3f}_{}'.format(self.Level_energies[final_index], self.DaughterNuclide),
                        'alpha': self.float0(transition['CC'].strip()),
                        'Hf':  "{:#.4g}".format(level_hf),
                        'gamma': 1
                    })
                    else:
                        # add information to the list if it's non-gamma transition
                        self.Transitions.append({
                        'energy': tran_energy,
                        'intensity': self.float0(transition['I'].strip()),
                        'InitialState': '{:.3f}_{}'.format(self.ParentEnergy, self.ParentNuclide),
                        'FinalState': '{:.3f}_{}'.format(self.Level_energies[init_index], self.DaughterNuclide),
                        'alpha': 0,
                        'Hf': self.HalfLife,
                        'gamma': 0
                    })
                except Exception as err:
                    print("Error occurs at transition parsing during transition setup:\n {}"\
                    .format(err))
                    print(transition)

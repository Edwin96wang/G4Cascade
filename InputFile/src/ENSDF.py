from src.ENSDF_Record import ENSDF_Record as Record, ID_positions as id_pos
import re

class ENSDF:
    def __init__(self, lines) -> None:
        self.lines = []
        self.levels =[]
        self.ID = {}
        self.Parent = {}
        self.Norm = {}

        # specify the transitions that are needed
        self.transition_type = ["Beta", "EC", "Alpha", "Gamma"]
        
        #select only lines that have ENSD_Record
        self.select(lines)

        # parsing each with the corresponding ENSD_record -> self.levels is filled
        self.scan()

    def select(self, lines):
        for line in lines:
            flag = line[id_pos]
            self.lines.extend([(k, line) for k, v in Record.items() \
                 if bool(re.match(v["ID_pa"], flag))])

    def scan(self):
        for line in self.lines:

            #create a new level
            if line[0] == 'Level':
                self.levels.append(level_ENSDF(line[1]))

            # get the information of daughter nuclide
            elif line[0] == 'ID':
                self.scan_ID(line)

            # get the information of parent nuclide
            elif line[0] == 'Parent':
                self.scan_Parent(line)

            # get the values of normalization factors
            elif line[0] == 'Norm':
                self.scan_Norm(line)

            # add transitions to the level just added
            if self.levels and line[0] in self.transition_type:
                self.levels[-1].Add_Transition(line)

    def scan_ID(self, line):
        self.ID = dict((k, line[1][v]) for k, v in\
                 Record[line[0]].items() if k != "ID_pa")

    def scan_Parent(self, line):
        self.Parent = dict((k, line[1][v]) for k, v in\
                 Record[line[0]].items() if k != "ID_pa")

    def scan_Norm(self, line):
        self.Norm = dict((k, line[1][v]) for k, v in\
                 Record[line[0]].items() if k != "ID_pa")

class level_ENSDF:
    def __init__(self, line) -> None:
        self.transitions = []

        # information of initial state
        self.par = self.level_parse(line)


    def level_parse(self, line):
        if not re.match(Record["Level"]["ID_pa"], line[id_pos]):
            raise Exception("wrong line being parsed to level")

        return dict((k, line[v]) for k, v in Record["Level"].items() if k != "ID_pa")

    # add transition to the inital state
    def Add_Transition(self, line):
        if not re.match(Record[line[0]]["ID_pa"], line[1][id_pos]):
            raise Exception("wrong line being parsed to "+line[0])

        if line[1][5] == " ":
            self.transitions.append({"Trans_type": line[0]})
            self.transitions[-1].update(dict((k, line[1][v]) for k, v in\
                 Record[line[0]].items() if k != "ID_pa"))
            self.transitions[-1]["Extra"] = []
        elif self.transitions:
            self.transitions[-1]["Extra"].append([line[1][5], line[1][9:]])
        else:
            raise Exception("Somthing wrong during transition parse")
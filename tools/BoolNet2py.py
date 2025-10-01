import re
from rpy2.robjects import r
from rpy2.robjects.packages import importr

class BoolNetNetworkAnalyzer:
    def __init__(self, sbml_path):
        self.boolnet = importr('BoolNet')
        self.sbml_path = sbml_path
        self.net_loaded = False
        self.genes = []
        self.functions = None

    def load_network(self):
        r(f'net <- loadSBML("{self.sbml_path}")')
        self.net_loaded = True
        print(r('net'))
        self.genes = list(r('net$genes'))
        self.functions = r('net$functions')

    def print_genes(self):
        print("Genes:", self.genes)

    def print_functions(self):
        print("Transition functions:")
        if self.functions != r('NULL'):
            names = list(self.functions.names)
            values = list(self.functions)
            for k, v in zip(names, values):
                print(f"{k}: {v}")
        else:
            print("No Boolean functions found in the network.")

    def get_inputs(self):
        inputs = list(r('setdiff(net$genes, unlist(lapply(net$functions, all.vars)))'))
        print("Inputs (cues/stimuli):", inputs)
        return inputs

    def get_outputs(self):
        outputs = list(r('setdiff(net$genes, names(net$functions))'))
        print("Outputs (readouts):", outputs)
        return outputs

    def get_inhibitors(self):
        inhibitors = set()
        if self.functions != r('NULL'):
            for rule in list(self.functions):
                found = re.findall(r'!([A-Za-z0-9_]+)', rule)
                inhibitors.update(found)
        print("Inhibitors (appear as !X in rules):", sorted(inhibitors))
        return sorted(inhibitors)

    def save_bnet(self, output_path="data/apoptosis.bnet"):
        r(f'saveNetwork(net, "{output_path}")')
        print(f"Network saved as {output_path}")

# Example usage:
if __name__ == "__main__":
    analyzer = BoolNetNetworkAnalyzer("data/apoptosis.xml")
    analyzer.load_network()
    analyzer.print_genes()
    analyzer.print_functions()
    analyzer.get_inputs()
    analyzer.get_outputs()
    analyzer.get_inhibitors()
    analyzer.save_bnet()
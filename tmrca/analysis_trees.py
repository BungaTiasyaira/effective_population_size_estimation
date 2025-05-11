import tskit
from collections import defaultdict

class Mutation:

    def __init__(self, generation, haplotypes):
        self.generation = generation
        self.haplotypes = haplotypes
        self.color = '#FFFFFF'

    def get_generation(self):
        return self.generation

    def get_color(self):
        return self.color
    
    def set_color(self, color):
        self.color = color

    def __repr__(self):
        return f"Mutation {self.generation}"

class Haplotype:

    def __init__(self, idx, mutation):
        self.idx = idx
        self.mutation = mutation
        self.parent = None

    def set_parent(self, parent):
        self.parent = parent

    def get_color(self):
        if self.mutation is not None:
            return self.mutation.get_color()
        
        return "#808080"

    def __repr__(self):
        return f"Haplotype {self.idx} {self.get_color()}"

    def get_age(self):
        return 0

class Node:

    def __init__(self, age, link_color):
        self.left = None
        self.right = None
        self.parent = None
        self.age = age
        self.link_color = link_color

    def set_left(self, left):
        if self.left is not None:
            self.left.set_parent(None)
        if left is not None:
            left.set_parent(self)
        self.left = left

    
    def set_right(self, right):
        if self.right is not None:
            self.right.set_parent(None)
        if right is not None:
            right.set_parent(self)
        self.right = right

    def set_parent(self, parent):
        self.parent = parent

    def get_color(self):
        return self.link_color

    def get_age(self):
        return self.age

    def __repr__(self):
        return f"Node ({self.age} old) {self.get_color()}"

class Cluster:

    def __init__(self, root):
        self.root = root

    def get_haplotypes(self):
        haplotypes = []
        def visit(node):
            if isinstance(node, Haplotype):
                haplotypes.append(node)
                return node

            visit(node.left)
            visit(node.right)

        visit(self.root)

        return haplotypes

    def get_haplotypes_by_mutation(self):
        haplotypes = self.get_haplotypes()

        haplotypes_by_mutation = defaultdict(list)
        for haplotype in haplotypes:
            haplotypes_by_mutation[haplotype.mutation].append(haplotype)

        return dict(haplotypes_by_mutation)

    def get_mrca_to_node(self, other):
        def visit_down(node):
            if node == other:
                return True
            elif isinstance(node, Haplotype):
                return False

            return visit_down(node.left) or visit_down(node.right)

        def visit_up(node, prev):
            if (node.left == prev and visit_down(node.right)) or (node.right == prev and visit_down(node.left)):
                return node
            
            return visit_up(node.parent, node)

        return visit_up(self.root.parent, self.root)

            

    def __repr__(self):
        return f"Cluster {self.root.get_age()} {[haplotype.idx for haplotype in self.get_haplotypes()]}"

def get_mapping_of_haplotype_to_mutation(trees_path, txt_path):

    # === Load your .trees file ===
    ts = tskit.load(trees_path)  # <-- Replace with your file

    # === Load included individual IDs from a file ===
    with open(txt_path, "r") as f:
        included_individuals = [int(line.strip()) for line in f if line.strip().isdigit()]

    # === Define target position ===
    target_position = 35000

    # === Find the site at that position ===
    site = next((s for s in ts.sites() if int(s.position) == target_position), None)
    if site is None:
        print(f"No mutation found at position {target_position}")
        exit()

    tree = ts.at(site.position)

    # === Maps for mutation states and origins ===
    state_map = {}
    origin_node_map = {}

    ancestral_state = site.ancestral_state

    # Hack to ignore mutations of type 0
    broken_type_0_mutations = {}
    for mut in site.mutations:
        if all([m['slim_time'] < 5000 for m in mut.metadata['mutation_list']]):
            broken_type_0_mutations[mut.id] = mut

    for mut in site.mutations:
        if mut.id in broken_type_0_mutations:
            continue
        state_map[mut.node] = mut.derived_state
        origin_node_map[mut.node] = mut.node
        if mut.parent != tskit.NULL and mut.parent not in broken_type_0_mutations:
            parent_mut_node = ts.mutation(mut.parent).node
            origin_node_map[mut.node] = origin_node_map.get(parent_mut_node, parent_mut_node)

    # === Function to trace mutation state and origin ===
    def get_state_origin_and_time(node_id):
        current = node_id
        while current != tskit.NULL and current not in broken_type_0_mutations:
            if current in state_map:
                origin_node = origin_node_map[current]
                origin_time = ts.node(origin_node).time
                return state_map[current], origin_node, origin_time, origin_node_map[current]
            current = tree.parent(current)
        return ancestral_state, None, None, None

    # === Dictionary to group by generation and parent mutation ===
    generation_to_parents = defaultdict(lambda: defaultdict(list))

    # === Process and print haplotype info ===
    print(f"Mutation states at position {target_position}:\n")

    # Create a dictionary to map individual ID to order in the file
    individual_order = {ind_id: idx + 1 for idx, ind_id in enumerate(included_individuals)}

    for ind_id in included_individuals:
        try:
            individual = ts.individual(ind_id)
        except ValueError:
            print(f"Individual {ind_id} not found.")
            continue

        for node_id in individual.nodes:
            state, origin_node, origin_time, parent_id = get_state_origin_and_time(node_id)

            # Group by generation and parent mutation ID
            if origin_time is not None:
                generation_to_parents[origin_time][parent_id].append((state, node_id, ind_id))

            # Optional: print to screen
            print(f"Haplotype: {node_id}")
            print(f"From individual: {ind_id}")
            print(f"Mutation state of position {target_position}: {state}")
            print(f"Node where this mutation appeared first: {origin_node}")
            print(f"Generation where it appeared: {origin_time}")
            print(f"Parent mutation ID: {parent_id}")
            print("-" * 60)

    # Generate mapping
    haplotype_to_mutation = {}

    for generation in sorted(generation_to_parents.keys()):
        print(f"Generation where it appeared: {generation}\n")
        for parent_id, haplotypes in generation_to_parents[generation].items():
            print(f"Parent ID: {parent_id}\n")
            this_generation_haplotypes = [(individual_order[hap[2]] - 1) * 2 + (1 if  hap[1] % 2 == 1 else 0) for hap in haplotypes]
            mutation = Mutation(generation, this_generation_haplotypes)
            for haplotype in this_generation_haplotypes:
                haplotype_to_mutation[int(haplotype)] = mutation
            #out.write(", ".join([f"{hap[1]}({hap[2]})({individual_order[hap[2]]})({(individual_order[hap[2]] - 1) * 2 + (1 if  hap[1] % 2 == 1 else 0)})" for hap in haplotypes]) + "\n")
            print(this_generation_haplotypes)
        print("------")

    return haplotype_to_mutation
    

import tskit
from collections import defaultdict

def get_mapping_of_haplotype_to_generation(trees_path, txt_path):

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

    for mut in site.mutations:
        state_map[mut.node] = mut.derived_state
        origin_node_map[mut.node] = mut.node
        if mut.parent != tskit.NULL:
            parent_mut_node = ts.mutation(mut.parent).node
            origin_node_map[mut.node] = origin_node_map.get(parent_mut_node, parent_mut_node)

    # === Function to trace mutation state and origin ===
    def get_state_origin_and_time(node_id):
        current = node_id
        while current != tskit.NULL:
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
    haplotype_to_generation = {}

    for generation in sorted(generation_to_parents.keys()):
        print(f"Generation where it appeared: {generation}\n")
        for parent_id, haplotypes in generation_to_parents[generation].items():
            print(f"Parent ID: {parent_id}\n")
            this_generation_haplotypes = [(individual_order[hap[2]] - 1) * 2 + (1 if  hap[1] % 2 == 1 else 0) for hap in haplotypes]
            for haplotype in this_generation_haplotypes:
                haplotype_to_generation[int(haplotype)] = generation
            #out.write(", ".join([f"{hap[1]}({hap[2]})({individual_order[hap[2]]})({(individual_order[hap[2]] - 1) * 2 + (1 if  hap[1] % 2 == 1 else 0)})" for hap in haplotypes]) + "\n")
            print(this_generation_haplotypes)
        print("------")

    return haplotype_to_generation
    

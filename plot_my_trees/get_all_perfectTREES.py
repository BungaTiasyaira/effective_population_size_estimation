import os

output_file = "plot_my_trees/all_tree_paths.txt"

with open(output_file, "w") as out:
    for n in sorted(range(1, 11)):  # Ensure directories are processed in order
        dir_path = f"{n}/TREES"
        if not os.path.exists(dir_path):
            continue
        # Sort filenames inside each directory
        for filename in sorted(os.listdir(dir_path)):
            if filename.endswith(".trees"):
                full_path = os.path.join(dir_path, filename)
                out.write(full_path + "\n")

print(f"Saved sorted .trees paths to {output_file}")

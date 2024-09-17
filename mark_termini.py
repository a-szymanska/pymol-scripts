import numpy as np

def mark_termini(term):
    if term == "C":
        get_res = lambda arr: int(np.max(arr))
        not_c_val = -np.inf
        color = "tv_red"
    elif term == "N":
        get_res = lambda arr: int(np.min(arr))
        not_c_val = np.inf
        color = "limon"
    else:
        print("Chose C or N terminus")
    
    chains = cmd.get_chains("all")
    print(f"Chains found:\n{chains}")
    for name in chains:
        cmd.select(f"chain_sele", f"chain {name}")
        model = cmd.get_model("chain_sele")
        residues = [int(a.resi) if a.name[0] == 'C' else not_c_val for a in model.atom]
        term_residue = get_res(residues)
        cmd.select("terminus", f"chain {name} and resi {term_residue} and name CA")
        cmd.show("spheres", "terminus")
        cmd.color(f"{color}", "terminus")

    cmd.delete("selection chain_sele")
    cmd.delete("selection terminus")
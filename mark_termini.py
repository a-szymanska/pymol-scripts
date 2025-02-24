import numpy as np
import re

def mark_termini(name="all", term="C"):
    if term == "C":
        get_res = lambda arr: int(np.max(arr))
        not_c_val = -1
        color = "tv_red"
    elif term == "N":
        get_res = lambda arr: int(np.min(arr))
        not_c_val = 1
        color = "limon"
    else:
        print("Chose C or N terminus")
    
    chains = cmd.get_chains(name)
    for ch_name in chains:
        cmd.select(f"chain_sele", f"chain {ch_name} and name CA")
        model = cmd.get_model("chain_sele")
        not_c_val *= len(model.atom)
        residues = [int(re.findall(r'\d+', a.resi)[0]) if a.name[0] == 'C' else not_c_val for a in model.atom]
        if not len(residues):
            continue
        term_residue = get_res(residues)
        cmd.select("terminus", f"chain {ch_name} and resi {term_residue} and name CA")
        cmd.show("spheres", "terminus")
        cmd.color(f"{color}", "terminus")

    cmd.delete("selection chain_sele")
    cmd.delete("selection terminus")
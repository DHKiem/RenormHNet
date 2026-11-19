#!/usr/bin python3
'''
# Author: DHK
this code has two mode. 
jxbinder config --element Mn --cutoff 15 
jxbinder
'''

from collections import defaultdict
import numpy as np
import pandas as pd

import os
import re, sys, argparse
try:
    import tomllib            # Python 3.11+
except ModuleNotFoundError:
    import tomli as tomllib   # Python 3.9–3.10

def extract_spin_block(text):
    """Extract the data block following 'Total spin moment (muB)'."""
    pattern_spin = r"Total spin moment\s*\(muB\)(?:.*?\n){3}(.*?)(?:\n\s*\n|\Z)"
    match = re.search(pattern_spin, text, re.DOTALL)
    if not match:
        print("❌ Error: 'Total spin moment (muB)' block not found.")
        sys.exit(1)
    return match.group(1).strip()

def encode_elements(elements):
    element_to_int = {}
    encoded = []
    current_id = 1
    for e in elements:
        if e not in element_to_int:
            element_to_int[e] = current_id
            current_id += 1
        encoded.append(element_to_int[e])

    return encoded, element_to_int


#def get_info_from_out(filename, atom_indices, cutoff, jxposition):
def get_info_from_out(filename, target_element_name, cutoff, jxposition):
    """Return list of spin moments for given atom indices."""
    with open(filename, "r") as f:
        text = f.read()

    F = open("sym_param.toml", "w")

    match_systemname = re.search(r"System\.Name\s+(\S+)", text)

    if match_systemname:
        systemname = match_systemname.group(1)
        print(f"System name: {systemname}")
    else:
        systemname = "openmx"
        print(f"System name is not found. Set by default: {systemname}")

    F.write("system_name= \""+str(systemname)+"\"\n")
    match_unit = re.search(r"Atoms\.UnitVectors\.Unit\s+(\S+)", text)

    if match_unit:
        unit = match_unit.group(1)
        if unit != "Ang":
            print(f"⚠️ Warning: Atoms.UnitVectors.Unit is '{unit}', expected 'Ang'")
        else:
            print("✅ Unit is Ang")
    else:
        print("Atoms.UnitVectors.Unit is not found.")
        #sys.exit(1)

    # <Atoms.UnitVectors ... Atoms.UnitVectors> extracting
    pattern_lat_vec = r"<Atoms\.UnitVectors\s*(.*?)\s*Atoms\.UnitVectors>"

    match = re.search(pattern_lat_vec, text, re.DOTALL)

    if match:
        block = match.group(1)
        numbers = re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", block)
        numbers = list(map(float, numbers))
        vectors = [numbers[i:i+3] for i in range(0, len(numbers), 3)]
        print(vectors)
        F.write("lattice= "+str(vectors)+"\n")
    else:
        print("Atoms.UnitVectors not found.")


    data_block = extract_spin_block(text)

    spin_data = {}
    for line in data_block.splitlines():
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            idx = int(parts[0])
            spin = float(parts[-1])
            spin_data[idx] = spin
        except ValueError:
            continue

    spin_values = []


    FLOAT = r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?"
    ROW_RE = re.compile(
        rf"^\s*(\d+)\s+([A-Za-z][a-z]?)\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*$",
        re.MULTILINE
    )

    hdr = re.search(r"Fractional\s+coordinates\s+of\s+the\s+final\s+structure", text, re.IGNORECASE)
    if not hdr:
        raise RuntimeError("Header 'Fractional coordinates of the final structure' not found.")

    after = text[hdr.end():]
    xyz = []
    elements = []
    for m in ROW_RE.finditer(after):
        x, y, z = map(float, (m.group(3), m.group(4), m.group(5)))
        elements.append(m.group(2))
        xyz.append([x, y, z])
    if not xyz:
        raise RuntimeError("No coordinate rows found after the header.")
    print(elements)
    print(target_element_name)

    atom_indices = [i for i, x in enumerate(elements,1) if x in target_element_name]
    
    for idx in atom_indices:
        if idx in spin_data:
            spin_values.append(spin_data[idx])
        else:
            print(f"⚠️ Warning: Atom index {idx} not found.")
            spin_values.append(None)
    print(f"✅ Spin moments for atoms {spin_values}")
    F.write("atom1 = "+str(atom_indices)+"\n")
    F.write("atom2 = "+str(atom_indices)+"\n")
    upspins = []
    dnspins = []
    for atom,spin in zip(atom_indices,spin_values):
        if spin > 0 :
            upspins.append(atom)
        elif spin < 0 :
            dnspins.append(atom)
    F.write("upspins = "+str(upspins)+"\n")
    F.write("dnspins = "+str(dnspins)+"\n")



    numbers, mapping = encode_elements(elements)
    print(f"positions={xyz}")
    F.write("positions= "+str(xyz)+"\n")
    F.write("numbers = "+str(numbers)+"\n")
    #F.write("file_positions = [\"jx2.col.spin_0.0/\"]")
    F.write(f"file_positions = [\"{jxposition}\"]")
    

    if cutoff > 0.0:
        F.write("\ncutoffradius = "+str(cutoff)+"\n")

    F.close()


def binder():
    #############   INPUT START ###########

    with open("sym_param.toml", "rb") as f:
        params = tomllib.load(f)

    system_name = params.get("system_name", "openmx")
    lattice = np.array(params["lattice"])
    positions = params["positions"]
    atom1 = params["atom1"]
    atom2 = params["atom2"]
    upspins = params["upspins"]
    dnspins = params["dnspins"]
    numbers = params["numbers"]
    file_positions = params["file_positions"]
    outputfilename = params.get("output","Spin_Matrix.csv")
    cutoff = params.get("cutoffradius", None)

    #############   INPUT END ###########

    ### Input print
    print("Lattice information: \n", lattice)
    print("Atomic positions (fractional): \n", positions)

    # calculate distance from Rx, Ry, Rz
    def calculate_distance(row):
        a1 = int(row['Atom1']-1)
        a2 = int(row['Atom2']-1)
        atom1position = positions[a1][0]*lattice[0] + positions[a1][1]*lattice[1] + positions[a1][2]*lattice[2]
        atom2position = positions[a2][0]*lattice[0] + positions[a2][1]*lattice[1] + positions[a2][2]*lattice[2]
        position_vec = row['Rx'] *lattice[0] + row['Ry'] *lattice[1] + row['Rz'] *lattice[2] + (atom2position - atom1position)
        distance = np.linalg.norm(position_vec)
        return distance

    ### file load Jx calculation
    df = pd.DataFrame([])
    # file_positions[0] == "jx2.col.spin_0.0"
    # jxprefix = "jx2.col.spin"  # chempot = "0.0"
    jxdir = file_positions[0].rstrip("/")
    jxprefix, chempot = jxdir.rsplit("_", 1)
    # system_name == "openmx"
    # co_prefix = "jx2.col.spin_openmx_atomij_"
    co_prefix = jxprefix+"_"+system_name+"_atomij_"
    #["jx2.col.spin_0.0/"]
    #jx2.col.spin_openmx_atomij_3_3_[all_all]_ChemPdelta_0.0.csv
    for a1 in atom1:
      for a2 in atom2:
        jxfile = os.path.join(".", jxdir, co_prefix+str(a1)+"_"+str(a2)+"_[all_all]_ChemPdelta_"+chempot+".csv")
        jxfile21 = os.path.join(".", jxdir, co_prefix+str(a2)+"_"+str(a1)+"_[all_all]_ChemPdelta_"+chempot+".csv")
        fm = ((a1 in upspins) == (a2 in upspins))*2 -1 
        if os.path.exists(jxfile):
            Jfile = pd.read_csv(jxfile)
            if cutoff is not None:
                Jfile = Jfile[Jfile['Distance'] <= cutoff]

            df2 = pd.DataFrame({
                            'Rx': Jfile.Rx,
                            'Ry': Jfile.Ry,
                            'Rz': Jfile.Rz,
                            'Atom1': a1,
                            'Atom2': a2,
                            'J': Jfile.JmeV*fm,
                            })


        elif os.path.exists(jxfile21):
            Jfile = pd.read_csv(jxfile21)
            if cutoff is not None:
                Jfile = Jfile[Jfile['Distance'] <= cutoff]
            df2 = pd.DataFrame({
                            'Rx': -Jfile.Rx,
                            'Ry': -Jfile.Ry,
                            'Rz': -Jfile.Rz,
                            'Atom1': a1,
                            'Atom2': a2,
                            'J': Jfile.JmeV*fm,
                            })   
        else:
            assert False, f"File for atom pair {a1}, {a2} not found."


        df = pd.concat([df,df2], ignore_index=True)

    df['distance'] = df.apply(calculate_distance, axis=1)
    df = df[df['distance'] <= cutoff] if cutoff is not None else df
    df = df.sort_values(by=['distance','J'])
    df.to_csv(outputfilename, index=False)
    print(f"{outputfilename} has been saved.\n")

def main():
    parser = argparse.ArgumentParser(prog="jxbinder",
        description="Bind Jx outputs into Spin_Matrix.csv")
    sub = parser.add_subparsers(dest="mode")

    p_cfg = sub.add_parser("config", help="Generate sym_param.toml from OpenMX output")
    #p_cfg.add_argument("--atom", required=True, help="Comma-separated atom indices, e.g. 1,2,3,4")
    p_cfg.add_argument("--element", required=True, help="Comma-separated element names, e.g. Mn")
    p_cfg.add_argument("--file", default="openmx.out", help="Input file name (default: openmx.out)")
    p_cfg.add_argument("--cutoff", type=float, default=0.0, help="Cutoff value (default: 0.0)")
    p_cfg.add_argument("--jxprefix", default="jx2.col.spin_0.0/", help="default=\"jx2.col.spin_0.0/\"")

    args = parser.parse_args()

    if args.mode == "config":
        #atom_indices = [int(x) for x in args.atom.split(",") if x.strip().isdigit()]
        #get_info_from_out(args.file, atom_indices, args.cutoff, args.jxprefix)
        target_elements = [str(x) for x in args.element.split(",") ]
        get_info_from_out(args.file, target_elements, args.cutoff, args.jxprefix)
    else:
        binder()


if __name__ == "__main__":
    main()

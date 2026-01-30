from qiskit import QuantumCircuit
from qiskit.transpiler import CouplingMap, PassManager
from qiskit.transpiler.passes import (
    SabreLayout, SabreSwap,
    FullAncillaAllocation, EnlargeWithAncilla, ApplyLayout
)
import time 
import csv
import os

def load_qasm(path:str) -> QuantumCircuit:
    return QuantumCircuit.from_qasm_file(path)

def load_architecture(arch_path: str) -> CouplingMap:
    """
    .arch format:
      line 1: N (number of qubits)
      next lines: "u v" directed edges 
    """
    with open(arch_path, "r") as f:
        first = f.readline()
        n_qubits = int(first.strip())
        edges = []
        for line in f:
            line = line.strip()
            parts = line.split()
            u = int(parts[0])
            v = int(parts[1])
            edges.append((u, v))

    return n_qubits, CouplingMap(edges)

def count_swaps(qc: QuantumCircuit) -> int:
    return sum(1 for ci in qc.data if ci.operation.name == "swap")

def compile_with_lightsabre(qasm_path: str, arch_path: str, heuristic: str, seed:int) -> QuantumCircuit:
    qc = load_qasm(qasm_path)
    _, cmap = load_architecture(arch_path)

    # LightSABRE routing (heuristic: basic/lookahead/decay)
    sr = SabreSwap(coupling_map=cmap, heuristic=heuristic, seed=seed)

    # SabreLayout for choice initial layout (not modify circuit)
    sl = SabreLayout(coupling_map=cmap, routing_pass=sr, seed=seed)

    pm = PassManager([
        sl,
        FullAncillaAllocation(cmap),
        EnlargeWithAncilla(),
        ApplyLayout(),
        sr,
    ])

    return pm.run(qc)

def compile_with_budget(qasm_path: str, arch_path: str, heuristic: str, seed: int, budget: float):
    start = time.perf_counter()
    deadline = start + budget

    swaps = depth = size = None
    last_iter_s = None

    while True:
        now = time.perf_counter()

        # If we have an estimate, don't start another run if it likely won't fit.
        if last_iter_s is not None and now + last_iter_s >= deadline:
            break

        iter_start = time.perf_counter()
        final = compile_with_lightsabre(qasm_path, arch_path, heuristic, seed)
        iter_end = time.perf_counter()

        last_iter_s = iter_end - iter_start

        swaps = count_swaps(final)
        depth = final.depth()
        size = final.size()

        if iter_end >= deadline:
            break

    real_time = time.perf_counter() - start
    return swaps, depth, size, seed, real_time

def save_results(
    file_out: str,
    filename: str,
    architecture_name: str,
    seed: int,
    swaps_in: int,
    swaps_out: int,
    depth_in: int,
    depth_out: int,
    size_in: int,
    size_out: int,
    time_budget: float,
    real_time: float,
    heuristic: str,
):
    fieldnames = [
        "seed",
        "instance_file",
        "architecture",
        "heuristic",
        "time_budget",
        "real_time",
        "swaps_in",
        "swaps_out",
        "depth_in",
        "depth_out",
        "size_in",
        "size_out",
    ]

    row = {
        "seed": seed,
        "instance_file": filename,
        "architecture": architecture_name,
        "heuristic": heuristic,
        "time_budget": time_budget,
        "real_time": real_time,
        "swaps_in": swaps_in,
        "swaps_out": swaps_out,
        "depth_in": depth_in,
        "depth_out": depth_out,
        "size_in": size_in,
        "size_out": size_out,
    }

    write_header = (not os.path.exists(file_out)) or (os.path.getsize(file_out) == 0)

    with open(file_out, mode="a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    qasm = "revlib/examples"   # can be a file OR a folder
    arch = "revlib/architectures/ibmq_tokyo.arch"
    heuristics = ["basic", "lookahead", "decay"]
    seconds_list = [10]
    seeds = [1, 2, 3, 4, 5]
    file_out = "results_sabre.csv"

    # Build list of qasm files (file or folder)
    if os.path.isdir(qasm):
        qasm_files = [
            os.path.join(root, f)
            for root, _, files in os.walk(qasm)
            for f in files
            if f.lower().endswith(".qasm")
        ]
    else:
        qasm_files = [qasm]

    for qasm_file in qasm_files:
        initial = load_qasm(qasm_file)
        swaps_in = count_swaps(initial)
        depth_in = initial.depth()
        size_in = initial.size()

        print("\n====================")
        print("QASM:", qasm_file)
        print("=== INITIAL ===")
        print("depth:", depth_in)
        print("size :", size_in)
        print("swaps:", swaps_in)

        for h in heuristics:
            print(f"\n### Heuristic = {h} ###")
            for budget in seconds_list:
                for s in seeds:
                    swaps_out, depth_out, size_out, seed, real_time = compile_with_budget(
                        qasm_path=qasm_file,
                        arch_path=arch,
                        heuristic=h,
                        seed=s,
                        budget=budget
                    )

                    print(f"heuristic={h} budget={budget}s -> swaps={swaps_out}, depth={depth_out}, size={size_out}, seed={seed}, elapsed={real_time:.3f}s")

                    save_results(
                        file_out=file_out,
                        filename=qasm_file,
                        architecture_name=os.path.basename(arch),
                        seed=seed,
                        swaps_in=swaps_in,
                        swaps_out=swaps_out,
                        depth_in=depth_in,
                        depth_out=depth_out,
                        size_in=size_in,
                        size_out=size_out,
                        time_budget=budget,
                        real_time=real_time,
                        heuristic=h,
                    )



from qiskit import QuantumCircuit
from qiskit.transpiler import CouplingMap, PassManager
from qiskit.transpiler.passes import (
    SabreLayout, SabreSwap,
    FullAncillaAllocation, EnlargeWithAncilla, ApplyLayout
)
import time 

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

def compile_with_budget(qasm_path: str, arch_path: str, heuristic: str, seed: int, budget: int):
    start = time.perf_counter()
    
    while True:
        final = compile_with_lightsabre(qasm_path, arch_path, heuristic, seed)
        swaps = count_swaps(final)
        depth = final.depth()
        size = final.size()
 
        if time.perf_counter() - start >= budget:
            break

    return swaps, depth, size, seed


if __name__ == "__main__":
    qasm = "revlib/examples/urf5_280.qasm"
    arch = "revlib/architectures/ibmq_tokyo.arch"
    heuristics = ["basic", "lookahead", "decay"]
    seconds_list = [5]
    seed = [1,2,3,4]
    
    initial = load_qasm(qasm)
    print("=== INITIAL ===")
    print("depth:", initial.depth())
    print("size :", initial.size())
    print("swaps:", count_swaps(initial))
    
    for h in heuristics:
        print(f"\n### Heuristic = {h} ###")
        for budget in seconds_list:
            for s in seed:
                swaps, depth, size, seed_used = compile_with_budget(qasm_path=qasm, arch_path=arch, heuristic=h, seed=s, budget=budget)
                print(f"heuristic={h} budget={budget}s -> swaps={swaps}, depth={depth}, size={size}, seed={seed_used}")
        

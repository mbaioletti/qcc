import pandas as pd, numpy as np, subprocess, os
from pathlib import Path
from itertools import permutations
from concurrent.futures import ThreadPoolExecutor, as_completed

base_dir = Path("../revlib/examples/")
df = pd.read_csv("../revlib/dati2.csv")

MASTER_CSV = Path("esp_all.csv")       # final csv
TMP_DIR = Path("tmp_results")          # folder for temp results
TMP_DIR.mkdir(exist_ok=True)
MAX_WORKERS = 10

# map filename and column to dict
file2n = (df.assign(filename=df["name"].astype(str),num_qubits=df["num_qubits"].astype(int)).set_index("name")["num_qubits"].to_dict())

def generate_init_states(n: int):
    if n <= 4:
        for p in permutations(range(n)):
            yield ";".join(map(str, p))        
    else:
        seen = set()
        while len(seen) < 30:
            full = np.random.permutation(n)    
            arr = tuple(full[:])            
            if arr not in seen:
                seen.add(arr)
                yield ";".join(map(str, arr))

def execute(qasm_path: Path, init_state: str, seed: int, n: int, idx: int) -> Path:
    tmp_csv = TMP_DIR / f"esp_{qasm_path.stem}_n{n}_i{idx}_seed{seed}.csv"
    args = [
        "./compile_qasm",
        "../revlib/architectures/ibmq_tokyo.arch",
        str(qasm_path),
        "-timeout","10","-min_pr","0","-max_pr","0",
        "-objf","depth",
        "-res", str(tmp_csv),          
        "-seed", str(seed),
        "-init_state", init_state
    ]
    subprocess.run(args, check=True)
    return tmp_csv

def append_to_master(master_csv: Path, part_csvs: list[Path]):
    write_header = not master_csv.exists()
    with open(master_csv, "a", newline="") as _:
        pass 
    for p in part_csvs:
        if not p.exists() or p.stat().st_size == 0:
            continue
        dfp = pd.read_csv(p)
        dfp.to_csv(master_csv, mode="a", index=False, header=write_header)
        write_header = False 

def clean_partials(paths: list[Path]):
    for p in paths:
        try:
            p.unlink(missing_ok=True)
        except Exception:
            pass

if __name__ == "__main__":
    if MASTER_CSV.exists():
        MASTER_CSV.unlink()

    for qasm_path in sorted(base_dir.glob("*.qasm")):
        name = qasm_path.name
        
        n = file2n[name]
        print(f"==> {name} | n={n}")

        init_idx = 0
        for init_state in generate_init_states(n):
            init_idx += 1

            # run seed in parallel MAX_WORKERS
            with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
                futs = [ex.submit(execute, qasm_path, init_state, seed, n, init_idx)
                        for seed in range(1, 11)]
                parts = [f.result() for f in as_completed(futs)]

            append_to_master(MASTER_CSV, parts)
            clean_partials(parts)

        print(f"[OK] Complete {name}")

    print(f"Master CSV create in: {MASTER_CSV.resolve()}")

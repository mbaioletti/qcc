from pathlib import Path
import subprocess

COMPILER_BIN = "./compile_qasm"
ARCH_PATH = Path("../revlib/architectures/ibmq_tokyo.arch")
QASM_DIR = Path("../qcc_single/time-exp")   
RESULT_CSV = "time-test.csv"           
TIMEOUT = 60

def run_one(qasm_path: Path) -> int:
    cmd = [
        COMPILER_BIN,
        str(ARCH_PATH),
        str(qasm_path),
        "-timeout", str(TIMEOUT),
        "-trace", RESULT_CSV,
    ]
    print("->", " ".join(cmd), flush=True)
    res = subprocess.run(cmd)
    return res.returncode
          
if __name__ == "__main__":
    qasm_files = sorted(QASM_DIR.glob("*.qasm"))
    
    print(f"Found {len(qasm_files)} file .qasm in {QASM_DIR}")
    for i, q in enumerate(qasm_files, 1):
        print(f"\n=== ({i}/{len(qasm_files)}) {q.name} ===")
        run_one(q)
    
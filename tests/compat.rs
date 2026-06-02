use std::path::PathBuf;
use std::process::{Command, Stdio};

fn golden(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

/// (cli flags, skbio rooted=..., skbio weight=...) for each variant we check.
const VARIANTS: &[(&[&str], &str, &str)] = &[
    (&[], "None", "False"),
    (&["--rooted", "--weight", "1"], "True", "True"),
    (&["--unrooted", "--weight", "1"], "False", "True"),
    (&["--rooted", "--weight", "0.25"], "True", "0.25"),
    (&["--unrooted", "--weight", "0.5"], "False", "0.5"),
];

fn run_binary(table: &PathBuf, tree: &PathBuf, flags: &[&str]) -> Vec<(String, f64)> {
    let bin = env!("CARGO_BIN_EXE_rsomics-phydiv");
    let out = Command::new(bin)
        .arg(table)
        .arg("--tree")
        .arg(tree)
        .args(flags)
        .args(["-p", "12"])
        .output()
        .expect("run rsomics-phydiv");
    assert!(out.status.success(), "binary failed: {:?}", out);
    parse(&String::from_utf8(out.stdout).unwrap())
}

fn parse(text: &str) -> Vec<(String, f64)> {
    text.lines()
        .skip(1)
        .filter(|l| !l.is_empty())
        .map(|l| {
            let (s, v) = l.split_once('\t').unwrap();
            (s.to_string(), v.parse().unwrap())
        })
        .collect()
}

#[test]
fn matches_committed_skbio_golden() {
    let table = golden("table.tsv");
    let tree = golden("tree.nwk");
    for (flags, rooted, weight) in VARIANTS {
        let got = run_binary(&table, &tree, flags);
        let want = parse(
            &std::fs::read_to_string(golden(&format!("expected_{rooted}_{weight}.tsv"))).unwrap(),
        );
        assert_eq!(got.len(), want.len(), "{flags:?}");
        for ((gs, gv), (ws, wv)) in got.iter().zip(&want) {
            assert_eq!(gs, ws);
            assert!(
                (gv - wv).abs() < 1e-9,
                "{flags:?} {gs}: ours {gv} vs skbio {wv}"
            );
        }
    }
}

#[test]
fn matches_live_skbio() {
    let python = std::env::var("SKBIO_PYTHON").unwrap_or_else(|_| "python3".to_string());
    let probe = Command::new(&python).args(["-c", "import skbio"]).status();
    match probe {
        Ok(s) if s.success() => {}
        _ => {
            eprintln!("SKIP matches_live_skbio: scikit-bio not importable via '{python}'");
            return;
        }
    }

    let table = golden("table.tsv");
    let tree = golden("tree.nwk");
    for (flags, rooted, weight) in VARIANTS {
        let want = oracle(&python, &table, &tree, rooted, weight);
        let got = run_binary(&table, &tree, flags);
        assert_eq!(got.len(), want.len(), "{flags:?}");
        for ((gs, gv), (ws, wv)) in got.iter().zip(&want) {
            assert_eq!(gs, ws);
            assert!(
                (gv - wv).abs() < 1e-9,
                "{flags:?} {gs}: ours {gv} vs skbio {wv}"
            );
        }
    }
}

fn oracle(
    python: &str,
    table: &PathBuf,
    tree: &PathBuf,
    rooted: &str,
    weight: &str,
) -> Vec<(String, f64)> {
    let script = r#"
import sys
import numpy as np
from skbio import TreeNode
from skbio.diversity.alpha import phydiv
table, treepath, rooted, weight = sys.argv[1:5]
rooted = {"None": None, "True": True, "False": False}[rooted]
weight = {"False": False, "True": True}.get(weight, weight)
if isinstance(weight, str):
    weight = float(weight)
tree = TreeNode.read([open(treepath).read().strip()])
lines = [l for l in open(table).read().splitlines() if l and not l.startswith('#')]
samples = lines[0].split('\t')[1:]
names, rows = [], []
for l in lines[1:]:
    f = l.split('\t')
    names.append(f[0]); rows.append([int(x) for x in f[1:]])
counts = np.array(rows)
print("sample\tphydiv")
for j, s in enumerate(samples):
    v = phydiv(counts[:, j], taxa=names, tree=tree, rooted=rooted, weight=weight)
    print(f"{s}\t{v:.12f}")
"#;
    let mut child = Command::new(python)
        .args(["-c", script])
        .arg(table)
        .arg(tree)
        .arg(rooted)
        .arg(weight)
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .spawn()
        .expect("spawn skbio");
    drop(child.stdin.take());
    let out = child.wait_with_output().expect("skbio output");
    assert!(out.status.success(), "skbio failed: {:?}", out);
    parse(&String::from_utf8(out.stdout).unwrap())
}

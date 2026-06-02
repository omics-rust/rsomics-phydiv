use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_phydiv::{Config, Rooted, Weight, run};
use rsomics_phylo_tree::Tree;

fn balanced_newick(n_tips: usize) -> String {
    fn rec(lo: usize, hi: usize, out: &mut String) {
        if hi - lo == 1 {
            out.push_str(&format!("T{lo}:1.0"));
            return;
        }
        let mid = (lo + hi) / 2;
        out.push('(');
        rec(lo, mid, out);
        out.push(',');
        rec(mid, hi, out);
        out.push_str("):1.0");
    }
    let mut s = String::new();
    rec(0, n_tips, &mut s);
    s.push_str("root;");
    s
}

fn table(n_tips: usize, n_samples: usize) -> String {
    let mut s = String::from("feature");
    for j in 0..n_samples {
        s.push_str(&format!("\tS{j}"));
    }
    s.push('\n');
    for i in 0..n_tips {
        s.push_str(&format!("T{i}"));
        for j in 0..n_samples {
            let c = (i * 7 + j * 3) % 5;
            s.push('\t');
            s.push_str(&c.to_string());
        }
        s.push('\n');
    }
    s
}

fn bench(c: &mut Criterion) {
    let tree = Tree::from_newick(&balanced_newick(5000)).unwrap();
    let tbl = table(5000, 500);
    let cfg = Config {
        delim: '\t',
        rooted: Rooted::Unrooted,
        weight: Weight::Theta(0.5),
        precision: 6,
    };
    c.bench_function("phydiv_unrooted_theta_5000tips_500samples", |b| {
        b.iter(|| {
            let mut out = Vec::new();
            run(
                std::io::Cursor::new(black_box(&tbl)),
                &mut out,
                black_box(&tree),
                &cfg,
            )
            .unwrap();
            black_box(out);
        });
    });
}

criterion_group!(benches, bench);
criterion_main!(benches);
